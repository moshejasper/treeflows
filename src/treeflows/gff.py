from __future__ import annotations # for type checking...

#import numpy as np
import pandas as pd
import subprocess
import functools
from pathlib import Path
#import gzip
#from random import Random
#from cyvcf2 import VCF

from bisect import bisect_right
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple, Union


# going to try and get a full features database in... 

class GFF_Feature:

    def __init__(self, gff, featureline):
        """Placeholder for a parsed GFF feature record."""
        pass

class GFF_Attributes:
    def __init__(self, feature, attr_txt):
        """Placeholder for parsed GFF attribute text into a structured form."""
        pass

class GFF:

    def __init__(self, gff_file):
        """Load a GFF file into a simple header string and feature table (prototype)."""
        self.header = ""
        self.feature_table = pd.DataFrame(columns=["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"])

        with open(gff_file, 'r') as gff:
            for ll in gff:
                if ll.startswith("#"): 
                    self.header += ll # keep \n presences
                else:
                    print(ll)
                    self.feature_table = pd.concat([self.feature_table, pd.Series(ll.strip().split('\t'))], axis=1)
                    print(self.feature_table)

        print(self.header)
        print(self.feature_table)

    pass


def ncbi_fetch_gff(to_file, target_acc="NC_035107.1,NC_035108.1,NC_035109.1"):
    """Download a GFF3 from NCBI's sviewer endpoint using `wget`.

    Args:
        to_file: Output filename to write.
        target_acc: Comma-separated list of NCBI accessions/contigs to request.

    Returns:
        None. Writes `to_file`.
    """

    subprocess.run(f'wget -O {to_file} "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id={target_acc}"', shell=True)

def gff_extract_gene(gene, infile, outfile):
    """Extract records for a gene from a GFF3 file (placeholder/incomplete).

    Args:
        gene: Gene identifier to match.
        infile: Input GFF3 filename.
        outfile: Output filename to write matching lines to.

    Returns:
        None.

    Notes:
        This function currently sets up an iteration skeleton but does not yet
        implement the full extraction logic.
    """

    with open(infile) as ifi, open(outfile, 'w') as ofi:

        active = False

        for ll in ifi:
            if ll.startswith('#'):
                ofi.write(ll)
                continue


awks="""
awk -F'\t' -v target="LOC5567355" '
function has_gene_field(attr,   n,i,kv,v) {
    n = split(attr, kv, /[;]+/)                  # split col6 on ;
    for (i=1; i<=n; i++) {
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", kv[i])
        if (kv[i] ~ /^gene=/) {
            v = kv[i]; sub(/^gene=/, "", v)      # strip key
            # strip optional quotes
            if (v ~ /^'\''.*'\''$/) v = substr(v, 2, length(v)-2)
            else if (v ~ /^".*"$/) v = substr(v, 2, length(v)-2)
            if (v == target) return 1
        }
    }
    return 0
}
/^#/ { print; next }                             # 1) comments
$3=="gene" && has_gene_field($9) {               # 2) anchor line
    print; seen=1; next
}
seen && has_gene_field($9) {                     # 3) subsequent matches
    print; next
}
' input.tsv > filtered.tsv

"""

def get_max_overlap(infile, feature):
    """Merge overlapping intervals of a given feature type from a GFF-like file.

    Args:
        infile: Input file path (tab-delimited GFF/GFF3 expected).
        feature: Feature type to extract (column 3 / index 2).

    Returns:
        List of merged `(start, end)` integer intervals.
    """

    intervals = []

    with open(infile) as ifi:
        for ll in ifi:
            if ll.startswith('#'):
                continue
            lparsed = ll.strip().split('\t')
            if lparsed[2] != feature:
                continue
            intervals.append((int(lparsed[3]), int(lparsed[4])))

    print(intervals)
    
    intervals = sorted(intervals, key=lambda x: (x[0], x[1]))
    print(intervals)
    merged = []

    if not intervals:
        return merged
    
    cur_left, cur_right = intervals[0]
    for left, right in intervals[1:]:
        overlaps = left <= cur_right + 1
        if overlaps:
            cur_right = max(cur_right, right)
        else:
            merged.append((cur_left, cur_right))
            cur_left, cur_right = left, right

    merged.append((cur_left, cur_right))
    print(merged)
    return merged



###  GFF3 MASKS


Interval = Tuple[int, int]          # 1-based, inclusive
MaskDict = Dict[str, List[Interval]]
PathLike = Union[str, Path]


class GenomeMask:
    """
    Interval mask for genomic coordinates.

    Internal coordinate convention:
        - 1-based, inclusive intervals: (start, end)

    Intended use:
        - Build once from a file or interval dictionary
        - Query repeatedly with pos_in_mask()

    Constructors:
        - GenomeMask.from_gff3(...)
        - GenomeMask.from_bed(...)
        - GenomeMask.from_repeatmasker_out(...)
        - GenomeMask.from_intervals(...)

    Notes:
        - GFF3 is read as 1-based inclusive
        - BED is converted from 0-based half-open to 1-based inclusive
        - RepeatMasker .out is treated as 1-based inclusive
    """

    def __init__(self, masks: Mapping[str, Sequence[Interval]]):
        self._masks: MaskDict = self._normalize_and_merge(masks)
        self._starts: Dict[str, List[int]] = self._build_starts(self._masks)

    # ----------------------------
    # Public alternate constructors
    # ----------------------------

    @classmethod
    def from_gff3(
        cls,
        gff3_file: PathLike,
        feature_types: Optional[Union[str, Iterable[str]]] = None,
        applications: Optional[Union[str, Iterable[str]]] = ("RepeatMasker",),
        chromosomes: Optional[Union[str, Iterable[str]]] = None,
        pos_min: Optional[int] = None,
        pos_max: Optional[int] = None,
    ) -> "GenomeMask":
        masks = cls._load_gff3_mask(
            gff3_file=gff3_file,
            feature_types=feature_types,
            applications=applications,
            chromosomes=chromosomes,
            pos_min=pos_min,
            pos_max=pos_max,
        )
        return cls(masks)

    @classmethod
    def from_bed(
        cls,
        bed_file: PathLike,
        chromosomes: Optional[Union[str, Iterable[str]]] = None,
        pos_min: Optional[int] = None,
        pos_max: Optional[int] = None,
    ) -> "GenomeMask":
        masks = cls._load_bed_mask(
            bed_file=bed_file,
            chromosomes=chromosomes,
            pos_min=pos_min,
            pos_max=pos_max,
        )
        return cls(masks)

    @classmethod
    def from_repeatmasker_out(
        cls,
        out_file: PathLike,
        chromosomes: Optional[Union[str, Iterable[str]]] = None,
        pos_min: Optional[int] = None,
        pos_max: Optional[int] = None,
    ) -> "GenomeMask":
        masks = cls._load_repeatmasker_out_mask(
            out_file=out_file,
            chromosomes=chromosomes,
            pos_min=pos_min,
            pos_max=pos_max,
        )
        return cls(masks)

    @classmethod
    def from_intervals(
        cls,
        intervals: Mapping[str, Sequence[Interval]],
        chromosomes: Optional[Union[str, Iterable[str]]] = None,
        pos_min: Optional[int] = None,
        pos_max: Optional[int] = None,
    ) -> "GenomeMask":
        chrom_set = cls._as_set(chromosomes)
        out: MaskDict = {}

        for chrom, ivals in intervals.items():
            if chrom_set is not None and chrom not in chrom_set:
                continue
            for start, end in ivals:
                clipped = cls._clip_interval(start, end, pos_min, pos_max)
                if clipped is None:
                    continue
                out.setdefault(chrom, []).append(clipped)

        return cls(out)

    # ----------------------------
    # Public accessors
    # ----------------------------

    @property
    def masks(self) -> MaskDict:
        return self._masks

    @property
    def starts(self) -> Dict[str, List[int]]:
        return self._starts

    def intervals(self, chrom: str) -> List[Interval]:
        return self._masks.get(chrom, [])

    def chromosomes(self) -> List[str]:
        return sorted(self._masks.keys())

    def pos_in_mask(self, chrom: str, pos: int) -> bool:
        """
        Test whether a 1-based position falls inside any masked interval.
        """
        intervals = self._masks.get(chrom)
        if not intervals:
            return False

        i = bisect_right(self._starts[chrom], pos) - 1
        if i < 0:
            return False

        return pos <= intervals[i][1]

    def span_in_mask(self, chrom: str, start: int, end: int) -> bool:
        """
        Test whether a 1-based inclusive span overlaps any masked interval.
        """
        intervals = self._masks.get(chrom)
        if not intervals:
            return False

        i = bisect_right(self._starts[chrom], end) - 1
        if i < 0:
            return False

        s0, e0 = intervals[i]
        if e0 >= start:
            return True

        j = i + 1
        return j < len(intervals) and intervals[j][0] <= end

    def fraction_masked(self, chrom: str, start: int, end: int) -> float:
        """
        Fraction of a 1-based inclusive interval covered by the mask.
        """
        if end < start:
            raise ValueError("end must be >= start")

        intervals = self._masks.get(chrom)
        if not intervals:
            return 0.0

        total_len = end - start + 1
        covered = 0

        i = bisect_right(self._starts[chrom], start) - 1
        if i < 0:
            i = 0

        while i < len(intervals):
            s, e = intervals[i]
            if s > end:
                break
            ov_start = max(start, s)
            ov_end = min(end, e)
            if ov_start <= ov_end:
                covered += (ov_end - ov_start + 1)
            i += 1

        return covered / total_len

    def to_bed(self, bed_file: PathLike) -> None:
        """
        Write mask as BED (0-based, half-open).
        """
        with open(bed_file, "w") as out:
            for chrom in sorted(self._masks):
                for start, end in self._masks[chrom]:
                    out.write(f"{chrom}\t{start - 1}\t{end}\n")

    # ----------------------------
    # Internal loaders
    # ----------------------------

    @classmethod
    def _load_gff3_mask(
        cls,
        gff3_file: PathLike,
        feature_types: Optional[Union[str, Iterable[str]]] = None,
        applications: Optional[Union[str, Iterable[str]]] = ("RepeatMasker",),
        chromosomes: Optional[Union[str, Iterable[str]]] = None,
        pos_min: Optional[int] = None,
        pos_max: Optional[int] = None,
    ) -> MaskDict:
        feature_set = cls._as_set(feature_types)
        app_set = cls._as_set(applications)
        chrom_set = cls._as_set(chromosomes)

        masks: MaskDict = {}

        with open(gff3_file) as f:
            for line in f:
                if not line or line.startswith("#"):
                    continue

                fields = line.rstrip().split("\t")
                if len(fields) < 5:
                    continue

                chrom = fields[0]
                source = fields[1]
                feature_type = fields[2]
                start = int(fields[3])   # GFF3: 1-based inclusive
                end = int(fields[4])

                if chrom_set is not None and chrom not in chrom_set:
                    continue
                if app_set is not None and source not in app_set:
                    continue
                if feature_set is not None and feature_type not in feature_set:
                    continue

                clipped = cls._clip_interval(start, end, pos_min, pos_max)
                if clipped is None:
                    continue

                masks.setdefault(chrom, []).append(clipped)

        return masks

    @classmethod
    def _load_bed_mask(
        cls,
        bed_file: PathLike,
        chromosomes: Optional[Union[str, Iterable[str]]] = None,
        pos_min: Optional[int] = None,
        pos_max: Optional[int] = None,
    ) -> MaskDict:
        chrom_set = cls._as_set(chromosomes)
        masks: MaskDict = {}

        with open(bed_file) as f:
            for line in f:
                if not line or line.startswith(("#", "track", "browser")):
                    continue

                fields = line.rstrip().split("\t")
                if len(fields) < 3:
                    continue

                chrom = fields[0]
                if chrom_set is not None and chrom not in chrom_set:
                    continue

                bed_start = int(fields[1])   # BED: 0-based
                bed_end = int(fields[2])     # BED: end-exclusive

                start = bed_start + 1        # convert to 1-based inclusive
                end = bed_end

                clipped = cls._clip_interval(start, end, pos_min, pos_max)
                if clipped is None:
                    continue

                masks.setdefault(chrom, []).append(clipped)

        return masks

    @classmethod
    def _load_repeatmasker_out_mask(
        cls,
        out_file: PathLike,
        chromosomes: Optional[Union[str, Iterable[str]]] = None,
        pos_min: Optional[int] = None,
        pos_max: Optional[int] = None,
    ) -> MaskDict:
        """
        Read RepeatMasker .out file.

        Expected columns after the header include:
            SW score, perc div., perc del., perc ins., query sequence,
            begin, end, (left), ...

        We extract:
            chrom = column 5
            start = column 6
            end   = column 7

        These are treated as 1-based inclusive.
        """
        chrom_set = cls._as_set(chromosomes)
        masks: MaskDict = {}

        with open(out_file) as f:
            for line in f:
                striped = line.strip()
                if not striped:
                    continue

                # Skip the standard header lines
                if (
                    striped.startswith("SW")
                    or striped.startswith("score")
                    or striped.startswith("perc")
                ):
                    continue

                fields = striped.split()
                if len(fields) < 7:
                    continue

                # Most RepeatMasker .out records begin with numeric SW score
                if not fields[0].isdigit():
                    continue

                chrom = fields[4]
                if chrom_set is not None and chrom not in chrom_set:
                    continue

                start = int(fields[5])
                end = int(fields[6])

                clipped = cls._clip_interval(start, end, pos_min, pos_max)
                if clipped is None:
                    continue

                masks.setdefault(chrom, []).append(clipped)

        return masks

    # ----------------------------
    # Internal utilities
    # ----------------------------

    @staticmethod
    def _as_set(x: Optional[Union[str, Iterable[str]]]) -> Optional[set]:
        if x is None:
            return None
        if isinstance(x, str):
            return {x}
        return set(x)

    @staticmethod
    def _clip_interval(
        start: int,
        end: int,
        pos_min: Optional[int],
        pos_max: Optional[int],
    ) -> Optional[Interval]:
        if end < start:
            return None

        if pos_min is not None and end < pos_min:
            return None
        if pos_max is not None and start > pos_max:
            return None

        if pos_min is not None:
            start = max(start, pos_min)
        if pos_max is not None:
            end = min(end, pos_max)

        if end < start:
            return None

        return (start, end)

    @staticmethod
    def _normalize_and_merge(masks: Mapping[str, Sequence[Interval]]) -> MaskDict:
        out: MaskDict = {}

        for chrom, intervals in masks.items():
            clean = [(int(s), int(e)) for s, e in intervals if e >= s]
            if not clean:
                continue

            clean.sort()
            merged: List[List[int]] = []

            for s, e in clean:
                if not merged or s > merged[-1][1] + 1:
                    merged.append([s, e])
                else:
                    merged[-1][1] = max(merged[-1][1], e)

            out[chrom] = [(s, e) for s, e in merged]

        return out

    @staticmethod
    def _build_starts(masks: Mapping[str, Sequence[Interval]]) -> Dict[str, List[int]]:
        return {chrom: [s for s, _ in ivals] for chrom, ivals in masks.items()}


# ----------------------------
# Example usage
# ----------------------------

if __name__ == "__main__":
    # From GFF3
    mask1 = GenomeMask.from_gff3(
        "repeats.gff3",
        applications="RepeatMasker",
        feature_types=None,
    )

    # From BED
    mask2 = GenomeMask.from_bed("repeats.bed")

    # From RepeatMasker .out
    mask3 = GenomeMask.from_repeatmasker_out("genome.fa.out")

    # From an in-memory dictionary
    mask4 = GenomeMask.from_intervals({
        "chr1": [(100, 200), (180, 250), (1000, 1100)],
        "chr2": [(50, 60)],
    })

    print(mask1.pos_in_mask("chr1", 150))
    print(mask4.pos_in_mask("chr1", 900))
    print(mask4.span_in_mask("chr1", 240, 260))
    print(mask4.fraction_masked("chr1", 150, 300))