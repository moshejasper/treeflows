import pyfaidx
import gzip
from treeflows.Spath import Mpath
from treeflows.bgzip import BgzfWriter
import numpy as np


def extract_region(infasta, outfasta, seqname, start, end, outwidth=40):
    '''Extracts region to new fasta file (inclusive). Note bp will be adjusted. need correction later. 
    
    VCF REFS HERE ARE 1-BASED AND INCLUSIVE'''

    # note. selection coords are 0-based, but then sequence names are 1-based (but not slicing)... 
    # start and stop will be 1-based here (per fasta type). And inclusive!!! Means we need conversion. 

    genome = pyfaidx.Fasta(infasta)

    with open(outfasta, 'w') as ofi:

        rawregion = genome[seqname][start-1:end]
        ofi.write(f'>{seqname}:{rawregion.start}-{rawregion.end}\n')
        chunks = len(rawregion) // outwidth # NOTE: CODE WAS BADLY WRITTEN AND IS NOW BENG REWRITTEN...
        for i in range(chunks): # this slicing is based on a string, not original sequence. 
            ofi.write(f'{rawregion[i*outwidth:(i+1)*outwidth]}\n')
        ofi.write(f'{rawregion[chunks*outwidth:]}')


def exclude_region(infasta, outfasta, blacklist, outwidth=40): # must have set this to 80 somewhere
    '''Extracts all regions but one to a new fasta file'''

    genome = pyfaidx.Fasta(infasta)

    with open(outfasta, 'w') as ofi: 

        for scaffold in genome.keys():
            if not scaffold in blacklist:
                scaf = genome[scaffold]
                ofi.write(f'>{scaf.name}\n')
                seglength = len(scaf)
                chunks = seglength // outwidth # NOTE: CODE HERE WAS BADLY WRITTEN AND HAS JUST BEEN FIXED... 
                for i in range(chunks):
                    ofi.write(f'{scaf[i*outwidth:(i+1)*outwidth]}\n')
                ofi.write(f'{scaf[chunks*outwidth:]}\n')




### This section will be for handling mappability fasta map interpolation and aggregation. 

def fasta_interpolate_mappability(infasta, outfasta, period=5, outwidth=120):
    '''Interpolates mappability fasta map to a new fasta file with a different period.'''

    genome = pyfaidx.Fasta(infasta)

    if outwidth % period:
        raise ValueError('Output width must be a multiple of the period.')
    if not period % 2:
        raise ValueError('Period must be odd.')

    with gzip.open(outfasta, 'wt') as ofi:

        for scaffold in genome.keys():
            scaf = genome[scaffold]
            ofi.write(f'>{scaf.name}\n')
            seglength = len(scaf)
            chunks = seglength // outwidth # NOTE: CODE HERE WAS BADLY WRITTEN AND HAS JUST BEEN FIXED... 
            for i in range(chunks):
                ofi.write(f'{fasta_interpolate_line(scaf, i*outwidth, outwidth, period)}\n')
                #ofi.write(f'{scaf[i*outwidth:(i+1)*outwidth]}\n')
            ofi.write(f'{fasta_interpolate_line(scaf, i*outwidth, outwidth, period)}\n')
            #ofi.write(f'{scaf[chunks*outwidth:]}\n')

def fasta_interpolate_line(chrom, startpos, outwidth=80, period=5):
    '''Interpolates a line of mappability fasta map to a new line with a different period.'''

    hemi = period // 2 + 1
    line = chrom[startpos:(startpos+outwidth + 1 if startpos+outwidth < len(chrom) else len(chrom))]

    outline = ''
    for i in range(0, len(line), period):
        if not i % period:
            outline += str(line[i])*hemi
    
    return outline[:outwidth]


def fasta_mappability_compress(
    infasta,
    outfasta=None,
    period=5,
    outwidth=120,
    preserve_widths=False,
    read_chunk_bases=5_000_000,
    write_buffer_chars=1_000_000,
    bgzf_compresslevel=1,
):
    """Compress mappability FASTA to a smaller-period FASTA."""

    infasta = Mpath(infasta)
    if ".fa" not in infasta.suffixes:
        if len(infasta.suffixes) > 1:
            raise ValueError("Input fasta file must have a .fa extension.")
        elif infasta.with_suffix(".fa").exists():
            infasta = infasta.with_suffix(".fa")
        elif infasta.with_suffix(".fa.gz").exists():
            infasta = infasta.with_suffix(".fa.gz")
        else:
            raise ValueError("Input fasta file must have a .fa extension.")

    if outfasta is None:
        outfasta = infasta.with_suffixes(".compressed.fa.gz")

    if period % 2 == 0:
        raise ValueError("Period must be odd.")

    if preserve_widths:
        if outwidth % period:
            raise ValueError(
                "Output width must be a multiple of the period if preserve_widths is True."
            )
        outwidth //= period

    genome = pyfaidx.Fasta(infasta)

    with BgzfWriter(outfasta, compresslevel=bgzf_compresslevel, encoding="utf-8") as ofi:
        for scaffold in genome.keys():
            scaf = genome[scaffold]
            seglength = len(scaf)

            ofi.write(
                f">{scaf.name}\ttype=mappability_compressed"
                f"\tscaf_length={seglength}\tscaf_period={period}\n"
            )

            write_buf = []
            write_buf_len = 0
            current_line = []
            current_line_len = 0

            for block_start in range(0, seglength, read_chunk_bases):
                block_end = min(block_start + read_chunk_bases, seglength)
                block = str(scaf[block_start:block_end])

                offset = (-block_start) % period
                compressed = block[offset::period]

                j = 0
                while j < len(compressed):
                    need = outwidth - current_line_len
                    piece = compressed[j:j + need]
                    current_line.append(piece)
                    current_line_len += len(piece)
                    j += len(piece)

                    if current_line_len == outwidth:
                        line = "".join(current_line)
                        write_buf.append(line)
                        write_buf.append("\n")
                        write_buf_len += len(line) + 1
                        current_line.clear()
                        current_line_len = 0

                        if write_buf_len >= write_buffer_chars:
                            ofi.write("".join(write_buf))
                            write_buf.clear()
                            write_buf_len = 0

            if current_line:
                line = "".join(current_line)
                write_buf.append(line)
                write_buf.append("\n")

            if write_buf:
                ofi.write("".join(write_buf))



def fasta_mappability_aggregate(
    *args,
    outfasta=None,
    consensus_threshold=0.75,
    outwidth=120,
    chromosomes=None,
    read_chunk_bases=5_000_000,
    write_buffer_chars=1_000_000,
    bgzf_compresslevel=1,
):
    """Aggregate multiple mappability FASTAs to a single consensus FASTA.

    Consensus rule:
      - C if frac(C) >= threshold
      - else A if frac(A or C) >= threshold
      - else N
    """

    if len(args) == 1 and isinstance(args[0], (list, tuple)):
        args = args[0]

    if len(args) == 0:
        raise ValueError("At least one input fasta must be provided.")
    
    if not (0 < consensus_threshold <= 1):
        raise ValueError("consensus_threshold must be in (0, 1].")

    # normalize input paths properly
    infastas = []
    for arg in args:
        infasta = Mpath(arg)
        if ".fa" not in infasta.suffixes:
            if len(infasta.suffixes) > 1:
                raise ValueError("Input fasta file must have a .fa extension.")
            elif infasta.with_suffix(".fa").exists():
                infasta = infasta.with_suffix(".fa")
            elif infasta.with_suffix(".fa.gz").exists():
                infasta = infasta.with_suffix(".fa.gz")
            else:
                raise ValueError("Input fasta file must have a .fa extension.")
        infastas.append(infasta)

    if outfasta is None:
        outfasta = infastas[0].with_suffixes(".aggregated.fa.gz")
    outfasta = Mpath(outfasta)
    if len(outfasta.suffixes) == 0:
        outfasta = outfasta.with_suffix(".aggregated.fa.gz")

    genomes = [pyfaidx.Fasta(infasta) for infasta in infastas]

    ref_keys = list(genomes[0].keys())
    ref_key_set = set(ref_keys)

    if chromosomes is not None:
        ref_keys = [key for key in ref_keys if key in chromosomes]
        ref_key_set = set(ref_keys)

    for gi, genome in enumerate(genomes[1:], start=1):
        genome_keys = list(genome.keys())
        genome_key_set = set(genome_keys)

        if (chromosomes is not None and not ref_key_set.issubset(genome_key_set)) or (chromosomes is None and genome_key_set != ref_key_set):
            missing = sorted(ref_key_set - genome_key_set)
            extra = sorted(genome_key_set - ref_key_set) if chromosomes is None else []
            raise ValueError(
                f"Genome {gi} has different scaffold names from genome 0.\n"
                f"Missing from genome {gi}: {missing[:10]}\n"
                f"Extra in genome {gi}: {extra[:10]}"
            )
            

    nfiles = len(genomes)
    c_thresh = int(np.ceil(consensus_threshold * nfiles))
    ac_thresh = int(np.ceil(consensus_threshold * nfiles))

    # choose compact integer type for counts
    if nfiles <= np.iinfo(np.uint8).max:
        count_dtype = np.uint8
    elif nfiles <= np.iinfo(np.uint16).max:
        count_dtype = np.uint16
    else:
        count_dtype = np.uint32

    with BgzfWriter(outfasta, compresslevel=bgzf_compresslevel, encoding="utf-8") as ofi:
        for scaffold in ref_keys:
            scafs = [genome[scaffold] for genome in genomes]
            seglength = len(scafs[0])

            for si, scaf in enumerate(scafs[1:], start=1):
                if len(scaf) != seglength:
                    raise ValueError(
                        f"Scaffold {scaffold!r} has different lengths across inputs "
                        f"(genome 0: {seglength}, genome {si}: {len(scaf)})."
                    )

            ofi.write(f">{genomes[0][scaffold].long_name}\n")

            write_buf = []
            write_buf_len = 0
            current_line = bytearray()

            for block_start in range(0, seglength, read_chunk_bases):
                block_end = min(block_start + read_chunk_bases, seglength)
                block_len = block_end - block_start

                c_counts = np.zeros(block_len, dtype=count_dtype)
                ac_counts = np.zeros(block_len, dtype=count_dtype)

                for scaf in scafs:
                    # Read chunk once as bytes
                    block = str(scaf[block_start:block_end]).upper().encode("ascii")
                    arr = np.frombuffer(block, dtype=np.uint8)

                    # Optional validation. Remove if you want max speed.
                    bad = ~((arr == 65) | (arr == 67) | (arr == 78))  # A/C/N
                    if bad.any():
                        bad_pos = int(np.flatnonzero(bad)[0])
                        bad_char = chr(arr[bad_pos])
                        raise ValueError(
                            f"Unexpected character {bad_char!r} in scaffold {scaffold!r} "
                            f"at position {block_start + bad_pos}."
                        )

                    is_c = (arr == 67)
                    is_ac = (arr != 78)   # since only A/C/N allowed, non-N means A or C

                    c_counts += is_c
                    ac_counts += is_ac

                # Build consensus block as bytes
                out = np.full(block_len, 78, dtype=np.uint8)      # 'N'
                mask_a = ac_counts >= ac_thresh
                out[mask_a] = 65                                  # 'A'
                mask_c = c_counts >= c_thresh
                out[mask_c] = 67                                  # 'C'

                out_bytes = out.tobytes()

                # Wrap output lines
                j = 0
                while j < block_len:
                    need = outwidth - len(current_line)
                    piece = out_bytes[j:j + need]
                    current_line.extend(piece)
                    j += len(piece)

                    if len(current_line) == outwidth:
                        write_buf.append(bytes(current_line))
                        write_buf.append(b"\n")
                        write_buf_len += outwidth + 1
                        current_line.clear()

                        if write_buf_len >= write_buffer_chars:
                            ofi.write(b"".join(write_buf))
                            write_buf.clear()
                            write_buf_len = 0

            if current_line:
                write_buf.append(bytes(current_line))
                write_buf.append(b"\n")

            if write_buf:
                ofi.write(b"".join(write_buf))

def fasta_mappability_decompress(
    infasta,
    outfasta=None,
    outwidth=120,
    read_chunk_bases=5_000_000,
    write_buffer_chars=1_000_000,
    bgzf_compresslevel=1,
):
    """
    Decompress mappability FASTA to original scaffold length.

    Assumes each compressed character represents the midpoint of a period-sized
    interpolation block:
      - first compressed base expands to hemi = period // 2 + 1
      - subsequent compressed bases expand to period
      - final block is truncated to scaf_length if needed

    Guarantees output length == scaf_length from header.
    Pads with 'N' if decompressed sequence is too short.
    """

    infasta = Mpath(infasta)
    if ".fa" not in infasta.suffixes:
        if len(infasta.suffixes) > 1:
            raise ValueError("Input fasta file must have a .fa extension.")
        elif infasta.with_suffix(".fa").exists():
            infasta = infasta.with_suffix(".fa")
        elif infasta.with_suffix(".fa.gz").exists():
            infasta = infasta.with_suffix(".fa.gz")
        else:
            raise ValueError("Input fasta file must have a .fa extension.")

    if outfasta is None:
        outfasta = infasta.with_suffixes(".decompressed.fa.gz")

    genome = pyfaidx.Fasta(infasta)

    with BgzfWriter(outfasta, compresslevel=bgzf_compresslevel, encoding="utf-8") as ofi:
        for scaffold in genome.keys():
            scaf = genome[scaffold]
            header = scaf.long_name.split("\t")
            header_new = header[0]
            scaf_length = int(header[2].split("scaf_length=")[1])
            period = int(header[3].split("scaf_period=")[1])
            hemi = period // 2 + 1

            ofi.write(f">{header_new}\n")

            write_buf = []
            write_buf_len = 0
            current_line = []

            produced = 0
            comp_index = 0

            def push_seq(seq):
                nonlocal write_buf, write_buf_len, current_line

                j = 0
                while j < len(seq):
                    curr_len = sum(len(x) for x in current_line)
                    need = outwidth - curr_len
                    piece = seq[j:j + need]
                    current_line.append(piece)
                    j += len(piece)

                    curr_len += len(piece)
                    if curr_len == outwidth:
                        line = "".join(current_line)
                        write_buf.append(line)
                        write_buf.append("\n")
                        write_buf_len += len(line) + 1
                        current_line.clear()

                        if write_buf_len >= write_buffer_chars:
                            ofi.write("".join(write_buf))
                            write_buf.clear()
                            write_buf_len = 0

            for block_start in range(0, len(scaf), read_chunk_bases):
                block_end = min(block_start + read_chunk_bases, len(scaf))
                block = str(scaf[block_start:block_end])

                for ch in block:
                    nrep = hemi if comp_index == 0 else period

                    remaining = scaf_length - produced
                    if remaining <= 0:
                        break
                    if nrep > remaining:
                        nrep = remaining

                    push_seq(ch * nrep)
                    produced += nrep
                    comp_index += 1

                if produced >= scaf_length:
                    break

            if produced < scaf_length:
                npad = scaf_length - produced
                push_seq("N" * npad)
                produced += npad

            if current_line:
                line = "".join(current_line)
                write_buf.append(line)
                write_buf.append("\n")

            if write_buf:
                ofi.write("".join(write_buf))

            if produced != scaf_length:
                raise RuntimeError(
                    f"Decompressed scaffold {header_new!r} has length {produced}, "
                    f"expected {scaf_length}."
                )