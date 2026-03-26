import numpy as np
from pathlib import Path
import gzip
from random import Random
from cyvcf2 import VCF
import functools

_INFO_VERSION = '0.1'
_INFO_SOURCE = 'vcf-wrangler'


class VarParsed:
    '''A companion class to the Variant class from cyvcf2. 
    Call with Variant class to access a set of properties and functions calculated on it. 

    Attributes: 
        var             (Variant)   The variant class it operates on
        gts             (uint16)    Flattened numpy array of (non-missing) genotypes (no longer mapped to individuals)
        v_a             (int)       Indices of (unique) nonmissing alleles (int). Usable with `vcf_alleles`
        v_c             (int)       Counts of nonmissing alleles (matches v_a). (int)
        allelicity      (int)       Number of nonmissing allele types at site (int)
        is_nonallelic   (bool)      Whether site contains any allelic information (allelicity == 0)
        is_monoallelic  (bool)      Whether site is monoallelic/invariant (allelicity == 1)
        is_biallelic    (bool)      Whether site is biallelic (allelicity == 2)
        is_multiallelic (bool)      Whether site is multiallelic (allelicity > 2)
        is_nonsnp       (bool)      Whether site is indel, etc. (i.e. not SNP or invariant)
        depths          (int32)     Flattened numpy array of site depths per individual (or 0 if no sites)
        has_data         (bool)      Whether site has genotyped data (i.e. at least 1 individual with depth > 3)
        vcf_alleles     (str)       List of VCF allele identities from combined REF & ALT fields (e.g. 'A' 'C' 'G' 'T' '*'). Corresponds to VCF but not nec. to data. 
                                        Often the basis of slices from e.g. v_a
        has_ref          (bool)      Whether REF allele is present in actual samples
        geno_alleles    (str)       Like `vcf_alleles`, but only contains alleles that are acturally present in the data. Otherwise preserves order. 
        geno_maj        (str)       Allele that is highest frequency in the data. None if data missing or nonallelic
        geno_min        (str)       Allele that is at lowest (non-zero) frequency in the data. Only non-None for biallelic sites at present. 
        geno_2nd        (str)       Allele that is at 2nd highest frequency in the data. Only not None for bi-multi-allelic sites. 
        an              (int)       Total count of alleles present in the data [currently not working adequately with cyvcf2 Writer]
        dp              (int)       Total read depth for individuals included in the data [currently not working adequately with cyvcr2 Writer]
        has_del          (bool)      Whether site contains a deleted or ambiguous marker (* or N) in `geno_alleles` (i.e alleles present in the data)
    '''
    def __init__(self, var, sortkey=False):
        """Wrap a cyvcf2 Variant with convenience accessors and derived stats.

        Args:
            var: cyvcf2 `Variant` instance.
            sortkey: Optional sort key / index array for sample ordering.
        """
        self.var = var
        self._sortkey=sortkey # either False or an np. sort array... 
        self._depthmodded=False
        self._allele_depths = None

    ### MAJOR STARTUP FUNTIONS... 
    # @functools.cached_property
    # def _sortkey_2(self):
    #     sk2 = np.empty(self._sortkey.size*2, dtype=self._sortkey.dtype)
    #     sk2[0::2]=self._sortkey*2
    #     sk2[1::2]=self._sortkey*2+1
    #     return sk2
    @functools.cached_property
    def gts(self):
        """Flattened array of non-missing genotype allele codes (dtype uint16)."""
        return np.array([h for g in self.var.genotypes for h in g[:2] if h >= 0], dtype=np.uint16)
    @functools.cached_property
    def missing_frac(self):
        """Fraction of allele calls that are missing (-1) across samples."""
        return np.mean([h < 0 for g in self.var.genotypes for h in g[:2]])
    @functools.cached_property
    def _gts(self):
        """Tuple of (unique_allele_codes, counts) for `gts`."""
        return np.unique_counts(self.gts)
    @property
    def v_a(self):
        """Array of unique allele codes present in the genotype calls."""
        return self._gts[0]
    @property
    def v_c(self):
        """Array of counts corresponding to `v_a`."""
        return self._gts[1]
    @property
    def allelicity(self):
        """Number of allele types present in genotype calls."""
        return len(self.v_a)
    @property
    def is_nonallelic(self):
        """True if no allele calls are present (allelicity == 0)."""
        return self.allelicity == 0
    @property
    def is_monoallelic(self):
        """True if site is invariant in called genotypes (allelicity == 1)."""
        return self.allelicity == 1
    @property
    def is_biallelic(self):
        """True if exactly two alleles are present (allelicity == 2)."""
        return self.allelicity == 2
    @property
    def is_multiallelic(self):
        """True if more than two alleles are present (allelicity > 2)."""
        return self.allelicity > 2
    @property
    def ref_is_sv(self):
        """True if the REF allele length suggests a structural variant-like record."""
        return len(self.REF) > 1
    @property
    def is_ordered(self):
        """True if allele codes are consistent with VCF allele ordering."""
        return all(self.v_a < self.allelicity)
    @property
    def is_nonsnp(self):
        """True if the record is not a simple SNP (MNP/INDEL/SV or long REF)."""
        v = self.var
        return v.is_mnp or v.is_indel or v.is_sv or self.ref_is_sv
    @property
    def _depths(self):
        """Raw DP format field array as provided by cyvcf2 (may be None)."""
        return self.var.format('DP', int)
    @property
    def allele_depths(self):
        """Raw AD format field array as provided by cyvcf2 (or overridden)."""
        if self._depthmodded:
            return self._allele_depths
        return self.var.format('AD', int)
    @functools.cached_property
    def depths(self):
        """Flattened per-sample depth array (int32), defaulting to [0] if missing."""
        return self._depths.flatten().astype(np.int32) if self._depths is not None else np.array([0], dtype=np.int32)
    @property
    def depths_from_alleles(self):
        """Per-sample total depth computed by summing AD across alleles."""
        return np.array(self.allele_depths, dtype=np.int32).sum(axis=1) if self.allele_depths is not None else np.array([0], dtype=np.int32)
    @functools.cache
    def low_depth_frac(self, depth=3):
        """Fraction of samples with depth < `depth` using DP."""
        return np.mean([d < depth for d in self.depths])
    def low_depth_frac_alleles(self, depth=3):
        """Fraction of samples with depth < `depth` using summed AD."""
        return np.mean([d < depth for d in self.depths_from_alleles])
    @property
    def has_data(self): #fix
        """True if any sample has depth > 3."""
        return any(self.depths > 3)
    @functools.cached_property
    def vcf_alleles(self):
        """List of allele strings from REF + ALT fields (VCF ordering)."""
        return [self.var.REF] + self.var.ALT
    @property
    def has_ref(self):
        """True if the REF allele is present among genotype calls."""
        return 0 in self.v_a
    @functools.cached_property
    def geno_alleles(self):
        """List of allele strings actually present in genotype calls (ordered by `v_a`)."""
        return [self.vcf_alleles[u] for u in self.v_a]
    @property
    def recode(self):
        """True if genotype-present alleles differ from the VCF allele list length."""
        return len(self.geno_alleles) != len(self.vcf_alleles)
    @property
    def geno_maj(self):
        """Most frequent allele among genotype calls (or None if nonallelic)."""
        return None if self.is_nonallelic else self.geno_alleles[np.argmax(self.v_c)]
    @property
    def geno_min(self):
        """Least frequent allele among genotype calls (biallelic only)."""
        return self.geno_alleles[np.argmin(self.v_c)] if self.is_biallelic else None
    @property
    def geno_2nd(self):
        """Second-most frequent allele among genotype calls (allelicity > 1)."""
        return self.geno_alleles[np.argsort(self.v_c)[-2]] if self.allelicity > 1 else None
    @property
    def geno_2nd_frac(self):
        """Frequency of the second-most frequent allele (allelicity > 1)."""
        return self.v_c[np.argsort(self.v_c)[-2]] / self.an if self.allelicity > 1 else 0
    @property
    def geno_2nd_count(self):
        """Count of the second-most frequent allele (allelicity > 1)."""
        return self.v_c[np.argsort(self.v_c)[-2]] if self.allelicity > 1 else 0 # note this doesn't work with polarization... 
    @property
    def an(self):
        """Total number of called allele copies (2 per non-missing diploid sample)."""
        return len(self.gts)
    @property
    def dp(self):
        """Total depth across samples (sum of `depths`)."""
        return np.sum(self.depths)
    @property
    def has_del(self):
        """True if a deletion/ambiguous marker ('*' or 'N') is present in `geno_alleles`."""
        return any(g in ['N', '*'] for g in self.geno_alleles)
    @property
    def has_singleton(self):
        """True if any allele is present exactly once among called allele copies."""
        return any([c == 1 for c in self.v_c]) # base this on allele counts
    @property
    def has_doubleton(self):
        """True if any allele is present exactly twice among called allele copies."""
        return any([c == 2 for c in self.v_c]) # base this on allele counts
    @property
    def has_tripleton(self):
        """True if any allele is present exactly three times among called allele copies."""
        return any([c == 3 for c in self.v_c]) # base this on allele counts
    
    def has_lowcount(self, thresh=3):
        """This does not include zero. Lowcount, not nocount (returns if has counts between 1 & thres (inclusive))"""
        return any([c <= thresh and c > 0 for c in self.v_c])
    @property
    def INFO(self):
        """INFO dict-like accessor for the underlying cyvcf2 Variant."""
        return self.var.INFO

    ## going to try and put in a bunch of properties here
    @property
    def genotypes(self):
        """Raw cyvcf2 genotype triples for each sample."""
        return self.var.genotypes
    # @genotypes.setter # example of setter object for reference
    # def genotypes(self, value):
    #     self.var.genotypes = value
    @property
    def POS(self):
        """1-based position from the underlying VCF record."""
        return self.var.POS
    @property
    def CHROM(self):
        """Chromosome/contig name from the underlying VCF record."""
        return self.var.CHROM
    @property
    def REF(self):
        """Reference allele string from the underlying VCF record."""
        return self.var.REF
    @property
    def ALT(self):
        """ALT allele list from the underlying VCF record."""
        return self.var.ALT
    # @property
    # def INFO(self):
    #     return self.var.INFO

    def __str__(self):
        """String representation of the wrapped cyvcf2 Variant record."""
        return str(self.var) # for now. need to think about this longer term... 
    
    def get_genocounts(self):
        '''Returns tuple of (a) integer indexes (to `vcf_alleles`) of existing alleles, and (b) the (integer) counts of the same alleles'''
        return self.v_a, self.v_c
    
    def get_geno_ref_idx(self):
        '''Returns index of REF in `vcf_alleles` if REF present in data, otherwise index to the first actually existing allele (int)'''
        return self.v_a[0]
    
    def get_geno_ref(self):
        '''Returns REF if REF present in data, otherwise the 0th index of array of existing alleles (str) '''
        return self.geno_alleles[0]
    
    def get_geno_alt_idx(self):
        '''Returns index to highest-freq ALT allele in `vcf_alleles` if REF exists; otherwise, highest other lowest-indexed (int)'''
        return np.argmax(self.v_c[1:]) + 1
    
    def get_geno_alt(self):
        '''Returns highest-frequency ALT allele in dataset, if REF allele present; otherwise, highest ALT allele excluding the lowest-indexed (which is treated as REF)'''
        return self.vcf_alleles[self.get_geno_alt_idx()] 
    
    def get_missing_fraction(self):
        '''Returns a float showing percentage of missing data at site across current samples: 0 is none missing 1 is all missing'''
    
    def check_ref_present(self):
        '''Returns True if REF allele is present in the data, False otherwise. Ultimately derives from genotype values rather than depths, probs, etc.'''
        return self.v_a[0] == 0
    
    def passes_filter(self, monoallelic=False, multiallelic=False, nonsnp=False, biallelic=True, badref=False):
        '''Quick filter to get rid of common allele types. Set `True` to allow datatype through, `False` to prevent. 
        Returns True if data passes filter (and can be used), False if it fails'''

        if self.is_nonallelic:
            return False
        if not badref and self.ref_is_sv:
            return False
        if monoallelic and self.is_monoallelic:
            return True
        if not nonsnp and self.is_nonsnp:
            return False
        if multiallelic and self.is_multiallelic:
            return True
        if nonsnp and self.is_nonsnp:
            return True
        if biallelic and self.is_biallelic:
            return True
        return False
    
    def genotypes_replace(self, var):
        '''Replaces genotypes in current variant with those from another, provided they exactly match chromosome & position. Needs to be another VarParsed Object
        Returns True if operation succeeds, False otherwise (need to fix this? )
        Parameter: var (other VarParsed object to replace from); careful: whether to remove any chances of inappropriate allele subsitution. 
            This point matters where alleles have been removed via missingness after file divergence'''
        if not self.CHROM == var.CHROM:
            exit(f'Chromosomes are not equal! {self.CHROM}, {var.CHROM}')
        if not self.POS == var.POS:
            exit(f'Positions are not equal! {self.POS}, {var.POS}')
        if not self.vcf_alleles == var.vcf_alleles:
            gmap = self._allelemap(var.geno_alleles, self.geno_alleles) # remake ??? 
            try: 
                ngen = [[gmap[j[0]]] + [gmap[j[1]]] + [j[2]] for j in var.genotypes]
            except IndexError:
                return False
            # testarray = np.unique_counts(np.array([h for g in ngen for h in g[:2] if g[0] >= 0], dtype=np.uint16))
            # if len(testarray[0]) != len(self.v_a):
            #     return False
            if any([i not in self.v_a for i in np.unique_counts(np.array([h for g in ngen for h in g[:2] if h >= 0], dtype=np.uint16))[0]]):
                return False
            self.var.genotypes = ngen
            return True
        # if any([i not in self.v_a for i in np.unique_counts(np.array([h for g in var.genotypes for h in g[:2] if h >= 0], dtype=np.uint16))[0]]):
        #     return False
        self.var.genotypes = var.genotypes
        return True
    
    def _allelemap(self, va1, va2):
        '''make index from first to second... '''
        v1 = np.array(va1) # the big one
        v2 = np.array(va2) # the small one
        vin = np.isin(v1, v2)
        return np.concatenate([np.where(vin)[0], np.where(np.invert(vin))[0]])
    
    def singleton_to_miss(self):
        """Function that (a) checks for singleton alleles, and (b) converts them to missing"""
        MISSING_VAL = -1
        if not self.has_singleton:
            return False
        genos = self.var.genotypes
        single_idx = self.v_a[np.where(self.v_c == 1)] # gives us npy array of int numbers corresponding to singletons
        self.var.genotypes = [[MISSING_VAL if a[0] in single_idx else a[0]] + [MISSING_VAL if a[1] in single_idx else a[1]] + [a[2]] for a in genos]
        adepths = self.var.format('AD', int)
        mask = [x in single_idx for x in self.geno_alleles] # should work
        adepths = [[d[i] if m else 0] for d in adepths for i, m in enumerate(mask)]
        newvar = VarParsed(self.var)
        newvar._allele_depths = adepths
        newvar._depthmodded = True
        return newvar
    
    def lowcount_to_miss(self, c=2):
        """Function that (a) checks for doubleton/singleton alleles, and (b) converts them to missing"""
        MISSING_VAL = -1
        if not self.has_lowcount(c):
            return False
        genos = self.var.genotypes
        multi_idx = tuple(self.v_a[(self.v_c > 0) & (self.v_c <= c)])
        self.var.genotypes = [[MISSING_VAL if a[0] in multi_idx else a[0]] + [MISSING_VAL if a[1] in multi_idx else a[1]] + [a[2]] for a in genos]
        adepths = self.var.format('AD', int)
        mask = [x in multi_idx for x in self.geno_alleles] # should work
        adepths = [[d[i] if m else 0] for d in adepths for i, m in enumerate(mask)]
        newvar = VarParsed(self.var)
        newvar._allele_depths = adepths
        newvar._depthmodded = True
        return newvar
    

    
    def make_trimmed_fields(self):
        '''produce the first part of the genotype fields'''
        return f"""{self.var.CHROM}\t{self.var.POS}\t{'.' if self.var.ID is None else self.var.ID}\t{self.get_geno_ref()}\t{','.join(self.geno_alleles[1:]) if len(self.geno_alleles) > 1 else '.'}\t
{'.' if self.var.QUAL is None else self.var.QUAL}\t
{'.' if self.var.FILTER is None else self.var.FILTER}\t{';'.join('='.join(map(str, inner)) for inner in self.var.INFO)}\t{':'.join(self.var.FORMAT)}"""


class VCFHeader:
    """
    A function that loads, modifies, and (possibly) writes VCF headers. gz is possible but no .bcf format. 
    """

    def __init__(self, input_path, contigs=list()):
        '''
        Initialize the header object from file path
        '''
        self.source = 'VCFHeader-v0.1'
        input_path = Path(input_path)
        if not '.vcf' in input_path.suffixes and not '.bcf' in input_path.suffixes: exit('File is not a vcf! (does not contain appropriate suffix)') # what about .bcf
        zipped = True if '.gz' in input_path.suffixes else False
        zmode = 'rb' if zipped else 'r'

        self._contigs = contigs
        _use_contigs = True if len(contigs) > 0 else False
        self._n = 0

        # load the relevant file
        with gzip.open(input_path, zmode) as ifi:
            hdict = self._headerdict_make()

            for n, line in enumerate(ifi):
                self._n = n
                ll = line.decode()
                if ll.startswith('##'):
                    hdict = self._headerdict_addline(hdict, ll, n)
                else:
                    valdict = {'id': self.source, 'num': n, 'line': f'##source={self.source}\n'}
                    hdict['source'][self.source] = valdict
                    self.titleline = ll
                    break
            self.hdict = hdict

    def add_info_to_header(self, infodict: dict, throw = True, overwrite = False, gentle=False):
        '''replicates the function in VCF file format'''
        if overwrite and gentle:
            exit('`overwrite` and `gentle` cannot both be set to True!')
        if not throw and not overwrite and not gentle:
            exit('No id-handling flags set (throw, overwrite, gentle)!')
        throw = not (overwrite or gentle)

        infostr = f'##INFO=<ID={infodict['ID']},Number={infodict['Number']},Type={infodict['Type']}'
        infostr += f',Description="{infodict['Description']}"'
        infostr += f',Source="{infodict['Source'] if 'Source' in infodict.keys() else _INFO_SOURCE}"'
        infostr += f',Version="{infodict['Version'] if 'Version' in infodict.keys() else _INFO_VERSION}">\n'
        #print(infostr)
        self._n += 1
        self._headerdict_addline(self.hdict, infostr, self._n, throw=throw, overwrite = overwrite)


    def write_header(self, outfile):
        """Write the reconstructed header (including title line) to a gzipped file."""
        with gzip.open(Path(outfile), 'wb') as ofi:
            self._write_header(ofi, True)

    def _write_header(self, f, title=False):
        """Write header categories to an open binary stream."""
        '''file stream'''
        for cat in ['fileformat', 'filedate', 'source', 'reference', 'INFO', 'FILTER', 'FORMAT', 'ALT', 'assembly', 'contig', 'SAMPLE', 'PEDIGREE', 'misc']:
            self._writeheadcat(f, cat)
        if title:
            self._write_titleline(f)

    def _write_titleline(self, f):
        """Write the `#CHROM ...` title line to an open binary stream."""
        '''file stream'''
        f.write(self.titleline.encode())
            

    def _headerdict_make(self):
            """Create an empty categorized header dictionary."""
            return {
                'fileformat': dict(), 'filedate': dict(), 'source': dict(), 'reference': dict(),
                'INFO': dict(), 'FILTER': dict(), 'FORMAT':dict(), 'ALT': dict(), 'assembly': dict(), 'contig': dict(),
                'SAMPLE': dict(), 'PEDIGREE': dict(), 'misc': dict()
                }
    
    def _headerdict_addline(self, hdict, ll, n, throw = False, overwrite = False):
        """Parse a header line and insert it into the categorized header dict."""
        lll = ll[2:-1] # trimming - this whole code needs a much better long-term solution.
        lparse = lll.split('=<')
        if len(lparse) == 1: # i.e. the split failed & the format is irregular
            lparse = lll.split('=')
            cat = lparse[0]
            id = lparse[1]
            if not cat in hdict.keys():
                cat = 'misc'
            if not id in hdict[cat]:
                valdict = {'id': id, 'num': n, 'line': ll}
                hdict[cat][id] = valdict
                return hdict
            else:
                return hdict

        cat = lparse[0]
        if not cat in hdict.keys():
            cat = 'misc'
        val = lparse[1]
        id = val.split(',')[0].split('=')[1]
        if not id in hdict[cat]: #
            valdict = {'id': id, 'num': n, 'line': ll}
            hdict[cat][id] = valdict
        elif throw:
            exit(f'ID {id} is already in header!')
        elif overwrite:
            n = hdict[cat][id]['num']
            valdict = {'id': id, 'num': n, 'line': ll}
            hdict[cat][id] = valdict
        return hdict
    
    def _writeheadcat(self, f, category):
        """Write a single header category to an open binary stream."""
        for val in self.hdict[category]:
            f.write(self.hdict[category][val]['line'].encode())


class VCFReader:

    def __init__(self, input_path):
        """Open a VCF and expose iteration yielding `VarParsed` objects."""
        self.header = VCFHeader(input_path)
        self.vcf = VCF(input_path)
        self.samples = np.array(self.vcf.samples, dtype=np.str_)
        self._sortkey = False
        self.idv_mask = np.full(len(self.samples), True, dtype=np.bool) # make numpy?
        self.idx = np.array([n for n in range(len(self.idv_mask))], dtype=np.int32)

    @property
    def ordered_samples(self):
        '''Prints based on sample order. Samples remain ordered by original vcf order'''
        if self._sortkey:
            return self.samples[self._sortkey]
        return self.samples

    def __call__(self, region=None):
        """Enabling region access on indexed files"""
        self._iter = iter(self.vcf(region))
        return self

    def __iter__(self):
        """Iterate over all records in the VCF."""
        self._iter = iter(self.vcf)
        return self
    
    def __next__(self):
        """Return the next record as a `VarParsed` wrapper."""
        raw = next(self._iter)
        return VarParsed(raw, self._sortkey)
    
    def mask_update(self):
        """Placeholder for updating sample/site masks (not yet implemented)."""
        pass

    def set_samples(self, samplelist):
        """Restrict the VCF view to a subset of samples (order preserved by cyvcf2)."""
        if any(s not in self.vcf.samples for s in samplelist):
            exit('sample not present in VCF!')
        self.vcf.set_samples([str(s) for s in samplelist])
        self.samples = np.array(self.vcf.samples, dtype=np.str_)

    def subsample(self, n):
        """Randomly subsample `n` samples and set them as the active sample set."""
        if n > len(self.vcf.samples):
            exit(f'Too large a sample set! max: {len(self.vcf.samples)}')
        samplelist = Random().sample([m for m in self.vcf.samples], n)
        self.set_samples(samplelist)

    def set_samples_ordered(self, samplelist):
        """Set samples and also attach a sort key to recover requested ordering."""
        self.set_samples(samplelist)
        self._sortkey = np.argsort([np.where(np.array(samplelist, dtype=np.str_) == i)[0][0] for i in self.samples]) # this injects the sortkey for appropriate functions


    def add_info_to_header(self, infodict: dict, throw = True, overwrite = False, gentle=False):
        """Add an INFO header line to both the cached header and cyvcf2 VCF handle."""
        self.header.add_info_to_header(infodict, throw, overwrite, gentle)
        self.vcf.add_info_to_header(infodict)

class VCFWriter:

    def __init__(self, output_path, template_vcf=None, vcf_header=None, samplelist=None):
        """Create a gzipped VCF writer compatible with `VCFReader`/`VarParsed`."""
        
        # TODO: samplelist here not properly aware of auto-sorting feature of VCF file substructure 
        
        if template_vcf is None and (vcf_header is None or samplelist is None):
            exit('need template files!')
        self.header = vcf_header if vcf_header is not None else template_vcf.header
        self._path = Path(output_path)
        self._fh = None
        self._header_written = False
        self.samples = samplelist if samplelist is not None else template_vcf.samples

    def _ensure_open(self):
        """Open the gzip file if it isn't already open"""
        if self._fh is None:
            self._fh = gzip.open(self._path, "wb")

    def __enter__(self):
        """Context-manager enter: open file handle and return self."""
        self._ensure_open()
        return self
    
    def __exit__(self, ex_type, exc, tb):
        """Context-manager exit: close the writer (do not suppress exceptions)."""
        self.close()
        return False # do not suppress exceptions


    def add_info_to_header(self, infodict: dict, throw = True, overwrite = False, gentle=False):
        """Add an INFO header line before the header has been written."""
        if self._header_written:
            exit('Header has already been written. Cannot add new info-line!')
        self.header.add_info_to_header(infodict, throw, overwrite, gentle)

    def write_header(self, samplelist = None):
        """Write the VCF header and `#CHROM ...` title line."""
        if samplelist is None:
            samplelist = self.samples
        if self._header_written:
            return 0
        self._ensure_open()
        self.header._write_header(self._fh)
        self._fh.write(('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(samplelist) + '\n').encode())
        self._header_written = True

    def write_record(self, var):
        """Write a `VarParsed` record (or compatible object with `.var`) to the file."""
        self._ensure_open()
        if not self._header_written:
            self.write_header()
        self._fh.write(str(var.var).encode())

    def write_alt_record(self, var):
        """Placeholder for writing alternate record representations (not implemented)."""
        self._ensure_open()
        if not self._header_written:
            self.write_header()
        self._fh.write(''.encode())

    def rectest(self):
        """Placeholder test helper (not functional in current form)."""
        for _ in range(10):
            self.write_record()
        self.close()

    def close(self):
        """Close only if open"""
        if self._fh is not None:
            self._fh.close()
            self._fh = None

class VCFCompare:
    '''Class for running loosely through two VCF files (of similar attributes)
    Currently requires files to have same individuals & only processes 'intercept' method
    Returns a tuple of two matching VarParsed objects (vcf1, vcf2) when iterated. 
    Currently vulnerable to multiple lines in vcf with same position. 

    Parameters:
        samplemode  - which samples to include (*intersect*, left, right, exclude, all)
        sitemode    - which sites to included (*intersect*, left, right, exclude, all)
        actionmode  - which actions to take (*compare*, replace)
        priority    - which site to favour in replacement (*None*, left, right)
        exclusive   - whether to *only* include priority sites (*False*, True)
    '''
    def __init__(self, vcf1, vcf2, samplemode='intersect', sitemode='intersect', actionmode='compare', priority=None, exclusive=False):
        """Create an iterator over two VCFs aligned by position.

        Args:
            vcf1: Left VCF path.
            vcf2: Right VCF path.
            samplemode: Sample handling mode (`intersect`, `left`, `right`, `exclude`, `all`).
            sitemode: Site handling mode (`intersect`, `left`, `right`, `exclude`, `all`).
            actionmode: Either `compare` (yield pairs) or `replace` (yield chosen site).
            priority: Replacement priority (`left`/`right`) when `actionmode='replace'`.
            exclusive: If True, drop sites missing on the priority side.
        """
        self.vcf1 = VCFReader(vcf1)
        self.vcf2 = VCFReader(vcf2)
        if not actionmode in ['compare', 'replace']:
            exit('action mode must be one of: compare, replace')
        self.actionmode = actionmode
        if not samplemode in ['intersect', 'left', 'right', 'exclude', 'all']:
            exit('sample mode must be one of: intersect, left, right, exclude') # exclude means only look at sites in 1 and not other... 
        self.samplemode = samplemode
        if not sitemode in ['intersect', 'left', 'right', 'exclude', 'all']:
            exit('site mode must be one of: intersect, left, right, exclude')
        self.sitemode = sitemode
        if actionmode == "replace":
            # try to sort out priorities
            if priority is None:
                if sitemode in ['intersect', 'all']:
                    exit('priority must be set for intersect & all sitemodes if actionmode is replace')
                if sitemode == "exclude" and exclusive:
                    exit('if exclusive is True and sitemode is exclude, priority must be explicitly set')
                elif sitemode == "exclusive":
                    priority = None
                else:
                    priority = sitemode

        if not priority in [None, 'left', 'right']:
            exit('priority must be one of: None, left, right')
        self.priority = priority
        self.exclusive = exclusive
        # set initial samples based on overlap... 
        if samplemode == 'intersect':
            sample_overlap = [s for s in self.vcf1.samples if s in self.vcf2.samples]
            self.set_samples(sample_overlap)
        elif samplemode == 'all': # poor code
            sample_union = list({s for s in self.vcf1.samples}.union({t for t in self.vcf2.samples}))
            self.samples = np.array(sample_union, dtype=np.str_)
        self.idv_mask = np.full(len(self.samples), True, dtype=np.bool) # make numpy?
        self.idx = np.array([n for n in range(len(self.idv_mask))], dtype=np.int32)
        self.idx = self.idx[:5]
        self._get1 = True
        self._get2 = True
        self._stop = 0
        self._stop1 = False
        self._stop2 = False

    def __iter__(self):
        """Initialize iterators over both VCFs and prime the first records."""
        self._iter1 = iter(self.vcf1)
        self._iter2 = iter(self.vcf2)
        self.raw1 = self._safe_next1()
        self.raw2 = self._safe_next2()
        return self
    
    # --- helpers --------------------

    def _safe_next1(self):
        """Return the next record from vcf1, or None when exhausted."""
        try:
            return next(self._iter1)
        except StopIteration:
            return None
        
    def _safe_next2(self):
        """Return the next record from vcf2, or None when exhausted."""
        try:
            return next(self._iter2)
        except StopIteration:
            return None
        
    def _advance1(self):
        """Advance the vcf1 cursor by one record."""
        self.raw1 = self._safe_next1()

    def _advance2(self):
        """Advance the vcf2 cursor by one record."""
        self.raw2 = self._safe_next2()
        
    
    def _next_inner(self): # TODO: fix this so that one loop can run when other has broken
        """Core merge logic yielding the next aligned record pair."""

        if self.sitemode == "intersect":
            # must have both sides available
            while True:
                if self.raw1 is None or self.raw2 is None:
                    raise StopIteration
                if self.raw1.POS < self.raw2.POS:
                    self._advance1()
                    continue
                if self.raw2.POS < self.raw1.POS:
                    self._advance2()
                    continue
                # equal
                out = (self.raw1, self.raw2)
                self._advance1()
                self._advance2()
                return out
        
        elif self.sitemode == "left":
            # stop when LEFT ends
            while True:
                if self.raw1 is None:
                    raise StopIteration
                if self.raw2 is None or self.raw1.POS < self.raw2.POS:
                    out = (self.raw1, None)
                    self._advance1()
                    return out
                if self.raw2.POS < self.raw1.POS:
                    self._advance2()
                    continue
                # equal
                out = (self.raw1, self.raw2)
                self._advance1()
                self._advance2()
                return out
                    
        elif self.sitemode == "right":
            # stop when RIGHT ends
            while True:
                if self.raw2 is None:
                    raise StopIteration
                if self.raw1 is None or self.raw2.POS < self.raw1.POS:
                    out = (None, self.raw2)
                    self._advance2()
                    return out
                if self.raw1.POS < self.raw2.POS:
                    self._advance1()
                    continue
                # equal
                out = (self.raw1, self.raw2)
                self._advance1()
                self._advance2()
                return out
                
        elif self.sitemode == "all":
            # union: stop only when BOTH end
            while True:
                if self.raw1 is None and self.raw2 is None:
                    raise StopIteration
                if self.raw2 is None or (self.raw1 is not None and self.raw1.POS < self.raw2.POS):
                    out = (self.raw1, None)
                    self._advance1()
                    return out
                if self.raw1 is None or self.raw2.POS < self.raw1.POS:
                    out = (None, self.raw2)
                    self._advance2()
                    return out
                # equal
                out = (self.raw1, self.raw2)
                self._advance1()
                self._advance2()
                return out
                
        elif self.sitemode == "exclude":
            # symmetric difference: emit only non-equal; stop when BOTH end
            while True:
                if self.raw2 is None and self.raw2 is None:
                    raise StopIteration
                # if side ended, emit the other side straight through
                if self.raw2 is None:
                    out = (self.raw1, None)
                    self._advance1()
                    return out
                if self.raw1 is None:
                    out = (None, self.raw2)
                    self._advance2()
                    return out
                # both present
                if self.raw1.POS < self.raw2.POS:
                    out = (self.raw1, None)
                    self._advance1()
                    return out
                if self.raw2.POS < self.raw1.POS:
                    out = (None, self.raw2)
                    self._advance2()
                    return out
                # equal: skip both (exclude)
                self._advance1()
                self._advance2()
                # and continue Loop to find next non-equal
        else:
            raise ValueError(f"Unknown sitemode: {self.sitemode}")
        
    def __next__(self):
        """Yield the next item based on `actionmode` (pair or chosen record)."""
        while True:
            try:
                left, right = self._next_inner()
            except StopIteration:
                # Inner exhausted - just propagate the normal StopIteration
                raise

            if self.actionmode == "compare":
                return left, right
            elif self.actionmode == "replace":
                # we need priority here
                if self.priority == "left":
                    if self.exclusive and left is None:
                        continue
                    if left is None:
                        return right
                    return left
                elif self.priority == "right":
                    if self.exclusive and right is None:
                        continue
                    if right is None:
                        return left
                    return right
                elif not self.exclusive and self.sitemode == "exclude":
                    if right is None:
                        return left
                    return right
                else:
                    raise ValueError(f"Unusable priority combo")


    def set_samples(self, samplelist):
        """Restrict both VCFs to the same subset of samples."""
        if any(s not in self.vcf1.samples for s in samplelist):
            exit('sample not present in VCF!')
        self.vcf1.set_samples([str(s) for s in samplelist])
        self.vcf2.set_samples([str(s) for s in samplelist])
        self.samples = np.array(self.vcf1.samples, dtype=np.str_)

    def subsample(self, n):
        '''Subsample individuals in downstream (without replacement).
        only for tests! (uses first vcf)'''
        if n > len(self.vcf1.samples):
            exit(f'Too large a sample set! max: {len(self.vcf1.samples)}')
        samplelist = Random().sample([m for m in self.vcf1.samples], n)
        self.set_samples(samplelist)

    def _add_info_to_headers(self, infodict: dict, throw = True, overwrite = False, gentle=False):
        """Add an INFO header line to both underlying VCFs."""
        # probably shouldn't use... but available if needed
        self.vcf1.add_info_to_header(infodict, throw, overwrite, gentle)
        self.vcf2.add_info_to_header(infodict, throw, overwrite, gentle)
