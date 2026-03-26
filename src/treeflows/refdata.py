import pandas as pd
import random

# from treeflows.aegypti_utilities import FileConfig


# myconfig = FileConfig()

# aeg_base = pd.read_csv(myconfig.refdir / "aegfiles_filtered.csv").sort_values("id")
# asiapac = aeg_base[~aeg_base["continent"].isin(["North America", "South America"]) | aeg_base["pop_short"].isin(["ASNC", "STLC"])]

# should have the dataset for the other samples now... 


def _ref2node(indexlist):
    """Convert row indices to diploid node indices.

    Each sample index `i` maps to the pair `(2*i, 2*i+1)`.

    Args:
        indexlist: Iterable of integer sample indices.

    Returns:
        Flat list of node indices (length `2 * len(indexlist)`).
    """
    return [k for j in [[i*2, i*2+1] for i in indexlist] for k in j]
# assume we have a refdata object... 

class RefData(pd.DataFrame):

    @property
    def _constructor(self):
        return type(self)
    
    def __finalize__(self, other, method=None, **kwargs):
        # copy over custom attrs here if you have any
        return self


    def check_index_intact(self):
        """Raise if the DataFrame index is not 0..n-1 in order."""
        if not [f for f in self.index]==[f for f in range(self.shape[0])]:
            raise Exception("Initial dataframe is not indexed in order!")
    
    def _getsubframe(self, column, valuelist=None):
        """Base method to extract subframe"""
        if valuelist is None:
            return self
        else:
            if type(valuelist) == str:
                valuelist = [valuelist]
            return self[self[column].isin(valuelist)]
        
    def _getsub(self, column, valuelist=None):
        """Base method to just extract column"""
        return self._getsubframe(column, valuelist)[column]
        
    def get_values(self, column, valuelist=None):
        """Return values in `column` optionally filtered by membership in `valuelist`."""
        return [val for val in self._getsub(column, valuelist)]
    
    def get_values_col(self, refcol, valuelist=None, valcol=None):
        """Return values from `valcol` after filtering rows by `refcol`/`valuelist`."""
        temp = self._getsubframe(refcol, valuelist)
        return [val for val in temp[valcol]]
    
    def get_indices(self, column, valuelist=None):
        """Return row indices for rows where `column` is in `valuelist`."""
        return [i for i in self._getsub(column, valuelist).index] # doesn't preserve order
    
    def get_group_indices(self, column, valuelist=None):
        """Return a list of index lists grouped by unique values in `column`."""
        vals = self[column].unique() if valuelist is None else [*dict.fromkeys([valuelist])] if type(valuelist)==str else [*dict.fromkeys(valuelist)] # preserve order
        return [[i for i in self._getsub(column, [val]).index] for val in vals]
    
    def get_nodes(self, column, valuelist=None):
        """Return diploid node indices for rows matching `column`/`valuelist`."""
        return _ref2node(self.get_indices(column, valuelist))
    
    def get_nodes_unique(self, column, valuelist=None):
        """Return one randomly-chosen node per diploid sample for matched rows."""
        nodes = self.get_nodes(column, valuelist)
        return [nodes[i + random.randint(0, 1)] for i in range(len(nodes) // 2)]
    
    def get_group_nodes(self, column, valuelist=None):
        """Return diploid node index lists grouped by unique values in `column`."""
        gindex = self.get_group_indices(column, valuelist)
        return [_ref2node(group) for group in gindex]
    
    

class AegData(RefData):
    # def __init__(self):
    #     super().__init__()
    #     # conditions here... 

    def get_pop_short_idv(self, idv):
        """Return population short code for individual ID `idv`."""
        return self.get_values_col("id", idv, "pop_short")[0]

    def get_idvs_pop(self, poplist=None):
        """Return individual IDs for populations in `poplist` (column `pop`)."""
        return self.get_values_col("pop", poplist, "id")
    def get_idx_pop(self, poplist=None):
        """Return row indices for populations in `poplist` (column `pop`)."""
        return self.get_indices("pop", poplist)
    def get_idx_pop_short(self, poplist=None):
        """Return row indices for populations in `poplist` (column `pop_short`)."""
        return self.get_indices("pop_short", poplist)
    def get_nodes_pop(self, poplist=None):
        """Return diploid node indices for populations in `poplist` (column `pop`)."""
        return self.get_nodes("pop", poplist)
    def get_idvs_country(self, poplist=None):
        """Return individual IDs for countries in `poplist` (column `country`)."""
        return self.get_values_col("country", poplist, "id")
    def get_idx_country(self, poplist=None):
        """Return row indices for countries in `poplist` (column `country`)."""
        return self.get_indices("country", poplist)
    def get_nodes_country(self, poplist=None):
        """Return diploid node indices for countries in `poplist` (column `country`)."""
        return self.get_nodes("country", poplist)
    def get_idvs_continent(self, poplist=None):
        """Return individual IDs for continents in `poplist` (column `continent`)."""
        return self.get_values_col("continent", poplist, "id")
    def get_idx_continent(self, poplist=None):
        """Return row indices for continents in `poplist` (column `continent`)."""
        return self.get_indices("continent", poplist)
    def get_nodes_continent(self, poplist=None):
        """Return diploid node indices for continents in `poplist` (column `continent`)."""
        return self.get_nodes("continent", poplist)
    
    def get_group_idx_pop(self, poplist=None):
        """Return row indices grouped by population (column `pop`)."""
        return self.get_group_indices("pop", poplist)
    def get_group_nodes_pop(self, poplist=None):
        """Return diploid node indices grouped by population (column `pop`)."""
        return self.get_group_nodes("pop", poplist)
    def get_group_idx_country(self, poplist=None):
        """Return row indices grouped by country (column `country`)."""
        return self.get_group_indices("country", poplist)
    def get_group_nodes_country(self, poplist=None):
        """Return diploid node indices grouped by country (column `country`)."""
        return self.get_group_nodes("country", poplist)
    def get_group_idx_continent(self, poplist=None):
        """Return row indices grouped by continent (column `continent`)."""
        return self.get_group_indices("continent", poplist)
    def get_group_nodes_continent(self, poplist=None):
        """Return diploid node indices grouped by continent (column `continent`)."""
        return self.get_group_nodes("continent", poplist)
    
    # spectial utilities
    def get_idx_idv(self, idvlist=None):
        """Return row indices for individuals in `idvlist` (column `id`)."""
        return self.get_indices("id", idvlist)
    def get_nodes_idv(self, idvlist=None):
        """Return diploid node indices for individuals in `idvlist` (column `id`)."""
        return self.get_nodes("id", idvlist)
    # bigger deal...
    def get_group_idx_idv(self, idvlists):
        """Return row indices for multiple individual lists."""
        return [self.get_idx_idv(idvlist) for idvlist in idvlists]
    def get_group_nodes_idv(self, idvlists):
        """Return diploid node indices for multiple individual lists."""
        return [self.get_nodes_idv(idvlist) for idvlist in idvlists]
    
    def get_group_pop_idx(self, idxlists):
        """Map index lists to population labels (column `pop`) preserving group structure."""
        temp = self["pop"]
        return [[t for t in temp.iloc[idxlist]] for idxlist in idxlists]
    
    def get_nodes_idx(self, idxlist):
        """Convert a list of row indices into diploid node indices."""
        return _ref2node(idxlist)
    def get_nodes_idx_unique(self, idxlist):
        """Pick one node per diploid sample for the provided row indices."""
        nodes = _ref2node(idxlist)
        # we want to randomly sample from either first or second node, as there could be genotyping biases
        return [nodes[i + random.randint(0, 1)] for i in range(len(nodes) // 2)]
    def get_group_nodes_idx(self, idxlists):
        """Convert a list of index lists into a list of diploid node index lists."""
        return [_ref2node(idxlist) for idxlist in idxlists]

def load_refdata(ref_csv):
    """Load a reference CSV as a `RefData` dataframe subclass."""
    df = pd.read_csv(ref_csv)
    return RefData(df)

def load_aegdata(ref_csv):
    """Load a reference CSV as an `AegData` dataframe subclass."""
    df = pd.read_csv(ref_csv)
    return AegData(df)
