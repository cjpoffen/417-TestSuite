
from mulrfmodel import *
from compbio import phylo
from rasmus import treelib
import unittest

#mul = MulRFModel(extra = None)
#gene2species = phylo.read_gene2species("../../../examples/config/fungi.smap")
#stree = treelib.read_tree('../../../examples/config/fungi.stree')
#treelib.draw_tree(stree,  minlen=5, maxlen=5)  

#print gene2species("smik_13")

class TestMulRFCost(unittest.TestCase):

    def testWithBinaryTres(self):
        mul = MulRFModel(extra = None)
        gene2species = phylo.read_gene2species("../../../examples/test/24Hits.smap")
        stree = treelib.read_tree('../../../examples/test/24Hits.stree')
        gtree = treelib.read_tree('../../../examples/test/24Hits.gtree')
        
        mul.stree = stree
        mul.gene2species = gene2species
        
        self.assertEqual(mul.compute_cost(gtree), 4)

    def test_Null_Trees(self):
        mul = MulRFModel(extra = None)
        stree = treelib.read_tree('../../../examples/test/EmptyTree.stree')
        gtree = treelib.read_tree('../../../examples/test/EmptyTree.stree')
        gene2species = phylo.read_gene2species("../../../examples/test/24Hits.smap")
        mul.stree = stree
        mul.gene2species = gene2species
        
        with self.assertRaises(AttributeError):
            mul.compute_cost(gtree)
            
    def test_null_gtree(self):
        mul = MulRFModel(extra = None)
        gene2species = phylo.read_gene2species("../../../examples/test/24Hits.smap")
        stree = treelib.read_tree('../../../examples/test/24Hits.stree')
        gtree = treelib.read_tree('../../../examples/test/EmptyTree.stree')
        
        mul.stree = stree
        mul.gene2species = gene2species
        
        with self.assertRaises(AttributeError):
            mul.compute_cost(gtree)
        
   
if __name__ == '__main__':
    unittest.main()