
from mulrfmodel import *
from deepcoalescencemodel import *
from compbio import phylo
from rasmus import treelib
import unittest

#mul = MulRFModel(extra = None)
#gene2species = phylo.read_gene2species("../../../examples/config/fungi.smap")
#stree = treelib.read_tree('../../../examples/config/fungi.stree')
#treelib.draw_tree(stree,  minlen=5, maxlen=5)  

#print gene2species("smik_13")

class TestMulRFCost(unittest.TestCase):

    def test_with_binary_trees(self):
        mul = MulRFModel(extra = None)
        gene2species = phylo.read_gene2species("../../../examples/test/24Hits.smap")
        stree = treelib.read_tree('../../../examples/test/24Hits.stree')
        gtree = treelib.read_tree('../../../examples/test/24Hits.gtree')
        
        mul.stree = stree
        mul.gene2species = gene2species
        
        self.assertEqual(mul.compute_cost(gtree), 4)

    def test_null_trees(self):
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
        
    def test_non_binary_gene_tree(self):
        mul = MulRFModel(extra = None)
        gene2species = phylo.read_gene2species("../../../examples/test/nonBinaryGene.smap")
        stree = treelib.read_tree('../../../examples/test/nonBinaryGene.stree')
        gtree = treelib.read_tree('../../../examples/test/nonBinaryGene.gtree')
        
        mul.stree = stree
        mul.gene2species = gene2species
        self.assertEqual(mul.compute_cost(gtree), 7)
   
    def test_nonbinary_trees(self):
        mul = MulRFModel(extra = None)
        gene2species = phylo.read_gene2species("../../../examples/test/nonBinaryAll.smap")
        stree = treelib.read_tree('../../../examples/test/nonBinaryAll.stree')
        gtree = treelib.read_tree('../../../examples/test/nonBinaryAll.gtree')
        
        mul.stree = stree
        mul.gene2species = gene2species
        self.assertEqual(mul.compute_cost(gtree), 6)
      
    def test_non_binary_gtree_error(self):
        mul = MulRFModel(extra = None)
        gene2species = phylo.read_gene2species("../../../examples/test/nonBinaryAll.smap")
        stree = treelib.read_tree('../../../examples/test/nonBinaryAll.stree')
        gtree = treelib.read_tree('../../../examples/test/nonBinaryAll.gtree')
        
        
        with self.assertRaises(Exception):
            mul.optimize_model(gtree, stree, gene2species)
            
    def test_non_binary_stree_error(self):
        mul = MulRFModel(extra = None)
        gene2species = phylo.read_gene2species("../../../examples/test/nonBinaryAll.smap")
        stree = treelib.read_tree('../../../examples/test/nonBinaryAll.stree')
        gtree = treelib.read_tree('../../../examples/test/24Hits.gtree')
        
        
        with self.assertRaises(Exception):
            mul.optimize_model(gtree, stree, gene2species)
        
    def test_smap_error(self):
        mul = MulRFModel(extra = None)
        gene2species = phylo.read_gene2species("../../../examples/test/nonBinaryAll.smap")
        stree = treelib.read_tree('../../../examples/test/24Hits.stree')
        gtree = treelib.read_tree('../../../examples/test/24Hits.gtree')
        
        
        with self.assertRaises(Exception):
            mul.optimize_model(gtree, stree, None)

    def test_single_node(self):
        
        sNode = treelib.TreeNode()
        gNode = treelib.TreeNode()
        
        mul = MulRFModel(extra = None)
        gene2species = phylo.read_gene2species("../../../examples/test/nonBinaryAll.smap")
        
        mul.stree = sNode
        mul.gene2species = gene2species
        
        with self.assertRaises(AttributeError):
            mul.compute_cost(gNode)
        
    def test_deep(self):
        deep = DeepCoalescenceModel(extra = None)
        gene2species = phylo.read_gene2species("../../../examples/test/24Hits.smap")
        stree = treelib.read_tree('../../../examples/test/test1.stree')
        gtree = treelib.read_tree('../../../examples/test/test1.gtree')
        
        deep.stree = stree
        deep.gene2species = gene2species
        
        self.assertEqual(deep.compute_cost(gtree), 2)
    
        
if __name__ == '__main__':
    unittest.main()