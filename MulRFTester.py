
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
    def test(self):
        mul = MulRFModel(extra = None)
        gene2species = phylo.read_gene2species("../../../examples/test/24Hits.smap")
        stree = treelib.read_tree('../../../examples/test/24Hits.stree')
        gtree = treelib.read_tree('../../../examples/test/24Hits.gtree')
        
        mul.stree = stree
        mul.gene2species = gene2species
        
        self.assertEqual(mul.compute_cost(gtree), 4)



class TestStringMethods(unittest.TestCase):

  def test_upper(self):
      self.assertEqual('foo'.upper(), 'FOO')

  def test_isupper(self):
      self.assertTrue('FOO'.isupper())
      self.assertFalse('Foo'.isupper())

  def test_split(self):
      s = 'hello world'
      self.assertEqual(s.split(), ['hello', 'world'])
      # check that s.split fails when the separator is not a string
      with self.assertRaises(TypeError):
          s.split(2)

          
class TestStringMethod(unittest.TestCase):

  def test_upper(self):
      self.assertEqual('foo'.upper(), 'FOO')

  def test_isupper(self):
      self.assertTrue('FOO'.isupper())
      self.assertFalse('Foo'.isupper())

  def test_split(self):
      s = 'hello world'
      self.assertEqual(s.split(), ['hello', 'world'])
      # check that s.split fails when the separator is not a string
      with self.assertRaises(TypeError):
          s.split(2)
          
if __name__ == '__main__':
    unittest.main()