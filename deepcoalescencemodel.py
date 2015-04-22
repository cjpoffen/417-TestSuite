#
# Python module for robinson foulds cost
#

# treefix libraries
from treefix.models import CostModel

# python libraries
import optparse

# rasmus libraries
from rasmus import treelib

# compbio libraries
from compbio import phylo

#debug
import StringIO
#=============================================================================

class DeepCoalescenceModel(CostModel):
    """Computes Robinson-Foulds or MulRF costs"""

    def __init__(self, extra):
        """Initializes the model"""
        CostModel.__init__(self, extra)

        self.VERSION = "1.0.1"
        self.mincost = 0
        self.count = 0
        self.log = open('matched.txt', 'w')
        
    def optimize_model(self, gtree, stree, gene2species):
        """Optimizes the model"""
        CostModel.optimize_model(self, gtree, stree, gene2species)
        
        # ensure gtree and stree are both rooted and binary
        if not (treelib.is_rooted(gtree) and treelib.is_binary(gtree)):
            raise Exception("gene tree must be rooted and binary")
        if not (treelib.is_rooted(stree) and treelib.is_binary(stree)):
            raise Exception("species tree must be rooted and binary")
        try:
            junk = phylo.reconcile(gtree, stree, gene2species)
        except:
            raise Exception("problem mapping gene tree to species tree")
    
    def compute_cost(self, gtree):
        """Returns the duplication-loss cost"""
        #recon = phylo.reconcile(gtree, self.stree, self.gene2species)            
        
        recon = {}
        streeCopy = self.stree.copy()
        
        def recon_dup(stree, gtree, recon):
            #Recon plus find duplications in gene tree, for each gene tree leaf we map it to a stree leaf, we also keep track of duplications in the gene tree 
            leaves = gtree.leaves()
            for node in leaves:
                sNode = stree.nodes[self.gene2species(node.name)]
                recon[node] = sNode
                
        def printLca(tree, lca):

            for node in lca:
                if  len(node.leaves()) > 1 and (node is not tree.root):
                    print >> self.log, node
                    print >> self.log, lca[node]
        
        #LCA method 
        def lca(node, lca_dict):
            """Creates a dictionary of (node, lca) pairs from given tree"""
            if node.is_leaf():
                lca_dict[node] = []
                lca_dict[node].append(node)
            
            else:
                lca_dict[node] = []
                append = lca_dict[node].append
                for child in node.children:
                    lca(child, lca_dict)
                    for x in lca_dict[child]:
                        append(x)

        stree_lca_dict = {}
        gtree_lca_dict = {}
        
        geneToSpeciesMap = {}
        
        recon_dup(streeCopy, gtree, recon)
        lca(gtree.root, gtree_lca_dict)
        lca(streeCopy.root, stree_lca_dict)
                       
        #create a mapping from gene nodes to species nodes
        for gNodeLca in gtree_lca_dict:
            if not gNodeLca.is_leaf():
                gNode = gtree_lca_dict[gNodeLca][0] #Get the first node in the lca 
                sNode = recon[gNode]
                
                gene2speciesLca = []
                for node in gtree_lca_dict[gNodeLca]:
                    gene2speciesLca.append(recon[node])
                
                found = False
                    
                while not found:
                    sNode = sNode.parent
                    for node in gene2speciesLca:
                        if node not in stree_lca_dict[sNode]:
                            break
                    else:
                        found = True
                        geneToSpeciesMap[gNodeLca] = sNode     
        
        #count distance
        count = 0
        for node in geneToSpeciesMap:
            for child in node.children:
                
                if child.is_leaf():
                    snode = recon[child]
                else:
                    snode = geneToSpeciesMap[child]
                    
                while snode != geneToSpeciesMap[node]:
                    count +=1
                    snode = snode.parent
                    
        #assuming the number of edges equals the number of nodes - 1 
        def countChildern(node):
            count = len(node.children)
            for child in node.children:
                count += countChildern(child)
                
            return count
                
                
        count = count - countChildern(self.stree.root)
        
        return count
     

#cherry yum diddly dip