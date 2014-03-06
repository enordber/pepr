package edu.vt.vbi.ci.pepr.tree;

import java.util.BitSet;
import java.util.HashMap;
import java.util.HashSet;

import edu.vt.vbi.ci.pepr.alignment.SequenceAlignment;
import edu.vt.vbi.ci.util.ExtendedBitSet;

public interface TreeI {

	public abstract int getNodeCount();

	public abstract int getBranchCount();

	/**
	 * Returns an array containing all TreeNodes in this Tree.
	 * @return
	 */
	public abstract TreeNodeI[] getNodes();

	public abstract TreeNodeI[] getLeafNodes();

	public abstract TreeNode[] getInternalNodes();

	/**
	 * Chooses the branch to cut yielding the largest representative
	 * clade in the tree. Returns the set of leaf nodes in this
	 * subtree. If this cut produces two representative subtrees,
	 * both sets of nodes are returned.
	 * 
	 * @param idToGenomeMap
	 * @return
	 */
	public abstract TreeNodeI[][] getLargestRepresentativeSubset(
			HashMap idToGenomeMap);

	public abstract BitSet[] getNodeDefinitions(HashMap idToGenomeMap);

	/**
	 * Returns a list of HashSets, where each HashSet
	 * contains the genome names for a single clade.
	 * 
	 * TODO need to restore initial rooting after this is done.
	 *       As it is written, this method will alter the 
	 *       rooting of the tree.
	 *       
	 * @param idToGenomeMap
	 * @return
	 */
	public abstract HashSet[] getCladeTaxa(HashMap idToGenomeMap);

	/**
	 * Returns a HashSet of HashSet pairs. Each pair represents
	 * one bipatition in the tree. The HashSets contain the
	 * names of the genomes represented in the clade. The pairs
	 * are in the form of HashSet[], with the smaller 
	 */
	public abstract HashSet getBipartitionTaxa(HashMap idToGenomeMap);

	/**
	 * Indicates if the tree is currently rooted.
	 * 
	 * @return
	 */
	public abstract boolean isRooted();

	/**
	 * Roots the tree at the given location along
	 * the given branch. rootPoint must be between
	 * 0.0 and 1.0 inclusive.
	 * 
	 * @param b
	 * @param rootPoint
	 */
	public abstract void root(TreeBranch b, double rootPoint);

	/**
	 * Roots the tree at the center of the
	 * given branch.
	 * 
	 * @param b
	 * @return the new root node
	 */
	public abstract void root(TreeBranch b);

	/**
	 * Returns the current 'root branch'. If the 
	 * tree is not rooted, this returns null.
	 * 
	 * @return
	 */
	public abstract TreeBranch getRootBranch();

	public abstract double getTotalDistance();

	/**
	 * Collapses internal nodes until the number of leaves
	 * is no greater than maxLeaves.
	 * 
	 * @param maxLeaves
	 */
	public abstract void collapseToMaximumLeaves(int maxTips);

	public abstract int getLeafCount();

	/**
	 * Returns the distance to the deepest tip/leaf.
	 * TODO his is currently checking distance to all leaves,
	 * but it need to be modified to only go until it
	 * reaches either a leaf or a collapsed node.
	 * 
	 * @return
	 */
	public abstract double getMaxTipDistance();

	/**
	 * Returns TreeNodes for all tips, in an 
	 * order that would produce a 'laderized up'
	 * tree. That is, nodes with more leaves are
	 * listed before nodes with fewer leaves.
	 * @return
	 */
	public abstract TreeNodeI[] getTipsLaderizeUp();

	/**
	 * Returns TreeNodes for all tips, in an 
	 * order that would produce a 'laderized down'
	 * tree. That is, nodes with fewer leaves are
	 * listed before nodes with more leaves.
	 * @return
	 */
	public abstract TreeNodeI[] getTipsLaderizeDown();

	/**
	 * Returns TreeNodes for all tips in an order that
	 * would produce a balanced tree. That is, nodes
	 * alternate between more leaves up and more leaves
	 * down.
	 * 
	 * @param startUp indicates how to start. If true, the 
	 *                first node will place the largest subnode
	 *                up.
	 * @return
	 */
	public abstract TreeNodeI[] getTipsBalance(boolean startUp);

	/**
	 * Calculates and return the Robinson-Foulds distance
	 * between this Tree and the given Tree. The 
	 * Robinson-Foulds distance is the size of the union
	 * of the clades in the two trees minus the size of
	 * the intersection of the clades in the two trees.
	 * 
	 * @param otherTree
	 * @return
	 */
	public abstract int getRobinsonFouldsDistance(TreeI otherTree);

	public abstract String getTreeString(boolean includeDistance,
			boolean includeSupport);

	/**
	 * This returns a topology-only String of the tree,
	 * with a trifurcation at the highest level. This is 
	 * being done specifically top support the use of 
	 * MrBayes, which requires such a tree for it's 
	 * 'usertree', if a usertree is provided.
	 * @return
	 */
	public abstract String getTreeStringTrifurcating();

	public abstract HashSet getCladeSet();

	public abstract Bipartition[] getBipartitions(HashMap idToTaxonMap,
			String[] distinctTaxa);
	
	/**
	 * Sets the alignment used to construct this Tree.
	 * @param alignment
	 */
	public abstract void setAlignment(SequenceAlignment alignment);

	/**
	 * Returns the SequenceAlignment used to construct this Tree, if
	 * the SequenceAlignment is available. 
	 * 
	 * @return
	 */
	public abstract SequenceAlignment getAlignment();

	public abstract TreeBranch[] getBranches();

	public abstract double[][] getNodeDistanceMatrix(TreeNodeI[] nodes);
}