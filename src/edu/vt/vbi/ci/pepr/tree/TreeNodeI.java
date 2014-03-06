package edu.vt.vbi.ci.pepr.tree;

import java.util.HashMap;

public interface TreeNodeI {

	/**
	 * Returns an array of all direct child nodes
	 * of this node from all branches except
	 * fromBranch.
	 * 
	 * @return
	 */
	public abstract TreeNodeI[] getChildNodes(TreeBranch fromBranch);

	/**
	 * Returns an array of all direct child nodes
	 * of this node.
	 * 
	 * @return
	 */
	public abstract TreeNodeI[] getChildNodes();

	/**
	 * Returns an array of all direct child branches
	 * of this node, sorted according to the
	 * number of leaves.
	 * 
	 * @param ascending if true, nodes are sorted in ascending
	 *                  order of number of leaves. If false,
	 *                  nodes are sorted in descending order
	 *                  of number of leaves.
	 * @return
	 */
	public abstract TreeBranch[] getChildBranches(TreeBranch fromBranch,
			boolean ascending);

	/**
	 * Returns an array of all direct child branches
	 * of this node, sorted according to the
	 * number of leaves.
	 * 
	 * @param ascending if true, nodes are sorted in ascending
	 *                  order of number of leaves. If false,
	 *                  nodes are sorted in descending order
	 *                  of number of leaves.
	 * @return
	 */
	public abstract TreeBranch[] getChildBranches(boolean ascending);

	/**
	 * Returns the newick format snippet for this node.
	 * @return
	 */
	public abstract String getOriginalNodeString();

	/**
	 * Return the label for this node. This differs from
	 * toString() in that it does not include branch
	 * length or support information.
	 * @return
	 */
	public abstract String getNodeLabel(TreeBranch fromBranch);

	public abstract boolean isLeaf();

	/**
	 * Returns an array of all branches for this node,
	 * except the given branch.
	 * 
	 * @param aBranch
	 * @return
	 */
	public abstract TreeBranch[] getOtherBranches(TreeBranch aBranch);

	public abstract TreeNodeI[] getDescendantLeaves(TreeBranch fromBranch);

	/**
	 * Returns the total distance between this node and
	 * the current root in the tree.
	 * 
	 * @return
	 */
	public abstract double getDistanceFromRoot();

	public abstract double getDistanceToRoot(TreeBranch fromBranch);

	public abstract TreeBranch getParentBranch();

	/**
	 * Sets the parent branch for this node, resulting in the
	 * orientation of the remainder of the tree beyond this node.
	 * @return
	 */
	public abstract void setParentBranch(TreeBranch parentBranch);

	public abstract TreeNodeI getParentNode();

	/**
	 * Sets the parent node for this node, resulting in the
	 * orientation of the remainder if the tree beyond this node.
	 * @param parentNode
	 */
	public abstract void setParentNode(TreeNodeI parentNode);

	/**
	 * Set the collapsed value for this Node.
	 * @param collapsed
	 */
	public abstract void setCollapsed(boolean collapsed);
	
	/**
	 * Returns the current collapsed state for this Node.
	 * @return
	 */
	public abstract boolean isCollapsed();
}