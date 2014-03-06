package edu.vt.vbi.ci.pepr.tree;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;

public class TreeNode implements TreeNodeI, Comparable {
	private static final String commaChar = ",";

	/*
	 * used for sorting nodes based on the number of leaves
	 * in ascending order
	 */
	private TreeBranchSizeAscendingComparator treeNodeSizeAscendingComparator 
	= new TreeBranchSizeAscendingComparator();

	/*
	 * used for sorting nodes based on the number of leaves
	 * in descending order
	 */
	private TreeBranchSizeDescendingComparator treeNodeSizeDescendingComparator
	= new TreeBranchSizeDescendingComparator();

	private HashSet branches;
	private TreeBranch parentBranch;
	//private ArrayList childNodes;
	private String originalNodeString;
	private String currentNodeString;
	private boolean nodeStringNeedsRebuild = true;
	private TreeNodeI parentNode;
	//private boolean isLeaf = false;
	private boolean isRepresentative = false;
	private boolean isRoot = false;
	private boolean isCollapsed = false;
	private String collapsedNodeString;
	private double distanceForCollapsed = Double.NaN;
	private double collapsedDistance;
	private String nodeLabel;
	
	/**
	 * Create a TreeNode without an original String.
	 * This is used for creating "pseudo-nodes", such 
	 * as those used for unrooted trees.
	 *
	 */
	public TreeNode() {
		
	}
	
	public TreeNode(String nodeString) {
		this.originalNodeString = nodeString;
	}

	void addBranch(TreeBranch branch) {
		if(branches == null) {
			branches = new HashSet();
		}
	    branches.add(branch);
	}

	void removeBranch(TreeBranch branch) {
		branches.remove(branch);
	}


	/**
	 * Returns all branches, including the parent branch
	 * if one is currently set. If a parent branch is included,
	 * it will be the last entry in the array (largest index).
	 * 
	 * @return
	 */
	TreeBranch[] getAllBranches() {
		TreeBranch[] r;
//		if(parentBranch == null) {
			r = new TreeBranch[branches.size()];
//		}else {
//			r = new TreeBranch[branches.size()+1];
//			r[r.length -1] = parentBranch;
//		}
		branches.toArray(r);
		return r;
	}

	/* (non-Javadoc)
	 * @see edu.vt.vbi.ci.nongraphic.tree.TreeNodeI#getChildNodes(edu.vt.vbi.ci.nongraphic.tree.TreeBranch)
	 */
	public TreeNodeI[] getChildNodes(TreeBranch fromBranch) {
		TreeNodeI[] r = null;//new TreeNode[childNodes.size()];
		TreeBranch[] cb = getOtherBranches(fromBranch);
		r = new TreeNodeI[cb.length];
		for(int i = 0; i < cb.length; i++) {
			r[i] = cb[i].getOtherNode(this);
		}
		return r;
	}
	
	/* (non-Javadoc)
	 * @see edu.vt.vbi.ci.nongraphic.tree.TreeNodeI#getChildNodes()
	 */
	public TreeNodeI[] getChildNodes() {
		TreeNodeI[] r = null;//new TreeNode[childNodes.size()];
		TreeBranch[] cb = getOtherBranches(parentBranch);
		r = new TreeNodeI[cb.length];
		for(int i = 0; i < cb.length; i++) {
			r[i] = cb[i].getOtherNode(this);
		}
		return r;
	}

	
	/* (non-Javadoc)
	 * @see edu.vt.vbi.ci.nongraphic.tree.TreeNodeI#getChildBranches(edu.vt.vbi.ci.nongraphic.tree.TreeBranch, boolean)
	 */
	public TreeBranch[] getChildBranches(TreeBranch fromBranch, boolean ascending) {
		TreeBranch[] r = getOtherBranches(fromBranch);
		if(ascending) {
			Arrays.sort(r, treeNodeSizeAscendingComparator);
		} else {
			Arrays.sort(r, treeNodeSizeDescendingComparator);
		}
		return r;
	}

	/* (non-Javadoc)
	 * @see edu.vt.vbi.ci.nongraphic.tree.TreeNodeI#getChildBranches(boolean)
	 */
	public TreeBranch[] getChildBranches(boolean ascending) {
		TreeBranch[] r = getOtherBranches(parentBranch);
		if(ascending) {
			Arrays.sort(r, treeNodeSizeAscendingComparator);
		} else {
			Arrays.sort(r, treeNodeSizeDescendingComparator);
		}
		return r;
	}

	
	/**
	 * Returns the number of nodes at or 
	 * below this node.
	 * 
	 * @return
	 */
	int getLeafCount(TreeBranch fromBranch) 
	{
		int r = 0;
		if(isLeaf()) {
			r = 1;
		} else {
			TreeBranch[] childBranches = getOtherBranches(fromBranch);
			for(int i = 0; i < childBranches.length; i++) {
				r += childBranches[i].getOtherNode(this).getLeafCount(childBranches[i]);
			}
		}

		return r;
	}

	HashSet getGenomeSet(TreeBranch fromBranch, HashMap idToGenomeMap) {
		HashSet genomeSet = new HashSet();
		if(isLeaf()) {
			Object genomeName = idToGenomeMap.get(originalNodeString);
			if(genomeName == null) {
				genomeName = originalNodeString;
			}
			genomeSet.add(genomeName);
		} else {
			TreeBranch[] branches = getOtherBranches(fromBranch);
			for(int i = 0; i < branches.length; i++) {
				TreeNode node = branches[i].getOtherNode(this);
				genomeSet.addAll(node.getGenomeSet(branches[i], idToGenomeMap));
			}
		}

		return genomeSet;
	}

	/* (non-Javadoc)
	 * @see edu.vt.vbi.ci.nongraphic.tree.TreeNodeI#getOriginalNodeString()
	 */
	public String getOriginalNodeString() {
		return originalNodeString;
	}

	private void buildCurrentNodeString(TreeBranch fromBranch, boolean includeDistance, boolean includeSupport) 
	{
		nodeStringNeedsRebuild = false;
		if(isLeaf()) {
			currentNodeString = originalNodeString;
		} else {
			StringBuffer sb = new StringBuffer();
			TreeBranch[] children = getOtherBranches(fromBranch);
			if(children.length > 0) {
				sb.append("(");
				int last = children.length - 1;
				for(int i = 0; i < last; i++) {
					sb.append(children[i].toString(this, includeDistance, includeSupport));
					sb.append(commaChar);
				}
				sb.append(children[last].toString(this, includeDistance, includeSupport));
				sb.append(")");
			}
			currentNodeString = sb.toString();
		}
	}

	/* (non-Javadoc)
	 * @see edu.vt.vbi.ci.nongraphic.tree.TreeNodeI#getNodeLabel(edu.vt.vbi.ci.nongraphic.tree.TreeBranch)
	 */
	public String getNodeLabel(TreeBranch fromBranch) {
		if(nodeLabel == null) {
			//construct the node label String
			if(isLeaf()) {
				nodeLabel = originalNodeString;
				int colonIndex = nodeLabel.indexOf(':');
				if(colonIndex > -1) {
					nodeLabel = nodeLabel.substring(0, colonIndex);
				}
			} else {
				StringBuffer sb = new StringBuffer();
				sb.append("(");
				TreeBranch[] branches = getOtherBranches(fromBranch);
				int last = branches.length - 1;
				for(int i = 0; i < last; i++) {
					sb.append(branches[i].getOtherNode(this).getNodeLabel(branches[i]));
					sb.append(commaChar);
				}
				sb.append(branches[last].getOtherNode(this).getNodeLabel(branches[last]));
				sb.append(")");
				nodeLabel = sb.toString();
			}

		}
		return nodeLabel;
	}
	
	String getCollapsedNodeString(TreeBranch fromBranch) {
		if(collapsedNodeString == null) {
			determineDefaultCollpasedNodeString(fromBranch);
		}
		return collapsedNodeString;
	}

	void setCollapsedNodeString(String collapsedNodeString) {
		this.collapsedNodeString = collapsedNodeString;
	}

	public String toString(TreeBranch fromBranch, boolean includeDistance, boolean includeSupport) {
		String r = null;
		if(isCollapsed()) {
			r = getCollapsedNodeString(fromBranch);
		} else {
		    buildCurrentNodeString(fromBranch, includeDistance, includeSupport);
			r = currentNodeString;
		}
		return r;
	}

	/* (non-Javadoc)
	 * @see edu.vt.vbi.ci.nongraphic.tree.TreeNodeI#isLeaf()
	 */
	public boolean isLeaf() {
		boolean r = false;
		if(branches == null || branches.size() == 1) {
			r = true;
		}
		return r;
	}

	String getLeafList(TreeBranch fromBranch) {
		String r = null;
		if(isLeaf()) {
			r = originalNodeString;
		} else {
			TreeBranch[] c = getOtherBranches(fromBranch);
			StringBuffer sb = new StringBuffer();
			sb.append(c[0].getLeafList(this));
			for(int i = 1; i < c.length; i++) {
				sb.append(", ");
				sb.append(c[i].getLeafList(this));
			}
			r = sb.toString();
		}
		return r;
	}

	/* (non-Javadoc)
	 * @see edu.vt.vbi.ci.nongraphic.tree.TreeNodeI#getOtherBranches(edu.vt.vbi.ci.nongraphic.tree.TreeBranch)
	 */
	public TreeBranch[] getOtherBranches(TreeBranch aBranch) {
		TreeBranch[] r = null;
		boolean didContainBranch = branches.remove(aBranch);
		r = new TreeBranch[branches.size()];
		branches.toArray(r);
		if(didContainBranch) {
			branches.add(aBranch);
		}
		return r;
	}
	
	/**
	 * Returns all descendant nodes of all branches from
	 * this node except fromBranch.
	 * @param fromBranch
	 * @return
	 */
	Collection getDescendantNodes(TreeBranch fromBranch) {
		ArrayList r = null;
		TreeBranch[] branches = getOtherBranches(fromBranch);
		//TreeNode[] cn = getChildNodes();
		r = new ArrayList();
		for(int i = 0; i < branches.length; i++) {
			TreeNode node = branches[i].getOtherNode(this);
			r.add(node);
			r.addAll(node.getDescendantNodes(branches[i]));
		}
		return r;
	}

	/* (non-Javadoc)
	 * @see edu.vt.vbi.ci.nongraphic.tree.TreeNodeI#getDescendantLeaves(edu.vt.vbi.ci.nongraphic.tree.TreeBranch)
	 */
	public TreeNodeI[] getDescendantLeaves(TreeBranch fromBranch) {
		TreeNodeI[] r = null;
		Collection dn = getDescendantNodes(fromBranch);
		TreeNodeI[] nodes = new TreeNodeI[dn.size()];
		dn.toArray(nodes);
		ArrayList leaves = new ArrayList();
		for(int i = 0; i < nodes.length; i++) {
			if(nodes[i].isLeaf()) {
				leaves.add(nodes[i]);
			}
		}
		r = new TreeNodeI[leaves.size()];
		leaves.toArray(r);
		return r;
	}

	void setAsRoot() {
		if(parentBranch != null) {
			addBranch(parentBranch);
			parentBranch = null;
		}
		isRoot = true;
	}

	boolean isRoot() {
		return isRoot;
	}

	double getTotalDistance(TreeBranch fromBranch) {
		double r = 0;
		TreeBranch[] childBranches = getOtherBranches(fromBranch);
		for(int i = 0; i < childBranches.length; i++) {
			r += childBranches[i].getTotalDistance(this);
		}
		return r;
	}

	/**
	 * Returns the current distance from this node
	 * to all tips that are descendants of this node.
	 * This takes into account the current collapse
	 * state of this node and all subnodes, and therefore
	 * is dependant on these states for this node and
	 * all descendant nodes.
	 * 
	 * @return
	 */
	double getCurrentDistance() {
		double r = 0;

		return r;
	}


	/**
	 * Returns the distance associated with this node
	 * when it is collapsed. The value returned by this
	 * method is independant of the current collapsed 
	 * state of the node. The value is the length of the
	 * parent pranch leading to this node plus the mean of the
	 * getDistanceForCollapsed() values of any child
	 * nodes. In the case of leaf nodes, the value returned
	 * is just the length of the parent branch leading to this node.
	 * @return
	 */
	double getDistanceForCollapsed(TreeBranch fromBranch) {
		if(Double.isNaN(distanceForCollapsed)) {
			distanceForCollapsed = parentBranch.getBranchLength();
			double childContribution = 0;
			TreeBranch[] branches = getOtherBranches(fromBranch);
			for(int i = 0; i < branches.length; i++) {
				childContribution += 
					branches[i].getOtherNode(this).getDistanceForCollapsed(branches[i]);
			}

			childContribution /= branches.length;
			if(!Double.isNaN(childContribution)) {
				distanceForCollapsed += childContribution;
			}
		}
		return distanceForCollapsed;
	}

	/**
	 * Returns the mean distance between this node and
	 * it's descendant leaves.
	 * 
	 * @return
	 */
	double getMeanChildDistance() {
		double r = 0;
		TreeBranch[] branches =getAllBranches();
		r = Double.MAX_VALUE;
		for(int i = 0; i < branches.length; i++) {
			double meanDistance = getMeanLeafDistance(branches[i]);
			r = Math.min(r, meanDistance);
		}
		
		return r;
	}

	/**
	 * Returns the mean distance from this node
	 * to all leaves. This is basically the 
	 * total distance for the node divided by the
	 * number of child branches the node has.
	 * 
	 * @return
	 */
	double getMeanLeafDistance(TreeBranch fromBranch) {
		double r = 0;
		r = getTotalDistance(fromBranch) / (branches.size() - 1);
		return r;
	}

	/**
	 * Returns the distance between this node and the
	 * given descendant node. If the given node is not
	 * actually a descendant of this node, Double.NaN is
	 * returned.
	 * 
	 * @param descendant
	 * @return
	 */
	double getDistanceToDescendant(TreeNodeI descendant, TreeBranch fromBranch) {
		double r = Double.NaN;
		if(descendant == this) {
			r = 0;
		} else {
			TreeBranch[] childBranches = getChildBranches(fromBranch, true);

			for(int i = 0; i < childBranches.length && Double.isNaN(r); i++) {
				r = childBranches[i].getOtherNode(this).
				       getDistanceToDescendant(descendant, childBranches[i]);
				if(!Double.isNaN(r)) {
					r += childBranches[i].getBranchLength();
				}
			}
		}
		return r;
	}
	
	/* (non-Javadoc)
	 * @see edu.vt.vbi.ci.nongraphic.tree.TreeNodeI#getDistanceFromRoot()
	 */
	public double getDistanceFromRoot() {
		double r = Double.NaN;
		
		r = getDistanceToRoot(null);
		return r;
	}

	
	/* (non-Javadoc)
	 * @see edu.vt.vbi.ci.nongraphic.tree.TreeNodeI#getDistanceToRoot(edu.vt.vbi.ci.nongraphic.tree.TreeBranch)
	 */
	public double getDistanceToRoot(TreeBranch fromBranch) {
		double r = Double.NaN;

		//try each branch, until the root branch is encountered.
		//If a tip is reached, return Double.NaN, indicating
		//that the root is not along that path
		
		//of all of the branches, only one should return a value
		//other than Double.NaN. Once that happens, return it

		TreeBranch[] branches = getOtherBranches(fromBranch);
		
		for(int i = 0; i < branches.length && Double.isNaN(r); i++) {
			r = branches[i].getDistanceFromRoot(this);
		}			
		
		return r;
	}
	/**
	 * Creates the default String to be used for this node
	 * in the collapsed state. This is based on the descendant
	 * leaf node with the shortest distance from this node.
	 */
	private void determineDefaultCollpasedNodeString(TreeBranch fromBranch) {
		TreeNodeI[] leaves = getDescendantLeaves(fromBranch);
		double minDistance = Double.MAX_VALUE;
		double meanDistance = 0;
		int indexWithMin = -1;
		for(int i = 0; i < leaves.length; i++) {
			double distance = getDistanceToDescendant(leaves[i], fromBranch);
			meanDistance += distance;
			if(distance < minDistance) {
				minDistance = distance;
				indexWithMin = i;
			}
		}
		meanDistance /= leaves.length;

		if(indexWithMin > -1) {
			TreeNodeI closestLeaf = leaves[indexWithMin];
			StringBuffer sb = new StringBuffer();
			sb.append(closestLeaf.getNodeLabel(null));
			sb.append("_and_");
			sb.append(leaves.length-1);
			sb.append("_others");
			collapsedNodeString = sb.toString();
		}
	}

	public boolean isCollapsed() {
		return isCollapsed;
	}

	/**
	 * Returns the current number of tips below
	 * this node. A tip is any uncollapsed leaf
	 * or collapsed node. This is the number of 
	 * tips that will be displayed in the current
	 * collapse state.
	 * @return
	 */
	int getTipCount(TreeBranch fromBranch) {
		int r = 0;
		if(isLeaf() || isCollapsed()) {
			r = 1;
		} else {
			TreeBranch[] branches = getOtherBranches(fromBranch); 
			for(int i = 0; i < branches.length; i++) {
				r += branches[i].getOtherNode(this).getTipCount(branches[i]);
			}
		}
		return r;
	}

	public void setCollapsed(boolean isCollapsed) {
		this.isCollapsed = isCollapsed;
	}

	/**
	 * Returns the tip nodes below this node in an
	 * order based on the tree structure, with the
	 * largest child node always coming before the
	 * smaller child node(s). The largest node is 
	 * defined as the one with the greatest number
	 * of leaves, not the greatest number of tips.
	 * This is done to avoid rotating a node simply
	 * as a result of collapsing a subnode.
	 * 
	 * @return
	 */
	ArrayList getTipsLaderizeUp(TreeBranch fromBranch) {
		ArrayList r = new ArrayList();
		TreeBranch[] branches = getOtherBranches(fromBranch);

		if(branches.length > 0) {
			Arrays.sort(branches, treeNodeSizeDescendingComparator);
			TreeNode[] children = new TreeNode[branches.length];
			for(int i = 0; i < branches.length; i++) {
				children[i] = branches[i].getOtherNode(this);
			}			
			for(int i = 0; i < children.length; i++) {
				if(children[i].isTip()) {
					r.add(children[i]);
				} else {
					r.addAll(children[i].getTipsLaderizeUp(branches[i]));
				}
			}
		}

		return r;
	}

	/**
	 * Returns the tip nodes below this node in an
	 * order based on the tree structure, with the
	 * smallest child node always coming before the
	 * larger child node(s). The largest node is 
	 * defined as the one with the greatest number
	 * of leaves, not the greatest number of tips.
	 * This is done to avoid rotating a node simply
	 * as a result of collapsing a subnode.
	 * 
	 * @return
	 */
	ArrayList getTipsLaderizeDown(TreeBranch fromBranch) {
		ArrayList r = new ArrayList();
		TreeBranch[] branches = getOtherBranches(fromBranch);
		if(branches.length > 0) {
			Arrays.sort(branches, treeNodeSizeAscendingComparator);
			for(int i = 0; i < branches.length; i++) {
				TreeNode childNode = branches[i].getOtherNode(this);
				if(childNode.isTip()) {
					r.add(childNode);
				} else {
					r.addAll(childNode.getTipsLaderizeDown(branches[i]));
				}
			}
		}
		return r;
	}

	/**
	 * Returns the tip nodes below this node in an
	 * order based on the tree structure, alternating
	 * between largest node first and smallest node 
	 * first. The result is a tree that should appear
	 * balanced, rather than leaning one direction
	 * The largest node is defined as the one with 
	 * the greatest number of leaves, not the greatest
	 * number of tips. This is done to avoid rotating 
	 * a node simply as a result of collapsing a subnode.
	 * 
	 * @param up true is this node should arrange its
	 *           results with larger nodes before smaller
	 *           nodes. Any child nodes will use the opposite
	 *           of this value.
	 * @return
	 */

	ArrayList getTipsBalance(boolean up, TreeBranch fromBranch) {
		ArrayList r = new ArrayList();
		TreeBranch[] branches = getOtherBranches(fromBranch);
		if(branches.length > 0) {
			if(up) {
				Arrays.sort(branches, treeNodeSizeDescendingComparator);
			} else {
			    Arrays.sort(branches, treeNodeSizeAscendingComparator);
			}
			for(int i = 0; i < branches.length; i++) {
				TreeNode childNode = branches[i].getOtherNode(this);
				if(childNode.isTip()) {
					r.add(childNode);
				} else {
					r.addAll(childNode.getTipsBalance(!up, branches[i]));
				}
			}
		}

		return r;
	}
	
	/**
	 * Returns true if this node is a tip, meaning
	 * it is either collapsed or it is a leaf.
	 * @return
	 */
	private boolean isTip() {
		boolean r = isCollapsed() || isLeaf();
		return r;
	}

	/**
	 * Compares TreeBranches based on the number of leaves
	 * they contain. Used to sort in ascending order.
	 * 
	 * @author enordber
	 *
	 */
	private class TreeBranchSizeAscendingComparator implements Comparator {

		public int compare(Object arg0, Object arg1) {
			TreeBranch branchA = (TreeBranch)arg0;
			TreeBranch branchB = (TreeBranch) arg1;
			TreeNode nodeA = branchA.getOtherNode(TreeNode.this);
			TreeNode nodeB = branchB.getOtherNode(TreeNode.this);
			int leafCountA = nodeA.getLeafCount(branchA);
			int leafCountB = nodeB.getLeafCount(branchB);
			int r = leafCountA - leafCountB;
			return r;
		}

	}

	/**
	 * Compares TreeBranches based on the number of leaves
	 * they contain. Used to sort in descending order.
	 * 
	 * @author enordber
	 *
	 */
	private class TreeBranchSizeDescendingComparator implements Comparator {

		public int compare(Object arg0, Object arg1) {
			TreeBranch branchA = (TreeBranch) arg0;
			TreeBranch branchB = (TreeBranch) arg1;
			TreeNode nodeA = branchA.getOtherNode(TreeNode.this);
			TreeNode nodeB = branchB.getOtherNode(TreeNode.this);
			int leafCountA = nodeA.getLeafCount(branchA);
			int leafCountB = nodeB.getLeafCount(branchB);
			int r = leafCountB - leafCountA;
			if(r == 0) {
				//number of leaves is the same, so decide based on total 
				//branch length
				double distanceA = nodeA.getTotalDistance(branchA);
				double distanceB = nodeB.getTotalDistance(branchB);
				if(distanceB < distanceA) {
					r = -1;
				} else if(distanceA < distanceB){
					r = 1;
				}
				
				if(leafCountA == 1 && leafCountB == 1) {
					//sort alphabetically
					r = nodeA.getNodeLabel(branchA).
					compareToIgnoreCase(nodeB.getNodeLabel(branchB));
				}
			}
			return r;
		}

	}

	/* (non-Javadoc)
	 * @see edu.vt.vbi.ci.nongraphic.tree.TreeNodeI#getParentBranch()
	 */
	public TreeBranch getParentBranch() {
		return parentBranch;
	}

	/* (non-Javadoc)
	 * @see edu.vt.vbi.ci.nongraphic.tree.TreeNodeI#setParentBranch(edu.vt.vbi.ci.nongraphic.tree.TreeBranch)
	 */
	public void setParentBranch(TreeBranch parentBranch) {
		this.parentBranch = parentBranch;
		TreeNodeI newParentNode = parentBranch.getOtherNode(this);
		this.setParentNode(newParentNode);
		TreeNodeI[] newChildNodes = getChildNodes();
		for(int i = 0; i < newChildNodes.length; i++) {
			newChildNodes[i].setParentNode(this);
		}
	}

	/* (non-Javadoc)
	 * @see edu.vt.vbi.ci.nongraphic.tree.TreeNodeI#getParentNode()
	 */
	public TreeNodeI getParentNode() {
		return parentNode;
	}

	/* (non-Javadoc)
	 * @see edu.vt.vbi.ci.nongraphic.tree.TreeNodeI#setParentNode(edu.vt.vbi.ci.nongraphic.tree.TreeNodeI)
	 */
	public void setParentNode(TreeNodeI parentNode) {
		this.parentNode = parentNode;
		TreeBranch[] branches = getAllBranches();
		for(int i = 0; i < branches.length; i++) {
			if(branches[i] == null) {
			}
			TreeNodeI otherNode = branches[i].getOtherNode(this);
			if(otherNode != parentNode) {
				otherNode.setParentBranch(branches[i]);
			}
		}
	}

	/**
	 * Compares based on the originalNodeString, alphabetically
	 */
	public int compareTo(Object arg0) {
		int r = this.getOriginalNodeString().compareTo(((TreeNode)arg0).
				getOriginalNodeString());
		return r;
	}
}
