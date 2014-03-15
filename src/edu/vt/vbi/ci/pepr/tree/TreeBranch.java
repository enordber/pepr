package edu.vt.vbi.ci.pepr.tree;

import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;

public class TreeBranch {
	
	private TreeNode nodeA;
	private TreeNode nodeB;
	private double branchLength = 1;
	private double branchSupport = Double.NaN;
	private boolean isRoot = false;
	
	/**
	 * The location along the branch where the root is
	 * located, if this is the root branch. This should
	 * be a value between 0.0 (at nodeA) and 1.0 (at nodeB).
	 */
	private double rootPoint;
	
	/*
	 * Newick format String of subtree represented by this
	 * branch.
	 */
	private String branchString;
	private boolean branchStringNeedsRebuild = true;

	/**
	 * @param nodeA
	 * @param nodeB
	 */
	public TreeBranch(TreeNode nodeA, TreeNode nodeB) {
		this.nodeA = nodeA;
		this.nodeB = nodeB;
		
		nodeA.addBranch(this);
		nodeB.addBranch(this);
	}

	/**
	 * 
	 * @param nodeA
	 * @param nodeB
	 * @param length
	 */
	public TreeBranch(TreeNode nodeA, TreeNode nodeB, 
			double length, double support) {
		this(nodeA, nodeB);
//		System.out.println("new TreeBranch() length: " + length);
		if(!Double.isNaN(length)) {
		    branchLength = length;
		}
		branchSupport = support;
//		System.out.println("branchLength: " + branchLength + "\tbranchSupport: "
//				+ branchSupport);
	}
	
	public double getBranchLength() {
		return branchLength;
	}

	public void setBranchLength(double branchLength) {
		System.out.println("TreeBranch.setBranchLength() from "
				+ this.branchLength + " to " + branchLength);
		this.branchLength = branchLength;
	}

	public double getBranchSupport() {
		return branchSupport;
	}

	public void setBranchSupport(double branchSupport) {
		this.branchSupport = branchSupport;
	}
	
	/**
	 * Returns the node for this Branch that is not the 
	 * given node. If the given node is the parent node, 
	 * then the child node is returned, If the given
	 * node is the child node, then the parent node is
	 * returned. If the given node is neither the parent
	 * or child node for this branch, then null is 
	 * returned.
	 * 
	 * @param aNode
	 * @return
	 */
	public TreeNode getOtherNode(TreeNodeI aNode) {
		TreeNode r = null;
		if(aNode == nodeA) {
			r = nodeB;
		} else if(aNode == nodeB) {
			r = nodeA;
		}
		return r;
	}

	/**
	 * Returns an aray containing both nodes for
	 * this Branch.
	 * 
	 * @return
	 */
	public TreeNode[] getNodes() {
		return new TreeNode[]{nodeA, nodeB};
	}

	/**
	 * Returns true if this branch is currently the root branch.
	 * @return
	 */
	public boolean isRoot() {
		return isRoot;
	}
	
	/**
	 * Makes this branch serve as the root for the tree.
	 * @param rootPoint Location along branch for root.
	 *                  This should be a value between
	 *                   0.0 (at nodeA) and 1.0 (at nodeB).
	 *
	 */
	public void setAsRoot(double rootPoint) {
		isRoot = true;
		branchStringNeedsRebuild = true;
		this.rootPoint = rootPoint;
		//buildBranchString();
	}
	
	/**
	 * Sets the isRoot value to false for this branch.
	 * If this Branch had been the root, it will no longer
	 * be the root after calling this method.
	 *
	 */
	void setAsNotRoot() {
		isRoot = false;
		branchStringNeedsRebuild = true;
	}
	
	
	public String getLeafList(TreeNodeI fromNode) {
		String r = null;
		if(fromNode == nodeA) {
		    r = nodeB.getLeafList(this);
		} else if(fromNode == nodeB) {
			r = nodeA.getLeafList(this);
		}
		return r;
	}

     public Collection getDescendantBranches(TreeNodeI fromNode) {
    	 HashSet r = new HashSet();
    	 TreeNodeI childNode = getOtherNode(fromNode);
    	 if(childNode != null) {
    		 TreeBranch[] childBranches = childNode.getOtherBranches(this);
    		 for(int i = 0; i < childBranches.length; i++) {
    			 r.add(childBranches[i]);
    			 r.addAll(childBranches[i].getDescendantBranches(childNode));
    		 }
    	 }
    	 return r;
     }
     
     private String buildBranchString(TreeNodeI fromNode, boolean includeDistance, boolean includeSupport) {
    	 //System.out.println(">TreeBranch.buildBranchString() fromNode: " +
    	//		 fromNode);
    	 branchStringNeedsRebuild = false;
    	 StringBuffer sb = new StringBuffer();
    	 if(isRoot || fromNode == null) {
    		 String nodeAString = nodeA.toString(this, includeDistance, includeSupport);
    		 String nodeBString = nodeB.toString(this, includeDistance, includeSupport);
    		 
//    		 System.out.println("A: " + nodeAString + "\t\tB: " + nodeBString);
    		 sb.append("(");
    		 sb.append(nodeAString);
    		 sb.append(",");
    		 sb.append(nodeBString);
    		 sb.append(")");

    	 } else {
    		 String childNodeString = getOtherNode(fromNode).toString(this, includeDistance, includeSupport);
    		 sb.append(childNodeString);
    	 }
    	 //add support, if it is requested and known
    	 if(includeSupport && !Double.isNaN(branchSupport)) {
    		 sb.append(branchSupport);
    	 }
    	 //add length, if it is requested and known
    	 if(includeDistance && !Double.isNaN(branchLength)) {
    		 sb.append(":");
    		 sb.append(branchLength);
    	 }
    	 branchString = sb.toString();
    	 //System.out.println("<TreeBranch.buildBranchString(): " + branchString);
    	 return branchString;
     }
     
     public String toString(TreeNodeI fromNode, boolean includeDistance, boolean includeSupport) {
        	 return buildBranchString(fromNode, includeDistance, includeSupport);
     }
     
     public String getToplogyString(TreeNodeI fromNode) {
    	     return buildBranchString(fromNode, false, false);
     }
     
     public double getTotalDistance(TreeNode fromNode) {
    	     double r = 0;
    	     TreeNode otherNode = getOtherNode(fromNode);
    	     r = this.branchLength + otherNode.getTotalDistance(this);
    	     if(isRoot) {
    	    	     r += fromNode.getTotalDistance(this);
    	     }
    	     return r;
     }
     
     /**
      * Returns the sets of leaf nodes that occur on either
      * side of this branch. That is, if you cut the tree 
      * at this branch, it would split it into two sets of 
      * leaves. The returned value contains these two sets
      * of leaves as a 2 dimension String array.
      * 
      * @return
      */
     String[][] getBipartition() {
    	 String[][] r = new String[2][];
    	 TreeNodeI[] childNodes = getNodes();
    	 if(childNodes[0].isLeaf()) {
    		 r[0] = new String[]{childNodes[0].getOriginalNodeString()};
    	 } else {
    		 TreeNodeI[] descendants = childNodes[0].getDescendantLeaves(this);
    		 r[0] = new String[descendants.length];
    		 for(int i = 0; i < descendants.length; i++) {
    			 r[0][i] = descendants[i].getOriginalNodeString();
    		 }
    	 }
    	 if(childNodes[1].isLeaf()) {
    		 r[1] = new String[]{childNodes[1].getOriginalNodeString()};
    	 } else {
    		 TreeNodeI[] descendants = childNodes[1].getDescendantLeaves(this);
    		 r[1] = new String[descendants.length];
    		 for(int i = 0; i < descendants.length; i++) {
    			 r[1][i] = descendants[i].getOriginalNodeString();
    		 }
    	 }
    	 return r;
     }
     
     HashSet[] getBipartitionLeafSets() {
    	 HashSet[] r = new HashSet[2];
    	 String[][] biparts = getBipartition();
    	 for(int i = 0; i < 2; i++) {
    		 r[i] = new HashSet();
    		 for(int j = 0; j < biparts[i].length; j++) {
    			 r[i].add(biparts[i][j]);
    		 }
    	 }
    	 return r;
     }
     
     Bipartition getBipartition(HashMap idToTaxonMap, String[] taxonNames) {
    	 Bipartition r = null;
    	 
    	 return r;
     }

	public double getDistanceFromRoot(TreeNodeI fromNode) {
		double r = branchLength;
		if(this.isRoot()) {
			if(fromNode == nodeA) {
				r *= rootPoint;
			} else if(fromNode == nodeB) {
				r *= (1.0 - rootPoint);
			}
			
		} else {
			TreeNodeI otherNode = getOtherNode(fromNode);
			r += otherNode.getDistanceToRoot(this);
		}
		return r;
	}
}
