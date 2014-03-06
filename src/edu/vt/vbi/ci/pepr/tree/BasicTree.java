package edu.vt.vbi.ci.pepr.tree;

import java.util.Arrays;

import edu.vt.vbi.ci.util.ExtendedBitSet;


/**
 * A simple representation of a phylogenetic tree.
 * @author enordber
 *
 */
public class BasicTree {

	public static final int NO_PARENT_INDICATOR = -1;

	private static String openParen = "(";
	private static String closeParen = ")";
	private static String comma = ",";
	private static String colon = ":";
	private static String openSquare = "[";
	private static String closeSquare = "]";
	private static char openParenChar = '(';
	private static char closeParenChar = ')';
	private static char commaChar = ',';
	private static char colonChar = ':';

	/*
	 * The original String representation of the tree/
	 */
	private String originalTreeString;

	/*
	 * All leaves (terminal nodes) in the tree.
	 */
	private String[] leaves;

	/*
	 * The original String representations for each node.
	 */
	private String[] nodeStrings;

	/*
	 * For each node, this contains the index of the parent node. 
	 * Nodes without parent nodes (such as a root) have a value of 
	 * NO_PARENT_INDICATOR.
	 * This array is the core data structure, and completely defines the
	 * tree topology when combined with the leaves array.
	 */
	private int[] nodeParentPointers;

	/*
	 * For each node, this contains the indices for child nodes. Terminal
	 * nodes (without children) have zero-length int[].
	 */
	private int[][] nodeChildPointers;

	/*
	 * Lengths for branches. The index here is the same as in 
	 * nodeParentPointers. A branch is defined by the more distal node.
	 * That is, the child node defines the branch, because each child node
	 * has only one parent node, but a parent node has multiple child nodes.
	 */
	private double[] branchLengths;

	/*
	 * Lengths for branches as cladogram.
	 */
	private double[] cladogramBranchLengths;

	/*
	 * Support values for branches. These are stored as Strings because there 
	 * may be multiple support values provided. The index here is the same as
	 * in nodeParentPointers. A branch is defined by the more distal node.
	 * That is, the child node define the branch, because each child node
	 * has only one parent node, but a parent node has multiple child nodes.
	 */
	private String[] branchSupports;

	/*
	 * Is the tree currently rooted?
	 */
	private boolean rooted;

	/*
	 * The index of the root node. This may be a valid index, even if
	 * the tree is not currently rooted. A previously rooted tree that gets
	 * unrooted will have rooted=false, but will have a valid index for 
	 * rootIndex.
	 */
	private int rootIndex = -1;


	public BasicTree(String treeString) {
		boolean looksLikeNewick = false;
		boolean looksLikeJSON = false;
		if(treeString.startsWith("(")) {
			looksLikeNewick = true;
		} else if(treeString.startsWith("{")) {
			looksLikeJSON = true;
		}
		if(looksLikeNewick) {
			parseNewickTreeString(treeString);
		} else if(looksLikeJSON) {
			parseJSONTreeString(treeString);
		}
	}

	/**
	 * No-argument constructor. Currently this is only used by the clone()
	 * method, so it is private.
	 */
	private BasicTree() {

	}

	private void parseJSONTreeString(String treeString) {

	}

	/**
	 * Parses the treeString to populate nodeParentPointers
	 * and leaves. This code is initially copied from
	 * TreeParser.parseTreeString().
	 * @param treeString
	 */
	private void parseNewickTreeString(String treeString) {
		treeString = treeString.trim();
		//remove ';' from the end, if it is there
		if(treeString.endsWith(";")) {
			treeString = treeString.substring(0,
					treeString.length()  - 1);
		}

		//get the tree string as a character array
		char[] treeChars = treeString.toCharArray();

		//count how many open parentheses and how many
		//close parentheses
		int openCount = 0;
		int closeCount = 0;
		for(int i = 0; i < treeChars.length; i++) {
			if(treeChars[i] == openParenChar) {
				openCount++;
			} else if(treeChars[i] == closeParenChar) {
				closeCount++;
			}
		}

		if(openCount != closeCount) {
			System.out.println("Unmatched parentheses in tree string. open: " 
					+ openCount
					+ ", close: " + closeCount);
		} else {

			//Find leaves
			//leaves are flanked by commas, or by (, or by ,)
			String commaString = ",";
			String[] commaSplits = treeString.split(commaString);
			leaves = new String[commaSplits.length];
			double[] leafBranchLengths = new double[leaves.length];
			for(int i = 0; i < commaSplits.length; i++) {
				boolean startsWithOpenParen =
						commaSplits[i].startsWith(openParen);
				boolean endsWithCloseParen =
						commaSplits[i].endsWith(closeParen);
				if(startsWithOpenParen && !endsWithCloseParen) {
					leaves[i] = 
							commaSplits[i].substring(commaSplits[i]
									.lastIndexOf(openParen) + 1);
				} else if(endsWithCloseParen && !startsWithOpenParen) {
					leaves[i] = 
							commaSplits[i].substring(0, 
									commaSplits[i].indexOf(closeParen));
				} else {
					leaves[i] = commaSplits[i];
					if(leaves[i].startsWith("(")) {
						System.out.println("BasicTree.parseNewickTreeString() leaf still has '(': " + leaves[i]);
					}
				}

				//if leaf contains branch distance info, remove it
				int indexOfColon = leaves[i].indexOf(colon);
				if(indexOfColon > -1) {
					String lengthString 
					= leaves[i].substring(indexOfColon+1);
					int indexOfParen = lengthString.indexOf(closeParen);
					if(indexOfParen > -1) {
						lengthString = lengthString.substring(0, 
								indexOfParen);
					}

					try {
						double leafBranchLength 
						= Double.parseDouble(lengthString);
						leaves[i] = leaves[i].substring(0, indexOfColon);
						leafBranchLengths[i] = leafBranchLength;
					} catch(NumberFormatException nfe) {
						nfe.printStackTrace();
						System.out.println("problem parsing length: " + lengthString + " for leaf: " + leaves[i]);
					}
				}

			}

			//for each leaf, find the index in the tree where the 
			//leaf starts
			int[][] leafNodePairs = new int[leaves.length][2];
			for(int i = 0; i < leaves.length; i++) {
				//The following line is causing a bug when there 
				//is one leaf name that is a substring of another
				//leaf name. If the leaf with the name that is a 
				//substring is after the leaf with the superstring
				//name, then the indexOf() search returns the
				//index for the wrong leaf.
				leafNodePairs[i][0] = treeString.indexOf(leaves[i]);
				leafNodePairs[i][1] = leafNodePairs[i][0] 
						+ leaves[i].length() - 1;
				//make sure correct index of the leaf string was found. 
				//If it is the correct index, then the next character
				//should be a comma, colon, or close parenthesis. 
				//This is to correct for the above mentioned bug.
				char followingChar = treeChars[leafNodePairs[i][1]+1];
				while(followingChar != commaChar
						&& followingChar != closeParenChar
						&& followingChar != colonChar) {
					//the index for a superstring of the leaf name was 
					//found.
					//Check for the next occurrence
					leafNodePairs[i][0] 
							= treeString.indexOf(leaves[i], 
									leafNodePairs[i][1]+1);
					leafNodePairs[i][1] = leafNodePairs[i][0] 
							+ leaves[i].length() - 1;
					followingChar = treeChars[leafNodePairs[i][1]+1];
				}
			}

			//tracks the locations open parentheses that have not yet
			//been paired with a closing parenthesis
			int[] openStack = new int[openCount];
			int[][] internalNodePairs = new int[openCount][2];

			//index of the next pair to be set
			int pairPointer = 0;

			//index of the current position in opens - next write
			//position
			int openStackPointer = 0;

			for(int i = 0; i < treeChars.length; i++) {
				if(treeChars[i] == openParenChar) {
					openStack[openStackPointer] = i;
					openStackPointer++;
				} else if(treeChars[i] == closeParenChar) {
					openStackPointer--;
					int openForPair = openStack[openStackPointer];
					int closeForPair = i;
					internalNodePairs[pairPointer][0] = openForPair;
					internalNodePairs[pairPointer][1] = closeForPair;
					pairPointer++;
				}
			}

			//nodePairs is a concatenation of leafNodePairs and
			//internalNodePairs
			int[][] nodePairs = 
					new int[leafNodePairs.length+internalNodePairs.length][];
			System.arraycopy(leafNodePairs, 0, nodePairs, 0, 
					leafNodePairs.length);
			System.arraycopy(internalNodePairs, 0, nodePairs, 
					leafNodePairs.length, internalNodePairs.length);

			//			for each node in nodePairs, gives the index for the 
			//nodes parent node (also an index in nodePairs)
			nodeParentPointers = new int[nodePairs.length];
			//support ond length are stored with the same index as
			//the child node of the branch. This is a good way to
			//do it because when the branches are created (see below)
			//they are done starting with the child node index
			branchLengths = new double[nodePairs.length];
			Arrays.fill(branchLengths, 0);
			branchSupports = new String[nodePairs.length];

			Arrays.fill(nodeParentPointers, NO_PARENT_INDICATOR);
			for(int i = 0; i < nodePairs.length; i++) {
				int nodeStart = nodePairs[i][0];
				int nodeEnd = nodePairs[i][1];

				//determine the branch length and support value for
				//the branch leading to this node, if these values
				//are provided.
				int firstCandidateIndex = nodeEnd+1;
				int nextCommaIndex = treeString.indexOf(comma, 
						firstCandidateIndex);
				int nextCloseParenIndex = treeString.indexOf(closeParen,
						firstCandidateIndex);
				int nextColonIndex = treeString.indexOf(colon,
						firstCandidateIndex);
				int treeEnd = treeString.length();
				//			String nodeString 
				//			= treeString.substring(nodeStart, nodeEnd);
				double branchLength = Double.NaN;
				//			double branchSupport = Double.NaN;
				int valueEnd = treeEnd;
				if(nextCommaIndex > -1) {
					valueEnd = Math.min(valueEnd, nextCommaIndex);
				}
				if(nextCloseParenIndex > -1) {
					valueEnd = Math.min(valueEnd, nextCloseParenIndex);
				}
				int supportEnd = -1;
				if(nextColonIndex >= 0) {
					supportEnd = Math.min(valueEnd, nextColonIndex);
				}
				if(supportEnd > firstCandidateIndex) {
					//there is a branch support value
					String branchSupportString = 
							treeString.substring(firstCandidateIndex,
									supportEnd);
					branchSupports[i] = branchSupportString;
				}
				if(nextColonIndex > nodeEnd && nextColonIndex < valueEnd) {
					//there is a length
					String branchLengthString = 
							treeString.substring(nextColonIndex+1, valueEnd);
					//see if support value is present, between square 
					//brackets after the branch length. This is an alternate 
					//way of encoding the support values in newick format
					int openSquareIndex = branchLengthString.indexOf(openSquare);
					if(openSquareIndex > -1) {
						int closeSquareIndex = branchLengthString.indexOf(closeSquare, openSquareIndex);
						String branchSupportString = branchLengthString.substring(openSquareIndex+1, closeSquareIndex);
						branchSupports[i] = branchSupportString;
						branchLengthString = branchLengthString.substring(0, openSquareIndex);
					}
					try {
						branchLength = Double.parseDouble(branchLengthString);
					} catch(NumberFormatException nfe) {
						System.out.println("problem parsing branch length: " 
								+ branchLengthString);
						System.out.println("using 0 for branch length. here is stack trace:");
						branchLength = 0;
						nfe.printStackTrace();
					}
				}
				if(!Double.isNaN(branchLength)) {
					branchLengths[i] = branchLength;
				}

				//look until a node is found such that this node
				//is between the start and end of that node. 
				//Because the internal nodes are parsed in a way
				//that finds the most internal nodes first,
				//the first node found containing another node 
				//should be the correct one.
				for(int j = 0; j < nodePairs.length; j++) {
					if(nodePairs[j][0] < nodeStart && 
							nodePairs[j][1] > nodeEnd) {
						nodeParentPointers[i] = j;
						j = nodePairs.length; 
						//break out of the loop here
					}
				}
			}

			//at this point, the tree topology can be completely represented by
			// the String[] leaves and the int[] nodeParentPointers
			//The indices in leaves map directly to parentPointers.
			//Additional indices in parentPointers correspond to 
			//internal nodes. Branch lengths and branch supports
			//are tracked in parallel arrays (if they are present
			//in the input tree string)

			nodeStrings = new String[nodeParentPointers.length];
			//populate nodeStrings
			for(int i = 0; i < nodeParentPointers.length; i++) {
				String nodeString = 
						treeString.substring(nodePairs[i][0], 
								nodePairs[i][1] + 1);
				nodeStrings[i] = nodeString;
			}
		}

		populateChildPointers();

		//determine if tree is rooted. If rooted there will be one entry
		//in nodeParentPointer that == NO_PARENT_INDICATOR
		//This node will also have exactly 2 child nodes
		int rootedAt = -1;
		int noParentsFound = 0;
		for(int i = 0; i < nodeParentPointers.length; i++) {
			if(nodeParentPointers[i] == NO_PARENT_INDICATOR) {
				rootedAt = i;
				noParentsFound++;
			}
		}

		if(noParentsFound == 1) {
			//			rootIndex = rootedAt;
			if(nodeChildPointers[rootedAt].length == 2) {
				rooted = true;
				rootIndex = rootedAt;
			} else {
				rooted = false;
			}
		}
	}

	private void populateChildPointers() {
		int[] nodeChildCounters = new int[nodeParentPointers.length];

		//count how many child nodes there are for each parent node
		for(int i = 0; i < nodeParentPointers.length; i++) {
			if(nodeParentPointers[i] > -1 &&
					nodeParentPointers[i] < nodeParentPointers.length) {
				nodeChildCounters[nodeParentPointers[i]]++;
			}
		}

		//if highest-level node has two child nodes, then consider the
		//tree t be rooted. Otherwise, it is unrooted.
		if(nodeChildCounters[nodeChildCounters.length-1] > 2) {
			rooted = false;
		}
		else {
			rooted = true;
		}

		//allocate nodeChildPointers arrays
		nodeChildPointers = new int[nodeChildCounters.length][];
		for(int i = 0 ; i < nodeChildPointers.length; i++) {
			nodeChildPointers[i] = new int[nodeChildCounters[i]];
		}

		//populate nodeChildPointers
		for(int i = 0; i < nodeParentPointers.length; i++) {
			int nodeIndex = nodeParentPointers[i];
			if(nodeIndex > -1 && nodeIndex < nodeParentPointers.length) {
				int[] childArray = nodeChildPointers[nodeIndex];
				int childIndex = childArray.length - nodeChildCounters[nodeIndex];
				childArray[childIndex] = i;
				nodeChildCounters[nodeIndex]--;
			}

		}
	}

	public String getTreeString(boolean includeLengths, boolean includeSupports) {
		String r = null;
		int topLevel = getTopLevelNode();
		if(topLevel == -1) {
			System.out.println("no top level nodes were found in this tree");
		} else {
			r = getNodeString(topLevel, includeLengths, includeSupports);
			r = r + ";";
		}
		return r;
	}

	public int getTopLevelNode() {
		int r = -1;
		//find out how many top-level nodes there are
		int topLevelNodes = 0;
		for(int i = 0; i < nodeParentPointers.length; i++) {
			if(nodeParentPointers[i] == NO_PARENT_INDICATOR &&
					nodeChildPointers[i].length > 0) {
				topLevelNodes++;
				r = i;
			}
		}
		if(r == -1) {
			System.out.println("no top level nodes were found in this tree");
		} 

		return r;
	}

	String getNodeString(int node, boolean includeLengths, boolean includeSupports) {
		String r = null;
		StringBuffer sb = new StringBuffer();
		if(nodeChildPointers[node].length == 0) {
			//this is a leaf node, return the nodeString
			sb.append(nodeStrings[node]);
		} else {
			//this is an internal node. Get the strings for its children
			String[] childStrings = new String[nodeChildPointers[node].length];
			for(int i = 0; i < childStrings.length; i++) {
				childStrings[i] = 
						getNodeString(nodeChildPointers[node][i], includeLengths, includeSupports);
			}
			if(childStrings.length > 1) {
				sb.append("(");
			}
			sb.append(childStrings[0]);
			if(includeLengths) {
				sb.append(":");
				sb.append(branchLengths[nodeChildPointers[node][0]]);
			}
			for(int i = 1; i < childStrings.length; i++) {
				sb.append(",");
				sb.append(childStrings[i]);
				if(includeLengths) {
					sb.append(":");
					sb.append(branchLengths[nodeChildPointers[node][i]]);
				}
			}
			if(childStrings.length > 1) {
				sb.append(")");
			}
			if(includeSupports) {
				if(branchSupports[node] != null) {
					sb.append(branchSupports[node]);
				}
			}
		}
		r = sb.toString();
		return r;
	}

	public String getTreeTopologyString() {
		String r = null;
		//collect top-level nodes

		//find out how many top-level nodes there are
		int topLevelNodes = 0;
		int topLevel = -1;
		for(int i = 0; i < nodeParentPointers.length; i++) {
			if(nodeParentPointers[i] == NO_PARENT_INDICATOR) {
				topLevelNodes++;
				topLevel = i;
			}
		}
		r = getNodeTopologyString(topLevel);
		r = r + ";";
		//		System.out.println("topLevelNodes: " + topLevelNodes);
		return r;
	}

	private String getNodeTopologyString(int node) {
		String r = null;
		if(nodeChildPointers[node].length == 0) {
			//this is a leaf node, return the nodeString
			r = nodeStrings[node];
		} else {
			//this is an internal node. Get the strings for its children
			String[] childStrings = new String[nodeChildPointers[node].length];
			for(int i = 0; i < childStrings.length; i++) {
				childStrings[i] = 
						getNodeTopologyString(nodeChildPointers[node][i]);
			}
			StringBuffer sb = new StringBuffer();
			sb.append("(");
			sb.append(childStrings[0]);
			for(int i = 1; i < childStrings.length; i++) {
				sb.append(",");
				sb.append(childStrings[i]);
			}
			sb.append(")");
			r = sb.toString();
		}
		return r;
	}

	public String[] getNodeStrings() {
		return nodeStrings;
	}

	public String[] getLeaves() {
		return leaves;
	}

	public int[] getNodeParentPointers() {
		return nodeParentPointers;
	}

	public int[][] getNodeChildPointers() {
		return nodeChildPointers;
	}

	public double[] getBranchLengths() {
		return branchLengths;
	}

	public void setBranchLengths(double[] branchLengths) {
		this.branchLengths = branchLengths;
	}

	public String[] getBranchSupportStrings() {
		return branchSupports;
	}

	public void setBranchSupportStrings(String[] supports) {
		branchSupports = supports;
	}

	public boolean isRooted() {
		boolean r = false;
		int topLevel = getTopLevelNode();

		r = rooted && topLevel >= 0 && nodeChildPointers[topLevel].length == 2;
		if(r) {
			rootIndex = topLevel;
		}
		return r;
	}

	public int getRootIndex() {
		return rootIndex;
	}

	/**
	 * Returns true if every internal node in this tree 
	 * has exactly 2 child nodes. If any multifurcations 
	 * are found, false is returned.
	 * 
	 * @return
	 */
	public boolean isBifurcating() {
		boolean r = true;
		int allowedMulti = 1; 
		if(isRooted()) {
			allowedMulti = 0;
		}
		int multi = 0;
		for(int i = 0; multi <= allowedMulti && i < nodeChildPointers.length; i++) {
			//multifurcating nodes will have 3 or mode child pointers.
			//leaf nodes have zero child pointers
			if(nodeChildPointers[i].length > 2) {
				multi++;
			}
		}
		r = multi <= allowedMulti;
		return r;
	}

	public Object clone() {
		BasicTree r = new BasicTree();
		r.branchLengths = new double[branchLengths.length];
		System.arraycopy(branchLengths, 0, r.branchLengths, 0, branchLengths.length);

		r.branchSupports = new String[branchSupports.length];
		System.arraycopy(branchSupports, 0, r.branchSupports, 0, branchSupports.length);

		r.leaves = new String[leaves.length];
		System.arraycopy(leaves, 0, r.leaves, 0, r.leaves.length);

		r.nodeStrings = new String[nodeStrings.length];
		System.arraycopy(nodeStrings, 0, r.nodeStrings, 0, nodeStrings.length);

		r.nodeParentPointers = new int[nodeParentPointers.length];
		System.arraycopy(nodeParentPointers, 0, r.nodeParentPointers, 0, nodeParentPointers.length);

		r.nodeChildPointers = new int[nodeChildPointers.length][];
		for(int i = 0; i < nodeChildPointers.length; i++) {
			r.nodeChildPointers[i] = new int[nodeChildPointers[i].length];
			System.arraycopy(nodeChildPointers[i], 0, r.nodeChildPointers[i], 0, nodeChildPointers[i].length);
		}

		r.rooted = rooted;
		r.rootIndex = rootIndex;

		return r;
	}

	/**
	 * Unroots the tree if it is rooted.
	 */
	void unroot() {
		if(isRooted()) {
//			System.out.println("BasicTree.unroot() tree is rooted");
			//find the root and remove it, making one child of the root
			//a child of the other child of the root
			int root = 0;
			for(int i = 0; i < nodeParentPointers.length; i++) {
				if(nodeParentPointers[i] == NO_PARENT_INDICATOR) {
					root = i;
					break;
				}
			}
			root = getTopLevelNode();
//			System.out.println("BasicTree.unroot() tree root is: " + root);

			int[] rootChildren = nodeChildPointers[root];
			//root Children should have length 2, otherwise root 
			//isn't really a root
			if(rootChildren.length != 2) {
				System.out.println("BasicTree.unroot() wrong number of rootChildren: " + rootChildren.length);
				System.exit(0);
			}
			//We'll make rootChildren[0] the parent of rootChildren[1].
			//This will make rootChildren[0] at least a trifurcation
			int newParent = rootChildren[0];
			int newChild = rootChildren[1];
			if(nodeChildPointers[newParent].length < 2) {
//				System.out.println("newParent only has one child, so reverse order");
				newParent = rootChildren[1];
				newChild = rootChildren[0];
			}

			nodeParentPointers[newChild] = newParent;
			int[] newChildren = new int[nodeChildPointers[newParent].length+1];
			System.arraycopy(nodeChildPointers[newParent], 0, newChildren,
					0, newChildren.length-1);
			newChildren[newChildren.length-1] = newChild;
			nodeChildPointers[newParent] = newChildren;

			//remove child pointers from defunct root
			nodeChildPointers[root] = new int[0];

			//the new parent now has no parent
			nodeParentPointers[newParent] = NO_PARENT_INDICATOR;

			//adjust branch length for newChild by adding the previous 
			//branch length for newParent
			branchLengths[newChild] += branchLengths[newParent];
			branchLengths[newParent] = 0;
			rooted = false;
		}
	}

	/**
	 * Root tree between nodes A and B. Any node that was a parent of A becomes
	 * a child of A. This change propagates until a node is reached with no parent.
	 * So, if D has no parent, and D is the parent of C, and C is the parent
	 * of A, then A becomes the parent of C, and C becomes the parent of D. The
	 * Root becomes the parent of A and C.
	 * 
	 * rootPoint should be a value between 0 and 1.0 indicating where 
	 * along the branch from nodeA to nodeB the root will be placed.
	 * rootPoint = 0   ->  root is at nodeA
	 * rootPoint = 1.0 ->  root is at nodeB
	 * rootPoint = 0.5 ->  root is in the middle of the nodeA-nodeB branch
	 * 
	 */
	public void rootBetweenNodes(int nodeA, int nodeB, double rootPoint) {
		//determine if nodeA is the parent of nodeB, or if nodeB is the
		//parent of nodeA. If neither of these is the case, then the tree
		//cannot be rooted between these two nodes, so exit without doing
		//anything (maybe change this to returning some error code)
		boolean nodeAIsParent = nodeParentPointers[nodeB] == nodeA;
		boolean nodeBIsParent = nodeParentPointers[nodeA] == nodeB;

		if(nodeAIsParent || nodeBIsParent) {
			int parentNode, childNode;
			if(nodeAIsParent) {
				parentNode = nodeA;
				childNode = nodeB;
			} else {
				rootPoint = 1d - rootPoint; //keep rootPoint relative to parentNode
				parentNode = nodeB;
				childNode = nodeA;
			}

			//ensure that the tree is unrooted. If it is not already unrooted,
			//then unroot it.
			unroot();

			//If a root node does not already exist, create a new node, 
			//which will be the root. Its parent will be 
			//NO_PARENT_INDICATOR. Its children will be nodeA and nodeB
			//for the node that was the parent, convert parents to children
			if(rootIndex == -1) {
				int[] newNodeParentPointers = 
						new int[nodeParentPointers.length+1];
				System.arraycopy(nodeParentPointers, 0, 
						newNodeParentPointers, 0, nodeParentPointers.length);
				newNodeParentPointers[nodeParentPointers.length] = 
						NO_PARENT_INDICATOR;
				nodeParentPointers = newNodeParentPointers;
				rootIndex = nodeChildPointers.length;

				int[][] newNodeChildPointers = 
						new int[nodeChildPointers.length+1][];
				System.arraycopy(nodeChildPointers, 0, newNodeChildPointers, 
						0, nodeChildPointers.length);
				newNodeChildPointers[nodeChildPointers.length] = 
						new int[]{parentNode,childNode};
				nodeChildPointers = newNodeChildPointers;

				//increase branchLengths array
				double[] newBranchLengths = new double[branchLengths.length+1];
				System.arraycopy(branchLengths, 0, newBranchLengths, 0, 
						branchLengths.length);
				newBranchLengths[branchLengths.length] = 0;
				branchLengths = newBranchLengths;

				//increase branchSupports array
				String[] newBranchSupports = new String[branchSupports.length+1];
				System.arraycopy(branchSupports, 0, newBranchSupports, 0, branchSupports.length);
				newBranchSupports[branchSupports.length] = "";
				branchSupports = newBranchSupports;
			}

			//add branches from the root to each node (nodeA and nodeB). 
			//These branches come from splitting the branch connecting nodeA
			//and nodeB to each other.
			double branchLength = branchLengths[childNode];
			double parentLength = branchLength * rootPoint;
			double childLength = branchLength * (1-rootPoint);
			String branchSupport = branchSupports[childNode];

			int grandparentNode = nodeParentPointers[parentNode];
			if(grandparentNode >=0 ) {
				convertParentsToChildren(grandparentNode, parentNode);
			}

			//remove childNode as a child of parentNode
			removeChildNode(parentNode, childNode);

			nodeParentPointers[parentNode] = rootIndex;
			nodeParentPointers[childNode] = rootIndex;
			branchLengths[parentNode] = parentLength;
			branchLengths[childNode] = childLength;
			branchSupports[parentNode] = branchSupport;
			nodeChildPointers[rootIndex] = new int[]{parentNode, childNode};

			rooted = true;
		}
	}

	private void removeChildNode(int parentNode, int childNode) {
		//remove childNode as a child of parentNode
		int[] newChildrenOfParent = new int[nodeChildPointers[parentNode].length-1];
		int newChildIndex = 0;
		for(int i = 0; i < nodeChildPointers[parentNode].length; i++) {
			if(nodeChildPointers[parentNode][i] != childNode) {
				newChildrenOfParent[newChildIndex] = nodeChildPointers[parentNode][i];
				newChildIndex++;
			}
		}

		nodeChildPointers[parentNode] = newChildrenOfParent;
		nodeParentPointers[childNode] = NO_PARENT_INDICATOR;

	}

	public int getNodeCount() {
		return nodeParentPointers.length;
	}

	/**
	 * This is a recursive method to reverse the direction of parent-child
	 * relationships
	 * @param node
	 */
	public void convertParentsToChildren(int parentNode, int childNode) {
		if(nodeParentPointers[parentNode] == NO_PARENT_INDICATOR) {
			//this is the end of the recursion. We have reached the node 
			//that has no parent
			nodeParentPointers[parentNode] = childNode;

		} else {
			//recursively call this method with the parentNode as the child
			//and the parentNode's parent as the parent.
			convertParentsToChildren(nodeParentPointers[parentNode], parentNode);

		}

		//remove childNode from the list of children of parentNode
		int[] newChildPointers = 
				new int[nodeChildPointers[parentNode].length-1];
		int index = 0;
		for(int i = 0; i < nodeChildPointers[parentNode].length; i++) {
			if(nodeChildPointers[parentNode][i] != childNode) {
				newChildPointers[index++] = nodeChildPointers[parentNode][i];
			}
		}
		nodeChildPointers[parentNode] = newChildPointers;

		//add parentNode as a child of childNode
		newChildPointers = new int[nodeChildPointers[childNode].length+1];
		System.arraycopy(nodeChildPointers[childNode], 0, newChildPointers, 0, newChildPointers.length-1);
		newChildPointers[newChildPointers.length-1] = parentNode;
		nodeChildPointers[childNode] = newChildPointers;

		//replace parent of parentNode with childNode
		nodeParentPointers[parentNode] = childNode;

		//the branch that used to lead to childNode now leads to parentNode
		//so update branchLengths
		branchLengths[parentNode] = branchLengths[childNode];

	}

	double[] getCladogramBranchLengths() {
		return cladogramBranchLengths;
	}

	void setCladogramBranchLengths(double[] cladogramBranchLengths) {
		this.cladogramBranchLengths = cladogramBranchLengths;
	}

	public void removeTaxon(String taxon) {
		//find index of the taxon
		int taxonIndex = -1;
		for(int i = 0; i < leaves.length; i++) {
			if(taxon.equalsIgnoreCase(leaves[i]) || TreeUtils.compressTaxonNameForComparison(taxon).equals(TreeUtils.compressTaxonNameForComparison(leaves[i]))) {
				taxonIndex = i;
			}
		}

		if(taxonIndex > -1) {
			int parent = nodeParentPointers[taxonIndex];

			if(parent == NO_PARENT_INDICATOR) {
				//the targeted taxon node has no parent node, so it is already essentially removed from the tree. 
				//all taxon nodes are leaf nodes, so they have no child nodes. This one has no parent node either,
				//so nothing in the tree points to it, so it will never be encountered in tree traversal.
			} else {
				//remove taxonIndex as a child node of parent
				removeChildNode(parent, taxonIndex);

				//if parent only has one child left, then parent is no longer a needed
				//node, so remove it and add the child directly to the grandparent
				if(nodeChildPointers[parent].length == 1) {
					int remainingChild = nodeChildPointers[parent][0];
					if(nodeParentPointers[parent] == NO_PARENT_INDICATOR) {
						System.out.println("making remaining child the new root");
						//remove parent and make the single child the new root.
						nodeParentPointers[remainingChild] = NO_PARENT_INDICATOR;
						nodeChildPointers[parent] = new int[0];
						rootIndex = remainingChild;
						rooted = true;
					} else {
						int grandparent = nodeParentPointers[parent];
						
						double branchLength = branchLengths[parent];
						branchLengths[remainingChild] += branchLength;
						//find which child of grandparent parent is
						int parentAsChildIndex = -1;
						for(int i = 0; i < nodeChildPointers[grandparent].length; i++) {
							if(nodeChildPointers[grandparent][i] == parent) {
								parentAsChildIndex = i;
								break;
							}
						}
						//replace parent with remainingChild as a child of grandparent
						nodeChildPointers[grandparent][parentAsChildIndex] = remainingChild;
						//point remainingChild to its new parent (formerly grandparent)
						nodeParentPointers[remainingChild] = grandparent;
					}
				} else {
					if(!isRooted()) {
						//just removed the only multifurcation in an unrooted tree
						//turn it into a rooted tree with parent as the root node
						rooted = true;
						rootIndex = parent;
					}
				}
			}
		}
	}

	private int[] getDescendantTips(int node) {
		int[] r = null;
		if(nodeChildPointers[node].length == 0) {
			//this is the end of the recursion - we have reached a tip
			r = new int[]{node};
		} else {
			r = new int[0];
			int[] descendants;
			int[] children = nodeChildPointers[node];
			for(int i = 0; i < children.length; i++) {
				descendants = getDescendantTips(children[i]);
				int[] holder = new int[r.length+descendants.length];
				System.arraycopy(r, 0, holder, 0, r.length);
				System.arraycopy(descendants, 0, holder, r.length, descendants.length);
				r = holder;
			}
		}
		return r;
	}

	/**
	 * Replaces the subtree below the specified node with the given tree.
	 * The newSubTree should be rooted, and the largest child node of the root
	 * will be swapped for the replaceBelowNode.
	 * 
	 * @param replaceBelowNode
	 * @param newSubTree
	 */
	public BasicTree replaceSubtreeBelow(int replaceBelowNode, BasicTree newSubTree) {
		BasicTree r = null;

		//create new arrays large enough to hold the current tree and the 
		//new subTree. Every node index in newSubTree will be offset by the 
		//number of nodes in the current tree

		int newSubRoot = newSubTree.getRootIndex();
		int[] rootChildren = newSubTree.nodeChildPointers[newSubRoot];
		int maxChild = -1;
		int maxSize = -1;
		int minSize = Integer.MAX_VALUE;
		int minChild = -1;
		int descendantsOfRoot = newSubTree.getDescendantTips(newSubRoot).length;
		for(int i = 0; i < rootChildren.length; i++) {
			int descendantCount = newSubTree.getDescendantTips(rootChildren[i]).length;
			if(descendantCount > maxSize) {
				maxSize = descendantCount;
				maxChild = rootChildren[i];
			}
			if(descendantCount < minSize) {
				minSize = descendantCount;
				minChild = rootChildren[i];
			}
		}

		int offset = nodeParentPointers.length;
		int newSubtreeNodeCount = newSubTree.nodeParentPointers.length;
		int newArraySize = offset + newSubtreeNodeCount;

		String[] newLeaves = new String[newArraySize];
		String[] newNodeStrings = new String[newArraySize];
		int[] newNodeParentPointers = new int[newArraySize];
		int[][] newNodeChildPointers = new int[newArraySize][];
		double[] newBranchLengths = new double[newArraySize];
		String[] newBranchSupports = new String[newArraySize];

		System.arraycopy(leaves, 0, newLeaves, 0, leaves.length);
		System.arraycopy(newSubTree.leaves, 0, newLeaves, offset, 
				newSubTree.leaves.length);

		System.arraycopy(nodeStrings, 0, newNodeStrings, 0, nodeStrings.length);
		System.arraycopy(newSubTree.nodeStrings, 0, newNodeStrings, offset, 
				newSubTree.nodeStrings.length);

		System.arraycopy(nodeParentPointers, 0, newNodeParentPointers, 0, nodeParentPointers.length);
		System.arraycopy(newSubTree.nodeParentPointers, 0, newNodeParentPointers, offset, 
				newSubTree.nodeParentPointers.length);
		for(int i = offset; i < newNodeParentPointers.length; i++) {
			newNodeParentPointers[i] +=offset;
		}

		System.arraycopy(nodeChildPointers, 0, newNodeChildPointers, 0, nodeChildPointers.length);
		System.arraycopy(newSubTree.nodeChildPointers, 0, newNodeChildPointers, offset, 
				newSubTree.nodeChildPointers.length);
		for(int i = offset; i < newNodeChildPointers.length; i++) {
			for(int j = 0; j < newNodeChildPointers[i].length; j++) {
				newNodeChildPointers[i][j] += offset;
			}
		}

		System.arraycopy(branchLengths, 0, newBranchLengths, 0, branchLengths.length);
		System.arraycopy(newSubTree.branchLengths, 0, newBranchLengths, offset, 
				newSubTree.branchLengths.length);

		//use the old branch length to the replaceNode instead of the
		//new branch length
		newBranchLengths[maxChild+offset] = newBranchLengths[replaceBelowNode];

		System.arraycopy(branchSupports, 0, newBranchSupports, 0, branchSupports.length);
		System.arraycopy(newSubTree.branchSupports, 0, newBranchSupports, offset, 
				newSubTree.branchSupports.length);


		//have the parent node that points to replaceBelowNode change so it 
		//points to the 'maxChild' node in the subtree
		int parent = nodeParentPointers[replaceBelowNode];
		int[] childrenOfParent = newNodeChildPointers[parent];
		int[] newChildrenOfParent = new int[childrenOfParent.length];
		for(int i = 0; i < childrenOfParent.length; i++) {
			if(childrenOfParent[i] == replaceBelowNode) {
				newChildrenOfParent[i] = maxChild + offset;
			} else {
				newChildrenOfParent[i] = childrenOfParent[i];
			}
		}
		newNodeChildPointers[parent] = newChildrenOfParent;

		newNodeChildPointers[replaceBelowNode] = newNodeChildPointers[maxChild+offset];

		r = (BasicTree) this.clone();
		r.branchLengths = newBranchLengths;
		r.branchSupports = newBranchSupports;
		r.leaves = newLeaves;
		r.nodeStrings = newNodeStrings;
		r.nodeParentPointers = newNodeParentPointers;
		r.nodeChildPointers = newNodeChildPointers;

		r = new BasicTree(r.getTreeString(true, true));

		return r;
	}

	public double[][] getDistanceMatrix() {
		double[][] r = null;

		int[] nodeParentPointers = this.getNodeParentPointers();
		String[] leaves = this.getLeaves();
		double[] branchLengths = this.getBranchLengths();

		ExtendedBitSet[] nodeAncestors = new ExtendedBitSet[leaves.length];
		double[][] ancestorDistances = new double[nodeAncestors.length][nodeParentPointers.length];
		for(int i = 0; i < nodeAncestors.length; i++) {
			nodeAncestors[i] = new ExtendedBitSet();
			int parent = nodeParentPointers[i];
			double cumulativeDistance = branchLengths[i];

			while(parent >= 0) {
				ancestorDistances[i][parent] = cumulativeDistance;
				nodeAncestors[i].set(parent);
				cumulativeDistance += branchLengths[parent];
				parent = nodeParentPointers[parent];
			}
		}

		//get the distances between each pair of leaves, by finding their
		//closest common ancestor and adding together the ancestorDistances
		r = new double[leaves.length][leaves.length];
		for(int i = 0; i < leaves.length; i++) {
			for(int j = i+1; j < leaves.length; j++) {
				//get common ancestors
				ExtendedBitSet commonAncestors = nodeAncestors[i].getAnd(nodeAncestors[j]);
				int nearestCommonAncestor = commonAncestors.nextSetBit(0);
				double distance = ancestorDistances[i][nearestCommonAncestor] + ancestorDistances[j][nearestCommonAncestor];
				r[i][j] = distance;
				r[j][i] = distance;
			}
		}

		return r;
	}

	public double getTotalTreeLength() {
		double r = 0;
		double[] branchLengths = getBranchLengths();
		for(int i= 0; i < branchLengths.length; i++) {
			if(!Double.isNaN(branchLengths[i])) {
				r += branchLengths[i];
			}
		}
		return r;
	}
	
	public String getTreeJSON() {
		String r = null;
		int topNode = getTopLevelNode();
		r = getNodeJSON(topNode, 0);
		return r;
	}
	
	public String getNodeJSON(int node, double distanceFromRoot) {
		String r = null;
		
		StringBuffer sb = new StringBuffer();
		sb.append("{");
		sb.append("\"name\":");
		if(nodeChildPointers[node] == null || nodeChildPointers[node].length == 0) {
			//this is a leaf node
			sb.append("\"");
			sb.append(leaves[node]); 
			sb.append("\"");
		} else {
		    //this is an internal node
			sb.append("\"");
			sb.append(node); //use node number for name of internal nodes
			sb.append("\"");
			sb.append(",\"children\":[");
			sb.append(getNodeJSON(nodeChildPointers[node][0], distanceFromRoot + branchLengths[nodeChildPointers[node][0]]));
			for(int i = 1; i < nodeChildPointers[node].length; i++) {
				sb.append(",");
				sb.append(getNodeJSON(nodeChildPointers[node][i], distanceFromRoot + branchLengths[nodeChildPointers[node][0]]));
			}
			sb.append("]");
		}
		sb.append(", \"toRoot\":");
		sb.append(distanceFromRoot);
		//if there is metadata, add it
		sb.append("}");
		r = sb.toString();
		return r;
	}
	
	public boolean containsTaxon(String taxon) {
		boolean r = false;
		String[] leaves = getLeaves();
		if(leaves == null) {
			System.out.println("Freak Out - leaves is null for this tree: " + originalTreeString);
		}
		for(int i = 0; !r && i < leaves.length; i++) {
			r = taxon.equalsIgnoreCase(leaves[i]);
		}
		return r;
	}

}
