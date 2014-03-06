package edu.vt.vbi.ci.pepr.tree;

import java.util.ArrayList;
import java.util.Arrays;

public class TreeParser {
	
	private static String openParen = "(";
	private static String closeParen = ")";
	private static String openParenRegX = "\\(";
	private static String closeParenRegX = "\\)";
	private static String comma = ",";
	private static String colon = ":";
	private static char openParenChar = '(';
	private static char closeParenChar = ')';
	private static char commaChar = ',';
	private static char colonChar = ':';

	/**
	 * Parses a newick format tree string and returns
	 * a Tree object.
	 * 
	 * @param treeString
	 * @return
	 */
	public static TreeI parseTreeString(String treeString) {
		TreeI r = null;
		
		BasicTree basicTree = new BasicTree(treeString);
		String[] leaves = basicTree.getLeaves();
		int[] parentPointers = basicTree.getNodeParentPointers();
		double[] branchLengths = basicTree.getBranchLengths();
		String[] branchSupportStrings = basicTree.getBranchSupportStrings();
		double[] branchSupports = new double[branchSupportStrings.length];
		for(int branchIndex = 0; branchIndex < branchSupports.length; branchIndex++) {
//			System.out.println("branchSupportStrings " + branchIndex + ": " +branchSupportStrings[branchIndex]);
    		if(branchSupportStrings[branchIndex] == null) {
    			branchSupports[branchIndex] = Double.NaN;
    		} else {
			branchSupports[branchIndex] = 
				Double.parseDouble(branchSupportStrings[branchIndex]);
			}
		}
		
		//at this point, the tree can be completely represented by
		// the String[] leaves and the int[] parentPointers
		//The indices in leaves map directly to parentPointers.
		//Additional indices in parentPointers correspand to 
		//internal nodes. Branch lengths and branch supports
		//are tracked in parallel arrays (if they are present
		//in the input tree string)
		TreeNode[] nodes = new TreeNode[parentPointers.length];
		ArrayList branches = new ArrayList(parentPointers.length);
		
		String[] nodeStrings = basicTree.getNodeStrings();
		//first create TreeNodes
		for(int i = 0; i < parentPointers.length; i++) {
			nodes[i] = new TreeNode(nodeStrings[i]);
		}
		
		//now that TreeNodes have been created, create 
		//TreeBranches to connect the TreeNodes
		int rootNodeIndex = -1;
		//find the root node index
		for(int i = 0; 
		    i < nodes.length && rootNodeIndex == -1;
		    i++) {
			if(parentPointers[i] == BasicTree.NO_PARENT_INDICATOR) {
				rootNodeIndex = i;
			}
		}
		
		//while connecting the nodes, count the
		//number of top-level nodes.
		int topLevelNodeCount = 0;
		TreeNode rootFlankA  = null;
		TreeNode rootFlankB = null;
		double rootLengthA = Double.NaN;
		double rootLengthB = Double.NaN;
		double rootSupportA = Double.NaN;
		double rootSupportB = Double.NaN;
		
		for(int i = 0; i < nodes.length; i++) {
			//see if this will be a root flanking node
			if(parentPointers[i] == rootNodeIndex) {
				topLevelNodeCount++;
//				System.out.println("topLevelNode " + topLevelNodeCount 
//						+ " encountered: " + 
//						nodes[i].getOriginalNodeString());
				if(rootFlankA == null) {
					rootFlankA = nodes[i];
					rootLengthA = branchLengths[i];
					rootSupportA = branchSupports[i];
				} else if(rootFlankB == null){
					rootFlankB = nodes[i];
					rootLengthB = branchLengths[i];
					rootSupportB = branchSupports[i];
				} else {
					
				}
			} else if(parentPointers[i] == BasicTree.NO_PARENT_INDICATOR) {
				//skip this - it's the 'root node', which I'm not
				//using anymore. This should ultimately be removed
			} else {
				TreeNode nodeA = nodes[i];
				TreeNode nodeB = nodes[parentPointers[i]];
				branches.add(new TreeBranch(nodeA, nodeB, 
						branchLengths[i], branchSupports[i]));
			}
		}
		
		//collect the top-level nodes, branch lengths, and
		//support values
		TreeNode[] topLevelNodes = new TreeNode[topLevelNodeCount];
		double[] topLevelLengths = new double[topLevelNodeCount];
		double[] topLevelSupports = new double[topLevelNodeCount];
		int topLevelNodeIndex = 0;
		for(int i = 0; i < nodes.length; i++) {
			if(parentPointers[i] == rootNodeIndex) {
				topLevelNodes[topLevelNodeIndex] = nodes[i];
				topLevelLengths[topLevelNodeIndex] 
				                = branchLengths[i];
				topLevelSupports[topLevelNodeIndex] 
				                 = branchSupports[i];
				topLevelNodeIndex++;
			}
		}

		//Trees may have either two ar three top-level nodes. 
		//In the case of two, connect them by a branch, which 
		//will be the Root Branch for the Tree. In the case 
		//of three, create a new internal node, and create 
		//three branches to this node. One of these branches
		//will be chosen as the Root Branch, with the root 
		//point = 1.0 (all the way at the new internal node).
		//This is basically, pseudo-rooting the tree.
		TreeBranch rootBranch = null;
		if(topLevelNodeCount == 2) {
			//Create root branch connecting the two flanking nodes
			//root branch length will be the sum of the two parts
			double rootBranchLength = rootLengthA + rootLengthB;
			//root point will be the proportion of the branch length
			//assigned to rootFlankA
			double rootPoint = rootLengthA / rootBranchLength;

			//root branch support will mean the mean of the two parts
			double rootBranchSupport 
			            = (rootSupportA + rootSupportB) / 2;
			rootBranch = new TreeBranch(rootFlankA, rootFlankB,
					rootBranchLength, rootBranchSupport);
			rootBranch.setAsRoot(rootPoint);
		} else if(topLevelNodeCount == 3) {
			TreeNode iNode = new TreeNode();
			TreeBranch rootBranch1 
			           = new TreeBranch(iNode, topLevelNodes[0], 
					topLevelLengths[0], topLevelSupports[0]);
			TreeBranch rootBranch2 
			           = new TreeBranch(iNode, topLevelNodes[1],
			         topLevelLengths[1], topLevelSupports[1]);
			rootBranch = new TreeBranch(topLevelNodes[2], iNode,
					topLevelLengths[2], topLevelSupports[2]);
			rootBranch.setAsRoot(1.0);
		} else {
			System.out.println("TreeParser.parseTreeString() tree has "
					+ topLevelNodeCount + " top level nodes. " +
							"Only 2 and 3 are valid.");
		}
		try {
		 r = new Tree(rootBranch);
		} catch(NullPointerException npe) {
			npe.printStackTrace();
		}

		return r;
	}
	
	/**
	 * Parses a newick format tree string and returns
	 * a Tree object.
	 * 
	 * @param treeString
	 * @return
	 */
	public static TreeI parseTreeString2(String treeString) {
		TreeI r = null;
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
			System.out.println("Unmatched parentheses. open: " 
					+ openCount
					+ ", close: " + closeCount);
			System.exit(0);
		}

		//Find leaves
		//leaves are flanked by commas, or by (, or by ,)
		String commaString = ",";
		String[] commaSplits = treeString.split(commaString);
		String[] leaves = new String[commaSplits.length];
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
				double leafBranchLength 
						= Double.parseDouble(lengthString);
				leaves[i] = leaves[i].substring(0, indexOfColon);
				leafBranchLengths[i] = leafBranchLength;
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
				//Check for the next ocurrance
				leafNodePairs[i][0] 
				                 = treeString.indexOf(leaves[i], 
				                		 leafNodePairs[i][1]);
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
				
		String[] nodeStrings = new String[internalNodePairs.length];
		for(int i = 0; i < internalNodePairs.length; i++) {
			nodeStrings[i] = 
				treeString.substring(internalNodePairs[i][0],
						internalNodePairs[i][1]+1);
		}
		
		//nodePairs is a concatenation of leafNodePairs and
		//internalNodePairs
		int[][] nodePairs = 
			new int[leafNodePairs.length+internalNodePairs.length][];
		System.arraycopy(leafNodePairs, 0, nodePairs, 0, 
				leafNodePairs.length);
		System.arraycopy(internalNodePairs, 0, nodePairs, 
				leafNodePairs.length, internalNodePairs.length);
		
//		for each node in nodePairs, gives the index for the 
		//nodes parent node (also an index in nodePairs)
		int[] parentPointers = new int[nodePairs.length];
		//support ond length are stored with the same index as
		//the child node of the branch. This si a good way to
		//do it because when the branches are created (see below)
		//they are done starting with the child node index
		double[] branchLengths = new double[nodePairs.length];
		Arrays.fill(branchLengths, 1);
		double[] branchSupports = new double[nodePairs.length];
		
		int rootIndicator = -1;
		Arrays.fill(parentPointers, rootIndicator);
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
			String nodeString 
			        = treeString.substring(nodeStart, nodeEnd);
			double branchLength = Double.NaN;
			double branchSupport = Double.NaN;
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
				branchSupport 
				      = Double.parseDouble(branchSupportString);
			}
			branchSupports[i] = branchSupport;
			if(nextColonIndex > nodeEnd && nextColonIndex < valueEnd) {
				//there is a length
				String branchLengthString = 
					treeString.substring(nextColonIndex+1, valueEnd);
				branchLength = Double.parseDouble(branchLengthString);
			}
			if(!Double.isNaN(branchLength)) {
			    branchLengths[i] = branchLength;
			}

			//look until a node is found such that this node
			//is between the start and end of that node. 
			//Because the internal nodes are parsed in a way
			//that finds the most internal nodes first,
			//the first node found countaining another node 
			//should be the correct one.
			for(int j = 0; j < nodePairs.length; j++) {
				if(nodePairs[j][0] < nodeStart && 
						nodePairs[j][1] > nodeEnd) {
					parentPointers[i] = j;
					j = nodePairs.length; 
					//break out of the loop here
				}
			}
		}
		
		//at this point, the tree can be completely represented by
		// the String[] leaves and the int[] parentPointers
		//The indices in leaves map directly to parentPointers.
		//Additional indices in parentPointers correspond to 
		//internal nodes. Branch lengths and branch supports
		//are tracked in parallel arrays (if they are present
		//in the input tree string)

		TreeNode[] nodes = new TreeNode[parentPointers.length];
		ArrayList branches = new ArrayList(parentPointers.length);
		
		//first create TreeNodes
		for(int i = 0; i < parentPointers.length; i++) {
			String nodeString = 
				treeString.substring(nodePairs[i][0], 
						nodePairs[i][1] + 1);
			nodes[i] = new TreeNode(nodeString);
		}
		
		//now that TreeNodes have been created, create 
		//TreeBranches to connect the TreeNodes
		int rootNodeIndex = -1;
		//find the root node index
		for(int i = 0; 
		    i < nodes.length && rootNodeIndex == -1;
		    i++) {
			if(parentPointers[i] == rootIndicator) {
				rootNodeIndex = i;
			}
		}
		
		//while connecting the nodes, count the
		//number of top-level nodes.
		int topLevelNodeCount = 0;
		TreeNode rootFlankA  = null;
		TreeNode rootFlankB = null;
		double rootLengthA = Double.NaN;
		double rootLengthB = Double.NaN;
		double rootSupportA = Double.NaN;
		double rootSupportB = Double.NaN;
		
		for(int i = 0; i < nodes.length; i++) {
			//see if this will be a root flanking node
			if(parentPointers[i] == rootNodeIndex) {
				topLevelNodeCount++;
//				System.out.println("topLevelNode " + topLevelNodeCount 
//						+ " encountered: " + 
//						nodes[i].getOriginalNodeString());
				if(rootFlankA == null) {
					rootFlankA = nodes[i];
					rootLengthA = branchLengths[i];
					rootSupportA = branchSupports[i];
				} else if(rootFlankB == null){
					rootFlankB = nodes[i];
					rootLengthB = branchLengths[i];
					rootSupportB = branchSupports[i];
				} else {
					
				}
			} else if(parentPointers[i] == rootIndicator) {
				//skip this - it's the 'root node', which I'm not
				//using anymore. This should ultimately be removed
			} else {
				TreeNode nodeA = nodes[i];
				TreeNode nodeB = nodes[parentPointers[i]];
				branches.add(new TreeBranch(nodeA, nodeB, 
						branchLengths[i], branchSupports[i]));
			}
		}
		
		//collect the top-level nodes, branch lengths, and
		//support values
		TreeNode[] topLevelNodes = new TreeNode[topLevelNodeCount];
		double[] topLevelLengths = new double[topLevelNodeCount];
		double[] topLevelSupports = new double[topLevelNodeCount];
		int topLevelNodeIndex = 0;
		for(int i = 0; i < nodes.length; i++) {
			if(parentPointers[i] == rootNodeIndex) {
				topLevelNodes[topLevelNodeIndex] = nodes[i];
				topLevelLengths[topLevelNodeIndex] 
				                = branchLengths[i];
				topLevelSupports[topLevelNodeIndex] 
				                 = branchSupports[i];
				topLevelNodeIndex++;
			}
		}

		//Trees may have either two ar three top-level nodes. 
		//In the case of two, connect them by a branch, which 
		//will be the Root Branch for the Tree. In the case 
		//of three, create a new internal node, and create 
		//three branches to this node. One of these branches
		//will be chosen as the Root Branch, with the root 
		//point = 1.0 (all the way at the new internal node).
		//This is basically, pseudo-rooting the tree.
		TreeBranch rootBranch = null;
		if(topLevelNodeCount == 2) {
			//Create root branch connecting the two flanking nodes
			//root branch length will be the sum of the two parts
			double rootBranchLength = rootLengthA + rootLengthB;
			//root point will be the proportion of the branch length
			//assigned to rootFlankA
			double rootPoint = rootLengthA / rootBranchLength;

			//root branch support will mean the mean of the two parts
			double rootBranchSupport 
			            = (rootSupportA + rootSupportB) / 2;
			rootBranch = new TreeBranch(rootFlankA, rootFlankB,
					rootBranchLength, rootBranchSupport);
			rootBranch.setAsRoot(rootPoint);
		} else if(topLevelNodeCount == 3) {
			TreeNode iNode = new TreeNode();
			TreeBranch rootBranch1 
			           = new TreeBranch(iNode, topLevelNodes[0], 
					topLevelLengths[0], topLevelSupports[0]);
			TreeBranch rootBranch2 
			           = new TreeBranch(iNode, topLevelNodes[1],
			         topLevelLengths[1], topLevelSupports[1]);
			rootBranch = new TreeBranch(topLevelNodes[2], iNode,
					topLevelLengths[2], topLevelSupports[2]);
			rootBranch.setAsRoot(1.0);
		} else {
			System.out.println("TreeParser.parseTreeString() tree has "
					+ topLevelNodeCount + " top level nodes. " +
							"Only 2 and 3 are valid.");
			for(int i = 0; i < topLevelNodeCount; i++) {
				System.out.println("topLevelLengths " + i + ": " 
						+ topLevelLengths[i]);
				System.out.println("topLevelSupports " + i + ": " +
						topLevelSupports[i]);
			}
		}
		 r = new Tree(rootBranch);

		return r;
	}
}
