package edu.vt.vbi.ci.pepr.tree;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.apache.log4j.Logger;

import edu.vt.vbi.ci.util.ExtendedBitSet;

/**
 * AdvancedTree may contain one or more BasicTrees. It supports more 
 * advanced operations than BasicTree. 
 * 
 * @author enordber
 *
 */
public class AdvancedTree {

	public static final int EXPANDED = 1;
	public static final int COLLAPSED = 2;

	private BasicTree[] treeStack;
	private TreeState[] treeStateStack;
	private BasicTree currentTree;
	
	//index in treeStateStack for currentTree
	private int currentTreePointer; 
	
	private double[] nodeDistancesFromRoot;
	private int[] preorderTraversalSequence;
	private int[] nodeDescendantLeafCounts;
	private Logger logger = Logger.getLogger(getClass());

	/*
	 * current options for nodeState are EXPANDED and COLLAPSED
	 */
	private int[] nodeStates;

	/*
	 * Number of times node has been refined. Works as true/false with 1 and 0
	 */
	private int[] nodeRefined;
	private float[][] nodeCoordinates;
	private float[][] branchXCoordinates;
	private float[] tipDistancesFromPrevious;
	private int[] leaves;
	private int[] tips;
	private int[] branchSupports;
	private HashMap<String,String>[] nodeMetadata;

	/*
	 * treeDisplay mode is either TreeState.CLADOGRAM or TreeState.PHYLOGRAM
	 */
	private int treeDisplayMode = TreeState.PHYLOGRAM;

	private boolean usingLogBranchLengths = false;

	private ChangeListener[] changeListeners = new ChangeListener[0];
	private ChangeEvent changeEvent = new ChangeEvent(this);

	public AdvancedTree(String treeString) {
		this(new BasicTree(treeString));

		boolean looksLikeNewick = false;
		boolean looksLikeJSON = false;
		if(treeString.startsWith("(")) {
			looksLikeNewick = true;
		} else if(treeString.startsWith("{")) {
			looksLikeJSON = true;
		}

		if(looksLikeJSON) {
			//add additional info, such as refined nodes and collapsed nodes
		}

	}

	public AdvancedTree(BasicTree startingTree) {
		treeStack = new BasicTree[1];
		treeStack[0] = startingTree;
		currentTreePointer = 0;
		currentTree = startingTree;
		calcPreorderTraversalSequence();
		calcNodeDesendantLeaves();
	}

	/**
	 * AdvancedTree keeps a history of previous tree states. When the tree,
	 * by re-rooting, or removing nodes, etc., a copy of the old tree is kept.
	 * For large trees that are being manipulated a lot, this can use a big
	 * chunk of memory. This clearHistory() method is provided to allow
	 * the user to avoid keeping all of the old trees and using all of that 
	 * memory.
	 * This makes undo() and redo() mostly useless.
	 */
	public void clearHistory() {
		treeStack = new BasicTree[]{currentTree};
		currentTreePointer = 0;
	}
	public double getMaxDistanceFromRoot() {
		double r = 0;
		if(nodeDistancesFromRoot == null) {
			calcDistancesFromRoot();
		}
		for(int i = 0; i < nodeDistancesFromRoot.length; i++) {
			r = Math.max(r, nodeDistancesFromRoot[i]);
		}
		return r;
	}

	private double[] getNodeDistancesFromRoot() {
		if(nodeDistancesFromRoot == null) {
			calcDistancesFromRoot();
		}
		return nodeDistancesFromRoot;
	}

	/**
	 * Does preorder traversal of the tree, adding up the distances
	 * to each node from the root.
	 */
	private void calcDistancesFromRoot() {
		int[] nodeParentPointers = currentTree.getNodeParentPointers();
		double[] branchLengths = currentTree.getBranchLengths();
		if(treeDisplayMode == TreeState.CLADOGRAM) {
			branchLengths = currentTree.getCladogramBranchLengths();
		}

		nodeDistancesFromRoot = new double[branchLengths.length];
		nodeDistancesFromRoot[preorderTraversalSequence[0]] = 0;
		for(int i = 1; i < preorderTraversalSequence.length; i++) {
			int nodeIndex = preorderTraversalSequence[i];
			int parentIndex = nodeParentPointers[nodeIndex];
			//total distance from root equals distance from root to parent
			//plus distance from parent to this node
			nodeDistancesFromRoot[nodeIndex] = 
				nodeDistancesFromRoot[parentIndex] + branchLengths[nodeIndex];
		}
	}

	/**
	 * Does postorder traversal of the tree, adding up the number of
	 * descendants (leaves) for each node.
	 */
	private void calcNodeDesendantLeaves() {
		int[][] nodeChildPointers = currentTree.getNodeChildPointers();
		nodeDescendantLeafCounts = new int[nodeChildPointers.length];

		//for postorder traversal, go through preorderTraversalSequence backwards
		for(int i = preorderTraversalSequence.length-1; i >= 0; i--) {
			int node = preorderTraversalSequence[i];
			int[] children = nodeChildPointers[node];
			if(children.length == 0) {
				//if this node is a leaf, count it as having 1 descendant (itself).
				//yes, it's weird, but it works out right.
				nodeDescendantLeafCounts[node] = 1;
			}
			for(int j = 0; j < children.length; j++) {
				nodeDescendantLeafCounts[node] += 
					nodeDescendantLeafCounts[children[j]];
			}
		}
	}

	/**
	 * Determines the order to visit nodes for preorder traversal.
	 * This is done once here, with the result stored in the 
	 * preorderTraversalSequence array, to speed up future traversals. 
	 */
	private void calcPreorderTraversalSequence() {
		int[][] nodeChildPointers = currentTree.getNodeChildPointers();

		int[] nodeStates = getNodeStates();
		//only traverse EXPANDED nodes (actually, nodes without a COLLAPSED
		//ancestor node. The highest=level COLLAPSED node is traversed,
		//and represents the entire collapsed clade)
		preorderTraversalSequence = new int[nodeChildPointers.length];
		int topIndex = currentTree.getTopLevelNode();

		int traversalIndex = 0;

		//element 0 in the stack has the pointer to the current stack top
		int[] nodeStack = new int[preorderTraversalSequence.length];
		nodeStack[0] = 1;
		nodeStack[1] = topIndex;
		int currentNode;
		while(nodeStack[0] > 0) {
			currentNode = nodeStack[nodeStack[0]--];
			preorderTraversalSequence[traversalIndex++] = currentNode;
			if(nodeStates[currentNode] != COLLAPSED) {
				//don't add descendants of COLLAPSED nodes to the stack
				for(int i = nodeChildPointers[currentNode].length-1; i >= 0; i--) {
					nodeStack[++nodeStack[0]] = nodeChildPointers[currentNode][i];
				}
			}
		}

		//if traversalIndex != preorderTraversalSequence.length, it is because 
		//some nodes were skipped due to being COLLAPSED. Shrink 
		//preorderTraversalSequence
		if(traversalIndex < preorderTraversalSequence.length) {
			int[] trimmed = new int[traversalIndex];
			System.arraycopy(preorderTraversalSequence, 0, trimmed,
					0, trimmed.length);
			preorderTraversalSequence = trimmed;
		}

	}

	/**
	 * Sort child nodes so that the one with the most descendant leaves 
	 * is first.
	 */
	public void ladderizeUp() {
		BasicTree ladderizeTree = (BasicTree) currentTree.clone();
		int[][] nodeChildPointers = ladderizeTree.getNodeChildPointers();
		Comparator ladderizeUpComp = 
			new LadderizeSortComparator(ladderizeTree.getNodeStrings(), true);
		for(int i = 0; i < nodeChildPointers.length; i++) {
			if(nodeChildPointers[i].length > 1) {
				Integer[] sortableChildPointers = 
					new Integer[nodeChildPointers[i].length];
				for(int j = 0; j < sortableChildPointers.length; j++) {
					sortableChildPointers[j] = 
						new Integer(nodeChildPointers[i][j]);
				}

				Arrays.sort(sortableChildPointers, ladderizeUpComp);

				for(int j = 0; j < sortableChildPointers.length; j++) {
					nodeChildPointers[i][j] = sortableChildPointers[j].intValue();
				}
			}
		}

		setCurrentTree(ladderizeTree);
	}

	public void unroot() {
		BasicTree unrootTree = (BasicTree)currentTree.clone();
		unrootTree.unroot();
		setCurrentTree(unrootTree);
	}

	public String[] getLeafLabels() {
		String[] r = null;
		String[] leafLabels = currentTree.getLeaves(); 
		leaves = getLeaves();

		r = new String[leaves.length];
		for(int i = 0; i < r.length; i++) {
			r[i] = leafLabels[leaves[i]];
		}
		return r;
	}

	public String[] getTipLabels() {
		String[] r = null;
		String[] tipLabels = currentTree.getLeaves(); 
		int[] tips = getTips();
		int[] nodeStates = getNodeStates();

		r = new String[tips.length];
		for(int i = 0; i < r.length; i++) {
			if(nodeStates[tips[i]] == COLLAPSED) {
				r[i] = "Collapsed Node";
			} else {
				r[i] = tipLabels[tips[i]];
			}
		}
		return r;
	}

	public int[] getLeaves() {
		if(leaves == null) {
			String[] leafLabels = currentTree.getLeaves(); 
			leaves = new int[leafLabels.length];
			int leafIndex = 0;
			for(int i = 0; i < preorderTraversalSequence.length 
			&& leafIndex < leaves.length; i++) {
				int node = preorderTraversalSequence[i];
				if(nodeDescendantLeafCounts[node] == 1 && node < leafLabels.length) {
					//this is a leaf
					leaves[leafIndex++] = node; 
				}
			}

		}
		return leaves;
	}

	/**
	 * Returns an array of indices for the tree tips. Tips are either
	 * leaves or collapsed nodes. Leaves that are descendants of
	 * collapsed nodes are not included.
	 * @return
	 */
	public int[] getTips() {
		if(tips == null) {
			int[] nodeStates = getNodeStates();
			int[] nodeParentPointers = currentTree.getNodeParentPointers();
			String[] leafLabels = currentTree.getLeaves(); 
			tips = new int[leafLabels.length];
			int tipIndex = 0;
			for(int i = 0; i < preorderTraversalSequence.length 
			&& tipIndex < tips.length; i++) {
				int node = preorderTraversalSequence[i];
				if(nodeDescendantLeafCounts[node] == 1 
						|| nodeStates[node] == COLLAPSED) {
					//see if any ancestor nodes are COLLAPSED, if not, then
					//this node is a tip
					int parent = nodeParentPointers[node];
					boolean collapsedParentFound = false;
					while(parent >=0 && !collapsedParentFound) {
						collapsedParentFound = nodeStates[parent] == COLLAPSED;
						parent = nodeParentPointers[parent];
					}
					if(!collapsedParentFound) {
						tips[tipIndex++] = node; 
					}
				}
			}

			//if tipIndex is < tips.length, then there are fewer tips than 
			//leaves (becuase of collapsing). Shrink tips accordingly
			if(tipIndex < tips.length) {
				int[] shrinkTips = new int[tipIndex];
				System.arraycopy(tips, 0, shrinkTips, 0, tipIndex);
				tips = shrinkTips;
			}
		}
		return tips;
	}

	public double[] getLeafDistancesFromRoot() {
		double[] r =null;
		int[] leaves = getLeaves();
		r = new double[leaves.length];
		for(int i = 0; i < r.length; i++) {
			r[i] = nodeDistancesFromRoot[leaves[i]];
		}
		return r;
	}

	public double[] getTipDistancesFromRoot() {
		double[] r =null;
		int[] tips = getTips();
		r = new double[tips.length];
		for(int i = 0; i < r.length; i++) {
			r[i] = nodeDistancesFromRoot[tips[i]];
		}
		return r;
	}


	/**
	 * Returns the distance between the two specified nodes.
	 * 
	 * @param nodeA
	 * @param nodeB
	 * @return
	 */
	public double getDistanceBetweenNodes(int nodeA, int nodeB) {
		double r = Double.NaN;
		//make a clone of the tree
		BasicTree cloneTree = (BasicTree)currentTree.clone();

		//root the tree at nodeA
		int nodeAParent = cloneTree.getNodeParentPointers()[nodeA];
		cloneTree.rootBetweenNodes(nodeA, nodeAParent, 0);

		AdvancedTree advancedClone = new AdvancedTree(cloneTree);		
		//get the distance to nodeB from the root
		r = advancedClone.getNodeDistancesFromRoot()[nodeB];
		return r;
	}

	/**
	 * Returns a list of y-coordinate distances. These are the distance from
	 * the leaf above to the current leaf. These will typically all be 1.0,
	 * but will be smaller durring a collapsing or expanding operation.
	 * @return
	 */
	private float[] getLeafDistanceFromPrevious() {
		float[] r = null;
		int[] leaves = getLeaves();
		r = new float[leaves.length];
		for(int i = 0; i < r.length; i++) {
			r[i] = 1.0f;
		}

		return r;
	}


	/**
	 * Returns a list of y-coordinate distances. These are the distance from
	 * the tip above to the current tip. These will typically all be 1.0,
	 * but will be smaller during a collapsing or expanding operation.
	 * @return
	 */
	private float[] getTipDistanceFromPrevious() {
		if(tipDistancesFromPrevious == null) {
			int[] tips = getTips();
			tipDistancesFromPrevious = new float[tips.length];
			for(int i = 0; i < tipDistancesFromPrevious.length; i++) {
				tipDistancesFromPrevious[i] = 1.0f;
			}
		}
		return tipDistancesFromPrevious;
	}

	/**
	 * Returns a set of coordinates for the y positions for nodes in the tree.
	 * This is a 2-dimension array, with the index of the first dimension
	 * being the node index. Each 2nd dimension array has 2 values. 
	 * [0] = y coordinate of center of node
	 * [1] = half the height of the vertical node line 
	 * The units for these are 1 per leaf in the tree, and must be scaled
	 * to be visually useful.
	 * 
	 * @return
	 */
	public float[][] getNodeCoordinates() {
		if(nodeCoordinates == null) {
			nodeCoordinates = new float[currentTree.getNodeCount()][2];
			int[][] nodeChildPointers = currentTree.getNodeChildPointers();
			//use postorder traversal to determine y coordinates. The scale of
			//y coordinates is 1 per leaf. The top leaf is 0 and the bottom leaf
			//is leaves.length - 1
			int[] tips = getTips();
			float[] tipYDistances = getTipDistanceFromPrevious();
			nodeCoordinates[tips[0]][0] = 0;
			nodeCoordinates[tips[0]][1] = 0;

			for(int i = 1; i < tips.length; i++) {
				nodeCoordinates[tips[i]][0] = 
					nodeCoordinates[tips[i-1]][0] + tipYDistances[i];
				nodeCoordinates[tips[i]][1] = 0;
			}

			int[] nodeStates = getNodeStates();
			for(int i = preorderTraversalSequence.length-1; i >= 0; i--) {
				int node = preorderTraversalSequence[i];
				int[] children = nodeChildPointers[node];
				if(children.length > 0 && nodeStates[node] != COLLAPSED) {
					//this is not a leaf. set the y coordinate to half way between
					//the smallest and the largest y coords for the children.
					//Children should already be sorted from smallest to largest
					float maxY = 0; 
					float minY = Float.MAX_VALUE;
					for(int j = 0; j < children.length; j++) {
						float childY = nodeCoordinates[children[j]][0];
						maxY = Math.max(childY, maxY);
						minY = Math.min(childY, minY);
					}
					float halfNode = (maxY-minY)/2; 
					nodeCoordinates[node][0] = minY + halfNode;
					nodeCoordinates[node][1] = halfNode;

				}
			}
		}
		return nodeCoordinates;
	}

	public int getNodeCount() {
		return currentTree.getNodeCount();
	}

	public int getParentNode(int node) {
		int r = -1;
		r = currentTree.getNodeParentPointers()[node];
		return r;
	}

	public int[] getBranchSupports() {
		if(branchSupports == null) {
			int defaultBranchSupport = 100;
			String[] supportStrings = currentTree.getBranchSupportStrings();
			branchSupports = new int[supportStrings.length];
			for(int i = 0; i < branchSupports.length; i++) {
				if(supportStrings[i] == null || supportStrings[i].length() ==0) {
					branchSupports[i] = defaultBranchSupport;
				} else {
					branchSupports[i] = Integer.parseInt(supportStrings[i]);
				}
			}
		}
		return branchSupports;
	}

	public float[][] getBranchXCoordinates() {
		if(branchXCoordinates == null) {
			double[] nodeDistancesFromRoot = getNodeDistancesFromRoot();
			branchXCoordinates = new float[currentTree.getNodeCount()][2];
			double[] branchLengths = currentTree.getBranchLengths();
			if(treeDisplayMode == TreeState.CLADOGRAM) {
				branchLengths = currentTree.getCladogramBranchLengths();
			}
			for(int i = 1; i < preorderTraversalSequence.length; i++) {
				int node = preorderTraversalSequence[i];
				branchXCoordinates[node][1] = (float) nodeDistancesFromRoot[node];
				branchXCoordinates[node][0] = 
					(float) (branchXCoordinates[node][1] - branchLengths[node]);
			}
		}
		return branchXCoordinates;
	}

	/**
	 * collapses the specified node. The clade rooted at this
	 * node will appear as a single tip.
	 * @param node
	 */
	public void collapseNode(int node) {
		System.out.println("AdvancedTree.collapseNode() " + node);

		//first animate the collapse by decreasing the spacing between leaves
		//that are descendants of the collapsing node
		int[] descendantTips = getDescendantTips(node);

		//now actually collapse the node
		int[] nodeStates = getNodeStates();
		nodeStates[node] = COLLAPSED;

		TreeState treeState = new TreeState();
		treeState.setTree(currentTree);
		treeState.setNodeStates(nodeStates);
		treeState.setGeneratingAction(TreeState.COLLAPSE_NODE, new int[]{node});
		setNewTreeState(treeState);
	}

	private boolean isTip(int node) {
		boolean r = false;
		if(node >= 0) {
			r = nodeStates[node] == COLLAPSED;
			if(!r) {
				r = currentTree.getNodeChildPointers()[node].length == 0;
			}
		}
		return r;
	}

	private int[] getDescendantTips(int node) {
		int[] r = null;
		if(isTip(node)) {
			//this is the end of the recursion - we have reached a tip
			r = new int[]{node};
		} else {
			r = new int[0];
			int[] descendants;
			int[] children = currentTree.getNodeChildPointers()[node];
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

	public String[] getDescendantLeaves(int node) {
		String[] r = null;
		String[] labels = currentTree.getLeaves();
		int[] indices = getDescendantTips(node);
		ArrayList<String> leafList = new ArrayList<String>();
		for(int i = 0; i < indices.length; i++) {
			if(labels.length > indices[i]) {
				try{
					leafList.add(labels[indices[i]]);
				} catch(ArrayIndexOutOfBoundsException e) {
					logger.error("array index out of bounds. trying index "
							+ indices[i] + " in labels array with length " + labels.length);
				}
			}
		}
		r = new String[leafList.size()];
		leafList.toArray(r);
		return r;
	}

	private int[] getNodeStates() {
		if(nodeStates == null) {
			//begin with all nodes EXPANDED
			nodeStates = new int[currentTree.getNodeCount()];
			Arrays.fill(nodeStates, EXPANDED);
		} else if(currentTree.getNodeCount() != nodeStates.length) {
			int[] newStates = new int[currentTree.getNodeCount()];
			Arrays.fill(newStates, EXPANDED);
			int copyCount = Math.min(nodeStates.length, newStates.length);
			System.arraycopy(nodeStates, 0, newStates, 0, copyCount);
			nodeStates = newStates;
		}
		return nodeStates;
	}

	public int[] getNodeRefined() {
		if(nodeRefined == null) {
			nodeRefined = new int[currentTree.getNodeCount()];
		}
		return nodeRefined;
	}

	public int getNodeRefined(int index) {
		if(nodeRefined == null) {
			nodeRefined = new int[currentTree.getNodeCount()];
		}
		return nodeRefined[index];
	}

	public void setNodeRefined(int index, int value) {
		if(nodeRefined == null) {
			nodeRefined = new int[currentTree.getNodeCount()];
		}
		nodeRefined[index] = value;
	}

	public void setNodeRefined(int[] refined) {
		nodeRefined = new int[refined.length];
		System.arraycopy(refined, 0, nodeRefined, 0, refined.length);
	}

	/**
	 * Provides a set of outgroup taxa. Causes the tree to attempt
	 * to root between the outgroup set and the ingroup set 
	 * (the remainder of the taxa). This may not be possible, so the root
	 * will be selected to maximize the separation between the two sets.
	 */
	public void setOutGroup(String[] outGroupLabels) {
		String[] allTaxa = getLeafLabels();
		String[] sortedTaxa = TreeUtils.compressTaxonNamesForComparison(allTaxa);
		Arrays.sort(sortedTaxa);
		
		//convert outGroupLabels to compressed names for comparison.
		String[] compressedOutgroup = TreeUtils.compressTaxonNamesForComparison(outGroupLabels);
		outGroupLabels = compressedOutgroup;

		BitSet outgroupSet = new BitSet(sortedTaxa.length);
		for(int i = 0; i < outGroupLabels.length; i++) {
			int index = Arrays.binarySearch(sortedTaxa, outGroupLabels[i]);
			if(index >= 0) {
				outgroupSet.set(index);
			}
		}

		//re-root at some branch leading to an ingroup member
		String[] nodeLabels = currentTree.getNodeStrings();
		//compress nodLabels for comparison 
		nodeLabels = TreeUtils.compressTaxonNamesForComparison(nodeLabels);
		int[] leaves = getLeaves();
		int rootBefore = -1;
		for(int i = 0; i < leaves.length && rootBefore < 0; i++) {
			//get the index for this leaf
			int index = Arrays.binarySearch(sortedTaxa, nodeLabels[leaves[i]]);
			//if this index is not in outgroupSet, then use this 
			//as the leaf for rooting
			if(!outgroupSet.get(index)) {
				rootBefore = leaves[i];
			}
		}
		int rootAfter = currentTree.getNodeParentPointers()[rootBefore];
		//create the new tree to be the rooted version. This one will not
		//be stored long term in the treeStack, because it is only being used to
		//determine where to root for the outgroup. It will be temporarily 
		//stored in the stack, to allow the use of AdvancedTree methods 
		BasicTree findOutgroupTree = (BasicTree) currentTree.clone();
		findOutgroupTree.rootBetweenNodes(rootBefore, rootAfter, 0.5);
		setCurrentTree(findOutgroupTree);
		calcPreorderTraversalSequence();

		int deepestNodeWithMax = getMostEnrichedNode(outGroupLabels);
		
		int parentNode = currentTree.getNodeParentPointers()[deepestNodeWithMax];

		BasicTree outgroupRootedTree = (BasicTree) currentTree.clone();
		outgroupRootedTree.rootBetweenNodes(deepestNodeWithMax, parentNode, 0.5);

		//move currentTreePointer back one, so temporary tree created in this
		//method won't be retained in the treeStack. 
		currentTreePointer--;
		TreeState treeState = new TreeState();
		treeState.setTree(outgroupRootedTree);
		treeState.setNodeStates(getNodeStates());
		treeState.setGeneratingAction(TreeState.REROOT, 
				new int[]{deepestNodeWithMax, parentNode});
		setNewTreeState(treeState);
	}

	private int getMostEnrichedNode(String[] labels) {
		String[] allTaxa = getLeafLabels();
		String[] sortedTaxa = TreeUtils.compressTaxonNamesForComparison(allTaxa);
		Arrays.sort(sortedTaxa);
		
		BitSet outgroupSet = new BitSet(sortedTaxa.length);
		for(int i = 0; i < labels.length; i++) {
			int index = Arrays.binarySearch(sortedTaxa, labels[i]);
			if(index >= 0) {
				outgroupSet.set(index);
			}
		}

		//determine the taxa set for each node, using BitSets
		int[][] nodeChildPointers = currentTree.getNodeChildPointers();
		BitSet[] nodeTaxaSets = new BitSet[nodeChildPointers.length];

		//initialize all of the bit sets
		for(int i = 0; i < nodeTaxaSets.length; i++) {
			nodeTaxaSets[i] = new BitSet(sortedTaxa.length);
		}

		String[] nodeLabels = TreeUtils.compressTaxonNamesForComparison(currentTree.getNodeStrings());
		//first create the leaf sets, which will each have one bit set
		for(int i = 0; i < leaves.length; i++) {
			int leafIndex = leaves[i];
			int taxonIndex = Arrays.binarySearch(sortedTaxa, nodeLabels[leafIndex]);
			nodeTaxaSets[leafIndex].set(taxonIndex);
		}

		//now build the taxon sets for the internal nodes by ORing the child 
		//node sets. Do this with a postorder traversal of the tree, so that 
		//every node visited will already have it's child nodes visited prior
		for(int i = preorderTraversalSequence.length-1; i >= 0; i--){
			int thisNode = preorderTraversalSequence[i];
			if(thisNode < nodeTaxaSets.length) {
				int[] childNodes = nodeChildPointers[thisNode];
				for(int j = 0; j < childNodes.length; j++) {
					int childNode = childNodes[j];
					BitSet nodeTaxaSet = nodeTaxaSets[childNode];
					nodeTaxaSets[thisNode].or(nodeTaxaSet);
				}
			} else {
				logger.warn("AdvancedTree.setOutgroup() node is " + 
						thisNode + " but nodeTaxaSets.length is only " +
						nodeTaxaSets.length);
			}
		}

		//for each node, determine how many of the taxa are in the outgroup set
		int[] outMinusInCounts = new int[nodeTaxaSets.length];
		BitSet testerSet = new BitSet(sortedTaxa.length);
		for(int i = 0; i < outMinusInCounts.length; i++) {
			testerSet.clear();
			testerSet.or(outgroupSet);
			testerSet.and(nodeTaxaSets[i]);
			int outgroupCount = testerSet.cardinality();
			int ingroupCount = nodeDescendantLeafCounts[i] - outgroupCount;
			outMinusInCounts[i] = outgroupCount - ingroupCount;
		}

		//traverse tree preorder, keeping the last node with the maximum out-in
		int maxValue = outMinusInCounts[preorderTraversalSequence[0]];
		int deepestNodeWithMax = preorderTraversalSequence[0];
		for(int i =1; i < preorderTraversalSequence.length; i++) {
			int node = preorderTraversalSequence[i];
			if(outMinusInCounts[node] >= maxValue) {
				maxValue = outMinusInCounts[node];
				deepestNodeWithMax = node;
			}
		}

		return deepestNodeWithMax;
	}

	private void setCurrentTree(BasicTree tree) {
		//create a new TreeState with this tree and the current nodeStates
		TreeState treeState = new TreeState();
		treeState.setTree(tree);
		treeState.setNodeStates(getNodeStates());
		setNewTreeState(treeState);
	}

	public void makePhylogram() {
		BasicTree phylogramTree = (BasicTree)currentTree.clone();
		TreeState newState = new TreeState();
		newState.setTree(phylogramTree);
		newState.setGeneratingAction(TreeState.PHYLOGRAM);
		setNewTreeState(newState);
	}

	public void makeCladogram() {
		BasicTree cladogramTree = (BasicTree)currentTree.clone();
		TreeState newState = new TreeState();
		newState.setTree(cladogramTree);
		newState.setGeneratingAction(TreeState.CLADOGRAM);

		if(currentTree.getCladogramBranchLengths() != null) {
			//nothing special needs to be done
		} else {
			//if currentTree does not have cladogram branch lengths, 
			//then create them.

			//traverse the tree, setting branch lengths for cladogram
			//display. Each node has a depth from the tip, which is
			//the maximum number of nodes along the path to the tips.
			//This represents the number of branch length units to be 
			//split up along each path from the node to the tips.

			//do postorder traversal to get the max distance to a tip for each node
			int[][] nodeChildPointers = currentTree.getNodeChildPointers();
			int[] maxStepsToTip = new int[preorderTraversalSequence.length];
			for(int i = preorderTraversalSequence.length-1; i >= 0; i--) {
				int node = preorderTraversalSequence[i];
				maxStepsToTip[node] = 0;
				int[] childNodes = nodeChildPointers[node];
				for(int j = 0; j < childNodes.length; j++) {
					maxStepsToTip[node] = 
						Math.max(maxStepsToTip[node], maxStepsToTip[childNodes[j]]);
				}
				//add a step for this node
				maxStepsToTip[node]++;
			}

			//do a preorder traversal dividing up the amount of distance available
			//betwen branches
			int maxSteps = maxStepsToTip[preorderTraversalSequence[0]];

			int[] remainingLength = new int[preorderTraversalSequence.length];
			double[] cladogramBranchLengths = new double[preorderTraversalSequence.length];
			remainingLength[preorderTraversalSequence[0]] = maxSteps;
			int[] nodeParentPointers = currentTree.getNodeParentPointers();
			int root = currentTree.getRootIndex();
			for(int i = 1; i < preorderTraversalSequence.length; i++) {
				int node = preorderTraversalSequence[i];
				int[] childNodes = nodeChildPointers[node];
				cladogramBranchLengths[node] = 
					remainingLength[nodeParentPointers[node]] - maxStepsToTip[node];
				remainingLength[node] = 
					(int) (remainingLength[nodeParentPointers[node]] - cladogramBranchLengths[node]);
			}

			//set cladogram branch lengths for cladogramTree
			cladogramTree.setCladogramBranchLengths(cladogramBranchLengths);
		}
		//make cladogram tree the current tree
		setNewTreeState(newState);
	}

	private void setNewTreeState(TreeState newState) {
		if(treeStateStack == null) {
			currentTreePointer = -1;
			treeStateStack = new TreeState[1];
		}

		currentTreePointer++;
		if(treeStateStack.length <= currentTreePointer) {
			//increase stack size to make room for the new state Object
			TreeState[] newStack = new TreeState[treeStateStack.length+1];
			System.arraycopy(treeStateStack, 0, newStack, 0, treeStateStack.length);
			treeStateStack = newStack;
		}

		//put the new state in the stack
		treeStateStack[currentTreePointer] = newState;

		//point currentTree at the tree in the new state
		currentTree = newState.getTree();
		if(newState.getGeneratingAction() == TreeState.CLADOGRAM) {
			treeDisplayMode = TreeState.CLADOGRAM;
		} else if(newState.generatingAction == TreeState.PHYLOGRAM) {
			treeDisplayMode = TreeState.PHYLOGRAM;
		}
		refresh();
		fireChangeEvent();
	}

	/**
	 * Several arrays are calculated and stored to speed up various operations.
	 * Some actions on the tree invalidate these and they need to be
	 * recalculated. This method will force all of these things to be
	 * recalculated, either immediately or upon the next time they are needed.
	 */
	private void refresh() {
		calcPreorderTraversalSequence();
		calcNodeDesendantLeaves();
		branchXCoordinates = null;
		nodeCoordinates = null;
		nodeDistancesFromRoot = null;
		leaves = null;
		tips = null;
		tipDistancesFromPrevious = null;
		branchSupports = null;
		nodeRefined = null;
	}


	public String getTreeString(boolean includeLengths, boolean includeSupports) {
		String r = currentTree.getTreeString(includeLengths, includeSupports);
		return r;
	}

	public void setUseLogBranchLengths(boolean useLogBranchLengths) {
		if(useLogBranchLengths != usingLogBranchLengths) {
			//changing status of using log branch lengths
			double logAdd = 1.01;
			double base = 10;
			double logBase = 1; //Math.log(base);

			if(useLogBranchLengths) {
				//clone tree
				BasicTree logTree = (BasicTree) currentTree.clone();

				//get branch lengths of tree and do log
				double[] branchLengths = logTree.getBranchLengths();
				for(int i = 0; i < branchLengths.length; i++) {
					branchLengths[i] = Math.log(branchLengths[i]+logAdd)/logBase;
				}

				//set branch lengths to log lengths
				logTree.setBranchLengths(branchLengths);

				//make this the current tree
				TreeState newState = new TreeState();
				newState.setTree(logTree);
				newState.setGeneratingAction(TreeState.LOG_BRANCH);
				setNewTreeState(newState);
			} else {
				//clone tree
				BasicTree expTree = (BasicTree) currentTree.clone();

				//do exponent of branch lengths
				double[] branchLengths = expTree.getBranchLengths();
				for(int i = 0; i < branchLengths.length; i++) {
					branchLengths[i] = Math.pow(branchLengths[i], base) - logAdd;
				}

				//set branch lengths to exponent lengths
				expTree.setBranchLengths(branchLengths);

				//make this the current tree
				TreeState newState = new TreeState();
				newState.setTree(expTree);
				newState.setGeneratingAction(TreeState.EXP_BRANCH);
			}
		}
	}

	private void fireChangeEvent() {
		for(int i = 0; i < changeListeners.length; i++) {
			changeListeners[i].stateChanged(changeEvent);
		}
	}

	/**
	 * Adds a listener to be notified of certain changes in this tree.
	 * Changes triggering an event include any change of the current tree
	 * to a new tree. This happens upon rerooting, ladderizing, node 
	 * expand/collapse.
	 * @param cl
	 */
	public void addChangeListener(ChangeListener cl) {
		ChangeListener[] newChangeListeners = new ChangeListener[changeListeners.length+1];
		System.arraycopy(changeListeners, 0, newChangeListeners, 0, changeListeners.length);
		newChangeListeners[changeListeners.length] = cl;
		changeListeners = newChangeListeners;
	}

	public void removeChangeListener(ChangeListener cl) {
		//see if cl is in the changeListsners list
		//if it is in the list, then remove it
		int clIndex = -1;
		for(int i = 0; i < changeListeners.length; i++) {
			if(changeListeners[i] == cl) {
				clIndex = i;
				break;
			}
		}

		if(clIndex >=0) {
			ChangeListener[] newChangeListeners = new ChangeListener[changeListeners.length-1];
			int index = 0;
			for(int i = 0; i < changeListeners.length; i++) {
				if(i != clIndex) {
					newChangeListeners[index++] = changeListeners[i];
				}
			}
			changeListeners = newChangeListeners;
		}
	}

	/**
	 * If the currentTree is not the first tree in the treeStack, then 
	 * undo() will set the currentTree to the previous tree in the stack.
	 * 
	 * @return true if an undo operation was performed. False if there
	 *         are no previous trees in the stack 
	 */
	public boolean undo() {
		boolean r = false;
		if(currentTreePointer > 0) {
			System.out.println("AdvancedTree.undo() changing current tree from "
					+ currentTreePointer + " to " + (currentTreePointer-1));
			currentTreePointer--;
			currentTree = treeStateStack[currentTreePointer].getTree();
			int[] redoStates = treeStateStack[currentTreePointer].getNodeStates();
			nodeStates = new int[redoStates.length];
			System.arraycopy(redoStates, 0, nodeStates, 0, redoStates.length);

			refresh();
			r = true;
			fireChangeEvent();
		}
		return r;
	}

	/**
	 * If the current tree is not the final tree in the treeStack, then 
	 * redo() will set the currentTree to the next tree in the stack.
	 * 
	 * @return true if a redo operation was performed. False if there
	 *         are no subsequent trees in the stack
	 */
	public boolean redo() {
		boolean r = false;
		if(treeStateStack.length > currentTreePointer+1) {
			currentTreePointer++;
			currentTree = treeStateStack[currentTreePointer].getTree();
			nodeStates = treeStateStack[currentTreePointer].getNodeStates();
			refresh();
			r = true;
			fireChangeEvent();
		}
		return r;
	}

	public void removeTaxon(String taxon) {
		BasicTree cloneTree = (BasicTree) currentTree.clone();
		cloneTree.removeTaxon(taxon);
		TreeState newState = new TreeState();
		newState.setTree(cloneTree);
		setNewTreeState(newState);
	}


	public void setBranchSupportValues(int[] branchSupportCounts) {
		String[] bsStrings = new String[branchSupportCounts.length];
		for(int i = 0; i < bsStrings.length; i++) {
			bsStrings[i] = "" + branchSupportCounts[i];
		}
		currentTree.setBranchSupportStrings(bsStrings);
		refresh();
	}

	public int[] getMeanDescendantSupportValues() {
		int[] r = new int[preorderTraversalSequence.length];

		int[] descendantSupportSums = new int[r.length];
		int[] nodeDescendantCounts = new int[r.length];
		int[] branchSupports = getBranchSupports();

		int[][] nodeChildPointers = currentTree.getNodeChildPointers();
		//do post-order traversal
		for(int i = preorderTraversalSequence.length-1; i >= 0; i--) {
			int node = preorderTraversalSequence[i];
			//get child nodes
			int[] childNodes = nodeChildPointers[node];
			for(int j = 0; j < childNodes.length; j++) {
				int child = childNodes[j];
				//add support values for each child
				descendantSupportSums[node] += descendantSupportSums[child];
				descendantSupportSums[node] += branchSupports[child];

				//add number of child nodes to descendant counts
				nodeDescendantCounts[node] += nodeDescendantCounts[child];
				if(!isTip(node)) {
					nodeDescendantCounts[node]++;
				}
			}

			if(nodeDescendantCounts[node] == 0) {
				//this is a terminal node. consider the support to be 0.
				//the terminal nodes are ignored for this calculation
				r[node] = 0;
			} else {
				r[node] = 
					(int) Math.floor((double)descendantSupportSums[node] / nodeDescendantCounts[node]);
			}
		}

		return r;
	}

	/**
	 * Finds the node containing the given set of leaves. If there is
	 * no node in the tree containing exactly this set of leaves, then
	 * -1 is returned.
	 * 
	 * @param leaves
	 * @return
	 */
	public int getMatchingNode(String[] leaves) {
		int r = -1;

		HashSet leafSet = new HashSet();
		for(int i = 0; i < leaves.length; i++) {
			leafSet.add(leaves[i]);
		}

		for(int i = 0; r < 0 && i < preorderTraversalSequence.length; i++) {
			int node = preorderTraversalSequence[i];
			String[] nodeDescendants = getDescendantLeaves(node);
			if(nodeDescendants.length == leaves.length) {
				//correct number of leaves
				//check to see if each node descendant is one of the 
				//target leaves. If a single descendant is found that is not
				//one of the target leave, then this is not the matching node
				boolean allMatch = true;
				for(int j = 0; allMatch && j < nodeDescendants.length; j++) {
					if(!leafSet.contains(nodeDescendants[j])) {
						allMatch = false;
					}
				}
				if(allMatch) {
					//this is the matching node
					r = node;
				}
			}
		}
		return r;
	}

	public BasicTree getBasicTree() {
		return currentTree;
	}

	public int[] getPreorderTraversalSequence() {
		return preorderTraversalSequence;
	}

	/**
	 * Tries to replace a node in this tree with the given subtree. 
	 * The subtree should be rooted, with the largest side being the
	 * portion to be added to this tree. The tree is searched for the largest
	 * node for which every leaf is contained within the given subtree.
	 *  
	 * @param subtree
	 * @return
	 */
	public AdvancedTree replaceNode(BasicTree subtree) {
		AdvancedTree r= null;
		int subtreeRoot = subtree.getRootIndex();
		int[] rootChildren = subtree.getNodeChildPointers()[subtreeRoot];
		HashSet subtreeMembers = new HashSet();
		String[][] rootChildLeaves = new String[rootChildren.length][];
		AdvancedTree subtreeAdvanced = new AdvancedTree(subtree);
		for(int i = 0; i < rootChildLeaves.length; i++) {
			rootChildLeaves[i] = subtreeAdvanced.getDescendantLeaves(rootChildren[i]);
		}
		String[] subtreeIngroup = null;
		if(rootChildLeaves[0].length > rootChildLeaves[1].length) {
			subtreeIngroup = rootChildLeaves[0];
		} else {
			subtreeIngroup = rootChildLeaves[1];
		}

		for(int i = 0; i < subtreeIngroup.length; i++) {
			subtreeMembers.add(subtreeIngroup[i]);
		}


		//subtreeMembers now contains all ingroup members from the subtree to
		//be added. Find that largest node in this tree for which all 
		//descendants are found in subtreeMembers.
		int replaceNode = -1;
		for(int i = 0; replaceNode < 0 && i < preorderTraversalSequence.length; i++) {
			int node = preorderTraversalSequence[i];
			//get leaves for node
			String[] nodeLeaves = getDescendantLeaves(node);
			//if nodeLeaves has more members than subtreeMembers, then this
			//can't be the best node
			if(nodeLeaves.length <= subtreeMembers.size()) {
				//see if all nodeLeaves are in subtreeMembers
				boolean allFound = true;
				for(int j = 0; j < nodeLeaves.length; j++) {
					if(!subtreeMembers.contains(nodeLeaves[j])) {
						allFound = false;
					}
				}
				if(allFound) {
					replaceNode = node;
				}
			}
		}

		BasicTree thisClone = (BasicTree)currentTree.clone();
		BasicTree replacedBasic = thisClone.replaceSubtreeBelow(replaceNode, subtree);
		r = new AdvancedTree(replacedBasic);

		return r;
	}

	/**
	 * Returns a set of Bipartitions defined by this tree. The index
	 * values for the Bipartitions are based on the sorted list of
	 * leaf names, which is not necessarily the order returned by 
	 * AdvancedTree.getLeafLabels(). 
	 * 
	 * @return
	 */
	public Bipartition[] getBipartitions() {
		Bipartition[] r = null;
		ArrayList<Bipartition> bipartList = new ArrayList<Bipartition>();
		String[] leaves = getLeafLabels();

		//sort the leaves, to allow using binary search to more quickly
		//find the indices.
		Arrays.sort(leaves);
		int[] preorderTraversalSequence = getPreorderTraversalSequence();

		int nodeCount = getNodeCount();
		for(int i = 0; i < nodeCount; i++) {
			String[] descendants = 
				getDescendantLeaves(i);
			if(descendants.length > 1 && descendants.length < leaves.length -1) {
				ExtendedBitSet ebs = new ExtendedBitSet();
				//set the bit for each leaf that is present in the list of
				//descendants for this node. This represents one side
				//of the bipartition based on the branch leading to this node.
				for(int j = 0; j < descendants.length; j++) {
					int index = Arrays.binarySearch(leaves, descendants[j]);
					if(index > -1) {
						ebs.set(index);
					}
				}
				Bipartition bipart = new Bipartition(ebs, leaves.length);
				bipartList.add(bipart);
			}
		}
		r = new Bipartition[bipartList.size()];
		r = bipartList.toArray(r);
		return r;
	}

	/**
	 * Returns null for trivial Bipartitions.
	 * 
	 * @param node
	 * @return
	 */
	public Bipartition getBipartitionForNode(int node) {
		Bipartition r = null;
		String[] leaves = getLeafLabels();
		//sort the leaves, to allow using binary search to more quickly
		//find the indices.
		Arrays.sort(leaves);
		
		String[] descendants = 
			getDescendantLeaves(node);
		if(descendants.length > 1 && descendants.length < leaves.length -1) {
			ExtendedBitSet ebs = new ExtendedBitSet();
			//set the bit for each leaf that is present in the list of
			//descendants for this node. This represents one side
			//of the bipartition based on the branch leading to this node.
			for(int j = 0; j < descendants.length; j++) {
				int index = Arrays.binarySearch(leaves, descendants[j]);
				if(index > -1) {
					ebs.set(index);
				}
			}
			r = new Bipartition(ebs, leaves.length);
		}
		return r;
	}
	
	/**
	 * Returns a set of Bipartitions defined by this tree. The index
	 * values for the Bipartitions are based on the order in the provided 
	 * array of taxonNames. The list of taxon names does not have to be identical to 
	 * the taxon names in the tree, but there must be some overlap to get
	 * any meaningful Bipartitions returned.
	 * 
	 * @return
	 */
	public Bipartition[] getBipartitions(String[] taxonNames) {
		Bipartition[] r = null;
		
		//determine the set of taxa in taxonNames that are participating in this tree
		ExtendedBitSet participatingTaxa = new ExtendedBitSet();
		String[] leaves = getLeafLabels();
		Arrays.sort(leaves);
		for(int i = 0; i < taxonNames.length; i++) {
			int index = Arrays.binarySearch(leaves, taxonNames[i]);
			if(index >= 0) {
				participatingTaxa.set(i);
			}
		}
		
		ArrayList<Bipartition> bipartList = new ArrayList<Bipartition>();

		int[] preorderTraversalSequence = getPreorderTraversalSequence();

		for(int i = 0; i < preorderTraversalSequence.length; i++) {
			String[] descendants = 
				getDescendantLeaves(preorderTraversalSequence[i]);
			if(descendants.length > 1 && descendants.length < leaves.length -1) {
				ExtendedBitSet ebs = new ExtendedBitSet();
				//set the bit for each leaf that is present in the list of
				//descendants for this node. This represents one side
				//of the bipartition based on the branch leading to this node.

				Arrays.sort(descendants);
				for(int j = 0; j < taxonNames.length; j++) {
					int index = Arrays.binarySearch(descendants, taxonNames[j]);
					if(index >= 0) {
						ebs.set(j);
					}
				}
				Bipartition bipart = new Bipartition(ebs, taxonNames.length);
				bipart.setParticipatingTaxonSet(participatingTaxa);
				bipartList.add(bipart);
			}
		}
		r = new Bipartition[bipartList.size()];
		r = bipartList.toArray(r);
		return r;
	}


	private class LadderizeSortComparator implements Comparator {
		boolean up = true;
		String[] nodeLabels;

		public LadderizeSortComparator(String[] labels, boolean ladderizeUp) {
			this.up = ladderizeUp;
			this.nodeLabels = labels;
		}

		public int compare(Object objA, Object objB) {
			int r = -1;
			Integer nodeA = (Integer)objA;
			Integer nodeB = (Integer)objB;

			int nodeADesc = nodeDescendantLeafCounts[nodeA.intValue()];
			int nodeBDesc = nodeDescendantLeafCounts[nodeB.intValue()];
			if(nodeADesc == 1 && nodeBDesc == 1) {
				//these are both leaves. Sort alphabetically by node label
				//				r = nodeLabels[nodeB.intValue()].compareTo(nodeLabels[nodeA.intValue()]);
			} else {
				r = nodeBDesc - nodeADesc; 

				if(up) {
					//r = -r;
				}
			}
			return r;
		}

	}

	/**
	 * This stores the state information for the tree, and is used
	 * for undo/redo functionality.
	 * 
	 * @author enordber
	 *
	 */
	private class TreeState {
		/*
		 * Action constants
		 */
		private static final int EXPAND_NODE = 1;
		private static final int COLLAPSE_NODE = 2;
		private static final int UNROOT = 4;
		private static final int REROOT = 8;
		private static final int LOG_BRANCH = 16;
		private static final int EXP_BRANCH = 32;
		private static final int CLADOGRAM = 64;
		private static final int PHYLOGRAM = 128;

		private BasicTree tree;
		private int[] nodeStates;

		/*
		 * The action that resulted in this version of the tree being created
		 * Possible values are EXPAND_NODE, COLLAPSE_NODE, REROOT, UNROOT
		 */
		private int generatingAction;

		/*
		 * for collapse/expand operations, this will have one value: the node
		 * that was expanded or collapsed. 
		 * For rooting operations this will contain the two nodes on either 
		 * side the root.
		 * For unrooting operations this will be a zero-length array
		 */
		private int[] nodesInvolved;

		public TreeState() {

		}

		private BasicTree getTree() {
			return tree;
		}

		private void setTree(BasicTree tree) {
			this.tree = tree;
		}

		private int[] getNodeStates() {
			return nodeStates;
		}

		private void setNodeStates(int[] nodeStates) {
			this.nodeStates = new int[nodeStates.length];
			System.arraycopy(nodeStates, 0, this.nodeStates, 0, nodeStates.length);
		}

		private void setGeneratingAction(int action, int[] nodesInvolved) {
			generatingAction = action;
			this.nodesInvolved = nodesInvolved;
		}

		private int getGeneratingAction() {
			return generatingAction;
		}

		private void setGeneratingAction(int generatingAction) {
			this.generatingAction = generatingAction;
		}

		private int[] getNodesInvolved() {
			return nodesInvolved;
		}

		private void setNodesInvolved(int[] nodesInvolved) {
			this.nodesInvolved = nodesInvolved;
		}


	}

	public int[] getChildNodes(int node) {
		int[] r = new int[0];
		int[][] nodeChildPointers = currentTree.getNodeChildPointers();
		if(node < nodeChildPointers.length && node >= 0) {
			r = new int[nodeChildPointers[node].length];
			System.arraycopy(nodeChildPointers[node], 0, r, 0, r.length);
		}
		return r;
	}

	public int getRobinsonFouldsDistance(AdvancedTree otherTree) {
		int r = 0;
		String[] taxa = this.getLeafLabels();
		Bipartition[] thisBiparts = this.getBipartitions(taxa);
		Bipartition[] otherBiparts = otherTree.getBipartitions(taxa);
		HashSet<Bipartition> thisSet = new HashSet<Bipartition>();
		for(int i = 0; i < thisBiparts.length; i++) {
			thisSet.add(thisBiparts[i]);
		}
		
		HashSet<Bipartition> otherSet = new HashSet<Bipartition>();
		for(int i = 0; i < otherBiparts.length; i++) {
			otherSet.add(otherBiparts[i]);
		}
		
		HashSet<Bipartition> union = new HashSet<Bipartition>(thisSet);
		union.addAll(otherSet);
		
		HashSet<Bipartition> intersection = new HashSet<Bipartition>(thisSet);
		intersection.retainAll(otherSet);

		r = (union.size() - intersection.size()) / 2;
		return r;
	}

	/**
	 * Associates a key-value pair of metadata with the indicated node. 
	 * Providing null as the value will remove the key from the metadata
	 * set for the node
	 * @param node
	 * @param key
	 * @param value
	 */
	public void addNodeMetadata(int node, String key, String value) {
		if(nodeMetadata == null) {
			nodeMetadata = new HashMap[getLeaves().length];
		}
		if(nodeMetadata[node] == null) {
			nodeMetadata[node] = new HashMap<String, String>();
		}
		if(value == null) {
			nodeMetadata[node].remove(key);
		} else {
			nodeMetadata[node].put(key, value);
		}
		
	}
	
	public void addNodeMetadata(String nodeLabel, String key, String value) {
		int index = getLeafIndex(nodeLabel);
		if(index >= 0) {
			addNodeMetadata(index, key, value);
		}
	} 
	
	private int getLeafIndex(String label) {
		int r = -1;
	    String[] leafLabels = getLeafLabels();
	    for(int i = 0; i < leafLabels.length && r < 0; i++) {
	    	if(leafLabels[i].equals(label)) {
	    		r = i;
	    	}
	    }
        return r;
	}
	
	public String getTreeJSON() {
		String r = null;
		int topNode = currentTree.getTopLevelNode();
		r = getNodeJSON(topNode, 0);
		return r;
	}
	
	public String getNodeJSON(int node, double distanceFromRoot) {
		String r = null;
		int[][] nodeChildPointers = currentTree.getNodeChildPointers();
		AdvancedTree phylogram = new AdvancedTree(this.currentTree);
		phylogram.makePhylogram();
		AdvancedTree cladogram = new AdvancedTree(this.currentTree);
		cladogram.makeCladogram();
		
		double[] branchLengths = currentTree.getBranchLengths();
		double[] phylogramNodeX = phylogram.getNodeDistancesFromRoot();
		double[] cladogramNodeX = cladogram.getNodeDistancesFromRoot();
		int[] branchSupports = getBranchSupports();

		if(distanceFromRoot != phylogramNodeX[node]) {
			System.out.println("AdvancedTree.getNodeJSON(): " + distanceFromRoot + " != " + phylogramNodeX[node]);
		}
		float[][] nodeY = getNodeCoordinates();
		String[] tipLabels = currentTree.getLeaves();
		
		StringBuffer sb = new StringBuffer();
		sb.append("{");
		sb.append("\'px\':");
		sb.append(phylogramNodeX[node]);
		sb.append(",");
		sb.append("cx:");
		sb.append(cladogramNodeX[node]);
		sb.append(",");
		sb.append("py:");
		sb.append(nodeY[node][0]);
		sb.append(",");
		sb.append("s:");
		sb.append(branchSupports[node]);
		sb.append(",");
		sb.append("toRoot:");
		sb.append(distanceFromRoot);
		sb.append(",");
        //if there is metadata for this node, add it
		if(nodeMetadata != null && nodeMetadata[node] != null && nodeMetadata[node].size() > 0) {
			Object[] keys = nodeMetadata[node].keySet().toArray();
			sb.append("md:");
			sb.append("{");
			sb.append(keys[0]);
			sb.append(":");
			sb.append(nodeMetadata[node].get(keys[0]));
			for(int i = 1; i < keys.length; i++) {
				sb.append(",");
				sb.append(keys[i]);
				sb.append(":");
				sb.append(nodeMetadata[node].get(keys[i]));				
			}
			sb.append("}");
			sb.append(",");
		}
		sb.append("n:");
		if(nodeChildPointers[node] == null || nodeChildPointers[node].length == 0) {
			//this is a leaf node
			sb.append("'");
			tipLabels[node] = tipLabels[node].replaceAll("'", "");
//			tipLabels[node] = tipLabels[node].replaceAll("_", " ");
			sb.append(tipLabels[node]);
//			sb.append(node);
			sb.append("'");
		} else {
		    //this is an internal node
			sb.append("'");
			sb.append(node); //use node number for name of internal nodes
			sb.append("'");
			sb.append(",c:[");
			sb.append(getNodeJSON(nodeChildPointers[node][0], distanceFromRoot + branchLengths[nodeChildPointers[node][0]]));
			for(int i = 1; i < nodeChildPointers[node].length; i++) {
				sb.append(",");
				sb.append(getNodeJSON(nodeChildPointers[node][i], distanceFromRoot + branchLengths[nodeChildPointers[node][i]]));
			}
			sb.append("]");
		}
		sb.append("}");
		r = sb.toString();
		return r;
	}

	public AdvancedTree clone() {
		AdvancedTree r = null;
		r = new AdvancedTree((BasicTree)this.getBasicTree().clone());
		return r;
	}
}
