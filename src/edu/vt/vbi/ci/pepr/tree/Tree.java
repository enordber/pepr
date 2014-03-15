package edu.vt.vbi.ci.pepr.tree;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;

import edu.vt.vbi.ci.pepr.alignment.SequenceAlignment;
import edu.vt.vbi.ci.util.ExtendedBitSet;

public class Tree implements TreeI {

	private TreeNodeI[] nodes;
	/*
	 * The nodeToIndexMap is to speed up searches for the index of a
	 * node in the nodes array. A linear search could be used, but the
	 * hash lookup will improve speed. This is not initialized and
	 * populated unless getNodeIndex() is called.
	 */
	private HashMap nodeToIndexMap;

	/*
	 * The alignment used to construct this tree.
	 */
	private SequenceAlignment alignment;
	private HashSet branchSet; //list of TreeBranches
	private TreeBranch[] branches;
	private TreeBranch currentRootBranch;
	private boolean isRooted;
	private HashSet bipartitionTaxa = null;

	/*
	 * nodeDistanceMatrix contains the distances between all pairs of
	 * nodes in the tree. It is populated by createNodeDistanceMatrix().
	 * The index values are the same as in the nodes array.
	 */
	private double[][] nodeDistanceMatrix;

	/**
	 * Location along root branch to place root, if no
	 * other valid root point is given.
	 */
	private double defaultRootPoint = 0.5;

	public Tree(TreeBranch startingBranch) {
		currentRootBranch = startingBranch;
		//currentRootBranch.setAsRoot(0.5);

		//get list of all TreeBranches
		branchSet = new HashSet();
		TreeNodeI[] rootBranchNodes = currentRootBranch.getNodes();
		TreeBranch[] branches0 = rootBranchNodes[0].getOtherBranches(currentRootBranch);
		TreeBranch[] branches1 = rootBranchNodes[1].getOtherBranches(currentRootBranch);
		for(int i = 0; i < branches0.length; i++) {
			branchSet.add(branches0[i]);
			branchSet.addAll(branches0[i].getDescendantBranches(rootBranchNodes[0]));
		}
		for(int i = 0; i < branches1.length; i++) {
			branchSet.add(branches1[i]);
			branchSet.addAll(branches1[i].getDescendantBranches(rootBranchNodes[1]));
		}
		branchSet.add(currentRootBranch);

		//get array of all TreeNodes
		ArrayList nodeList = new ArrayList();
		TreeNode[] rootNodes = currentRootBranch.getNodes();
		nodeList.add(rootNodes[0]);
		nodeList.addAll(rootNodes[0].getDescendantNodes(currentRootBranch));
		nodeList.add(rootNodes[1]);
		nodeList.addAll(rootNodes[1].getDescendantNodes(currentRootBranch));
		nodes = new TreeNodeI[nodeList.size()];
		nodeList.toArray(nodes);

	}


	/* (non-Javadoc)
	 * @see edu.vt.vbi.ci.nongraphic.tree.TreeI#getNodeCount()
	 */
	public int getNodeCount() {
		return nodes.length;
	}

	/* (non-Javadoc)
	 * @see edu.vt.vbi.ci.nongraphic.tree.TreeI#getBranchCount()
	 */
	public int getBranchCount() {
		return branchSet.size();
	}

	/* (non-Javadoc)
	 * @see edu.vt.vbi.ci.nongraphic.tree.TreeI#getNodes()
	 */
	public TreeNodeI[] getNodes() {
		TreeNodeI[] r = new TreeNodeI[nodes.length];
		System.arraycopy(nodes, 0, r, 0, nodes.length);
		return r;
	}

	/* (non-Javadoc)
	 * @see edu.vt.vbi.ci.nongraphic.tree.TreeI#getLeafNodes()
	 */
	public TreeNodeI[] getLeafNodes() {
		TreeNodeI[] r = null;
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

	public double[][] getLeafDistanceMatrix() {
		double[][] r = null;
//		HashMap nodeToIndex = new HashMap();
		TreeNodeI[] leaves = getLeafNodes();
//
//		for(int i = 0; i < leaves.length; i++) {
//			nodeToIndex.put(leaves[i], new Integer(i));
//		}
//
//		r = new double[leaves.length][leaves.length];
//		for(int i = 0; i < r.length; i++) {
//			TreeNodeI leafi = leaves[i];
//			HashMap nodeToDistanceMap = new HashMap();
//		}

		r = getNodeDistanceMatrix(leaves);
		return r;
	}

	/* (non-Javadoc)
	 * @see edu.vt.vbi.ci.nongraphic.tree.TreeI#getInternalNodes()
	 */
	public TreeNode[] getInternalNodes() {
		TreeNode[] r;
		ArrayList internalNodes = new ArrayList(nodes.length/2 + 1);
		for(int i = 0; i < nodes.length; i++) {
			if(!nodes[i].isLeaf()) {
				internalNodes.add(nodes[i]);
			}
		}
		r = new TreeNode[internalNodes.size()];
		internalNodes.toArray(r);
		return r;
	}

	/* (non-Javadoc)
	 * @see edu.vt.vbi.ci.nongraphic.tree.TreeI#getLargestRepresentativeSubset(java.util.HashMap)
	 */
	public TreeNodeI[][] getLargestRepresentativeSubset(HashMap idToGenomeMap) {
		TreeNodeI[][] r = null;

		//try setting each branch as root successively
		TreeBranch[] branches = new TreeBranch[branchSet.size()];
		branchSet.toArray(branches);
		int[][] branchCutSizes = new int[branches.length][2];
		boolean[][] branchCutRepresentative = new boolean[branches.length][2];

		for(int i = 0; i < branches.length; i++) {
			TreeNode[] rootNodes = branches[i].getNodes();
			branches[i].setAsRoot(0.5);
			int genomeSet0 = rootNodes[0].getGenomeSet(branches[i], idToGenomeMap).size();
			int genomeSet1 = rootNodes[1].getGenomeSet(branches[i], idToGenomeMap).size();
			int leafCount0 = rootNodes[0].getLeafCount(branches[i]);
			int leafCount1 = rootNodes[1].getLeafCount(branches[i]);
			branchCutSizes[i][0] = leafCount0;
			branchCutSizes[i][1] = leafCount1;
			if(genomeSet0 == leafCount0) {
				branchCutRepresentative[i][0] = true;
			} else {
				branchCutRepresentative[i][0] = false;
			}
			if(genomeSet1 == leafCount1) {
				branchCutRepresentative[i][1] = true;
			} else {
				branchCutRepresentative[i][1] = false;
			}

		}

		//find best branch to cut
		int sizeOfLargestRepresentative = 0;
		int indexOfBranchForLargestRepresentative = -1;
		for(int i = 0; i < branchCutSizes.length; i++) {
			if(branchCutRepresentative[i][0]) {
				if(branchCutSizes[i][0] > sizeOfLargestRepresentative) {
					sizeOfLargestRepresentative = branchCutSizes[i][0];
					indexOfBranchForLargestRepresentative = i;
				}
			}
			if(branchCutRepresentative[i][1]) {
				if(branchCutSizes[i][1] > sizeOfLargestRepresentative) {
					sizeOfLargestRepresentative = branchCutSizes[i][1];
					indexOfBranchForLargestRepresentative = i;
				}
			}
		}

		branches[indexOfBranchForLargestRepresentative].setAsRoot(0.5);

		TreeNodeI[] candidateNodes = 
			branches[indexOfBranchForLargestRepresentative].getNodes();
		TreeNodeI[] leafNodes0;
		if(candidateNodes[0].isLeaf()) {
			leafNodes0 = new TreeNodeI[]{candidateNodes[0]};
		}else {
			leafNodes0 = candidateNodes[0].getDescendantLeaves(branches[indexOfBranchForLargestRepresentative]);
		}
		TreeNodeI[] leafNodes1;
		if(candidateNodes[1].isLeaf()) {
			leafNodes1 = new TreeNodeI[]{candidateNodes[1]};
		}else {
			leafNodes1 = candidateNodes[1].getDescendantLeaves(branches[indexOfBranchForLargestRepresentative]);
		}

		if(branchCutRepresentative[indexOfBranchForLargestRepresentative][0] &&
				branchCutRepresentative[indexOfBranchForLargestRepresentative][1]) {
			if(leafNodes0.length >= leafNodes1.length) {
				r = new TreeNodeI[][]{leafNodes0, leafNodes1};
			} else {
				r = new TreeNodeI[][]{leafNodes1, leafNodes0};
			}

		} else {
			if(branchCutRepresentative[indexOfBranchForLargestRepresentative][0]) {
				r = new TreeNodeI[][]{leafNodes0};
			} else if(branchCutRepresentative[indexOfBranchForLargestRepresentative][1]){
				r = new TreeNodeI[][]{leafNodes1};
			}
		}
		return r;
	}

	/* (non-Javadoc)
	 * @see edu.vt.vbi.ci.nongraphic.tree.TreeI#getNodeDefinitions(java.util.HashMap)
	 */
	public BitSet[] getNodeDefinitions(HashMap idToGenomeMap) {
		// TODO Auto-generated method stub
		return null;
	}

	public TreeBranch[] getBranches() {
		TreeBranch[] branches = new TreeBranch[branchSet.size()];
		branchSet.toArray(branches);
		return branches;
	}

	/* (non-Javadoc)
	 * @see edu.vt.vbi.ci.nongraphic.tree.TreeI#getCladeTaxa(java.util.HashMap)
	 */
	public HashSet[] getCladeTaxa(HashMap idToGenomeMap) {
		TreeBranch[] branches = new TreeBranch[branchSet.size()];
		branchSet.toArray(branches);
		HashSet[] r = new HashSet[branches.length*2];

		for(int i = 0; i < branches.length; i++) {
			TreeNode[] branchNodes = branches[i].getNodes();
			branches[i].setAsRoot(0.5);
			r[i*2] = branchNodes[0].getGenomeSet(branches[i], idToGenomeMap);
			r[(i*2)+1] = branchNodes[1].getGenomeSet(branches[i], idToGenomeMap);
		}
		return r;
	}

	/**
	 * Returns the set of bipartitions of taxa represented by this Tree.
	 * The returned HashSet contains HashSet[] elements. Each element is a
	 * pair of HashSets containing the names of the taxa represented by
	 * each side of the branches on the tree. The idToGenomeMap should
	 * contain an entry for each leaf name in the tree, mapping it to
	 * the taxon name that should be used for that name.
	 */
	public HashSet getBipartitionTaxa(HashMap idToGenomeMap) {
		if(bipartitionTaxa == null) {
			//populate bipartitionTaxa
			bipartitionTaxa = new HashSet();
			TreeBranch[] branches = new TreeBranch[branchSet.size()];
			branchSet.toArray(branches);
			BipartHashSetComparator comparator = new BipartHashSetComparator();

			for(int i = 0; i < branches.length; i++) {
				TreeNode[] branchNodes = branches[i].getNodes();
				root(branches[i]);
				HashSet[] bipartPair = new HashSet[2];
				bipartPair[0] = branchNodes[0].getGenomeSet(branches[i], idToGenomeMap);
				bipartPair[1] = branchNodes[1].getGenomeSet(branches[i], idToGenomeMap);
				Arrays.sort(bipartPair, comparator);
				bipartitionTaxa.add(bipartPair);
			}

		}
		return bipartitionTaxa;
	}

	/**
	 * Returns an ExtendedBitSet for each clade representing the presence
	 * and absence of each taxon in that clade. This is different from a
	 * bipartition, in that a taxon may appear on both sides of a branch.
	 * Each side of the branch is represented by a different ExtendedBitSet
	 * in the result. Multiple equivalent ExtendedBitSets may be present. 
	 * The taxa are represented in the result sets in the same order as given 
	 * in the taxonList. There may be taxa in the taxonList that are not in
	 * this Tree, in which case they will show as absent in all of the 
	 * result sets. 
	 * @return
	 */
	public Bipartition[] getBipartitions(HashMap idToTaxonMap, String[] taxonNames) {
		Bipartition[] r = null;
		HashMap taxonToIndexMap = new HashMap();
		for(int i = 0; i < taxonNames.length; i++) {
			taxonToIndexMap.put(taxonNames[i], new Integer(i));
		}

		HashSet bipartSet = getBipartitionTaxa(idToTaxonMap);
		HashSet[][] biparts = new HashSet[bipartSet.size()][];
		bipartSet.toArray(biparts);
		ArrayList bpList = new ArrayList();
		HashSet sidesSeen = new HashSet();
		for(int i = 0; i < biparts.length; i++) {
			if(biparts[i][0].size() > 1 && biparts[i][1].size() > 1) { 
				String[] sideA = new String[biparts[i][0].size()];
				biparts[i][0].toArray(sideA);
				String[] sideB = new String[biparts[i][1].size()];
				biparts[i][1].toArray(sideB);

				ExtendedBitSet bsa = new ExtendedBitSet();
				for(int j = 0; j < sideA.length; j++) {
					bsa.set(((Integer)taxonToIndexMap.get(sideA[j])).intValue());
				}

				ExtendedBitSet bsb = new ExtendedBitSet();
				for(int j = 0; j < sideB.length; j++) {
					bsb.set(((Integer)taxonToIndexMap.get(sideB[j])).intValue());
				}
				if(!sidesSeen.contains(bsa) && !sidesSeen.contains(bsb)) { 
					ExtendedBitSet participatingTaxa = bsa.getOr(bsb);
					Bipartition bp = new Bipartition(bsa, bsb, taxonNames.length);
					bp.setParticipatingTaxonSet(participatingTaxa);
					bpList.add(bp);
					sidesSeen.add(bsa);
					sidesSeen.add(bsb);
				} 
			}
		}

		r = new Bipartition[bpList.size()];
		bpList.toArray(r);
		return r;
	}

	/* (non-Javadoc)
	 * @see edu.vt.vbi.ci.nongraphic.tree.TreeI#isRooted()
	 */
	public boolean isRooted() {
		return isRooted;
	}


	/* (non-Javadoc)
	 * @see edu.vt.vbi.ci.nongraphic.tree.TreeI#root(edu.vt.vbi.ci.nongraphic.tree.TreeBranch, double)
	 */
	public void root(TreeBranch b, double rootPoint) {
		//set the currentRootBranch to no longer be root
		if(currentRootBranch != null) {
			currentRootBranch.setAsNotRoot();
		}
		//make sure rootPoint is valid. if not, defaultRootPoint
		if(rootPoint < 0.0 || rootPoint > 1.0) {
			rootPoint = defaultRootPoint;
		}
		//set the new root branch as root
		b.setAsRoot(rootPoint);
		currentRootBranch = b;
	}

	/* (non-Javadoc)
	 * @see edu.vt.vbi.ci.nongraphic.tree.TreeI#root(edu.vt.vbi.ci.nongraphic.tree.TreeBranch)
	 */
	public void root(TreeBranch b) {
		root(b, defaultRootPoint);
	}

	/* (non-Javadoc)
	 * @see edu.vt.vbi.ci.nongraphic.tree.TreeI#getRootBranch()
	 */
	public TreeBranch getRootBranch() {
		return currentRootBranch;
	}

	/* (non-Javadoc)
	 * @see edu.vt.vbi.ci.nongraphic.tree.TreeI#getTotalDistance()
	 */
	public double getTotalDistance() {
		double r = 0;
		if(isRooted()) {
			//r = currentRootNode.getTotalDistance();
		} else {
			r = currentRootBranch.getTotalDistance(null);
		}
		return r;
	}

	/* (non-Javadoc)
	 * @see edu.vt.vbi.ci.nongraphic.tree.TreeI#collapseToMaximumLeaves(int)
	 */
	public void collapseToMaximumLeaves(int maxTips) {
		Comparator collapseComparator = new CollapseComparator();
		TreeNode[] iNodes = getInternalNodes();
		Arrays.sort(iNodes, collapseComparator);
		System.out.println("Tree.collapseToMaximumLeaves() tips before collapsing:"
				+ getTipCount());		
		//while there are still too many tips, find the node
		//with the smallest mean child distance that has only
		//tips as children and collapse it
		while(getTipCount() > maxTips) {
			int nextCollapseIndex = 0;
			boolean nodeToCollapseNotFound = true;
			while(nodeToCollapseNotFound) {
				if(!iNodes[nextCollapseIndex].isCollapsed() &&
						iNodes[nextCollapseIndex].getTipCount(null) == 
							iNodes[nextCollapseIndex].getChildNodes(null).length) {
					nodeToCollapseNotFound = false;
				} else {
					nextCollapseIndex++;
				}
			}
			System.out.println("collapse " 
					+ getTipCount() + ": "+ iNodes[nextCollapseIndex].toString());
			iNodes[nextCollapseIndex].setCollapsed(true);
		}
	}

	/* (non-Javadoc)
	 * @see edu.vt.vbi.ci.nongraphic.tree.TreeI#getLeafCount()
	 */
	public int getLeafCount() {
		int r = 0;
		if(isRooted()) {
			//r = currentRootNode.getLeafCount();
		} else {
			TreeNode[] nodes = currentRootBranch.getNodes();
			r += nodes[0].getLeafCount(currentRootBranch);
			r += nodes[1].getLeafCount(currentRootBranch); 
		}
		return r;
	}

	/* (non-Javadoc)
	 * @see edu.vt.vbi.ci.nongraphic.tree.TreeI#getMaxTipDistance()
	 */
	public double getMaxTipDistance() {
		double maxDistance = -1;
		TreeNodeI[] leaves = getLeafNodes();
		if(isRooted()) {
//			for(int i = 0; i < leaves.length; i++) {
//			double distance = 
//			currentRootNode.getDistanceToDescendant(leaves[i]);
//			maxDistance = Math.max(maxDistance, distance);
//			}
		} else {
			for(int i = 0; i < leaves.length; i++) {
				double distance = leaves[i].getDistanceFromRoot();
				if(!Double.isNaN(distance)) {
					maxDistance = Math.max(distance, maxDistance);
				}
				System.out.print("");
			}
		}
		return maxDistance;
	}

	/**
	 * Returns the number of tips, where a tip
	 * is either a leaf with no nodes collapsed 
	 * between the leaf and the root, or a collapsed
	 * node, regardless of the number of leaves
	 * below that node.
	 * @return
	 */
	int getTipCount() {
		int r = -1;
		if(isRooted()) {
			//r = currentRootNode.getTipCount();
		} else {
			TreeNode[] nodes = currentRootBranch.getNodes();
			r = nodes[0].getTipCount(currentRootBranch) + 
			nodes[1].getTipCount(currentRootBranch);
		}
		return r;
	}

	/**
	 * Returns a HashSet of HashSets. Each lower-level
	 * HashSet contains the leaves from a single clade
	 * in the tree.
	 * @return
	 */
	public HashSet getCladeSet() {
		HashSet r = new HashSet();
		TreeBranch[] branches = getBranches();
		for(int i = 0; i < branches.length; i++) {
			String[][] biparts = branches[i].getBipartition();
			HashSet s1 = new HashSet();
			Arrays.sort(biparts[0]);
			for(int j = 0; j < biparts[0].length; j++) {
				s1.add(biparts[0][j].intern());
			}
			HashSet s2 = new HashSet();
			Arrays.sort(biparts[1]);
			for(int j = 0; j < biparts[1].length; j++) {
				s2.add(biparts[1][j].intern());
			}
			r.add(s1);
			r.add(s2);
		}
		return r;
	}

	/* (non-Javadoc)
	 * @see edu.vt.vbi.ci.nongraphic.tree.TreeI#getTipsLaderizeUp()
	 */
	public TreeNodeI[] getTipsLaderizeUp() {
		TreeNodeI[] r = null;
		ArrayList tips = new ArrayList();
		TreeNode[] nodes = currentRootBranch.getNodes();
		if(nodes[0].getLeafCount(currentRootBranch) > nodes[1] .getLeafCount(currentRootBranch)) {
			tips.addAll(nodes[0].getTipsLaderizeUp(currentRootBranch));
			tips.addAll(nodes[1].getTipsLaderizeUp(currentRootBranch));
		} else {
			tips.addAll(nodes[1].getTipsLaderizeUp(currentRootBranch));
			tips.addAll(nodes[0].getTipsLaderizeUp(currentRootBranch));
		}
		r = new TreeNodeI[tips.size()];
		tips.toArray(r);

		return r;
	}

	/* (non-Javadoc)
	 * @see edu.vt.vbi.ci.nongraphic.tree.TreeI#getTipsLaderizeDown()
	 */
	public TreeNodeI[] getTipsLaderizeDown() {
		TreeNodeI[] r = null;
		ArrayList tips = new ArrayList();
		TreeNode[] nodes = currentRootBranch.getNodes();
		if(nodes[0].getLeafCount(currentRootBranch) > nodes[1].getLeafCount(currentRootBranch)) {
			tips.addAll(nodes[1].getTipsLaderizeDown(currentRootBranch));
			tips.addAll(nodes[0].getTipsLaderizeDown(currentRootBranch));
		} else {
			tips.addAll(nodes[0].getTipsLaderizeDown(currentRootBranch));
			tips.addAll(nodes[1].getTipsLaderizeDown(currentRootBranch));
		}
		r = new TreeNodeI[tips.size()];
		tips.toArray(r);

		return r;
	}

	/* (non-Javadoc)
	 * @see edu.vt.vbi.ci.nongraphic.tree.TreeI#getTipsBalance(boolean)
	 */
	public TreeNodeI[] getTipsBalance(boolean startUp) {
		TreeNodeI[] r = null;
		ArrayList tips = new ArrayList();
		TreeNode[] nodes = currentRootBranch.getNodes();
		if(nodes[0].getLeafCount(currentRootBranch) > nodes[1].getLeafCount(currentRootBranch) && startUp) {
			tips.addAll(nodes[0].getTipsBalance(!startUp, currentRootBranch));
			tips.addAll(nodes[1].getTipsBalance(!startUp, currentRootBranch));
		} else {
			tips.addAll(nodes[1].getTipsBalance(!startUp, currentRootBranch));
			tips.addAll(nodes[0].getTipsBalance(!startUp, currentRootBranch));
		}
		r = new TreeNodeI[tips.size()];
		tips.toArray(r);

		return r;
	}

	/* (non-Javadoc)
	 * @see edu.vt.vbi.ci.nongraphic.tree.TreeI#getRobinsonFouldsDistance(edu.vt.vbi.ci.nongraphic.tree.Tree)
	 */
	public int getRobinsonFouldsDistance(TreeI otherTree) {
		int r = -1;
		HashSet thisTreeClades = this.getCladeSet();
		HashSet otherTreeClades = otherTree.getCladeSet();

		HashSet unionSet = new HashSet(thisTreeClades);
		unionSet.addAll(otherTreeClades);

		HashSet intersectionSet = new HashSet(thisTreeClades);
		intersectionSet.retainAll(otherTreeClades);

		r = unionSet.size() - intersectionSet.size();
		r /= 4;
		return r;
	}

	/**
	 * Returns a matrix of distances for the given list of nodes. The order
	 * of nodes in the returned matrix is the same as the order in the 
	 * input array.
	 * 
	 * @param nodes
	 * @return
	 */
	public double[][] getNodeDistanceMatrix(TreeNodeI[] nodes) {
		double[][] r = null;
		if(nodeDistanceMatrix == null) {
			createNodeDistanceMatrix();
		}

		int[] indexMap = new int[nodes.length];
		for(int i = 0; i < indexMap.length; i++) {
			indexMap[i] = getNodeIndex(nodes[i]);
		}

		r = new double[nodes.length][nodes.length];
		for(int i = 0; i < r.length; i++) {
			for(int j = 0; j < r[i].length; j++) {
				r[i][j] = nodeDistanceMatrix[indexMap[i]][indexMap[j]];
			}
		}

		return r;
	}

	/**
	 * Initializes and populates the nodeDistanceMatrix with the distances
	 * for all pairs of nodes in the tree.
	 */
	private void createNodeDistanceMatrix() {
		nodeDistanceMatrix = new double[nodes.length][nodes.length];
		if(branches == null) {
			branches = getBranches();
		}
		//start with values for the node pairs connected by a branch
		for(int i = 0; i < branches.length; i++) {
			//get the node indices for this branch
			TreeNodeI[] branchNodes = branches[i].getNodes();
			int nodeAIndex = getNodeIndex(branchNodes[0]);
			int nodeBIndex = getNodeIndex(branchNodes[1]);

			//get the length for this branch
			double branchLength = branches[i].getBranchLength();
			nodeDistanceMatrix[nodeAIndex][nodeBIndex] = branchLength;
			nodeDistanceMatrix[nodeBIndex][nodeAIndex] = branchLength;
		}

		//fill in the remaining cells in the matrix based on 
		//these initial pairs
		boolean emptyCellsRemain = true;
		while(emptyCellsRemain) {
			emptyCellsRemain = false;

			//check each pair. if not filled, see if it can be filled this
			//round.
			for(int i = 0; i < nodes.length; i++) {
				for(int j = i+1; j < nodes.length; j++) {
					if(nodeDistanceMatrix[i][j] == 0) {
						double minDistance = 0;
						for(int k = 0; k < nodes.length; k++) {
							double distance = 0;
							if(nodeDistanceMatrix[i][k] > 0 
									&& nodeDistanceMatrix[j][k] > 0) {
								distance = nodeDistanceMatrix[i][k]
								              + nodeDistanceMatrix[j][k];
							}
							if(distance > 0) {
								if(minDistance == 0) {
									minDistance = distance;
								} else {
								    minDistance = Math.min(distance, minDistance);
								}
							}
						}
						nodeDistanceMatrix[i][j] = minDistance;
						nodeDistanceMatrix[j][i] = minDistance;
						emptyCellsRemain |= nodeDistanceMatrix[i][j] == 0;
					}
				}
			}
		}
	}

	private int getNodeIndex(TreeNodeI node) {
		int r = -1;
		if(nodeToIndexMap == null) {
			//populate the nodeToIndexMap
			nodeToIndexMap = new HashMap();
			for(int i = 0; i < nodes.length; i++) {
				nodeToIndexMap.put(nodes[i], new Integer(i));
			}
		}

		Integer index = (Integer)nodeToIndexMap.get(node);
		if(index != null) {
			r = index.intValue();
		}
		return r;
	}

	public String toString() {
		String r = null;
		r = currentRootBranch.toString(null, true, true) + ";";
		return r;
	}

	/* (non-Javadoc)
	 * @see edu.vt.vbi.ci.nongraphic.tree.TreeI#getTreeString(boolean, boolean)
	 */
	public String getTreeString(boolean includeDistance, boolean includeSupport) {
		String r = null;
		r = currentRootBranch.toString(null, includeDistance, includeSupport) + ";";
		return r;
	}

	/* (non-Javadoc)
	 * @see edu.vt.vbi.ci.nongraphic.tree.TreeI#getTreeStringTrifurcating()
	 */
	public String getTreeStringTrifurcating() {
		String r = null;

		//determine which node flanking the root branch will be 
		//removed
		TreeNodeI splitNode;
		TreeNode nonSplitNode;
		TreeNode[] splitCandidates = currentRootBranch.getNodes();
		//split the one with the greatest number of leaves
		if(splitCandidates[0].getLeafCount(currentRootBranch) >
		splitCandidates[1].getLeafCount(currentRootBranch)) {
			splitNode = splitCandidates[0];
			nonSplitNode = splitCandidates[1];
		} else {
			splitNode = splitCandidates[1];
			nonSplitNode = splitCandidates[0];
		}

		//TODO there could be a bug here is the splitCandidate node is
		//multifurcating
		TreeBranch[] splitBranches = splitNode.getOtherBranches(currentRootBranch);
		if(splitBranches.length != 2) {
			System.out.println("Tree.getTreeStringTrifurcating() wrong number" +
					"of splitBranches. Should be 2, is " + splitBranches.length);
		}
		String part1 = nonSplitNode.toString(currentRootBranch, false, false);
		String part2 = splitBranches[0].toString(splitNode, false, false);
		String part3 = splitBranches[1].toString(splitNode, false, false);

		r = "(" + part1 + "," + part2 + "," + part3 + ");";
		return r;
	}

	private static class CollapseComparator implements Comparator {

		public int compare(Object arg0, Object arg1) {
			int r;
			TreeNode nodeA = (TreeNode) arg0;
			TreeNode nodeB = (TreeNode) arg1;
			double distanceA = nodeA.getMeanChildDistance();
			double distanceB = nodeB.getMeanChildDistance();
			double diff = distanceA - distanceB;
			if(diff < 0) {
				r = -1;
			} else if(diff > 0) {
				r = 1;
			} else {
				r = 0;
			}
			return r;
		}};

		/**
		 * Compares HashSets on the basis of size. If the sizes are
		 * equal, comparison is made on the basis of the hash codes.
		 * 
		 * @author enordber
		 *
		 */
		private static class BipartHashSetComparator implements Comparator{

			public int compare(Object arg0, Object arg1) {
				int r;

				HashSet s0 = (HashSet)arg0;
				HashSet s1 = (HashSet)arg1;
				r = s0.size() - s1.size();
				if(r == 0) {
					r = s0.hashCode() - s1.hashCode();
				}
				if(r == 0) {
					System.out.println("Tree.BipartHashSetComparator.compare() " +
					"comparator still finds 0");
				}
				return r;
			}

		}


		/**
		 * Returns the alignment used to construct this tree,
		 * if it is available. If the alignment has not been set 
		 * (using setAlignment()) this method returns null
		 */
		public SequenceAlignment getAlignment() {
			return alignment;
		}


		public void setAlignment(SequenceAlignment alignment) {
			this.alignment = alignment;
		}

}
