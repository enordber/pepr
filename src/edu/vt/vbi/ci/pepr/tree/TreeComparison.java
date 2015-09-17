package edu.vt.vbi.ci.pepr.tree;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

import edu.vt.vbi.ci.pepr.alignment.SequenceAlignment;
import edu.vt.vbi.ci.pepr.alignment.SequenceAlignmentParser;
import edu.vt.vbi.ci.pepr.stats.StatisticsUtilities;
import edu.vt.vbi.ci.pepr.tree.pipeline.PhyloPipeline;
import edu.vt.vbi.ci.util.CommandLineProperties;
import edu.vt.vbi.ci.util.CommandResults;
import edu.vt.vbi.ci.util.ExecUtilities;
import edu.vt.vbi.ci.util.ExtendedBitSet;
import edu.vt.vbi.ci.util.HandyConstants;
import edu.vt.vbi.ci.util.file.FastaSequenceFile;
import edu.vt.vbi.ci.util.file.TextFile;

public class TreeComparison {

	private TreeI standardTree;
	private TreeI[] comparisonTrees;
	private int[] treeDistancesFromStandard;

	public static void main(String[] args) {
		TreeComparison tc = new TreeComparison(args);
	}

	public TreeComparison(String[] args) {
		CommandLineProperties clp = new CommandLineProperties(args);
		PhyloPipeline.setCommandPaths();
		int processors = Runtime.getRuntime().availableProcessors();
		String procString = clp.getValues(HandyConstants.MAX_CONCURRENT_PROCESS_PARAM, ""+processors)[0];
		processors = Integer.parseInt(procString);

		String matrix = clp.getValues(HandyConstants.ML_MATRIX, "PROTGAMMALGF")[0];
		String alignmentFileName = clp.getValues(HandyConstants.ALIGNMENT_FILE_COMMAND, null)[0];
		if(alignmentFileName == null) {
			alignmentFileName = clp.getValues(HandyConstants.ALIGN_FILE, null)[0];
		}

		String[] treeFileNames = clp.getValues(HandyConstants.TREE_FILE);
		boolean doRF = clp.getValues(HandyConstants.RF, HandyConstants.TRUE)[0].equalsIgnoreCase(HandyConstants.TRUE);
		boolean doConsel = clp.getValues(HandyConstants.CONSEL, HandyConstants.TRUE)[0].equalsIgnoreCase(HandyConstants.TRUE);
		boolean doTreeDist = clp.getValues(HandyConstants.TREE_DIST, HandyConstants.TRUE)[0].equalsIgnoreCase(HandyConstants.TRUE);
		boolean useLengths = clp.getValues(HandyConstants.USE_LENGTHS, HandyConstants.TRUE)[0].equalsIgnoreCase(HandyConstants.TRUE);
		boolean compareSupports = clp.getValues(HandyConstants.SUPPORT_REPS, HandyConstants.FALSE)[0].equalsIgnoreCase(HandyConstants.TRUE);


		//if outgroup files are given, then remove outgroup taxa from the trees and the alignment before proceeding
		String[] outgroupFileNames = clp.getValues(HandyConstants.OUTGROUP);
		if(outgroupFileNames == null) {
			outgroupFileNames = new String[0];
		}


		try {
			TextFile[] treeFiles = new TextFile[treeFileNames.length];
			for(int i = 0; i < treeFiles.length; i++) {
				System.out.println("treeFile " + i + ": " + treeFileNames[i]);
				treeFiles[i] = new TextFile(treeFileNames[i]);
			}

			ArrayList treeList = new ArrayList(treeFiles.length);
			for(int i = 0; i < treeFiles.length;i ++) {
				String[] treeStrings = treeFiles[i].getAllLines();
				for(int j = 0; j < treeStrings.length; j++) {
					if(treeStrings[j].length() > 0) {
						BasicTree tree = new BasicTree(treeStrings[j]);
						//						TreeI tree = TreeParser.parseTreeString(treeStrings[j]);
						treeList.add(tree);
					}
				}
			}
			BasicTree[] initialTrees = new BasicTree[treeList.size()];
			treeList.toArray(initialTrees);

			String[] outgroupTaxa = new String[0];
			try {
				FastaSequenceFile[] outgroupFiles = 
					loadSequenceFiles(outgroupFileNames);
				outgroupTaxa = getDistinctTaxa(outgroupFiles);
				Arrays.sort(outgroupTaxa);
				for(int i = 0; i < outgroupTaxa.length; i++) {
					System.out.println("" + outgroupTaxa[i]);
					//remove each outgroup taxon from each tree
					for(int j = 0; j < initialTrees.length; j++) {
						initialTrees[j].removeTaxon(outgroupTaxa[i]);
						initialTrees[j] = new BasicTree(initialTrees[j].getTreeString(useLengths, true));
					}
				}
			} catch(IOException ioe) {
				ioe.printStackTrace();
			}

			boolean commonOnly = true;
			if(commonOnly) {
				HashSet[] treeTaxa = new HashSet[initialTrees.length];
				for(int i = 0; i < initialTrees.length; i++) {
					treeTaxa[i] = new HashSet();
					String[] taxa = initialTrees[i].getLeaves();
					for(int j = 0; j < taxa.length; j++) {
						treeTaxa[i].add(taxa[j]);
					}
				}

				HashSet commonTaxaSet = new HashSet(treeTaxa[0]);
				for(int i = 1; i < treeTaxa.length; i++) {
					commonTaxaSet.retainAll(treeTaxa[i]);
				}

				for(int i = 0; i < initialTrees.length; i++) {
					String[] taxa = initialTrees[i].getLeaves();
					for(int j = 0; j < taxa.length; j++) {
						if(!commonTaxaSet.contains(taxa[j])) {
							System.out.println("tree " + i +
									" removing taxon: " + taxa[j]);
							initialTrees[i].removeTaxon(taxa[j]);
							initialTrees[i] = new BasicTree(initialTrees[i].getTreeString(useLengths, true));					}
					}
				}

			}
			TreeI[] trees = new TreeI[initialTrees.length];
			AdvancedTree[] advancedTrees = new AdvancedTree[initialTrees.length];
			//System.out.println("parsing " + advancedTrees.length + " trees");
			for(int i = 0; i < initialTrees.length; i++) {
				trees[i] = 
					TreeParser.parseTreeString(initialTrees[i].getTreeString(useLengths, true));
				advancedTrees[i] = new AdvancedTree(initialTrees[i]);
			}
			System.out.println("done parsing");

			if(compareSupports) {
				compareBranchSupports(new AdvancedTree(initialTrees[0]),
						new AdvancedTree(initialTrees[1]));
			}

			if(doRF) {
				compareAllvsAll(trees);
//				compareAllvsAll(advancedTrees);
			}

			if(doConsel) {

				TextFile alignmentFile = new TextFile(alignmentFileName);
				String alignmentString = alignmentFile.toString();
				SequenceAlignment alignment = null;
				if(alignmentString.startsWith(">")) {
					//Fasta alignment
					alignment = SequenceAlignmentParser.parseFastaAlignmentFile(alignmentFileName);
				} else {
					//Phylip alignment
					alignment = SequenceAlignmentParser.parsePhylipAlignment(alignmentString);
				}
				alignment.setName(alignmentFile.getFile().getName());

				if(outgroupTaxa != null && outgroupTaxa.length > 0) {
					//remove each outgroup sequence from the alignment
					for(int i = 0; i < outgroupTaxa.length; i++) {
						int index =  -1;
						if(alignment.hasTaxonNames()) {
							index = alignment.getTaxonIndex(outgroupTaxa[i]);
						} else {
							//find index for sequence name
							index = alignment.getSequenceIndex(outgroupTaxa[i]);
						}
						alignment.removeSequence(index);
					}
				}
				double gapOrMissing = alignment.getProportionGapOrMissing();
				int totalChars = alignment.getNTax() * alignment.getLength();
				int nonGapOrMissing = (int) (totalChars * (1-gapOrMissing));
				System.out.println("Alignment: " + alignment.getName() + 
						"\n taxa: " + alignment.getNTax() + 
						"\n length: "+ 
						alignment.getLength() + 
						"\n total characters: " + totalChars +
						"\n proportion of missing characters or gaps: " + 
						alignment.getProportionGapOrMissing() + 
						"\n number of non-gap non-missing characters: " + nonGapOrMissing);

				String[] conselResults = runConsel(alignment, initialTrees, procString, matrix);
				System.out.println("Results:");
				for(int i = 0; i < conselResults.length; i++) {
					System.out.println(conselResults[i]);
				}
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	public TreeComparison() {
		String standardTreeFileName = "/Users/enordber/vbi/brucella/brucella_full_tree.nwk";
		standardTreeFileName = "/Users/enordber/vbi/brucella/full_tree_2";
		standardTreeFileName = "/Users/enordber/vbi/brucella/fullTree100k/brucellaFullTree100k.nwk";

		String comparisonTreeFileName = 
			"/Users/enordber/vbi/brucella/individualTrees/brucellaIndividualTreeLog";
		comparisonTreeFileName = 
			"/Users/enordber/vbi/fastTree/enteroTrees_2.nwk";

		String fileType = HandyConstants.TAB_DELIMITED;

		String mrTfileName = 
			"/Users/enordber/vbi/brucella/Brucella_full.nex.run1.t";

		//		comparisonTreeFileName = mrTfileName;
		//		fileType = HandyConstants.MR_T_FILE;

		try {
			if(fileType.equals(HandyConstants.TAB_DELIMITED)) {
				loadComparisonTreesFromFile(comparisonTreeFileName, 2);
			} else if(fileType.equals(HandyConstants.MR_T_FILE)) {
				loadComparisonTreesFromMrTFile(comparisonTreeFileName);
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

//				compareAllvsAll(comparisonTrees);

		compareSlidingWindow(comparisonTrees, 10);

		//		int[] rollingRF = rollingComparison(comparisonTrees);

		//		for(int i = 0; i < rollingRF.length; i++) {
		//			System.out.println("" + i + ": " + rollingRF[i]);
		//		}
	}

	private void compareSlidingWindow(TreeI[] trees, int windowSize) {

		for(int i = windowSize; i < trees.length; i++) {
			int min = Integer.MAX_VALUE;
			int max = 0;
			double sum = 0;
			int count = 0;
			if(trees[i] != null) {
				for(int j = i - windowSize; j < i ; j++) {
					if(trees[j] != null) {
						int distance = trees[i].getRobinsonFouldsDistance(trees[j]);
						min = Math.min(min, distance);
						max = Math.max(max, distance);
						sum += distance;
						count++;
						System.out.println(i + " vs " + (j) + ": " + distance);
					}
				}
				double mean = sum / count;
				System.out.println(i + "(" + count + ")  min: " + min + " max: " + max 
						+ " mean: " + mean);
			}
		}

	}

	private int[][] compareAllvsAll(TreeI[] trees) {
		int[][] r = new int[trees.length][trees.length];

		for(int i = 0; i < trees.length; i++) {
			for(int j = i+1; j < trees.length; j++) {
				if(trees[i]!= null && trees[j]!= null) {
					int rfDistance = getRFDistance(trees[i], trees[j]);
					r[i][j] = rfDistance;
					r[j][i] = rfDistance;

					//					System.out.println(i + " vs " + j + ": " + rfDistance);
					double treeDist = getBranchDistance(trees[i], trees[j]);
					double discDist = getDiscrepantBranchDistance(trees[i], trees[j]);
					System.out.println("Tree Distance: " + i + " vs " + j +
							"\t RF: " + rfDistance + "\tBranch Dist: " +
							treeDist + "\t\t" + discDist);

				}	
			}
		}
		return r;
	}

	/**
	 * Does Robinson-Foulds distance comparison of all pairs
	 * of trees.
	 * Assumes all trees have the same set of leaves.
	 * 
	 * @param trees
	 * @return
	 */
	private int[][] compareAllvsAll(AdvancedTree[] trees) {
		int[][] r = new int[trees.length][trees.length];

		HashSet[] treeBipartSets = new HashSet[trees.length];
		String[] taxonNames = trees[0].getLeafLabels();
		
		System.out.println("Get Bipartitions for trees");
		for(int i = 0; i < treeBipartSets.length; i++) {
			Bipartition[] treeBiparts = trees[i].getBipartitions(taxonNames);
			HashSet hs = new HashSet();
			for(int j = 0; j < treeBiparts.length; j++) {
				hs.add(treeBiparts[j]);
			}
			treeBipartSets[i] = hs;
		}
		System.out.println("done getting Bipartitions");
		
		for(int i = 0; i < r.length; i++) {
			for(int j = i+1; j < r.length; j++) {
				HashSet union = new HashSet(treeBipartSets[i]);
				union.addAll(treeBipartSets[j]);
				
				HashSet intersection = new HashSet(treeBipartSets[i]);
				intersection.retainAll(treeBipartSets[j]);
				
				r[i][j] = (union.size() - intersection.size()) /2;
				System.out.println("Tree Distance: " + i + " vs " + j +
						"\t RF: " + r[i][j]); 
			}
		}
		
		return r;
	}

	/**
	 * Returns a set of Bipartitions defined by the given tree. The index
	 * values for the Bipartitions are based on the sorted list of
	 * leaf names, which is not necessarily the order returned by 
	 * AdvancedTree.getLeafLabels()
	 * 
	 * @param tree
	 * @return
	 */
	private Bipartition[] getTreeBipartitions(AdvancedTree tree) {
		Bipartition[] r = null;
		ArrayList<Bipartition> bipartList = new ArrayList<Bipartition>();
		String[] leaves = tree.getLeafLabels();
		
		//sort the leaves, to allow using binary search to more quickly
		//find the indices.
		Arrays.sort(leaves);
		int[] preorderTraversalSequence = tree.getPreorderTraversalSequence();
		
		for(int i = 0; i < preorderTraversalSequence.length; i++) {
			String[] descendants = 
				tree.getDescendantLeaves(preorderTraversalSequence[i]);
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

	private void compareBranchSupports(AdvancedTree treeA, AdvancedTree treeB) {
		
		String[] tempOG = new String[]{treeA.getLeafLabels()[0]};
		treeA.setOutGroup(tempOG);
		treeB.setOutGroup(tempOG);
		
		int[] treeASupports = treeA.getBranchSupports();
		int[] treeBSupports = treeB.getBranchSupports();

		int[][] branchSupports = new int[treeASupports.length][2];
		for(int i = 0; i < treeASupports.length; i++) {
			String[] treeALeaves = treeA.getDescendantLeaves(i);
			int treeBNode = treeB.getMatchingNode(treeALeaves);
			int bSupport = -1;
			if(treeBNode >= 0) {
				bSupport = treeBSupports[treeBNode];
			}

			branchSupports[i][0] = treeASupports[i];
			branchSupports[i][1] = bSupport;
		}

		ArrayList pairStringList = new ArrayList(branchSupports.length);
		HashMap<String,Integer> pairToCountMap = new HashMap<String,Integer>();

		for(int i = 0; i < branchSupports.length;i++) {
			String pairString = "" + branchSupports[i][0] + ", " +
			branchSupports[i][1];
//		System.out.println(pairString);
		Integer count = pairToCountMap.get(pairString);
			if(count == null) {
				count = new Integer(0);
			}
			count = new Integer(count.intValue() + 1);
			pairToCountMap.put(pairString, count);

		}

		branchSupports = new int[treeBSupports.length][2];
		for(int i = 0; i < treeBSupports.length; i++) {
			String[] treeBLeaves = treeB.getDescendantLeaves(i);
			int treeANode = treeA.getMatchingNode(treeBLeaves);
			int aSupport = -1;
			if(treeANode >= 0) {
				aSupport = treeASupports[treeANode];
			}

			branchSupports[i][0] = aSupport;
			branchSupports[i][1] = treeBSupports[i];
		}

		for(int i = 0; i < branchSupports.length;i++) {
			if(branchSupports[i][0] < 0) {
				String pairString = "" + branchSupports[i][0] + ", " +
    				branchSupports[i][1];
				
				Integer count = pairToCountMap.get(pairString);
				if(count == null) {
					count = new Integer(0);
				}
				count = new Integer(count.intValue() + 1);
				pairToCountMap.put(pairString, count);
			}
		}
		
		String[] distinctPairs = new String[pairToCountMap.size()];
		pairToCountMap.keySet().toArray(distinctPairs);
		Arrays.sort(distinctPairs);

		System.out.print("[");
		for(int i = 0; i < distinctPairs.length; i++) {
			int count = (pairToCountMap.get(distinctPairs[i])).intValue();
			System.out.print("[" + distinctPairs[i] + ", " + count + "]");
			if(i+1 < distinctPairs.length) {
				System.out.print(",");
			}
		}
		System.out.print("]");
		System.out.println();
	}


	private int[] rollingComparison(TreeI[] trees) {
		int[] r = null;
		r = new int[trees.length-1];
		for(int i = 0; i < r.length; i++) {
			r[i] = trees[i+1].getRobinsonFouldsDistance(trees[i]);
			System.out.println("rf between " + i + " and " + (i+1) + ": " 
					+ r[i]);
		}
		return r;
	}

	/**
	 * load trees from a Mr Bayes .t file.
	 * @throws IOException 
	 */
	private void loadComparisonTreesFromMrTFile(String fileName) throws IOException {
		TextFile tf = new TextFile(fileName);
		int firstTreeLine = 0;

		//skip lines until the first line with a tree 
		//tree lines begin with several spaces, the 'tree rep'
		boolean lookingForStart = true;
		while(lookingForStart) {
			if(tf.getLine(firstTreeLine).trim().startsWith("tree rep")) {
				lookingForStart = false;
			} else {
				firstTreeLine++;
			}
		}

		//get the tree String from each line, as long as the lines have trees
		char equalsChar = '=';
		ArrayList treeStringList = 
			new ArrayList(tf.getLineCount()-firstTreeLine);
		for(int i = firstTreeLine; i < tf.getLineCount(); i++) {
			String line = tf.getLine(i);
			int equalsIndex = line.indexOf(equalsChar);
			if(equalsIndex >= 0) {
				line = line.substring(equalsIndex+1).trim();
				treeStringList.add(line);
				System.out.println(line);
			}
		}

		String[] treeStrings = new String[treeStringList.size()];
		treeStringList.toArray(treeStrings);

		System.out.println("TreeComparison.loadComparisonTreesFromMrTFile() " +
				"begin parsing " + treeStrings.length + " trees...");
		comparisonTrees = new TreeI[treeStrings.length];
		for(int i = 0; i < comparisonTrees.length; i++) {
			comparisonTrees[i] = TreeParser.parseTreeString(treeStrings[i]);
			if(i%100 == 0){
				System.out.println();
			}
			System.out.print(".");
		}
		System.out.println("TreeComparison.loadComparisonTreesFromMrTFile() " +
		"done parsing trees");

	}

	private void loadStandardTreeFromFile(String fileName, int row, int col) throws IOException {
		//		System.out.println("TreeComparison.loadStandardTreeFromFile()");
		TextFile tf = new TextFile(fileName);
		String delimiter = "\t";
		String treeString = tf.getLine(row).split(delimiter)[col];
		standardTree = TreeParser.parseTreeString(treeString);
		System.out.println("standardTree leaves: " + standardTree.getLeafCount());
	}

	private void loadComparisonTreesFromFile(String fileName, int col) throws IOException {
		//		System.out.println("TreeComparison.loadComparisonTreesFromFile()");
		TextFile tf = new TextFile(fileName);
		String delimiter = "\t";
		int treeCount = tf.getLineCount()-1; //first line is column headers, not tree
		if(tf.getLine(treeCount-1).length() < 5) {
			treeCount--;
		}
		comparisonTrees = new TreeI[treeCount];
		for(int i = 1; i < comparisonTrees.length; i++) {
			System.out.println("TreeComparison.loadComparisonTreesFromFile() tree " + i);
			String line = tf.getLine(i);
			System.out.println(line);
			String[] fields = line.split(delimiter);
			if(fields.length > col) {
				String treeString = fields[col];
				try 
				{
					TreeI tree = TreeParser.parseTreeString(treeString);
					comparisonTrees[i-1] = tree;
				} catch(NullPointerException npe) {
					System.out.println("problem parsing tree " + i + 
							": " + treeString);
				}
			}
		}

	}

	private void calculateTreeDistances() {
		treeDistancesFromStandard = new int[comparisonTrees.length];
		for(int i = 0; i < treeDistancesFromStandard.length; i++) {
			if(comparisonTrees[i] != null) {
				treeDistancesFromStandard[i] = 
					getRFDistance(standardTree, comparisonTrees[i]);
			}
		}
	}

	static int getRFDistance(TreeI tree1, TreeI tree2) {
		int r = -1;
		r = tree1.getRobinsonFouldsDistance(tree2);
		return r;
	}

	private void printTreeDistances() {
		//		System.out.println("standardTree leaves: " + standardTree.getLeafCount());
		int sum = 0;
		int count = 0;
		int min = Integer.MAX_VALUE;
		int max = 0;
		for(int i = 0; i < treeDistancesFromStandard.length; i++) {
			if(comparisonTrees[i] != null) {
				int distance = treeDistancesFromStandard[i];
				System.out.println("distance between standard tree and tree " +
						i + ": " + distance + " leaves: "
						+ comparisonTrees[i].getLeafCount());
				sum += distance;
				count++;
				min = Math.min(distance, min);
				max = Math.max(distance, max);
			}
		}

		System.out.println("mean: " + ((double)sum/count)
				+ "  min: " + min + "  max: " + max);
	}

	/**
	 * Custom method for calculating a distance between two trees.
	 * Square root of the sum of squares of the differences of the lengths of
	 * the branches. 
	 * 
	 * This is similar to the Branch Score Distance from the Phylip
	 * program treeDist, except the scores here are normalized to between
	 * 0.0 and 1.0 by dividing by the maximum possible distance.
	 * 
	 * @return
	 */
	static double getBranchDistance(TreeI tree1, TreeI tree2) {
		double r = 0;
		boolean normalizeBranchLengths = true;
		double sumOfSquaresTree1 = 0;
		double sumOfSquaresTree2 = 0;
		HashMap tree1SetsToBranches = new HashMap();
		HashMap tree2SetsToBranches = new HashMap();
		HashSet allBranches = new HashSet();
		HashSet checkedBranches = new HashSet();

		TreeBranch[] tree1Branches = tree1.getBranches();
		double[] tree1BranchLengths = new double[tree1Branches.length];
		double[] tree1BranchSupports = new double[tree1Branches.length];
		for(int i= 0; i < tree1Branches.length; i++) {
			tree1BranchLengths[i] = tree1Branches[i].getBranchLength();
			tree1BranchSupports[i] = tree1Branches[i].getBranchSupport();
			HashSet[] leafSets = tree1Branches[i].getBipartitionLeafSets();
			for(int j = 0; j < leafSets.length; j++) {
				tree1SetsToBranches.put(leafSets[j], tree1Branches[i]);
			}
		}

		TreeBranch[] tree2Branches = tree2.getBranches();
		double[] tree2BranchLengths = new double[tree2Branches.length];
		for(int i= 0; i < tree2Branches.length; i++) {
			tree2BranchLengths[i] = tree2Branches[i].getBranchLength();
			HashSet[] leafSets = tree2Branches[i].getBipartitionLeafSets();
			for(int j = 0; j < leafSets.length; j++) {
				tree2SetsToBranches.put(leafSets[j], tree2Branches[i]);
			}
		}

		//normalize branch lengths
		//		if(normalizeBranchLengths) {
		double tree1MaxLength = StatisticsUtilities.getMaximum(tree1BranchLengths);
		//			System.out.println("tree1 maxLength: " + tree1MaxLength);
		for(int i = 0; i < tree1BranchLengths.length; i++) {
			tree1BranchLengths[i] /= tree1MaxLength;
		}

		double tree2MaxLength = StatisticsUtilities.getMaximum(tree2BranchLengths);
		//			System.out.println("tree2 maxLength: " + tree2MaxLength);
		for(int i = 0; i < tree2BranchLengths.length; i++) {
			tree2BranchLengths[i] /= tree2MaxLength;
		}
		//		}

		//check each tree1 branch for a counterpart in tree2. Add square of the 
		//difference in the branch lengths
		ArrayList<Double> removedBranchSupports = new ArrayList<Double>();
		for(int i = 0; i < tree1Branches.length; i++) {
			double counterpartLength = 0;
			TreeBranch tree1Branch = tree1Branches[i];
			HashSet[] leafSets = tree1Branch.getBipartitionLeafSets();
			double tree1Support = tree1Branch.getBranchSupport();
			TreeBranch counterpartBranch = (TreeBranch) tree2SetsToBranches.get(leafSets[1]);
			if(counterpartBranch != null) {
				counterpartLength = counterpartBranch.getBranchLength()/tree2MaxLength ;
				if((tree1Support == tree1Support || tree1BranchSupports[i] == tree1BranchSupports[i]) && tree1Support != tree1BranchSupports[i]) {
					System.out.println("branch supports for branch " + i + " disagree. " + tree1BranchSupports[i] + " -- " + tree1Support);
				}
				double tree2Support = counterpartBranch.getBranchSupport();
				if(tree1Support == tree1Support && tree1Support + tree2Support < 200) {
					HashSet smallerLeafSet = leafSets[0];
					if(leafSets[1].size() < leafSets[0].size()) {
						smallerLeafSet = leafSets[1];
					}

					HashSet[] counterpartLeafSets = counterpartBranch.getBipartitionLeafSets();
					HashSet smallerCounterpartLeafSet = counterpartLeafSets[0];
					if(counterpartLeafSets[1].size() < counterpartLeafSets[0].size()) {
						smallerCounterpartLeafSet = counterpartLeafSets[1];
					}
//					System.out.println(smallerLeafSet);
//					System.out.println(smallerCounterpartLeafSet);
//					System.out.println("counterpart found for tree1Branches " + i + " size: " + smallerLeafSet.size() + ", " 
//						 + tree1Support + " --> " + tree2Support);
				}
			} else {
				System.out.println("Branch removed. Initial support: " + tree1Support);
				removedBranchSupports.add(tree1Support);
			}
			double branchLength = tree1BranchLengths[i];
			sumOfSquaresTree1 += branchLength*branchLength;
			double difference = branchLength - counterpartLength;
			r += difference*difference;

			checkedBranches.add(counterpartBranch);
		}
		double[] removedSupports = new double[removedBranchSupports.size()];
		double sum = 0;
		for(int i = 0; i < removedSupports.length; i++) {
			removedSupports[i] = removedBranchSupports.get(i);
			sum += removedSupports[i];
		}
		double meanRemovedSupport = sum / removedSupports.length;
		Arrays.sort(removedSupports);
		for(int i = 0; i< removedSupports.length; i++) {
			System.out.println(removedSupports[i]);
		}
		double medianRemovedSupport = removedSupports[removedSupports.length/2];
		System.out.println("median support value of removed branches: " + medianRemovedSupport);
		System.out.println("mean support value of removed branches: " + meanRemovedSupport);
		
		ArrayList<Double> tree1NonNanSupports = new ArrayList<Double>();
		sum = 0;
		for(int i = 0; i < tree1BranchSupports.length; i++) {
			if(tree1BranchSupports[i] == tree1BranchSupports[i]) {
				tree1NonNanSupports.add(tree1BranchSupports[i]);
				sum += tree1BranchSupports[i];
			}
		}
		Collections.sort(tree1NonNanSupports);
		double medianAll = tree1NonNanSupports.get(tree1NonNanSupports.size()/2);
		double meanAll = sum / tree1NonNanSupports.size();
		System.out.println("median support value of all branches in tree1: " + medianAll);
		System.out.println("mean support value of all branches in tree1: " + meanAll);

		//add any tree2 branches that have not already been handled. Because 
		//these have not been handled, it means they don't have a counterpart
		//in tree1. So add the squares of these branch lengths, with no need
		//to determine a difference with counterpart branches
		for(int i = 0; i < tree2Branches.length; i++) {
			double branchLength = tree2BranchLengths[i];
			sumOfSquaresTree2 += branchLength*branchLength;
			if(!checkedBranches.contains(tree2Branches[i])) {
				r += branchLength*branchLength;
			}
		}

		r = Math.sqrt(r) / Math.sqrt(sumOfSquaresTree1 + sumOfSquaresTree2);
		return r;
	}

	/**
	 * Custom method for calculating a distance between two trees.
	 * Square root of the sum of squares of the differences of the lengths of
	 * the branches. This is the same as getBranchDistance() excpet it 
	 * only includes branches found in one tree and not the other..
	 * @return
	 */
	private double getDiscrepantBranchDistance(TreeI tree1, TreeI tree2) {
		double r = 0;
		double sumOfSquaresTree1 = 0;
		double sumOfSquaresTree2 = 0;
		HashMap tree1SetsToBranches = new HashMap();
		HashMap tree2SetsToBranches = new HashMap();
		HashSet allBranches = new HashSet();
		HashSet checkedBranches = new HashSet();

		TreeBranch[] tree1Branches = tree1.getBranches();
		for(int i= 0; i < tree1Branches.length; i++) {
			HashSet[] leafSets = tree1Branches[i].getBipartitionLeafSets();
			for(int j = 0; j < leafSets.length; j++) {
				tree1SetsToBranches.put(leafSets[j], tree1Branches[i]);
			}
		}

		TreeBranch[] tree2Branches = tree2.getBranches();
		for(int i= 0; i < tree2Branches.length; i++) {
			HashSet[] leafSets = tree2Branches[i].getBipartitionLeafSets();
			for(int j = 0; j < leafSets.length; j++) {
				tree2SetsToBranches.put(leafSets[j], tree2Branches[i]);
			}
		}

		//check each tree1 branch for a counterpart in tree2. Add square of the 
		//difference in the branch lengths
		for(int i = 0; i < tree1Branches.length; i++) {
			double counterpartLength = 0;
			HashSet[] leafSets = tree1Branches[i].getBipartitionLeafSets();
			double branchLength = tree1Branches[i].getBranchLength();
			sumOfSquaresTree1 += branchLength*branchLength;
			TreeBranch counterpart = (TreeBranch) tree2SetsToBranches.get(leafSets[0]);
			if(counterpart == null) {
				double difference = branchLength - counterpartLength;
				r += difference*difference;
			}			
			checkedBranches.add(counterpart);
		}

		//add any tree2 branches that have not already been handled. Because 
		//these have not been handled, it means they don't have a counterpart
		//in tree1. So add the squares of these branch lengths, with no need
		//to determine a difference with counterpart branches
		for(int i = 0; i < tree2Branches.length; i++) {
			double branchLength = tree2Branches[i].getBranchLength();
			sumOfSquaresTree2 += branchLength*branchLength;
			if(!checkedBranches.contains(tree2Branches[i])) {
				r += branchLength*branchLength;
			}
		}

		r = Math.sqrt(r) / Math.sqrt(sumOfSquaresTree1 + sumOfSquaresTree2);
		return r;
	}
	
	public String[] runConsel(SequenceAlignment alignment, BasicTree[] trees, 
			String processors, String matrix) {
		String[] r = null;
		try {
			File workingDir = new File(System.getProperty("user.dir"));
			//Write alignment to a file in phylip format
			String alignmentFileName = alignment.getName() + ".phy";
			FileWriter alignmentWriter = new FileWriter(alignmentFileName);
			alignmentWriter.write(alignment.getAlignmentAsExtendedPhylipUsingSequenceNames()+ "\n");
			alignmentWriter.flush();
			alignmentWriter.close();

			//Write trees to file 
			File treeFile = File.createTempFile("trees", ".nwk", workingDir);
			String treeFileName = treeFile.getPath();
			FileWriter treeFileWriter = new FileWriter(treeFile);
			for(int i = 0; i < trees.length; i++) {
				//				System.out.println("trees " + i + ": " + trees[i]);
				treeFileWriter.write(trees[i].getTreeString(true, true) + "\n");
			}
			treeFileWriter.flush();
			treeFileWriter.close();

			//Generate a temporary file name to use for raxml run name
			File raxmlRunFile = File.createTempFile("run", "", workingDir);
			String raxmlRunName = raxmlRunFile.getName();
			raxmlRunFile.delete();

			//run consel with the alignment and the trees
			//run raxml to get per-site log likelihoods
			String raxmlPath = ExecUtilities.getCommandPath("raxmlHPC-PTHREADS");
			String raxmlCommand = raxmlPath + " -f g -m "
			+ matrix + " -T " + processors + " -z " + treeFileName 
			+ " -s " + alignmentFileName + " -n " + raxmlRunName;
			System.out.println(raxmlCommand);
			ExecUtilities.exec(raxmlCommand);

			//Raxml result file needs to be renamed so the makermt result files
			//will have unique names.
			String raxmlOutFileName = "RAxML_perSiteLLs." + raxmlRunName;
			String renamedRaxmlOutFile = alignmentFileName + "_perSiteLLs." +
			raxmlRunName;
			File raxmlOut = new File(raxmlOutFileName);
			raxmlOut.renameTo(new File(renamedRaxmlOutFile));

			//run makermt on per site likelihood file
			String makermtPath = ExecUtilities.getCommandPath("makermt");
			String makermtCommnd = makermtPath + 
			" -b 10" + //-b 10 multiplies the number of bootstraps by 10, for less error
			" --puzzle " + renamedRaxmlOutFile;
			System.out.println(makermtCommnd);
			ExecUtilities.exec(makermtCommnd);			
			String makermtOutFileName = alignmentFileName + "_perSiteLLs";

			//run consel on makermt result file
			String conselPath = ExecUtilities.getCommandPath("consel");
			String conselCommand = conselPath + " " + makermtOutFileName;
			System.out.println(conselCommand);
			ExecUtilities.exec(conselCommand);

			//run catpv on consel result file
			String catpvPath = ExecUtilities.getCommandPath("catpv");
			String catpvCommand = catpvPath + " -v " + makermtOutFileName;
			System.out.println(catpvCommand);
			CommandResults results = ExecUtilities.exec(catpvCommand);

			String[] resultTable = results.getStdout();
			r = resultTable;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return r;
	}

	private String[] getDistinctTaxa(FastaSequenceFile[] inputSequenceFiles) {
		HashSet uniqueTaxa = new HashSet();

		for(int i = 0; i < inputSequenceFiles.length; i++) {
			String[] taxa = inputSequenceFiles[i].getTaxa();
			for(int j = 0; j < taxa.length; j++) {
				uniqueTaxa.add(taxa[j]);
			}
		}

		String[] taxa = new String[uniqueTaxa.size()];
		uniqueTaxa.toArray(taxa);
		return taxa;
	}

	private FastaSequenceFile[] loadSequenceFiles(String[] sequenceFileNames) 
	throws IOException {
		FastaSequenceFile[] r = new FastaSequenceFile[sequenceFileNames.length];
		for(int i = 0; i < r.length; i++) {
			r[i] = new FastaSequenceFile(sequenceFileNames[i]);
		}
		return r;
	}
}
