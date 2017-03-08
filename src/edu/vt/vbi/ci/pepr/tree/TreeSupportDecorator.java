package edu.vt.vbi.ci.pepr.tree;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

import edu.vt.vbi.ci.util.CommandLineProperties;
import edu.vt.vbi.ci.util.ExtendedBitSet;
import edu.vt.vbi.ci.util.file.TextFile;

public class TreeSupportDecorator {

	private static final String MAIN = "main";
	private static final String SUPPORT = "support";
	private static final String NONE = "none";



	public static void main(String[] args) {
		CommandLineProperties clp = new CommandLineProperties(args);
		String mainTreeFileName = clp.getValues(MAIN, NONE)[0];
		String[] supportTreeFileNames = clp.getValues(SUPPORT, NONE);

		//load tree strings from files
		try {
			System.out.println("main tree file: " + mainTreeFileName);
			String mainTree = new TextFile(mainTreeFileName).toString().trim();

			ArrayList<String> supportTrees = new ArrayList<String>();
			for(int i = 0; i < supportTreeFileNames.length; i++) {
				TextFile supportTreeFile = new TextFile(supportTreeFileNames[i]);
				String[] treeLines = supportTreeFile.getAllLines();
				for(String treeLine: treeLines) {
					if(treeLine.trim().length() > 0) {
						supportTrees.add(treeLine.trim());
					}
				}
			}

			String[] supports = new String[supportTrees.size()];
			supportTrees.toArray(supports);
			String decoratedTree = addSupportValues(mainTree, supports);
			System.out.println(decoratedTree);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	/**
	 * Adds support values to the mainTree based on the presence of
	 * bipartitions in the supportTrees. All tree should have the 
	 * same set of taxa.
	 * 
	 * @param mainTree
	 * @param supportTrees
	 * @return
	 */
	public static TreeI calculateSupports(TreeI mainTree, TreeI[] supportTrees) {
		TreeI r = null;

		TreeNodeI[] leaves = mainTree.getLeafNodes();
		String[] taxonNames = new String[leaves.length];
		for(int i = 0; i < taxonNames.length; i++) {
			taxonNames[i] = leaves[i].getNodeLabel(null);
		}

		//load taxon names into a HashMap. This is needed by the 
		//Tree.getBipartitions() method.
		HashMap taxaMap = new HashMap();
		for(int i = 0; i < taxonNames.length; i++) {
			taxaMap.put(taxonNames[i], taxonNames[i]);
		}

		//get the bipartitions that are present in the main tree

		//go through each support tree and count determine which main tree 
		//biparts are present in the support tree

		HashSet mainTreeBiparts = mainTree.getBipartitionTaxa(taxaMap);

		return r;
	}

	public static String addSupportValues(String main, String[] supports) {
		String r = null;

		//create BasicTrees for all tree strings
		AdvancedTree mainTree = new AdvancedTree(main);
		mainTree.unroot();

		System.out.println("support trees: " + supports.length);
		AdvancedTree[] supportTrees = new AdvancedTree[supports.length];
		for(int i = 0; i < supportTrees.length; i++) {
			supportTrees[i] = new AdvancedTree(supports[i]);
			supportTrees[i].unroot();
		}

		//create a list of taxa contained in the main tree. The index in this
		//list will be used for the bipartition BitSets
		String[] taxa = mainTree.getLeafLabels();
		Arrays.sort(taxa);

		//Get all Bipartitions from all trees
		//main tree bipartitions are stored in an array with the same indices
		//as the corresponding nodes in the BasicTree
		int nodeCount = mainTree.getNodeCount();
		Bipartition[] mainTreeBiparts = new Bipartition[nodeCount];
		int taxonCount = taxa.length;
		String[][] nodeLeaves = new String[nodeCount][];
		for(int i = 0; i < nodeCount; i++) {
			nodeLeaves[i] = mainTree.getDescendantLeaves(i);
			ExtendedBitSet ebs = new ExtendedBitSet();
			for(int j = 0; j < nodeLeaves[i].length; j++) {
				int index = Arrays.binarySearch(taxa, nodeLeaves[i][j]);
				if(index > -1) {
					ebs.set(index);
				}
			}
			mainTreeBiparts[i] = new Bipartition(ebs, taxonCount);
		}
		//done getting mainTree bipartitions

		//get support tree bipartitions as a BipartitionSet
		BipartitionSet supportTreeBipartSet = new BipartitionSet();
		for(int i = 0; i < supportTrees.length; i++) {
			nodeCount = supportTrees[i].getNodeCount();
			nodeLeaves = new String[nodeCount][];
			for(int j = 0; j < nodeCount; j++) {
				nodeLeaves[j] = supportTrees[i].getDescendantLeaves(j);
				ExtendedBitSet ebs = new ExtendedBitSet();
				for(int k = 0; k < nodeLeaves[j].length; k++) {
					int index = Arrays.binarySearch(taxa, nodeLeaves[j][k]);
					if(index > -1) {
						ebs.set(index);
					}
				}
				Bipartition bipart = new Bipartition(ebs, taxonCount);
				supportTreeBipartSet.add(bipart);
			}	
		}

		//For each bipartiion in the main tree, find out how many support
		//trees contain the same bipartition
		int[] bipartSupportCounts = new int[mainTreeBiparts.length];
		for(int i = 0; i < bipartSupportCounts.length; i++) {
			bipartSupportCounts[i] = supportTreeBipartSet.getCount(mainTreeBiparts[i]);
			if(bipartSupportCounts[i] < 100) {
				System.out.println("main tree bipartition:");
				System.out.println(mainTreeBiparts[i].getString());
				System.out.println("support: " + bipartSupportCounts[i]);
			}
		}

		//add support values to the main tree
		mainTree.setBranchSupportValues(bipartSupportCounts);

		//get String of main tree with support values added
		r = mainTree.getTreeString(true, true);

		return r;
	}

}
