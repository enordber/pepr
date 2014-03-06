package edu.vt.vbi.ci.pepr.tree;

import java.io.IOException;
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
			String mainTree = new TextFile(mainTreeFileName).toString().trim();
			
			String[] supportTrees = new String[supportTreeFileNames.length];
			for(int i = 0; i < supportTrees.length; i++) {
				supportTrees[i] = new TextFile(supportTreeFileNames[i]).toString().trim();
			}
			
			String decoratedTree = addSupportValues(mainTree, supportTrees);
			System.out.println(decoratedTree);
		} catch (IOException e) {
			// TODO Auto-generated catch block
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
//			if(mainTreeBiparts[i].getSmallerSide().cardinality() >= 1) {
//				System.out.println(mainTreeBiparts[i].getString()
//						+ "\t" + bipartSupportCounts[i] + " of " + supportTrees.length);
//			}
		}

		//add support values to the main tree
		mainTree.setBranchSupportValues(bipartSupportCounts);

		//get String of main tree with support values added
		r = mainTree.getTreeString(true, true);

		return r;
	}

}
