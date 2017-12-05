package edu.vt.vbi.ci.pepr.tree;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

import org.apache.log4j.Logger;

import edu.vt.vbi.ci.pepr.tree.pipeline.PhyloPipeline;
import edu.vt.vbi.ci.util.CommandLineProperties;
import edu.vt.vbi.ci.util.HandyConstants;
import edu.vt.vbi.ci.util.PEPRTracker;
import edu.vt.vbi.ci.util.file.FastaSequenceFile;
import edu.vt.vbi.ci.util.file.TextFile;

public class PhylogeneticTreeRefiner {

	private HashMap taxonToFileName;
	private CommandLineProperties clp;
	private HashSet refinedSubsets;
	private String runName;
	private int refineRound = 0;
	private static Logger logger = Logger.getLogger(PhylogeneticTreeRefiner.class);
	private String mostRefinedTreeString;


	public static void main(String[] args) {
		CommandLineProperties clp = new CommandLineProperties(args);

		String treeFileName = clp.getValues(HandyConstants.TREE_FILE)[0];
		String[] inputSequenceFileNames = 
			clp.getValues(HandyConstants.GENOME_FILE);
		if(inputSequenceFileNames == null) {
			logger.error("No input sequence files were provided. " +
					"Please provide sequence files with the command -"
					+ HandyConstants.GENOME_FILE);
		}

		String[] outgroupFileNames = clp.getValues(HandyConstants.OUTGROUP);
		FastaSequenceFile[] inputSequenceFiles = null;
		FastaSequenceFile[] outgroupSequenceFiles = null;
		String initialTree = null;

		//See if any 'track' is specified, and load the relevant properties.
		//Tracks can be a starting point and be modified by other properties.
		String track = clp.getValues(HandyConstants.TRACK, HandyConstants.FALSE)[0];
		if(!track.equals(HandyConstants.FALSE)) {
			String[] trackProperties = PhyloPipeline.getTrackProperties(track);
			clp.addArgs(trackProperties);
		}

		try {
			inputSequenceFiles = loadSequenceFiles(inputSequenceFileNames);
			if(outgroupFileNames != null) {
				outgroupSequenceFiles = loadSequenceFiles(outgroupFileNames);
			} else {
				outgroupSequenceFiles = new FastaSequenceFile[0];
			}
			
			TextFile treeFile = new TextFile(treeFileName);
			initialTree = treeFile.toString();
			PhylogeneticTreeRefiner refiner = 
				new PhylogeneticTreeRefiner(initialTree, clp, 
						inputSequenceFiles, outgroupSequenceFiles);
		} catch (IOException e) {
			logger.error(e);
		}
	}

	public PhylogeneticTreeRefiner(String initialTree, 
			CommandLineProperties clp, FastaSequenceFile[] ingroupSequences, 
			FastaSequenceFile[] outgroupSequences) {
		this.clp = clp;

		runName = clp.getValues(HandyConstants.RUN_NAME, runName)[0];
		refine(initialTree, clp, ingroupSequences, outgroupSequences);
	}

	private void refine(String initialTree, 
			CommandLineProperties clp, FastaSequenceFile[] ingroupSequences, 
			FastaSequenceFile[] outgroupSequences) {
		logger.info("PhylogeneticTreeRefiner.refine()");
		int threadCount = Integer.parseInt(clp.getValues(HandyConstants.MAX_CONCURRENT_PROCESS_PARAM, "1")[0]);
		logger.info("initialTree string: " + initialTree.trim());

		setMostRefinedTreeString(initialTree);
		int requiredSpeciesForUniqueSpeciesFilter = 5;
		
		//create taxonToFileName map
		taxonToFileName = new HashMap();

		refinedSubsets = new HashSet();

		for(int i = 0; i < ingroupSequences.length; i++) {
			String fileName = ingroupSequences[i].getFullName();
			String[] taxa = ingroupSequences[i].getTaxa();
			for(int j = 0; j < taxa.length; j++) {
				taxonToFileName.put(taxa[j], fileName);
			}
		}

		HashSet outgroupTaxa = new HashSet();
		for(int i = 0; i < outgroupSequences.length; i++) {
			String fileName = outgroupSequences[i].getFullName();
			String[] taxa = outgroupSequences[i].getTaxa();
			for(int j = 0; j < taxa.length; j++) {
				taxonToFileName.put(taxa[j], fileName);
				outgroupTaxa.add(taxa[j]);
			}
		}

		String[] outgroupList = new String[outgroupTaxa.size()];
		outgroupTaxa.toArray(outgroupList);
		AdvancedTree tree = new AdvancedTree(initialTree);
		tree.setOutGroup(outgroupList);

		int refineCutoff = 
			Integer.parseInt(clp.getValues(HandyConstants.REFINE_CUTOFF, "100")[0]);

		AdvancedTree refiningTree = tree;

		//get index for next node to be refined
		int nextRefineIndex = getNextIndexToRefine(refiningTree, refineCutoff);
		logger.info("nextRefineIndex: " + nextRefineIndex);

		int refiningRound = 0;
		while(nextRefineIndex > -1) {
			refiningRound++;
			logger.info("tree refining round: " + refiningRound);
			//get ingroup sequences for subtree
			String[] descendants = refiningTree.getDescendantLeaves(nextRefineIndex);

			HashSet refineIngroup = new HashSet();
			String[] genomeFiles = new String[descendants.length];
			logger.info("next node to refine is " + nextRefineIndex + 
					" and it has " + descendants.length + " descendants:");
			for(int i = 0; i < descendants.length; i++) {
				refineIngroup.add(descendants[i]);
				genomeFiles[i] = (String)taxonToFileName.get(descendants[i]);
				logger.info("\t" + descendants[i] + "\t" + genomeFiles[i]);
			}

			int uniqueSpeciesCount = getUniqueSpeciesCountForGenomeFiles(genomeFiles);
			String useUniqueSpeciesFilter = HandyConstants.FALSE;
			if(uniqueSpeciesCount >= requiredSpeciesForUniqueSpeciesFilter) {
				useUniqueSpeciesFilter = HandyConstants.TRUE;
			}

			//get outgroup sequence for subtree
			HashSet outgroupSet = new HashSet();
			String[] inAndOut = null;
			int parentNode = refiningTree.getParentNode(nextRefineIndex);
			if(parentNode == -1) {
				logger.info("refine() parentNode is " +
				"-1, so there is no outgroup");
				inAndOut = refiningTree.getLeafLabels();
			} else {
				inAndOut = refiningTree.getDescendantLeaves(parentNode);
				logger.info("the parent of the nextRefineNode is " +
						parentNode + " and it has " + (inAndOut.length-descendants.length) + 
				" additional descendants:");
				for(int i = 0; i < inAndOut.length; i++) {
					if(!refineIngroup.contains(inAndOut[i])) {
						outgroupSet.add(inAndOut[i]);
						logger.info(inAndOut[i]);
					}
				}
			}

			String[] refineOutgroupList = new String[outgroupSet.size()];
			outgroupSet.toArray(refineOutgroupList);
			String[] outgroupFiles = new String[refineOutgroupList.length];
//			String[] refineOutgroupFiles = new String[refineOutgroupList.length];
			logger.info("outgroup pool:");
			for(int i = 0; i < refineOutgroupList.length; i++) {
				outgroupFiles[i] = (String) taxonToFileName.get(refineOutgroupList[i]);
				logger.info("\t" + refineOutgroupList[i] + "\t" + outgroupFiles[i]);
			}

			//build subtree
			//assemble commands for PhylogenomicPipeline
			CommandLineProperties pipelineCLP = new CommandLineProperties(clp.getArgs());
			//remove values that don't apply to the refinement run
			pipelineCLP.remove(HandyConstants.GENOME_FILE);
			pipelineCLP.remove(HandyConstants.OUTGROUP);
			pipelineCLP.remove(HandyConstants.RUN_NAME);
			pipelineCLP.remove(HandyConstants.REFINE);
			pipelineCLP.remove(HandyConstants.UNIQUE_SPECIES);
			pipelineCLP.remove(HandyConstants.MAX_TAXA);
			pipelineCLP.remove(HandyConstants.MIN_TAXA);
			pipelineCLP.remove(HandyConstants.TARGET_NTAX);
			pipelineCLP.remove(HandyConstants.OUTGROUP_COUNT);

			String[] genomeFileArgs = new String[genomeFiles.length+1];
			genomeFileArgs[0] = "-" + HandyConstants.GENOME_FILE;
			System.arraycopy(genomeFiles, 0, genomeFileArgs, 1, genomeFiles.length);
			pipelineCLP.addArgs(genomeFileArgs);

			String[] outgroupFileArgs = new String[outgroupFiles.length+1];
			outgroupFileArgs[0] = "-" + HandyConstants.OUTGROUP;
			System.arraycopy(outgroupFiles, 0, outgroupFileArgs, 1, outgroupFiles.length);
			pipelineCLP.addArgs(outgroupFileArgs);
			
			//if a track was specified, remove it. Any track parameters
			//will have been added, so there's no need to include the track.
			//Including the track may end up overriding other parameters.
			pipelineCLP.remove(HandyConstants.TRACK);
			
			logger.info("PhylogeneticTreeRefiner run pipeline with unique_species set to "
					+ useUniqueSpeciesFilter);
			int outgroupCount = Math.min(outgroupFiles.length, 2);
			//additional pipeline args
			pipelineCLP.addArgs(
					new String[]{
							"-" + HandyConstants.OUTGROUP_COUNT, 
							"" + outgroupCount,
							"-" + HandyConstants.UNIQUE_SPECIES,
							useUniqueSpeciesFilter,
							"-" + HandyConstants.RUN_NAME,
							runName + "_refine_sub" + refiningRound,
							"-" + HandyConstants.REFINE,
							HandyConstants.FALSE,
							"-" + HandyConstants.SUBTREE,
							HandyConstants.TRUE
					}
			);

			PhyloPipeline pipeline = new PhyloPipeline(pipelineCLP.getArgs());
			String refinedTreeString = pipeline.getTree();
			AdvancedTree refinedSubTree = new AdvancedTree(refinedTreeString);
			refinedSubTree.setOutGroup(refineOutgroupList);
			String rootedSubtreeString = refinedSubTree.getTreeString(true, true);
			PEPRTracker.setTree(rootedSubtreeString);
			
			logger.info("refined subtree: " + refinedTreeString);

			//replace subtree in the full tree
			logger.info("replacing tree node " + nextRefineIndex + 
			" with the new subtree");
			logger.info("refiningTree leaf count: " + refiningTree.getLeafLabels().length);
			logger.info("refiningTree before replacing node: " 
					+ refiningTree.getTreeString(true, true));

			refiningTree = refiningTree.replaceNode(refinedSubTree.getBasicTree());
			refiningTree.unroot();
			refiningTree.setOutGroup(outgroupList);

			String refinedFullTreeString = refiningTree.getTreeString(true, true);
			setMostRefinedTreeString(refinedFullTreeString);
			PEPRTracker.setFullTree(refinedFullTreeString);
			
			logger.info("refiningTree with node " + nextRefineIndex + 
					" refined: ");
			logger.info(refinedFullTreeString);
			logger.info("refiningTree leaf count: " + 
					refiningTree.getLeafLabels().length);

			//write refined tree to file
			try {
				FileWriter fw = new FileWriter(runName + "_refine_" + refiningRound + ".nwk");
				fw.write(refiningTree.getTreeString(true, true));
				fw.flush();
				fw.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

			refiningTree.setOutGroup(outgroupList);
			nextRefineIndex = getNextIndexToRefine(refiningTree, refineCutoff);
			logger.info("nextRefineIndex: " + nextRefineIndex);
		}
	}

	private int getUniqueSpeciesCountForGenomeFiles(String[] fileNames) {
		int r = -1;
		//create FastaSequenceFiles 
		FastaSequenceFile[] sequenceFiles = new FastaSequenceFile[fileNames.length];
		try {
			for(int i = 0; i < sequenceFiles.length; i++) {
				sequenceFiles[i] = new FastaSequenceFile(fileNames[i]);
			}
		} catch(IOException ioe) {
			ioe.printStackTrace();
		}

		//use PhyloPipeline method to filter for unique species
		FastaSequenceFile[] filtered = 
			PhyloPipeline.filterOutDuplicateSpecies(sequenceFiles);
		
		//return the number of filtered files
		r = filtered.length;
		return r;
	}

	private int getNextIndexToRefine(AdvancedTree tree, int refineCutoff) {
		int r = -1;
		
		int minLeavesToRefine = 3; //don't try refining nodes with fewer than this many leaves
		int[] preorderTraversalSequence = tree.getPreorderTraversalSequence();
		int[] meanSupports = tree.getMeanDescendantSupportValues();
		int[] branchSupports = tree.getBranchSupports();
		boolean[] alreadyRefined = new boolean[meanSupports.length];

		//start at 2 instead of 0 because 0 is the root and the parent of 1
		//is the root, which would end up trying to refine the full tree
		//and there's no point in trying to refine the full tree.
		//preorderTraversalSequence[2] is the first node that can be
		//refined without trying to refine the full tree
		for(int i = 2; r == -1 && i < preorderTraversalSequence.length; i++) {
			r = preorderTraversalSequence[i];
			int descendantLeaves = tree.getDescendantLeaves(r).length;
			
			if(!alreadyRefined[r] &&
					meanSupports[r] < refineCutoff && 
					branchSupports[r] >= refineCutoff &&
					descendantLeaves >= minLeavesToRefine) {
				//see if all child nodes have at least refineCutoff support. If they do, then skip
				//this node.
				int[] childNodes = tree.getChildNodes(r);
				boolean allChildrenFullSupport = true;
				for(int j = 0; allChildrenFullSupport && j < childNodes.length; j++) {
					if(branchSupports[childNodes[j]] < refineCutoff) {
						allChildrenFullSupport = false;
					}
				}

				if(allChildrenFullSupport) {
					r = -1;
				} else {
					HashSet ingroup = new HashSet();
					String[] descendants = tree.getDescendantLeaves(r);
					ingroup = new HashSet();
					String[] genomeFiles = new String[descendants.length];
					for(int j = 0; j < descendants.length; j++) {
						ingroup.add(descendants[j]);
						genomeFiles[j] = (String)taxonToFileName.get(descendants[j]);
					}
					if(refinedSubsets.contains(ingroup)) {
						//this group has already been refined
						r = -1;
					} else if(ingroup.size() <= 2) {
						refinedSubsets.add(ingroup);
						alreadyRefined[r] = true;
						r = -1;
					} else {
						refinedSubsets.add(ingroup);
						alreadyRefined[r] = true;
					}

				}
			} else {
				r = -1;
			}
		}
		return r;
	}

	private static FastaSequenceFile[] loadSequenceFiles(String[] sequenceFileNames) throws IOException {
		FastaSequenceFile[] r = new FastaSequenceFile[sequenceFileNames.length];
		for(int i = 0; i < r.length; i++) {
			r[i] = new FastaSequenceFile(sequenceFileNames[i]);
		}
		return r;
	}

	public String getMostRefinedTreeString() {
		return mostRefinedTreeString;
	}

	private void setMostRefinedTreeString(String mostRefinedTreeString) {
		this.mostRefinedTreeString = mostRefinedTreeString;
		PEPRTracker.setFullTree(mostRefinedTreeString);
	}

}
