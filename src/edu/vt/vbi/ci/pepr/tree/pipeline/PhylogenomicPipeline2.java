package edu.vt.vbi.ci.pepr.tree.pipeline;

import java.io.File;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;

import org.apache.log4j.Logger;

import edu.vt.vbi.ci.pepr.alignment.AlignmentUtilities;
import edu.vt.vbi.ci.pepr.alignment.ConcatenatedSequenceAlignment;
import edu.vt.vbi.ci.pepr.alignment.MSAConcatenator;
import edu.vt.vbi.ci.pepr.alignment.MSATrimmer;
import edu.vt.vbi.ci.pepr.alignment.MultipleSequenceAligner;
import edu.vt.vbi.ci.pepr.alignment.SequenceAlignment;
import edu.vt.vbi.ci.pepr.alignment.SequenceAlignmentParser;
import edu.vt.vbi.ci.pepr.stats.StatisticsUtilities;
import edu.vt.vbi.ci.pepr.tree.BasicTree;
import edu.vt.vbi.ci.pepr.tree.Bipartition;
import edu.vt.vbi.ci.pepr.tree.BipartitionSet;
import edu.vt.vbi.ci.pepr.tree.FastTreeRunner;
import edu.vt.vbi.ci.pepr.tree.RAxMLRunner;
import edu.vt.vbi.ci.pepr.tree.TreeBuilder;
import edu.vt.vbi.ci.pepr.tree.TreeBuilderComparator;
import edu.vt.vbi.ci.pepr.tree.TreeSupportDecorator;
import edu.vt.vbi.ci.util.CommandLineProperties;
import edu.vt.vbi.ci.util.HandyConstants;
import edu.vt.vbi.ci.util.PEPRTracker;
import edu.vt.vbi.ci.util.RandomSetUtils;
import edu.vt.vbi.ci.util.file.FastaSequenceSet;
import edu.vt.vbi.ci.util.file.FastaUtilities;
import edu.vt.vbi.ci.util.file.SequenceSetProvider;


/**
 * Rewrite/refactoring of PhylogenomicPipeline. The original class
 * has gotten to be bloated and clunky and full of junk.
 * @author enordber
 *
 */
public class PhylogenomicPipeline2 {

	private static HashMap commands;
	private Logger logger = Logger.getLogger(getClass());

	private CommandLineProperties commandLineProperties;

	private SequenceSetProvider inputSequenceSetProvider;
	/*
	 * cache for alignments, to avoid duplicating work
	 */
	private HashMap sequenceSetToAlignment;

	private MultipleSequenceAligner multipleSequenceAligner;
	private MSATrimmer msaTrimmer;
	private boolean gblocksTrim;
	private boolean uniformTrim;

	/*
	 * smallest number of taxa for a sequence set to be included
	 */
	private int minTaxa;

	/*
	 * largest number of taxa for a sequence set to be included
	 */
	private int maxTaxa;

	/*
	 * verbose indicates if a lot of information should be printed out while
	 * the pipeline is running.
	 */
	private boolean verbose = true;
	private int verboseLevel = 2;

	private String finalTree;

	//optional constraint tree to be used for tree building steps
	private String constraintTree;

	private String runName;

	//This value can be set by a command line argument, or by checking the 
	//system (at least on linux). By default, it is large enough 
	//that it will not override the specified number of tree threads
	//unless it is set.
	private long availableRAMForTrees = Long.MAX_VALUE;

	static{
		defineHelp();
	}


	public static void main(String[] args) {
		PhylogenomicPipeline2 pp = new PhylogenomicPipeline2(args);
	}

	public PhylogenomicPipeline2(String[] args) {
		commandLineProperties = new CommandLineProperties(args);

		boolean printHelp = args == null || args.length == 0 ||
				commandLineProperties.getValues(HandyConstants.HELP, 
						HandyConstants.FALSE)[0].equals(HandyConstants.TRUE) ||
						commandLineProperties.getValues("h", 
								HandyConstants.FALSE)[0].equals(HandyConstants.TRUE) ||
								commandLineProperties.getValues("?", 
										HandyConstants.FALSE)[0].equals(HandyConstants.TRUE);
		if(printHelp) {
			printHelp();
			System.exit(0);
		}

		String availableMemory = commandLineProperties.getValues(HandyConstants.AVAILABLE_RAM, null)[0];

		if(availableMemory != null) {
			setAvailableRAMForTrees(availableMemory);
		}

		constraintTree = 
				commandLineProperties.getValues(HandyConstants.CONSTRAINT_TREE, null)[0];
		runName = commandLineProperties.getValues(HandyConstants.RUN_NAME, 
				"PhylogenomicPipeline_" + System.currentTimeMillis())[0];

		//load sequence files
		int minTaxa = Integer.parseInt(commandLineProperties.
				getValues(HandyConstants.MIN_TAXA, "0")[0]); 
		if(verboseLevel > 0) {
			System.out.println("min taxa: " + minTaxa);
		}

		int maxTaxa = Integer.parseInt(commandLineProperties.
				getValues(HandyConstants.MAX_TAXA, "100000000")[0]); 
		if(verboseLevel > 0) {
			System.out.println("maxTaxa: " + maxTaxa);
		}

		//targetTaxa defaults to the value for maxTaxa.
		int targetTaxa = Integer.parseInt(commandLineProperties.
				getValues(HandyConstants.TARGET_NTAX, ""+ maxTaxa)[0]); 
		if(verboseLevel > 0) {
			System.out.println("targetTaxa: " + targetTaxa);
		}

		int targetMinGeneCount = Integer.parseInt(commandLineProperties.
				getValues(HandyConstants.TARGET_MIN_GENE_COUNT, "0")[0]); 
		if(verboseLevel > 0) {
			System.out.println("targetMinGeneCount: " + targetMinGeneCount);
		}

		String[] sequenceFileNames = 
				commandLineProperties.getValues(HandyConstants.FILE);
		boolean inputsAligned = false;
		if(sequenceFileNames == null || sequenceFileNames.length == 0) {
			sequenceFileNames = 
					commandLineProperties.getValues(HandyConstants.ALIGN_FILE);
			if(sequenceFileNames != null && sequenceFileNames.length > 0) {
				//use aligned sequence families as input.
				System.out.println("loading pre-aligned input sequence sets");
				inputsAligned = true;
			}
		}

		if(sequenceFileNames == null) {
			//see if a directory name was provided instead of individual
			//file names. If there is a directory name, the load all
			//files in the directory that end with the .faa extension
			String directoryName = 
					commandLineProperties.getValues(HandyConstants.DIRECTORY, 
							null)[0];
			if(directoryName != null) {
				File directory = new File(directoryName);
				sequenceFileNames = directory.list(new FaaFilter());

				//add the directory name to the file name
				for(int i = 0; i < sequenceFileNames.length; i++) {
					sequenceFileNames[i] = 
							directoryName + File.separator + sequenceFileNames[i];
				}
			}
		}

		//filter sequence files based on being representative
		boolean repOnly = commandLineProperties.getValues(HandyConstants.REPRESETATIVE_ONLY,
				HandyConstants.FALSE)[0].equalsIgnoreCase(HandyConstants.TRUE);

		//how many genes should be filtered in the congruence filtering step
		double percentileToRemove = 
				Double.parseDouble(commandLineProperties.getValues(HandyConstants.TRIM_PERCENT, "0.1")[0]);

		SequenceSetProvider inputSequenceSetProvider = 
				createInputSequenceSetProvider(sequenceFileNames, minTaxa, maxTaxa, targetTaxa, targetMinGeneCount, repOnly);

		setInputSequenceSetProvider(inputSequenceSetProvider);

		//determine alignment trimming parameters
		setGblocksTrim(commandLineProperties.getValues(HandyConstants.GBLOCKS_TRIM, 
				HandyConstants.FALSE)[0].equalsIgnoreCase(HandyConstants.TRUE));
		setUniformTrim(commandLineProperties.getValues(HandyConstants.UNIFORM_TRIM, 
				HandyConstants.FALSE)[0].equalsIgnoreCase(HandyConstants.TRUE));

		if(inputsAligned) {
			loadInputAlignments(inputSequenceSetProvider.getAllSequenceSets());
		}

		//see if sequences should be pre-aligned
		boolean prealign = 
				commandLineProperties.getValues(HandyConstants.PREALIGN, 
						HandyConstants.FALSE)[0].equalsIgnoreCase(HandyConstants.TRUE);

		//find out how many threads to use for pre-aligning
		int alignmentThreadCount = Integer.parseInt(commandLineProperties.
				getValues(HandyConstants.ALIGN_THREADS, "1")[0]);
		//limit the number of threads to the number of available processors
		alignmentThreadCount = Math.min(alignmentThreadCount, 
				Runtime.getRuntime().availableProcessors());

		//determine number of total number of threads for tree building
		int treeThreads = Integer.parseInt(
				commandLineProperties.getValues(HandyConstants.TREE_THREADS, "1")[0]);

		//See if matrix evaluation should be done
		boolean doMatrixEvaluation = 
				!commandLineProperties.getValues(
						HandyConstants.MATRIX_EVALUATION, 
						HandyConstants.FALSE)[0].equals(HandyConstants.FALSE);

		boolean filterAlignmentsForCongruence = 
				!commandLineProperties.getValues(
						HandyConstants.CONGRUENCE_FILTER,
						HandyConstants.FALSE)[0].equals(HandyConstants.FALSE);


		if(prealign && !inputsAligned) {
			alignAll(inputSequenceSetProvider, alignmentThreadCount);			
		}

		if(filterAlignmentsForCongruence) {
			inputSequenceSetProvider = 
					filterForCongruence(inputSequenceSetProvider, percentileToRemove);
		}
		String mlMatrix = 
				commandLineProperties.getValues(
						HandyConstants.ML_MATRIX, "PROTGAMMAWAG")[0];

		if(doMatrixEvaluation) {
			//get list of matrices o be evaluated
			//if only one value for MATRIX_EVALUATION is present, and it it TRUE,
			//then evaluate a pre-selected list of matrices.
			String[] matrixList = 
					commandLineProperties.getValues(HandyConstants.MATRIX_EVALUATION);
			if(matrixList.length == 1 && matrixList[0].equals(HandyConstants.TRUE)) {
				//evaluate pre-selected list
				matrixList = new String[]{
						"PROTGAMMALG", 
						"PROTGAMMALGF",
						"PROTGAMMAWAG", 
						"PROTGAMMAWAGF",
						"PROTGAMMADAYHOFF", 
						"PROTGAMMADAYHOFFF", 
						"PROTGAMMADCMUT", 
						"PROTGAMMADCMUTF", 
						"PROTGAMMAJTT", 
						"PROTGAMMAJTTF", 
						"PROTGAMMAMTREV", 
						"PROTGAMMAMTREVF", 
						"PROTGAMMARTREV", 
						"PROTGAMMARTREVF", 
						"PROTGAMMACPREV", 
						"PROTGAMMACPREVF", 
						"PROTGAMMAVT", 
						"PROTGAMMAVTF", 
						"PROTGAMMABLOSUM62", 
						"PROTGAMMABLOSUM62F", 
						"PROTGAMMAMTMAM", 
						"PROTGAMMAMTMAMF", 
						"PROTGAMMAGTR"
				};
			}
			//if multiple values are given for MATRIX_EVALUATION, then they
			//should each be names of matrices to be evaluated
			//evaluateSubstitutionMatrices(inputSequenceSetProvider, matrixList, treeThreads);
			SequenceAlignment alignment = getConcatenatedAlignment(inputSequenceSetProvider);
			mlMatrix = evaluateSubstitutionMatricesForAlignment(alignment, 
					matrixList, treeThreads);
			if(verboseLevel > 0) {
				System.out.println("Preferred matrix is " + mlMatrix);
			}
		}

		boolean concatenated = 
				commandLineProperties.getValues(
						HandyConstants.CONCATENATED, 
						HandyConstants.FALSE)[0].
						equalsIgnoreCase(HandyConstants.TRUE);

		boolean geneWiseJackknife = 
				commandLineProperties.getValues(
						HandyConstants.GENE_WISE_JACKKNIFE, 
						HandyConstants.FALSE)[0].
						equalsIgnoreCase(HandyConstants.TRUE);

		boolean single = 
				commandLineProperties.getValues(
						HandyConstants.SINGLE, 
						HandyConstants.FALSE)[0].
						equalsIgnoreCase(HandyConstants.TRUE);

		//determine the number of threads to use for building the concatenated tree.
		//This comes out of the treeThreads amount, and is not in addition to it.
		//By default, the full treeThreads amount is used, but a smaller number can be given
		//which is especially useful for machines that don't support pthreads
		int fullTreeThreads = 
				Integer.parseInt(commandLineProperties.getValues(HandyConstants.FULL_TREE_THREADS, ""+treeThreads)[0]);

		int bootstrapReps = Integer.parseInt(
				commandLineProperties.getValues(HandyConstants.SUPPORT_REPS, "100")[0]);

		boolean nj = commandLineProperties.
				getValues(HandyConstants.NEIGHBOR_JOINING,
						HandyConstants.FALSE)[0].
						equals(HandyConstants.TRUE);

		String concatenatedTreeMethod = 
				commandLineProperties.getValues(
						HandyConstants.FULL_TREE_METHOD,
						HandyConstants.MAXIMUM_LIKELIHOOD)[0];

		String supportTreeMethod = 
				commandLineProperties.getValues(
						HandyConstants.SUPPORT_TREE_METHOD, 
						HandyConstants.FAST_TREE)[0];

		String singleTreeMethod = 
				commandLineProperties.getValues(
						HandyConstants.SINGLE_TREE_METHOD,
						HandyConstants.FALSE)[0];

		String outputFileName = commandLineProperties.
				getValues(HandyConstants.OUT, 
						"phyloPipeline.out")[0];

		if(single) {
			//if single was specified, get method, using parsimony as default
			singleTreeMethod = 
					commandLineProperties.getValues(
							HandyConstants.SINGLE_TREE_METHOD,
							HandyConstants.PARSIMONY)[0];
		} else {
			//even if -single was not specified
			// assume single trees should be built, 
			//if a single tree method is given
			single = !singleTreeMethod.equalsIgnoreCase(HandyConstants.FALSE);
		}


		boolean doToolCompare = commandLineProperties.getValues(HandyConstants.TOOL_COMPARE, HandyConstants.FALSE)[0].equalsIgnoreCase(HandyConstants.TRUE);
		//Parameters have been collected. 
		//Start the work now.
		if(doToolCompare) {
			//			compareTreeBuildingTools(inputSequenceSetProvider, mlMatrix, treeThreads);
			buildTreesFromSamples(inputSequenceSetProvider);
		}

		String finalTree = null;

		if(nj) {
			finalTree = buildConcatenatedNeighborJoiningTree(inputSequenceSetProvider);
		}

		if(concatenated && geneWiseJackknife) {
			if(verboseLevel > 0) {
				System.out.println("do concatenated and gene-wise jackknife");
			}
			finalTree = buildConcatenatedTreeWithGeneWiseJackKnifeSupport(
					inputSequenceSetProvider, concatenatedTreeMethod, 
					fullTreeThreads, bootstrapReps,
					supportTreeMethod, treeThreads, mlMatrix);
		} else if(concatenated && single) {
			if(verboseLevel > 0) {
				System.out.println("do concatenated and single");
			}
			buildConcatenatedAndSingleTrees(inputSequenceSetProvider, 
					concatenatedTreeMethod, 
					fullTreeThreads, 
					singleTreeMethod, treeThreads);
		} else if(concatenated) {
			if(verboseLevel > 0) {
				System.out.println("do concatenated");
			}
			ConcatenatedSequenceAlignment concatenatedAlignment = 
					getConcatenatedAlignment(inputSequenceSetProvider);
			finalTree = buildConcatenatedTree(concatenatedAlignment, fullTreeThreads, 
					bootstrapReps, concatenatedTreeMethod, mlMatrix);
		} else if(geneWiseJackknife){
			if(verboseLevel > 0) {
				System.out.println("do gene-wise jacknife");
			}
			buildGeneWiseJackknifeTrees(
					inputSequenceSetProvider,
					bootstrapReps, 
					supportTreeMethod, treeThreads,
					outputFileName, mlMatrix);
		} else if(single) {
			if(verboseLevel > 0) {
				System.out.println("do single");
			}
		}

		setFinalTree(finalTree);
	}

	private void setFinalTree(String tree) {
		finalTree = tree;
	}

	public String getFinalTree() {
		return finalTree;
	}

	private SequenceSetProvider filterForCongruence(
			SequenceSetProvider provider, double percentileToRemove) {
		SequenceSetProvider r = null;
		//		System.out.println(">PhylogenomicPipeline2.filterForCongruence()");
		System.out.println("filtering out most incongruent " + 
				" sequence sets, based on bipartition support");
		FastaSequenceSet[] sequenceSets = provider.getAllSequenceSets();
		ConcatenatedSequenceAlignment cat = getConcatenatedAlignment(provider);
		int topN = cat.getNTax() * 4;
		BipartitionSet fullBipartSet = 
				new BipartitionSet(cat.getBipartitionsForColumns(10000), topN);
//				new BipartitionSet(cat.getBipartitionsForColumns(10000));
		fullBipartSet.setTaxa(cat.getTaxa());
//		fullBipartSet = new BipartitionSet(fullBipartSet.getNonTrivialBipartitions(), topN);
//		fullBipartSet.setTaxa(cat.getTaxa());
		System.out.println("PhylogenomicPipeline2.filterForCongruence() topN (" + topN + ") bipartitions:");
		fullBipartSet.printBipartitionsAndCounts();
//		fullBipartSet.printNonTrivialBipartitionsAndCounts();

		SequenceAlignment[] alignments = cat.getAlignments();
		Bipartition[][] alignmentBiparts = 
				new Bipartition[alignments.length][];
		int[] alignmentMeanCosts = new int[alignments.length];

		long sumOfMeans = 0;
		long fullCostSum = 0;
		for(int i = 0; i < alignmentBiparts.length; i++) {
			alignmentBiparts[i] = cat.getBipartitionsForAlignment(i);
			int costSum = 0;
			for(int j = 0; j < alignmentBiparts[i].length; j++) {
				int bipartCost = fullBipartSet.getCost(alignmentBiparts[i][j]); 
				//System.out.println("cost for bipart " + i + ", " + j + ": " + bipartCost);
				costSum += bipartCost;
				fullCostSum += bipartCost;
				if(fullCostSum < 0) {
					//					System.out.println("negative fullCostSum: " + fullCostSum);
				}
			}

			alignmentMeanCosts[i] = costSum/alignments[i].getLength();
			sumOfMeans += alignmentMeanCosts[i];
			//			System.out.println("" + i + " length: " + alignments[i].getLength() 
			//					+ " fullCost: " + costSum + " meanCost: "+ alignmentMeanCosts[i]);
		}

		double[] sortMeans = new double[alignmentMeanCosts.length];
		for(int i = 0; i < sortMeans.length; i++) {
			sortMeans[i] = alignmentMeanCosts[i];
		}
		Arrays.sort(sortMeans);
		double meanOfMeans = StatisticsUtilities.getMean(sortMeans);
		double sdOfMeans = StatisticsUtilities.getStandardDeviation(sortMeans, meanOfMeans);
		//		System.out.println("distribution of bipart costs:");
		//		StatisticsUtilities.printDistribution(alignmentMeanCosts);

		System.out.println("bipart costs: " + meanOfMeans + " +/- " + sdOfMeans);
		//		int percentileIndex = (int) (sortMeans.length - (sortMeans.length*percentileToRemove));
		//		double maxCost = sortMeans[percentileIndex];
		double maxCost = meanOfMeans + (sdOfMeans*2);
//		maxCost = meanOfMeans + sdOfMeans;
//		maxCost = meanOfMeans;
		System.out.println("remove any alignments with meanCost > " +
				maxCost);

		//get the alignments being kept 
		ArrayList keepSequencesList = new ArrayList(alignmentMeanCosts.length);
		for(int i = 0;i < alignmentMeanCosts.length; i++) {
			if(alignmentMeanCosts[i] <= maxCost) {
				keepSequencesList.add(sequenceSets[i]);
			} else {
				System.out.println("discarding set " + sequenceSets[i].getName() 
						+ " cost: " + alignmentMeanCosts[i]);
			}
		}
		FastaSequenceSet[] keepSequences = new FastaSequenceSet[keepSequencesList.size()];
		keepSequencesList.toArray(keepSequences);
		System.out.println("Discarded: " + (sequenceSets.length - keepSequences.length) + "\tKept: " + keepSequences.length);

		r = new SequenceSetProviderImpl(keepSequences);

		return r;
	}

	private void loadInputAlignments(FastaSequenceSet[] sequenceSets) {
		sequenceSetToAlignment = new HashMap();
		if(msaTrimmer == null) {
			msaTrimmer = new MSATrimmer();
		}
		for(int i = 0; i < sequenceSets.length; i++) {
			try {
				SequenceAlignment alignment = 
						SequenceAlignmentParser.parseFastaAlignment(sequenceSets[i]);

				alignment.setName(sequenceSets[i].getName());
				alignment.setTitles(sequenceSets[i].getTitles());
				alignment.setTaxa(FastaUtilities.
						getTaxaFromTitles(sequenceSets[i].getTitles()));
				if(isGblocksTrim()) {
					alignment = msaTrimmer.trim(alignment);
				}
				sequenceSetToAlignment.put(sequenceSets[i], alignment);

			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}

	private void buildConcatenatedAndSingleTrees(
			SequenceSetProvider inputSequenceSetProvider2,
			String concatenatedTreeMethod, int fullTreeThreads,
			String singleTreeMethod, int treeThreads) {
		// TODO Auto-generated method stub

	}

	/**
	 * Creates a SequenceSetProvider with the specified constraints.
	 * 
	 * @param inputSequences
	 * @param minTaxa The fewest taxa a SequenceSet can contain
	 * @param maxTaxa The most taxa a SequenceSet can contain
	 * @param targetTaxa The ideal taxon count for SequenceSets. Sets with fewer
	 *                   taxa will be removed if doing so does not drop the
	 *                   number of sets below targetMinGeneCount.
	 * @param targetMinGeneCount The number of SequenceSets that should be kept,
	 *                           if possible. Any number with >= targetTaxa
	 *                           will be kept.
	 * @param repOnly Indicates if SequenceSets that are non-representative 
	 *                should be removed. Representative sets include a maximum 
	 *                of one sequence per taxon.
	 * @return
	 */
	private SequenceSetProviderImpl createInputSequenceSetProvider(
			String[] sequenceFileNames, int minTaxa, int maxTaxa, 
			int targetTaxa, int targetMinGeneCount, boolean repOnly) {
		if(verboseLevel > 0) {
			System.out.println(">createInputSequenceSetProvider() " + 
					sequenceFileNames.length + " " + minTaxa + " " + maxTaxa + " "
					+ targetTaxa + " " + targetMinGeneCount + " " + repOnly);
		}
		SequenceSetProviderImpl r =
				new SequenceSetProviderImpl(sequenceFileNames, minTaxa, maxTaxa);

		if(repOnly) {
			r.removeNonRepresentative();
		}

		if(r.getSequenceSetCount() > targetMinGeneCount) {
			//Remove SequenceSets with fewer than targetTaxa taxa, as long as there
			//are still >= targetMinGeneCount SequenceSets remaining
			for(int taxonCount = minTaxa+1; taxonCount <= targetTaxa; taxonCount++) {
				if(r.getSequenceSetCountForMinTaxa(taxonCount) >= targetMinGeneCount) {
					if(verboseLevel > 0) {
						System.out.println("PhylogenomicPipeline2.createInputSequenceSet() " +
								"removing sets smaller than " + taxonCount + 
								". There are " + r.getSequenceSetCount() + 
								" sequence sets before removing");
					}
					r.removeSetsWithFewerTaxaThan(taxonCount);
					if(verboseLevel > 0) {
						System.out.println("there are " + r.getSequenceSetCount() +
								" sequence sets after removing.");
					}
				}
			}
		}
		if(verboseLevel > 0) {
			System.out.println("<createInputSequenceSetProvider() " +
					r.getSequenceSetCount());
		}

		return r;
	}

	private static void defineHelp() {
		commands = new HashMap();

		commands.put(HandyConstants.MIN_TAXA, 
				"Miniumum number of taxa in a sequence set (smaller sets are filtered out)");
		commands.put(HandyConstants.MAX_TAXA, 
				"Maximum number of taxa in a sequence set (larger sets are filtered out)");
		commands.put(HandyConstants.REPRESETATIVE_ONLY, 
				"Only use representative sequence sets (maximum of one member per taxon)");
		commands.put(HandyConstants.ALIGN_THREADS, 
				"Number of threads (processors or processor cores) to use for multiple sequence alignment stage");
		commands.put(HandyConstants.PTHREADS, 
				"Number of threads to use for raxml tree-building stage");
		commands.put(HandyConstants.FILE,
				"\tName of input sequence set files (any number of file names may follow)");
		commands.put(HandyConstants.ALIGN_FILE, 
				"\tName of input alignment file (any number of file names may follow)/n" +
						"\t\tFor now, this should not be used in combination with " + HandyConstants.FILE);
		commands.put(HandyConstants.PREALIGN, 
				"Perform all alignments first, then begin tree-building stage");
		commands.put(HandyConstants.GBLOCKS_TRIM, 
				"Use Gblocks for automated trimming of alignments");
		commands.put(HandyConstants.UNIFORM_TRIM, 
				"Remove any uniform columns from alignments (these are not phylogenetically informative, but may affect the absolute branch lengths)");
		commands.put(HandyConstants.PARSIMONY, 
				"Build a parsimony tree without branch lengths, rather than a maximum-likelihood tree.");
		commands.put(HandyConstants.PARSIMONY_BL, 
				"Build a parsimony tree with branch lengths, rather than a maximum-likelihood tree.");
		commands.put(HandyConstants.CONCATENATED, 
				"Build a single tree from the concatenation of all alignments");
		commands.put(HandyConstants.SUPPORT_REPS, 
				"Number of bootstrap replicates to perform");
		commands.put(HandyConstants.HELP, 
				"\tPrint this help.");
		commands.put(HandyConstants.DIRECTORY, 
				"\tName of directory containing input faa files. All file from this directory will be loaded. Only applies if -" + HandyConstants.FILE + " is not specified");
		commands.put(HandyConstants.GENE_WISE_JACKKNIFE, 
				"For tree support values, build trees from a subset of the genes used for the full tree. ");
		commands.put(HandyConstants.TREE_THREADS, 
				"Number of processes to use for tree building step(s)");
		commands.put(HandyConstants.FULL_TREE_METHOD, 
				"Method used to build full concatenated tree. Default is Parsimony with Maximum Likelihood branch lengths (\"" + HandyConstants.PARSIMONY_BL + "\"). The other options are Maximum Likelihood (\"" + HandyConstants.MAXIMUM_LIKELIHOOD + "\") and FastTree (\"" + HandyConstants.FAST_TREE + "\")");
		commands.put(HandyConstants.SUPPORT_TREE_METHOD, 
				"Method used for building support trees (trees used for branch support values fo the full tree). Default is Parsimony (\"" + HandyConstants.PARSIMONY + "\"). The other current option is Maximum Likelihood (\"" + HandyConstants.MAXIMUM_LIKELIHOOD + "\")");
		commands.put(HandyConstants.TARGET_NTAX, 
				"\tThe ideal number of taxa per sequence set. By default, this is set to the value for " + HandyConstants.MAX_TAXA);
		commands.put(HandyConstants.TARGET_MIN_GENE_COUNT, 
				"\tThe ideal minimum number of genes (sequence sets) to be used. The value for " + HandyConstants.MIN_TAXA + " takes precedence over this.");
		commands.put(HandyConstants.FULL_TREE_THREADS, 
				"the number of threads to use for building the concatenated tree. This comes out of the treeThreads amount, and is not in addition to it. By default, the full treeThreads amount is used, but a smaller number can be given which is especially useful for machines that don't support pthreads");
		commands.put(HandyConstants.SINGLE, "\tBuild a tree for each input sequence set meeting the membership criteria.");
		commands.put(HandyConstants.SINGLE_TREE_METHOD, "Method used to build single-gene trees.");
		commands.put(HandyConstants.MATRIX_EVALUATION, "\tPerform matrix evaluation to determine which amino acid substitution matrix performs best for this dataset.");
		//		commands.put(HandyConstants.HELP, "\tPrint this help.");
		//		commands.put(HandyConstants.HELP, "\tPrint this help.");
	}

	/**
	 * Prints a help message, including a list of available commands.
	 */
	private static void printHelp() {
		Object[] commandNames = commands.keySet().toArray();

		Arrays.sort(commandNames);

		for(int i = 0; i < commandNames.length; i++) {
			Object commandDesc = commands.get(commandNames[i]);
			String line = "\t-" + commandNames[i] + "\t\t" + commandDesc;
			System.out.println(line);
		}

	}

	private void setInputSequenceSetProvider(SequenceSetProvider provider) {
		inputSequenceSetProvider = provider;
	}

	private SequenceSetProvider getInputSequenceSetProvider() {
		return inputSequenceSetProvider;
	}

	/**
	 * Returns the constraint tree, as a String, if a constraint tree
	 * has been provided. If no constraint tree has been provided
	 * then this method returns null.
	 */
	private String getConstraintTree() {
		return constraintTree;
	}

	/**
	 * Performs all alignments, using the specified number
	 * of Threads. 
	 */
	private void alignAll(SequenceSetProvider provider, int threads) {
		System.out.println("Aligning "
				+ provider.getSequenceSetCount() + " sequence sets using "
				+ threads + " threads");
		Iterator sequenceIterator = 
				provider.getSynchronizedIterator();

		//create and start all Threads
		Thread[] alignmentThreads = new Thread[threads];
		for(int i = 0; i < alignmentThreads.length; i++) {
			alignmentThreads[i] = 
					new Thread(new AlignmentRunnable(sequenceIterator));
			alignmentThreads[i].start();
		}

		//wait for all Threads to finish
		for(int i = 0; i < alignmentThreads.length; i++) {
			try {
				alignmentThreads[i].join();
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		logger.info("done aligning");
	}

	/**
	 * Returns a SequenceAlignment for the given FastaSequenceSet. If 
	 * an alignment file exists for this sequence set (determined by the
	 * set name and the derived alignment file name) then the alignment
	 * is taken from the file. If no alignment file exists, the sequences
	 * are aligned with the MultipleSequenceAligner.
	 * 
	 * @param sequenceSet
	 * @return
	 */
	private SequenceAlignment getAlignment(FastaSequenceSet sequenceSet) {
		SequenceAlignment r = null;

		boolean reuseAlignmentFiles = false;
		if(sequenceSetToAlignment == null) {
			sequenceSetToAlignment = new HashMap();
		}

		synchronized(sequenceSetToAlignment) {
			r = (SequenceAlignment) sequenceSetToAlignment.get(sequenceSet);
		}
		if(r == null) {
			//see if alignment file already exist. If it does, load the 
			//alignment from the file. If not, use the multiple sequence aligner
			//to perform the alignment.

			String alignmentFileName = new String(System.getProperty("user.dir") +
					File.separator + sequenceSet.getName() 
					+ HandyConstants.ALIGNMENT_FILE_SUFFIX);

			File alignmentFile = new File(alignmentFileName);
			if(reuseAlignmentFiles && alignmentFile.exists()) {
				try {
					r = SequenceAlignmentParser.
							parseFastaAlignmentFile(alignmentFileName);
					r.setName(sequenceSet.getName() + 
							HandyConstants.ALIGNMENT_FILE_SUFFIX);
				} catch (IOException e) {
					e.printStackTrace();
				} catch(ArithmeticException ae) {
					//this happens when the alignment file is empty. In this
					//case, the returned alignment is empty
				}
			} else {
				if(multipleSequenceAligner == null) {
					multipleSequenceAligner = new MultipleSequenceAligner();
				}
//				if(alignmentFile.exists()) {
//					logger.info("Realigning sequence set " + sequenceSet.getName());
//				}
				r = multipleSequenceAligner.getMSA(sequenceSet);
			}
			if(r != null) {
				if(isGblocksTrim()) {
					if(msaTrimmer == null) {
						msaTrimmer = new MSATrimmer();
					}
					//trim alignment
					r = msaTrimmer.trim(r);
				}

				if(isUniformTrim()) {
					if(msaTrimmer == null) {
						msaTrimmer = new MSATrimmer();
					}
					//trim out uninformative columns
					if(r != null) {
						r = msaTrimmer.trimUniformColumns(r);
					}
				}
				synchronized(sequenceSetToAlignment) {
					sequenceSetToAlignment.put(sequenceSet, r);
				}
			}
//			System.out.println("alignmentFileName: " + alignmentFileName);
//			new File(alignmentFileName).deleteOnExit();
		} 

		return r;
	}

	private boolean isUniformTrim() {
		return uniformTrim;
	}

	private void setUniformTrim(boolean trim) {
		uniformTrim = trim;
	}

	private boolean isGblocksTrim() {
		return gblocksTrim;
	}

	private void setGblocksTrim(boolean trim) {
		gblocksTrim = trim;
	}

	/**
	 * Builds a tree from the concatenation of all alignments.
	 * 
	 * @param provider
	 * @param threads number of thread to use for tree building. This is 
	 *                currently passed to RAxML as the pthreads value.
	 * @param bootstraps number of bootstrap replicates 
	 * @param concatenatedTreeMethod tree-building method to use
	 * @param mlMatrix likelihood matrix to use if maximum likelihood 
	 * 				   tree-building is used
	 * @return
	 */
	private String buildConcatenatedTree(ConcatenatedSequenceAlignment concatenatedAlignment, 
			int threads, int bootstraps, String concatenatedTreeMethod, 
			String mlMatrix) {
		String r = null;
		String concatenatedAlignmentName = runName + 
				System.currentTimeMillis();
		concatenatedAlignment.setName(concatenatedAlignmentName);

		logger.info("concatenated alignment: " +
				concatenatedAlignmentName);
		logger.info("\ttaxa: " + concatenatedAlignment.getNTax());
		logger.info("\tgenes: " + concatenatedAlignment.getAlignments().length);
		logger.info("\tcharacters: " + 
				concatenatedAlignment.getLength());

		PEPRTracker.setTaxonNumber(concatenatedAlignment.getNTax());
		PEPRTracker.setGeneNumber(concatenatedAlignment.getAlignments().length);
		PEPRTracker.setAlignedPositions(concatenatedAlignment.getLength());

		PhylogeneticTreeBuilder treeBuilder = createPhylogeneticTreeBuilder();
		if(concatenatedTreeMethod.equals(HandyConstants.MAXIMUM_LIKELIHOOD)) {
			treeBuilder.setTreeBuildingMethod(HandyConstants.MAXIMUM_LIKELIHOOD);
			treeBuilder.setMLMatrix(mlMatrix);
		}

		treeBuilder.useRaxmlBranchLengths(true);
		treeBuilder.setRunName(runName);
		treeBuilder.setTreeBuildingMethod(concatenatedTreeMethod);
		treeBuilder.setBootstrapReps(bootstraps);
		treeBuilder.setAlignment(concatenatedAlignment);
		String conTree = getConstraintTree();
		if(conTree != null) {
			treeBuilder.setConstraintTree(conTree);
		}
		treeBuilder.setProcesses(threads);

		//write alignment to file
		String alignmentFileName = runName + "_align.phy";
		String alignmentString = concatenatedAlignment.getAlignmentAsExtendedPhylipUsingTaxonNames();
		try {
			FileWriter fw = new FileWriter(alignmentFileName);
			fw.write(alignmentString);
			fw.flush();
			fw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		treeBuilder.run();

		r = treeBuilder.getTreeString();
		writeConsensusTreeFile(r);

		return r;
	}

	/**
	 * Creates a file containing the given tree for the given alignment.
	 * The file is named based on the run name.
	 * 
	 * @param alignment
	 * @param tree
	 */
	private void writeConsensusTreeFile(String tree) {
		String newickTreeFileName = runName + ".nwk";
		logger.info("tree: " + tree);
		System.out.println("writing tree to file: " + newickTreeFileName);
		try {
			FileWriter fw = new FileWriter(newickTreeFileName);
			fw.write(tree);
			fw.write("\n");
			fw.flush();
			fw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	/**
	 * Creates and returns an instance of PhylogeneticTreeBuilder.
	 * 
	 * @return PhylogeneticTreeBuilder
	 */
	private PhylogeneticTreeBuilder createPhylogeneticTreeBuilder() {
		PhylogeneticTreeBuilder r = new PhylogeneticTreeBuilder();

		//see if building parsimony tree only
		boolean parsimonyOnly = 
				commandLineProperties.getValues(HandyConstants.PARSIMONY, 
						HandyConstants.FALSE)[0].equals(HandyConstants.TRUE); 
		boolean parsimonyWithBL = 
				commandLineProperties.getValues(HandyConstants.PARSIMONY_BL, 
						HandyConstants.FALSE)[0].equals(HandyConstants.TRUE); 
		if(parsimonyOnly) {
			logger.info("build parsimony tree only");
			r.setTreeBuildingMethod(HandyConstants.PARSIMONY);
		} else if(parsimonyWithBL) {
			logger.info("build parsimony tree with branch lengths");
			r.setTreeBuildingMethod(HandyConstants.PARSIMONY_BL);
		} else {
			logger.info("build maximum likelihood tree");
			r.setTreeBuildingMethod(HandyConstants.MAXIMUM_LIKELIHOOD);
		}
		return r;
	}

	/**
	 * Selects a subset of the sequence sets and creates a concatenated
	 * sequence alignment from this subset. the number of sequences to be
	 * selected is specified by n. The n sequences are selected at random, and
	 * no sequence will appear more than once. n must not be larger than the
	 * total number of sequences available in the SequenceSetProvider.
	 * 
	 * This method is synchronized because there is some problem with 
	 * simultaneous access to the same alignment files and synchronizing
	 * within the getAlignment() method was creating threadlock problems.
	 * This becomes a temporary bottleneck in beginning individual 
	 * tree-building processes, but it isn't much of an issue compared with the
	 * time required for the tree building step.
	 * 
	 * @param sequences
	 * @return
	 */
	private synchronized ConcatenatedSequenceAlignment getConcatenatedAlignmentOfSubset(SequenceSetProvider provider, int n) {
		ConcatenatedSequenceAlignment r = null;

		if(verboseLevel > 1) {
			logger.info("getConcatenated alignment of " + n + " sequences");
		}
		int[] indices = RandomSetUtils.getRandomSet(n, 0,
				provider.getSequenceSetCount()-1, false);
		SequenceAlignment[] alignments = new SequenceAlignment[indices.length];
		for(int i = 0; i < alignments.length; i++) {
			FastaSequenceSet sequenceSet = provider.getSequenceSet(indices[i]); 
			alignments[i] = getAlignment(sequenceSet);
			if(alignments[i] != null) {
				alignments[i].setTaxa(FastaUtilities.getTaxaFromTitles(sequenceSet.getTitles()));
			}
		}
		r = MSAConcatenator.concatenate(alignments);
		return r;
	}

	/**
	 * Build a single tree from the concatenated alignments of all sequence
	 * sets, and build trees for jackknife replicates.
	 * 
	 * @param provider	Source for the Sequence Sets
	 * @param concatenatedTreeMethod Tree-building method for the full Tree
	 * @param fullTreeThreads Number of Threads to use to build the full tree.
	 * @param reps Number of jackknife trees to build
	 * @param supportTreeMethod Tree-building method for jacknife trees
	 * @param supportTreeThreads Number of Threads to use for support trees. The 
	 *                           fullTreeThreads comes out of this number. 
	 *                           Once the fullTree Threads finish, those threads
	 *                           become available to help build support trees.
	 * @return
	 */
	private String buildConcatenatedTreeWithGeneWiseJackKnifeSupport(
			SequenceSetProvider provider, String concatenatedTreeMethod, int fullTreeThreads, 
			int reps, String supportTreeMethod, int supportTreeThreads, 
			final String mlMatrix) {
		String r = null;
		String fullTreeString;
		final String finalConcatenatedTreeMethod = concatenatedTreeMethod;
		final SequenceSetProvider finalProvider = provider;
		final int finalThreadCount = fullTreeThreads;

		logger.info("building tree from concatenation of " + 
				provider.getSequenceSetCount() + " alignment using " 
				+ fullTreeThreads + " threads with method " + concatenatedTreeMethod);

		final ConcatenatedSequenceAlignment concatenatedAlignment = 
				getConcatenatedAlignment(provider);

		//estimate amount of RAM needed for each support tree
		//This is based on a formula derived from measurements of different
		//sized alignments used with FastTree
		//RAM in bytes =~ 100*L * N + 2000*L
		//where L is the alignment length and N is the number of taxa.
		//For this calculation, L will be half of the complete alignment,
		//since support trees are built from half of the sequences
		long L = concatenatedAlignment.getLength();
		long N = concatenatedAlignment.getNTax();
		long supportTreeRAM = 100l * (L/2) * N + 2000 * (L/2);
		long availableRAM = getAvailableRAMForTrees();

		int maxThreadsForRAM = (int) Math.floor(availableRAM/supportTreeRAM);

		logger.info("Estimated bytes of RAM needed for each tree-building thread: " +
				supportTreeRAM);
		logger.info("based on estimated RAM requirement for tree building with this alignment size, " +
				"the amount of available RAM will support " + maxThreadsForRAM 
				+ " parallel tree-building Threads.");

		if(maxThreadsForRAM == 0) {
			logger.info("There may not be enough RAM to support even a single tree-building thread, " +
					"but we will try it anyway");
			maxThreadsForRAM = 1;
		}

		supportTreeThreads = Math.min(supportTreeThreads, maxThreadsForRAM);

		Thread fullTreeThread = new Thread(){
			private String result;

			public void run() {
				result = buildConcatenatedTree(concatenatedAlignment, finalThreadCount, 0, 
						finalConcatenatedTreeMethod, mlMatrix);
			}

			public String toString() {
				return result;
			}
		};

		//start the full tree thread
		//		System.out.println("PhylogenomicPipeline2.buildConcatenatedTreeWithGeneWiseJackKnife() start full tree thread");
		fullTreeThread.start();

		//Have at least two threads wait for the full tree thread to finish.
		//This is to avoid using too much RAM, since the full tree will require
		//roughly twice as much RAM as each support tree
		int threadsToWait = 2;

		//threadsToWait might need to be increased for RAxML 
		//why is it * 80 * 8 instead of * 640? - because it's from raxml documentation:
		//http://sco.h-its.org/exelixis/web/software/raxml/:
		//"MEM(AA+GAMMA) = (n-2) * m * (80 * 8) bytes"  
		//where m is number of patterns, so this should overestimate
		long raxmlBytes = (N-2) * L * 80 * 8; 

		//Figure out if there is enough RAM to run the RAxML on the alignment 
		//(which is either used for the full tree, or for adding branch lengths
		//if FastTree is used for the full tree). If there is not enough RAM
		//available, hold back additional support tree threads until the 
		//full tree is finished, so the RAM allocated for those Threads will be
		//available for RAxML
		long spareBytes = availableRAM - (supportTreeRAM * supportTreeThreads);
		while(spareBytes < raxmlBytes && threadsToWait < supportTreeThreads) {
			spareBytes += supportTreeRAM;
			threadsToWait++;
		}

		logger.info("Estimated bytes of RAM needed for RAxML: " + 
				raxmlBytes);
		logger.info("support tree threads to wait for full tree: " + 
				threadsToWait);

		int jackknifeSize = provider.getSequenceSetCount() /2;
		int jackknifeReps = reps; 
		//start and wait for the jack knife threads
		String[] jackknifeTrees = getGeneWiseJackKnifeTrees(provider, 
				jackknifeSize, jackknifeReps, supportTreeThreads, 
				fullTreeThread, threadsToWait, supportTreeMethod, mlMatrix);

		//get the full tree string
		try {
			fullTreeThread.join();
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		fullTreeString = fullTreeThread.toString();

		//get full concatenated alignment
		//		SequenceAlignment concatenatedAlignment = 
		//			getConcatenatedAlignment(provider);
		//
		//		concatenatedAlignment.setName(runName+"_concatenated_with_support_" +
		//				System.currentTimeMillis()%10000000);

		r = TreeSupportDecorator.addSupportValues(fullTreeString, jackknifeTrees);
		//		System.out.println("tree: " + r);
		//write support trees to a file
		String supportTreeFileName = runName + ".sup";
		try {
			FileWriter fw = new FileWriter(supportTreeFileName);
			for(int i = 0; i < jackknifeTrees.length; i++) {
				fw.write(jackknifeTrees[i]);
				fw.write("\n");
			}
			fw.flush();
			fw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		writeConsensusTreeFile(r);

		return r;
	}

	private long getAvailableRAMForTrees() {
		return availableRAMForTrees;
	}

	/**
	 * Sets the amount of RAM available for tree-building (i.e. RAM not 
	 * being used by the Java process or other running processes).
	 * This is used to limit the number of parallel tree-building threads
	 * if the number of threads would otherwise cause memory problems.
	 * @param ram
	 */
	private void setAvailableRamForTrees(long ram) {
		if(availableRAMForTrees > 0) {
			availableRAMForTrees = ram;
		}
	}

	/**
	 * Parses the available memory string and sets the memory.
	 * @param memString
	 */
	private void setAvailableRAMForTrees(String memString) {
		long factor = 1;
		if(memString.endsWith("k") || memString.endsWith("K")) {
			//value is in kilobytes
			factor = 1024;
			memString = memString.substring(0, memString.length()-1);
		} else if(memString.endsWith("m") || memString.endsWith("M")) {
			//value is in megabytes
			factor = 1024 * 1024;
			memString = memString.substring(0, memString.length()-1);
		} else if(memString.endsWith("g") || memString.endsWith("G")) {
			//value is in gigabytes
			factor = 1024 * 1024 * 1024;
			memString = memString.substring(0, memString.length()-1);
		} else {
			//value is in bytes
		}

		try {
			long memoryValue = Long.parseLong(memString);
			long memoryBytes = memoryValue * factor;
			setAvailableRamForTrees(memoryBytes);
		} catch(NumberFormatException nfe) {
			System.out.println("Not able to parse the available ram parameter");
			nfe.printStackTrace();
		}

	}

	private void buildGeneWiseJackknifeTrees(
			SequenceSetProvider provider, int reps,
			String treeMethod, int treeThreads, String outputFileName, String mlMatrix) {
		logger.info("PhylogenomicPipeline2.buildGeneWiseJackknifeTrees() "
				+ reps + " trees using " + treeThreads + " threads and "
				+ treeMethod + " as tree-building method");

		int jackknifeSize = provider.getSequenceSetCount() /2;
		//start and wait for the jack knife threads
		String[] jackknifeTrees = getGeneWiseJackKnifeTrees(provider, 
				jackknifeSize, reps, treeThreads, 
				null, 0, treeMethod, mlMatrix);

		try {
			FileWriter fw  = new FileWriter(outputFileName);
			for(int i = 0; i < jackknifeTrees.length; i++) {
				fw.write(jackknifeTrees[i]);
				fw.write("\n");
			}
			fw.flush();
			fw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	/**
	 * Builds a single tree for each replicate, using a subset of the 
	 * sequence sets. For example, 100 replicates, each containing half of
	 * the sequence sets, can be used instead of bootstrap replicates
	 * to calculate support values for a single tree built from all of the
	 * sequences.
	 *
	 * @param size number of sequence sets per tree
	 * @param reps number of trees
	 * @param threadCount number of parallel threads to run
	 * @param waitFor A Thread to wait for before allowing one of the gene
	 *                subset tree building threads to proceed. This is to allow
	 *                maximum effective use of all available processes. The 
	 *                waitFor Thread is expected to be running the full 
	 *                concatenated tree. When it finishes, that processor 
	 *                becomes available for building subset trees.
	 * @param threadsToWait Number of support tree Threads that should wait
	 *                      for the waitFor Thread to finish before starting.
	 *                      This will probably be the number of Thread (or
	 *                      processes) being used to run the full tree. (For 
	 *                      example, the -T parameter with RAxML)
	 * @param mlMatrix 
	 */
	private String[] getGeneWiseJackKnifeTrees(SequenceSetProvider provider,
			int size, int reps, int threadCount, Thread waitFor, 
			int threadsToWait, String supportTreeMethod, String mlMatrix) {
		logger.info("PhylogenomicPipleline2.getGeneWiseJackKnifeTrees() " +
				"reps: " + reps + " threadsToWait: " + threadsToWait);
		String[] r = null;
		Thread[] threads = new Thread[threadCount];
		GeneSubsetTreeRunnable[] runnables = 
				new GeneSubsetTreeRunnable[threadCount];
		int[] repTracker = new int[]{reps};
		for(int i = 0; i < threads.length; i++) {
			if(i < threadsToWait) {
				runnables[i] = new GeneSubsetTreeRunnable(provider, size, 
						repTracker, waitFor);				
			} else {
				runnables[i] = new GeneSubsetTreeRunnable(provider, size, 
						repTracker, null);
			}

			runnables[i].setTreeBuildingMethod(supportTreeMethod);
			runnables[i].setMatrix(mlMatrix);
			String conTree = getConstraintTree();
			if(conTree != null) {
				runnables[i].setConstraintTree(conTree);
			}
			threads[i] = new Thread(runnables[i]);
			threads[i].start();
		}

		logger.info(threads.length + " support tree threads started (" + threadsToWait + " of them will wait for the main tree to finish)");
		//wait for threads to finish
		for(int i = 0; i < threads.length; i++) {
			try {
				threads[i].join();
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}

		//get all tree strings 
		ArrayList treeStrings = new ArrayList(reps);
		for(int i = 0; i < runnables.length; i++) {
			treeStrings.addAll(runnables[i].getTreeStrings());
		}

		r = new String[treeStrings.size()];
		treeStrings.toArray(r);
		return r;
	}



	private String buildConcatenatedNeighborJoiningTree(SequenceSetProvider provider) {
		logger.info("PhylogenomicPipeline2.buildConcatenatedNeighborJoiningTree()");
		String r = null;
		SequenceAlignment concatenatedAlignment =
				getConcatenatedAlignment(provider);
		String[] taxa = concatenatedAlignment.getTaxa();
		int[][] scoreMatrix = AlignmentUtilities.
				loadScoringMatrixFromResource(AlignmentUtilities.BLOSUM_62);
		double[][] similarityMatrix = 
				concatenatedAlignment.getPairwiseScores(scoreMatrix);

		r = TreeBuilder.getNJTopology(taxa, similarityMatrix, TreeBuilder.SIMILARITY);
		logger.info("neighbor-joining tree: " + r);
		return r;
	}

	private ConcatenatedSequenceAlignment getConcatenatedAlignment(SequenceSetProvider provider) {
		if(verbose) {
			logger.info("concatenating all " + 
					provider.getSequenceSetCount() + " sequence alignments");
		}
		ConcatenatedSequenceAlignment r = null;
		FastaSequenceSet[] seqSets = provider.getAllSequenceSets();

		ArrayList alignmentList = new ArrayList(seqSets.length);


		for(int i = 0; i < seqSets.length; i++) {
			SequenceAlignment alignment = getAlignment(seqSets[i]);
			if(alignment != null && alignment.getLength() > 0) {
				alignmentList.add(alignment);
			}
		}
		SequenceAlignment[] alignments = new SequenceAlignment[alignmentList.size()];
		alignmentList.toArray(alignments);
		alignmentList = null;

		r = MSAConcatenator.concatenate(alignments);

		//create file describing the homolog sets that were used in the 
		//concatenated alignment
		String hsFileName = runName + ".hs";
		String[] taxa = r.getTaxa();
		try {
			FileWriter fw = new FileWriter(hsFileName);
			//each row is for one homolog set (alignments, here)
			//each column is for one taxon.
			String tab = "\t";
			String newline = "\n";
			for(int i = 0; i < alignments.length; i++) {
				StringBuffer rowBuff = new StringBuffer();
				for(int j = 0; j < taxa.length; j++) {
					String member = null;
					int index = alignments[i].getTaxonIndex(taxa[j]);
					if(index >= 0) {
						member = alignments[i].getSequenceTitle(index);
						if(member != null) {
							member = member.split(" ")[0];
						} else {
							logger.error("Problem with alignment " + i + 
									" for taxon " + j + 
									" '" + taxa[j] + "'. index in alignment is "
									+ index + 
									" . here are the taxa: " );
							String[] alignmentTaxa = alignments[i].getTaxa();
							for(int k = 0; k < alignmentTaxa.length; k++) {
								logger.error("" + k + ": " 
										+ alignmentTaxa[k]);
							}
							logger.error("here are the titles: ");
							String[] titles = alignments[i].getSequenceTitles();
							for(int k = 0; k < titles.length; k++) {
								logger.error("" + k + ": " 
										+ titles[k]);
							}
						}
					}
					if(member == null) {
						member = "?";
					}
					rowBuff.append(member);
					if(j + 1 < taxa.length) {
						rowBuff.append(tab);
					}
				}
				fw.write(rowBuff.toString() + newline);

			}
			fw.flush();
			fw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		if(verbose) {
			logger.info("Done concatenating alignment. " +
					"Total concatenated alignment length is " +
					r.getLength());
		}
		return r;
	}

	/**
	 * Evaluate several amino acid substitution matrices for this alignment.
	 * not sure what to do with results yet. A simple thing would be to return
	 * the name of the best-performing matrix. The problem with this is that the
	 * best matrix may not be significantly better than the second-best, and 
	 * some people might grumble about that.
	 * 
	 * @param alignment
	 */
	private String evaluateSubstitutionMatricesForAlignment(
			SequenceAlignment alignment, String[] matrices, int threads) {
		String r = null;
		logger.info(">PhylogenomicPipeline2.evaluateSubstitutionMatricesForAlignment()");
		//for the alignment, build a parsimony tree with raxml
		RAxMLRunner rr = new RAxMLRunner(threads);
		rr.setAlignment(alignment);
		rr.setParsimonyOnly(true);
		rr.setBootstrapReps(0);
		rr.run();
		String parsimonyTree = rr.getParsimonyTree();

		//determine the matrices to be tested
		//		String[] matrices = new String[]{
		//				"PROTGAMMALG", 
		//				"PROTGAMMALGF",
		//				"PROTGAMMAWAG", 
		//				"PROTGAMMAWAGF",
		//				"PROTGAMMADAYHOFF", 
		//				"PROTGAMMADAYHOFFF", 
		//				"PROTGAMMADCMUT", 
		//				"PROTGAMMADCMUTF", 
		//				"PROTGAMMAJTT", 
		//				"PROTGAMMAJTTF", 
		//				"PROTGAMMAMTREV", 
		//				"PROTGAMMAMTREVF", 
		//				"PROTGAMMARTREV", 
		//				"PROTGAMMARTREVF", 
		//				"PROTGAMMACPREV", 
		//				"PROTGAMMACPREVF", 
		//				"PROTGAMMAVT", 
		//				"PROTGAMMAVTF", 
		//				"PROTGAMMABLOSUM62", 
		//				"PROTGAMMABLOSUM62F", 
		//				"PROTGAMMAMTMAM", 
		//				"PROTGAMMAMTMAMF", 
		//				"PROTGAMMAGTR"};
		//		
		//use raxml to score the parsimony tree with each matrix,
		//and keep the scores.
		double[] scores = new double[matrices.length];
		for(int i = 0 ; i < scores.length; i++) {
			scores[i] = getTreeScore(parsimonyTree, alignment, matrices[i], threads);
			logger.info("score for " + matrices[i] + ": " + scores[i]);
		}

		//see which matrix did the best
		int maxIndex = 0;
		double maxScore = scores[0];
		logger.info("score for " + matrices[0] + ": " + scores[0]);
		for(int i = 1; i < scores.length; i++) {
			if(scores[i] > maxScore) {
				maxIndex = i;
				maxScore = scores[i];
			}
			logger.info("score for " + matrices[i] + ": " + scores[i]);
		}
		r = matrices[maxIndex];
		System.out.println("best score: " + scores[maxIndex] + 
				" from " + matrices[maxIndex]);
		return r;
	}

	private void compareTreeBuildingTools(SequenceSetProvider provider, String matrix, int threads) {
		FastaSequenceSet[] sequenceSets = provider.getAllSequenceSets();

		TreeBuilderComparator tbc = new TreeBuilderComparator();

		for(int i = 60; i < sequenceSets.length; i+=10) {
			logger.info("concatenate " + i);
			SequenceAlignment alignment = getConcatenatedAlignmentOfSubset(provider, i);
			logger.info("concatenated alignment length: " + alignment.getLength());
			tbc.compareFastTreeAndRaxmlForAlignment(alignment, matrix, threads);
		}
	}

	private void buildTreesFromSamples(SequenceSetProvider provider) {
		FastaSequenceSet[] sequenceSets = provider.getAllSequenceSets();

		logger.info("seqs\tlength\tree");
		for(int i = 70; i < sequenceSets.length; i+=10) {
			SequenceAlignment alignment = getConcatenatedAlignmentOfSubset(provider, i);
			FastTreeRunner ftr = new FastTreeRunner();
			ftr.setAlignment(alignment);
			ftr.run();
			String treeString = ftr.getResult();
			String resultLine = "" + i + "\t" + alignment.getLength() + "\t" 
					+ treeString;
			logger.info(resultLine);
		}
	}

	private double getTreeScore(String tree, SequenceAlignment alignment, 
			String matrix, int threads) {
		double r = Double.NaN;
		RAxMLRunner rr = new RAxMLRunner(threads);
		rr.setAlignment(alignment);
		rr.setMatrix(matrix);
		rr.setPerSiteLogLikelihoods(true);
		rr.setPerSiteLLTrees(new String[]{tree});
		rr.run();
		String[] resultLines = rr.getPerSiteLLResultFile(); 

		String[] scores = resultLines[1].split("\\s+");

		r = 0.0;
		for(int i = 1; i < scores.length; i++) {
			r += Double.parseDouble(scores[i]);
		}
		return r;
	}

	/*
	 * Inner classes follow.
	 */	

	/**
	 * A runnable class to perform sequence alignments until
	 * all sequence sets have been aligned.
	 * 
	 * @author enordber
	 *
	 */
	private class AlignmentRunnable implements Runnable {

		private Iterator sequenceIterator;

		public AlignmentRunnable(Iterator sequenceIterator) {
			this.sequenceIterator = sequenceIterator;
		}

		public void run() {
			logger.info(">AlignmentRunnable.run() " + Thread.currentThread().getName());
			try {
				FastaSequenceSet sequenceSet =
						(FastaSequenceSet) sequenceIterator.next();
				while(sequenceSet != null) {
					SequenceAlignment alignment = getAlignment(sequenceSet);
					sequenceSet = (FastaSequenceSet) sequenceIterator.next();
				}
			} catch(Exception e) {
				System.out.println("EXCEPTION from eric's stupid 'catch any Exception' block in AlignmentRunnable");
				e.printStackTrace();
			}
			logger.info("<AlignmentRunnable.run() " + Thread.currentThread().getName());
		}
	}

	private class FaaFilter implements FilenameFilter {

		public boolean accept(File dir, String name) {
			boolean r = name.endsWith(".faa");
			return r;
		}
	}

	private class GeneSubsetTreeRunnable implements Runnable {
		private int size;
		private int[] reps;
		private SequenceSetProvider provider;
		private Thread waitFor;
		private ArrayList treeStrings = new ArrayList();
		private String treeBuildingMethod = HandyConstants.MAXIMUM_LIKELIHOOD;
		private String constraintTree;
		private String mlMatrix;

		/**
		 * 
		 * @param provider source for the sequence sets
		 * @param size number of genes from provider to include in each subset tree
		 * @param reps element 0 in this array tracks the remaining number of
		 *             trees to be built. This is shared by all active
		 *             GeneSubsetTreeRunnables to allow them to decrement a
		 *             single counter to synchronize their activity
		 * @param waitFor An optional Thread that this instance should 
		 *                wait for (by calling join()) before beginning work.
		 */
		public GeneSubsetTreeRunnable(SequenceSetProvider provider, int size,
				int[] reps, Thread waitFor) {
			this.provider = provider;
			this.size = size;
			this.reps = reps;
			this.waitFor = waitFor;
		}

		public void setConstraintTree(String conTree) {
			constraintTree = conTree;
		}

		public void setTreeBuildingMethod(String method) {
			treeBuildingMethod = method;
		}
		
		public void setMatrix(String mlMatrix) {
			this.mlMatrix = mlMatrix;
		}

		public void run() {
			logger.info(">GeneSubsetTreeRunnable.run() " + Thread.currentThread().getName());
			if(waitFor != null) {
				//wait for the waitFor Thread to complete before continuing
				try {
					logger.info("GeneSubsetTreeRunnable waiting for " + waitFor.getName() + " to join");
					waitFor.join();
					logger.info("GeneSubsetTreeRunnable waitFor Thread is finished. Begining " + Thread.currentThread().getName());
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
			int building = reps[0]--;
			while(building > 0) {
				logger.info("building tree: " + building + "  on " + Thread.currentThread().getName());
				ConcatenatedSequenceAlignment alignment = null;
				try {
					while(alignment == null || alignment.getLength() ==0) {
						alignment = getConcatenatedAlignmentOfSubset(provider, size);
						if(alignment == null || alignment.getLength() == 0) {
							System.out.println("&*^%#&^%#&*^%#&^%#&^%#&^$%#&%^$#^%$#bad concatenated alignment returned");
						}
					}
				} catch(Exception e) {
					logger.error("Exception trying to get subset " +
							"alignment. Will try again. Here's the " +
							"exception, but it's being handled: ");
					e.printStackTrace();
				}
				//				}
				PhylogeneticTreeBuilder treeBuilder = 
						new PhylogeneticTreeBuilder();
				treeBuilder.setTreeBuildingMethod(treeBuildingMethod);
				treeBuilder.setMLMatrix(mlMatrix);
				treeBuilder.setBootstrapReps(0);
				treeBuilder.setAlignment(alignment);
				if(constraintTree != null) {
					treeBuilder.setConstraintTree(constraintTree);
				}
				treeBuilder.run();
				String treeString = treeBuilder.getTreeString();
				System.out.println(Thread.currentThread().getName() + " tree for " + building + ": " + treeString);
				treeStrings.add(treeString);
				building = reps[0]--;
			}
			logger.info("<GeneSubsetTreeRunnable.run() " + Thread.currentThread().getName());
		}

		private boolean treeIsBifurcating(String treeString) {
			boolean r = false;
			BasicTree tree = new BasicTree(treeString);
			r = tree.isBifurcating();
			return r;
		}

		private ArrayList getTreeStrings() {
			return treeStrings;
		}

	}
}
