package edu.vt.vbi.ci.pepr.tree.pipeline;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Pattern;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.json.simple.JSONObject;

import edu.vt.vbi.ci.pepr.alignment.BlastRunner;
import edu.vt.vbi.ci.pepr.alignment.BlatRunner;
import edu.vt.vbi.ci.pepr.stats.StatisticsUtilities;
import edu.vt.vbi.ci.pepr.tree.AdvancedTree;
import edu.vt.vbi.ci.pepr.tree.BasicTree;
import edu.vt.vbi.ci.pepr.tree.PhylogeneticTreeRefiner;
import edu.vt.vbi.ci.util.CommandLineProperties;
import edu.vt.vbi.ci.util.CommandResults;
import edu.vt.vbi.ci.util.ExecUtilities;
import edu.vt.vbi.ci.util.HandyConstants;
import edu.vt.vbi.ci.util.IntPair;
import edu.vt.vbi.ci.util.PEPRTracker;
import edu.vt.vbi.ci.util.SequenceSetExtractor;
import edu.vt.vbi.ci.util.file.FastaSequenceFile;
import edu.vt.vbi.ci.util.file.FastaUtilities;
import edu.vt.vbi.ci.util.file.TextFile;

/**
 * This class will contain as much as possible of the full 
 * phylogenomic pipeline. This is being created to handle 
 * things that are currently handled by a perl script called
 * phyloPipeline.pl. That script wraps several different steps
 * in the pipeline, starting from genome faa files and resulting
 * in a species tree. If such a script is kept, I want it to be
 * for doing very basic environment set up, then calling this class.
 * All the other steps in the pipeline will be orchestrated by 
 * this class.
 * 
 * @author enordber
 *
 */
public class PhyloPipeline {

	private static final int MINIMUM_GENOME_COUNT = 4;

	private static Logger logger;
	private static HashMap<String,String> commands;
	private static long startTime;	
	static{
		logger = Logger.getLogger("PEPR");
		defineHelp();
	}

	private String runName;
	private int verbose = 1;
	private String tree;
	private String[] selectedOutgroupGenomes;
	private boolean writeJSON = false;

	public static void main(String[] args) {
		startTime = System.currentTimeMillis();

		try {
			PhyloPipeline pp = new PhyloPipeline(args);
			pp.run();
		} catch(Exception e) {
			e.printStackTrace();
			if(logger != null) {
				logger.info("Exiting with an Exception: " + e.toString());
			}
		}
		exit();
	}

	private static void exit() {
		if(logger != null) {
			logger.info("Exiting");
		}
		long endTime = System.currentTimeMillis();

		long elapsedMillis = endTime -startTime;
		long elapsedSeconds = elapsedMillis/1000;
		long elapsedMinutes = elapsedSeconds /60;
		long remainderSeconds = elapsedSeconds%60;
		if(logger != null) {
			logger.info("elapsed time: " + elapsedMinutes + "m " 
					+ remainderSeconds + "s");
		}
		System.exit(0);
	}

	public PhyloPipeline(String[] args) {
		CommandLineProperties clp = new CommandLineProperties(args);
		//set runName
		runName = "pepr-" + System.currentTimeMillis();
		runName = clp.getValues(HandyConstants.RUN_NAME, runName)[0];
		PEPRTracker.newTree(runName);

		String logfile = clp.getValues("logfile", runName+".log")[0];
		//set log file name, based on run name, if it has not already been set
		System.setProperty("logfile.name",System.getProperty("logfile.name", logfile));
		//		logger = Logger.getLogger("PEPR");
		logger.info("Starting " + new Date());
		System.out.println("logger configured with log file name '" + System.getProperty("logfile.name"));
		System.out.println("logger: " + logger);

		boolean printHelp = args == null || args.length == 0 ||
				clp.getValues(HandyConstants.HELP, 
						HandyConstants.FALSE)[0].equals(HandyConstants.TRUE) ||
						clp.getValues("h", 
								HandyConstants.FALSE)[0].equals(HandyConstants.TRUE) ||
								clp.getValues("?", 
										HandyConstants.FALSE)[0].equals(HandyConstants.TRUE);


		boolean useInstalledThirdPartyBinaries = true;
		if(clp.getValues(HandyConstants.PATRIC, HandyConstants.FALSE)[0].equals(HandyConstants.TRUE)) {
			FastaUtilities.setStripPipeAndSuffix(false);
			writeJSON = true;
			useInstalledThirdPartyBinaries = false;
		}

		useInstalledThirdPartyBinaries = clp.getValues(HandyConstants.USE_INSTALLED_THIRD_PARTY_BINARIES, ""+useInstalledThirdPartyBinaries)[0].equalsIgnoreCase(HandyConstants.TRUE);

		if(useInstalledThirdPartyBinaries) {
			System.out.println("using pre-installed binaries for third-party tools");
			setCommandPaths();
		} else {
			System.out.println("using system binaries for third-party tools");
			//check for specific version of muscle, and use if it is present
			String muscleVersion = "muscle-3.6";
			String fullCommandPath = ExecUtilities.getCommandPath(muscleVersion);
			if(fullCommandPath != null && fullCommandPath.length() > 0) {
				logger.info("set path for '" + muscleVersion + "' to '" + fullCommandPath + "'" );
				System.setProperty("muscle", fullCommandPath);
			}

		}

		//check for required programs
		checkForRequiredPrograms();

		if(printHelp) {
			printHelp();
			exit();
		}

		logger.info("run name: " + runName);

		//extract the necessary parameters from the CommandLineProperties

		//determine if a properties file should be loaded
		String confFileName = clp.getValues(HandyConstants.CONF_FILE_COMMAND,
				HandyConstants.FALSE)[0];
		if(!confFileName.equalsIgnoreCase(HandyConstants.FALSE)) {
			//load the properties file and add any new properties. New
			//properties override any loaded from the file.
			try {
				logger.info("loading configuration file: " + confFileName);
				CommandLineProperties confCLP = 
						CommandLineProperties.loadFromFile(confFileName);
				clp.addArgs(confCLP);
			} catch (IOException e) {
				logger.log(Level.DEBUG, "problem loading configuration file: " + 
						confFileName);
				logger.log(Level.DEBUG, "Problem trying to load configuration file: " +
						confFileName);
				logger.log(Level.DEBUG,e.getMessage());
			}
			clp.remove("-"+HandyConstants.CONF_FILE_COMMAND);
		}

		CommandLineProperties initialCLP = clp;

		//See if any 'track' is specified, and load the relevant properties.
		//Tracks can be a starting point and be modified by other properties.
		String track = clp.getValues(HandyConstants.TRACK, HandyConstants.TRACK)[0];
		clp = new CommandLineProperties();
		if(!track.equals(HandyConstants.FALSE)) {
			String[] trackProperties = getTrackProperties(track);
			clp.addArgs(trackProperties);
		}

		clp.addArgs(initialCLP.getArgs());

		//write the CommandLineProperties to a log file that can be re-used
		try {
			writePropertiesToLogFile(clp);
		} catch (IOException e) {
			e.printStackTrace();
			logger.log(Level.DEBUG,e.getMessage());
		}

		//Determine number of Threads to use. By default, every thread parameter
		//will be set to the number of processors available.
		int maxThreads = Runtime.getRuntime().availableProcessors();
		String maxThreadsS = "" + maxThreads;

		//An upper limit, other than the number of available processors, may be
		//placed on all thread values with the MAX_CONCURRENT_PROCESSES option
		maxThreads = Integer.parseInt(
				clp.getValues(HandyConstants.MAX_CONCURRENT_PROCESS_PARAM, 
						maxThreadsS)[0]);
		maxThreadsS = "" + maxThreads;

		int alignThreads;
		int treeThreads;
		int fullTreeThreads; //default to treeThreads, but can be overridden
		int homologyThreads;
		int hmmThreads;
		int mclThreads;

		alignThreads = Integer.parseInt(
				clp.getValues(HandyConstants.ALIGN_THREADS, maxThreadsS)[0]);
		treeThreads = Integer.parseInt(
				clp.getValues(HandyConstants.TREE_THREADS, maxThreadsS)[0]);
		fullTreeThreads = Integer.parseInt(
				clp.getValues(HandyConstants.FULL_TREE_THREADS, ""+treeThreads)[0]);
		homologyThreads = Integer.parseInt(
				clp.getValues(HandyConstants.HOMOLOGY_THREADS, maxThreadsS)[0]);
		mclThreads = Integer.parseInt(
				clp.getValues(HandyConstants.MCL_THREADS, maxThreadsS)[0]);

		hmmThreads = Integer.parseInt(
				clp.getValues(HandyConstants.HMM_THREADS, maxThreadsS)[0]);


		//Examine input sequence files to see how many taxa there are.
		//If the files are not already exactly one file per taxon, then
		//create temporary sequence files with one file per taxon.
		String[] inputSequenceFileNames = 
				clp.getValues(HandyConstants.GENOME_FILE);
		if(inputSequenceFileNames == null) {
			logger.info("No input sequence files were provided. " +
					"Please provide sequence files with the command -"
					+ HandyConstants.GENOME_FILE);
			exit();
		}

		FastaSequenceFile[] inputSequenceFiles = null;

		try {
			inputSequenceFiles = loadSequenceFiles(inputSequenceFileNames);
		} catch (IOException e) {
			e.printStackTrace();
			logger.log(Level.DEBUG,e.getMessage());
		}

		PEPRTracker.setInputSequenceFiles(inputSequenceFiles);

		//find out if blast/blat searches should use only a single member
		//from each species. If true, then inputSequenceFiles will be filtered
		//to remove duplicate taxa from the same species. This allows HMMs
		//to be built from clusters without duplicate species. The family
		//members from the duplicate species are added in based on hmmsearch
		//results
		boolean uniqueSpecies = clp.getValues(HandyConstants.UNIQUE_SPECIES,
				HandyConstants.FALSE)[0].equalsIgnoreCase(HandyConstants.TRUE);
		if(uniqueSpecies) {
			logger.info("Use unique species for first round homology search: " 
					+ uniqueSpecies);
		}

		boolean uniqueGenus = clp.getValues(HandyConstants.UNIQUE_GENUS,
				HandyConstants.FALSE)[0].equalsIgnoreCase(HandyConstants.TRUE);

		FastaSequenceFile[] homologySearchSequenceFiles = inputSequenceFiles;
		if(uniqueSpecies) {
			homologySearchSequenceFiles = 
					filterOutDuplicateSpecies(homologySearchSequenceFiles);
			if(homologySearchSequenceFiles.length < MINIMUM_GENOME_COUNT) {
				System.out.println("There are not enough unique species to use the unique species filter, so all genomes will be used.");
				homologySearchSequenceFiles = inputSequenceFiles;
			}
			Arrays.sort(homologySearchSequenceFiles);
		} else if(uniqueGenus) {
			homologySearchSequenceFiles = 
					filterOutDuplicateGenera(homologySearchSequenceFiles);
			Arrays.sort(homologySearchSequenceFiles);			
		}
		boolean treeFromUnique = clp.getValues(HandyConstants.UNIQUE_ONLY,
				HandyConstants.FALSE)[0].equalsIgnoreCase(HandyConstants.TRUE);
		if(treeFromUnique) {
			logger.info("Using only unique species for tree.");
			inputSequenceFiles = homologySearchSequenceFiles;
		}

		//count taxa in input files
		String[] distinctTaxa = getDistinctTaxa(inputSequenceFiles);
		int taxonCount = distinctTaxa.length;
		logger.info(taxonCount + " distinct taxa found in " + 
				inputSequenceFiles.length + " sequence files.");

		//If called for, run Blast or Blat searches
		boolean runBlast = clp.getValues(HandyConstants.HOMOLOGY_SEARCH_METHOD, 
				HandyConstants.FALSE)[0].equals(HandyConstants.BLAST);
		boolean runBlat = clp.getValues(HandyConstants.HOMOLOGY_SEARCH_METHOD, 
				HandyConstants.FALSE)[0].equals(HandyConstants.BLAT);

		//Blast and Blat parameters that are being hard-coded for now, but may be
		//made accessible by command-line options eventually
		int hitsPerQuery = 1;
		float evalueCutoff = 0.1f;
		int minIdentity = 10; //for blat def 10
		int minScore = 15; //for blat def 15

		logger.info("runBlast: " + runBlast);

		TextFile hitPairFile = null;
		if(runBlast) {
			logger.info("run blast on " + homologySearchSequenceFiles.length + " sequence files");
			hitPairFile = runBlast(homologySearchSequenceFiles, 
					homologyThreads, hitsPerQuery, evalueCutoff);
		} else if(runBlat) {
			hitPairFile = runBlat(homologySearchSequenceFiles, 
					homologyThreads, hitsPerQuery, evalueCutoff, 
					minIdentity, minScore);
		} else {
			//check for a file of homology results
			String homologyFile = clp.getValues(HandyConstants.HOMOLOGY_SEARCH_METHOD, 
					HandyConstants.FALSE)[0];

			//If Blast/Blat result files are given, load results.			
			if(!homologyFile.equals(HandyConstants.FALSE)) {
				try {
					if(verbose > 0) {
						logger.info("loading homology file: " + homologyFile);
					}
					hitPairFile = new TextFile(homologyFile);
				} catch (IOException e) {
					e.printStackTrace();
					logger.log(Level.DEBUG,e.getMessage());
				}
			}
		}
		//If called for, filter Blast/Blat hits for bidirectional pairs
		String bidirectional = 
				clp.getValues(HandyConstants.BIDIRECTIONAL, 
						HandyConstants.FALSE)[0];
		boolean bidirectionalFilter = bidirectional.equals(HandyConstants.TRUE);
		if(bidirectionalFilter) {
			try {
				TextFile unfilteredHitPairFile = hitPairFile;
				hitPairFile = filterForBidirectional(hitPairFile, 11);
				//remove the unfiltered version of the file, which may 
				//be using a lot of disk space
				//unfilteredHitPairFile.getFile().delete();
			} catch (IOException e) {
				e.printStackTrace();
				logger.log(Level.DEBUG,e.getMessage());
			}
		} else {
			//Filter hit pair file to contain only ID1	ID2	score
			try {
				hitPairFile = filterHitPairFile(hitPairFile);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

		//Create input file for MCL: ID1	ID2	score

		//Run MCL
		String inflation = clp.getValues(HandyConstants.INFLATION, "1.5")[0];
		TextFile mclResult = runMCL(hitPairFile, inflation, mclThreads);

		if(verbose > 0) {
			try {
				checkMCLResults(mclResult);
			} catch (IOException e) {
				e.printStackTrace();
				logger.log(Level.DEBUG,e.getMessage());
			}
		}

		//Create the sequence set faa files based on the results from MCL
		//create the directory for the Homolog Group files
		String hgDirName = "hg_" + runName;
		File hgDir = new File(hgDirName);
		hgDir.mkdir();

		SequenceSetExtractor sse = 
				new SequenceSetExtractor(mclResult.getFile().getAbsolutePath(),
						inputSequenceFiles, hgDirName + "/set", "faa");

		float minTaxaMultiplier = 
				Float.parseFloat(clp.getValues(HandyConstants.MIN_TAXA_MULTIPLIER,
						"0.99")[0]);
		int maxTaxa = 
				Integer.parseInt(clp.getValues(HandyConstants.MAX_TAXA,
						"" + inputSequenceFiles.length)[0]);

		//if only unique species were used for homology search, adjust 
		//the maxTaxa value accordingly, by subtracting the number of genomes 
		//that were filtered out by the unique species filter
		int fullUniqueDifference = 
				inputSequenceFiles.length - homologySearchSequenceFiles.length;
		int hsMaxTaxa = maxTaxa;
		if(fullUniqueDifference > 0) {
			hsMaxTaxa -= fullUniqueDifference;
			logger.info("adjusting maxTaxa for homology set enhancement from "
					+ maxTaxa + " to "
					+ hsMaxTaxa + 
					" because of unique species filter");
		}

		int minTaxa = (int) (maxTaxa * minTaxaMultiplier);
		minTaxa = Integer.parseInt(clp.getValues(HandyConstants.MIN_TAXA,
				"" + minTaxa)[0]);


		//if a set of outgroup genomes have been provided, add them to the 
		//list of input sequence files. These sequences are not included in 
		//the initial homolog group creation, but they are included in the
		//enhanced groups. This is to speed up the homology search step and
		//to avoid the HGs being biased by outgoup taxa.

		String[] outgroupSequenceFileNames = clp.getValues(HandyConstants.OUTGROUP);
		FastaSequenceFile[] outgroupSequenceFiles = new FastaSequenceFile[0];
		if(outgroupSequenceFileNames != null) {
			logger.info("adding " + outgroupSequenceFileNames.length + 
					" outgroup files");
			Arrays.sort(outgroupSequenceFileNames);
			try {
				outgroupSequenceFiles = 
						loadSequenceFiles(outgroupSequenceFileNames);
			} catch (IOException e) {
				e.printStackTrace();
				logger.log(Level.DEBUG,e.getMessage());
			}
		}

		PEPRTracker.setOutgroupPoolSequenceFiles(outgroupSequenceFiles);
		int outgroupCount = Integer.parseInt(clp.getValues(HandyConstants.OUTGROUP_COUNT, ""+outgroupSequenceFiles.length)[0]);

		Arrays.sort(inputSequenceFiles);

		//In the future, this is the place to add additional Homolog groups
		//stuff, like using HMMs to look for group members that got missed.
		//The future is now.
		int enhancerMinTaxa = Math.max(3, hsMaxTaxa/2);

		//pass a String[][] to enhanceHomologGroups so the enhancer can fill 
		//it with the names of the retained outgroup sequences
		String[][] retainedOutgroupHolder = new String[1][];
		hgDirName = enhanceHomologGroups(hgDirName, inputSequenceFiles,
				alignThreads,
				hmmThreads, 
				enhancerMinTaxa, 
				(int)(hsMaxTaxa*1.5), outgroupSequenceFiles, outgroupCount,
				retainedOutgroupHolder);

		//now add the number of outgroup taxa to maxTaxa, but keep targetTaxa at
		//the current minTaxa value
		int targetNTax = maxTaxa;
		if(outgroupSequenceFileNames != null) {
			maxTaxa += outgroupCount;
		}

		for(int i = 0; i < retainedOutgroupHolder[0].length; i++) {
			logger.info("Selected Outgroup Genome: " + retainedOutgroupHolder[0][i]);
		}
		setSelectedOutgroupGenomes(retainedOutgroupHolder[0]);

		//Run the PhylogenomicPipeline2.jar. The naming here is going to get
		//confusing, and probably needs to be addressed. This is the part that
		//filters the Homolog Groups for representative groups and for 
		//certain numbers of taxa, then does alignment, alignment trimming,
		//and tree-building.

		long hgTime = System.currentTimeMillis();
		long elapsedMillis = hgTime -startTime;
		long elapsedSeconds = elapsedMillis/1000;
		long elapsedMinutes = elapsedSeconds /60;
		long remainderSeconds = elapsedSeconds%60;
		logger.info("Elapsed time for homolog group creation: " 
				+ elapsedMinutes + "m " + elapsedSeconds + "s");

		//add properties that aren't coming in directly from the 
		String[] phylogenomicPipelineOptions = new String[]{
				"-" + HandyConstants.DIRECTORY,
				hgDirName,
				"-" + HandyConstants.MIN_TAXA,
				"" + minTaxa,
				"-" + HandyConstants.MAX_TAXA,
				"" + maxTaxa,
				"-" + HandyConstants.TARGET_NTAX,
				"" + targetNTax,
				"-" + HandyConstants.TREE_THREADS,
				"" + treeThreads,
				"-" + HandyConstants.FULL_TREE_THREADS,
				"" + fullTreeThreads,
				"-" + HandyConstants.ALIGN_THREADS,
				"" + alignThreads,
				"-" + HandyConstants.MCL_THREADS,
				"" + mclThreads
		};
		clp.addArgs(phylogenomicPipelineOptions);
		String[] phyloArgs = clp.getArgs();

		PhylogenomicPipeline2 pgp = new PhylogenomicPipeline2(phyloArgs);

		String firstRoundTree = pgp.getFinalTree();
		setTree(firstRoundTree);

		//do tree refinement if appropriate
		boolean refine = 
				!clp.getValues(HandyConstants.REFINE, HandyConstants.FALSE)[0]
						.equalsIgnoreCase(HandyConstants.FALSE);
		logger.info("refine: " + refine);
		if(refine) {
			refineTree(firstRoundTree, clp, inputSequenceFiles, outgroupSequenceFiles);
		}

		PEPRTracker.setTree(getTree());

		try {
			printTreeAndOutgroupFile();
		} catch (IOException e2) {
			System.out.println("There was a problem writing the tree and outgroup (*.tog) output file: ");
			e2.printStackTrace();
		}

		try {
			printRootedFinalTree();
		} catch (IOException e1) {
			System.out.println("There was a problem writing the final rooted tree output files: ");
			e1.printStackTrace();
		}

		if(writeJSON) {
			try {
				printJSONTreeAndMetadata();
			} catch (IOException e) {
				System.out.println("There was a problem writing the PATRIC json output file: ");
				e.printStackTrace();
			}
		}
	}

	private void printJSONTreeAndMetadata() throws IOException {
		String newickTree = getTree();
		AdvancedTree tree = new AdvancedTree(newickTree);
		tree.setOutGroup(getSelectedOutgroupGenomes());
		String rootedNewickTree = tree.getTreeString(true, true);
		String jsonTree = newickToPATRICJSON(rootedNewickTree, getSelectedOutgroupGenomes(), "unknown", "unknown");
		String jsonFileName = runName + ".json";
		logger.info("writing json output to file " + jsonFileName);
		FileWriter fw = new FileWriter(jsonFileName);
		fw.write(jsonTree);
		fw.close();


	}

	private void printRootedFinalTree() throws IOException {
		String newickTree = getTree();
		AdvancedTree tree = new AdvancedTree(newickTree);
		tree.setOutGroup(getSelectedOutgroupGenomes());
		String rootedNewickTree = tree.getTreeString(true, true);
		String rootedNewickFileName = runName + "_final_rooted.nwk";

		FileWriter fw = new FileWriter(rootedNewickFileName);
		fw.write(rootedNewickTree);
		fw.close();

		String rootedJSONFileName = runName + "_final_rooted.json";
		String rootedJSONTree = tree.getTreeJSON();
		fw = new FileWriter(rootedJSONFileName);
		fw.write(rootedJSONTree);
		fw.close();
	}

	private static String newickToPATRICJSON(String newickTree, String[] outgroupGenomes, String taxonName, String taxonRank) {
		String r = null;
		String delimiter = "_@_";
		int nameIndex = 0;
		int idIndex = 1;
		String idOnlyNewick = newickTree;
		BasicTree tree = new BasicTree(newickTree);
		String[] leaves = tree.getLeaves();
		ArrayList<String> idList = new ArrayList<String>(leaves.length);
		HashMap<String,String> nameToId = new HashMap<String,String>();
		HashSet<String> outgroup = new HashSet<String>();
		for(String og: outgroupGenomes) {
			String[] parts = og.split(delimiter);
			outgroup.add(parts[nameIndex].replaceAll("_", " "));
		}

		JSONObject labelJSON = new JSONObject();
		//verify that leaves are in <genome_name>@<genome_id> format
		//extract genome name and genome id to make map
		//replace leaves with just genome id
		for(String leaf: leaves) {
			if(leaf.contains(delimiter)) {
				String[] parts = leaf.split(delimiter);
				String name = parts[nameIndex].replaceAll("_", " "); 
				String id = parts[idIndex];
				idOnlyNewick = idOnlyNewick.replaceFirst(leaf, parts[idIndex]);
				nameToId.put(name, id);
				labelJSON.put(id,name);
				idList.add(id);
			}
		}

		JSONObject outgroupJSON = new JSONObject();
		String[] names = outgroup.toArray(new String[0]);
		Arrays.sort(names);
		for(String name: names) {
			String id = nameToId.get(name);
			outgroupJSON.put(id,name);
		}

		JSONObject infoJSON = new JSONObject();
		infoJSON.put("outgroups", outgroupJSON);

		String genomeCount = "" + (leaves.length - outgroupGenomes.length);
		infoJSON.put("count", genomeCount);

		infoJSON.put("taxon_name", taxonName);
		infoJSON.put("taxon_rank", taxonRank);
		infoJSON.put("ids", idList);

		JSONObject fullJSON = new JSONObject();
		fullJSON.put("info", infoJSON);
		fullJSON.put("labels", labelJSON);
		fullJSON.put("tree", idOnlyNewick);
		//		System.out.println(idOnlyNewick);
		r = fullJSON.toJSONString();
		//		System.out.println(r);
		return r;
	}

	private void printTreeAndOutgroupFile() throws IOException {
		String fileName = runName + ".tog";
		String[] outgroups = getSelectedOutgroupGenomes();
		FileWriter fw = new FileWriter(fileName);
		fw.write("trees[\'" + runName + "\'] = \""+ getTree() + "\";\n");
		if(outgroups != null && outgroups.length > 0) {
			fw.write("outgroups[\'" + runName + "\'] = [\"" + outgroups[0] + "\"");
			for(int i = 1; i < outgroups.length; i++) {
				fw.write(", \"" + outgroups[i] + "\"");
			}
			fw.write("];\n");
		}
		fw.flush();
		fw.close();
	}

	private void setTree(String tree) {
		this.tree = tree;
	}

	public String getTree() {
		return tree;
	}

	private void refineTree(String initialTree, CommandLineProperties clp,
			FastaSequenceFile[] ingroupSequences, FastaSequenceFile[] outgroupSequences) {
		PhylogeneticTreeRefiner refiner = 
				new PhylogeneticTreeRefiner(initialTree, clp, 
						ingroupSequences, outgroupSequences);
		setTree(refiner.getMostRefinedTreeString());
	}

	/**
	 * Returns a filtered list of sequence files containing only one member per
	 * species. Duplicate species/sub-species/strains are removed.
	 * Species is determined by taking the first two tokens of the taxon name.
	 * Current implementation assumes one genome per file.
	 * @param homologySearchSequenceFiles
	 * @return
	 */
	public static FastaSequenceFile[] filterOutDuplicateSpecies(
			FastaSequenceFile[] homologySearchSequenceFiles) {
		FastaSequenceFile[] r = null;
		ArrayList keepFiles = new ArrayList(homologySearchSequenceFiles.length/2);
		//speciesKept stores each species name being kept
		HashMap speciesKept = new HashMap();
		Pattern spacePat = Pattern.compile("_");
		String space = " ";
		for(int i = 0; i < homologySearchSequenceFiles.length; i++) {
			String[] taxa =  homologySearchSequenceFiles[i].getTaxa();
			if(taxa != null && taxa.length > 0) {
				String taxon = taxa[0];
				String[] taxonTokens = spacePat.split(taxon);
				if(taxonTokens.length > 1) {
					String species = taxonTokens[0] + space + taxonTokens[1];
					if(!speciesKept.containsKey(species)) {
						keepFiles.add(homologySearchSequenceFiles[i]);
						speciesKept.put(species,homologySearchSequenceFiles[i]);
					}
					else {
						FastaSequenceFile previousSpeciesFile = 
								(FastaSequenceFile)speciesKept.get(species);
						int previousGeneCount = previousSpeciesFile.getSequenceCount();
						int currentGeneCount = homologySearchSequenceFiles[i].getSequenceCount();
						if(currentGeneCount > previousGeneCount) {
							speciesKept.put(species, homologySearchSequenceFiles[i]);
						}
					}
				}
			}
		}

		r = new FastaSequenceFile[speciesKept.size()];
		speciesKept.values().toArray(r);
		return r;
	}

	/**
	 * Returns a filtered list of sequence files containing only one member per
	 * species. Duplicate species/sub-species/strains are removed.
	 * Species is determined by taking the first two tokens of the taxon name.
	 * Current implementation assumes one genome per file.
	 * @param homologySearchSequenceFiles
	 * @return
	 */
	private FastaSequenceFile[] filterOutDuplicateGenera(
			FastaSequenceFile[] homologySearchSequenceFiles) {
		FastaSequenceFile[] r = null;
		ArrayList<FastaSequenceFile> keepFiles = new ArrayList(homologySearchSequenceFiles.length/2);
		//speciesKept stores each species name being kept
		HashMap<String, FastaSequenceFile> generaKept = new HashMap();
		Pattern spacePat = Pattern.compile("_");
		for(int i = 0; i < homologySearchSequenceFiles.length; i++) {
			String taxon = homologySearchSequenceFiles[i].getTaxa()[0];
			logger.info("checking taxon: " + taxon);
			String[] taxonTokens = spacePat.split(taxon);
			if(taxonTokens.length > 1) {
				String genus = taxonTokens[0];
				if(!generaKept.containsKey(genus)) {
					logger.info("new genus: " + genus + "\t" 
							+ taxon);
					keepFiles.add(homologySearchSequenceFiles[i]);
					generaKept.put(genus,homologySearchSequenceFiles[i]);
				}
				else {
					logger.info("duplicate genus: " + genus + 
							"\t" + taxon);
					FastaSequenceFile previouSpeciesFile = 
							(FastaSequenceFile)generaKept.get(genus);
					int previousGeneCount = previouSpeciesFile.getSequenceCount();
					int currentGeneCount = homologySearchSequenceFiles[i].getSequenceCount();
					if(currentGeneCount > previousGeneCount) {
						logger.info("replacing genus representative for "
								+ genus + " with " + taxon + ". " +
								currentGeneCount + " > " + previousGeneCount 
								+ " genes");
						generaKept.put(genus, homologySearchSequenceFiles[i]);
					}
				}
			}
		}

		r = new FastaSequenceFile[generaKept.size()];
		generaKept.values().toArray(r);

		logger.info("genomes before filtering: " + homologySearchSequenceFiles.length);
		logger.info("genomes after filtering: " + r.length);
		return r;
	}

	private String enhanceHomologGroups(String hgDirName,
			FastaSequenceFile[] inputSequenceFiles, int alignThreads, int hmmThreads, 
			int minTaxa, int maxTaxa, FastaSequenceFile[] outgroupGenomeSequencFiles,
			int outgroupGenomesToInclude, String[][] retainedOutgroupHolder) {
		logger.info(">PhyloPipeline.enhanceHomologGroups()");
		String r = hgDirName;
		HMMSetEnhancer hmmSetEnhancer = new HMMSetEnhancer(hgDirName, 
				inputSequenceFiles, alignThreads, hmmThreads, minTaxa, maxTaxa, 
				outgroupGenomeSequencFiles, outgroupGenomesToInclude);
		hmmSetEnhancer.run();
		r = hmmSetEnhancer.getResultDirectoryName();
		retainedOutgroupHolder[0] = hmmSetEnhancer.getRetainedOutgroups();
		logger.info("<PhyloPipeline.enhanceHomologGroups()");
		return r;
	}

	public static void setCommandPaths() {
		String jarPath = PhyloPipeline.class.getProtectionDomain().
				getCodeSource().getLocation().getPath();
		//expect executables to be one level up from this jar, then in bin.
		//i.e. ../bin

		logger.info("jar path: " + jarPath);

		String bin = "bin";
		String os = System.getProperty("os.name");
		logger.info("os: " + os);
		if(os.startsWith("Mac")) {
			bin = "bin_mac";
		} 
		File jarFile = new File(jarPath);
		File jarParentDir = jarFile.getParentFile();
		File jarGrandparentDir = jarParentDir.getParentFile();
		String binDirName = jarGrandparentDir.getAbsolutePath() + 
				"/" + bin + "/";

		logger.info("bin directory: " + binDirName);

		String[] commands = new String[]{
				"muscle",
				"Gblocks",
				"blastall",
				"blat",
				"formatdb",
				"FastTree",
				"FastTree_LG",
				"FastTree_WAG",
				"mcl",
				"raxmlHPC",
				"raxmlHPC-PTHREADS",
				"consel",
				"makermt",
				"catpv",
				"hmmbuild",
				"hmmsearch"
		};

		for(int i = 0; i < commands.length; i++) {
			String fullCommandPath = binDirName + commands[i];
			logger.info("set path for '" + commands[i] + "' to '" + fullCommandPath + "'" );
			System.out.println("set path for '" + commands[i] + "' to '" + fullCommandPath + "'" );
			System.setProperty(commands[i], fullCommandPath);
		}

		//check for specific version of muscle, and use if it is present
		String muscleVersion = "muscle-3.6";
		String fullCommandPath = ExecUtilities.getCommandPath(muscleVersion);
		if(fullCommandPath != null && fullCommandPath.length() > 0) {
			logger.info("set path for '" + muscleVersion + "' to '" + fullCommandPath + "'" );
			System.setProperty("muscle", fullCommandPath);
		}

	}

	private TextFile runMCL(TextFile hitPairFile, String inflation, int threads) {
		TextFile r = null;

		//get the name if the hitPairFile, which is the mcl input file
		String mclInputFileName = hitPairFile.getFile().getAbsolutePath();

		//determine the name of the mcl output file.
		String mclOutputFileName = mclInputFileName + "_mclOut";

		//construct mcl command
		String mclPath = ExecUtilities.getCommandPath("mcl");

		String mclCommand = mclPath + " " + mclInputFileName + " --abc "
				+ " -I " + inflation + " -te " + threads + " -o " +
				mclOutputFileName;

		if(verbose > 0) {
			System.out.println("mcl command: " + mclCommand);
		}
		CommandResults results = ExecUtilities.exec(mclCommand);
		try {
			r = new TextFile(mclOutputFileName);
		} catch (IOException e) {
			e.printStackTrace();
			logger.log(Level.DEBUG,e.getMessage());
		}
		return r;
	}

	private TextFile filterForBidirectional(TextFile hitPairFile, int scoreField) throws IOException {

		TextFile r = null;
		int id1Field = 0;
		int id2Field = 1;
		String tab = "\t";
		String nl = "\n";

		String bidirectionalFileName = 
				hitPairFile.getFile().getAbsolutePath() + "_bidir";
		FileWriter fw = new FileWriter(bidirectionalFileName);

		int lineCount = hitPairFile.getLineCount();
		HashSet allIds = new HashSet();
		if(verbose > 0) {
			logger.info("filterForBidirectional() lines in input file: "
					+ lineCount);
		}

		hitPairFile.openFile();
		Pattern tabPattern = Pattern.compile(tab);
		for(int i = 0; i < lineCount; i++) {
			String[] fields = tabPattern.split(hitPairFile.getLine(i));
			if(fields.length > 1) {
				String id1 = new String(fields[id1Field]);
				String id2 = new String(fields[id2Field]);

				allIds.add(id1);
				allIds.add(id2);
			}
		}

		String[] ids = new String[allIds.size()];
		allIds.toArray(ids);
		allIds = null;
		Arrays.sort(ids);
		HashMap pairToScore = new HashMap();

		for(int i = 0; i < lineCount; i++) {
			String[] fields = tabPattern.split(hitPairFile.getLine(i));
			if(fields.length > scoreField) {
				String id1 = new String(fields[id1Field]);
				String id2 = new String(fields[id2Field]);
				int index1 = Arrays.binarySearch(ids, id1);
				int index2 = Arrays.binarySearch(ids, id2);
				IntPair pair = new IntPair(index1, index2);
				Float storedScore = (Float) pairToScore.get(pair);
				if(storedScore == null) {
					//this pair has not been seen. Add it to the hash
					Float score = new Float(fields[scoreField]);
					pairToScore.put(pair, score);
				} else {
					//this pair has already been stored, so this is a bidirectional
					//hit. Print two lines to the output file for this pair.
					String line1 = id1 + tab + id2
							+ tab + fields[scoreField] + nl;
					String line2 = id2 + tab + id1 
							+ tab + storedScore.toString() + nl;
					fw.write(line1);
					fw.write(line2);
					fw.flush();
					pairToScore.remove(pair);
				}
			}
		}

		hitPairFile.closeFile();

		pairToScore = null;
		fw.close();

		r = new TextFile(bidirectionalFileName);
		if(verbose > 0) {
			logger.info("lines in filtered file: " + r.getLineCount());
		}
		return r;
	}

	private TextFile filterHitPairFile(TextFile hitPairFile) throws IOException {
		TextFile r = null;
		if(verbose > 0) {
			logger.info("PhyloPipeline.filterHitPairFile() " 
					+ hitPairFile.getLineCount() + " lines in file");
		}
		int id1Field = 0;
		int id2Field = 1;
		int scoreField = 11;
		String tab = "\t";
		String nl = "\n";

		int lineCount = hitPairFile.getLineCount();
		String filteredFileName = 
				hitPairFile.getFile().getAbsolutePath() + "_filtered";
		FileWriter fw = new FileWriter(filteredFileName);

		for(int i = 0; i < lineCount; i++) {
			String[] fields = hitPairFile.getLine(i).split(tab);
			if(fields.length > 1) {
				String filteredLine = fields[id1Field] + tab + fields[id2Field]
						+ tab + fields[scoreField] + nl;
				fw.write(filteredLine);
			}
		}

		fw.flush();
		fw.close();

		r = new TextFile(filteredFileName);
		if(verbose > 0) {
			logger.info("Done filtering for columns. " +
					r.getLineCount() + " lines in file");
		}
		return r;
	}

	/**
	 * Runs blast of all-vs-each for the given sequence files.
	 * Returns a TextFile representing the concatenated results of each of
	 * the independent blast runs.
	 * 
	 * @param inputSequenceFiles
	 * @return
	 */
	private TextFile runBlat(FastaSequenceFile[] inputSequenceFiles, 
			int threads, int hitsPerQuery, float evalue, int minIdentity,
			int minScore) {
		TextFile r = null;
		BlatRunner br = new BlatRunner();
		br.setSequenceSets(inputSequenceFiles);
		br.setQuerySequenceFiles(inputSequenceFiles);
		br.setThreadCount(threads);
		br.setHitsPerQuery(hitsPerQuery);
		br.setEvalueCutoff(evalue);
		br.setMinScore(minScore);
		br.setMinIdentity(minIdentity);
		br.setRnName(runName);
		br.run();
		r = br.getResults();

		return r;
	}

	/**
	 * Runs blat of all-vs-each for the given sequence files.
	 * Returns a TextFile representing the concatenated results of each of
	 * the independent blat runs.
	 * 
	 * @param inputSequenceFiles
	 * @return
	 */
	private TextFile runBlast(FastaSequenceFile[] inputSequenceFiles, 
			int threads, int hitsPerQuery, float evalue) {
		System.out.println("PhyloPipeline.runBlast() inputSequenceFiles: " + inputSequenceFiles.length);
		TextFile r = null;
		BlastRunner br = new BlastRunner();
		br.setSequenceSets(inputSequenceFiles);
		br.setQuerySequenceFiles(inputSequenceFiles);
		br.setThreadCount(threads);
		br.setHitsPerQuery(hitsPerQuery);
		br.setEvalueCutoff(evalue);
		br.setExtensionThreshold(18);
		br.setRunName(runName);
		br.run();
		r = br.getResults();

		return r;
	}

	private String[] getDistinctTaxa(FastaSequenceFile[] inputSequenceFiles) {
		HashSet<String> uniqueTaxa = new HashSet<String>();

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

	private FastaSequenceFile[] loadSequenceFiles(String[] sequenceFileNames) throws IOException {
		FastaSequenceFile[] r = new FastaSequenceFile[sequenceFileNames.length];
		for(int i = 0; i < r.length; i++) {
			r[i] = new FastaSequenceFile(sequenceFileNames[i]);
		}
		return r;
	}

	public static String[] getTrackProperties(String track) {
		String[] r = null;
		ArrayList<String> propertyLines = new ArrayList<String>();
		propertyLines.add("-" + HandyConstants.HOMOLOGY_SEARCH_METHOD);
		propertyLines.add(HandyConstants.BLAST);
		propertyLines.add("-" + HandyConstants.BIDIRECTIONAL);
		propertyLines.add(HandyConstants.TRUE);
		propertyLines.add("-" + HandyConstants.CONCATENATED);
		propertyLines.add(HandyConstants.TRUE);
		propertyLines.add("-" + HandyConstants.FULL_TREE_METHOD);
		propertyLines.add(HandyConstants.MAXIMUM_LIKELIHOOD);
		propertyLines.add("-" + HandyConstants.SUPPORT_TREE_METHOD);
		propertyLines.add(HandyConstants.FAST_TREE);
		propertyLines.add("-" + HandyConstants.GENE_WISE_JACKKNIFE);
		propertyLines.add(HandyConstants.TRUE);
		propertyLines.add("-" + HandyConstants.SUPPORT_REPS);
		propertyLines.add("100");
		propertyLines.add("-" + HandyConstants.UNIFORM_TRIM);
		propertyLines.add(HandyConstants.FALSE);
		propertyLines.add("-" + HandyConstants.GBLOCKS_TRIM);
		propertyLines.add(HandyConstants.TRUE);
		propertyLines.add("-" + HandyConstants.PREALIGN);
		propertyLines.add(HandyConstants.TRUE);
		propertyLines.add("-" + HandyConstants.TARGET_MIN_GENE_COUNT);
		propertyLines.add("9999");
		propertyLines.add("-" + HandyConstants.MIN_TAXA_MULTIPLIER);
		propertyLines.add("0.99");
		propertyLines.add("-" + HandyConstants.UNIQUE_SPECIES);
		propertyLines.add(HandyConstants.TRUE);
		propertyLines.add("-" + HandyConstants.CONGRUENCE_FILTER);
		propertyLines.add(HandyConstants.FALSE);
		propertyLines.add("-" + HandyConstants.REFINE);
		propertyLines.add(HandyConstants.TRUE);

		if(track.equals(HandyConstants.TRACK_BLAST_FAST)){
			r = new String[]{
					"-" + HandyConstants.FULL_TREE_METHOD,
					HandyConstants.FAST_TREE,					
					"-" + HandyConstants.FULL_TREE_THREADS,
					"1", 
			};
		}

		r = propertyLines.toArray(new String[0]);
		return r;
	}

	private static void defineHelp() {
		commands = new HashMap<String,String>();

		commands.put(HandyConstants.MIN_TAXA, 
				"Miniumum number of taxa in a sequence set (smaller sets are filtered out). Default value is the number of unique taxa detected.");
		commands.put(HandyConstants.MAX_TAXA, 
				"Maximum number of taxa in a sequence set (larger sets are filtered out)");
		commands.put(HandyConstants.REPRESETATIVE_ONLY, 
				"Only use representative sequence sets (maximum of one member per taxon)");
		commands.put(HandyConstants.ALIGN_THREADS, 
				"Number of threads (processors or processor cores) to use for multiple sequence alignment stage");
		commands.put(HandyConstants.PTHREADS, 
				"Number of threads to use for raxml tree-building stage");
		commands.put(HandyConstants.GENOME_FILE,
				"Name of input sequence set files (any number of file names may follow)");
		commands.put(HandyConstants.OUTGROUP,
				"Name of outgroup sequence set files (any number of file names may follow). Outgroup files are not included in the initial steps of homology group formation, but are added in at a later stage.");
		commands.put(HandyConstants.ALIGN_FILE, 
				"Name of input alignment file (any number of file names may follow)");
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
				"Number of support replicates to perform");
		commands.put(HandyConstants.HELP, 
				"\tPrint this help.");
		commands.put(HandyConstants.DIRECTORY, 
				"\tName of directory containing input faa files. All file from this directory will be loaded. Only applies if -" + HandyConstants.FILE + " is not specified");
		commands.put(HandyConstants.GENE_WISE_JACKKNIFE, 
				"For tree support values, build trees from a subset of the genes used for the full tree. ");
		commands.put(HandyConstants.TREE_THREADS, 
				"Number of processes to use for tree building step(s)");
		commands.put(HandyConstants.FULL_TREE_METHOD, 
				"Method used to build full concatenated tree. Default is Maximum Likelihood (\"" + HandyConstants.MAXIMUM_LIKELIHOOD + "\"). The other options are Parsimony with Maximum Likelihood branch lengths (\"" + HandyConstants.PARSIMONY_BL + "\") and FastTree (\"" + HandyConstants.FAST_TREE + "\")");
		commands.put(HandyConstants.SUPPORT_TREE_METHOD, 
				"Method used for building support trees (trees used for branch support values fo the full tree). Default is FastTree (\"" + HandyConstants.FAST_TREE + "\"). The other options are Maximum Likelihood (\"" + HandyConstants.MAXIMUM_LIKELIHOOD + "\") and Parsimony (\""+ HandyConstants.PARSIMONY + "\").");
		commands.put(HandyConstants.TARGET_NTAX, 
				"The ideal number of taxa per sequence set. By default, this is set to the value for " + HandyConstants.MAX_TAXA);
		commands.put(HandyConstants.TARGET_MIN_GENE_COUNT, 
				"The ideal minimum number of genes (sequence sets) to be used. The value for " + HandyConstants.MIN_TAXA + " takes precedence over this.");
		commands.put(HandyConstants.FULL_TREE_THREADS, 
				"The number of threads to use for building the concatenated tree. This comes out of the treeThreads amount, and is not in addition to it. By default, the full tree_threads amount is used.");
		commands.put(HandyConstants.SINGLE, 
				"\tBuild a tree for each input sequence set meeting the membership criteria.");
		commands.put(HandyConstants.SINGLE_TREE_METHOD, 
				"Method used to build single-gene trees.");
		commands.put(HandyConstants.MATRIX_EVALUATION, 
				"Perform matrix evaluation to determine which amino acid substitution matrix performs best for this dataset.");
		commands.put(HandyConstants.UNIQUE_SPECIES, 
				"For initial homology search step, only use one member from each species. Homology group seed clusters are built from these results and then expanded with hmmsearch.");
		commands.put(HandyConstants.CONGRUENCE_FILTER, 
				"Filter homolog sets based on potential phylopgenetic congruence before concatenating alignments to infer main tree. Congruence filter is based on conflicting ns in the data set.");
		commands.put(HandyConstants.RUN_NAME, 
				"Provide an optional run name for. Output files will contain this name, making it easier to track the results. If no run name is provided, one will be automatically generated.");
		commands.put(HandyConstants.MIN_TAXA_MULTIPLIER, 
				"Proportion of max_taxa to be used as a minimum taxa value. Sequence sets with fewer taxa are not included in tree building. Default is 0.8");
		//		commands.put(HandyConstants.TRACK, "\tSpecify a pre-defined track (set of options) to use. Tracks are identified by the sequence similarity search program used and the tree-building program used. Track options are:\n" +
		//				"\t\t\t\t\t\tblat_fast\n\t\t\t\t\t\tblast_fast\n\t\t\t\t\t\tblat_raxml\n\t\t\t\t\t\tblast_raxml");
	}

	/**
	 * Prints a help message, including a list of available commands.
	 */
	private static void printHelp() {
		logger.info("Printing Help");
		Object[] commandNames = commands.keySet().toArray();

		Arrays.sort(commandNames);

		for(int i = 0; i < commandNames.length; i++) {
			Object commandDesc = commands.get(commandNames[i]);
			String line = "\t-" + commandNames[i] + "\t\t\t" + commandDesc;
			System.out.println(line);
		}

	}
	private void checkForRequiredPrograms() {
		logger.info("Check for required programs...");
		String[] requiredPrograms = new String[]{
				"blastall",
				//				"blat",
				"formatdb",
				"FastTree",
				"FastTree_WAG",
				"Gblocks",
				"hmmbuild",
				"hmmsearch",
				"mcl",
				"muscle",
				"muscle-3.6",
				"raxmlHPC",
				"raxmlHPC-PTHREADS"
		};

		boolean[] programFound = new boolean[requiredPrograms.length];

		for(int i = 0; i < programFound.length; i++) {
			String path = ExecUtilities.getCommandPath(requiredPrograms[i]);
			programFound[i] = path != null;
			if(!programFound[i]) {
				logger.info("The required program '" + requiredPrograms[i] +
						"' was not found");
			} else {
				File programFile = new File(path);
				programFound[i] = programFile.exists();
				if(!programFound[i]) {
					logger.info("The required program '" + requiredPrograms[i] + 
							"' was not found at expected location: " + path);
				} else {
					programFile.setExecutable(true);
					programFound[i] = programFile.canExecute();
					if(programFound[i]) {
						logger.info("The required program '" + requiredPrograms[i] + 
								"' found at: " + path);
					} else {
						logger.info("The required program '" + requiredPrograms[i] + 
								"' was found at '" + path + "', but is not executable.");
					}
				}
			}
		}

		for(int i = 0; i < requiredPrograms.length; i++) {
			logger.info("required program available: " + 
					requiredPrograms[i] + "\t" + programFound[i]);
		}
	}

	private void writePropertiesToLogFile(CommandLineProperties clp) throws IOException {
		if(verbose > 0) {
			System.out.println("PhyloPipeline.writePropertiesToLogFile()");
		}
		String logFileName = runName + ".clp";

		if(verbose > 0) {
			System.out.println("log file: " + logFileName);
		}

		FileWriter fw = new FileWriter(logFileName);
		String[] clpArgs = clp.getArgs();
		for(int i = 0; i < clpArgs.length; i++) {
			fw.write(clpArgs[i] + "\n");
		}
		fw.flush();
		fw.close();
	}

	public void run() {

	}

	private static void checkMCLResults(TextFile mclResultFile) throws IOException {
		int[] fieldCounts = new int[mclResultFile.getLineCount()];
		String tab = "\t";
		for(int i = 0; i < fieldCounts.length; i++) {
			String[] fields = mclResultFile.getLine(i).split(tab);
			fieldCounts[i] = fields.length;
		}

		StatisticsUtilities.printDistribution(fieldCounts);
	}

	private String[] getSelectedOutgroupGenomes() {
		return selectedOutgroupGenomes;
	}

	private void setSelectedOutgroupGenomes(String[] selectedOutgroupGenomes) {
		this.selectedOutgroupGenomes = selectedOutgroupGenomes;
	}
}
