package edu.vt.vbi.ci.util;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Pattern;

import edu.vt.vbi.ci.pepr.stats.StatisticsUtilities;
import edu.vt.vbi.ci.pepr.tree.AdvancedTree;
import edu.vt.vbi.ci.pepr.tree.TreeBuilder;
import edu.vt.vbi.ci.util.file.FastaSequenceFile;
import edu.vt.vbi.ci.util.file.FastaUtilities;

public class NeighborMasher {
	/* -ingroup_tree_only */
	private static final String INGROUP_TREE_ONLY = "ingroup_tree_only";
	private boolean ingroupTreeOnly = false;

	/* -outgroup_sketch */
	private static final String OUTGROUP_SKETCH = "outgroup_sketch";
	private String outgroupSketchFileName = null;

	/* -auto_adjust */
	private static final String AUTO_ADJUST = "auto_adjust";
	private boolean autoAdjustParameters = false;

	/* -k */
	private static final String K = "k";
	private int mashKmerLength = 21;

	/* -s */
	private static final String S = "s";
	private int mashSketchSize = 50000;

	/* -outgroup_count */
	private static final String OUTGROUP_COUNT = "outgroup_count";
	private int outgroupCount = 3;

	/* -expand_ingroup */ 
	private static final String EXPAND_INGROUP = "expand_ingroup";
	private int expandedIngroupSize = 0;

	/* -run_name */
	private static final String RUN_NAME = "run_name";
	private String runName = "run_" + System.currentTimeMillis();

	/* -outgroup_only */
	private static final String OUTGROUP_ONLY = "outgroup_only";
	private boolean outgroupOnly = false;

	/* -ingroup */
	private static final String INGROUP = "ingroup";
	private String[] ingroupGenomeFileNames;

	/* -outgroup */
	private static final String OUTGROUP = "outgroup";
	private String[] outgroupGenomeFileNames;

	/* -mash */
	private static final String MASH = "mash";
	private String mashPath = ExecUtilities.getCommandPath(MASH);

	/* -p */
	private static final String MASH_PROCESSOR_COUNT = "p";
	private int mashProcessorCount = Runtime.getRuntime().availableProcessors();
	
	/* -h */
	private static final String HELP = "h";

	/*
	 * Maps from a genome name to a HashMap of other genome names mapped and their distances.
	 * queryToTargetDistances.get(genomeA).get(genomeB) == distance between genomeA and genomeB.
	 * These distances are symmetrical, and the map will contain both orders, genomeA --> genomeB, genomeB --> genomeA
	 */
	private HashMap<String,HashMap<String,Double>> genomeToGenomeDistances = new HashMap<String,HashMap<String,Double>>();

	private String[] ingroupTaxa;
	private HashMap<String,String> genomeFileNameToTaxonName = new HashMap<String,String>();

	public static void main(String[] args) {
		FastaUtilities.setStripPipeAndSuffix(false);
		NeighborMasher nm = new NeighborMasher();
		if(nm.readParameters(args)) {
			nm.run();
		}
	}

	private void run() {
		System.out.println(">NeighborMasher.run()");
		//		run mash ingroup vs ingroup
		try {
			runIngroupVsIngroupMash();
		} catch (IOException e1) {
			System.out.println("There was a problem trying to run mash on ingroup genomes:");
			e1.printStackTrace();
		}
		//		if(!outgroup_only) build unrooted ingroup tree
		if(!isOutgroupOnly()) {
			try {
				String ingroupTree = buildIngroupTree();
				String ingroupTreeFileName = getRunName() + "_k"  + getMashKmerLength() + "_s" + getMashSketchSize() + "_ingroup_nj.nwk";
				writeTreeFile(ingroupTree, ingroupTreeFileName);
			} catch (IOException e) {
				System.out.println("There was a problem trying to build the unrooted ingroup tree:");
				e.printStackTrace();
			}
		}

		//if(ingroup_tree_only) exit
		if(!isIngroupTreeOnly()) {
			//calculate ingroup distance mean and standard deviation, and use them to determine minimum outgroup distance
			double[] ingroupDistances = getIngroupPairDistances();
			double meanIngroupDistance = StatisticsUtilities.getMean(ingroupDistances);
			double ingroupDistanceStdDev = StatisticsUtilities.getStandardDeviation(ingroupDistances, meanIngroupDistance);
			double maximumIngroupDistance = StatisticsUtilities.getMax(ingroupDistances);
			System.out.println("ingroup distance: " + meanIngroupDistance + " +/- " + ingroupDistanceStdDev);
			System.out.println("maximum ingroup pair distance: " + maximumIngroupDistance);

			double minimumRequiredDistance = maximumIngroupDistance + (ingroupDistanceStdDev*2);
			System.out.println("minimum required outgroup distance: " + minimumRequiredDistance);

			//		build outgroup sketch file if necessary (if it doesn’t already exist)
			//		run mash outgroup vs ingroup
			OutgroupCandidate[] outgroupCandidates = getMinDistanceForEachOutgroupToIngroupSketch();
			//		select outgroup genomes
			OutgroupCandidate[] selectedOutgroupCandidates = selectOutgroupCandidates(outgroupCandidates, maximumIngroupDistance, ingroupDistanceStdDev, getOutgroupCount());
			//		write outgroup file
			try {
				writeOutgroupFile(selectedOutgroupCandidates, getRunName() + "_outgroup.txt");
			} catch (IOException e) {
				System.out.println("There was a problem writing the outgroup file:");
				e.printStackTrace();
			}
			//if(outgroup_only) exit
			if(!outgroupOnly) {
			//build rooted tree, including outgroup
				try {
					String rootedTree = buildRootedTree(selectedOutgroupCandidates);
					String rootedTreeFileName = getRunName() + "_k"  + getMashKmerLength() + "_s" + getMashSketchSize() + "_rooted_nj.nwk";
					
					//write rooted tree file
					writeTreeFile(rootedTree, rootedTreeFileName);

				} catch (IOException e) {
					System.out.println("There was a problem building the rooted tree:");
					e.printStackTrace();
				}
			}
			cleanup();
		}
		System.out.println("<NeighborMasher.run()");
	}

	private boolean readParameters(String[] args) {
		boolean r = true;
		CommandLineProperties clp = new CommandLineProperties(args);
		setIngroupTreeOnly(clp.getValues(INGROUP_TREE_ONLY, ""+ingroupTreeOnly)[0].equalsIgnoreCase(HandyConstants.TRUE));
		setOutgroupSketchFileName(clp.getValues(OUTGROUP_SKETCH, null)[0]);
		setAutoAdjustParameters(clp.getValues(AUTO_ADJUST, ""+autoAdjustParameters)[0].equalsIgnoreCase(HandyConstants.TRUE));
		setMashKmerLength(Integer.parseInt(clp.getValues(K, ""+mashKmerLength)[0]));
		setMashSketchSize(Integer.parseInt(clp.getValues(S, ""+mashSketchSize)[0]));
		setOutgroupCount(Integer.parseInt(clp.getValues(OUTGROUP_COUNT, ""+outgroupCount)[0]));
		setExpandedIngroupSize(Integer.parseInt(clp.getValues(EXPAND_INGROUP, ""+expandedIngroupSize)[0]));
		setRunName(clp.getValues(RUN_NAME, runName)[0]);
		setOutgroupOnly(clp.getValues(OUTGROUP_ONLY, ""+outgroupOnly)[0].equalsIgnoreCase(HandyConstants.TRUE));
		setIngroupGenomeFileNames(clp.getValues(INGROUP));
		setOutgroupGenomeFileNames(clp.getValues(OUTGROUP));
		setMashPath(clp.getValues(MASH, mashPath)[0]);
		setMashProcessorCount(Integer.parseInt(clp.getValues(MASH_PROCESSOR_COUNT, ""+mashProcessorCount)[0]));

		boolean printHelp = clp.getValues(HELP, HandyConstants.FALSE)[0].equalsIgnoreCase(HandyConstants.TRUE);
		if(printHelp) {
			printHelp();
			r = false;
		}
		
		if(ingroupGenomeFileNames == null || ingroupGenomeFileNames.length == 0) {
			System.out.println("No ingroup genomes have been specified. Specify ingroup genome file names with the option '-" + INGROUP + "'");
			r = false;
		}
		return r;
	}

	private void writeTreeFile(String tree, String fileName) throws IOException {
		FileWriter fw = new FileWriter(fileName);
		fw.write(tree);
		fw.flush();
		fw.close();
	}
	
	private void writeOutgroupFile(OutgroupCandidate[] selectedOutgroupCandidates, String fileName) throws IOException {
		FileWriter fw = new FileWriter(fileName);
		for(OutgroupCandidate candidate: selectedOutgroupCandidates) {
			String taxonName = candidate.getTaxonName();
			String[] parts = taxonName.split("_@_");
			fw.write(parts[1]);
			fw.write("\t");
			fw.write(parts[0]);
			fw.write("\t");
			fw.write(candidate.getGenomeFileName());
			fw.write("\n");
		}
		fw.flush();
		fw.close();
	}

	/**
	 * returns a previously calculated and cached distance for the two
	 * specified genomes. Returns null if there is no cached distance 
	 * for this pair. This method will not calculate a genome distance.
	 *  
	 * @param genomeA
	 * @param genomeB
	 * @return
	 */
	private Double getGenomeDistance(String genomeA, String genomeB) {
		Double r = null;
		r = genomeToGenomeDistances.get(genomeA).get(genomeB);
		return r;
	}

	/**
	 * Caches the distance for the specified pair of genomes. This distance
	 * can be retrieved by calling getGenomeDistance(genomeA, genomeB).
	 * 
	 * @param genomeA
	 * @param genomeB
	 * @param distance
	 */
	private void setGenomeDistance(String genomeA, String genomeB, Double distance) {
		HashMap<String,Double> genomeAMap = genomeToGenomeDistances.get(genomeA);
		if(genomeAMap == null) {
			genomeAMap = new HashMap<String,Double>();
			genomeToGenomeDistances.put(genomeA, genomeAMap);
		}
		genomeAMap.put(genomeB, distance);

		HashMap<String,Double> genomeBMap = genomeToGenomeDistances.get(genomeB);
		if(genomeBMap == null) {
			genomeBMap = new HashMap<String,Double>();
			genomeToGenomeDistances.put(genomeB, genomeBMap);
		}
		genomeBMap.put(genomeA, distance);
	}

	private OutgroupCandidate[] selectOutgroupCandidates(OutgroupCandidate[] outgroupCandidates, double maximumIngroupDistance, double ingroupDistanceStdDev, int count) {
		OutgroupCandidate[] r = new OutgroupCandidate[count];
		Arrays.sort(outgroupCandidates);
		double minimumRequiredDistance = maximumIngroupDistance + (ingroupDistanceStdDev*2);
		System.out.println("minimum required outgroup distance: " + minimumRequiredDistance);
		HashMap<Double,HashSet<OutgroupCandidate>> acceptableCandidatesByDistance = new HashMap<Double,HashSet<OutgroupCandidate>>();
		for(int i = 0; i < outgroupCandidates.length && acceptableCandidatesByDistance.size() < count; i++) {
			double distance = outgroupCandidates[i].getMinimumIngroupDistance();
			if(distance > minimumRequiredDistance && distance < 1.0) {
				HashSet<OutgroupCandidate> candidatesAtDistance = acceptableCandidatesByDistance.get(distance);
				if(candidatesAtDistance == null) {
					candidatesAtDistance = new HashSet<OutgroupCandidate>();
					acceptableCandidatesByDistance.put(distance, candidatesAtDistance);
				}
				candidatesAtDistance.add(outgroupCandidates[i]);
//				double zScore = (outgroupCandidates[i].getMinimumIngroupDistance() - maximumIngroupDistance) / ingroupDistanceStdDev;
//				try {
//					System.out.print(outgroupCandidates[i].getGenomeFileName() +  " " + getTaxonName(outgroupCandidates[i].getGenomeFileName()) + ":\t" + outgroupCandidates[i].getMinimumIngroupDistance());
//					System.out.println(" * z: " + zScore);
//				} catch (IOException e) {
//					e.printStackTrace();
//				}
			}
		}

		ArrayList<OutgroupCandidate> selectedCandidates = new ArrayList<OutgroupCandidate>();
		Double[] distanceKeys = acceptableCandidatesByDistance.keySet().toArray(new Double[0]);
		Arrays.sort(distanceKeys);
		for(Double distance: distanceKeys) {
			HashSet<OutgroupCandidate> candidatesAtDistance = acceptableCandidatesByDistance.get(distance);
			for(OutgroupCandidate candidate: candidatesAtDistance) {
				selectedCandidates.add(candidate);
				if(selectedCandidates.size() == count) {
					break;
				}
			}
			if(selectedCandidates.size() == count) {
				break;
			}
		}
		
		r = selectedCandidates.toArray(r);
		System.out.println("selectedCandidates (should be " + count + "): " + selectedCandidates.size());
		return r;
	}

	private String getIngroupSketchFileName() {
		String r = "ingroup_sketch_k" + getMashKmerLength() + "_s" + getMashSketchSize() + ".msh";
		return r;
	}

	private String getSelectedOutgroupSketchFileName() {
		String r = "selected_outgroup_sketch_k" + getMashKmerLength() + "_s" + getMashSketchSize() + ".msh";
		return r;
	}

	private String buildIngroupTree() throws IOException {
		String r = null;
		String[] taxa = getIngroupTaxa();
		r = buildNJTreeOfTaxa(taxa);
		return r;
	}
	
	private String buildRootedTree(OutgroupCandidate[] selectedOutgroupGenomes) throws IOException {
		String r = null;
		//calculate distances for each outgroup genome vs each ingroup genome
		
		//create sketch file for selected outgroup genomes
		String[] outgroupFileNames = new String[selectedOutgroupGenomes.length];
		for(int i = 0; i < outgroupFileNames.length; i++) {
			outgroupFileNames[i] = selectedOutgroupGenomes[i].getGenomeFileName();
		}
		String sketchCommand = getMashPath() + " sketch -p " + getMashProcessorCount()
				+ " -s " + getMashSketchSize() + " -k " + getMashKmerLength()
				+ " -o " + getSelectedOutgroupSketchFileName()
				+ " " + getTargetStringStartingAt(outgroupFileNames, 0);
		
		System.out.println("create mash sketch file for selected outgroup genomes...");
		System.out.println(sketchCommand);
		ExecUtilities.exec(sketchCommand);

		//run mash on selected outgroup vs ingroup 
		String mashCommand = getMashPath() + " dist -p " + getMashProcessorCount() + " -s " + getMashSketchSize()
				+ " " + getIngroupSketchFileName() + " " + getSelectedOutgroupSketchFileName();			

		System.out.println("running outgroup vs ingroup mash...");
		System.out.println(mashCommand);
		CommandResults results = ExecUtilities.exec(mashCommand);
		readMashResults(results);
		
		//run mash on outgroup vs outgroup
		mashCommand = getMashPath() + " dist -p " + getMashProcessorCount() + " -s " + getMashSketchSize()
				+ " " + getSelectedOutgroupSketchFileName() + " " + getSelectedOutgroupSketchFileName();			

		System.out.println("running outgroup vs outgroup mash...");
		System.out.println(mashCommand);
		results = ExecUtilities.exec(mashCommand);
		readMashResults(results);

		String[] outgroupTaxa = new String[outgroupFileNames.length];
		String[] ingroupGenomeFileNames = getIngroupGenomeFileNames();
		String[] taxonNames = new String[ingroupGenomeFileNames.length + outgroupFileNames.length];
		int index = 0;
		for(String genomeFileName: outgroupFileNames) {
			outgroupTaxa[index] = getTaxonName(genomeFileName);
			taxonNames[index++] = getTaxonName(genomeFileName);
		}
		for(String genomeFileName: ingroupGenomeFileNames) {
			taxonNames[index++] = getTaxonName(genomeFileName);
		}
		
		String treeString = buildNJTreeOfTaxa(taxonNames);
		
		//actually root the tree
		AdvancedTree tree = new AdvancedTree(treeString);
		tree.setOutGroup(outgroupTaxa);
		r = tree.getTreeString(true, true);

		return r;
	}
	
	/**
	 * Returns a neighbor-joining tree of the specified taxa.
	 * This requires that all genome pair distances have 
	 * already been caluculated and stored, and can be retrieved
	 * from getGenomeDistance().
	 * 
	 * @param taxa
	 * @return
	 */
	private String buildNJTreeOfTaxa(String[] taxa) {
		String r = null;
		double[][] distances = new double[taxa.length][taxa.length];
		for(int i = 0; i < distances.length; i++) {
			distances[i][i] = 0;
			for(int j = i+1; j < distances.length; j++) {
				try {
				double distance = getGenomeDistance(taxa[i], taxa[j]);
				distances[i][j] = distance;
				distances[j][i] = distance;
				}catch(NullPointerException npe) {
					System.out.println("null distance for " + taxa[i] + " vs " + taxa[j]);
				}
			}
		}

		r = TreeBuilder.getNJTreeString(taxa, distances, TreeBuilder.DISTANCE);

		return r;
	}
	
	private void runIngroupVsIngroupMash() throws IOException {
		System.out.println(">NeighborMasher.runIngroupVsIngroupMash()");
		String sketchFileName = getIngroupSketchFileName();
		System.out.println("creating mash sketch file " + sketchFileName + "...");
		String sketchCommand = getMashPath() + " sketch -p " + getMashProcessorCount()
				+ " -s " + getMashSketchSize() + " -k " + getMashKmerLength()
				+ " -o " + sketchFileName
				+ " " + getTargetStringStartingAt(getIngroupGenomeFileNames(), 0);

		System.out.println(sketchCommand);
		ExecUtilities.exec(sketchCommand);

		String mashCommand = getMashPath() + " dist -p " + getMashProcessorCount() + " -s " + getMashSketchSize()
				+ " " + sketchFileName + " " + sketchFileName;			

		System.out.println("running ingroup vs ingroup mash...");
		System.out.println(mashCommand);
		CommandResults results = ExecUtilities.exec(mashCommand);
		readMashResults(results);
		System.out.println("<NeighborMasher.runIngroupVsIngroupMash()");
	}

	private void readMashResults(CommandResults mashResults) throws IOException {
		String[] stdout = mashResults.getStdout();
		Pattern delimiter = Pattern.compile("\\t");
		int targetIndex = 0;
		int queryIndex = 1;
		int distanceIndex = 2;

		for(String line: stdout) {
			String[] fields = delimiter.split(line);
			String target = fields[targetIndex];
			String targetGenomeName = getTaxonName(target); //new FastaSequenceFile(target).getTaxa()[0];
			String query = fields[queryIndex];
			String queryGenomeName = getTaxonName(query); //new FastaSequenceFile(query).getTaxa()[0];
			Double distance = Double.valueOf(fields[distanceIndex]);
			setGenomeDistance(targetGenomeName, queryGenomeName, distance);
		}
	}
	
	private OutgroupCandidate[] getMinDistanceForEachOutgroupToIngroupSketch() {
		OutgroupCandidate[] r = null;
		HashMap<String,Double> outgroupCandidateToMinDistance = new HashMap<String,Double>();
		String sketchFileName = getOutgroupSketchFileName();
		if(sketchFileName == null) {
			//create outgroup sketch file
			createOutgroupSketchFile();
			sketchFileName = getOutgroupSketchFileName();
		}

		String ingroupSketchFileName = getIngroupSketchFileName();
		String mashCommand = getMashPath() + " dist -p " + getMashProcessorCount() + " -s " + getMashSketchSize() 
				+ " " + ingroupSketchFileName + " " + sketchFileName;			

		System.out.println("run mash on outgroup vs ingroup...");
		System.out.println(mashCommand);
		CommandResults results = ExecUtilities.exec(mashCommand);
		String[] stdout = results.getStdout();

		Pattern pattern = Pattern.compile("\\t");
		int queryIndex = 1;
		int distanceFieldIndex = 2;
		for(String line: stdout) {
			String[] fields = pattern.split(line);
			String queryGenome = fields[queryIndex];
			double distance = Double.parseDouble(fields[distanceFieldIndex]);
			Double minDistance = outgroupCandidateToMinDistance.get(queryGenome);
			if(minDistance == null || distance < minDistance) {
				outgroupCandidateToMinDistance.put(queryGenome, distance);
			}			
		}

		r = new OutgroupCandidate[outgroupCandidateToMinDistance.size()];
		int index = 0;
		for(String genomeFileName: outgroupCandidateToMinDistance.keySet()) {
			Double minDistance = outgroupCandidateToMinDistance.get(genomeFileName);
			try {
				r[index++] = new OutgroupCandidate(genomeFileName, getTaxonName(genomeFileName), minDistance);
			} catch (IOException e) {
				System.out.println("There was a problem trying to read the taxon name form the genome file, " + genomeFileName);
				e.printStackTrace();
			}
		}

		return r;
	}

	private void createOutgroupSketchFile() {
		String sketchFileName = "outgroup_sketch_k" + getMashKmerLength() + "_s" + getMashSketchSize() + ".msh";
		System.out.println("creating mash sketch file " + sketchFileName + "...");
		String sketchCommand = mashPath + " sketch -p " + getMashProcessorCount()
				+ " -s " + getMashSketchSize() + " -k " + getMashKmerLength()
				+ " -o " + sketchFileName 
				+ " " + getTargetStringStartingAt(getOutgroupGenomeFileNames(), 0);
		System.out.println(sketchCommand);
		ExecUtilities.exec(sketchCommand);
		setOutgroupSketchFileName(sketchFileName);
	}

	private String getTargetStringStartingAt(String[] fileNames, int start) {
		String r = null;
		StringBuffer sb = new StringBuffer(fileNames[start]);
		for(int i = start+1; i < fileNames.length; i++) {
			sb.append(" ");
			sb.append(fileNames[i]);
		}
		r = sb.toString();
		return r;
	}

	private double[] getIngroupPairDistances() {
		double[] r = null;
		try {
			String[] taxa = getIngroupTaxa();
			r = new double[taxa.length * (taxa.length-1) / 2];
			int index = 0;
			for(int i = 0; i < taxa.length; i++) {
				for(int j = i+1; j < taxa.length; j++) {
					r[index++] = getGenomeDistance(taxa[i], taxa[j]);
				}
			}
		} catch (IOException e) {
			System.out.println("There was a problem reading the taxon name from an ingroup genome file:");
			e.printStackTrace();
		}
		return r;
	}

	private boolean isIngroupTreeOnly() {
		return ingroupTreeOnly;
	}

	private void setIngroupTreeOnly(boolean ingroupTreeOnly) {
		this.ingroupTreeOnly = ingroupTreeOnly;
	}

	private String getOutgroupSketchFileName() {
		return outgroupSketchFileName;
	}

	private void setOutgroupSketchFileName(String outgroupSketchFileName) {
		this.outgroupSketchFileName = outgroupSketchFileName;
	}

	private boolean isAutoAdjustParameters() {
		return autoAdjustParameters;
	}

	private void setAutoAdjustParameters(boolean autoAdjustParameters) {
		this.autoAdjustParameters = autoAdjustParameters;
	}

	private int getMashKmerLength() {
		return mashKmerLength;
	}

	private void setMashKmerLength(int mashKmerLength) {
		this.mashKmerLength = mashKmerLength;
	}

	private int getMashSketchSize() {
		return mashSketchSize;
	}

	private void setMashSketchSize(int mashSketchSize) {
		this.mashSketchSize = mashSketchSize;
	}

	private int getOutgroupCount() {
		return outgroupCount;
	}

	private void setOutgroupCount(int outgroupCount) {
		this.outgroupCount = outgroupCount;
	}

	private int getExpandedIngroupSize() {
		return expandedIngroupSize;
	}

	private void setExpandedIngroupSize(int expandedIngroupSize) {
		this.expandedIngroupSize = expandedIngroupSize;
	}

	private String getRunName() {
		return runName;
	}

	private void setRunName(String runName) {
		this.runName = runName;
	}

	private boolean isOutgroupOnly() {
		return outgroupOnly;
	}

	private void setOutgroupOnly(boolean outgroupOnly) {
		this.outgroupOnly = outgroupOnly;
	}

	private String[] getIngroupGenomeFileNames() {
		return ingroupGenomeFileNames;
	}

	private String[] getIngroupTaxa() throws IOException {
		if(ingroupTaxa == null) {
			String[] fileNames = getIngroupGenomeFileNames();
			ingroupTaxa = new String[fileNames.length];
			for(int i = 0; i < ingroupTaxa.length; i++) {
				ingroupTaxa[i] = getTaxonName(fileNames[i]);
			}
		}
		return ingroupTaxa;
	}

	private void setIngroupGenomeFileNames(String[] ingroupGenomeFileNames) {
		this.ingroupGenomeFileNames = ingroupGenomeFileNames;
	}

	private String[] getOutgroupGenomeFileNames() {
		return outgroupGenomeFileNames;
	}

	private void setOutgroupGenomeFileNames(String[] outgroupGenomeFileNames) {
		this.outgroupGenomeFileNames = outgroupGenomeFileNames;
	}

	private String getMashPath() {
		return mashPath;
	}

	private void setMashPath(String mashPath) {
		this.mashPath = mashPath;
	}

	private int getMashProcessorCount() {
		return mashProcessorCount;
	}

	private void setMashProcessorCount(int mashProcessorCount) {
		this.mashProcessorCount = mashProcessorCount;
	}

	private String getTaxonName(String genomeFileName) throws IOException {
		String r = genomeFileNameToTaxonName.get(genomeFileName);
		if(r == null) {
			r = new FastaSequenceFile(genomeFileName).getTaxa()[0];
			genomeFileNameToTaxonName.put(genomeFileName, r);
		}
		return r;
	}

	private void cleanup() {
		File ingroupSketchFile = new File(getIngroupSketchFileName());
		ingroupSketchFile.delete();
		
		File selectedOutgroupSketchFile = new File(getSelectedOutgroupSketchFileName());
		selectedOutgroupSketchFile.delete();
	}
	
	private static void printHelp() {
		HashMap<String,String> commandDescriptions = new HashMap<String,String>();
		commandDescriptions.put("ingroup_tree_only","Exit after writing ingroup tree file. This tree can be viewed with midpoint rooting, rather than outgroup. [false]");
		commandDescriptions.put("outgroup_sketch","Name of sketch file to use for outgroup.");
		commandDescriptions.put("k","Mash kmer size [21]");
		commandDescriptions.put("s","Mash sketch size [100000]");
		commandDescriptions.put("outgroup_count","Number of outgroup candidates to report. [3]");
		commandDescriptions.put("run_name","Name for this run. Used as base name for some output files. [timestamp]");
		commandDescriptions.put("outgroup_only","Don’t make any trees. Just output the selected outgroup genomes. [false]");
		commandDescriptions.put("ingroup","List of .fna files for ingroup genomes. Required.");
		commandDescriptions.put("outgroup","List of .fna files for outgroup candidate genomes. Required unless -outgroup_sketch is provided or -ingroup_tree_only.");
		commandDescriptions.put("mash","Path to mash executable. Required.");
		commandDescriptions.put("h","Print this help information, then exit.");
		commandDescriptions.put("p","Number of processors to use. [all available]");
//		commandDescriptions.put("","");
//		commandDescriptions.put("","");
				
		ArrayList<String> commands = new ArrayList<String>(commandDescriptions.keySet());
		Collections.sort(commands);
		int commandPrintLength = 17;
		
		for(String command: commands) {
			String description = commandDescriptions.get(command);
			while(command.length() < commandPrintLength) {
				command = command + " ";
			}
			System.out.println("\t-" + command + "\t" + description);
		}
	}
	
	private static class OutgroupCandidate implements Comparable<OutgroupCandidate> {
		private String genomeFileName;
		private String taxonName;
		private double minimumIngroupDistance;

		private OutgroupCandidate(String genomeFileName, String taxonName, double minimumIngroupDistance) {
			this.genomeFileName = genomeFileName;
			this.taxonName = taxonName;
			this.minimumIngroupDistance = minimumIngroupDistance;
		}

		@Override
		public int compareTo(OutgroupCandidate o) {
			int r = 0;
			if(o.minimumIngroupDistance > this.minimumIngroupDistance) {
				r = -1;
			} else if(o.minimumIngroupDistance < this.minimumIngroupDistance) {
				r = 1;
			}
			return r;
		}

		private String getGenomeFileName() {
			return genomeFileName;
		}

		private double getMinimumIngroupDistance() {
			return minimumIngroupDistance;
		}

		private String getTaxonName() {
			return taxonName;
		}
	}

}
