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
	private int mashSketchSize = 100000;

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
	 * These distances are symmetrical, and the map will contain both directions, genomeA --> genomeB, genomeB --> genomeA
	 */
	private HashMap<String,HashMap<String,Double>> genomeToGenomeDistances = new HashMap<String,HashMap<String,Double>>();

	private String[] ingroupTaxa;
	private HashMap<String,String> genomeFileNameToTaxonName = new HashMap<String,String>();
	private HashMap<String,String> taxonNameToGenomeFileName = new HashMap<String,String>();

	public static void main(String[] args) {
		FastaUtilities.setStripPipeAndSuffix(false);
		NeighborMasher nm = new NeighborMasher();
		if(nm.readParameters(args)) {
			nm.run();
		} else {
			nm.printHelp();
		}
	}

	private void run() {
		System.out.println(">NeighborMasher.run()");
		//run mash ingroup vs ingroup
		try {
			if(expandedIngroupSize > getIngroupGenomeFileNames().length) {
//				expandIngroup(expandedIngroupSize);
				slowlyExpandIngroup(expandedIngroupSize);
				System.out.println("build tree with expanded ingroup...");
			} 
			runIngroupVsIngroupMash(); //doing this a second time can be avoided by mashing new ingroup members vs each other
		} catch (IOException e1) {
			System.out.println("There was a problem trying to run mash on ingroup genomes:");
			e1.printStackTrace();
		}
		//if(!outgroup_only) build unrooted ingroup tree
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

			//build outgroup sketch file if necessary (if it doesn’t already exist)
			//run mash outgroup vs ingroup
			GenomeCandidate[] outgroupCandidates = getMinDistanceForEachOutgroupToIngroupSketch();
			//select outgroup genomes
			GenomeCandidate[] selectedOutgroupCandidates = selectOutgroupCandidates(outgroupCandidates, maximumIngroupDistance, ingroupDistanceStdDev, getOutgroupCount());
			try {
				//write ingroup file
				writeIngroupFile(getRunName() + "_ingroup.txt");
				//write outgroup file
				writeOutgroupFile(selectedOutgroupCandidates, getRunName() + "_outgroup.txt");
			} catch (IOException e) {
				System.out.println("There was a problem writing the outgroup file:");
				e.printStackTrace();
			}
			//if(outgroup_only) exit
			if(!outgroupOnly) {
				//build rooted tree, including outgroup
				try {
					printSomething(selectedOutgroupCandidates);
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

	/**
	 * Expands the ingroup in steps, at most doubling in each step.
	 * This is to try to keep the ingroup balanced and tight.
	 * @param targetSize
	 * @throws IOException
	 */
	private void slowlyExpandIngroup(int targetSize) throws IOException {
		String[] ingroupTaxa = getIngroupTaxa();
		while(ingroupTaxa.length < targetSize) {
			int nextSize = targetSize;
			if(ingroupTaxa.length * 2 < targetSize) {
				nextSize = ingroupTaxa.length * 2;
			}
			System.out.println("NeighborMasher.slowlyExpandIngroup() current size: " + ingroupTaxa.length + " targetSize: " + targetSize + " nextSize: " + nextSize);
			expandIngroup(nextSize);
			addNearestNeighborsToIngroup();
			ingroupTaxa = getIngroupTaxa();
		}
	}

	private void printSomething(GenomeCandidate[] selectedOutgroupCandidates) {
		String[] ingroup = getIngroupTaxa();
		String[] ingroupIds = new String[ingroup.length];
		Pattern pattern = Pattern.compile("_@_");
		for(int i = 0; i < ingroupIds.length; i++) {
			String[] parts = pattern.split(ingroup[i]);
			ingroupIds[i] = parts[1];
		}
		String[] outgroupIds = new String[selectedOutgroupCandidates.length];
		for(int i = 0; i < outgroupIds.length; i++) {
			String[] parts = pattern.split(selectedOutgroupCandidates[i].getTaxonName());
			outgroupIds[i] = parts[1];
		}

		StringBuffer sb = new StringBuffer();
		sb.append("-i ");
		for(String id: ingroupIds) {
			sb.append(id);
			sb.append(",");
		}

		sb.append(" -o ");
		for(String id: outgroupIds) {
			sb.append(id);
			sb.append(",");
		}

		System.out.println(sb);

	}

	private boolean readParameters(String[] args) {
		boolean r = true;
		CommandLineProperties clp = new CommandLineProperties(args);
		boolean printHelp = clp.getValues(HELP, HandyConstants.FALSE)[0].equalsIgnoreCase(HandyConstants.TRUE);
		if(printHelp) {
			r = false;
		} else {
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
			try {
				setMashProcessorCount(Integer.parseInt(clp.getValues(MASH_PROCESSOR_COUNT, ""+mashProcessorCount)[0]));
			} catch(NumberFormatException nfe) {}
			if(ingroupGenomeFileNames == null || ingroupGenomeFileNames.length == 0) {
				System.out.println("No ingroup genomes have been specified. Specify ingroup genome file names with the option '-" + INGROUP + "'");
				r = false;
			}
		}
		return r;
	}

	private void expandIngroup(int targetIngroupSize) throws IOException {
		System.out.println(">NeighborMasher.expandOutgroup() targetSize: " + targetIngroupSize + " initial size: " + getIngroupGenomeFileNames().length);
		//if outgroup sketch was provided, load genome file names from sketch
		String[] outgroupGenomeNamesFromSketch = getGenomeFileNamesFromSketch(getOutgroupSketchFileName());
		runIngroupVsIngroupMash();
		runOutgroupVsIngroupMash();

		GenomeCandidate[] expandedIngroupCandidates = new GenomeCandidate[outgroupGenomeNamesFromSketch.length];
		//for each outgroup genome, calculate the median distance from the ingroup genomes
		for(int i = 0; i < expandedIngroupCandidates.length; i++) {
			String outgroupGenomeFileName = outgroupGenomeNamesFromSketch[i];
			double ingroupDistance = getDistanceFromIngroup(outgroupGenomeFileName);
			expandedIngroupCandidates[i] = new GenomeCandidate(outgroupGenomeFileName, getTaxonName(outgroupGenomeFileName), ingroupDistance);
		}

		Arrays.sort(expandedIngroupCandidates);

		ArrayList<String> expandedIngroupGenomeFileNames = new ArrayList<String>(targetIngroupSize);
		for(String ingroupFileName: getIngroupGenomeFileNames()) {
			expandedIngroupGenomeFileNames.add(ingroupFileName);
		}

		HashSet<String> originalIngroupTaxa = new HashSet<String>();
		for(String taxon: getIngroupTaxa()) {
			originalIngroupTaxa.add(taxon);
		}
		System.out.println("candidates for expanding ingroup:");
		System.out.println("Genome file name\tTaxon name\tDistance from ingroup genomes");
		for(int i = 0; i < expandedIngroupCandidates.length && expandedIngroupGenomeFileNames.size() < targetIngroupSize; i++) {
			double distance = expandedIngroupCandidates[i].getDistance();
			if(distance < 1.0 && !originalIngroupTaxa.contains(expandedIngroupCandidates[i].getTaxonName())) {
				//distance == 1.0 means the genome is too distant / there are no shared kmers, so it's not a good candidate
				expandedIngroupGenomeFileNames.add(expandedIngroupCandidates[i].getGenomeFileName());
				System.out.println(expandedIngroupCandidates[i].getGenomeFileName() + "\t" + expandedIngroupCandidates[i].getTaxonName() + "\t" + expandedIngroupCandidates[i].getDistance());
			}
		}
		String[] newIngroup = expandedIngroupGenomeFileNames.toArray(new String[0]);
		setIngroupGenomeFileNames(newIngroup);
		ingroupTaxa = null; //so taxa will be re-done on next call to getIngroupTaxa()
		System.out.println("<NeighborMasher.expandOutgroup() new size: " + getIngroupGenomeFileNames().length);
	}

	private void addNearestNeighborsToIngroup() {
		System.out.println(">NeighborMasher.addNearestNeighborsToIngroup() " + getIngroupTaxa().length);
		String[] initialIngroup = getIngroupTaxa();
		HashSet<String> expandedIngroup = new HashSet<String>();
		for(String ingroupTaxon: initialIngroup) {
			expandedIngroup.add(ingroupTaxon);
		}
		//of each ingroup taxon, find the nearest genome and add it to the ingroup
		for(String ingroupTaxon: initialIngroup) {
			HashMap<String,Double> genomeDistances = genomeToGenomeDistances.get(ingroupTaxon);
			Double minDistance = Double.POSITIVE_INFINITY;
			String taxonAtMinDistance = null;
			for(String taxon: genomeDistances.keySet()) {
				Double distance = genomeDistances.get(taxon);
				if(distance > 0 && distance < minDistance) {
					minDistance = distance;
					taxonAtMinDistance = taxon;
				}
			}
			if(taxonAtMinDistance != null) {
				if(expandedIngroup.add(taxonAtMinDistance)) {
					System.out.println("adding neighbor of " + ingroupTaxon + ": " + taxonAtMinDistance);					
				}
			}
		}
		ingroupTaxa = expandedIngroup.toArray(new String[0]);
		String[] expandedIngroupGenomeFileNames = new String[ingroupTaxa.length];
		for(int i = 0; i < expandedIngroupGenomeFileNames.length; i++)  {
			expandedIngroupGenomeFileNames[i] = getGenomeFileName(ingroupTaxa[i]);
		}
		setIngroupGenomeFileNames(expandedIngroupGenomeFileNames);
		System.out.println("<NeighborMasher.addNearestNeighborsToIngroup() " + getIngroupTaxa().length);
	}
	
	private double getDistanceFromIngroup(String outgroupGenomeFileName) {
		double r = -1;
		String[] ingroupTaxa = getIngroupTaxa();
		double[] distances = new double[ingroupTaxa.length];
		Arrays.fill(distances, 1.0);
		try {
			String outgroupTaxon = getTaxonName(outgroupGenomeFileName);
			for(int i = 0; i < distances.length; i++) {
				Double distance = getGenomeDistance(outgroupTaxon, ingroupTaxa[i]);
				if(distance != null) {
					distances[i] = distance;
				}
			}
		} catch (IOException e) {
			System.out.println("There was a problem trying to read the taxon name for the genome file, '" + outgroupGenomeFileName + "'");
			e.printStackTrace();
		}

		Arrays.sort(distances);
//		r = distances[distances.length/2];
		r = distances[distances.length-1];
		return r;
	}

	private String[] getGenomeFileNamesFromSketch(String sketchFileName) {
		String[] r = null;
		String mashCommand = getMashPath() + " info " + sketchFileName;
		System.out.println("read outgroup genome file names from sketch file...");
		System.out.println(mashCommand);
		CommandResults results = ExecUtilities.exec(mashCommand);
		String[] stdout = results.getStdout();

		ArrayList<String> fileNames = new ArrayList<String>();
		int expectedFields = 5;
		int genomeFileNameField = 3;
		Pattern delimiter = Pattern.compile("\\s+");
		for(String line: stdout) {
			String[] fields = delimiter.split(line);
			if(fields.length == expectedFields && !fields[genomeFileNameField].startsWith("[")) {
				String genomeFileName = fields[genomeFileNameField];
				fileNames.add(genomeFileName);
			}
		}
		r = fileNames.toArray(new String[0]);
		return r;
	}

	private void writeTreeFile(String tree, String fileName) throws IOException {
		FileWriter fw = new FileWriter(fileName);
		fw.write(tree);
		fw.flush();
		fw.close();
	}

	private void writeOutgroupFile(GenomeCandidate[] selectedOutgroupCandidates, String fileName) throws IOException {
		System.out.println("NeighborMasher.writeOutgroupFile() " + selectedOutgroupCandidates);
		System.out.println(selectedOutgroupCandidates.length);
		FileWriter fw = new FileWriter(fileName);
		for(GenomeCandidate candidate: selectedOutgroupCandidates) {
			String taxonName = candidate.getTaxonName();
			String[] parts = taxonName.split("_@_");
			fw.write(parts[1]);
//			fw.write("\t");
//			fw.write(parts[0]);
//			fw.write("\t");
//			fw.write(candidate.getGenomeFileName());
			fw.write("\n");
		}
		fw.flush();
		fw.close();
	}
	
	private void writeIngroupFile(String fileName) {
		System.out.println("NeighborMasher.writeIngroupFile() " + fileName);
		try {
			FileWriter fw = new FileWriter(fileName);
			String[] ingroupFileNames = getIngroupGenomeFileNames();
			for(String ingroupFileName: ingroupFileNames) {
				FastaSequenceFile fsf = new FastaSequenceFile(ingroupFileName);
				String taxonName = fsf.getTaxa()[0];
				String[] parts = taxonName.split("_@_");
				fw.write(parts[1]);
//				fw.write("\t");
//				fw.write(parts[0]);
//				fw.write("\t");
//				fw.write(ingroupFileName);
				fw.write("\n");
			}
			fw.flush();
			fw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
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
		Double r = 1.0;
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
		//System.out.println(">NeighborMasher.setGenomeDistance() " + genomeA + " vs " + genomeB + ": " + distance);
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

	private GenomeCandidate[] selectOutgroupCandidates(GenomeCandidate[] outgroupCandidates, double maximumIngroupDistance, double ingroupDistanceStdDev, int count) {
		GenomeCandidate[] r = new GenomeCandidate[0];
		Arrays.sort(outgroupCandidates);
		double minimumRequiredDistance = maximumIngroupDistance + (ingroupDistanceStdDev*2);
		System.out.println("minimum required outgroup distance: " + minimumRequiredDistance);
		HashMap<Double,HashSet<GenomeCandidate>> acceptableCandidatesByDistance = new HashMap<Double,HashSet<GenomeCandidate>>();
		HashMap<Double,HashSet<GenomeCandidate>> tooNearCandidatesByDistance = new HashMap<Double,HashSet<GenomeCandidate>>();
		HashSet<GenomeCandidate> tooFarCandidates = new HashSet<GenomeCandidate>();

		for(int i = 0; i < outgroupCandidates.length && acceptableCandidatesByDistance.size() < count; i++) {
			double distance = outgroupCandidates[i].getDistance();
			//if(distance > minimumRequiredDistance && distance < 1.0) {
			if(distance < 1.0) {
				if(distance > minimumRequiredDistance) {
					HashSet<GenomeCandidate> candidatesAtDistance = acceptableCandidatesByDistance.get(distance);
					if(candidatesAtDistance == null) {
						candidatesAtDistance = new HashSet<GenomeCandidate>();
						acceptableCandidatesByDistance.put(distance, candidatesAtDistance);
					}
					candidatesAtDistance.add(outgroupCandidates[i]);
				} else {
					//this candidate to too close to the ingroup. keep it in case there are no good options
					HashSet<GenomeCandidate> candidatesAtDistance = tooNearCandidatesByDistance.get(distance);
					if(candidatesAtDistance == null) {
						candidatesAtDistance = new HashSet<GenomeCandidate>();
						tooNearCandidatesByDistance.put(distance, candidatesAtDistance);
					}
					candidatesAtDistance.add(outgroupCandidates[i]);

				}
			} else {
				tooFarCandidates.add(outgroupCandidates[i]);
			}
		}

		ArrayList<GenomeCandidate> selectedCandidates = new ArrayList<GenomeCandidate>();
		Double[] distanceKeys = acceptableCandidatesByDistance.keySet().toArray(new Double[0]);
		Arrays.sort(distanceKeys);
		//		if(distanceKeys.length == 0 || distanceKeys[0] < minimumRequiredDistance) {
		//			System.out.println("There are no good outgroup candidates with these settings.");
		//			if(distanceKeys.length > 0) {
		//				System.out.println("minimum required distance from ingroup is " + minimumRequiredDistance + " and the best match is at " + distanceKeys[0]);
		//			}
		//		}
		for(Double distance: distanceKeys) {
			HashSet<GenomeCandidate> candidatesAtDistance = acceptableCandidatesByDistance.get(distance);
			for(GenomeCandidate candidate: candidatesAtDistance) {
				selectedCandidates.add(candidate);
				if(selectedCandidates.size() == count) {
					break;
				}
			}
			if(selectedCandidates.size() == count) {
				break;
			}
		}

		if(selectedCandidates.size() == 0) {
			System.out.println("There are no good outgroup candidates with these settings.");
			System.out.println("too far candidates: " + tooFarCandidates.size());
			if(tooNearCandidatesByDistance.size() > 0) {
				System.out.println("There are some genomes that are too near to the ingroup. The farthest of these will be selected.");
				distanceKeys = tooNearCandidatesByDistance.keySet().toArray(new Double[0]);
				Arrays.sort(distanceKeys);
				for(int i = distanceKeys.length-1; i >= 0; i--) {
					HashSet<GenomeCandidate> candidatesAtDistance = tooNearCandidatesByDistance.get(distanceKeys[i]);
					for(GenomeCandidate candidate: candidatesAtDistance) {
						System.out.println(candidate + "\t" + distanceKeys[i]);
						selectedCandidates.add(candidate);
//						if(selectedCandidates.size() == count) {
							break;
//						}
					}
//					if(selectedCandidates.size() == count) {
						break;
//					}					
				}
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

	private String buildRootedTree(GenomeCandidate[] selectedOutgroupGenomes) throws IOException {
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
			String targetGenomeName = getTaxonName(target);
			String query = fields[queryIndex];
			String queryGenomeName = getTaxonName(query);
			Double distance = Double.valueOf(fields[distanceIndex]);
			setGenomeDistance(targetGenomeName, queryGenomeName, distance);
		}
	}

	private void runOutgroupVsIngroupMash() throws IOException {
		System.out.println(">NeighborMasher.runOutgroupVsIngroupMash()");
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
		readMashResults(results);
		System.out.println("<NeighborMasher.runOutgroupVsIngroupMash()");
	}

	private GenomeCandidate[] getMinDistanceForEachOutgroupToIngroupSketch() {
		GenomeCandidate[] r = null;
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

		r = new GenomeCandidate[outgroupCandidateToMinDistance.size()];
		int index = 0;
		for(String genomeFileName: outgroupCandidateToMinDistance.keySet()) {
			Double minDistance = outgroupCandidateToMinDistance.get(genomeFileName);
			try {
				r[index++] = new GenomeCandidate(genomeFileName, getTaxonName(genomeFileName), minDistance);
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
		String[] taxa = getIngroupTaxa();
		r = new double[taxa.length * (taxa.length-1) / 2];
		int index = 0;
		for(int i = 0; i < taxa.length; i++) {
			for(int j = i+1; j < taxa.length; j++) {
				r[index++] = getGenomeDistance(taxa[i], taxa[j]);
//				System.out.println("\t" + taxa[i] + " vs " + taxa[j] + ": " + r[index-1]);
			}
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

	private String[] getIngroupTaxa() {
		if(ingroupTaxa == null) {
			String[] fileNames = getIngroupGenomeFileNames();
			ingroupTaxa = new String[fileNames.length];
			for(int i = 0; i < ingroupTaxa.length; i++) {
				try {
					ingroupTaxa[i] = getTaxonName(fileNames[i]);
				} catch (IOException e) {
					System.out.println("There was a problem trying to get the taxon name from the genome file '" + fileNames[i] + "'");
					e.printStackTrace();
				}
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
		this.mashProcessorCount = Math.min(mashProcessorCount, Runtime.getRuntime().availableProcessors());
	}

	private String getTaxonName(String genomeFileName) throws IOException {
		String r = genomeFileNameToTaxonName.get(genomeFileName);
		if(r == null) {
			r = new FastaSequenceFile(genomeFileName).getTaxa()[0];
			genomeFileNameToTaxonName.put(genomeFileName, r);
			taxonNameToGenomeFileName.put(r, genomeFileName);
		}
		return r;
	}
	
	private String getGenomeFileName(String taxon) {
		String r = null;
		r = taxonNameToGenomeFileName.get(taxon);
		return r;
	}

	private void cleanup() {
		File ingroupSketchFile = new File(getIngroupSketchFileName());
		ingroupSketchFile.delete();

		File selectedOutgroupSketchFile = new File(getSelectedOutgroupSketchFileName());
		selectedOutgroupSketchFile.delete();
	}

	private void printHelp() {
		HashMap<String,String> commandDescriptions = new HashMap<String,String>();
		commandDescriptions.put("ingroup_tree_only","Exit after writing ingroup tree file. This tree can be viewed with midpoint rooting, rather than outgroup. [false]");
		commandDescriptions.put("outgroup_sketch","Name of mash sketch file to use for outgroup. Required unless -outgroup is provided or -ingroup_tree_only.");
		commandDescriptions.put("k","Mash kmer size [" + mashKmerLength + "]");
		commandDescriptions.put("s","Mash sketch size [" + mashSketchSize +"]");
		commandDescriptions.put("outgroup_count","Number of outgroup candidates to report. [" + outgroupCount + "]");
		commandDescriptions.put("run_name","Name for this run. Used as base name for some output files. [timestamp]");
		commandDescriptions.put("outgroup_only","Don’t make any trees. Just output the selected outgroup genomes. [false]");
		commandDescriptions.put("ingroup","List of .fna files for ingroup genomes. Required.");
		commandDescriptions.put("outgroup","List of .fna files for outgroup candidate genomes. Required unless -outgroup_sketch is provided or -ingroup_tree_only.");
		commandDescriptions.put("mash","Path to mash executable. Required.");
		commandDescriptions.put("h","Print this help information, then exit.");
		commandDescriptions.put("p","Number of processors to use. [all available]");
		commandDescriptions.put("expand_ingroup","Add closest genomes from outgroup pool to ingroup. value is target size for ingroup after expansion. Uses same k and s values for mash. [0]");
		//commandDescriptions.put("","");

		ArrayList<String> commands = new ArrayList<String>(commandDescriptions.keySet());
		Collections.sort(commands);
		int commandPrintLength = 17;

		System.out.println("Options:");
		for(String command: commands) {
			String description = commandDescriptions.get(command);
			while(command.length() < commandPrintLength) {
				command = command + " ";
			}
			System.out.println("\t-" + command + "\t" + description);
		}
	}

	private static class GenomeCandidate implements Comparable<GenomeCandidate> {
		private String genomeFileName;
		private String taxonName;
		private double distance;

		private GenomeCandidate(String genomeFileName, String taxonName, double distance) {
			this.genomeFileName = genomeFileName;
			this.taxonName = taxonName;
			this.distance = distance;
		}

		@Override
		public int compareTo(GenomeCandidate o) {
			int r = 0;
			if(o.distance > this.distance) {
				r = -1;
			} else if(o.distance < this.distance) {
				r = 1;
			}
			return r;
		}

		private String getGenomeFileName() {
			return genomeFileName;
		}

		private double getDistance() {
			return distance;
		}

		private String getTaxonName() {
			return taxonName;
		}

		@Override
		public String toString() {
			String r = null;
			StringBuffer sb = new StringBuffer();
			String taxonName = getTaxonName();
			String[] parts = taxonName.split("_@_");
			sb.append(parts[1]);
			sb.append("\t");
			sb.append(parts[0]);
			sb.append("\t");
			sb.append(getGenomeFileName());
			r = sb.toString();
			return r;
		}
	}

}
