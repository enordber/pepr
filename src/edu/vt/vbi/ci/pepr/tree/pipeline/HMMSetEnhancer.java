package edu.vt.vbi.ci.pepr.tree.pipeline;

import java.io.File;
import java.io.FileFilter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

import org.apache.log4j.Logger;

import edu.vt.vbi.ci.pepr.alignment.MultipleSequenceAligner;
import edu.vt.vbi.ci.pepr.alignment.SequenceAlignment;
import edu.vt.vbi.ci.util.CommandLineProperties;
import edu.vt.vbi.ci.util.ExecUtilities;
import edu.vt.vbi.ci.util.HandyConstants;
import edu.vt.vbi.ci.util.SequenceSetExtractor;
import edu.vt.vbi.ci.util.file.FastaSequenceFile;
import edu.vt.vbi.ci.util.file.TextFile;

public class HMMSetEnhancer {

	private int alignThreads;
	private int hmmThreads;
	private FastaSequenceFile[] genomeSequenceFiles;
	private FastaSequenceFile[] outgroupGenomeSequenceFiles;
	private int outgroupsToInclude;
	private String[] retainedOutgroups;
	private FastaSequenceFile[] initialHGFiles;
	private String[] hmmFileNames;
	private String inputDirName;
	private String outputDirName;
	private int currentSet = -1; //used by hmmbuild threads
	private int currentGenome = -1; //used by hmmsearch threads
	private int minTaxa;
	private int maxTaxa;
	private Logger logger = Logger.getLogger(getClass());

	public static void main(String[] args) {
		PhyloPipeline.setCommandPaths();
		CommandLineProperties clp = new CommandLineProperties(args);
		String hgDirName = clp.getValues(HandyConstants.DIRECTORY)[0];
		String[] genomeFiles = clp.getValues(HandyConstants.GENOME_FILE);
		int threads = Integer.parseInt(clp.getValues(HandyConstants.MAX_CONCURRENT_PROCESS_PARAM, "1")[0]);
		int minTaxa = Integer.parseInt(clp.getValues(HandyConstants.MIN_TAXA, "0")[0]);
		int maxTaxa = Integer.parseInt(clp.getValues(HandyConstants.MAX_TAXA, "999999999")[0]);
		int hmmThreads = threads;

		FastaSequenceFile[] genomeSequenceFiles = 
				new FastaSequenceFile[genomeFiles.length];
		for(int i = 0; i < genomeSequenceFiles.length; i++) {
			try {
				genomeSequenceFiles[i] = new FastaSequenceFile(genomeFiles[i]);
			} catch (IOException e) {
				e.printStackTrace();
			}
		}

		HMMSetEnhancer hmmse = new HMMSetEnhancer(hgDirName, genomeSequenceFiles, 
				threads, hmmThreads, minTaxa, maxTaxa, new FastaSequenceFile[0], 0);
		hmmse.run();
	}

	public HMMSetEnhancer(String hgDirName, FastaSequenceFile[] ingroupGenomeSequenceFiles, 
			int alignThreads, int hmmThreads, int minTaxa, int maxTaxa, 
			FastaSequenceFile[] outgroupGenomeSequenceFiles, int outGroupsToInclude) {
		logger.info("HMMSetEnhancer sequence files: " + 
				ingroupGenomeSequenceFiles.length + " minTaxa: " + minTaxa 
				+ " maxTaxa: " + maxTaxa + 
				" outgroup files: " + outgroupGenomeSequenceFiles.length + 
				" outgroupsToInclude: " + outGroupsToInclude + ". alignThreads: "
				+ alignThreads + " hmmThreads: " + hmmThreads);
		this.alignThreads = alignThreads;
		this.hmmThreads = hmmThreads;
		this.genomeSequenceFiles = ingroupGenomeSequenceFiles;
		this.outgroupGenomeSequenceFiles = outgroupGenomeSequenceFiles;
		this.outgroupsToInclude = outGroupsToInclude;
		inputDirName = hgDirName;
		outputDirName = "hmm_" + hgDirName;
		this.minTaxa = minTaxa;
		this.maxTaxa = maxTaxa;
	}

	public void run() {
		//load initialHGFiles
		try {
			initialHGFiles = loadSequenceFilesFromDirectory(inputDirName);
			//sort initialHGFiles so larger files come first
			Arrays.sort(initialHGFiles);

			hmmFileNames = new String[initialHGFiles.length];
			//create worker threads to create hmms
			Thread[] hmmbuildThreads = new Thread[alignThreads];
			logger.info("Building HMMs for " + hmmFileNames.length 
					+ " sequences sets using "+ hmmbuildThreads.length + 
					" threads.");
			for(int i = 0; i < hmmbuildThreads.length; i++) {
				hmmbuildThreads[i] = new Thread(new HMMBuildRunner());
				hmmbuildThreads[i].setName("Build Thread " + i);
				hmmbuildThreads[i].start();
			}

			//wait for worker threads to finish
			for(int i = 0; i < hmmbuildThreads.length; i++) {
				try {
					hmmbuildThreads[i].join();
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}

			logger.info("Done building HMMs");

			//concatenate all hmms into a single hmm file
			String hmmFileName = null;
			File hmmFile = null;
			FileWriter hmmWriter = null;
			hmmFile = File.createTempFile("hmm_", ".hmm", new File(System.getProperty("user.dir")));
			hmmFileName = hmmFile.getAbsolutePath();
			hmmWriter = new FileWriter(hmmFile);

			for(int i = 0; i < hmmFileNames.length; i++) {
				try{
					TextFile tf = new TextFile(hmmFileNames[i]);
					hmmWriter.write(tf.toString());
					tf.getFile().delete(); //done with individual hmm file
				} catch(IOException ioe) {
					//problem loading this file. //for now, just skip it
				}
			}
			hmmWriter.flush();
			hmmWriter.close();

			//prepend outgroupGenomeSequenceFiles to genomesequenceFiles
			FastaSequenceFile[] allGenomes = new FastaSequenceFile[outgroupGenomeSequenceFiles.length + genomeSequenceFiles.length];
			System.arraycopy(outgroupGenomeSequenceFiles, 0, allGenomes, 0, outgroupGenomeSequenceFiles.length);
			System.arraycopy(genomeSequenceFiles, 0, allGenomes, outgroupGenomeSequenceFiles.length, genomeSequenceFiles.length);
			genomeSequenceFiles = allGenomes;

			//run hmmsearch on each input genome file and collect result
			TextFile[] hmmsearchResultFiles = new TextFile[genomeSequenceFiles.length];
			HMMResult[][] hmmResultHolder = new HMMResult[genomeSequenceFiles.length][];
			//create hmmsearch threads
			Thread[] hmmsearchThreads = new Thread[hmmThreads];
			for(int i = 0; i < hmmsearchThreads.length; i++) {
				hmmsearchThreads[i] = 
						new Thread(new HMMSearchRunner(hmmFileName, 
								hmmsearchResultFiles, hmmResultHolder));
				hmmsearchThreads[i].start();
			}

			//wait for hmmsearch threads to finish
			for(int i = 0; i < hmmsearchThreads.length; i++) {
				try{
					hmmsearchThreads[i].join();
				} catch(InterruptedException ie) {
					ie.printStackTrace();
				}
			}

			logger.info("hmmsearch results have been collected");
			//get hmm score sum for each genome
			double[] hmmScoreSums = new double[hmmResultHolder.length];
			for(int i = 0; i < hmmScoreSums.length; i++) {
				hmmScoreSums[i] = 0;
				for(int j = 0; j < hmmResultHolder[i].length; j++) {
					hmmScoreSums[i] += hmmResultHolder[i][j].getScore();
				}
				logger.info("hmm score sum for genome " + i + " " + 
						genomeSequenceFiles[i].getTaxa()[0] + ": "
						+ hmmScoreSums[i]);
				if(hmmScoreSums[i] < 1) {
					System.out.println("no hmm scores for genome " + i + ", " + genomeSequenceFiles[i].getTaxa()[0]);
					System.out.println("\thmmResultHolder[i].length: " + hmmResultHolder[i].length);
					for(HMMResult result: hmmResultHolder[i]) {
						System.out.println("\t" + result);
					}
				}
			}
			//determine which outgroup taxa to keep 
			boolean[] keepOutgroup = new boolean[outgroupGenomeSequenceFiles.length];
			ArrayList retainedOutgroupTaxa = new ArrayList();
			for(int i = 0; i < outgroupsToInclude; i++) {
				double maxScore = 0;
				int maxScoreIndex = -1;
				for(int j = 0; j < outgroupGenomeSequenceFiles.length; j++) {
					if(!keepOutgroup[j]) {
						if(hmmScoreSums[j] > maxScore) {
							maxScore = hmmScoreSums[i];
							maxScoreIndex = j;
						}
					}
				}
				if(maxScoreIndex > -1) {
					keepOutgroup[maxScoreIndex] = true;
					retainedOutgroupTaxa.add(
							genomeSequenceFiles[maxScoreIndex].getTaxa()[0]);
				}
			}

			String[] retainedOutgroups = new String[retainedOutgroupTaxa.size()];
			retainedOutgroupTaxa.toArray(retainedOutgroups);
			setRetainedOutgroups(retainedOutgroups);

			//print out the outgroups being kept
			logger.info("selected outgroup genomes: " + outgroupsToInclude);
			for(int i = 0; i < keepOutgroup.length; i++) {
				if(keepOutgroup[i]) {
					logger.info("keeping outgroup genome " + i + ": " 
							+ genomeSequenceFiles[i].getTaxa()[0] + 
							" hmm score sum: " + hmmScoreSums[i]);
				}
			}

			for(int i = 0; i < keepOutgroup.length; i++) {
				if(!keepOutgroup[i]) {
					//this will skip over this genome, because it will have no 
					//hmm results
					hmmResultHolder[i] = new HMMResult[0];
				}
			}

			HashMap setNameToMemberList = new HashMap();

			//combine hmmResultHolder arrays into a single HMMResult[]
			int resultCount = 0;
			for(int i = 0; i < hmmResultHolder.length; i++) {
				resultCount += hmmResultHolder[i].length;
			}

			HMMResult[] hmmResults = new HMMResult[resultCount];
			int startIndex = 0;
			for(int i = 0; i < hmmResultHolder.length; i++) {
				System.arraycopy(hmmResultHolder[i], 0, hmmResults, startIndex, hmmResultHolder[i].length);
				startIndex += hmmResultHolder[i].length;
				hmmResultHolder[i] = null;
			}

			Arrays.sort(hmmResults);

			//collect entries for each set until a duplicate from a 
			//genome is found
			String previousSetName = "";
			if(hmmResults.length > 0){
				hmmResults[0].getSet();
			}
			ArrayList workingSet = new ArrayList();
			boolean[] genomeSeen = new boolean[genomeSequenceFiles.length];
			float[] genomeScore = new float[genomeSequenceFiles.length];
			boolean collectingForSet = true;
			for(int i = 0; i < hmmResults.length; i++) {
				if(hmmResults[i].getSet() != previousSetName) {
					//finish up previous set and prepare for the next set
					setNameToMemberList.put(previousSetName, workingSet);

					workingSet = new ArrayList();
					Arrays.fill(genomeSeen, false);
					Arrays.fill(genomeScore, 0);
					collectingForSet = true;
				}

				int genomeIndex = hmmResults[i].getGenomeIndex();
				if(genomeSeen[genomeIndex] || !collectingForSet) {
					if(hmmResults[i].getScore() == genomeScore[genomeIndex]) {
						//if this hit for this genome has the same score as the
						//previous hit for this genome, just skip it.
						//This is either an exact duplicate gene, which will
						//not alter the tree at all, or there is a problem in the
						//sequence file with duplicate entries. The duplicate
						//entry files cause the homolog sets to be prematurely
						//truncated and this is to correct for that problem.	
						logger.info("ignoring duplicate member from genome " 
								+ genomeSequenceFiles[genomeIndex].getFileName() + ". hmmResults[" + i + "]:.getId(): " +
								hmmResults[i].getId() + " with score: " + hmmResults[i].getScore() + " for set " + hmmResults[i].getSet());
					} else {
						//duplicate genome seen - stop collecting for this set
						collectingForSet = false;
					}
				} else {
					workingSet.add(hmmResults[i].getId());
				}
				genomeSeen[genomeIndex] = true;
				genomeScore[genomeIndex] = hmmResults[i].getScore();
				previousSetName = hmmResults[i].getSet();
			}

			//print set definitions to a set file
			File setFile = File.createTempFile("sets", ".out", new File(System.getProperty("user.dir")));
			FileWriter setWriter = new FileWriter(setFile);
			String[] sets = new String[setNameToMemberList.size()];
			setNameToMemberList.keySet().toArray(sets);
			Arrays.sort(sets);
			String tab = "\t";
			for(int i = 0; i < sets.length; i++) {
				ArrayList memberList = (ArrayList)setNameToMemberList.get(sets[i]);
				String[] members = new String[memberList.size()];
				memberList.toArray(members);
				for(int j = 0; j < members.length; j++) {
					setWriter.write(members[j]);
					if(j+1 < members.length) {
						setWriter.write(tab);
					}
				}
				setWriter.write("\n");
			}
			setWriter.flush();
			setWriter.close();


			logger.info("creating new homolog set files in " + outputDirName);
			//create directory for new homolog set sequence files
			File outputDir = new File(outputDirName);
			if(!outputDir.exists()) {
				outputDir.mkdir();
			}

			//create new homolog set sequence files
			SequenceSetExtractor sse = 
					new SequenceSetExtractor(setFile.getAbsolutePath(),
							genomeSequenceFiles, outputDirName + "/set", "faa");
			logger.info("done creating new homolog set files");

		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	private FastaSequenceFile[] loadSequenceFilesFromDirectory(String dirName) throws IOException {
		File dir = new File(dirName);
		FileFilter faaFilter = new FileFilter(){
			public boolean accept(File pathname) {
				return pathname.getName().endsWith("faa");
			}
		};
		File[] files = dir.listFiles(faaFilter);

		ArrayList fileList = new ArrayList(files.length);

		for(int i = 0; i < files.length; i++) {
			if(files[i].exists() && files[i].length() > 0) {
				FastaSequenceFile fsf = new FastaSequenceFile(files[i].getAbsolutePath());
				if(fsf.getSequenceCount() >= minTaxa && fsf.getSequenceCount() <= maxTaxa) {
					fileList.add(fsf);
				}
			}
		}
		initialHGFiles = new FastaSequenceFile[fileList.size()];
		fileList.toArray(initialHGFiles);

		return initialHGFiles;
	}

	public String getResultDirectoryName() {
		return outputDirName;
	}

	private synchronized int getNextSetIndex() {
		int r = ++currentSet;
		if(r >= initialHGFiles.length) {
			r = -1;
		}
		System.out.print(".");
		if(r >= 0 && (r+1) % 10 == 0) {
			System.out.print(" ");
			if((r+1 )%100 == 0) {
				System.out.println("\t" + (r+1));
			}
		}

		return r;
	}

	private synchronized int getNextGenomeIndex() {
		int r = ++currentGenome;

		return r;
	}

	private int getGenomeIndexForID(String id) {
		int r = -1;
		for(int i =0; i < genomeSequenceFiles.length; i++) {
			int indexInFile = genomeSequenceFiles[i].getIndexOfSequence(id);
			if(indexInFile >= 0) {
				r = i;
				break;
			}
		}
		return r;
	}


	private void setRetainedOutgroups(String[] outgroups) {
		this.retainedOutgroups = outgroups;
	}

	public String[] getRetainedOutgroups() {
		return this.retainedOutgroups;
	}

	private class HMMResult implements Comparable{
		private String id;
		private String set;
		float score;
		int genomeIndex;

		public HMMResult(String resultLine) {
			//parse the result line and extract the important parts
			String space = "\\s+";
			int idField = 0;
			int setField = 2;
			int scoreField = 5;

			String[] fields = resultLine.split(space);
			this.id = fields[idField].intern();
			this.set = fields[setField].intern();
			this.score = Float.parseFloat(fields[scoreField]);
			this.genomeIndex = getGenomeIndexForID(this.id);
		}

		public int compareTo(Object arg0) {
			int r = 0;
			HMMResult o = (HMMResult) arg0;
			r = this.set.compareTo(o.set);
			if(r == 0) {
				r = (int) (o.score - this.score);
			}
			return r;
		}

		public String toString() {
			String r = null;
			r = set + "\t" + genomeIndex + "\t" + score + "\t" + id;
			return r;
		}

		public String getId() {
			return id;
		}

		public String getSet() {
			return set;
		}

		public float getScore() {
			return score;
		}

		public int getGenomeIndex() {
			return genomeIndex;
		}
	}

	private class HMMBuildRunner implements Runnable {
		public void run() {
			String hmmbuildPath = ExecUtilities.getCommandPath("hmmbuild");

			//get sequence set file to work on
			int setIndex = getNextSetIndex();
			while(setIndex >= 0 && setIndex < initialHGFiles.length) {
				//align sequence set
				MultipleSequenceAligner msa = new MultipleSequenceAligner();
				SequenceAlignment alignment = 
						msa.getMSA(initialHGFiles[setIndex]);

				//run hmmbuild on alignment to create an hmm
				//create alignment file
				//
				String alignmentFileName = alignment.getName();
				try {
					FileWriter fw = new FileWriter(alignmentFileName);
					fw.write(alignment.getAlignmentAsFasta());
					fw.flush();
					fw.close();
				} catch (IOException e) {
					e.printStackTrace();
				}

				//determine name of hmm file
				String hmmFileName = alignmentFileName + ".hmm";

				//call hmmbuild
				//hmmer calls fasta "afa"
				String hmmbuildCmd = hmmbuildPath + " --informat afa " +
						hmmFileName + " " + alignmentFileName;

				ExecUtilities.exec(hmmbuildCmd);

				//delete alignment now
				new File(alignmentFileName).delete();

				hmmFileNames[setIndex] = hmmFileName;
				setIndex = getNextSetIndex();
			}
		}
	}

	private class HMMSearchRunner implements Runnable {

		private String hmmFileName;
		private TextFile[] resultHolder;
		private HMMResult[][] hmmResultHolder; 
		boolean keepHMMSearchFiles = false;

		public HMMSearchRunner(String hmmfileName, TextFile[] resultHolder,
				HMMResult[][] hmmResultHolder) {
			this.hmmFileName = hmmfileName;
			this.resultHolder = resultHolder;
			this.hmmResultHolder = hmmResultHolder;
		}

		public void run() {
			int genomeIndex = getNextGenomeIndex();
			while(genomeIndex < genomeSequenceFiles.length) {
				try {
					String outfileName = "";
					File f = File.createTempFile("hmmsearch_" + genomeIndex + "_", ".out", new File(System.getProperty("user.dir")));
					outfileName = f.getAbsolutePath();
					String hmmsearchPath = ExecUtilities.getCommandPath("hmmsearch");

					// "-o /dev/null" redirects main output so it doesn't go 
					//to stdout. I am using the tblout results, so I don't
					//need the stdout output
					String hmmsearchCmd = hmmsearchPath + " --tblout " + outfileName +
							" -o /dev/null --noali -E 1e-40 --cpu 2 " 
							+ hmmFileName + " " 
							+ genomeSequenceFiles[genomeIndex].getFile().getPath();

					ExecUtilities.exec(hmmsearchCmd);
					TextFile resultFile =  new TextFile(outfileName);
					if(!keepHMMSearchFiles) {
						resultFile.getFile().deleteOnExit();
					}
					resultHolder[genomeIndex] = resultFile;
					HMMResult[] hmmResults = getHMMResults(resultFile);
					hmmResultHolder[genomeIndex] = hmmResults;
					genomeIndex = getNextGenomeIndex();
				} catch(IOException ioe) {
					logger.error("IOException on genomeIndex " + 
							genomeIndex + ": " + genomeSequenceFiles[genomeIndex].getFileName());
					logger.error(ioe);
				}
			}
		}

		private HMMResult[] getHMMResults(TextFile resultFile) throws IOException {
			HMMResult[] r = null;
			String[] resultLines = resultFile.getAllLines();
			System.out.println("HMMSetEnhancer.getHMMResults() resultsFile: " + resultFile.getFile().getName() + " lines: " + resultLines.length);
			ArrayList<HMMResult> hmmResults = new ArrayList<HMMResult>();
			for(int i = 0; i < resultLines.length; i++) {
				HMMResult hmmResult = parseHMMResult(resultLines[i]);
				if(hmmResult != null) {
					hmmResults.add(hmmResult);
				}
			}

			r = new HMMResult[hmmResults.size()];
			hmmResults.toArray(r);
			System.out.println("results from " + resultFile.getFile().getName() + ": " + r.length);
			return r;
		}

		private HMMResult parseHMMResult(String resultLine) {
			HMMResult r = null;
			if(resultLine.length() > 1 && !resultLine.startsWith("#")) {
				try {
					r = new HMMResult(resultLine);
				} catch(Exception e) {
					System.out.println("HMMSetEnhancer.parseHMMResult() problem parsing hmmresult line:\n" + resultLine + " :: message: " + e.getMessage());
					e.printStackTrace();
				}
			}
			return r;
		}

	}
}
