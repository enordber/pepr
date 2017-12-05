package edu.vt.vbi.ci.pepr.alignment;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

import edu.vt.vbi.ci.util.CommandLineProperties;
import edu.vt.vbi.ci.util.CommandResults;
import edu.vt.vbi.ci.util.ExecUtilities;
import edu.vt.vbi.ci.util.HandyConstants;
import edu.vt.vbi.ci.util.file.FastaSequenceFile;
import edu.vt.vbi.ci.util.file.TextFile;

public class BlastRunner {

	private FastaSequenceFile[] sequenceSets;
	private FastaSequenceFile[] querySequenceFiles;
	private FastaSequenceFile concatenatedSequenceSet;
	private FastaSequenceFile[] querySets;
	private TextFile[][] setResults;
	private TextFile allResults;
	private String resultFileName = null;
	private int threadCount = 1;
	private String remoteWorkingDir = null;
	private int hitsPerQuery = 1;
	private double evalueCutoff;
	private String blastd;
	private int wordSize;
	private int extensionThreshold;
	private boolean keepFiles = false;
	private String runName;
	private boolean blastn = false;

	HashSet<String> formattedTargets;
	private int nextTargetIndex = 0;
	private int[] nextTargetIndices;

	private static int defaultThreads = Runtime.getRuntime().availableProcessors();
	private static double defaultEvalueCutoff = 1e-10;
	private static int defaultHitsPerQuery = 1;

	public static void main(String[] args) {

		BlastRunner br = new BlastRunner(args);
		br.run();

	}

	public BlastRunner(String[] args) {
		CommandLineProperties commandLineProperties = 
				new CommandLineProperties(args);

		String[] infileNames = 
				commandLineProperties.getValues(HandyConstants.FILE); 
		if(infileNames == null || infileNames.length == 0) {
			System.out.println("please provide input file names with option '-"
					+ HandyConstants.FILE + "'");
			printHelp();
			System.exit(0);
		}
		FastaSequenceFile[] sequenceFiles =
				new FastaSequenceFile[infileNames.length];
		for(int i = 0; i < sequenceFiles.length; i++) {
			try {
				System.out.println("loading sequence file: " + infileNames[i]);
				sequenceFiles[i] = new FastaSequenceFile(infileNames[i]);
				System.out.println("done loading file with " +
						sequenceFiles[i].getSequenceCount() + " sequences");
			} catch (IOException e) {
				e.printStackTrace();
			}
		}

		//check for separate list of query files
		String[] queryFileNames =
				commandLineProperties.getValues(HandyConstants.QUERY);
		if(queryFileNames == null) {
			queryFileNames = infileNames;
		}
		FastaSequenceFile[] queryFiles = 
				new FastaSequenceFile[queryFileNames.length];
		for(int i = 0; i < queryFiles.length; i++) {
			try {
				queryFiles[i] = new FastaSequenceFile(queryFileNames[i]);
			} catch(IOException e) {
				e.printStackTrace();
			}
		}

		//get name of preformatted blast database, if provided
		String blastd = commandLineProperties.getValues(HandyConstants.BLASTD, null)[0];

		//get blast word size (-W)
		String wordSizeParam = 
				commandLineProperties.getValues(HandyConstants.WORD_SIZE, "0")[0];
		int wordSize = Integer.parseInt(wordSizeParam);

		//get Extension Threshold (-f)
		String etParam = 
				commandLineProperties.getValues(HandyConstants.EXT_THRESH, "11")[0];
		int extThresh = Integer.parseInt(etParam);

		String threadParam = 
				commandLineProperties.getValues("blast_threads", 
						""+defaultThreads)[0];
		int blastThreads = Integer.parseInt(threadParam);

		String hitsPerQueryParam = 
				commandLineProperties.getValues("hits_per_query", 
						""+ defaultHitsPerQuery)[0];
		int blastHitsPerQuery = Integer.parseInt(hitsPerQueryParam);

		String evalueParam = commandLineProperties.getValues("evalue", ""+defaultEvalueCutoff)[0];
		double evalue = Double.parseDouble(evalueParam);
		String resultFileName = 
				commandLineProperties.getValues(HandyConstants.OUT, null)[0];

		boolean blastn = !commandLineProperties.getValues(HandyConstants.BLASTN, HandyConstants.FALSE)[0].equals(HandyConstants.FALSE);

		setSequenceSets(sequenceFiles);
		setQuerySequenceFiles(queryFiles);
		setBlastd(blastd);
		setWordSize(wordSize);
		setExtensionThreshold(extThresh);
		setResultFileName(resultFileName);
		setThreadCount(blastThreads);
		setHitsPerQuery(blastHitsPerQuery);
		setEvalueCutoff(evalue);
		setBlastn(blastn);
	}

	public void setExtensionThreshold(int extThresh) {
		this.extensionThreshold = extThresh;
	}

	private void setWordSize(int w) {
		this.wordSize = w;
	}

	private void setBlastd(String blastd) {
		this.blastd = blastd;
	}

	private void setBlastn(boolean bn) {
		blastn = bn;
	}

	private static void printHelp() {
		System.out.println("options:");
		System.out.println("-"+ HandyConstants.FILE + "\t\tinput file name(s)");
		System.out.println("-blast_threads\tnumber of parallel process to run." +
				" Default: " + defaultThreads);
		System.out.println("-hits_per_query" +
				"\tmaximum number of hits to record for each query. " +
				"Default: " + defaultHitsPerQuery);
		System.out.println("-evalue \te-value cutoff. Default: " + 
				defaultEvalueCutoff);
	}
	public BlastRunner() {

	}

	public void setSequenceSets(FastaSequenceFile[] sequenceSets)  {
		this.sequenceSets = sequenceSets;
	}

	public void setQuerySequenceFiles(FastaSequenceFile[] queries) {
		this.querySequenceFiles = queries;
	}

	public void setEvalueCutoff(double evalue) {
		evalueCutoff = evalue;
	}

	/**
	 * Sets the number of individual processes to used for performing blast.
	 * Individual runs are not split up, but multiple runs can be performed
	 * in parallel.
	 */
	public void setThreadCount(int threads) {
		this.threadCount = threads;
	}

	public void setHitsPerQuery(int hits) {
		hitsPerQuery = hits;
	}

	private void createConcatenatedSet() {
		try {
			File tempFile = File.createTempFile("cat", ".faa", new File(System.getProperty("user.dir")));
			BufferedWriter bw = new BufferedWriter( new FileWriter(tempFile));
			if(!keepFiles) {
				tempFile.deleteOnExit();
			}
			//write each line from each sequenceSet to the concatenated file
			for(int i = 0; i < querySequenceFiles.length; i++) {
				int fileLineCount = querySequenceFiles[i].getLineCount();
				querySequenceFiles[i].openFile();
				for(int j = 0; j < fileLineCount; j++) {
					String line = querySequenceFiles[i].getLine(j);
					if(line.length() > 1) {
						bw.write(line);
						bw.write("\n");
					}
				}
				querySequenceFiles[i].closeFile();
			}
			bw.flush();
			bw.close();

			concatenatedSequenceSet = new FastaSequenceFile(tempFile.getName());
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	/**
	 * Creates the concatenated files to be used as blast queries. The
	 * number of files created depends upon the number of input sequence
	 * files and the number of threads to be used. If the number of input
	 * files is greater than the number of threads, then only a single
	 * concatenated file is created. If a larger number of threads is to
	 * be used, then the concatenated query file is split into enough
	 * files to fully use the number of threads. That is:
	 * (number of query files) * (number of target files) >= (number of threads)
	 */
	private void createConcatenatedQuerySets() {
		int numberOfQueryFiles = (int) Math.ceil((double)threadCount/querySequenceFiles.length);
		//System.out.println("BlatRunner.createConcatenatedQuerySets() number " +
		//		"of query files to create: " + numberOfQueryFiles);

		//if only one query set is needed, use the concatenatedSequenceSet
		if(numberOfQueryFiles == 1) {
			querySets = new FastaSequenceFile[]{concatenatedSequenceSet};
		} else {
			querySets = new FastaSequenceFile[numberOfQueryFiles];
			//determine start and stop sequences for each query file
			int sequencesPerFile = 
					(int) Math.ceil((double)concatenatedSequenceSet.
							getSequenceCount()/ numberOfQueryFiles);
			int[][] startStops = new int[numberOfQueryFiles][2];
			for(int i = 0; i < startStops.length; i++) {
				startStops[i][0] = sequencesPerFile * i;
				startStops[i][1] = sequencesPerFile * (i+1);
			}
			startStops[startStops.length-1][1] =
					concatenatedSequenceSet.getSequenceCount();

			try {
				File directory = new File(System.getProperty("user.dir"));
				concatenatedSequenceSet.openFile();
				for(int i = 0; i < querySets.length; i++) {
					File f = File.createTempFile("cat", ".faa", directory);
					f.deleteOnExit();
					FileWriter fw = new FileWriter(f);
					for(int j = startStops[i][0]; j < startStops[i][1]; j++) {
						String[] seq = concatenatedSequenceSet.getSequence(j);
						for(int k = 0; k < seq.length; k++) {
							fw.write(seq[k]);
							fw.write("\n");
						}
					}
					fw.flush();
					fw.close();
					querySets[i] = new FastaSequenceFile(f.getAbsolutePath());

				}
				concatenatedSequenceSet.closeFile();

			} catch(IOException ioe) {
				ioe.printStackTrace();
			}
		}
		nextTargetIndices = new int[querySets.length];
	}	

	/**
	 * Creates the concatenated files to be used as blast queries. The
	 * number of files created depends upon the number of input sequence
	 * files and the number of threads to be used. If the number of input
	 * files is greater than the number of threads, then only a single
	 * concatenated file is created. If a larger number of threads is to
	 * be used, then the concatenated query file is split into enough
	 * files to fully use the number of threads. That is:
	 * (number of query files) * (number of target files) >= (number of threads)
	 * 
	 * This method splits into files attempting to have roughly equal total
	 * sequence length in each file, as opposed to equal numbers of sequences.
	 * This should allow a more balanced load for each Blast run.
	 */
	private void createBalancedConcatenatedQuerySets() {
		int numberOfQueryFiles = (int) Math.ceil((double)threadCount/querySequenceFiles.length);

		//if only one query set is needed, use the concatenatedSequenceSet
		if(numberOfQueryFiles == 1) {
			querySets = new FastaSequenceFile[]{concatenatedSequenceSet};
		} else {
			querySets = new FastaSequenceFile[numberOfQueryFiles];
			//determine start and stop sequences for each query file
			int linesPerFile = 
					(int) Math.ceil((double)concatenatedSequenceSet.
							getLineCount()/ numberOfQueryFiles);

			int sequenceCount = concatenatedSequenceSet.getSequenceCount();
			int nextSeqIndex = 0;
			int overCount = 0;
			try {
				File directory = new File(System.getProperty("user.dir"));
				concatenatedSequenceSet.openFile();
				for(int i = 0; i < querySets.length; i++) {
					File f = File.createTempFile("cat", ".faa", directory);
					if(!keepFiles) {
						f.deleteOnExit();
					}
					FileWriter fw = new FileWriter(f);
					//add sequences to this file until it has at 
					//least linesPerFile lines
					int linesInFile = 0;
					while(linesInFile < linesPerFile && nextSeqIndex < sequenceCount) {
						String[] seq = 
								concatenatedSequenceSet.getSequence(nextSeqIndex);
						boolean addNextSequence = false;
						if(linesInFile + seq.length < linesPerFile) {
							addNextSequence = true;
						} else if(overCount > 0) {
							addNextSequence = false;
							overCount--;
							break; //break from while and go to next query file
						} else {
							addNextSequence = true;
							overCount++;
						}
						if(addNextSequence) {
							nextSeqIndex++;
							linesInFile += seq.length;
							for(int k = 0; k < seq.length; k++) {
								fw.write(seq[k]);
								fw.write("\n");
							}
						}
					}

					//make sure any remaining sequences get put in the last file
					if(i == querySets.length -1) {
						//this is the last file - add any remaining sequences
						while(nextSeqIndex < sequenceCount) {
							String[] seq = 
									concatenatedSequenceSet.getSequence(nextSeqIndex);
							nextSeqIndex++;
							linesInFile += seq.length;
							for(int k = 0; k < seq.length; k++) {
								fw.write(seq[k]);
								fw.write("\n");
							}
						}
					}
					fw.flush();
					fw.close();
					querySets[i] = new FastaSequenceFile(f.getAbsolutePath());
				}
				concatenatedSequenceSet.closeFile();

			} catch(IOException ioe) {
				ioe.printStackTrace();
			}
		}
		nextTargetIndices = new int[querySets.length];
	}

	public void run() {
		//create the concatenated sequence set file.
		System.out.println("run() createConcatenatedSet:");
		createConcatenatedSet();
		//		createConcatenatedQuerySets();
		System.out.println("createBalancedConcatenatedQuerySets:");
		createBalancedConcatenatedQuerySets();
		//for now, this assumes a single RemoteHost is being used for all blasts

		String workingDirectory = System.getProperty("user.dir");

		formattedTargets = new HashSet<String>();;

		//create and start the appropriate number of worker Threads
		Thread[] threads = new Thread[threadCount];
		System.out.println("create blast runner threads: " + threadCount);
		for(int i = 0; i < threads.length; i++) {
			int querySetIndex = i%querySets.length;
			FastaSequenceFile queryForThread = querySets[querySetIndex];
			Runnable br = null;
			if(blastn) {
				br = new BlastnRunnable(queryForThread.getFileName(), querySetIndex, 
						blastd, 0, 4);				
			} else {
				br = new BlastRunnable(queryForThread.getFileName(), querySetIndex, 
						blastd, wordSize, extensionThreshold);
			}
			threads[i] = new Thread(br);
			threads[i].start();
		}

		//wait for the Threads to finish
		for(int i = 0; i < threads.length; i++) {
			try {
				threads[i].join();
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

		TextFile result = getResults();
	}

	public void setRemoteWorkingDirectory(String dir) {
		remoteWorkingDir = dir;
	}

	private void setResultFileName(String resultFileName) {
		this.resultFileName = resultFileName;
	}

	private File getResultFile() throws IOException {
		File r = null;
		if(resultFileName == null) {
			if(runName == null) {
				r = File.createTempFile("all_blast_", ".out", new File(System.getProperty("user.dir")));
			} else {
				r = new File(runName + "_all_blast.out");
			}
		} else {
			r = new File(resultFileName);
		}
		return r;
	}

	public TextFile getResults() {
		TextFile r = null;
		if(allResults == null) {
			try {
				File allResultsFile = getResultFile();
				System.out.println("combined result file: " + allResultsFile.getAbsolutePath());
				FileWriter fw = new FileWriter(allResultsFile);

				for (int i = 0; i < setResults.length; i++) {
					for (int j = 0; j < setResults[i].length; j++) {
						if(setResults[i][j] != null) {
							setResults[i][j].openFile();
							int lineCount = setResults[i][j].getLineCount();
							for (int k = 0; k < lineCount; k++) {
								String line = setResults[i][j].getLine(k);
								if (line.length() > 1) {
									fw.write(line);
									fw.write("\n");
								}
							}
							setResults[i][j].closeFile();

							if(!keepFiles) {
								setResults[i][j].getFile().delete();
							}
						}
					}
				}
				fw.flush();
				fw.close();
				allResults = new TextFile(allResultsFile.getAbsolutePath());
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		r = allResults;
		return r;
	}

	private void setBlastOutput(TextFile resultFile, int querySetIndex, 
			int targetIndex) {
		if(setResults == null) {
			synchronized (this) {
				setResults = new TextFile[querySets.length][sequenceSets.length];
			}
		}

		setResults[querySetIndex][targetIndex] = resultFile;
	}

	private synchronized int getNextTargetIndex() {
		int r = nextTargetIndex;
		nextTargetIndex++;
		return r;
	}

	private synchronized int getNextTargetIndex(int querySetIndex) {
		int r = nextTargetIndices[querySetIndex];
		nextTargetIndices[querySetIndex]++;
		return r;
	}

	private class BlastRunnable implements Runnable {

		private String formatdbPath;
		private String blastallPath;
		private String concatenatedSequenceName;
		private String outputFormat = " -m 8 ";

		//when word size of 0 is given, blast uses the default word size
		//this is 3 for blastp
		private int wordSize = 0;
		private int querySetIndex;
		//host index may change after multiple RemoteHost support is added
		private int hostIndex = 0;
		private boolean usePreformattedDB = false;
		private String blastd = "";
		private int extensionThreshold;


		public BlastRunnable(String concatenatedSequenceName, 
				int querySetIndex, String blastd, int wordSize, int extThresh) {
			this.concatenatedSequenceName = concatenatedSequenceName;
			this.querySetIndex = querySetIndex;
			System.out.println("new BlastRunnable " + querySetIndex + ", "
					+ concatenatedSequenceName);
			setBlastd(blastd);
			this.wordSize = wordSize;
			this.extensionThreshold = extThresh;
		}

		private void setBlastd(String bd) {
			blastd = bd;
			if(bd != null && bd.length() > 0) {
				usePreformattedDB = true;
			}
		}

		public void run() {

			System.out.println("BlastRunner.BlastRunnable.run() sequenceSets: " + sequenceSets.length);
			blastallPath = ExecUtilities.getCommandPath("blastall");
			formatdbPath = ExecUtilities.getCommandPath("formatdb");

			try {
				int index = getNextTargetIndex(querySetIndex);
				while (index < sequenceSets.length) {
					synchronized(sequenceSets[index].getName()) {
						//synchronize on the sequenceSet name so only one thread will 
						//run formatdb on this sequence set. By the time the second thread 
						//gets a lock on the name, it will have been formatted by the first thread
						if(!formattedTargets.contains(sequenceSets[index].getName())) {
							//format sequence set with formatdb, if not already done 
							String formatdbCommand = formatdbPath + " -i "
									+ sequenceSets[index].getFullName();
							System.out.println("formatdb command: " + formatdbCommand);
							CommandResults formatdbResult = ExecUtilities.exec(formatdbCommand);
							String[] res = formatdbResult.getStderr();
							for(int i = 0; i < res.length; i++) {
								System.out.println(res[i]);
							}
							formattedTargets.add(sequenceSets[index].getName());
						}
					}

					File forOutfileName = File.createTempFile("blast", ".out");
					String outputName = forOutfileName.getName();
					forOutfileName.delete();

					String dbName = sequenceSets[index].getFullName();
					if(usePreformattedDB) {
						dbName = blastd;
					}
					String blastCommand = blastallPath + " -p blastp " +
							" -d " + dbName + 
							" -i " + concatenatedSequenceName +
							" -e " + evalueCutoff + 
							" -W " + wordSize + 
							" -f " + extensionThreshold +
							" -v " + hitsPerQuery + " -b " +
							hitsPerQuery + outputFormat + " -o " + outputName;

					System.out.println("blastCommand: " + blastCommand);
					ExecUtilities.exec(blastCommand);

					TextFile resultFile = null;
					try {
						resultFile = new TextFile(outputName); 							
						setBlastOutput(resultFile, querySetIndex, index);

					} catch(IOException ioe) {
						System.out.println("failed to get " + outputName);
					}
					index = getNextTargetIndex(querySetIndex);
				}
			} catch (IOException ioe) {
				ioe.printStackTrace();
			}
		}

	}

	private class BlastnRunnable implements Runnable {

		private String formatdbPath;
		private String blastallPath;
		private String concatenatedSequenceName;
		private String outputFormat = " -F F -m 8 ";

		//when word size of 0 is given, blast uses the default word size
		//this is 3 for blastp
		private int wordSize = 0;
		private int querySetIndex;
		//host index may change after multiple RemoteHost support is added
		private boolean usePreformattedDB = false;
		private String blastd = "";
		private int extensionThreshold;


		public BlastnRunnable(String concatenatedSequenceName, 
				int querySetIndex, String blastd, int wordSize, int extThresh) {
			this.concatenatedSequenceName = concatenatedSequenceName;
			this.querySetIndex = querySetIndex;
			setBlastd(blastd);
			this.wordSize = wordSize;
			this.extensionThreshold = extThresh;
		}

		private void setBlastd(String bd) {
			blastd = bd;
			if(bd != null && bd.length() > 0) {
				usePreformattedDB = true;
			}
		}

		public void run() {

			blastallPath = ExecUtilities.getCommandPath("blastall");
			formatdbPath = ExecUtilities.getCommandPath("formatdb");

			try {
				int index = getNextTargetIndex(querySetIndex);
				while (index < sequenceSets.length) {
					
					synchronized(sequenceSets[index].getName()) {
						//synchronize on the sequenceSet name so only one thread will 
						//run formatdb on this sequence set. By the time the second thread 
						//gets a lock on the name, it will have been formatted by the first thread
						if(!formattedTargets.contains(sequenceSets[index].getName())) {							//format sequence set with formatdb, if not already done 
							String formatdbCommand = formatdbPath + " -p F -i "
									+ sequenceSets[index].getFullName();
							System.out.println("formatdb command: " + formatdbCommand);
							CommandResults formatdbResult = ExecUtilities.exec(formatdbCommand);
							String[] res = formatdbResult.getStderr();
							for(int i = 0; i < res.length; i++) {
								System.out.println(res[i]);
							}
							formattedTargets.add(sequenceSets[index].getName());
						}
					}
					
					File forOutfileName = File.createTempFile("blast", ".out");
					String outputName = forOutfileName.getName();

					String dbName = sequenceSets[index].getFullName();
					if(usePreformattedDB) {
						dbName = blastd;
					}
					String blastCommand = blastallPath + " -p blastn " +
							" -d " + dbName + 
							" -i " + concatenatedSequenceName +
							" -e " + evalueCutoff + 
							" -W " + wordSize + 
							" -f " + extensionThreshold +
							" -v " + hitsPerQuery + " -b " +
							hitsPerQuery + outputFormat + " -o " + outputName;

					System.out.println("blastCommand: " + blastCommand);
					CommandResults blastCR =   ExecUtilities.exec(blastCommand);
					String[] stdout = blastCR.getStdout();
					String[] stderr = blastCR.getStderr();

					for(int i = 0; i < stdout.length;i++) {
						System.out.println(stdout[i]);
					}

					for(int i = 0; i < stderr.length;i++) {
						System.out.println(stderr[i]);
					}

					TextFile resultFile = null;
					try {
						resultFile = new TextFile(outputName); 							
						setBlastOutput(resultFile, querySetIndex, index);

					} catch(IOException ioe) {
						System.out.println("failed to get " + outputName);
					}
					index = getNextTargetIndex(querySetIndex);
				}
			} catch (IOException ioe) {
				ioe.printStackTrace();
			}
		}

	}

	public void setRunName(String name) {
		runName = name;
	}


}
