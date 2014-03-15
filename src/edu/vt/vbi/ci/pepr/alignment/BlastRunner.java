package edu.vt.vbi.ci.pepr.alignment;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import edu.vt.vbi.ci.util.CommandLineProperties;
import edu.vt.vbi.ci.util.CommandResults;
import edu.vt.vbi.ci.util.HandyConstants;
import edu.vt.vbi.ci.util.RemoteHost;
import edu.vt.vbi.ci.util.file.FastaSequenceFile;
import edu.vt.vbi.ci.util.file.TextFile;

public class BlastRunner {

	private FastaSequenceFile[] sequenceSets;
	private FastaSequenceFile[] querySequenceFiles;
	private FastaSequenceFile concatenatedSequenceSet;
	private FastaSequenceFile[] querySets;
	private RemoteHost host = RemoteHost.getLocalHost();
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
	
	//blast params
//	private String 

	//useSW determines if smith-waterman alignment scores should be calculated
	//for all hit pairs before filtering for the top hits
	private boolean useSW = false;

	//first dimension is index of the RemoteHost. second dimension is
	//index of the sequenceSet. if true, it means this sequence set has
	//been formatted using formatdb on the host
	boolean targetFormattedOnHost[][];
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
				// TODO Auto-generated catch block
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

		//
		//
		RemoteHost remoteHost = RemoteHost.getLocalHost();

		String hostName = 
			commandLineProperties.getValues(HandyConstants.HOST,
			"localhost")[0];
		if(!hostName.equals("localhost")) {

			String user = commandLineProperties.getValues(HandyConstants.USER, 
					System.getProperty("user.name"))[0];
			String keyPath = System.getProperty("user.home") + "/.ssh/id_rsa";
			keyPath = 
				commandLineProperties.getValues(HandyConstants.RSA_PATH_PARAM, 
						keyPath)[0];

			remoteHost = new RemoteHost(hostName, user, keyPath);
			remoteHost.disableStrictHostKeyChecking();
			System.out.println("using remote host: " + hostName);
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

		setRemoteHost(remoteHost);
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

	public void setRemoteHost(RemoteHost host) {
		this.host = host;
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
			File tempFile = File.createTempFile("cat", ".faa");
			String concatenatedFileName = tempFile.getName();
			FileWriter fw = new FileWriter(concatenatedFileName);
			tempFile.deleteOnExit();
			//write each line from each sequenceSet to the concatenated file
			for(int i = 0; i < querySequenceFiles.length; i++) {
				int fileLineCount = querySequenceFiles[i].getLineCount();
				querySequenceFiles[i].openFile();
				for(int j = 0; j < fileLineCount; j++) {
					String line = querySequenceFiles[i].getLine(j);
					if(line .length() > 1) {
						fw.write(line);
						fw.write("\n");
					}
				}
				querySequenceFiles[i].closeFile();
			}
			fw.flush();
			fw.close();

			concatenatedSequenceSet = new FastaSequenceFile(concatenatedFileName);
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
//		System.out.println("query sets: ");
//		for(int i = 0; i < querySets.length; i++) {
//			System.out.println(querySets[i].getFileName());
//		}

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
//		System.out.println("BlastRunner.createBalancedConcatenatedQuerySets() number " +
//				"of query files to create: " + numberOfQueryFiles);

		//if only one query set is needed, use the concatenatedSequenceSet
		if(numberOfQueryFiles == 1) {
			querySets = new FastaSequenceFile[]{concatenatedSequenceSet};
		} else {
			querySets = new FastaSequenceFile[numberOfQueryFiles];
			//determine start and stop sequences for each query file
			int linesPerFile = 
				(int) Math.ceil((double)concatenatedSequenceSet.
						getLineCount()/ numberOfQueryFiles);

//			System.out.println("lines per file >= " + linesPerFile);
			int sequenceCount = concatenatedSequenceSet.getSequenceCount();
			int nextSeqIndex = 0;
			int overCount = 0;
			try {
				File directory = new File(System.getProperty("user.dir"));
				concatenatedSequenceSet.openFile();
				for(int i = 0; i < querySets.length; i++) {
					File f = File.createTempFile("cat", ".faa", directory);
					f.deleteOnExit();
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
//							System.out.println("under: " + overCount);
							overCount--;
							break; //break from while and go to next query file
						} else {
							addNextSequence = true;
//							System.out.println("over: " + overCount);
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
//					System.out.println("lines in file " + i + ": " + linesInFile);
					querySets[i] = new FastaSequenceFile(f.getAbsolutePath());

				}
				concatenatedSequenceSet.closeFile();

			} catch(IOException ioe) {
				ioe.printStackTrace();
			}
		}
		nextTargetIndices = new int[querySets.length];
//		System.out.println("query sets: ");
//		for(int i = 0; i < querySets.length; i++) {
//			System.out.println(querySets[i].getFileName());
//		}
//
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
		//if host is not local host, then transfer the files
		//to the host
		if(!host.isLocalHost()) {
			if(remoteWorkingDir == null) {
				String pwdCommand = host.getCommandPath("pwd");
				remoteWorkingDir =
					host.executeCommand(pwdCommand).getStdout()[0];
			}

			host.setRemoteWorkingDirectory(remoteWorkingDir);
			for(int i = 0; i < sequenceSets.length; i++) {
				String localFileName = sequenceSets[i].getFullName();
				String remoteFileName = sequenceSets[i].getFileName();
				host.copyToRemoteWorkingDirectory(localFileName, remoteFileName);
				workingDirectory = remoteWorkingDir;
			}

			for(int i = 0; i < querySets.length; i++) {
				String localFileName = querySets[i].getFullName();
				String remoteFileName = querySets[i].getFileName();
				host.copyToRemoteWorkingDirectory(localFileName, remoteFileName);

			}
		}

		//create the tracking array for formatdb
		targetFormattedOnHost = new boolean[1][sequenceSets.length];

		//create and start the appropriate number of worker Threads
		Thread[] threads = new Thread[threadCount];
		System.out.println("create blast runner threads: " + threadCount);
		for(int i = 0; i < threads.length; i++) {
			int querySetIndex = i%querySets.length;
			FastaSequenceFile queryForThread = querySets[querySetIndex];
			Runnable br = null;
			if(blastn) {
				br = new BlastnRunnable(host,
						queryForThread.getFileName(), querySetIndex, 
						blastd, 0, 4);				
			} else {
			br = new BlastRunnable(host,
					queryForThread.getFileName(), querySetIndex, 
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

		private RemoteHost host;
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


		public BlastRunnable(RemoteHost host, String concatenatedSequenceName, 
				int querySetIndex, String blastd, int wordSize, int extThresh) {
			this.host = host;
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

			blastallPath = host.getCommandPath("blastall");
			formatdbPath = host.getCommandPath("formatdb");

			try {
				int index = getNextTargetIndex(querySetIndex);
				while (index < sequenceSets.length) {
					synchronized(host) {
						if(!targetFormattedOnHost[hostIndex][index]) {
							//format sequence set with formatdb, if not already done 
							String formatdbCommand = formatdbPath + " -i "
							+ sequenceSets[index].getFullName();
							System.out.println("formatdb command: " + formatdbCommand);
							CommandResults formatdbResult = host.executeCommand(formatdbCommand);
							String[] res = formatdbResult.getStderr();
							for(int i = 0; i < res.length; i++) {
								System.out.println(res[i]);
							}
							targetFormattedOnHost[hostIndex][index] = true;
						}
					}
					File forOutfileName = File.createTempFile("blast", ".out");
					String outputName = forOutfileName.getName();

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
					host.executeCommand(blastCommand);
					host.copyFromRemoteWorkingDirectory(outputName, outputName);

					TextFile resultFile = null;
					try {
						if(useSW) {
							//pre-filter to limit the number of pairs per query
//							resultFile = new TextFile(outputName); 							
//							outputName = filterForTopHits(resultFile, 10);
							resultFile = 
								getSWScoresForHits(outputName, concatenatedSequenceName);
						} else {
							resultFile = new TextFile(outputName); 							
						}
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

		private RemoteHost host;
		private String formatdbPath;
		private String blastallPath;
		private String concatenatedSequenceName;
		private String outputFormat = " -F F -m 8 ";

		//when word size of 0 is given, blast uses the default word size
		//this is 3 for blastp
		private int wordSize = 0;
		private int querySetIndex;
		//host index may change after multiple RemoteHost support is added
		private int hostIndex = 0;
		private boolean usePreformattedDB = false;
		private String blastd = "";
		private int extensionThreshold;


		public BlastnRunnable(RemoteHost host, String concatenatedSequenceName, 
				int querySetIndex, String blastd, int wordSize, int extThresh) {
			this.host = host;
			this.concatenatedSequenceName = concatenatedSequenceName;
			this.querySetIndex = querySetIndex;
//			System.out.println("new BlastRunnable " + querySetIndex + ", "
//					+ concatenatedSequenceName);
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

			blastallPath = host.getCommandPath("blastall");
			formatdbPath = host.getCommandPath("formatdb");

			try {
				int index = getNextTargetIndex(querySetIndex);
				while (index < sequenceSets.length) {
					synchronized(host) {
						if(!targetFormattedOnHost[hostIndex][index]) {
							//format sequence set with formatdb, if not already done 
							String formatdbCommand = formatdbPath + " -p F -i "
							+ sequenceSets[index].getFullName();
							System.out.println("formatdb command: " + formatdbCommand);
							CommandResults formatdbResult = host.executeCommand(formatdbCommand);
							String[] res = formatdbResult.getStderr();
							for(int i = 0; i < res.length; i++) {
								System.out.println(res[i]);
							}
							targetFormattedOnHost[hostIndex][index] = true;
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
					CommandResults blastCR =   host.executeCommand(blastCommand);
					String[] stdout = blastCR.getStdout();
					String[] stderr = blastCR.getStderr();
					
					for(int i = 0; i < stdout.length;i++) {
						System.out.println(stdout[i]);
					}
					
					for(int i = 0; i < stderr.length;i++) {
						System.out.println(stderr[i]);
					}
					
					host.copyFromRemoteWorkingDirectory(outputName, outputName);

					TextFile resultFile = null;
					try {
						if(useSW) {
							//pre-filter to limit the number of pairs per query
//							resultFile = new TextFile(outputName); 							
//							outputName = filterForTopHits(resultFile, 10);
							resultFile = 
								getSWScoresForHits(outputName, concatenatedSequenceName);
						} else {
							resultFile = new TextFile(outputName); 							
						}
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
	private TextFile getSWScoresForHits(String pairFileName, String seqFileName) {
		TextFile r = null;
//		System.out.println(">getSWScoresForHits() " + 
//				pairFileName + " " + seqFileName + " " + 
//				Thread.currentThread().getName());
		String outFileName = pairFileName + ".sw";
		System.out.println("outFileName: " + outFileName);
		String[] alignerCommands = new String[]{
				"-" + HandyConstants.FILE,
				pairFileName,
				"-" + HandyConstants.SEQ_FILE_NAME,
				seqFileName, 
				"-" + HandyConstants.SCORE_COL,
				"11",
				"-" + HandyConstants.M8,
				outFileName,
				"-" + HandyConstants.BANDED,
				HandyConstants.TRUE, 
				"-" + HandyConstants.BAND_WIDTH,
				"32"
		};
		AAPairAligner aligner = new AAPairAligner(alignerCommands);
		aligner.run();
		try {
			r = new TextFile(outFileName);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
//		System.out.println("<getSWScoresForHits() " + Thread.currentThread().getName());
		return r;
	}

	public void setRunName(String name) {
		runName = name;
	}


}
