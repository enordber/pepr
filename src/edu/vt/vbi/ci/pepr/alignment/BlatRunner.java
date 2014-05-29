package edu.vt.vbi.ci.pepr.alignment;

import java.io.Console;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.regex.Pattern;

import edu.vt.vbi.ci.util.CommandLineProperties;
<<<<<<< HEAD
import edu.vt.vbi.ci.util.ExecUtilities;
import edu.vt.vbi.ci.util.HandyConstants;
import edu.vt.vbi.ci.util.file.FastaSequenceFile;
import edu.vt.vbi.ci.util.file.TextFile;

public class BlatRunner {

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
	private int minScore = 30;
	private int minIdentity = 25;
	private int stepSize = 1;
	private boolean keepFiles = false;
	private String runName;

	private int nextTargetIndex = 0;
	private int[] nextTargetIndices;

	private static int defaultThreads = 1;
	private static double defaultEvalueCutoff = 1e-10;
	private static int defaultHitsPerQuery = 1;
	private String[] sortedIds;
	
	//useSW determines if smith-waterman alignment scores should be calculated
	//for all hit pairs before filtering for the top hits
	private boolean useSW = false;


	public static void main(String[] args) {
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
				sequenceFiles[i] = new FastaSequenceFile(infileNames[i]);
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

		String minScoreParam =
			commandLineProperties.getValues(HandyConstants.MIN_SCORE, "30")[0];
		int minScore = Integer.parseInt(minScoreParam);

		String minIdentityParam =
			commandLineProperties.getValues(HandyConstants.MIN_IDENTITY, "25")[0];
		int minIdentity = Integer.parseInt(minIdentityParam);

		String resultFileName = 
			commandLineProperties.getValues(HandyConstants.OUT, null)[0];

		BlatRunner br = new BlatRunner();
		br.setSequenceSets(sequenceFiles);
		br.setQuerySequenceFiles(queryFiles);
		br.setThreadCount(blastThreads);
		br.setHitsPerQuery(blastHitsPerQuery);
		br.setEvalueCutoff(evalue);
		br.setMinScore(minScore);
		br.setMinIdentity(minIdentity);
		br.setResultFileName(resultFileName);
		br.run();
		TextFile result = br.getResults();
		System.out.println("Results: " + result.getLineCount());

	}

	private void setResultFileName(String resultFileName) {
		this.resultFileName = resultFileName;
	}

	public void setMinIdentity(int minIdentity) {
		this.minIdentity = minIdentity;
	}

	public void setMinScore(int minScore) {
		this.minScore = minScore;
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
	public BlatRunner() {

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
			File tempFile = File.createTempFile("pepr_seqs", ".faa");
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

			//populate sortedIds array, which is used to provide an index number
			//for each sequence id
			sortedIds = 
				new String[concatenatedSequenceSet.getIDToIndexMap().size()];
			concatenatedSequenceSet.getIDToIndexMap().keySet().toArray(sortedIds);
			Arrays.sort(sortedIds);
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
	 * (number of query files) * (number of input files) >= (number of threads)
	 */
	private void createConcatenatedQuerySets() {
		int numberOfQueryFiles = (int) Math.ceil((double)threadCount/querySequenceFiles.length);
		System.out.println("BlatRunner.createConcatenatedQuerySets() number " +
				"of query files to create: " + numberOfQueryFiles);

		//if only one query set is needed, use the concatenatedSequenceSet
		if(numberOfQueryFiles == 1) {
			querySets = new FastaSequenceFile[]{concatenatedSequenceSet};
		} else {
			querySets = new FastaSequenceFile[numberOfQueryFiles];
			//determine start and stop sequences for each query file
			int sequencesPerFile = 
				concatenatedSequenceSet.getSequenceCount()/ numberOfQueryFiles;
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
		System.out.println("query sets: ");
		for(int i = 0; i < querySets.length; i++) {
			System.out.println(querySets[i].getFileName());
		}

	}

	public void run() {
		//create the concatenated sequence set file.
		createConcatenatedSet();
		createConcatenatedQuerySets();
		//for now, this assumes a single RemoteHost is being used for all blats

		String workingDirectory = System.getProperty("user.dir");

		//create and start the appropriate number of worker Threads
		Thread[] threads = new Thread[threadCount];
		System.out.println("create blat runner threads: " + threadCount);
		for(int i = 0; i < threads.length; i++) {
			int querySetIndex = i%querySets.length;
			FastaSequenceFile queryForThread = querySets[querySetIndex];
			threads[i] = new Thread(new BlatRunnable(
					queryForThread.getFileName(), querySetIndex));
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

	}

	public void setRemoteWorkingDirectory(String dir) {
		remoteWorkingDir = dir;
	}
	
	public void setRnName(String name) {
		runName = name;
	}

	private File getResultFile() throws IOException {
		File r = null;
		if(resultFileName == null) {
			if(runName == null) {
			    r = File.createTempFile("all_blat_", ".out", new File(System.getProperty("user.dir")));
			} else {
				r = new File(runName + "_all_blat.out");
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
					if(setResults[i] != null) {
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

	private void setBlatOutput(TextFile resultFile, int querySetIndex, 
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

	private class BlatRunnable implements Runnable {

		private String blatPath;
		private String concatenatedSequenceName;
		private String outputFormat = " -out=blast8 ";
		private int querySetIndex;
		//host index may change after multiple RemoteHost support is added
		private int hostIndex = 0;


		public BlatRunnable(String concatenatedSequenceName, int querySetIndex) {
			this.concatenatedSequenceName = concatenatedSequenceName;
			this.querySetIndex = querySetIndex;
			System.out.println("new BlatRunnable " + querySetIndex + ", "
					+ concatenatedSequenceName);
		}

		public void run() {

			blatPath = ExecUtilities.getCommandPath("blat");

			try {
				int index = getNextTargetIndex(querySetIndex);
				while (index < sequenceSets.length) {
					System.out.println("blat query set " + querySetIndex 
							+ " against target set " + index + 
							" of " + sequenceSets.length);
					File forOutfileName = File.createTempFile("blat", ".out");
					String outputName = forOutfileName.getName();

					String blatCommand = blatPath + 
					" " + sequenceSets[index].getFullName() + 
					" " + concatenatedSequenceName +
					" -prot " + outputFormat + outputName +
					" -minScore=" + minScore + 
					" -stepSize=" + stepSize + 
					" -minIdentity=" + minIdentity;

					System.out.println("blatCommand: " + blatCommand);
					ExecUtilities.exec(blatCommand);
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
						outputName = filterForTopHits(resultFile, hitsPerQuery);
						resultFile = new TextFile(outputName);
						setBlatOutput(resultFile, querySetIndex, index);
					} catch(IOException ioe) {
						System.out.println("failed to get " + outputName);
					}
					new File(forOutfileName.getName()).delete();
					forOutfileName.delete();
=======
import edu.vt.vbi.ci.util.HandyConstants;
import edu.vt.vbi.ci.util.RemoteHost;
import edu.vt.vbi.ci.util.file.FastaSequenceFile;
import edu.vt.vbi.ci.util.file.TextFile;

public class BlatRunner {

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
	private int minScore = 30;
	private int minIdentity = 25;
	private int stepSize = 1;
	private boolean keepFiles = false;
	private String runName;

	private int nextTargetIndex = 0;
	private int[] nextTargetIndices;

	private static int defaultThreads = 1;
	private static double defaultEvalueCutoff = 1e-10;
	private static int defaultHitsPerQuery = 1;
	private String[] sortedIds;
	
	//useSW determines if smith-waterman alignment scores should be calculated
	//for all hit pairs before filtering for the top hits
	private boolean useSW = false;


	public static void main(String[] args) {
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
				sequenceFiles[i] = new FastaSequenceFile(infileNames[i]);
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
		//		//verify host
		//		String[] results = remoteHost.executeCommand("ls").getStdout();
		//		for(int i = 0; i < results.length; i++) {
		//			System.out.println(results[i]);
		//		}

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

		String minScoreParam =
			commandLineProperties.getValues(HandyConstants.MIN_SCORE, "30")[0];
		int minScore = Integer.parseInt(minScoreParam);

		String minIdentityParam =
			commandLineProperties.getValues(HandyConstants.MIN_IDENTITY, "25")[0];
		int minIdentity = Integer.parseInt(minIdentityParam);

		String resultFileName = 
			commandLineProperties.getValues(HandyConstants.OUT, null)[0];

		BlatRunner br = new BlatRunner();
		br.setRemoteHost(remoteHost);
		br.setSequenceSets(sequenceFiles);
		br.setQuerySequenceFiles(queryFiles);
		br.setThreadCount(blastThreads);
		br.setHitsPerQuery(blastHitsPerQuery);
		br.setEvalueCutoff(evalue);
		br.setMinScore(minScore);
		br.setMinIdentity(minIdentity);
		br.setResultFileName(resultFileName);
		br.run();
		TextFile result = br.getResults();
		System.out.println("Results: " + result.getLineCount());

	}

	private void setResultFileName(String resultFileName) {
		this.resultFileName = resultFileName;
	}

	public void setMinIdentity(int minIdentity) {
		this.minIdentity = minIdentity;
	}

	public void setMinScore(int minScore) {
		this.minScore = minScore;
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
	public BlatRunner() {

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

			//populate sortedIds array, which is used to provide an index number
			//for each sequence id
			sortedIds = 
				new String[concatenatedSequenceSet.getIDToIndexMap().size()];
			concatenatedSequenceSet.getIDToIndexMap().keySet().toArray(sortedIds);
			Arrays.sort(sortedIds);
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
	 * (number of query files) * (number of input files) >= (number of threads)
	 */
	private void createConcatenatedQuerySets() {
		int numberOfQueryFiles = (int) Math.ceil((double)threadCount/querySequenceFiles.length);
		System.out.println("BlatRunner.createConcatenatedQuerySets() number " +
				"of query files to create: " + numberOfQueryFiles);

		//if only one query set is needed, use the concatenatedSequenceSet
		if(numberOfQueryFiles == 1) {
			querySets = new FastaSequenceFile[]{concatenatedSequenceSet};
		} else {
			querySets = new FastaSequenceFile[numberOfQueryFiles];
			//determine start and stop sequences for each query file
			int sequencesPerFile = 
				concatenatedSequenceSet.getSequenceCount()/ numberOfQueryFiles;
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
		System.out.println("query sets: ");
		for(int i = 0; i < querySets.length; i++) {
			System.out.println(querySets[i].getFileName());
		}

	}

	public void run() {
		//create the concatenated sequence set file.
		createConcatenatedSet();
		createConcatenatedQuerySets();
		//for now, this assumes a single RemoteHost is being used for all blats

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


		//create and start the appropriate number of worker Threads
		Thread[] threads = new Thread[threadCount];
		System.out.println("create blat runner threads: " + threadCount);
		for(int i = 0; i < threads.length; i++) {
			int querySetIndex = i%querySets.length;
			FastaSequenceFile queryForThread = querySets[querySetIndex];
			threads[i] = new Thread(new BlatRunnable(host,
					queryForThread.getFileName(), querySetIndex));
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

	}

	public void setRemoteWorkingDirectory(String dir) {
		remoteWorkingDir = dir;
	}
	
	public void setRnName(String name) {
		runName = name;
	}

	private File getResultFile() throws IOException {
		File r = null;
		if(resultFileName == null) {
			if(runName == null) {
			    r = File.createTempFile("all_blat_", ".out", new File(System.getProperty("user.dir")));
			} else {
				r = new File(runName + "_all_blat.out");
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
					if(setResults[i] != null) {
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

	private void setBlatOutput(TextFile resultFile, int querySetIndex, 
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

	private class BlatRunnable implements Runnable {

		private RemoteHost host;
		private String blatPath;
		private String concatenatedSequenceName;
		private String outputFormat = " -out=blast8 ";
		private int querySetIndex;
		//host index may change after multiple RemoteHost support is added
		private int hostIndex = 0;


		public BlatRunnable(RemoteHost host, String concatenatedSequenceName, int querySetIndex) {
			this.host = host;
			this.concatenatedSequenceName = concatenatedSequenceName;
			this.querySetIndex = querySetIndex;
			System.out.println("new BlatRunnable " + querySetIndex + ", "
					+ concatenatedSequenceName);
		}

		public void run() {

			blatPath = host.getCommandPath("blat");

			try {
				int index = getNextTargetIndex(querySetIndex);
				while (index < sequenceSets.length) {
					System.out.println("blat query set " + querySetIndex 
							+ " against target set " + index + 
							" of " + sequenceSets.length);
					File forOutfileName = File.createTempFile("blat", ".out");
					String outputName = forOutfileName.getName();

					String blatCommand = blatPath + 
					" " + sequenceSets[index].getFullName() + 
					" " + concatenatedSequenceName +
					" -prot " + outputFormat + outputName +
					" -minScore=" + minScore + 
					" -stepSize=" + stepSize + 
					" -minIdentity=" + minIdentity;

					System.out.println("blatCommand: " + blatCommand);
					host.executeCommand(blatCommand);
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
						outputName = filterForTopHits(resultFile, hitsPerQuery);
						resultFile = new TextFile(outputName);
						setBlatOutput(resultFile, querySetIndex, index);
					} catch(IOException ioe) {
						System.out.println("failed to get " + outputName);
					}
					new File(forOutfileName.getName()).delete();
					forOutfileName.delete();
					host.copyFromRemoteWorkingDirectory(outputName, outputName);

>>>>>>> refs/remotes/origin/master
					index = getNextTargetIndex(querySetIndex);
				}
			} catch (IOException ioe) {
				ioe.printStackTrace();
			}
		}

		private TextFile getSWScoresForHits(String pairFileName, String seqFileName) {
			TextFile r = null;
//			System.out.println(">getSWScoresForHits() " + 
//					pairFileName + " " + seqFileName + " " + 
//					Thread.currentThread().getName());
			String outFileName = pairFileName + ".sw";
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
//			System.out.println("<getSWScoresForHits() " + Thread.currentThread().getName());
			return r;
		}

		/**
		 * Filters the result file for the correct number of top hits
		 * per query per target genome. A new output file is written, and 
		 * the original output file is deleted.
		 * @param resultFile
		 * @throws IOException 
		 */
		private String filterForTopHits(TextFile resultFile, int topHits) throws IOException {
			File forFilteredFileName = 
				File.createTempFile("filtered_blat", ".out");
			forFilteredFileName.deleteOnExit();
			String r = forFilteredFileName.getName();
			FileWriter fw = new FileWriter(r);
//			System.out.println("BlatRunner.filterForTopHits() allowing "
//					+ topHits + " hits");
//			System.out.println("filtered file: " + 
//					forFilteredFileName.getAbsolutePath());
			int lineCount = resultFile.getLineCount();

			/*
			 * queryTargetHitCount tracks the hits from a query to
			 * the same target sequence.
			 */
			HashMap queryTargetHitCount = new HashMap();

			/*
			 * queryToHitCounts tracks the total number of hits for 
			 * each query.
			 */
			HashMap queryToHitCounts = new HashMap();

			Integer zero = new Integer(0);

			String tab = "\t";
			Pattern tabPattern = Pattern.compile(tab);
			int queryField = 0;
			int targetField = 1;
			resultFile.openFile();
			String previousQuery = "";
			for(int i = 0; i < lineCount; i++) {
				String line = resultFile.getLine(i);
				String[] fields = tabPattern.split(line);
				if(fields.length > 1) {
					if(!previousQuery.equals(fields[queryField])) {
						queryToHitCounts.remove(previousQuery);
						queryTargetHitCount.clear();
					}

					String key = fields[queryField] + "," + fields[targetField];
					Integer count = (Integer) queryTargetHitCount.get(key);
					if(count == null) {
						count = zero;
					}
					count = new Integer(count.intValue()+1);
					queryTargetHitCount.put(key, count);
					//if this is the first occurrence of this query-target pair,
					//allow it though into the filtered file.
					if(count.intValue() == 1) {
						//now check for the number of hits for this query
						Integer hitCount = 
							(Integer) queryToHitCounts.get(fields[queryField]);
						if(hitCount == null) {
							hitCount = zero;
						}
						hitCount = new Integer(hitCount.intValue()+1);
						queryToHitCounts.put(fields[queryField], hitCount);
						if(hitCount.intValue() <= topHits) {
							fw.write(line);
							fw.write("\n");
						}
					}

					previousQuery = fields[queryField];
				}
			}
			fw.flush();
			fw.close();
			resultFile.closeFile();
			return r;
		}
	}
}
