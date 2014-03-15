package edu.vt.vbi.ci.pepr.alignment;

import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.regex.Pattern;

import edu.vt.vbi.ci.util.CommandLineProperties;
import edu.vt.vbi.ci.util.HandyConstants;
import edu.vt.vbi.ci.util.StringPair;
import edu.vt.vbi.ci.util.file.FastaSequenceFile;
import edu.vt.vbi.ci.util.file.TextFile;

public class AAPairAligner implements Runnable{

	private static final int PAIR_1_COL = 0;
	private static final int PAIR_2_COL = 1;
	private int scoreCol = 2;
	private CommandLineProperties clp;
	private HashMap seqIDToSeqInts;
	private HashMap pairToScore;
	private int nextIndex = 0;
	private StringPair[] pairs;
	private int gapOpen = 11;
	private int gapExtend = 1;
	
	public static void main(String[] args) {
		AAPairAligner aapa = new AAPairAligner(args);
		aapa.run();
	}
	
	/**
	 * commands available: 
	 * HandyConstants.SEQ_FILE_NAME seq_file
	 * HandyConstants.FILE pair file
	 * HandyConstants.MAX_CONCURRENT_PROCESSES max_concurrent_processes
	 * HandyConstants.OUT out
	 * @param args
	 */
	public AAPairAligner(String[] args) {
		clp = new CommandLineProperties(args);
	}
	
	private void loadPairs(String pairFileName) throws IOException {
		String tab = "\\s+";
//		System.out.println("1create TextFile...");
		TextFile pairFile = new TextFile(pairFileName);
		pairFile.openFile();
		Integer zero = new Integer(0);
		int lineCount = pairFile.getLineCount();
//		lineCount = 8000000;
		int hashStartSize = lineCount + lineCount/4;
		pairToScore = new HashMap(hashStartSize);

		Pattern tabRE = Pattern.compile(tab);
//		System.out.println("loading lines in file: " + lineCount);
		for(int i = 0; i < lineCount; i++) {
//			if((i+1) %10000== 0) {
//				System.out.print(".");
//				if((i+1) % 100000 == 0) {
//					System.out.println(" " + (i+1));
//				}
//			}
			String[] fields = tabRE.split(pairFile.getLine(i));
			if(fields.length > PAIR_2_COL) {
				StringPair pair = 
//					new StringPair(fields[PAIR_1_COL], 
//							fields[PAIR_2_COL]);
				new StringPair(new String(fields[PAIR_1_COL]), 
						new String(fields[PAIR_2_COL]));
//				System.out.println("adding pair: " + pair);
				pairToScore.put(pair, zero);
			} else {
//				System.out.println("problem with line " + i + 
//						". fields.length: " + fields.length);
			}
		}
		pairFile.closeFile();
//		System.out.println();
//		System.out.println("loaded " + pairToScore.size() + " pairs");
//		System.out.println("create pairs array...");
		pairs = new StringPair[pairToScore.size()];
//		System.out.println("fill pairs array...");
		pairToScore.keySet().toArray(pairs);
		
	}

	private void loadSequenceFiles(String[] sequenceFiles) throws IOException {

		for(int i = 0; i < sequenceFiles.length; i++) {
			FastaSequenceFile fsf = new FastaSequenceFile(sequenceFiles[i]);
			loadSequences(fsf);
		}
	}

	private void loadSequences(FastaSequenceFile seq) {
		if(seqIDToSeqInts == null) {
			seqIDToSeqInts = new HashMap();
		}
		
		HashMap idToTokenMap = seq.getIDToIndexMap();
		String[] ids = new String[idToTokenMap.size()];
		idToTokenMap.keySet().toArray(ids);
		
		for(int i = 0; i < ids.length; i++) {
			String sequence = 
				seq.getPlainSequence(seq.getIndexOfSequence(ids[i]));
			int[] sequenceInts = 
				AlignmentUtilities.convertAminaAcidSequenceToInts(sequence);
			seqIDToSeqInts.put(ids[i], sequenceInts);
		}
	}	
	
	public void run() {
//		System.out.println("AAPairAligner.run()");
		String[] sequenceFiles = clp.getValues(HandyConstants.SEQ_FILE_NAME);
		String pairFileName = clp.getValues(HandyConstants.FILE, null)[0];
		try {
//			System.out.println("load sequence files...");
			loadSequenceFiles(sequenceFiles);
//			System.out.println("load pair file...");
			loadPairs(pairFileName);
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		scoreCol = 
			Integer.parseInt(clp.getValues(HandyConstants.SCORE_COL, "2")[0]);
			
		boolean banded = 
			clp.getValues(HandyConstants.BANDED, HandyConstants.FALSE)[0].
				equalsIgnoreCase(HandyConstants.TRUE);
		int bandWidth = 
			Integer.parseInt(clp.getValues(HandyConstants.BAND_WIDTH, "32")[0]);
			
		int threadCount = 
			Integer.parseInt(clp.getValues(
					HandyConstants.MAX_CONCURRENT_PROCESS_PARAM, "1")[0]);
		Thread[] threads = new Thread[threadCount];
//		System.out.println("banded: " + banded);
		if(banded) {
//			System.out.println("band width: " + bandWidth);
		}
//		System.out.println("create and start " + threadCount + " alignment threads");
		for(int i = 0; i < threads.length; i++) {
			threads[i] = new Thread(new AlignmentRunnable(banded, bandWidth));
			threads[i].start();
		}
//		System.out.println("waiting for threads to finish all " + pairs.length 
//				+ " alignments");
		
		//wait for threads to finish
		for(int i = 0; i < threads.length; i++) {
			try {
				threads[i].join();
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		
		boolean m8  = !clp.getValues(HandyConstants.M8, 
		HandyConstants.FALSE)[0].equalsIgnoreCase(HandyConstants.FALSE);
		
		if(m8) {
			try {
				replaceScores(pairFileName, pairToScore, clp.getValues(HandyConstants.M8, 
						"AAPairAligner.m8")[0]);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		} else {
			//print pairs and scores to file
			String outFileName = 
				clp.getValues(HandyConstants.OUT, "AAPairAligner.out")[0];
			try {
				printOutputFile(outFileName);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

		}

	}
	
	private void replaceScores(String pairFileName, HashMap pairToScore2,
			String outFileName) throws IOException {
		String tab = "\t";
//		System.out.println("create TextFile...");
		TextFile pairFile = new TextFile(pairFileName);
		pairFile.openFile();

//		System.out.println("AAPairAligner.replaceScores() writing to " + outFileName);
		FileWriter fw = new  FileWriter(outFileName);
		
		int lineCount = pairFile.getLineCount();
		Pattern tabRE = Pattern.compile("\t");
//		System.out.println("loading lines in file: " + lineCount);
		for(int i = 0; i < lineCount; i++) {
//			if((i+1) %10000== 0) {
//				System.out.print(".");
//				if((i+1) % 100000 == 0) {
//					System.out.println(" " + (i+1));
//				}
//			}
			String[] fields = tabRE.split(pairFile.getLine(i));
			if(fields.length > PAIR_2_COL) {
				StringPair pair = 
				new StringPair(new String(fields[PAIR_1_COL]), 
						new String(fields[PAIR_2_COL]));
				Integer score = (Integer)pairToScore.get(pair);
				if(score != null) {
					fields[scoreCol] = score.toString();
				} else {
					fields[scoreCol] = "No Score";
				}
				fw.write(fields[0]);
				for(int j = 1; j < fields.length; j++) {
					fw.write(tab);
					fw.write(fields[j]);
				}
				fw.write("\n");
			}
			
		}
		fw.flush();
		fw.close();
		pairFile.closeFile();
	}

	private void printOutputFile(String outFileName) throws IOException {
		FileWriter fw = new FileWriter(outFileName);
		String tab = "\t";
		String nl = "\n";
		for(int i = 0; i < pairs.length; i++) {
			StringPair pair = pairs[i];
			int score = ((Integer)pairToScore.get(pair)).intValue();
			String line1 = pair.getA() + tab + pair.getB() + tab + score + nl; 
			String line2 = pair.getB() + tab + pair.getA() + tab + score + nl; 
			fw.write(line1);
			fw.write(line2);
		}
		fw.flush();
		fw.close();
	}

	private class AlignmentRunnable implements Runnable 
	{
		private SWProteinAligner aligner = 
			new SWProteinAligner(gapOpen, gapExtend, SWProteinAligner.BLOSUM_62);
		private boolean banded = false;
		private int bandWidth = 32;

		public AlignmentRunnable(boolean banded, int bandWidth) {
			this.banded = banded;
			this.bandWidth = bandWidth;
		}
		
		public void run() {
			
			for(int i = nextIndex++; i < pairs.length; i = nextIndex++) {
//				if((i+1) %1000 == 0) {
//					System.out.print(".");
//					if((i+1) % 25000 == 0){
//						System.out.println(" " + (i+1));
//					}
//				}
				//get pair
				StringPair pair = pairs[i];
				//get sequence ints for the pair
				int[] seq1 = (int[]) seqIDToSeqInts.get(pair.getA());
				int[] seq2 = (int[]) seqIDToSeqInts.get(pair.getB());
				
				//get score for alignment
				Integer score = null;
				if(banded) {
					score = new Integer(aligner.getBandedAlignmentScore(seq1,
							seq2, bandWidth));
				} else {
					score = new Integer(aligner.getAlignmentScore(seq1, seq2));
				}
				
				//put alignment score in the pairToScore map
				pairToScore.put(pair, score);
			}
		}
		
	}
	
}
