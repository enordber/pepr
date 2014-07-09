/*
 * Created on Nov 17, 2006
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.vt.vbi.ci.util;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;

import edu.vt.vbi.ci.util.file.FastaSequenceFile;
import edu.vt.vbi.ci.util.file.FastaSequenceSet;
import edu.vt.vbi.ci.util.file.TextFile;

/**
 * This class reads a file describing sets of sequences
 * and extracts those sets from fasta sequence files.
 * A single file is created for each set of sequences.
 * 
 * @author enordber
 */
public class SequenceSetExtractor {

	/*
	 * Command line options.
	 */
	private static final String SET_FILE_CMD = "set";
	private static final String OG_FILE_CMD = "og";
	private static final String SEQUENCE_FILE_CMD = "seq";
	private static final String OUTPUT_BASE_CMD = "out";
	private static final String OUTPUT_EXT_CMD = "ext";

	/**
	 * A list of sequence files containing the sequences
	 * to be extracted.
	 */
	private FastaSequenceSet[] sequenceFiles;

	/**
	 * A list of Strings describing the sets of sequences
	 * to be extracted. Each element is a tab-delimited
	 * String, where each field is an identifier of
	 * a sequence in the set.
	 */
	private String[] groupDefinitions;

	/**
	 * The base of the names for output files. They
	 * will have a number appended, indicating the
	 * number of the sequence set, and will end in
	 * either .faa or .ffn.
	 */
	private String outputBaseName;

	/**
	 * The extension to be used for output sequence files.
	 * Will be either "faa" or "ffn". 
	 */
	private String outputFileExtension;

	private HashMap parameters;

	/*
	 * sequenceIndexMap enables quick searching for sequences based
	 * on id tokens.
	 * key: id token
	 * value: int[] 
	 *        [0]=index of sequence file in sequenceFiles
	 *        [1]=index of sequence in the file
	 */
	private HashMap sequenceIndexMap;

	public static void main(String[] args) {
		if(args.length == 0 || 
				args[0].startsWith("-h") || 
				args[0].startsWith("-?")) {
			printHelp();
			System.exit(0);

		}
		new SequenceSetExtractor(args);
	}

	/**
	 * 
	 */
	private static void printHelp() {
		String[] help = new String[]{
				"SequenceSetExtractor extracts sets of sequences from Fasta sequence files",
				"Command-line parameters: ",
				"\t" + OG_FILE_CMD + "\tname of OG file defining sequence sets",
				"\t" + SET_FILE_CMD + "\tname of file defining sequence sets",
				"\t" + SEQUENCE_FILE_CMD + "\tname of fasta sequence file(s)",
				"\t" + OUTPUT_BASE_CMD + "\tbase of output file name(s)",
				"\t" + OUTPUT_EXT_CMD + "\textension to use for output file(s)"
		};

		for(int i = 0; i < help.length; i++) {
			System.out.println(help[i]);
		}

	}

	public SequenceSetExtractor(String[] args) {
		CommandLineProperties clp = new CommandLineProperties(args);

		outputBaseName = clp.getValues(OUTPUT_BASE_CMD)[0];
		outputFileExtension = clp.getValues(OUTPUT_EXT_CMD)[0]; 
		String setFileName = clp.getValues(SET_FILE_CMD, "")[0];
		String ogFileName = clp.getValues(OG_FILE_CMD, "")[0];
		boolean useOGFile = false;
		if(setFileName.equals("")) {
			setFileName = ogFileName;
			useOGFile = true;
		}
		String[] seqFileNames = clp.getValues(SEQUENCE_FILE_CMD);

		groupDefinitions = getFileLines(setFileName);

		sequenceFiles = getFastaSequenceFiles(seqFileNames);

		int minSetSize = 
			Integer.parseInt(clp.getValues(HandyConstants.MIN, "0")[0]);
		int maxSetSize = 
			Integer.parseInt(clp.getValues(HandyConstants.MAX, 
			"99999999")[0]);

		if(useOGFile) {
			processOGSequenceSets(minSetSize, maxSetSize);
		} else {
			processSequenceSets(minSetSize, maxSetSize);
		}

	}

	public SequenceSetExtractor(String setFileName, 
			FastaSequenceSet[] sequenceFiles, String outputBaseName, 
			String outputFileExtension) {

		this.outputBaseName = outputBaseName;
		this.outputFileExtension = outputFileExtension;
		groupDefinitions = getFileLines(setFileName);
		this.sequenceFiles = sequenceFiles;

		int minSetSize = 0;
		int maxSetSize = Integer.MAX_VALUE;
		processSequenceSets(minSetSize, maxSetSize);
	}

	/**
	 * This is the method that does the actual work.
	 */
	private void processSequenceSets(int minSetSize, int maxSetSize) {
		try {
			String delimiter = "\t";
			//for each set, 
			for(int set = 0; set < groupDefinitions.length; set++) {
				//get the list of ids for the set by parsing the set line
				String[] setIds = groupDefinitions[set].split(delimiter);
				if(setIds.length >= minSetSize && setIds.length <= maxSetSize) {
					//create a file for output
					String fileName = outputBaseName + 
					"_" + set + "." + 
					outputFileExtension;
					FileWriter fw = new FileWriter(fileName);
					//search for each id in the set in the fasta files
					for(int id = 0; id < setIds.length; id++) {
						//check each FastaFile until a match is found, or
						//all files have been checked
						String[] seq = null;
						seq = getSequenceByQuickLookup(setIds[id]);
						for(int file = 0; 
						file < sequenceFiles.length && seq == null; 
						file++) {
							seq = sequenceFiles[file].getSequence(setIds[id]);
						}
						//				print the sequence to the output file
						if(seq != null) {
							for(int i = 0; i < seq.length; i++) {
								fw.write(seq[i]);
								fw.write("\n");
							}
						} 
					}

					fw.flush();
					fw.close();
				}
			}
		} catch(IOException ioe) {
			ioe.printStackTrace();
		}
	}

	/**
	 * This is the method that does the actual work.Assumes the first field
	 * will be the name of the OG.
	 */
	private void processOGSequenceSets(int minSetSize, int maxSetSize) {
		try {
			String delimiter = "\t";
			String openParen = "(";
			//for each set, 
			for(int set = 0; set < groupDefinitions.length; set++) {
				delimiter = "\t";
				//get the list of ids for the set by parsing the set line
				String[] setIds = groupDefinitions[set].split(delimiter);
				int endIndex = setIds[0].indexOf(openParen);
				if(endIndex < 0) {
					System.out.println("no openParen found in set " + set 
							+ " set name: " + setIds[0]);
				} else {
					String ogName = 
						setIds[0].substring(0, endIndex);

					delimiter = "\\s+";
					String memberString = setIds[1].trim();
					setIds = memberString.split(delimiter);
					if(setIds.length > 1 && setIds.length + 1 >= minSetSize &&
							setIds.length + 1 <= maxSetSize) {
						//create a file for output
						String fileName = outputBaseName + ogName +
						"_" + set + "." + 
						outputFileExtension;
						FileWriter fw = new FileWriter(fileName);
						//search for each id in the set in the fasta files
						for(int id = 0; id < setIds.length; id++) {
							//check each FastaFile until a match is found, or
							//all files have been checked
							endIndex = setIds[id].indexOf(openParen);
							if(endIndex < 0) {
								System.out.println("no openParen found in set " + set 
										+ " set name: " + ogName + " member " 
										+ id + ": " + setIds[id]);
							}
							String seqId  = 
								setIds[id].substring(0, endIndex).trim();
							String[] seq = null;
							seq = getSequenceByQuickLookup(seqId);
							for(int file = 0; 
							file < sequenceFiles.length && seq == null; 
							file++) {
								seq = sequenceFiles[file].getSequence(seqId);
							}
							//				print the sequence to the output file
							if(seq != null) {
								for(int i = 0; i < seq.length; i++) {
									fw.write(seq[i]);
									fw.write("\n");
								}
							}
						}

						fw.flush();
						fw.close();
					}
				}
			}
		} catch(IOException ioe) {
			ioe.printStackTrace();
		}
	}

	/**
	 * Returns the lines for the sequence with the given id, if it
	 * can be found by the quick lookup index in sequenceIndexMap.
	 * Returns null if the id is not found.
	 * 
	 * @param string
	 * @return
	 */
	private String[] getSequenceByQuickLookup(String seqId) {
		String[] r = null;
		if(sequenceIndexMap == null) {
			buildSequenceIndexMap();
		}
		int[] indices = (int[]) sequenceIndexMap.get(seqId);
		if(indices != null && indices.length == 2) {
			r = sequenceFiles[indices[0]].getSequence(indices[1]);
		}
		return r;
	}

	public static String[] getSequenceLines(String id, FastaSequenceSet[] fsfList) {
		String[] r = null;

		return r;
	}

	/**
	 * @param seqFileNames
	 * @return
	 */
	private FastaSequenceSet[] getFastaSequenceFiles(String[] seqFileNames) {
		FastaSequenceSet[] r = null;
		try {
			ArrayList fastaFileList = new ArrayList(seqFileNames.length);
			for(int i = 0; i < seqFileNames.length; i++) {
				FastaSequenceFile nextFile = new FastaSequenceFile(seqFileNames[i]);
				nextFile.open();
				fastaFileList.add(nextFile);
			}

			r = new FastaSequenceSet[fastaFileList.size()];
			fastaFileList.toArray(r);
		} catch(IOException ioe) {
			ioe.printStackTrace();
		}
		return r;
	}

	/**
	 * @param setFileName
	 * @return
	 */
	private String[] getFileLines(String setFileName) {
		String[] r = null;
		try {
			TextFile tf = new TextFile(setFileName);
			r = tf.getAllLines();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return r;
	}

	/**
	 * @param args
	 */
	private void readParameters(String[] args) {
		if(parameters == null) {
			parameters = new HashMap();
		}
		String activeCommand = null;
		for(int i = 0; i < args.length; i++) {
			if(args[i].equals(OUTPUT_BASE_CMD)) {
				activeCommand = OUTPUT_BASE_CMD;
				i++;
			} else if(args[i].equals(OUTPUT_EXT_CMD)) {
				activeCommand = OUTPUT_EXT_CMD;
				i++;
			} else if(args[i].equals(SEQUENCE_FILE_CMD)) {
				activeCommand = SEQUENCE_FILE_CMD;
				i++;
			} else if(args[i].equals(SET_FILE_CMD)) {
				activeCommand = SET_FILE_CMD;
				i++;
			}

			if(activeCommand == OUTPUT_BASE_CMD) {
				parameters.put(OUTPUT_BASE_CMD, args[i]);
			} else if(activeCommand == OUTPUT_EXT_CMD) {
				parameters.put(OUTPUT_EXT_CMD, args[i]);
			} else if(activeCommand == SEQUENCE_FILE_CMD) {
				//multiple sequence files can be used, so
				//store these in a Collection
				if(!parameters.containsKey(SEQUENCE_FILE_CMD)) {
					parameters.put(SEQUENCE_FILE_CMD, new ArrayList());
				}

				Collection c = (Collection)parameters.get(SEQUENCE_FILE_CMD);
				c.add(args[i]);
			} else if(activeCommand == SET_FILE_CMD) {
				parameters.put(SET_FILE_CMD, args[i]);
			}
		}
	}

	/**
	 * Reads all sequence titles from all sequence files and constructs
	 * a quick look-up index map. This map links from a sequence Id token
	 * to a sequence file index and sequence index within that file.
	 */
	private void buildSequenceIndexMap() {
		sequenceIndexMap = new HashMap();
		//go through each file
		for(int i = 0; i < sequenceFiles.length; i++) {
			//get titles for file
			String[] titles = sequenceFiles[i].getTitles();
			for(int j = 0; j < titles.length; j++) {
				//get token for the title
				String token = FastaSequenceFile.getIDToken(titles[j]);
				//add entry to index map
				sequenceIndexMap.put(token, new int[]{i,j});
			}
		}
	}
}
