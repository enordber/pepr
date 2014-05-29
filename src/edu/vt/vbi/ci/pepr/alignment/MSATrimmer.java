package edu.vt.vbi.ci.pepr.alignment;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

import edu.vt.vbi.ci.util.ExecUtilities;
import edu.vt.vbi.ci.util.HandyConstants;

/**
 * Trims a SequenceAlignment. The actual trimming method 
 * has not bee determined yet. For initial implementation, 
 * a dummy trim() method will be implemented that simply
 * returns the SequenceAlignment unaltered.
 * 
 * @author enordber
 *
 */
public class MSATrimmer {

	private boolean debug = false;

	private double portionOfSequencesForConservedSite = 0.8;

	private File workingDir;
	private boolean keepFiles = false;

	public MSATrimmer() {
		//check for debug System property
		String debugProp = System.getProperty(HandyConstants.DEBUG);
		if(debugProp != null && debugProp.equals(HandyConstants.TRUE)) {
			debug = true;
		}
	}
	/**
	 * @param alignment
	 * @return
	 */
	public SequenceAlignment trim(SequenceAlignment alignment) {
		SequenceAlignment r = null;

		try {
			r = trimWithGBlocks(alignment);
			r.setName(alignment.getName() + "_trim");
		} catch (IOException e) {
			// TODO Auto-generated catch block
//			e.printStackTrace();
		}
		return r;
	}

	/**
	 * Uses GBlocks to trim the alignment.
	 * 
	 * @param alignment
	 * @return
	 * @throws IOException 
	 */
	private SequenceAlignment trimWithGBlocks(SequenceAlignment alignment) 
	throws IOException {
		SequenceAlignment r = null;
		//write the alignment to a file, in a format that GBlocks
		//accepts as input. Replace sequence titles with index numbers for
		//the Gblocks input file. These will be replaced with the full 
		//title again after the alignment is trimmed.
		int sequenceCount = alignment.getNTax();
		if(debug) {
			System.out.println("MSATrimmer.trimWithGblocks() sequences in alignment: "
					+ sequenceCount);
		}
		String alignmentFileName = getWorkingDir().getAbsolutePath()
		+ System.getProperty("file.separator") + alignment.getName();
		File gblocksInputFile = new File(alignmentFileName);
		FileWriter fw = new FileWriter(gblocksInputFile);
		String[] taxa = alignment.getTaxa();
				
		for(int i = 0; i < sequenceCount; i++) {
			String seq = alignment.getSequenceString(i);
			//Gblocks limits sequence titles to 73 characters (74, including the '>')
			String taxon = taxa[i];
			if(taxon.length() > 73) {
				taxon = taxon.substring(0, 73);
			}
			fw.write(">" + taxon + "\n");
			//			fw.write(">" + i + "\n");
			fw.write(seq);
			fw.write("\n");
		}
		fw.flush();
		fw.close();

		int b1ParamValue = 
			(int)Math.ceil(taxa.length * portionOfSequencesForConservedSite);
		String gblocksB1Param = "-b1=" + b1ParamValue;
		String gblocksB3Param = "-b3=8";
		String gblocksB5Param = "-b5=h";
		//run GBlocks on the temporary alignment file
		String gblocksCommand = ExecUtilities.getCommandPath("Gblocks") + 
		" " + gblocksInputFile.getAbsolutePath() 
		+ " -o " + gblocksB1Param + " " + gblocksB3Param + " " + gblocksB5Param;
		if(debug) {
			System.out.println("MSATrimmer.trimWithGblocks() " + gblocksCommand);
		}
		ExecUtilities.exec(gblocksCommand);

		//read results from GBlocks output file
		String gblocksResultFileName = 
			gblocksInputFile.getAbsolutePath() + "-gb";
		r = SequenceAlignmentParser.
		parseFastaAlignmentFile(gblocksResultFileName);
		r.setTitles(alignment.getSequenceTitles());
		//		r.setTitles(alignment.getTaxa());
		r.setTaxa(alignment.getTaxa());

		if(!keepFiles) {
			//delete GBlocks input and output files on exit
			gblocksInputFile.deleteOnExit();
			File gblocksOutput = new File(gblocksResultFileName);
			gblocksOutput.deleteOnExit();
			gblocksOutput = new File(gblocksResultFileName + ".htm");
			gblocksOutput.deleteOnExit();
		}
		return r;
	}

	/**
	 * Create a new SequenceAlignment that contains only the 
	 * potentially phylogenetically informative columns of the
	 * original alignment. For now, this means any columns that
	 * are identical in all sequences are removed. Columns that
	 * contain only one character other than gaps are also
	 * considered uninformative and are removed. This may be 
	 * a bad assumption, as there may be methods that use gaps
	 * as informative.
	 * 
	 * @param alignment
	 * @return
	 */
	public SequenceAlignment trimUniformativeColumns(SequenceAlignment alignment) {
		SequenceAlignment r = null;
		int length = alignment.getLength();

		int uniformCount = 0;
		int variableCount = 0;

		//store variable columns, to be used below to create the trimmed
		//alignment. There may be a better default length for this 
		//ArrayList: maybe assume some percentage of variable columns?
		ArrayList variableColumns = new ArrayList(length);

		for(int i = 0; i < length; i++) {
			char[] column = alignment.getColumn(i);

			//find the first non-gap character
			int j = 0;
			while(j < column.length && column[j] == SequenceAlignment.GAP_CHAR) {
				j++;
			}

			boolean uniform = true;

			if(j < column.length) {
				char firstChar = column[j];
				j++;
				for(; j < column.length && uniform; j++) {
					uniform = column[j] == firstChar ||
					column[j] == SequenceAlignment.GAP_CHAR;
				}
			}


			if(uniform) {
				uniformCount++;
			} else {
				variableCount++;
				variableColumns.add(column);
			}
		}

		//create char[][] for trimmed alignment. This requires swapping
		//the dimensions of the stored variable column arrays
		char[][] trimmedMatrix = 
			new char[alignment.getNTax()][variableColumns.size()];

		char[][] variableColumnArrays = new char[variableColumns.size()][];
		variableColumns.toArray(variableColumnArrays);

		for(int i = 0; i < trimmedMatrix.length; i++) {
			for(int j = 0; j < trimmedMatrix[i].length; j++) {
				trimmedMatrix[i][j] = variableColumnArrays[j][i];
			}
		}
		r = new SequenceAlignment(trimmedMatrix, 
				alignment.getSequenceTitles(), alignment.getTaxa());
		r.setName(alignment.getName()+"_ui");
		if(debug) {
			System.out.println("MSATrimmer uniform: " + uniformCount + 
					" variable: " + variableColumns.size());
		}
		return r;
	}

	public SequenceAlignment trimUniformColumns(SequenceAlignment alignment) {
		SequenceAlignment r = null;
		int length = alignment.getLength();

		int uniformCount = 0;
		int variableCount = 0;

		//store variable columns, to be used below to create the trimmed
		//alignment. There may be a better default length for this 
		//ArrayList: maybe assume some percentage of variable columns?
		ArrayList variableColumns = new ArrayList(length);

		int[] steps = alignment.getMinimumStepsPerSite();

		for(int i = 0; i < length; i++) {

			boolean uniform = steps[i] == 0;

			if(uniform) {
				uniformCount++;
			} else {
				variableCount++;
				char[] column = alignment.getColumn(i);				
				variableColumns.add(column);
			}
		}

		//create char[][] for trimmed alignment. This requires swapping
		//the dimensions of the stored variable column arrays
		char[][] trimmedMatrix = 
			new char[alignment.getNTax()][variableColumns.size()];

		char[][] variableColumnArrays = new char[variableColumns.size()][];
		variableColumns.toArray(variableColumnArrays);

		for(int i = 0; i < trimmedMatrix.length; i++) {
			for(int j = 0; j < trimmedMatrix[i].length; j++) {
				trimmedMatrix[i][j] = variableColumnArrays[j][i];
			}
		}
		r = new SequenceAlignment(trimmedMatrix, 
				alignment.getSequenceTitles(), alignment.getTaxa());
		r.setName(alignment.getName()+"_ui");
		if(debug) {
			System.out.println("MSATrimmer uniform: " + uniformCount + 
					" variable: " + variableColumns.size());
		}
		return r;
	}

	/**
	 * Topologically uninformative columns include columns with
	 * only one character, columns with only one character other
	 * than gaps, and columns with two characters other than gaps
	 * where one of the characters occurs only once.
	 *  
	 * @param alignment
	 * @return
	 */
	public SequenceAlignment trimTopologicallyUninformativeColumns(
			SequenceAlignment alignment) {
		SequenceAlignment r = null;
		int length = alignment.getLength();

		int uninformativeCount = 0;
		int informativeCount = 0;

		//store variable columns, to be used below to create the trimmed
		//alignment. There may be a better default length for this 
		//ArrayList: maybe assume some percentage of variable columns?
		ArrayList informatieColumns = new ArrayList(length);

		for(int i = 0; i < length; i++) {
			char[] column = alignment.getColumn(i);

			//find the first non-gap character
			int j = 0;
			while(j < column.length && column[j] == SequenceAlignment.GAP_CHAR) {
				j++;
			}

			boolean uniform = true;

			if(j < column.length) {
				char firstChar = column[j];
				j++;
				for(; j < column.length && uniform; j++) {
					uniform = column[j] == firstChar ||
					column[j] == SequenceAlignment.GAP_CHAR;
				}

				if(uniform) {
					uninformativeCount++;
				} else {
					//the column is non-uniform. To be kept, there must be 
					//at least two characters that each occur at least two times
					HashSet multipleClassCharacters = new HashSet();
					int[] counts = new int[column.length];
					for(j = 0; j < column.length; j++) {
						if(column[j] != SequenceAlignment.GAP_CHAR) {
							for(int k = 0; k < column.length; k++) {
								if(column[j] == column[k]) {
									counts[j]++;
									if(counts[j] > 1) {
										multipleClassCharacters
										.add(new Character(column[j]));
										k = column.length;
										if(multipleClassCharacters.size() ==2) {
											informativeCount++;
											informatieColumns.add(column);
											j = column.length;
										}
									} else if(k == (column.length-1)) {
										//this is the last item in the column, 
										//and it has not been found to be 
										//informative, so count it as uninformative
										uninformativeCount++;
									}
								}
							}
						}
					}
				}
			}
		}

		//create char[][] for trimmed alignment. This requires swapping
		//the dimensions of the stored variable column arrays
		char[][] trimmedMatrix = 
			new char[alignment.getNTax()][informatieColumns.size()];

		char[][] variableColumnArrays = new char[informatieColumns.size()][];
		informatieColumns.toArray(variableColumnArrays);

		for(int i = 0; i < trimmedMatrix.length; i++) {
			for(int j = 0; j < trimmedMatrix[i].length; j++) {
				trimmedMatrix[i][j] = variableColumnArrays[j][i];
			}
		}
		r = new SequenceAlignment(trimmedMatrix, 
				alignment.getSequenceTitles(), alignment.getTaxa());
		if(debug) {
			System.out.println("MSATrimmer uninformative: " + uninformativeCount + 
					" informative: " + informatieColumns.size());
		}
		return r;
	}

	private File getWorkingDir() {
		if(workingDir == null) {
			workingDir = new File(System.getProperty("user.dir"));
		}
		return workingDir;
	}

	public void setWorkingDir(File dir) {
		workingDir = dir;
	}

}
