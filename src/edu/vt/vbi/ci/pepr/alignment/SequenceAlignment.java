package edu.vt.vbi.ci.pepr.alignment;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;

import edu.vt.vbi.ci.pepr.tree.Bipartition;
import edu.vt.vbi.ci.util.ExtendedBitSet;
import edu.vt.vbi.ci.util.IntList;
import edu.vt.vbi.ci.util.file.FastaUtilities;

/**
 * Represents an alignment of two or more sequences.
 * For now I intend this to cover nucleotide and amino
 * acid sequence alignment. If that becomes a problem,
 * I can split into separate classes for each type.
 * 
 * This representation is independent of how the alignment
 * was generated. This representation is independent of 
 * any file formats. Parsing from files and writing to 
 * files in various formats will be handled by other classes.
 * This representation is independent of graphical display
 * of the alignment. Formatting and display will be handled
 * by other classes. The purpose of this class is to serve
 * as a common currency for sequence alignment data.
 * @author enordber
 *
 */
public class SequenceAlignment {

	private boolean debug = false;
	private static Random random = new Random();

	/*
	 * stepsPerSite stores the number of steps required per site to
	 * fit to the/a most parsimonious tree, as determined by protpars.
	 */
	private int[] stepsPerSite;

	/*
	 * minimum steps per site is based on the number of different 
	 * states in each column of the alignment. This is the number of 
	 * bipartitions for the column -1.
	 */
	private int[] minimumStepsPerSite;

	public static final String GAP_STRING = "-";
	public static final char GAP_CHAR = '-';
	public static final char MISSING_CHAR = '?';

	private String name;

	/*
	 * Stores sequences in alignment. Each element
	 * is a char[] for one sequences. Each char[]
	 * should be the same length, with gap characters
	 * providing any padding needed.
	 */
	private char[][] alignedSequenceChars = new char[0][0];

	/*
	 * The titles, for each sequence in the 
	 * alignment.The index values are the same as
	 * in alignedSequenceChars
	 */
	private String[] sequenceTitles;

	/*
	 * The sequence name is the first token (until the first space)
	 * in the title (without the '>')
	 */
	private String[] sequenceNames;

	/*
	 * The taxon names or ids for each sequence in
	 * the alignment. The index values are the same
	 * as in alignedSequenceChars. These values can 
	 * be used to relate alignments of different
	 * sequences from a common set of taxa. (e.g.
	 * to support generating a concatenated alignment.)
	 */
	private String[] sequenceTaxa;

	private Bipartition[][] columnBipartitions;
	private int[][] columnCharacterClasses;

	public SequenceAlignment() {
		sequenceTitles = new String[0];
	}

	public SequenceAlignment(String[] alignedSequences, String[] titles, String[] taxa) {
		this(alignedSequences, titles);
		setTaxa(taxa);
	}

	public SequenceAlignment(char[][] alignedSequences, String[] titles, String[] taxa) {
		this(alignedSequences, titles);
		setTaxa(taxa);
	}

	public SequenceAlignment(String[] alignedSequences, String[] titles) {
		alignedSequenceChars = new char[alignedSequences.length][];
		for(int i = 0; i < alignedSequences.length; i++) {
			alignedSequenceChars[i] = alignedSequences[i].toCharArray();
		}
		setTitles(titles);
	}

	public SequenceAlignment(char[][] alignedSequences, String[] titles) {
		this();
		setAlignedSequences(alignedSequences);
		setTitles(titles);
	}

	protected void setAlignedSequences(char[][] alignedSequences) {
		alignedSequenceChars = new char[alignedSequences.length][];
		for(int i = 0; i < alignedSequenceChars.length; i++) {
			alignedSequenceChars[i] = new char[alignedSequences[i].length];
			System.arraycopy(alignedSequences[i], 0, alignedSequenceChars[i], 0,
					alignedSequences[i].length);
		}
	}

	public void setTitles(String[] titles) {
		sequenceTitles = new String[titles.length];
		System.arraycopy(titles, 0, sequenceTitles, 0, titles.length);

		//remove '>' from titles, if present
		for(int i = 0; i < sequenceTitles.length; i++) {
			if(sequenceTitles[i].startsWith(">")) {
				sequenceTitles[i] = sequenceTitles[i].substring(1);
			}
		}

		//set the sequence names
		sequenceNames = new String[sequenceTitles.length];
		String ws = "\\s+";
		for(int i = 0; i < sequenceNames.length; i++) {
			sequenceNames[i] = sequenceTitles[i].split(ws)[0];
		}
	}

	public void setTaxa(String[] taxa) {
		sequenceTaxa = new String[taxa.length];
		System.arraycopy(taxa, 0, sequenceTaxa, 0, taxa.length);
	}

	/**
	 * Adds a new sequence to the alignment. The sequence must
	 * already be aligned. That is, it must contain gap characters
	 * where appropriate. This class does not perform alignments.
	 * 
	 * @param sequence
	 */
	public void addSequence(String sequence, String title) {
		char[] seqChars = sequence.toCharArray();
		if(alignedSequenceChars == null) {
			alignedSequenceChars = new char[0][];
			if(debug) {
				System.out.println("SequenceAlignment.addSequence() initialize alignedSequenceChars");
			}
		}

		char[][] tempHolder = alignedSequenceChars;
		alignedSequenceChars = new char[alignedSequenceChars.length+1][];
		if(debug) {
		    System.out.println("SequenceAlignment.addSequence() expand alignedSequenceChars for new sequence");
		}
		for(int i = 0; i < tempHolder.length; i++) {
			alignedSequenceChars[i] = tempHolder[i];
		}
		alignedSequenceChars[tempHolder.length] = seqChars;

		padSequencesIfNeeded();

		//add new title, after removing leading '>', if present
		if(title.startsWith(">")) {
			title = title.substring(1);
		}
		String[] newTitles = new String[sequenceTitles.length+1];
		System.arraycopy(sequenceTitles, 0, newTitles, 0, sequenceTitles.length);
		newTitles[sequenceTitles.length] = title;
		sequenceTitles = newTitles;

		//add new sequence name
		if(sequenceNames == null) {
			sequenceNames = new String[0];
		}
		String newName = title.split("\\s+")[0];
		String[] newNames = new String[sequenceNames.length + 1];
		System.arraycopy(sequenceNames, 0, newNames, 0, sequenceNames.length);
		newNames[sequenceNames.length] = newName;
		sequenceNames = newNames;
	}

	/**
	 * Removes the specified sequence from this alignment.
	 * 
	 * @param index
	 */
	public void removeSequence(int index) {
		int secondLength = getNTax() - index -1;
		//remove index from alignedSequenceChars
		char[][] nASC = new char[alignedSequenceChars.length-1][];
		System.arraycopy(alignedSequenceChars, 0, nASC, 0, index);
		if(secondLength > 0) {
			System.arraycopy(alignedSequenceChars, index+1, nASC, index, secondLength);
		}
		alignedSequenceChars = nASC;

		//remove index from sequenceNames
		String[] nSN = new String[sequenceNames.length-1];
		System.arraycopy(sequenceNames, 0, nSN, 0, index);
		System.arraycopy(sequenceNames, index+1, nSN, index, secondLength);
		sequenceNames = nSN;

		//remove index from sequenceTaxa
		if(sequenceTaxa != null) {
			String[] nSTa = new String[sequenceTaxa.length-1];
			System.arraycopy(sequenceTaxa, 0, nSTa, 0, index);
			System.arraycopy(sequenceTaxa, index+1, nSTa, index, secondLength);
			sequenceTaxa = nSTa;
		}
		//remove index from sequenceTitles
		String[] nSTi = new String[sequenceTitles.length-1];
		System.arraycopy(sequenceTitles, 0, nSTi, 0, index);
		System.arraycopy(sequenceTitles, index+1, nSTi, index, secondLength);
		sequenceTitles = nSTi;
	}

	public String getSequenceString(int index) {
		String r = null;
		if(index < alignedSequenceChars.length) {
			r = new String(alignedSequenceChars[index]);
		} else {
			System.out.println("SequenceAlignment.getSequenceString() requested " +
					"sequence at index " + index + ", but this alignment " +
					"only has " + alignedSequenceChars.length + 
			" aligned sequences.");
		}
		return r;
	}

	/**
	 * Checks the lengths of sequences. If they are not
	 * already the same length, they will be padded
	 * at the end with gap characters to make them
	 * all the same length.
	 *
	 */
	private void padSequencesIfNeeded() {
		boolean allSameLength = true;
		int firstSequenceLength = alignedSequenceChars[0].length;

		for(int i = 0; i < alignedSequenceChars.length && allSameLength; i++) {
			allSameLength = firstSequenceLength == alignedSequenceChars[i].length;
		}

		if(!allSameLength) {
			//find the maximum length, then pad everything to that length
			int maxLength = firstSequenceLength;
			for(int i = 1; i < alignedSequenceChars.length; i++) {
				int thisLength = alignedSequenceChars[i].length;
				if(thisLength > maxLength) {
					maxLength = thisLength;
				}
			}

			for(int i = 0; i < alignedSequenceChars.length; i++) {
				if(alignedSequenceChars[i].length != maxLength) {
					char[] forPadding = new char[maxLength];
					Arrays.fill(forPadding, GAP_CHAR);
					System.arraycopy(alignedSequenceChars[i],
							0, forPadding, 0, 
							alignedSequenceChars[i].length);
					alignedSequenceChars[i] = forPadding;
				}
			}
		}
	}

	/**
	 * Returns the number of sequences in this alignment.
	 * 
	 * @return
	 */
	public int getNTax() {
		return alignedSequenceChars.length;
	}

	/**
	 * Returns the length of this alignment.
	 * 
	 * @return
	 */
	public int getLength() {
		int r = 0;
		if(alignedSequenceChars.length > 0 && alignedSequenceChars[0] != null) {
			r = alignedSequenceChars[0].length;
		}
		return r;
	}

	/**
	 * Returns the character at the given position in the 
	 * alingmnet matrix.
	 * 
	 * @param seq
	 * @param pos
	 * @return
	 */
	public char getChar(int seq, int pos) {
		char r = 0;
		r = alignedSequenceChars[seq][pos];
		return r;
	}

	public String[] getSequenceTitles() {
		return sequenceTitles;
	}

	public String getSequenceTitle(int index) {
		String r = null;
		if(sequenceTitles.length > index) {
			r = sequenceTitles[index];
		}
		return r;
	}

	/**
	 * Returns a HashMap that maps sequence names to taxon names.
	 * 
	 * @return
	 */
	public HashMap getSequenceNameToTaxonMap() {
		HashMap r = new HashMap();
		String[] taxa = getTaxa();
		for(int i = 0; i < sequenceNames.length; i++) {
			r.put(sequenceNames[i], taxa[i]);
		}

		return r; 
	}

	public String[] getDistinctTaxa() {
		String[] r = null;
		String[] taxa = getTaxa();
		HashSet uniqueSet = new HashSet((int) (taxa.length * 1.25));
		for(int i = 0; i < taxa.length; i++) {
			uniqueSet.add(taxa[i]);
		}
		r = new String[uniqueSet.size()];
		uniqueSet.toArray(r);
		return r;
	}

	public String[] getTaxa() {
		if(sequenceTaxa == null || sequenceTaxa.length != sequenceTitles.length) {
			//try to get taxon names from titles
			sequenceTaxa = FastaUtilities.getTaxaFromTitles(sequenceTitles);
		}
		if(sequenceTaxa != null) {
			if(alignedSequenceChars == null) {
				System.out.println("SequenceAlignment.getTaxa() alignedSequenceChars is null");
			} else if(sequenceTaxa.length != alignedSequenceChars.length) {
				System.out.println("SequenceAlignment.getTaxa() mismatch between" +
						" number of taxa (" + sequenceTaxa.length + ") and " +
						"number of sequences (" + 
						alignedSequenceChars.length + ")");
			}
		}
		return sequenceTaxa;
	}

	public int getTaxonIndex(String taxon) {
		int r = -1;
		for(int i = 0; i < sequenceTaxa.length && r < 0; i++) {
			if(sequenceTaxa[i] != null && sequenceTaxa[i].equals(taxon)) {
				r = i;
			}
		}
		return r;
	}

	public int getSequenceIndex(String name) {
		int r = -1;
		for(int i = 0; i < sequenceNames.length; i++) {
			if(name.equals(sequenceNames[i])) {
				r = i;
				break;
			}
		}
		return r;
	}

	/**
	 * Returns a fasta format version of the alignment.
	 * Currently, each sequence is on one line. I need to add 
	 * code to split them into some number of characters 
	 * per line. This is still valid fasta, it's just not ideal.
	 * @return
	 */
	public String getAlignmentAsFasta() {
		StringBuffer sb = new StringBuffer();
		for(int i = 0; i < alignedSequenceChars.length; i++) {
			sb.append(">");
			sb.append(getSequenceTitle(i));
			sb.append("\n");
			sb.append(alignedSequenceChars[i]);
			sb.append("\n");
		}

		return sb.toString();
	}

	public String getAlignmentAsClustalW() {
		int maxNameLength = 50;
		int sequenceCharactersPerLine = 60;

		StringBuffer sb = new StringBuffer();
		String[] sequenceNames = getSequenceTitles();
		String[] truncatedNames = new String[sequenceNames.length];
		int maxTruncatedNameLength = 0;
		int minTruncatedNameLength = sequenceNames[0].length();

		for(int i = 0; i < truncatedNames.length; i++) {
			int spaceIndex = sequenceNames[i].indexOf(" ");
			int truncIndex = maxNameLength;
			if(spaceIndex >= 0 && spaceIndex < truncIndex) {
				truncIndex = spaceIndex;
			}

			truncatedNames[i] = sequenceNames[i].substring(0, truncIndex);
			maxTruncatedNameLength = Math.max(maxTruncatedNameLength, 
					truncatedNames[i].length());
			minTruncatedNameLength = Math.min(minTruncatedNameLength,
					truncatedNames[i].length());
		}

		//pad all sequence names to the same final length
		int maxPadLength = maxTruncatedNameLength - minTruncatedNameLength;
		char[] padBuffer = new char[maxPadLength];
		Arrays.fill(padBuffer, ' ');
		for(int i = 0; i < truncatedNames.length; i++) {
			int thisPadLength = 
				maxTruncatedNameLength - truncatedNames[i].length();
			String padString = new String(padBuffer, 0, thisPadLength);
			truncatedNames[i] = truncatedNames[i] + padString;
		}


		String spacer = "\t";

		int blocks = getLength() / sequenceCharactersPerLine;

		if(debug) {
			System.out.println("alignment length: " + getLength());
			System.out.println("blocks: " + blocks);
		}

		sb.append("CLUSTAL W (1.81) multiple sequence alignment\n\n");

		for(int i = 0; i < blocks; i++) {
			int blockStart = sequenceCharactersPerLine*i;
			int blockEnd = sequenceCharactersPerLine*(i+1);
			for(int j = 0; j < sequenceNames.length; j++) {
				sb.append(truncatedNames[j]);
				sb.append(spacer);
				sb.append(getSequenceString(j).substring(blockStart, blockEnd));
				sb.append("\n");
			}
			sb.append("\n");
		}

		if(debug) {
			System.out.println("SequenceAlignment.getAlignmentAsClustalW()");
			System.out.println(sb.toString());
		}
		return sb.toString();
	}

	/**
	 * Returns the alignment in phylip format. The taxon
	 * names are used, rather than the full titles.
	 * @return
	 */
	public String getAlignmentAsExtendedPhylipUsingTaxonNames() {
		String r = "";
		StringBuffer sb = new StringBuffer();
		sb.append(getNTax());
		sb.append(" ");
		sb.append(getLength());
		sb.append("\n");

		//create fake taxon names
		String[] paddedTaxa = new String[sequenceTaxa.length];

		//determine length of padding string, so it is one longer than the 
		//longest taxon name
		int longestTaxon = 0;
		for(int i = 0; i < sequenceTaxa.length; i++) {
			longestTaxon = Math.max(longestTaxon, sequenceTaxa[i].length());
		}
		longestTaxon++;
		char[] padChars = new char[longestTaxon];
		Arrays.fill(padChars, ' ');
		String pad = new String(padChars);
		for(int i = 0; i < paddedTaxa.length; i++) {
			paddedTaxa[i] = sequenceTaxa[i];
			paddedTaxa[i] = paddedTaxa[i] + pad.substring(paddedTaxa[i].length());
		}
		for(int i = 0; i < getNTax(); i++) {
			sb.append(paddedTaxa[i]);
			sb.append(getSequenceString(i));
			sb.append("\n");
		}
		r = sb.toString();
		return r;

	}

	/**
	 * Returns the alignment in phylip format. Sequence titles are used.
	 * @return
	 */
	public String getAlignmentAsExtendedPhylipUsingSequenceNames() {

		String r = "";
		StringBuffer sb = new StringBuffer();
		sb.append(getNTax());
		sb.append(" ");
		sb.append(getLength());
		sb.append("\n");

		//create trimmed names. These are the first token in the name
		String[] trimmedNames = new String[sequenceTitles.length];
		for(int i = 0; i < trimmedNames.length; i++) {
			trimmedNames[i] = sequenceTitles[i].split("\\s+")[0];
		}

		//create padded sequence names
		String[] paddedNames = new String[trimmedNames.length];

		//determine length of padding string, so it is one longer than the 
		//longest taxon name
		int longestTaxon = 0;
		for(int i = 0; i < trimmedNames.length; i++) {
			longestTaxon = Math.max(longestTaxon, trimmedNames[i].length());
		}
		longestTaxon++;
		char[] padChars = new char[longestTaxon];
		Arrays.fill(padChars, ' ');
		String pad = new String(padChars);
		for(int i = 0; i < paddedNames.length; i++) {
			paddedNames[i] = trimmedNames[i];
			paddedNames[i] = paddedNames[i] + pad.substring(paddedNames[i].length());
		}
		for(int i = 0; i < getNTax(); i++) {
			sb.append(paddedNames[i]);
			sb.append(getSequenceString(i));
			sb.append("\n");
		}
		r = sb.toString();
		return r;

	}

	/**
	 * Returns the alignment in phylip format. 'Fake' taxon
	 * names are used, rather than the full titles, or the
	 * real taxon names.
	 * @return
	 */
	public String getAlignmentAsPhylip() 
	{
		String r = "";
		StringBuffer sb = new StringBuffer();
		sb.append(getNTax());
		sb.append(" ");
		sb.append(getLength());
		sb.append("\n");

		//create fake taxon names
		String[] fakeTaxa = new String[sequenceTaxa.length];
		String pad = "          ";
		for(int i = 0; i < fakeTaxa.length; i++) {
			fakeTaxa[i] = "Taxon_" + i;
			fakeTaxa[i] = fakeTaxa[i] + pad.substring(fakeTaxa[i].length());
		}
		for(int i = 0; i < getNTax(); i++) {
			sb.append(fakeTaxa[i]);
			sb.append(getSequenceString(i));
			sb.append("\n");
		}
		r = sb.toString();
		return r;
	}
	/**
	 * Returns a string representation of this alignment. The format
	 * is one line per sequence: title-tab-sequence-newline
	 */
	public String toString() {
		String r = null;
		StringBuffer sb = new StringBuffer();

		int ntax = getNTax();
		for(int i = 0; i < ntax; i++) {
			sb.append(getSequenceTitle(i));
			sb.append("\t");
			sb.append(getSequenceString(i));
			sb.append("\n");
		}

		r = sb.toString();
		return r;
	}

	public boolean hasTaxonNames() {
		boolean r = true;
		if(sequenceTaxa == null) {
			r = false;
		} else {
			for(int i = 0; r && i < sequenceTaxa.length; i++) {
				if(sequenceTaxa[i] == null || sequenceTaxa[i].length() == 0) {
					r = false;
				}
			}
		}
		return r;
	}

	public String getName() {
		if(name == null) {
			name = "unnamed_alignment";
		}
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	/**
	 * Returns a single column of the alignment.
	 * If the specified index is not valid, null is returned.
	 * 
	 * @param columnIndex
	 * @return
	 */
	public char[] getColumn(int columnIndex) {
		char[] r = null;
		if(columnIndex > -1 && columnIndex < getLength()) {
			r = new char[alignedSequenceChars.length];
			for(int i = 0; i < r.length; i++) {
				r[i] = alignedSequenceChars[i][columnIndex];
			}
		}
		return r;
	}

	public Bipartition[] getBipartitionsForColumn(int index) {
		Bipartition[] r = null;
		if(columnBipartitions == null) {
			determineColumnBipartitions();
		}
		r = columnBipartitions[index];

		return r;
	}

	int[] getMinimumStepsPerSite() {
		if(minimumStepsPerSite == null || 
				minimumStepsPerSite.length != getLength()) {
			minimumStepsPerSite = new int[getLength()];
			for(int i = 0; i < minimumStepsPerSite.length; i++) {
				minimumStepsPerSite[i] = getBipartitionsForColumn(i).length -1;
			}
		}
		return minimumStepsPerSite;
	}

	public int[] getColumnIndicesForBipartition(Bipartition b) {
		int[] r = null;
		IntList indexList = new IntList();
		return r;
	}

	/**
	 * Returns Bipartitions for all columns that have no more than
	 * maxState different character states.
	 * 
	 * @param maxState
	 * @return
	 */
	public Bipartition[] getBipartitionsForColumns(int maxState) {
		Bipartition[] r = null;
		ArrayList bipartitions = new ArrayList();
		int columnCount = getLength();

		for(int i = 0; i < columnCount; i++) {
			Bipartition[] columnBiparts = getBipartitionsForColumn(i);
			if(columnBiparts == null) {
			}
			if(columnBiparts != null && columnBiparts.length <= maxState) {
				for(int j = 0; j < columnBiparts.length; j++) {
					bipartitions.add(columnBiparts[j]);
				}
			}
		}
		r = new Bipartition[bipartitions.size()];
		bipartitions.toArray(r);
		return r;
	}

	private void determineColumnBipartitions() {
		columnBipartitions = new Bipartition[getLength()][];

		for(int i = 0; i < columnBipartitions.length; i++) {
			columnBipartitions[i] = 
				getBipartitionsForAlignmentColumn(getColumn(i));
		}
	}

	/**
	 * Calculates and stores the character state arrays for all columns
	 * in the alignment.
	 */
	private void determineColumnCharacterClasses() {
		columnCharacterClasses = new int[getLength()][];

		for(int i = 0; i < columnCharacterClasses.length; i++) {
			columnCharacterClasses[i] = 
				getCharacterClassAssignmentsForAlignmentColumn(getColumn(i));
		}
	}

	/**
	 * Returns an array with character class assignments for the specified
	 * column in the alignment. The returned array has an entry for each 
	 * sequence in the alignment. Sequences with the same character in the
	 * column will have the same value in the returned array. The values 
	 * are bit flags, in ascending order. The value for missing data is 0.
	 * @param index
	 * @return
	 */
	public int[] getCharacterClassAssignmentsForColumn(int index) {
		int[] r = null;
		if(columnCharacterClasses == null) {
			determineColumnCharacterClasses();
		}
		r = columnCharacterClasses[index];
		return r;
	}

	/**
	 * Returns an array the same length as the column array. Values in the
	 * array indicate the character class for that position in the column.
	 * The values are bit flags, so they will all be powers of 2. Missing data
	 * columns are given a value of 0. Gaps are treated the same as other 
	 * characters.
	 * 
	 * @param column
	 * @return
	 */
	private int[] getCharacterClassAssignmentsForAlignmentColumn(char[] column) {
		int[] r = new int[column.length];
		int currentClass = 1;
		int nextClass = 2;
		int index = 0;

		while(column[index] == SequenceAlignment.MISSING_CHAR ||
				column[index] == SequenceAlignment.GAP_CHAR) {
			index++;
		}
		r[index] = currentClass;

		index++;
		for(;index < r.length; index++) {
			//ignore missing characters
			if(column[index] != SequenceAlignment.MISSING_CHAR) {
				//set to -1 so we can tell later if no class was found 
				r[index] = -1;
				//see if this character already has a class
				for(int j = 0 ; j < index; j++) {
					if(column[index] == column[j]) {
						//this class is the same as for element j
						r[index] = r[j];
					}
				}
				if(r[index] == -1) {
					//no existing class was found for this character, 
					//so it defines a new class
					r[index] = nextClass;
					nextClass <<= 1;
				}
			}
		}

		return r;
	}

	private Bipartition[] getBipartitionsForAlignmentColumn(char[] column) {
		Bipartition[] r = null;

		int[] characterClassAssignments = new int[column.length];
		int currentClass = 0;
		int nextClass = 1;
		int index = 0;

		//start with class of -1 for everything
		Arrays.fill(characterClassAssignments, -1);

		while(index < column.length && 
				(column[index] == SequenceAlignment.MISSING_CHAR ||
						column[index] == SequenceAlignment.GAP_CHAR)) {
			index++;
		}

		if(index >= column.length) {
		} else {
			characterClassAssignments[index] = currentClass;
			char currentClassChar = column[index];

			ExtendedBitSet participatingTaxa = new ExtendedBitSet();
			index++;
			for(;index < characterClassAssignments.length; index++) {
				//ignore missing characters and paps
				if(column[index] != SequenceAlignment.GAP_CHAR && 
						column[index] != SequenceAlignment.MISSING_CHAR) {
					participatingTaxa.set(index);
					//see if this character already has a class
					for(int j = 0 ; j < index; j++) {
						if(column[index] == column[j]) {
							//this class is the same as for element j
							characterClassAssignments[index] = characterClassAssignments[j];
						}
					}
					if(characterClassAssignments[index] == -1) {
						//no existing class was found for this character, 
						//so it defines a new class
						characterClassAssignments[index] = nextClass;
						nextClass++;
					}
				}
			}

			//all elements in the column have now been assigned to a class.
			//create the Bipartitions representing each class. The total number
			//of bipartitions currentClass+1
			ExtendedBitSet[] bitsets = new ExtendedBitSet[nextClass];
			//initialize Bipartitions
			for(int i = 0; i < bitsets.length; i++) {
				bitsets[i] = new ExtendedBitSet();
			}

			for(int i = 0; i < characterClassAssignments.length; i++) {
				if(characterClassAssignments[i] >= 0) {
					bitsets[characterClassAssignments[i]].set(i);
				}
			}

			r = new Bipartition[bitsets.length];
			for(int i = 0; i < r.length; i++) {
				r[i] = new Bipartition(bitsets[i], column.length);
				r[i].setParticipatingTaxonSet(participatingTaxa);
			}

			HashSet duplicateRemovalSet= new HashSet();
			for(int i = 0; i < r.length; i++) {
				duplicateRemovalSet.add(r[i]);
			}

			if(duplicateRemovalSet.size() != r.length) {
				r = new Bipartition[duplicateRemovalSet.size()];
				duplicateRemovalSet.toArray(r);
			}

			if(false) { //debug printing stuff
				String columnString = new String(column);
				System.out.println(columnString);
				StringBuffer classBuffer = new StringBuffer();
				for(int i = 0; i < characterClassAssignments.length; i++) {
					if(characterClassAssignments[i] == -1) {
						classBuffer.append('-');
					} else {
						classBuffer.append(characterClassAssignments[i]);
					}
				}
				System.out.println(classBuffer);
				for(int i = 0; i < r.length; i++) {
					System.out.println(r[i].getString());
				}
			}
		}
		return r;
	}

	public boolean isRepresentative() {
		boolean r = true;
		HashSet taxonSet = new HashSet();
		String[] taxa = getTaxa();
		for(int i = 0; r && i < taxa.length; i++) {
			if(taxonSet.contains(taxa[i])) {
				r = false;
			} else {
				taxonSet.add(taxa[i]);
			}
		}
		return r;
	}

	/**
	 * Gets an alignment with the columns beginning at from (inclusive)
	 * and ending at to (exclusive)
	 * @param from
	 * @param to
	 * @return
	 */
	public SequenceAlignment getSubAlignment(int from, int to) {
		SequenceAlignment r = null;
		int length = to-from;
		char[][] columns = new char[length][];
		for(int i = 0; i < columns.length; i++) {
			columns[i] = getColumn(from+i);
		}

		//transpose the columns
		char[][] subAlignment = new char[getNTax()][columns.length];
		for(int i = 0; i < subAlignment.length; i++) {
			for(int j = 0; j < subAlignment[i].length; j++) {
				subAlignment[i][j] = columns[j][i];
			}
		}
		r = new SequenceAlignment(subAlignment, getSequenceTitles());
		r.setTaxa(getTaxa());
		String subName = getName()+"_" + from + "-" + to;
		r.setName(subName);
		return r;
	}

	/**
	 * Returns a version of this alignment with the columns randomly
	 * shuffled.
	 * @return
	 */
	public SequenceAlignment getShuffledAlignment() {
		SequenceAlignment r = null;
		int length = getLength();
		char[][] shuffledColumns = new char[length][];
		int[] newColumnOrder = new int[length];
		boolean[] columnAssigned = new boolean[length];

		for(int i = 0; i < length; i++) {
			int nextIndex = random.nextInt(length);
			while(columnAssigned[nextIndex]) {
				nextIndex = random.nextInt(length);
			}
			newColumnOrder[i] = nextIndex;
			columnAssigned[nextIndex] = true;
		}

		for(int i = 0; i < shuffledColumns.length; i++) {
			shuffledColumns[i] = getColumn(newColumnOrder[i]);
		}

		//now transpose the shuffled columns, 
		char[][] shuffledAlignment = new char[getNTax()][getLength()];
		for(int i = 0; i < shuffledAlignment.length; i++) {
			for(int j = 0; j < shuffledAlignment[i].length; j++) {
				shuffledAlignment[i][j] = shuffledColumns[j][i];
			}
		}

		r = new SequenceAlignment(shuffledAlignment, getSequenceTitles());
		r.setTaxa(getTaxa());
		r.setName(getName()+ "_shuffled");

		return r;
	}

	/**
	 * Generates a shuffled alignment using sampling with or 
	 * without replacement. Used with replacement, this method 
	 * generates alignments suitable for use in bootstrap
	 * analyses.
	 * @return
	 */
	public SequenceAlignment getSampleWithReplacement(boolean withReplacement) {
		SequenceAlignment r = null;
		int length = getLength();
		char[][] shuffledColumns = new char[length][];
		int[] newColumnOrder = new int[length];
		boolean[] columnAssigned = new boolean[length];

		for(int i = 0; i < length; i++) {
			int nextIndex = random.nextInt(length);
			while(!withReplacement && columnAssigned[nextIndex]) {
				nextIndex = random.nextInt(length);
			}
			newColumnOrder[i] = nextIndex;
			columnAssigned[nextIndex] = true;
		}

		for(int i = 0; i < shuffledColumns.length; i++) {
			shuffledColumns[i] = getColumn(newColumnOrder[i]);
		}

		//now transpose the shuffled columns, 
		char[][] shuffledAlignment = new char[getNTax()][getLength()];
		for(int i = 0; i < shuffledAlignment.length; i++) {
			for(int j = 0; j < shuffledAlignment[i].length; j++) {
				shuffledAlignment[i][j] = shuffledColumns[j][i];
			}
		}

		r = new SequenceAlignment(shuffledAlignment, getSequenceTitles());
		r.setTaxa(getTaxa());
		r.setName(getName()+ "_shuffled");

		return r;
	}

	public int[] getStepsPerSite() {
		return stepsPerSite;
	}

	public void setStepsPerSite(int[] stepsPerSite) {
		if(stepsPerSite.length != getLength()) {
			System.out.println("Wrong size of stepsPerSite array for "
					+ getName() + " should be " + getLength() +
					" but is: " + stepsPerSite.length);
		}
		this.stepsPerSite = stepsPerSite;
	}

	public double[][] getPairwiseScores(int[][] scoreMatrix) {
		double[][] r = null;
		//convert sequences to ints
		int[][] sequenceInts = new int[getNTax()][];
		padSequencesIfNeeded();
		for(int i = 0; i < sequenceInts.length; i++) {
			sequenceInts[i] = 
				AlignmentUtilities.
				convertAminaAcidSequenceToInts(getSequenceString(i));
		}

		//calculate score for each pair of sequences, using the scorMatirx
		r = new double[sequenceInts.length][sequenceInts.length];
		for(int i = 0; i < r.length; i++) {
			Arrays.fill(r[i], -1f);
			for(int j = 0; j < r[i].length; j++) {
				//calculate the score for the pair (i,j) by summing the
				//scores for each position
				for(int k = 0; k < sequenceInts[i].length; k++) {
					r[i][j] += 
						scoreMatrix[sequenceInts[i][k]][sequenceInts[j][k]];
				}
			}
		}
		return r;
	}

	/**
	 * Returns the proportion of positions in this alignment
	 * that are either Gap or Missing characters.
	 * @return
	 */
	public double getProportionGapOrMissing() {
		double r = 0;
		int totalChars = alignedSequenceChars.length * alignedSequenceChars[0].length;
		int gapOrMissing = 0;
		for(int i = 0; i < alignedSequenceChars.length; i++) {
			for(int j = 0; j < alignedSequenceChars[i].length; j++) {
				switch(alignedSequenceChars[i][j]) {
				case GAP_CHAR:
					gapOrMissing++;
					break;
				case MISSING_CHAR:
					gapOrMissing++;
					break;
				}
			}
		}
		r = (double)gapOrMissing/(double) totalChars;
		return r;
	}
}
