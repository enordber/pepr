package edu.vt.vbi.ci.util.file;

public interface FastaSequenceSet {

	public abstract int getSequenceCount();

	public abstract String[] getTitles();

	/**
	 * Returns the ith sequence as an array of 
	 * Strings. The first element in the array 
	 * is the title line.
	 * 
	 * @param i
	 * @return
	 */
	public abstract String[] getSequence(int i);

	/**
	 * Returns the sequence with the given index as a single
	 * String. This String will not contain the title line
	 * or any whitespace.
	 * 
	 * @param index
	 * @return
	 */
	public abstract String getPlainSequence(int index);

	/**
	 * Looks for a title that starts with the given
	 * search ID. If such a title is found, that sequence
	 * will be returned. If no sequence matches, null is
	 * returned. If multiple sequences match, the first
	 * one that is found is returned.
	 *  
	 * @param searchID
	 * @return
	 */
	public abstract String[] getSequence(String searchID);

	public abstract int getIndexOfSequence(String searchID);

	/**
	 * Returns true if this sequence file contains a title
	 * matching the given title. The entire title need not
	 * match, but the actual title must start with the given
	 * String.
	 * 
	 * @param title
	 * @return
	 */
	public abstract boolean containsTitle(String title);
	
	/**
	 * Returns the name of this sequence set.
	 * @return
	 */
	public abstract String getName();

	/**
	 * Returns true if each taxon name occurs only once.
	 * @return
	 */
	public abstract boolean isRepresentative();

	public abstract FastaSequenceSet getTaxonSubset(String[] lines);

	public abstract void setName(String name);

	public abstract int getDistinctTaxonCount();

	public abstract String getFullName();
	
	public abstract int getTotalSequenceLength();

	public abstract String[] getTaxa();
}