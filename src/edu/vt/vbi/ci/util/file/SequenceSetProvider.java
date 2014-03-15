package edu.vt.vbi.ci.util.file;

import java.util.Iterator;

/**
 * Implemented by classes that provide access to multiple
 * FastaSequenceSet. Examples include genome .faa files
 * and homology groups.
 * 
 * @author enordber
 *
 */
public interface SequenceSetProvider {

	public abstract int getSequenceSetCount();
	public abstract FastaSequenceSet[] getAllSequenceSets();
	public abstract FastaSequenceSet getSequenceSet(int index);
	public abstract Iterator getSynchronizedIterator();
	public abstract void removeNonRepresentative();
	public abstract void useTaxonIncludeFile(String taxonIncludeFileName);
}
