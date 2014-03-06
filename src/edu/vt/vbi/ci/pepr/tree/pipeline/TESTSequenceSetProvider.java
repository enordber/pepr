package edu.vt.vbi.ci.pepr.tree.pipeline;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Random;

import edu.vt.vbi.ci.util.file.FastaSequenceFile;
import edu.vt.vbi.ci.util.file.FastaSequenceSet;
import edu.vt.vbi.ci.util.file.SequenceSetProvider;
import edu.vt.vbi.ci.util.file.TextFile;

public class TESTSequenceSetProvider implements SequenceSetProvider {

	private boolean debug = false;

	private FastaSequenceSet[] sequenceSets;
	private String[] sequenceFileNames;

	public TESTSequenceSetProvider() {
	}


	/**
	 * creates a SequenceSetProvider with all the files in the specified
	 * directory.
	 * 
	 * @param sequenceFileDir
	 * @param minTax the smallest number of sequences in a set to be loaded.
	 *                Any sequence files with less than minTax sequences will 
	 *                not be loaded.
	 */
	public TESTSequenceSetProvider(String sequenceFileDir, int minTax) {
		//get a list of all files in the directory
		File directory = new File(sequenceFileDir);
		File[] fileList = directory.listFiles();
		try {
			ArrayList sequenceSetList = new ArrayList(fileList.length); 
			//create FastaSequenceSet for each file
			for(int i = 0; i < fileList.length; i++) {
				FastaSequenceFile fsf =
					new FastaSequenceFile(fileList[i].getAbsolutePath());
				if(fsf.getSequenceCount() >= minTax) {
					sequenceSetList.add(fsf);
				}
			}
			sequenceSets = new FastaSequenceSet[sequenceSetList.size()];
			sequenceSetList.toArray(sequenceSets);
			sequenceFileNames = new String[sequenceSets.length];
			for(int i = 0; i < sequenceFileNames.length; i++) {
				sequenceFileNames[i] = sequenceSets[i].getName();
			}

		} catch(IOException ioe) {
			ioe.printStackTrace();
		}
	}

	public TESTSequenceSetProvider(FastaSequenceSet[] sequenceSets) {
		this.sequenceSets = new FastaSequenceSet[sequenceSets.length];
		System.arraycopy(sequenceSets, 0, this.sequenceSets, 0, sequenceSets.length);

		sequenceFileNames = new String[this.sequenceSets.length];
		for(int i = 0; i < sequenceFileNames.length; i++) {
			sequenceFileNames[i] = this.sequenceSets[i].getName();
		}
	}

	/**
	 * Creates a SequenceSetProvider with all the files specified
	 * in the given list.
	 * 
	 * @param fileNames
	 * @param minTax the smallest number of sequences in a set to be loaded.
	 *                Any sequence files with less than minTax sequences will 
	 *                not be loaded.
	 */
	public TESTSequenceSetProvider(String[] fileNames, int minTax) {
		try {
			ArrayList sequenceSetList = new ArrayList(fileNames.length);
			sequenceFileNames = fileNames;
			for(int i = 0; i < fileNames.length; i++) {
				FastaSequenceFile fsf = null;
				try {
					fsf = new FastaSequenceFile(fileNames[i]);
				} catch(ArithmeticException e) {

				}
				if(fsf != null && fsf.getSequenceCount() >= minTax) {
					sequenceSetList.add(fsf);
				}
				if(debug) {
					System.out.println("TESTSequenceSetProvider loaded" +
							" sequence file. sequence count: " +
							sequenceSets[i].getSequenceCount());
				}
			}
			sequenceSets = new FastaSequenceSet[sequenceSetList.size()];
			sequenceSetList.toArray(sequenceSets);
			sequenceFileNames = new String[sequenceSets.length];
			for(int i = 0; i < sequenceFileNames.length; i++) {
				sequenceFileNames[i] = sequenceSets[i].getFullName();
			}
		} catch(IOException ioe) {
			ioe.printStackTrace();
		}
	}

	/**
	 * Creates a SequenceSetProvider with all the files specified
	 * in the given list.
	 * 
	 * @param fileNames
	 * @param minTax the smallest number of sequences in a set to be loaded.
	 *                Any sequence files with less than minTax sequences will 
	 *                not be loaded.
	 *  @param maxTax the largest number of sequences in a set to be loaded.
	 *                Any sequence files with more than maxTax sequences will
	 *                not be loaded. 
	 */
	public TESTSequenceSetProvider(String[] fileNames, int minTax, int maxTax) {
		try {
			ArrayList sequenceSetList = new ArrayList(fileNames.length);
			sequenceFileNames = fileNames;
			for(int i = 0; i < fileNames.length; i++) {
				FastaSequenceFile fsf = null;
				try {
					fsf = new FastaSequenceFile(fileNames[i]);
				} catch(ArithmeticException e) {

				}
				if(fsf != null && fsf.getSequenceCount() >= minTax
						&& fsf.getSequenceCount() <= maxTax) {
					sequenceSetList.add(fsf);
				}
				if(debug) {
					System.out.println("TESTSequenceSetProvider loaded" +
							" sequence file. sequence count: " +
							fsf.getSequenceCount());
					System.out.println("retained sequence sets: " + sequenceSetList.size());
				}
			}
			sequenceSets = new FastaSequenceSet[sequenceSetList.size()];
			sequenceSetList.toArray(sequenceSets);
			sequenceFileNames = new String[sequenceSets.length];
			for(int i = 0; i < sequenceFileNames.length; i++) {
				sequenceFileNames[i] = sequenceSets[i].getFullName();
			}
		} catch(IOException ioe) {
			ioe.printStackTrace();
		}
	}

	public FastaSequenceSet[] getAllSequenceSets() {
		return sequenceSets;
	}

	/**
	 * Returns a random subset of the sequences in this set. The
	 * number of sequences in the subset is determined by the 
	 * subsetSize parameter. If subsetSize is greater than the 
	 * number of sequences in the set, then all sequences are
	 * returned.
	 * 
	 * @param subsetSize
	 * @return
	 */
	public FastaSequenceSet[] getRandomSubset(int subsetSize) {
		FastaSequenceSet[] r = null;
		subsetSize = Math.min(subsetSize, sequenceSets.length);
		boolean[] assigned = new boolean[sequenceSets.length];
		Random rand = new Random();
		for(int i = 0; i < subsetSize; i++) {
			int nextIndex = rand.nextInt(subsetSize);
			while(assigned[nextIndex]) {
				nextIndex = rand.nextInt(subsetSize);
			}
			assigned[nextIndex] = true;
		}

		return r;
	}

	public FastaSequenceSet getSequenceSet(int index) {
		return sequenceSets[index];
	}

	public int getSequenceSetCount() {
		return sequenceSets.length;
	}

	/**
	 * Returns an Iterator for the sequences in this sequence set. 
	 * The Iterators methods are synchronized.
	 */
	public Iterator getSynchronizedIterator() {
		Iterator r = new SynchronizedIterator();
		return r;
	}

	private class SynchronizedIterator implements Iterator {

		private int nextIndex = 0;

		public synchronized boolean hasNext() {
			boolean r = false;
			r = nextIndex < sequenceSets.length;
			return r;
		}

		public synchronized Object next() {
			Object r = null;
			if(nextIndex < sequenceSets.length) {
				r = sequenceSets[nextIndex];
				nextIndex++;
			}
			return r;
		}

		/**
		 * Not Implemented.
		 */
		public synchronized void remove() {
		}

	}

	public String[] getFileNames() {
		return sequenceFileNames;
	}

	/**
	 * Removes any SequenceSets that contain any taxon more than once
	 */
	public void removeNonRepresentative() {
		boolean[] keep = new boolean[sequenceSets.length];
		int keepCount = 0;

		for(int i = 0; i < sequenceSets.length; i++) {
			if(sequenceSets[i].isRepresentative()) {
				keep[i] = true;
				keepCount++;
			}
		}
		String[] newFileNames = new String[keepCount];
		FastaSequenceSet[] newSets = new FastaSequenceSet[keepCount];

		int nextIndex = 0;
		for(int i = 0; i < keep.length; i++) {
			if(keep[i]) {
				newFileNames[nextIndex] = sequenceFileNames[i];
				newSets[nextIndex] = sequenceSets[i];
				nextIndex++;
			}
		}

		sequenceFileNames = newFileNames;
		sequenceSets = newSets;
	}


	/**
	 * Loads the specified file. Each line in the file should contain the
	 * name of a toxon that is to be retained in the sequence sets
	 * in this provider. All taxa not in this file will be removed 
	 * from the sequence sets.
	 */
	public void useTaxonIncludeFile(String taxonIncludeFileName) {
		try {
			TextFile taxonFile = new TextFile(taxonIncludeFileName);
			String[] lines = taxonFile.getAllLines();
			ArrayList taxonList = new ArrayList();
			ArrayList reducedSets = new ArrayList();
			for(int i = 0; i < lines.length; i++) {
				if(lines[i].length() > 0) {
					taxonList.add(lines[i].trim());
				}
			}
			if(lines.length > taxonList.size()) {
				lines = new String[taxonList.size()];
				taxonList.toArray(lines);
			}
			for(int i = 0; i < sequenceSets.length; i++) {
				FastaSequenceSet reducedSet = 
					sequenceSets[i].getTaxonSubset(lines);
				if(reducedSet.getSequenceCount() > 1) {
					reducedSets.add(reducedSet);
				}
			}

			sequenceSets = new FastaSequenceSet[reducedSets.size()];
			reducedSets.toArray(sequenceSets);

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	/**
	 * Removes any SequenceSets that do not contain at least the specified
	 * number of taxa.
	 * 
	 * @param taxa
	 */
	public void removeSetsWithFewerTaxaThan(int taxa) {
		ArrayList nameList = new ArrayList(sequenceFileNames.length);
		ArrayList setList = new ArrayList(sequenceSets.length);

		//keep those with sufficient taxa
		for(int i= 0; i < sequenceSets.length; i++) {
			if(sequenceSets[i].getDistinctTaxonCount() >= taxa) {
				setList.add(sequenceSets[i]);
				nameList.add(sequenceFileNames[i]);
			}
		}

		sequenceFileNames = new String[nameList.size()];
		nameList.toArray(sequenceFileNames);

		sequenceSets = new FastaSequenceSet[setList.size()];
		setList.toArray(sequenceSets);
	}

	/**
	 * Returns the number of sequence sets that have at least minTaxa taxa.
	 * 
	 * @param minTaxa
	 * @return
	 */
	public int getSequenceSetCountForMinTaxa(int minTaxa) {
		int r = 0;
		int count = getSequenceSetCount();
		for(int i = 0; i < count; i++) {
			if(getSequenceSet(i).getDistinctTaxonCount() >= minTaxa) {
				r++;
			}
		}
		return r;
	}

}
