package edu.vt.vbi.ci.util.file;

import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import edu.vt.vbi.ci.util.IntList;

/**
 * @author ericnordberg
 */
public class FastaSequenceFile extends TextFile implements FastaSequenceSet,
Comparable{

	public static final byte ASCII_GT = '>';
	
	private String name;
	private String fullFileName;
	private int[] sequenceStarts;
	private String[] titles;
	private String[] taxa;
	
	/*
	 * Many sequences with a standardized format have a unique
	 * identifier token at the start of the title line. This token 
	 * occurs between the '>' and the first space. idToIndexMap
	 * maps these tokens to their index in the sequence, for those
	 * sequences in this file having a unique id token. This is
	 * useful for speeding up the search for a sequence based on this
	 * id.
	 */
	private HashMap idToIndexMap;
	
	/**
	 * @param fileName
	 * @throws IOException
	 */
	public FastaSequenceFile(String fileName) throws IOException {
		super(fileName);
		setFileName(fileName);
		findSequenceStarts();
	}

	private void findSequenceStarts() throws IOException {
		RandomAccessFile raFile = new RandomAccessFile(getFile(), "r");
		
		IntList gtIndices = new IntList();
		long[] newlineIndices = getNewlineIndices();
		//any greater-than characters '>' should occur immediately
		//after a newline, except the first one, which should be the
		//first character
		
		long byteBufferSize = raFile.length();
		
		//use a byte[] as a searchBuffer. Fill this array and
		//look through it for value. When done looking through
		//array, refill it from bb and look again. Do this until
		//the end of bb is reached.
		
		//use a buffer of size 128K, unless the complete file
		//is smaller than this
		long searchBufferSize = Math.min(1024*128, byteBufferSize);

		byte[] searchBuffer = new byte[(int)searchBufferSize];
		long fullBufferSearches = byteBufferSize / searchBufferSize;

		raFile.seek(0);
		byte candidate = raFile.readByte();
		if(candidate == ASCII_GT) {
			gtIndices.add(0);
		}
		raFile.seek(0);

		int newlineIndex = 0;
		long startOfBuffer;
		long endOfBuffer;
		long candidateAbsIndex;
		int candidateRelIndex;
		for(int i = 0; i < fullBufferSearches; i++) {
			raFile.readFully(searchBuffer);
			startOfBuffer  = i * searchBufferSize;
		    endOfBuffer = startOfBuffer + searchBufferSize;
			for(; newlineIndices[newlineIndex] + 1 < endOfBuffer; 
			newlineIndex++){
				candidateAbsIndex = newlineIndices[newlineIndex] + 1;
				candidateRelIndex = (int)(candidateAbsIndex - startOfBuffer);
				if(searchBuffer[candidateRelIndex] == ASCII_GT) {
					gtIndices.add(newlineIndex+1);
				}
			}
		}
		
		searchBuffer = null;
		
		//search the remaining bytes
		byte[] remainder = new byte[(int)(byteBufferSize % searchBufferSize)];
		raFile.readFully(remainder);
		startOfBuffer = fullBufferSearches * searchBufferSize;
		endOfBuffer = startOfBuffer + remainder.length;
		for(; newlineIndices[newlineIndex]+1 < endOfBuffer; newlineIndex++){
			candidateAbsIndex = (int)newlineIndices[newlineIndex] + 1;
			candidateRelIndex = (int)(candidateAbsIndex - startOfBuffer);
			if(remainder[candidateRelIndex] == ASCII_GT) {
				gtIndices.add(newlineIndex+1);
			}
		}

		sequenceStarts = gtIndices.getInts();
		raFile.close();
	}
	
	/* (non-Javadoc)
	 * @see util.file.FasaSequenceSet#getSequenceCount()
	 */
	public int getSequenceCount() {
		int r = -1;
		r = sequenceStarts.length;
		
		return r;
	}
	
	/* (non-Javadoc)
	 * @see util.file.FasaSequenceSet#getTitles()
	 * This is synchronized because it opens then closes the file, so 
	 * access by multiple Threads can result in one Thread closing the file 
	 * while another is still trying to read from the file.
	 */
	public synchronized String[] getTitles() {
		if(titles == null || !titles[0].startsWith(">")) {
			titles= new String[sequenceStarts.length];
			try {
				boolean wasClosed = openFile();
				for(int i = 0; i < sequenceStarts.length; i++) {
					titles[i] = getLine(sequenceStarts[i]);
				}
				if(wasClosed) {
				    closeFile();
				}
			} catch(IOException ioe) {
				ioe.printStackTrace();
			}
		} 
		return titles;
	}
	
	/**
	 * Returns a HashMap from unique ID to index for the
	 * sequence. Non-unique IDs are not in this map.
	 * 
	 * @return HashMap
	 * 			key: ID token (between '>' and first whitespace in title line)
	 * 			value: Integer for index of sequence with the ID
	 */
	public HashMap getIDToIndexMap() {
		if(idToIndexMap == null) {
			idToIndexMap = new HashMap();
			HashSet nonUniqueIDs = new HashSet();
			String[] titles = getTitles();
			
			for(int i = 0; i < titles.length; i++) {
				//get the id token
				String idToken = getIDToken(titles[i]);
				//System.out.println("id token: " + idToken);
				if(idToIndexMap.containsKey(idToken)) {
					//this is a duplicate id.
					//remove from map, and add to nonUniqueIDs
					idToIndexMap.remove(idToken);
					nonUniqueIDs.add(idToken);
				} else if(!nonUniqueIDs.contains(idToken)) {
					//this is the first time this id has been seen.
					//consider it unique, for now, and add it to the map
					Integer index = new Integer(i);
					idToIndexMap.put(idToken, index);
				}
			}
		}
		
		return idToIndexMap;
	}
	
	/**
	 * Convenience method to get the substring of a title
	 * line between '>' and the first whitespace character. Note
	 * that the first whitespace may be the newline at the end of
	 * the title line, in which case the idToken is the entire title
	 * line minus the '>'.
	 * @param title
	 * @return
	 */
	public static String getIDToken(String title) {
		String r = null;
		int spaceIndex = title.indexOf(" ", 1);
		if(spaceIndex < 0) {
			spaceIndex = title.length();
		}
		int startIndex = 0;
		if(title.startsWith(">")) {
			startIndex = 1;
		}
		r = new String(title.substring(startIndex, spaceIndex));
		return r;
	}

	/* (non-Javadoc)
	 * @see util.file.FasaSequenceSet#getSequence(int)
	 * This is synchronized because it opens then closes the file, so 
	 * access by multiple Threads can result in one Thread closing the file 
	 * while another is still trying to read from the file.
	 */
	public synchronized String[] getSequence(int i) {
		String[] r = null;
		int firstLine = sequenceStarts[i];
		int lastLine = -1; //the name lastLine is a bit misleading. this is
		                    //the index of the line just after the last line
		                    //of this sequence
		if(i + 1 < sequenceStarts.length) {
			lastLine = sequenceStarts[i+1];
		} else {
			lastLine = getLineCount();
		}
		
		r = new String[lastLine - firstLine];
		int currentIndex = 1;
		try {
			boolean wasClosed = openFile();
			//copy title line exactly - do not remove spaces
			r[0] = getLine(firstLine); 
			for(int line = firstLine+1; line < lastLine; line++) {
				r[currentIndex] = getLine(line);
				//remove any spaces in the sequence line
				r[currentIndex] = r[currentIndex].replaceAll(" ", "");
				currentIndex++;
			}
			if(wasClosed) {
			    closeFile();
			}
		} catch(IOException ioe) {
			ioe.printStackTrace();
		}
		return r;
	}
	
	/* (non-Javadoc)
	 * @see util.file.FasaSequenceSet#getPlainSequence(int)
	 */
	public String getPlainSequence(int index) {
		String sequence = null;
		String[] seqArray = getSequence(index);
		StringBuffer seqBuffer = new StringBuffer();
		for(int i = 1; i < seqArray.length; i++) {
			seqBuffer.append(seqArray[i]);
		}
		sequence = seqBuffer.toString();
		sequence = sequence.replaceAll(" ", "");
		return sequence;
	}

	/* (non-Javadoc)
	 * @see util.file.FasaSequenceSet#getSequence(java.lang.String)
	 */
	public String[] getSequence(String searchID) {
     	String[] r = null;
     	int foundIndex = -1;
     	
     	foundIndex = getIndexOfSequence(searchID);
     	
     	if(foundIndex > -1) {
     	    r = getSequence(foundIndex);
     	}
		return r;
	}
	
	/* (non-Javadoc)
	 * @see util.file.FasaSequenceSet#getIndexOfSequence(java.lang.String)
	 */
	public int getIndexOfSequence(String searchID) {
		int foundIndex = -1;

		String idToken = getIDToken(searchID);
		HashMap idToIndex = getIDToIndexMap();
		Integer index = (Integer)idToIndex.get(idToken);
		if(index != null) {
			foundIndex = index.intValue();
		} else {
			//rapid lookup failed, so do linear search through all titles
			int offset = 1;
			if(searchID.startsWith(">")) {
				offset = 0;
			}
			String[] titles = getTitles();
			for(int i = 0; i < titles.length; i++) {
				if(titles[i] != null && titles[i].startsWith(searchID, offset)) {
					foundIndex = i;
					break;
				} else if(titles[i] != null && titles[i].indexOf(searchID) > -1) {
					foundIndex = i;
					break; 
				}
			}
		}
		return foundIndex;
	}
	
	/* (non-Javadoc)
	 * @see util.file.FasaSequenceSet#containsTitle(java.lang.String)
	 */
	public boolean containsTitle(String title) {
		boolean r = false;
		int index = getIndexOfSequence(title);
		r = index > -1;
		return r;
	}
	
	public void open() throws IOException {
		openFile();
	}
	
	public void close() throws IOException {
		closeFile();
	}
	
	/**
	 * Returns the name for this sequence set.
	 * This is the file name (excluding directory and extension).
	 */
	public String getName() {
		return name;
	}
	
	/**
	 * Creates and sets the sequence set name based on the 
	 * file name. The set name is the file name excluding
	 * directory and extension.
	 * 
	 * @param fileName
	 */
	private void setFileName(String fileName) {
		int start = fileName.lastIndexOf(System.getProperty("file.separator"));
		int end = fileName.lastIndexOf(".");
		if(start < 0) {
			start = 0;
		} else {
			start++;
		}
		
		if(end < start) {
			end = fileName.length();
		}
		setName(fileName.substring(start, end));
		fullFileName = fileName;
	}
	
	public void setName(String name) {
		this.name = name;
	}
	
	public String getFullName() {
		return fullFileName;
	}

	public String getFileName() {
		String r = getFile().getName();
		return r;
	}
	/**
	 * Returns the contents as a fasta-format String.
	 */
	public String toString() {
		StringBuffer sb = new StringBuffer();
		int count = getSequenceCount();
		for(int i = 0; i < count; i++) {
			String[] seq = getSequence(i);
			for(int j = 0; j < seq.length; j++) {
			    sb.append(seq[j]);
			    sb.append('\n');
			}
		}

		return sb.toString();
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

	public String[] getTaxa() {
		if(taxa == null) {
			taxa = FastaUtilities.getTaxaFromTitles(getTitles());
		}
		return taxa;
	}
	
	public String[] getSequenceForTaxon(String taxon) {
		String[] r = null;
		//get the first index for a sequence with this taxon
		int index = -1;
		String[] taxa = getTaxa();
		for(int i = 0; i < taxa.length; i++) {
			if(taxon.equals(taxa[i])) {
				index = i;
				break;
			}
		}
		
		//get the sequence for this index
		r = getSequence(index);
		return r;
	}
	
	public int getDistinctTaxonCount() {
		int r = -1;
		String[] taxa = getTaxa();
		HashSet distinctTaxa = new HashSet();
		for(int i = 0; i < taxa.length; i++) {
			distinctTaxa.add(taxa[i]);
		}
		r = distinctTaxa.size();
		return r;
	}

	public FastaSequenceSet getTaxonSubset(String[] taxonList) {
		FastaSequenceSet r;
		HashSet retainTaxaSet = new HashSet();
		for(int i = 0; i < taxonList.length; i++) {
			retainTaxaSet.add(taxonList[i].trim());
		}
		boolean[] keep = new boolean[getTaxa().length];

		for(int i = 0; i < keep.length; i++) {
			if(retainTaxaSet.contains(taxa[i])) {
			    keep[i] = true;
			}
		}
		
		StringBuffer subSequence = new StringBuffer();
		for(int i = 0; i < taxa.length; i++) {
			if(keep[i]) {
				String[] linesToKeep = getSequence(i);
				for(int j = 0; j < linesToKeep.length; j++) {
					subSequence.append(linesToKeep[j]);
					subSequence.append("\n");
				}
			}
		}
		
		r = new FastaSequenceSetImpl(subSequence.toString());
		r.setName(this.getName());
		return r;
	}
	
	/**
	 * Returns the combined length of all sequences in this Sequence Set.
	 * There is a lot of room for efficiency improvements in the current
	 * implementation of this method.
	 * 
	 * @return
	 */
	public int getTotalSequenceLength() {
		int r = 0;
		int sequenceCount = getSequenceCount();
		for(int i = 0; i < sequenceCount; i++) {
			r += getPlainSequence(i).length();
		}
		
		return r;
	}

	public int compareTo(Object arg0) {
		int r = 0;
		FastaSequenceFile other = (FastaSequenceFile)arg0;
		r = (int) (other.getFile().length() - this.getFile().length());
		return r;
	}

}
