package edu.vt.vbi.ci.util.file;

import java.util.HashMap;
import java.util.HashSet;

public class FastaSequenceSetImpl implements FastaSequenceSet {
	private String originalSequenceString;
	private String name;
	private int[] titleIndices;
	private String[] titles;
	private String[] lines;
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
	

	public FastaSequenceSetImpl(String fastaSequence) {
		originalSequenceString = fastaSequence;
		lines = fastaSequence.split("\n");
		int sequenceCount = 0;
		for(int i = 0; i < lines.length; i++) {
			if(lines[i].startsWith(">")) {
				sequenceCount++;
			}
		}
		titleIndices = new int[sequenceCount];
		int indexTracker = 0;
		for(int i = 0; i < lines.length; i++) {
			if(lines[i].startsWith(">")) {
				titleIndices[indexTracker] = i;
				indexTracker++;
			}
		}
		
		titles = new String[indexTracker];
		for(int i = 0; i < titleIndices.length; i++) {
			titles[i] = lines[titleIndices[i]];
		}
	}
	
	public boolean containsTitle(String title) {
		boolean r = false;
		int index = getIndexOfSequence(title);
		r = index > -1;
		return r;
	}

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
				if(titles[i].startsWith(searchID, offset)) {
					foundIndex = i;
					break;
				} else if(titles[i].indexOf(searchID) > -1) {
					foundIndex = i;
					break; 
				}
			}
		}
		return foundIndex;
	}

	/**
	 * Returns the sequence as a single String, with no 
	 * newlines and with no title line.
	 */
	public String getPlainSequence(int index) {
		int[] firstAndLast = getLineStartAndEndForSequence(index);
		StringBuffer sb = new StringBuffer();
		
		for(int i = firstAndLast[0]+1; i < firstAndLast[1]; i++) {
			sb.append(lines[i]);
		}
		return sb.toString();
	}

	public String[] getSequence(int index) {
		String[] r = null;
		int[] firstAndLast = getLineStartAndEndForSequence(index);
		int lineCount = firstAndLast[1] - firstAndLast[0];
		r = new String[lineCount];
		int rIndex = 0;
		for(int i = firstAndLast[0]; i < firstAndLast[1]; i++) {
			r[rIndex] = lines[i];
			rIndex++;
		}
		return r;
	}

	public String[] getSequence(String searchID) {
    	String[] r = null;
     	int foundIndex = -1;
     	
     	foundIndex = getIndexOfSequence(searchID);
     	
     	if(foundIndex > -1) {
     	    r = getSequence(foundIndex);
     	}
		return r;
	}

	public int getSequenceCount() {
		return titles.length;
	}

	public String[] getTitles() {
		return titles;
	}
	
	private int[] getLineStartAndEndForSequence(int index) {
		int[] r = new int[2];
		
		r[0]= titleIndices[index];
		//the second value is the index of the first line after 
		//the desired sequence. If the desired sequence is the last
		//sequence in the set, then this value will be the number
		//of lines
		
		if(index < titleIndices.length-1) {
			r[1] = titleIndices[index+1];
		} else if(index == titleIndices.length - 1) {
			r[1] = lines.length;
		} else {
			System.out.println("**BAD LAST INDEX FOR: " + index + 
					" in FastaSequenceSetImpl.getLineStartAndEndForSequence()");
		}
		return r;
	}
	
	public String toString() {
		return originalSequenceString;
	}

	public void setName(String n) {
		name = n;
	}

	/**
	 * Returns the name of this sequence set. If no
	 * name has been set, this will return null
	 */
	public String getName() {
		return name;
	}
	
	public String getFullName() {
		return name;
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
	 * Internal convenience method to get the substring of a title
	 * line between '>' and the first whitespace character. Note
	 * that the first whitespace may be the newline at the end of
	 * the title line, in which case the idToken is the entire title
	 * line minus the '>'.
	 * @param title
	 * @return
	 */
	private String getIDToken(String title) {
		String r = null;
		int spaceIndex = title.indexOf(" ", 1);
		if(spaceIndex < 0) {
			spaceIndex = title.length();
		}
		int startIndex = 0;
		if(title.startsWith(">")) {
			startIndex = 1;
		}
		r = title.substring(startIndex, spaceIndex);

		return r;
	}

	public FastaSequenceSet getTaxonSubset(String[] lines) {
		// TODO Auto-generated method stub
		return null;
	}

	public int getTotalSequenceLength() {
		int r = 0;
		int count = getSequenceCount();
		
		for(int i = 0; i < count; i++) {
			r += getPlainSequence(i).length();
		}
		return r;
	}


}
