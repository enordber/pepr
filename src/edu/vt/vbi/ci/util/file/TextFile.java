package edu.vt.vbi.ci.util.file;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;

import edu.vt.vbi.ci.util.LongList;


/**
 * Represents a file containing ASCII text. Allows
 * access to the lines in the file by index, or to
 * all lines. Entire file contents may be obtained
 * as a single String, or as an array of Strings,
 * with one String per line.
 * 
 * @author ericnordberg
 */
public class TextFile {
	/*
	 * \n
	 */
	public final static byte ASCII_NL = 10;
	
	private boolean fileIsOpen = false;
	private RandomAccessFile raFile;
	private File file;
	
	/**
	 * Stores the index of each newline ('\n')
	 * in the file. This array will be in sorted
	 * (ascending) order.
	 */
	private long[] newlineIndices;
	
	public TextFile(File f) throws IOException {
		file = f;
		findNewlines();
	}
	
	public TextFile(String fileName) throws IOException {
		file = new File(fileName);
		findNewlines();
	}

	public int getLineCount() {
	    return newlineIndices.length;	
	}
	
	/**
	 * Returns the specified line from the file.
	 * @throws IOException
	 */
	public String getLine(int index) throws IOException {
		String r = null;
		boolean wasClosed =	openFile();

		char cr = 13;
		
		long startIndex;
		long endIndex;
		
		if(index == 0) {
			startIndex = 0;
		} else {
			//start one after the previous newline
			startIndex = newlineIndices[index-1] + 1;
		}
		endIndex = newlineIndices[index];
		
		//there will be a problem with the following if the length of the line
		//is too large for an int. This seems highly unlikely.
		byte[] lineBytes = new byte[(int)(endIndex - startIndex)];
		
		raFile.seek(startIndex);
		raFile.readFully(lineBytes);

		if(wasClosed) {
			closeFile();
		}

		r = new String(lineBytes);
		
		if(lineBytes.length > 0 && lineBytes[lineBytes.length-1] == cr) {
			r = r.substring(0, r.length()-1);
		}
		
		return r;
	}
	
	
	/**
	 * Gets the specified subset of lines from the file.
	 *  
	 * @param from inclusive
	 * @param to exclusive
	 * @return
	 * @throws IOException
	 */
	public String[] getLines(int from, int to) throws IOException {
	    String[] r = null;
	    boolean wasClosed = openFile();
	    r = new String[to - from];
	    int currentIndex = 0;
	    for(int i = from; i < to; i++) {
	    	    r[currentIndex] = getLine(i);
	    	    currentIndex++;
	    }
	    if(wasClosed) {
	        closeFile();
	    }
	    return r;
	}
	
	public String[] getAllLines() throws IOException {
		String[] r = null;
		boolean wasClosed = openFile();
		int lineCount = getLineCount();
		r = new String[lineCount];
		
		for(int i = 0; i < lineCount; i++) {
			r[i] = getLine(i);
		}
		if(wasClosed) {
		    closeFile();
		}
		return r;
	}
	
	/**
	 * Returns complete file contents as a single String,
	 * with newlines at the end of each line.
	 */
	public String toString() {
		StringBuffer sb = new StringBuffer();
		boolean wasClosed = !fileIsOpen;
		try {
			wasClosed = openFile();
			int lineCount = getLineCount();

			for(int i = 0; i < lineCount; i++) {
				sb.append(getLine(i));
				sb.append("\n");
			}
		} catch(IOException ioe) {
			ioe.printStackTrace();
		}
		if(wasClosed) {
			try {
				closeFile();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}			
		return sb.toString();

	}
	
	/**
	 * Finds the index of each newline (\n) in the file,
	 * and populates the nelineIndices array with these 
	 * values.
	 * 
	 * @throws IOException
	 */
	private void findNewlines() throws IOException {
		boolean wasClosed =	openFile();
		
		LongList nlIndices = new LongList();
		long byteBufferSize = raFile.length();
		
		//use a byte[] as a searchBuffer. Fill this array and
		//look through it for value. When done looking through
		//array, refill it from bb and look again. Do this until
		//the end of bb is reached.
		
		//use a buffer size of 128K, unless the file is
		//smaller than that
		long searchBufferSize = Math.min(1024*128, byteBufferSize);

		byte[] searchBuffer = new byte[(int) searchBufferSize];
		int fullBufferSearches = 0;
		if(searchBufferSize > 0) {
			fullBufferSearches = (int) (byteBufferSize / searchBufferSize);
		}

		raFile.seek(0);
		for(int i = 0; i < fullBufferSearches; i++) {
			raFile.readFully(searchBuffer);
			for(int j = 0; j < searchBufferSize; j++){
				if(searchBuffer[j] == ASCII_NL) {
					long index = ((i * searchBufferSize) + j);
					nlIndices.add(index);
				}
			}
		}
		
		searchBuffer = null;
		
		//search the remaining bytes
		if(searchBufferSize == 0) {
//			avoid division by zero
			searchBufferSize = 1;
		}
		byte[] remainder = new byte[(int) (byteBufferSize % searchBufferSize)];
		raFile.readFully(remainder);
		for(int i = 0; i < remainder.length; i++) {
			if(remainder[i] == ASCII_NL) {
				long index = (searchBufferSize * fullBufferSearches) + i;
				nlIndices.add(index);
			} 
		}
		
		//if file does not end with a newline, assume a newline should be
		//at the end of the file
		if(nlIndices.size() == 0 || 
				nlIndices.get(nlIndices.size()-1) != raFile.length()) {
			nlIndices.add(raFile.length());
		}
		
		newlineIndices = nlIndices.getLongs();
		if(wasClosed) {
			closeFile();
		}
	}

	/**
	 * Opens the file for reading, if it is not already open.
	 * Boolean return value indicates if state of file was
	 * changed by calling this method, That is, true is returned
	 * if the file was closed before calling.
	 * @return
	 * @throws IOException
	 */
	public boolean openFile() throws IOException {
		boolean r = false;
		if(!fileIsOpen) {
		    raFile = new RandomAccessFile(file, "r");
		    fileIsOpen = true;
		    r = true;
		}
		return r;
	}
	
	/**
	 * Closes the file, if it is open.
	 * Boolean return value indicates if state of file was
	 * changed by calling this method, That is, true is returned
	 * if the file was open before calling.
	 * @return
	 * @throws IOException
	 */
	public boolean closeFile() throws IOException {
		boolean r = false;
		if(fileIsOpen) {
			raFile.close();
			raFile = null;
			fileIsOpen = false;
			r = true;
		}
		return r;
	}
	
	public File getFile() {
		return file;
	}
	
	protected long[] getNewlineIndices() {
		return newlineIndices;
	}
}
