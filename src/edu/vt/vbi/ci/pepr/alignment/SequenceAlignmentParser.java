package edu.vt.vbi.ci.pepr.alignment;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import edu.vt.vbi.ci.util.file.FastaSequenceFile;
import edu.vt.vbi.ci.util.file.FastaSequenceSet;
import edu.vt.vbi.ci.util.file.FastaUtilities;
import edu.vt.vbi.ci.util.file.TextFile;


/**
 * Contains static methods for parsing various sequence
 * alignment file formats and generating a SequenceAlignment
 * Object.
 * @author enordber
 *
 */
public class SequenceAlignmentParser {
	
	public static SequenceAlignment parseNexusAlignmentFile(String filename) throws IOException {
		SequenceAlignment r = null;
		String whitespaceRegEx = "\\s+";
		TextFile tf = new TextFile(filename);
		ArrayList taxa = new ArrayList();
		HashMap taxaToSequence = new HashMap();
		String[] fileLines = tf.getAllLines();
		//move through lines until matrix begins
		int i = 0;
		for(;!fileLines[i].trim().equalsIgnoreCase("matrix"); i++) {}
		i++;
		//until semicolon is encountered, get sequences
		String taxon = null;
		String sequence = null;
		for(;fileLines[i].indexOf(';') < 0; i++) {
			String[] splits = fileLines[i].split(whitespaceRegEx);
			if(splits.length == 1) {
				sequence = splits[0];
			} else if(splits.length == 2) {
				taxon = splits[0];
				sequence = splits[1];
			} else if(splits.length == 3) {
				taxon = splits[1];
				sequence = splits[2];
			}
			
			if(!taxa.contains(taxon)) {
				taxa.add(taxon);
			}
			StringBuffer seqBuffer = (StringBuffer)taxaToSequence.get(taxon);
			if(seqBuffer == null) {
				seqBuffer = new StringBuffer();
			}
			seqBuffer.append(sequence);
			taxaToSequence.put(taxon, seqBuffer);
		}
		
		String[] taxaArray = new String[taxa.size()];
		taxa.toArray(taxaArray);
		
		String[] sequenceArray = new String[taxaArray.length];
		for(i = 0; i < taxaArray.length; i++) {
			sequenceArray[i] = taxaToSequence.get(taxaArray[i]).toString();
		}
		
		r = new SequenceAlignment(sequenceArray, taxaArray, taxaArray);
		return r;
	}
	
	public static SequenceAlignment parseClustalAlignment(String alignment) {
		SequenceAlignment r = null;
		String[] alignmentLines = alignment.split("\\n");
		ArrayList titleList = new ArrayList();

		HashMap titleToStringBuffer = new HashMap();
		//get the list of all sequence names (titles)
		for(int i = 0; i < alignmentLines.length; i++) {
			String[] splits = alignmentLines[i].split("\\s+");
			if(splits.length > 1 && !splits[0].toLowerCase().startsWith("clustal")) {
				//if this name is not known, add it to the list.
				//if the name is known, assume that means that all names
				//have been encountered, so exit the loop
				String title = splits[0];
				if(!titleList.contains(title)) {
					titleList.add(title);
				}
				if(!titleToStringBuffer.containsKey(title)) {
					titleToStringBuffer.put(title, new StringBuffer());
				}
				StringBuffer sb = 
					(StringBuffer) titleToStringBuffer.get(title);
				sb.append(splits[1]);
			}
		}
		
		String[] titles = new String[titleList.size()];
		titleList.toArray(titles);
		String[] sequences = new String[titles.length];
		for(int i = 0; i < titles.length; i++) {
			String sequence = titleToStringBuffer.get(titles[i]).toString();
			sequences[i] = sequence;
		}
		
		r = new SequenceAlignment(sequences, titles);

		return r;
	}
	
	public static SequenceAlignment parsePhylipAlignment(String alignment) {
		SequenceAlignment r = null;
		//remove the first line, containing the number of taxa and length
		//of alignment, then use parseClustalAlignment() to complete parsing
		
		int endOfFirstLine = alignment.indexOf('\n');
		String trimmedAlignment =
			alignment.substring(endOfFirstLine+1, alignment.length());
		r = parseClustalAlignment(trimmedAlignment);
		return r;
	}
	
	public static SequenceAlignment parseClustalAlignmentFile(String fileName) 
	       throws IOException {
		SequenceAlignment r = null;
		TextFile alignmentFile = new TextFile(fileName);
		String alignment = alignmentFile.toString();
		r = parseClustalAlignment(alignment);
		alignmentFile.closeFile();
		return r;
	}

	public static SequenceAlignment parseFastaAlignmentFile(String fileName) throws IOException {
		SequenceAlignment r = null;
		FastaSequenceFile alignmentFile = new FastaSequenceFile(fileName);
		r = parseFastaAlignment(alignmentFile);
		r.setName(fileName);
		alignmentFile.closeFile();
		return r;
	}

	public static SequenceAlignment parseFastaAlignment(FastaSequenceSet alignment) throws IOException {
		SequenceAlignment r = null;
		String[] titles = alignment.getTitles();
		String[] sequences = new String[titles.length];
		for(int i = 0; i < titles.length; i++) {
			sequences[i] = alignment.getPlainSequence(i);
		}
		
		String[] taxa = FastaUtilities.getTaxaFromTitles(titles);
		r = new SequenceAlignment(sequences, titles, taxa);
		return r;
	}	
}
