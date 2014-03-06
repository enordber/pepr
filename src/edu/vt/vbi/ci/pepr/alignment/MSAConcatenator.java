package edu.vt.vbi.ci.pepr.alignment;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;

import edu.vt.vbi.ci.util.file.TextFile;

/**
 * Takes a set of input SequenceAlignments and generates
 * a concatenated SequenceAlignment.
 * 
 * @author enordber
 *
 */
public class MSAConcatenator {

	private static boolean debug = false;

	/**
	 * This main() method is for testing purposes only
	 */
	public static void main(String[] args) {
		//assume each argument is the name of an alignment file
		//in fasta format. Get SequenceAlignment for each file, 
		//the concatenate them all

		SequenceAlignment[] alignments = new SequenceAlignment[args.length];
		//		System.out.println("load alignment files...");
		try {
			for(int i = 0; i < alignments.length; i++) {
				//				System.out.print("\r                                     \r");
				//				System.out.println(i + ": " + args[i]);
				String phylipAlignment = new TextFile(args[i]).toString();
				alignments[i] = 
					SequenceAlignmentParser.parseFastaAlignmentFile(args[i]);
//					SequenceAlignmentParser.parsePhylipAlignment(phylipAlignment);
			}
		} catch(IOException ioe) {
			ioe.printStackTrace();
		}
		//		System.out.println();
		//		System.out.println("all alignments have been loaded.");
		//		System.out.println("concatenate..");
		SequenceAlignment concatenated = concatenate(alignments);
		//		System.out.println("concatenation complete");
		String[] alignmentTaxa = alignments[1].getTaxa();
		String[] concatenatedTaxa = concatenated.getTaxa();

		int from = 0;
		int to = 50;
//				System.out.println("alignments[1] taxa:");
//				for(int i = 0; i < alignmentTaxa.length;i ++) {
//					System.out.println(i + ": \t" 
//							+ alignments[1].getSequenceString(i).substring(from,to)
//							+ "\t" + alignmentTaxa[i] );
//				}

				System.out.println();
				System.out.println("concatenated taxa:");
				for(int i = 0; i < concatenatedTaxa.length; i++) {
					System.out.println(i + ": \t" 
							+ concatenated.getSequenceString(i).
							substring(alignments[0].getLength()+from,
									alignments[0].getLength()+to)
									+ "\t" + concatenatedTaxa[i]);
				}

//		String concatenatedContents = concatenated.getAlignmentAsExtendedPhylipUsingTaxonNames();
//		System.out.println(concatenatedContents);
	}

	/**
	 * @param alignments
	 * @return
	 */
	public static ConcatenatedSequenceAlignment concatenate(SequenceAlignment[] alignments) {
		if(debug) {
			System.out.println("MSAConcatenator.concatenate() alignments: " 
					+ alignments.length);
		}
		
		//remove any null entries in alignments
		//In some cases, alignments trimmed by Gblock are being 
		//trimmed down to nothing, resulting in a null alignment.
		ArrayList alignmentList = new ArrayList(alignments.length);
		for(int i = 0; i < alignments.length; i++) {
			if(alignments[i] != null) {
				alignmentList.add(alignments[i]);
			}
		}
		
		alignments = new SequenceAlignment[alignmentList.size()];
		alignmentList.toArray(alignments);
		
		ConcatenatedSequenceAlignment r = new ConcatenatedSequenceAlignment();
		r.setSequenceAlignments(alignments);
		//concatenate alignments on the basis of common taxa
		//determine the union of the taxa in the various alignments
		HashSet taxonUnionSet = new HashSet();
		for(int i = 0; i < alignments.length; i++) {
			String[] taxa = alignments[i].getTaxa();
			for(int j = 0; j < taxa.length; j++) {
				taxonUnionSet.add(taxa[j]);
				if(debug) {
					System.out.println("MSAConcatenator.concatenate() adding to taxonUnionSet: "
							+ taxa[j]);
				}
			}
		}

		//At some point, a null pointer exception was thrown by the following 
		//sort operation. Try removing null from the taxonUnionSet in case
		//it somehow managed to get in there.
		taxonUnionSet.remove(null);

		String[] taxonUnion = new String[taxonUnionSet.size()];
		taxonUnionSet.toArray(taxonUnion);
		Arrays.sort(taxonUnion);

		//determine length of the concatenated alignment
		int concatenatedLength = 0;
		for(int i = 0; i < alignments.length; i++) {
			concatenatedLength += alignments[i].getLength();
		}

		//for each taxon in the union, create a concatenated sequence
		for(int i = 0; i < taxonUnion.length; i++) {
			StringBuffer sb = new StringBuffer();
			for(int j = 0; j < alignments.length; j++) {
				int taxonIndex = alignments[j].getTaxonIndex(taxonUnion[i]);
				if(debug) {
					System.out.println("MSAConcatenator.concatenate() add sequence for taxon: "
							+ taxonUnion[i]);
				}
				if(taxonIndex >= 0) {
					String taxonSequence = 
						alignments[j].getSequenceString(taxonIndex);
					if(taxonSequence != null) {
						sb.append(taxonSequence);
					} else {
						System.out.println("Problem with alignment " + j + 
								" for taxon " + i + 
								" '" + taxonUnion[i] + "'. index in alignment is "
								+ taxonIndex + 
						" . here are the taxa: " );
						String[] alignmentTaxa = alignments[j].getTaxa();
						for(int k = 0; k < alignmentTaxa.length; k++) {
							System.out.println("" + k + ": " 
									+ alignmentTaxa[k]);
						}
						System.out.println("here are the titles: ");
						String[] titles = alignments[j].getSequenceTitles();
						for(int k = 0; k < titles.length; k++) {
							System.out.println("" + k + ": " 
									+ titles[k]);
						}
						
						System.out.println("here is the alignment:");
						String phylipAlignment = alignments[j].getAlignmentAsExtendedPhylipUsingSequenceNames();
						System.out.println(phylipAlignment);
					}
				} else {
					//this taxon is not present in this alignment. Pad this 
					//region of the alignment with "missing" characters ('?')
					char[] gaps = new char[alignments[j].getLength()];
					Arrays.fill(gaps, SequenceAlignment.MISSING_CHAR);
					sb.append(gaps);
				}
			}

			//add each concatenated sequence to the concatenated alignment
			r.addSequence(sb.toString(), taxonUnion[i]);
		}
		r.setTaxa(taxonUnion);
		StringBuffer alignmentName =new StringBuffer("cat");
		for(int i = 0; i < alignments.length; i++) {
			alignmentName.append("_");
			alignmentName.append(alignments[i].getName());
		}
		int maximumNameLength = 64;
		String name = 
			alignmentName.substring(0, 
					Math.min(alignmentName.length(), maximumNameLength));

		r.setName(name);
		return r;
	}
}