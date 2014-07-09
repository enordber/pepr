package edu.vt.vbi.ci.pepr.alignment;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import edu.vt.vbi.ci.util.file.TextFile;

/**
 * Contains methods for sequence conversions and
 * scoring matrices.
 * @author enordber
 *
 */
public class AlignmentUtilities {
	
//	 one source says "lambda was estimated at 0.252 and K at 0.035"
	//for BLOUM62 matrix
	//To calculate a bit score from a raw score :
	//S' = (lambda x S - lnK) / ln2
	//To calculate E value from the bit score:
	//E = mn2^-S'
	
	//The S' needed to get a certain E value is:
	//S' = 1/log2(E/mn)
	public static final String BLOSUM_62 = "BLOSUM62";
	public static final String AMINO_ACID = "Amino Acid";
	public static final String NUCLEIC_ACID = "Nucleic Acid";
	// A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X

	//Amino Acid constants
	public static final int A = 0;
	public static final int R = 1;
	public static final int N = 2;
	public static final int D = 3;
	public static final int C = 4;
	public static final int Q = 5;
	public static final int E = 6;
	public static final int G = 7;
	public static final int H = 8;
	public static final int I = 9;
	public static final int L = 10;
	public static final int K = 11;
	public static final int M = 12;
	public static final int F = 13;
	public static final int P = 14;
	public static final int S = 15;
	public static final int T = 16;
	public static final int W = 17;
	public static final int Y = 18;
	public static final int V = 19;
	public static final int B = 20;
	public static final int Z = 21;
	public static final int X = 22;
	public static final int GAP = 23;

	//Amino Acid bit mask constants
	public static final int A_MASK = 0;
	public static final int R_MASK = 1 << 0;
	public static final int N_MASK = 1 << 1;
	public static final int D_MASK = 1 << 2;
	public static final int C_MASK = 1 << 3;
	public static final int Q_MASK = 1 << 4;
	public static final int E_MASK = 1 << 5;
	public static final int G_MASK = 1 << 6;
	public static final int H_MASK = 1 << 7;
	public static final int I_MASK = 1 << 8;
	public static final int L_MASK = 1 << 9;
	public static final int K_MASK = 1 << 10;
	public static final int M_MASK = 1 << 11;
	public static final int F_MASK = 1 << 12;
	public static final int P_MASK = 1 << 13;
	public static final int S_MASK = 1 << 14;
	public static final int T_MASK = 1 << 15;
	public static final int W_MASK = 1 << 16;
	public static final int Y_MASK = 1 << 17;
	public static final int V_MASK = 1 << 18;
	public static final int B_MASK = 1 << 19;
	public static final int Z_MASK = 1 << 20;
	public static final int X_MASK = 1 << 21;
	public static final int GAP_MASK = 1 << 22;

	private static final char[] aaTranslate = 
		new char[]{'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K',
		'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X'};
	
	/*
	 *Nucleotide constants 
	 */
	public static final int naA = 0;
	public static final int naT = 1;
	public static final int naG = 2;
	public static final int naC = 3;
	
	private static final char[] naTranslate = 
		new char[]{'A', 'T', 'G', 'C'};
	
	private static final char GAP_CHAR = '-';
	
	/**
	 * Returns an array for translating from int[] sequence
	 * to char sequence. Use the int in the sequence array
	 * as index to get the corresponding char.
	 * 
	 * @return
	 */
	public static char[] getNucleotideTranslation() {
		char[] r = new char[naTranslate.length];
		System.arraycopy(naTranslate, 0, r, 0, r.length);
		return r;
	}

	public static int[] convertNucleotideSequenceToInts(String sequence) {
		int[] r = null;
		int rIndex = 0;

		byte[] seqChars = sequence.getBytes();
		r = new int[seqChars.length];
		for(int i = 0; i < seqChars.length; i++) {
			switch(seqChars[i]) {
			case 'A':
			case 'a':
				r[rIndex] = naA;
				rIndex++;
				break;
			case 'T':
			case 't':
			case 'U':
			case 'u':
				r[rIndex] = naT;
				rIndex++;
				break;
			case 'G':
			case 'g':
				r[rIndex] = naG;
				rIndex++;
				break;
			case 'C':
			case 'c':
				r[rIndex] = naC;
				rIndex++;
				break;
			}
		}
		if(rIndex < r.length) {
			int[] compress = new int[rIndex];
			System.arraycopy(r, 0, compress, 0, rIndex);
			r = compress;
		}
		return r;
	}
	
	/**
	 * Returns an array for translating from int[] sequence
	 * to char sequence. Use the int in the sequence array
	 * as index to get the corresponding char.
	 * 
	 * @return
	 */
	public static char[] getAminoAcidTranslation() {
		char[] r = new char[aaTranslate.length];
		System.arraycopy(aaTranslate, 0, r, 0, r.length);
		return r;
	}
	
	/**
	 * Converts the given amino acid sequence to an int[].
	 * Sequence String case does not matter. Any non-sequence
	 * characters will be ignored. 
	 * 
	 * @param sequence
	 * @return
	 */
	public static int[] convertAminaAcidSequenceToInts(String sequence) {
		int[] r = null;
		int rIndex = 0;
		
		byte[] seqChars = sequence.getBytes();
		r = new int[seqChars.length];
		for(int i = 0; i < seqChars.length; i++) {
			switch(seqChars[i]) {
			case 'A' :
			case 'a' :
				r[rIndex] = A;
				rIndex++;
				break;
				
			case 'B' :
			case 'b' :
				r[rIndex] = B;
				rIndex++;
				break;
				
			case 'C' :
			case 'c' :
				r[rIndex] = C;
				rIndex++;
				break;
				
			case 'D' :
			case 'd' :
				r[rIndex] = D;
				rIndex++;
				break;
				
			case 'E' :
			case 'e' :
				r[rIndex] = E;
				rIndex++;
				break;
				
			case 'F' :
			case 'f' :
				r[rIndex] = F;
				rIndex++;
				break;
				
			case 'G' :
			case 'g' :
				r[rIndex] = G;
				rIndex++;
				break;
				
			case 'H' :
			case 'h' :
				r[rIndex] = H;
				rIndex++;
				break;
				
			case 'I' :
			case 'i' :
				r[rIndex] = I;
				rIndex++;
				break;
				
			case 'K' :
			case 'k' :
				r[rIndex] = K;
				rIndex++;
				break;
				
			case 'L' :
			case 'l' :
				r[rIndex] = L;
				rIndex++;
				break;
				
			case 'M' :
			case 'm' :
				r[rIndex] = M;
				rIndex++;
				break;
				
			case 'N' :
			case 'n' :
				r[rIndex] = N;
				rIndex++;
				break;
				
			case 'P' :
			case 'p' :
				r[rIndex] = P;
				rIndex++;
				break;
				
			case 'Q' :
			case 'q' :
				r[rIndex] = Q;
				rIndex++;
				break;
				
			case 'R' :
			case 'r' :
				r[rIndex] = R;
				rIndex++;
				break;
				
			case 'S' :
			case 's' :
				r[rIndex] = S;
				rIndex++;
				break;
				
			case 'T' :
			case 't' :
				r[rIndex] = T;
				rIndex++;
				break;
				
			case 'V' :
			case 'v' :
				r[rIndex] = V;
				rIndex++;
				break;
				
			case 'W' :
			case 'w' :
				r[rIndex] = W;
				rIndex++;
				break;
				
			case 'X' :
			case 'x' :
				r[rIndex] = X;
				rIndex++;
				break;
				
			case 'Y' :
			case 'y' :
				r[rIndex] = Y;
				rIndex++;
				break;
				
			case 'Z' :
			case 'z' :
				r[rIndex] = Z;
				rIndex++;
				break;
				
			case '-' :
				r[rIndex] = GAP;
				rIndex++;
				break;

				//There is not special int value for missing right now. 
				//The only place I am using this, I treat missing the same
				//as a gap. This can change if I ever treat missing differently
			case '?' :
				r[rIndex] = GAP;
				rIndex++;
				break;
}
		}
		
		if(rIndex < r.length) {
			int[] compress = new int[rIndex];
			System.arraycopy(r, 0, compress, 0, rIndex);
			r = compress;
		}
		
		return r;
	}
	
	public static int[][] loadScoringMatrixFromFile(String fileName) {
		int[][] r = null;
		try {
			TextFile file = new TextFile(fileName);
			
			//expect a 23x23 matrix of integers
			r = new int[24][24];
			
			String[] fileLines = file.getAllLines();
			int firstMatrixLine = 0;
			for(; !fileLines[firstMatrixLine].startsWith("A");firstMatrixLine++) {}
			
			for(int i = 0; i < r.length-1; i++) {
				String matrixLine = fileLines[firstMatrixLine+i];
				String[] cells = matrixLine.split("\\s+");
				for(int j = 0; j < r[i].length-1; j++) {
					int nextInt = Integer.parseInt(cells[j+1]);
					r[i][j] = nextInt;
				}
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
//			e.printStackTrace();
		}
		return r;
	}
	
	public static int[][] loadScoringMatrixFromResource(String resourceName) {
		int[][] r = null;
		try {
			InputStream is = ClassLoader.getSystemClassLoader().getResourceAsStream(resourceName);
			InputStreamReader isr = new InputStreamReader(is);
			BufferedReader br = new BufferedReader(isr);
			
			//expect a 23x23 matrix of integers
			r = new int[24][24];
			
			String matrixLine;
			while(!(matrixLine = br.readLine()).startsWith("A")) {}
			
			int i = 0;
			do {
				String[] cells = matrixLine.split("\\s+");
				for(int j = 0; j < r[i].length-1; j++) {
					int nextInt = Integer.parseInt(cells[j+1]);
					r[i][j] = nextInt;
				}
				i++;
			} while((matrixLine = br.readLine()) != null && i < r.length - 1);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return r;
	}
	
    /**
     * one source says "lambda was estimated at 0.252 and K at 0.035"
     * @param rawScore
     * @param lambda
     * @param K
     * @return
  	//for BLOUM62 matrix
	//To calculate a bit score from a raw score :
	//S' = (lambda x S - lnK) / ln2
	
	
	//The S' needed to get a certain E value is:
	//S' = 1/log2(E/mn)
	 * **/
	public static double convertRawScoreToBitScore(double rawScore, 
			double lambda, double K) {
		double r = Double.NaN;
		r = (lambda * rawScore - Math.log(K)) / Math.log(2);
		return r;
	}

	/**
	 * To calculate E value from the bit score:
	 * E = mn2^-S'
	 * 
	 * @return
	 */
	public static double getEValue(int m, int n, double bitScore) {
		double r = Double.NaN;
		bitScore = 0-bitScore;
		r = m * n * Math.pow(2, bitScore);
		return r;
	}
	
	/**
	 * Takes an aligned amino acid sequence and maps the gaps to an unaligned
	 * nucleotide sequence. This is used to convert from aligned amino acids
	 * to aligned nucleotides, keeping the codons intact.
	 * 
	 * This method works on one sequence at a time, so converting a multiple 
	 * sequence alignment will require multiple calls - one for each sequence
	 * in the alignment.
	 * 
	 * @param alignedAASeq
	 * @param unalignedNTSeq
	 * @return
	 */
	public static String mapAlignmentGapsToNTSeq(String alignedAASeq, String unalignedNTSeq) {
		String r= null;
		char[] unalignedNTChars = unalignedNTSeq.toCharArray();
		int ntIndex = 0;
		StringBuffer alignedNTBuffer = new StringBuffer(alignedAASeq.length()*3);
		char[] alignedAAChars = alignedAASeq.toCharArray();
		
		for(int i = 0; i < alignedAAChars.length; i++) {
			if(alignedAAChars[i] == GAP_CHAR) {
				//append three gaps, taking the place of the single
				//gap in the amino acid alignment
				alignedNTBuffer.append("---");
			} else {
				//append the next three nucleotides, which should be 
				//the next codon
				alignedNTBuffer.append(unalignedNTChars[ntIndex++]);
				alignedNTBuffer.append(unalignedNTChars[ntIndex++]);
				alignedNTBuffer.append(unalignedNTChars[ntIndex++]);
			}
		}
		r = alignedNTBuffer.toString();
		return r;
	}
}
