package edu.vt.vbi.ci.pepr.alignment;

import java.util.Arrays;

/**
 * This class performs Smith-Waterman sequence alignments.
 * This is intended for high-throughput applications
 * in which many alignments will be performed using the
 * same scoring matrix and (generally) the same gap 
 * penalties. 
 * @author enordber
 *
 */
public class SWProteinAligner implements SWAligner {
	public static final String BLOSUM_62 = "BLOSUM62";
	public static final String BLOSUM_100 = "BLOSUM_100";

	private static String blosum62FileName =
		"/Users/enordber/vbi/workspace/cigswg/BLOSUM62";
	private static String blosum100FileName =
		"/Users/enordber/vbi/workspace/cigswg/BLOSUM100";

	/**
	 * Cost of initiating a gap.
	 */
	private int gapOpenPenalty;

	/**
	 * Cost of extending a gap by one.
	 */
	private int gapExtensionPenalty;

	/**
	 * Scores for aligning amino acids.
	 */
	private int[][] scoringMatrix;
	
	/**
	 * Maps from int values to characters for amino acid sequences.
	 * Using lazy initialization.
	 */
	private static char[] translationTable;

	public SWProteinAligner(int gapOpenPenalty, 
			int gapExtensionPenalty, int[][] scoringMatrix) {
		this.gapOpenPenalty = Math.abs(gapOpenPenalty);
		this.gapExtensionPenalty = Math.abs(gapExtensionPenalty);
		this.scoringMatrix = scoringMatrix;
	}

	public SWProteinAligner(int gapOpenPenalty,
			int gapExtensionPenalty, String matrixName) {
		this.gapOpenPenalty = Math.abs(gapOpenPenalty);
		this.gapExtensionPenalty = Math.abs(gapExtensionPenalty);

		loadMatrix(matrixName);
	}

	private void loadMatrix(String matrixName) {
		if(matrixName.equals(BLOSUM_100)) {
			scoringMatrix = 
				AlignmentUtilities.loadScoringMatrixFromFile(blosum100FileName);
		} else if(matrixName.equals(BLOSUM_62)) {
			try {
				scoringMatrix = 
					AlignmentUtilities.loadScoringMatrixFromFile(
							blosum62FileName);
			} catch(Exception e) {
			}
			if(scoringMatrix == null) {
				scoringMatrix =
					AlignmentUtilities.loadScoringMatrixFromResource(
					"BLOSUM62");

			}
		}
	}

	/* (non-Javadoc)
	 * @see alignment.SWAligner#getAlignmentScore(int[], int[])
	 */
	public int getAlignmentScore(int[] seq1, int[] seq2) {
		int r = 0;

		//use index i for seq1 and index j for seq2

		/*
		 * The various row arrays are all 1 element longer than
		 * seq2. This is to allow for the first row and first column
		 * in the alignment matrix to be filled with zeros. These 
		 * represent a leading gap in either of the sequences, which
		 * does not incur a penalty in Smith-Waterman alignment.
		 * A result of including the extra row and column is that
		 * the index values in these tracking arrays are off by one
		 * compared with seq1 and seq2. This can make for some 
		 * confusion in the loop that calculates the score. Element j
		 * in seq2 will correspond to element j+1 in the tracking 
		 * arrays.
		 */
		int trackingArrayLength = seq2.length+1;

		/*
		 * previousRowScores tracks the score values for cells in the i-1
		 * row.
		 */
		int[] previousRowScores = new int[trackingArrayLength];

		/*
		 * currentRowScores tracks the scores on the ith row as
		 * they are being calculated. This becomes the previousRowScores
		 * when moving to the next row.  
		 */
		int[] currentRowScores = new int[trackingArrayLength];

		/*
		 * previousRowGaps track the occurance of gaps in the i-1
		 * row. This is used to determine if a new gap will be
		 * a gap open or a gap extension.
		 */
		boolean[] previousRowGaps = new boolean[trackingArrayLength];

		/*
		 * currentRowGaps tracks the gaps introduced in the 
		 * ith row. This becomes previousRowGaps when moving
		 * to the next row.
		 */
		boolean[] currentRowGaps = new boolean[trackingArrayLength];

		int matchScore; //align i, j
		int scoreFromPreviousRow; //insert/extend gap in seq1
		int scoreFromPreviousCol; //insert/extend gap in seq2

		//index in tracking arrays is offset by one (see above)
		int trackingIndex;

		//used to hold index for scoring matrix element in the flatScoringMatrix
		int matrixIndex;
		int matrixLength = scoringMatrix.length;
		int seq1Length = seq1.length;
		int seq2Length = seq2.length;
		int[] swapI;
		boolean[] swapB;

		for(int i = 0; i < seq1Length; i++){
			for(int j = 0; j < seq2Length; j++) {
//				System.out.println("check " + i + "," + j +
//						". seq elements are " + seq1[i] + "," + seq2[j] +
//						": score for match is " + scoringMatrix[seq1[i]][seq2[j]]);
				trackingIndex = j+1;

				//get value for i,j match 
				//matchScore is going diagonally in the DP matrix 
				matrixIndex = seq1[i]*matrixLength + seq2[j];
				matchScore = previousRowScores[j] + scoringMatrix[seq1[i]][seq2[j]];

				//scoreFromPreviousRow is going vertically in the DP matrix
				if(previousRowGaps[trackingIndex]) {
					//extend gap by 1
					scoreFromPreviousRow = previousRowScores[trackingIndex] - gapExtensionPenalty;
				} else {
					//open a gap
					scoreFromPreviousRow = previousRowScores[trackingIndex] - gapOpenPenalty;
				}
				
				//scoreFromPreviousCol is going horizontally in the DP matix
				if(currentRowGaps[j]) { 
					//extend gap by 1
					scoreFromPreviousCol = currentRowScores[j] - gapExtensionPenalty;
				} else {
					//open a gap
					scoreFromPreviousCol = currentRowScores[j] - gapOpenPenalty;
				}

				//determine which of the scores is the maximum
				currentRowScores[trackingIndex] = scoreFromPreviousRow;
				currentRowGaps[trackingIndex] = true;
				if(matchScore > currentRowScores[trackingIndex]) {
					currentRowScores[trackingIndex] = matchScore;
					currentRowGaps[trackingIndex] = false;
				}
				if(scoreFromPreviousCol > currentRowScores[trackingIndex]) {
					currentRowScores[trackingIndex] = scoreFromPreviousCol;
					currentRowGaps[trackingIndex] = true;
				}

				//SW alignment imposes a minimum score of 0
				if(currentRowScores[trackingIndex] < 0) {
					currentRowScores[trackingIndex] = 0;
				} else if(currentRowScores[trackingIndex] > r) {
					//SW if local alignment, so the maximum score 
					//anywhere in the matrix is what we want.
					r = currentRowScores[trackingIndex];
				}

				//System.out.println("value for this cell: " + 
				//		currentRowScores[trackingIndex]);
			}

			//prepare tracking arrays for next row, if there is a next row
			if(i+1 < seq1.length) {
				swapI = previousRowScores;
				swapB = previousRowGaps;
				previousRowScores = currentRowScores;
				previousRowGaps = currentRowGaps;
				currentRowScores = swapI;
				currentRowGaps = swapB;
			}
		}
		return r;
	}

	/* (non-Javadoc)
	 * @see alignment.SWAligner#getAlignment(int[], int[])
	 */
	public SequencePairAlignment getAlignment(int[] seq1, int[] seq2) {
		/*
		 * Tracks the highest score encountered
		 */
		int highScore = -1;

		/*
		 * The various row arrays are all 1 element longer than
		 * seq2. This is to allow for the first row and first column
		 * in the alignment matrix to be filled with zeros. These 
		 * represent a leading gap in either of the sequences, which
		 * does not incur a penalty in Smith-Waterman alignment.
		 * A result of including the extra row and column is that
		 * the index values in these tracking arrays are off by one
		 * compared with seq1 and seq2. This can make for some 
		 * confusion in the loop that calculates the score. Element j
		 * in seq2 will correspond to element j+1 in the tracking 
		 * arrays.
		 */

		//use index i for seq1 and index j for seq2

		/*
		 * index in seq1(+1) where highest scoring alignment ends
		 */
		int highI = -1;
		/*
		 * index in seq2(+1) where highest scoring alignment ends
		 */
		int highJ = -1;


		/*
		 * Dynamic programming matrix that stores all of the scores.
		 */
		int[][] dpMatrix = new int[seq1.length+1][seq2.length+1];

		/*
		 * Bactracking info
		 */
		final byte btVert = 1;
		final byte btDiag = 0;
		final byte btHor  = 2;
		byte[][] backTrackMatrix = new byte[seq1.length+1][seq2.length+1];

		/*
		 * previousRowScores tracks the score values for cells in the i-1
		 * row.
		 */
		int[] previousRowScores = dpMatrix[0];

		/*
		 * currentRowScores tracks the scores on the ith row as
		 * they are being calculated. This becomes the previousRowScores
		 * when moving to the next row.  
		 */
		int[] currentRowScores = dpMatrix[1];

		/*
		 * previousRowGaps track the occurance of gaps in the i-1
		 * row. This is used to determine if a new gap will be
		 * a gap open or a gap extension.
		 */
		byte[] previousRowBackTrack = backTrackMatrix[0];

		/*
		 * currentRowGaps tracks the gaps introduced in the 
		 * ith row. This becomes previousRowGaps when moving
		 * to the next row.
		 */
		byte[] currentRowBackTrack = backTrackMatrix[1];

		int matchScore; //align i, j
		int scoreFromPreviousRow; //insert/extend gap in seq1
		int scoreFromPreviousCol; //insert/extend gap in seq2

		//index in tracking arrays are offset by one (see above)
		int jTrackingIndex;
		int iTrackingIndex;

		int seq1i;
		byte btForCell;
		int bestForCell;

		for(int i = 0; i < seq1.length; i++){
			iTrackingIndex = i+1;
			seq1i = seq1[i];
			for(int j = 0; j < seq2.length; j++) {
				jTrackingIndex = j+1;
				//get score for i,j match in this cell
				matchScore = previousRowScores[j] + scoringMatrix[seq1i][seq2[j]];
				bestForCell = matchScore;
				btForCell = btDiag;

				//get score for vertical gap in this cell
				if(previousRowBackTrack[jTrackingIndex] == btVert) {
					//extend gap by 1
					scoreFromPreviousRow = previousRowScores[jTrackingIndex] - gapExtensionPenalty;
				} else {
					//open a gap
					scoreFromPreviousRow = previousRowScores[jTrackingIndex] - gapOpenPenalty;
				}
				if(scoreFromPreviousRow > bestForCell) {
					bestForCell = scoreFromPreviousRow;
					btForCell = btVert;
				}

				//get score for horizontal gap in this cell
				if(currentRowBackTrack[j] == btHor) { 
					//extend gap by 1
					scoreFromPreviousCol = currentRowScores[j] - gapExtensionPenalty;
				} else {
					//open a gap
					scoreFromPreviousCol = currentRowScores[j] - gapOpenPenalty;
				}
				if(scoreFromPreviousCol > bestForCell) {
					bestForCell = scoreFromPreviousCol;
					btForCell = btHor;
				}

				//SW alignment imposes a minimum score of 0
				if(bestForCell < 0) {
					bestForCell = 0;
				} else if(bestForCell > highScore) {
					//SW if local alignment, so the maximum score 
					//anywhere in the matrix is what we want.
					highScore = bestForCell;
					highI = iTrackingIndex;
					highJ = jTrackingIndex;
				}
				//set the score and backtrack for this cell
				currentRowScores[jTrackingIndex] = bestForCell;
				currentRowBackTrack[jTrackingIndex] = btForCell;				
			}

			//prepare tracking arrays for next row, if there is a next row
			if(iTrackingIndex < seq1.length) {
				previousRowScores = dpMatrix[iTrackingIndex];
				previousRowBackTrack = backTrackMatrix[iTrackingIndex];
				currentRowScores = dpMatrix[iTrackingIndex+1];
				currentRowBackTrack = backTrackMatrix[iTrackingIndex+1];
			}
		}

		//use backtracking matrix to build alignment
		//Alignment Strings will be build backwards, then reversed at 
		//the end. This is just easier, because I can use StringBuffer.append()
		//this way.
		//StringBuffer align1 = new StringBuffer();
		//StringBuffer align2 = new StringBuffer();
		int seq1Index = highI;
		int seq2Index = highJ;

		/*
		 *The alignment is stored in the following two arrays,
		 *seq1Align and seq2Align. seq1Align[i] contains the
		 *index in seq2 with which seq1[i] is matched. If 
		 *an element is not matched (opposite a gap) then
		 *-1 is used for the value. For this reason, the arrays
		 *are initially filled with -1, so leading and lagging
		 *portions (outside of the local alignment region) will
		 *be indicated as opposite gaps, rather than opposite
		 *element 0 of the other sequence. 
		 */
		int[] seq1Align = new int[seq1.length];
		int[] seq2Align = new int[seq2.length];
		Arrays.fill(seq1Align, -1);
		Arrays.fill(seq2Align, -1);
		//remember that all tracking index values are 1 more than the 
		//sequence index values (to allow for initial gaps in the alignment)
		//char gap = '-';
		int nextScore = dpMatrix[seq1Index][seq2Index];
		while(nextScore > 0) {
			switch(backTrackMatrix[seq1Index][seq2Index]) {
			case btHor:
				//align1.append(gap);
				//align2.append(translationTable[seq2[seq2Index-1]]);
				seq2Index--;	
				//seq2Align[seq2Index] = -1;
				break;
			case btDiag:
				//align1.append(translationTable[seq1[seq1Index-1]]);
				seq1Index--;
				//align2.append(translationTable[seq2[seq2Index-1]]);
				seq2Index--;
				seq1Align[seq1Index] = seq2Index;
				seq2Align[seq2Index] = seq1Index;
				break;
			case btVert:
				//align1.append(translationTable[seq1[seq1Index-1]]);
				seq1Index--;
				//align2.append(gap);
				break;
			}
			nextScore = dpMatrix[seq1Index][seq2Index];

		}

		//Create result Object
		if(translationTable == null) {
			translationTable = AlignmentUtilities.getAminoAcidTranslation();
		}
		char[] seq1Chars = new char[seq1.length];
		char[] seq2Chars = new char[seq2.length];
		for(int i = 0; i < seq1.length; i++) {
			seq1Chars[i] = translationTable[seq1[i]];
		}

		for(int i = 0; i < seq2.length; i++) {
			seq2Chars[i] = translationTable[seq2[i]];
		}
		SequencePairAlignment result = new SequencePairAlignment(seq1Chars,
				seq2Chars, seq1Align, seq2Align, highScore);
		return result;

	}

	/**
	 * Calculates the maximum score for an alignment between seq1 and seq2
	 * within a band of the specified length. There may be a better scoring
	 * alignment that extends outside of this band.
	 * 
	 * @param seq1
	 * @param seq2
	 * @param bandWidth
	 * @return
	 */
	public int getBandedAlignmentScore(int[] seq1, int[] seq2, int bandWidth) {
		int r = 0;

		//use index i for seq1 and index j for seq2

		/*
		 * The various row arrays are all 1 element longer than
		 * seq2. This is to allow for the first row and first column
		 * in the alignment matrix to be filled with zeros. These 
		 * represent a leading gap in either of the sequences, which
		 * does not incur a penalty in Smith-Waterman alignment.
		 * A result of including the extra row and column is that
		 * the index values in these tracking arrays are off by one
		 * compared with seq1 and seq2. This can make for some 
		 * confusion in the loop that calculates the score. Element j
		 * in seq2 will correspond to element j+1 in the tracking 
		 * arrays.
		 */
		int trackingArrayLength = seq2.length+1;

		/*
		 * previousRowScores tracks the score values for cells in the i-1
		 * row.
		 */
		int[] previousRowScores = new int[trackingArrayLength];

		/*
		 * currentRowScores tracks the scores on the ith row as
		 * they are being calculated. This becomes the previousRowScores
		 * when moving to the next row.  
		 */
		int[] currentRowScores = new int[trackingArrayLength];

		/*
		 * previousRowGaps track the occurance of gaps in the i-1
		 * row. This is used to determine if a new gap will be
		 * a gap open or a gap extension.
		 */
		boolean[] previousRowGaps = new boolean[trackingArrayLength];

		/*
		 * currentRowGaps tracks the gaps introduced in the 
		 * ith row. This becomes previousRowGaps when moving
		 * to the next row.
		 */
		boolean[] currentRowGaps = new boolean[trackingArrayLength];

		int matchScore; //align i, j
		int scoreFromPreviousRow; //insert/extend gap in seq1
		int scoreFromPreviousCol; //insert/extend gap in seq2

		//index in tracking arrays is offset by one (see above)
		int trackingIndex;

		//used to hold index for scoring matrix element in the flatScoringMatrix
		int matrixIndex;
		int matrixLength = scoringMatrix.length;
		int seq1Length = seq1.length;
		int seq2Length = seq2.length;
		int[] swapI;
		boolean[] swapB;
		int halfBand = bandWidth/2;

		for(int i = 0; i < seq1Length; i++){
			int bandStart = i-halfBand;
			int bandEnd = i+halfBand;
			if(bandStart < 0) {
				bandStart = 0;
			}
			if(bandEnd > seq2Length) {
				bandEnd = seq2Length;
			}
			
			for(int j = bandStart; j < bandEnd; j++) {
//				for(int j = 0; j < seq2Length; j++) {
//				System.out.println("check " + i + "," + j +
//						". seq elements are " + seq1[i] + "," + seq2[j] +
//						": score for match is " + scoringMatrix[seq1[i]][seq2[j]]);
				trackingIndex = j+1;

				//get value for i,j match 
				//matchScore is going diagonally in the DP matrix 
				matrixIndex = seq1[i]*matrixLength + seq2[j];
				matchScore = previousRowScores[j] + scoringMatrix[seq1[i]][seq2[j]];

				//scoreFromPreviousRow is going vertically in the DP matrix
				if(previousRowGaps[trackingIndex]) {
					//extend gap by 1
					scoreFromPreviousRow = previousRowScores[trackingIndex] - gapExtensionPenalty;
				} else {
					//open a gap
					scoreFromPreviousRow = previousRowScores[trackingIndex] - gapOpenPenalty;
				}
				
				//scoreFromPreviousCol is going horizontally in the DP matix
				if(currentRowGaps[j]) { 
					//extend gap by 1
					scoreFromPreviousCol = currentRowScores[j] - gapExtensionPenalty;
				} else {
					//open a gap
					scoreFromPreviousCol = currentRowScores[j] - gapOpenPenalty;
				}

				//determine which of the scores is the maximum
				currentRowScores[trackingIndex] = scoreFromPreviousRow;
				currentRowGaps[trackingIndex] = true;
				if(matchScore > currentRowScores[trackingIndex]) {
					currentRowScores[trackingIndex] = matchScore;
					currentRowGaps[trackingIndex] = false;
				}
				if(scoreFromPreviousCol > currentRowScores[trackingIndex]) {
					currentRowScores[trackingIndex] = scoreFromPreviousCol;
					currentRowGaps[trackingIndex] = true;
				}

				//SW alignment imposes a minimum score of 0
				if(currentRowScores[trackingIndex] < 0) {
					currentRowScores[trackingIndex] = 0;
				} else if(currentRowScores[trackingIndex] > r) {
					//SW if local alignment, so the maximum score 
					//anywhere in the matrix is what we want.
					r = currentRowScores[trackingIndex];
				}

				//System.out.println("value for this cell: " + 
				//		currentRowScores[trackingIndex]);
			}

			//prepare tracking arrays for next row, if there is a next row
			if(i+1 < seq1.length) {
				swapI = previousRowScores;
				swapB = previousRowGaps;
				previousRowScores = currentRowScores;
				previousRowGaps = currentRowGaps;
				currentRowScores = swapI;
				currentRowGaps = swapB;
			}
		}
		return r;
	}
}
