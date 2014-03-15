package edu.vt.vbi.ci.pepr.alignment;

public interface SWAligner {

	/**
	 * Retruns the score for the optimal alignment of 
	 * seq1 and seq2. This should be more efficient than 
	 * getting the actual alignment, because traceback
	 *  is not needed.
	 * 
	 * @param seq1
	 * @param seq2
	 * @return
	 */
	public abstract int getAlignmentScore(int[] seq1, int[] seq2);

	public abstract SequencePairAlignment getAlignment(int[] seq1, int[] seq2);

	public abstract int getBandedAlignmentScore(int[] seqI, int[] seqJ,
			int bandWidth);

}