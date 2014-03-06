package edu.vt.vbi.ci.pepr.alignment;

public class SequencePairAlignment {

	private char gap = '-';
	private char space = ' ';
	private char pipe = '|';
	private char[] seq1;
	private char[] seq2;
	private int[] seq1Align;
	private int[] seq2Align;
	private int score;
	private String resultString;
	
	/*
	 * Alignment start and end values use lazy calculation.
	 * They will be determined the first time any one 
	 * of them is requested.
	 */
	private int seq1AlignStart = -1;
	private int seq1AlignEnd = -1;
	private int seq2AlignStart = -1;
	private int seq2AlignEnd = -1;
	
	public SequencePairAlignment(char[] seq1, char[] seq2, 
			int[] seq1Align, int[] seq2Align, int score) {
		this.seq1 = seq1;
		this.seq2 = seq2;
		this.seq1Align = seq1Align;
		this.seq2Align = seq2Align;
		this.score = score;
		
		//calling the following methods for debugging purposes
		//only. This can be removed after debugging
		determineAlignmentStartsAndEnds();
		resultString = toString();
	}
	
	public int getScore() {
		return score;
	}
	
	/**
	 * Returns an array of three Strings.
	 * 0: Sequence 1 in the alignment
	 * 1: pipes or spaces
	 * 2: Sequence 2 in the alignment
	 * @return
	 */
	public String[] getAlignmentStrings() {
		String[] r = new String[3];
		//the String for seq1 is actuallly based on the content in 
		//seq2Align, and vise-versa
		StringBuffer seq1Result = new StringBuffer();
		StringBuffer seq2Result = new StringBuffer();
		StringBuffer pipes = new StringBuffer();
		char atSeq1; //char most recently added to seq1Result
		char atSeq2; //char most recently added to seq2Result

		//start and end positions for local alignment
		int seq1Start = getSeq1AlignStart();
		int seq1End = getSeq1AlignEnd();
		int seq2Start = getSeq2AlignStart();
		int seq2End = getSeq2AlignEnd();

		int seq1Index = seq1Start; //index in seq1
		int seq2Index = seq2Start; //index in seq2
		int seq1AIndex = seq1Start; //index in seq1Align
		int seq2AIndex = seq2Start; //index in seq2Align
		
		
		//add leading portion of seq2, if any exists
		if(seq1Align[seq1Start] > 0) {
			for(int i = seq2Start; i < seq1Align[seq1Start]; i++) {
				seq2Result.append(seq2[i]);
			}
		}
		//add leading portion of seq1, if any exists
		if(seq2Align[seq2Start] > 0) {
			for(int i = seq1Start; i < seq2Align[seq2Start]; i++) {
				seq1Result.append(seq1[i]);
			}
		}
		while(seq1AIndex <= seq1End && seq2AIndex <= seq2End) {
			if(seq1Align[seq1AIndex] < 0) {
				//this position of sequence 1 is opposite a gap,
				//so enter a gap in sequence 2
				atSeq2 = gap;
			} else {
				seq2Index = seq1Align[seq1AIndex];
				atSeq2 = seq2[seq2Index];
			}
			seq2Result.append(atSeq2);
			
			if(seq2Align[seq2AIndex] < 0) {
				//this position of sequence 2 is opposite a gap,
				//so enter a gap in sequence 1
				atSeq1 = gap;
			} else {
				seq1Index = seq2Align[seq2AIndex];
				atSeq1 = seq1[seq1Index];
			}
			seq1Result.append(atSeq1);

			seq1AIndex++;
			seq2AIndex++;
		}
		

		//add any remaining elements in seq1, 
		//and corresponding gaps in seq2
		while(++seq1Index < seq1End) {
			seq1Result.append(seq1[seq1Index]);
			if(seq1Align[seq1Index] < 0) {
				seq2Result.append(gap);
			}
		}
		
		//add any remaining elements in seq2,
		//and corresponding gaps in seq1
		while(++seq2Index < seq2End) {
			seq2Result.append(seq2[seq2Index]);
			if(seq2Align[seq2Index] < 0) {
				seq1Result.append(gap);
			}
		}
		
		//determine pipes - this could be made more efficient
		int length = Math.min(seq1Result.length(), seq2Result.length());
		for(int i = 0; i < length; i++) {
			if(seq1Result.charAt(i) == seq2Result.charAt(i)) {
				if(seq1Result.charAt(i) == gap) {
					//alignment has a gap in both sequences at this position,
					//so just remove both. This happend because of the
					//way a tie is broken in SW alignment. This should not
					//indicate any problem with the actual alignment

					seq1Result.deleteCharAt(i);
					seq2Result.deleteCharAt(i);
					length--;
					i--;
				} else {
					pipes.append(pipe);
				}
			} else {
				pipes.append(space);
			}
		}
		 
		r[0] = seq1Result.toString();
		r[1] = new String(pipes);
		r[2] = seq2Result.toString();
		return r;
	}
	
	private void determineAlignmentStartsAndEnds() {
		boolean start2Found = false;
		for(int i = 0; i < seq1Align.length; i++) {
			if(seq1Align[i] > -1) {
				seq2AlignEnd = seq1Align[i];
				if(!start2Found) {
					seq2AlignStart = seq1Align[i];
					start2Found = true;
				}
			}
		}
		
		boolean start1Found = false;
		for(int i = 0; i < seq2Align.length; i++) {
			if(seq2Align[i] > -1) {
				seq1AlignEnd = seq2Align[i];
				if(!start1Found) {
					seq1AlignStart = seq2Align[i];
					start1Found = true;
				}
			}
		}
	}

	public int getSeq1AlignEnd() {
		if(seq1AlignEnd == -1) {
			determineAlignmentStartsAndEnds();
		}
		return seq1AlignEnd;
	}

	public int getSeq1AlignStart() {
		if(seq1AlignStart == -1) {
			determineAlignmentStartsAndEnds();
		}
		return seq1AlignStart;
	}

	public int getSeq2AlignEnd() {
		if(seq2AlignEnd == -1) {
			determineAlignmentStartsAndEnds();
		}
		return seq2AlignEnd;
	}

	public int getSeq2AlignStart() {
		if(seq2AlignStart == -1) {
			determineAlignmentStartsAndEnds();
		}
		return seq2AlignStart;
	}
	
	public String toString() {
		StringBuffer sb = new StringBuffer();
		sb.append("score: ");
		sb.append(getScore());
		sb.append("\t");
		sb.append("seq1 ");
		sb.append(getSeq1AlignStart());
		sb.append(" to ");
		sb.append(getSeq1AlignEnd());
		sb.append(" of ");
		sb.append(seq1.length);
		sb.append("\t");
		sb.append("seq2 ");
		sb.append(getSeq2AlignStart());
		sb.append(" to ");
		sb.append(getSeq2AlignEnd());
		sb.append(" of ");
		sb.append(seq2.length);
		sb.append("\n");
		String[] alignLines = getAlignmentStrings();
		for(int i = 0; i < alignLines.length; i++) {
			sb.append(alignLines[i]);
			sb.append("\n");
		}
		return new String(sb);
	}
	
}
