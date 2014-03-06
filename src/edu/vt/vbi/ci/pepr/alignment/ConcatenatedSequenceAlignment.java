package edu.vt.vbi.ci.pepr.alignment;

import java.util.Arrays;
import java.util.Random;

import edu.vt.vbi.ci.pepr.stats.StatisticsUtilities;
import edu.vt.vbi.ci.pepr.tree.Bipartition;
import edu.vt.vbi.ci.pepr.tree.BipartitionSet;

/**
 * A SequenceAlignment that has been created by concatenating a set
 * of SequenceAlingments. This class keeps track of the original 
 * SequenceAlignments and where each of them starts and stops in 
 * the concatenated version.
 * 
 * @author enordber
 *
 */
public class ConcatenatedSequenceAlignment extends SequenceAlignment {

	private static Random random = new Random();
	private SequenceAlignment[] alignments;
	private int[][] startStops;

	private int[] stepsBeyondMinimum;
	private float[] stepsToMinimumRatios;

	void setSequenceAlignments(SequenceAlignment[] sa) {
		alignments = new SequenceAlignment[sa.length];
		System.arraycopy(sa, 0, alignments, 0, sa.length);

		startStops = new int[alignments.length][2];
		int nextStart = 0;
		for(int i = 0; i < startStops.length; i++) {
			if(alignments[i] != null) {
				startStops[i][0] = nextStart;
				startStops[i][1] = startStops[i][0] + alignments[i].getLength();
				nextStart = startStops[i][1];
			}
		}
	}

	public Bipartition[] getBipartitionsForAlignment(int alignmentIndex) {
		Bipartition[] r = null;
		//		System.out.println(">ConcatenatedSequenceAlignment.getBipartitionsForAlignment() "
		//				+ alignmentIndex);
		int alignmentStart = startStops[alignmentIndex][0];
		int alignmentStop = startStops[alignmentIndex][1];
		BipartitionSet bipartSet = new BipartitionSet();
		for(int i = alignmentStart; i < alignmentStop; i++) {
			Bipartition[] columnBiparts = getBipartitionsForColumn(i);
			if(columnBiparts != null) {
				for(int j = 0; j < columnBiparts.length; j++) {
					bipartSet.add(columnBiparts[j]);
				}
			}
		}
		r = bipartSet.getBipartitions();
		//		System.out.println("<ConcatenatedSequenceAlignment.getBipartitionsForAlignment() ");

		return r;
	}


	public int getStepsForAlignment(int alignmentIndex) {
		int r = 0;
		int[] steps = getStepsPerSite();
		int alignmentStart = startStops[alignmentIndex][0];
		int alignmentStop = startStops[alignmentIndex][1];

		for(int i = alignmentStart; i < alignmentStop; i++) {
			r += steps[i];
		}
		return r;
	}

	public int getTotalStepsBeyondMinimumForAlignment(int alignmentIndex) {
		int r = 0;
		int[] steps = getStepsBeyondMinimumPerSite();
		int alignmentStart = startStops[alignmentIndex][0];
		int alignmentStop = startStops[alignmentIndex][1];

		for(int i = alignmentStart; i < alignmentStop; i++) {
			r += steps[i];
		}

		return r;
	}

	public float getTotalStepsToMinimumRatioForAlignment(int alignmentIndex) {
		float r = 0;
		float[] steps = getStepsToMinumumRatioPerSite();
		int alignmentStart = startStops[alignmentIndex][0];
		int alignmentStop = startStops[alignmentIndex][1];

		for(int i = alignmentStart; i < alignmentStop; i++) {
			if(steps[i] == steps[i]) {
				r += steps[i];
			}
		}

		return r;
	}

	private int[] getStepsBeyondMinumumForAlignment(int alignmentIndex) {
		int[] r = null;
		int[] steps = getStepsBeyondMinimumPerSite();
		int alignmentStart = startStops[alignmentIndex][0];
		int alignmentStop = startStops[alignmentIndex][1];
		r = new int[alignmentStop-alignmentStart];
		System.arraycopy(steps, alignmentStart, r, 0, r.length);
		return r;
	}

	private float[] getStepsToMinumumRatioPerSite() {
		if(stepsToMinimumRatios == null) {
			int[] steps = getStepsPerSite();
			int[] minSteps = getMinimumStepsPerSite();

			stepsToMinimumRatios = new float[steps.length];
			for(int i = 0; i < stepsToMinimumRatios.length; i++) {
				stepsToMinimumRatios[i] = (float)steps[i]/(float)(minSteps[i]+1);
			}
		}
		return stepsToMinimumRatios;
	}

	private int[] getStepsBeyondMinimumPerSite() {
		if(stepsBeyondMinimum == null) {
			int[] steps = getStepsPerSite();
			int[] minSteps = getMinimumStepsPerSite();

			stepsBeyondMinimum = new int[steps.length];
			for(int i = 0; i < stepsBeyondMinimum.length; i++) {
				stepsBeyondMinimum[i] = steps[i] - minSteps[i];
			}
		}
		return stepsBeyondMinimum;
	}

	public int getThresholdStepsForAlignment(int alignmentIndex, 
			int reps, double alpha) {
		int r = 0;

		int[] steps = getStepsPerSite();
		int[] repSteps = new int[reps];
		int alignmentLength = alignments[alignmentIndex].getLength();

		int concatenatedLength = getLength();
		for(int i = 0; i < reps; i++) {
			//for randomization, to avoid re-use and to avoid using
			//columns from the target alignment
			boolean[] indexUsed = new boolean[concatenatedLength];
			int alignmentStart = startStops[alignmentIndex][0];
			int alignmentStop = startStops[alignmentIndex][1];
			for(int j = alignmentStart; j < alignmentStop; j++) {
				indexUsed[j] = true;
			}

			for(int j = 0; j < alignmentLength; j++) {
				int nextIndex = random.nextInt(concatenatedLength);
				while(indexUsed[nextIndex]) {
					nextIndex = random.nextInt(concatenatedLength);
				}
				indexUsed[nextIndex] = true;
				repSteps[i] += steps[nextIndex];
			}
		}

		int alphaSubtract = (int) Math.ceil(reps*alpha);

		Arrays.sort(repSteps);
		r = repSteps[reps-alphaSubtract];

		return r;
	}

	public int getThresholdStepsBeyondMinimumForAlignment(int alignmentIndex, 
			int reps, double alpha) {
		int r = 0;

		int[] steps = getStepsBeyondMinimumPerSite();
		int[] repSteps = new int[reps];
		int alignmentLength = alignments[alignmentIndex].getLength();

		int concatenatedLength = getLength();
		for(int i = 0; i < reps; i++) {
			//for randomization, to avoid re-use and to avoid using
			//columns from the target alignment
			boolean[] indexUsed = new boolean[concatenatedLength];
			int alignmentStart = startStops[alignmentIndex][0];
			int alignmentStop = startStops[alignmentIndex][1];
			for(int j = alignmentStart; j < alignmentStop; j++) {
				indexUsed[j] = true;
			}

			for(int j = 0; j < alignmentLength; j++) {
				int nextIndex = random.nextInt(concatenatedLength);
				while(indexUsed[nextIndex]) {
					nextIndex = random.nextInt(concatenatedLength);
				}
				indexUsed[nextIndex] = true;
				repSteps[i] += steps[nextIndex];
			}
		}

		int alphaSubtract = (int) Math.ceil(reps*alpha);

		Arrays.sort(repSteps);
		r = repSteps[reps-alphaSubtract];

		return r;
	}

	/**
	 * generates and returns a set of randomly shuffled arrays of the 
	 * stepsBeyondMinimumPerSite. These shuffled arrays can be used to get 
	 * estimates of the probability of finding a particular number of steps
	 * beyond the minimum required for a particular length of alignment.
	 * 
	 * @return
	 */
	public int[][] getShuffledStepsBeyondMinimumPerSite(int replicates, int length) {
		int[][] r = new int[replicates][length];
		int[] steps = getStepsBeyondMinimumPerSite();
		boolean[] unused = new boolean[steps.length];

		return r;
	}

	public int[] getThresholdStepsBeyondMinimumForLengths(int reps, 
			double alpha, int maxLength) {
		int[] r = new int[maxLength];
		int[] shuffledIndices = new int[maxLength];
		int concatenatedLength = getLength();

		return r;
	}

	public int getThresholdStepsBeyondMinimumForAlignment(int alignmentIndex, 
			int reps, double alpha, boolean[] alignmentMask) {
		int r = 0;

		int[] steps = getStepsBeyondMinimumPerSite();
		int[] repSteps = new int[reps];
		int alignmentLength = alignments[alignmentIndex].getLength();

		int concatenatedLength = getLength();
		boolean[] initialIndexMask = new boolean[concatenatedLength];
		for(int i = 0; i < alignments.length; i++) {
			if(i == alignmentIndex || alignmentMask[i]) {
				int alignmentStart = startStops[i][0];
				int alignmentStop = startStops[i][1];
				Arrays.fill(initialIndexMask, alignmentStart, alignmentStop, true);
				//				for(int j = alignmentStart; j < alignmentStop; j++) {
				//					initialIndexMask[j] = true;
				//				}
			}
		}

		//make sure enough columns remain unmasked to calculate a
		//value
		int multiple = 3;
		int unmaskedCount = 0;
		for(int i = 0; i < initialIndexMask.length; i++) {
			if(!initialIndexMask[i]) {
				unmaskedCount++;
			}
		}
		if(unmaskedCount < alignmentLength * multiple) {
			r = -1;
		} else {

			//get subset of steps that are available for this run.
			//This will exclude any steps that are currently masked
			int[] unmaskedSteps = new int[unmaskedCount];
			int unmaskedStepIndex = 0;
			for(int i = 0; i < steps.length; i++) {
				if(!initialIndexMask[i]) {
					unmaskedSteps[unmaskedStepIndex] = steps[i];
					unmaskedStepIndex++;
				}
			}


			//			boolean[] indexUsed = new boolean[unmaskedCount];

			for(int i = 0; i < reps; i++) {
				//for randomization, to avoid re-use
				//				indexUsed = new boolean[unmaskedCount];
				for(int j = 0; j < alignmentLength; j++) {
					int nextIndex = random.nextInt(unmaskedCount);
					//					while(indexUsed[nextIndex]) {
					//						nextIndex = random.nextInt(unmaskedCount);
					//					}
					//					indexUsed[nextIndex] = true;
					repSteps[i] += unmaskedSteps[nextIndex];
				}
			}

			int alphaSubtract = (int) Math.ceil(reps*alpha);

			Arrays.sort(repSteps);
			r = repSteps[reps-alphaSubtract];
		}
		return r;
	}

	public int getThresholdStepsBeyondMinimumForAlignmentOLD(int alignmentIndex, 
			int reps, double alpha, boolean[] alignmentMask) {
		int r = 0;

		int[] steps = getStepsBeyondMinimumPerSite();
		int[] repSteps = new int[reps];
		int alignmentLength = alignments[alignmentIndex].getLength();

		int concatenatedLength = getLength();
		boolean[] initialIndexMask = new boolean[concatenatedLength];
		for(int i = 0; i < alignments.length; i++) {
			if(i == alignmentIndex || alignmentMask[i]) {
				int alignmentStart = startStops[i][0];
				int alignmentStop = startStops[i][1];
				Arrays.fill(initialIndexMask, alignmentStart, alignmentStop, true);
				//				for(int j = alignmentStart; j < alignmentStop; j++) {
				//					initialIndexMask[j] = true;
				//				}

			}
		}

		//make sure enough columns remain unmasked to calculate a
		//value
		int multiple = 2;
		int unmaskedCount = 0;
		for(int i = 0; i < initialIndexMask.length; i++) {
			if(!initialIndexMask[i]) {
				unmaskedCount++;
			}
		}
		if(unmaskedCount < alignmentLength * multiple) {
			r = -1;
		} else {
			boolean[] indexUsed = new boolean[concatenatedLength];

			for(int i = 0; i < reps; i++) {
				//for randomization, to avoid re-use and to avoid using
				//columns from the target alignment
				System.arraycopy(initialIndexMask, 0, indexUsed, 0, indexUsed.length);
				for(int j = 0; j < alignmentLength; j++) {
					int nextIndex = random.nextInt(concatenatedLength);
					while(indexUsed[nextIndex]) {
						nextIndex = random.nextInt(concatenatedLength);
					}
					indexUsed[nextIndex] = true;
					repSteps[i] += steps[nextIndex];
				}
			}

			int alphaSubtract = (int) Math.ceil(reps*alpha);

			Arrays.sort(repSteps);
			r = repSteps[reps-alphaSubtract];
		}
		return r;
	}

	public float getThresholdStepsToMinimumRatioForAlignment(int alignmentIndex, 
			int reps, double alpha, boolean[] alignmentMask) {

		float r = 0;

		float[] stepRatios = getStepsToMinumumRatioPerSite();
		float[] repSteps = new float[reps];
		int alignmentLength = alignments[alignmentIndex].getLength();

		int concatenatedLength = getLength();
		boolean[] initialIndexMask = new boolean[concatenatedLength];
		for(int i = 0; i < alignments.length; i++) {
			if(i == alignmentIndex || alignmentMask[i]) {
				int alignmentStart = startStops[i][0];
				int alignmentStop = startStops[i][1];
				Arrays.fill(initialIndexMask, alignmentStart, alignmentStop, true);
				//				for(int j = alignmentStart; j < alignmentStop; j++) {
				//					initialIndexMask[j] = true;
				//				}

			}
		}

		//make sure enough columns remain unmasked to calculate a
		//value
		int multiple = 2;
		int unmaskedCount = 0;
		for(int i = 0; i < initialIndexMask.length; i++) {
			if(!initialIndexMask[i]) {
				unmaskedCount++;
			}
		}
		if(unmaskedCount < alignmentLength * multiple) {
			r = -1;
		} else {
			boolean[] indexUsed = new boolean[concatenatedLength];

			for(int i = 0; i < reps; i++) {
				//for randomization, to avoid re-use and to avoid using
				//columns from the target alignment
				System.arraycopy(initialIndexMask, 0, indexUsed, 0, indexUsed.length);
				for(int j = 0; j < alignmentLength; j++) {
					int nextIndex = random.nextInt(concatenatedLength);
					while(indexUsed[nextIndex]) {
						nextIndex = random.nextInt(concatenatedLength);
					}
					indexUsed[nextIndex] = true;
					repSteps[i] += stepRatios[nextIndex];
				}
			}

			int alphaSubtract = (int) Math.ceil(reps*alpha);

			Arrays.sort(repSteps);
			r = repSteps[reps-alphaSubtract];
		}
		return r;

	}

	public SequenceAlignment[] getAlignments() {
		return alignments;
	}

	public void testStepDistributions() {
		System.out.println("ConcatenatedAlignment.testStepDistributions()");
		int[] steps = getStepsBeyondMinimumPerSite();
		steps = getMinimumStepsPerSite();
		double[] stepsD = new double[steps.length];
		for(int i = 0; i < stepsD.length; i++) {
			stepsD[i] = steps[i];
		}
		double mean = StatisticsUtilities.getMean(stepsD);
		double sd = StatisticsUtilities.getStandardDeviation(stepsD, mean);
		System.out.println("steps for concatenated alignment: " 
				+ mean + " +/- " + sd);
		System.out.println("length of concatenated alignment: " + steps.length);
		int[] stepDistribution = StatisticsUtilities.getDistribution(steps);
		double[] stepPercent = new double[stepDistribution.length];
		for(int i = 0; i < stepPercent.length; i++) {
			stepPercent[i] = (double)stepDistribution[i]/steps.length;
		}
		//		for(int i = 0; i < stepDistribution.length; i++) {
		//			System.out.println(i + ": " + stepDistribution[i] + "\t"
		//					+ stepPercent[i]);
		//		}
		//		
		//		for(int i = 0; i < steps.length; i++) {
		//			if(steps[i] == 0) {
		//				System.out.println("column " + i + " requires 0 steps: ");
		//			}
		//		}
	}
}
