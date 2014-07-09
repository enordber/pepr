package edu.vt.vbi.ci.pepr.tree;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;

import edu.vt.vbi.ci.pepr.alignment.SequenceAlignment;
import edu.vt.vbi.ci.pepr.stats.StatisticsUtilities;

public class BipartitionSet {

	private String[] taxa;

	/*
	 * Key: Bipartition
	 * value: Integer - number of times this Bipartition has been added
	 */
	private HashMap bipartitionToCount;

	/*
	 *Key: Bipartition
	 *value: Double - ratio of non-conflicting columns to all columns for
	 *                this Bipartition 
	 */
	private HashMap bipartitionToSupport;

	/*
	 * Key: Bipartition
	 * value: Integer - cost for this Bipartition (number of changes needed in
	 *                  complete set to make consistent with this Bipartition)
	 */
	private HashMap bipartitionToCost;

	/*
	 * Tracks if the support values need to be calculated, either because
	 * they have not been calculated yet or because new Bipartitions
	 * have been added since they were last calculated.
	 */
	private boolean supportCalculationNeeded = true;

	private Bipartition[] bipartitions;

	/*
	 * key: a Bipartition
	 * value: ArrayList of all Bipartitions that are supertree-equivalent
	 *        to the key Bipartition
	 *        
	 * The specific Bipartition used as the key is whichever one happens
	 * to be added first.
	 */
	private HashMap equivalentBipartitionToIndividualBipartitionList;

	public BipartitionSet() {
		bipartitionToCount = new HashMap();
		equivalentBipartitionToIndividualBipartitionList = new HashMap();
		bipartitions = new Bipartition[0];
	}

	public BipartitionSet(Bipartition[] biparts) {
		this();
		setBipartitions(biparts, Integer.MAX_VALUE);
	}

	/**
	 * Create a BiparttitionSet with the given Bipartitions. Only the topN most
	 * frequent Bipartitions are kept. topN = biparts.length will keep all. It
	 * is possible that more than topN will be kept. A cutoff frequency is
	 * determined based on topN, and all bipartitions occurring at least that 
	 * many times are kept. The alternative is for some frequency some of the
	 * bipartitions will be kept and others will be removed. This seems less 
	 * sensible than keeping a few more than topN if they all have the same
	 * frequency.
	 * 
	 * @param biparts
	 * @param topN
	 */
	public BipartitionSet(Bipartition[] biparts, int topN) {
		this();
		setBipartitions(biparts, topN);
	}

	private void setBipartitions(Bipartition[] biparts, int topN) {
		//sort bipartitions
		//go through array keeping each unique bipart once, and keeping the 
		//count of the number of times it is seen
		int uniqueBipartCount = 0;
		Bipartition previousBipart = null;
		Arrays.sort(biparts);
		for(int i = 0; i < biparts.length; i++) {
			if(!biparts[i].equals(previousBipart)) {
				uniqueBipartCount++;
				previousBipart = biparts[i];
			}
		}

		bipartitions = new Bipartition[uniqueBipartCount];
		int[] bipartCounts = new int[uniqueBipartCount];
		int index = -1;
		for(int i = 0; i < biparts.length; i++) {
			if(!biparts[i].equals(previousBipart)) {
				index++;
				bipartitions[index] = biparts[i];
				previousBipart = biparts[i];
			}
			bipartCounts[index]++;
		}

		int[] countDistribution = StatisticsUtilities.getDistribution(bipartCounts);

		//determine count cutoff
		int countCutoff = 0;
		int included = 0;
		for(int i = countDistribution.length-1; i >= 0 && included < topN; i--) {
			countCutoff = i;
			included += countDistribution[i];
		}

		Bipartition[] topBiparts = new Bipartition[included];
		int[] topCounts = new int[included];
		index = -1;
		for(int i = 0; i < bipartitions.length; i++) {
			if(bipartCounts[i] >=countCutoff) {
				index++;
				topBiparts[index] = bipartitions[i];
				topCounts[index] = bipartCounts[i];
			}
		}

		bipartitions = topBiparts;
		bipartCounts = topCounts;

		bipartitionToCount = new HashMap();
		for(int i = 0; i < bipartitions.length; i++) {
			Integer count = new Integer(bipartCounts[i]);
			bipartitionToCount.put(bipartitions[i], count);
		}
	}

	public void add(Bipartition[] bipartitions) {
		if(bipartitions != null) {
			for(int i = 0; i < bipartitions.length; i++) {
				add(bipartitions[i]);
			}
		}
	}

	/**
	 * Add a Bipartition to the set. If this Bipartition has already
	 * been added its count is incremented.
	 * 
	 * @param bipartition
	 */
	public void add(Bipartition bipartition) {
		Integer count = (Integer) bipartitionToCount.get(bipartition);

		if(count == null) {
			count = new Integer(1);
			bipartitionToCount.put(bipartition, count);
		} else {
			count = new Integer(count.intValue()+1);
			bipartitionToCount.put(bipartition, count);
		}

		ArrayList equivalentList = 
			(ArrayList)equivalentBipartitionToIndividualBipartitionList.get(bipartition);
		if(equivalentList == null) {
			equivalentList = new ArrayList(1);
			equivalentBipartitionToIndividualBipartitionList.put(bipartition, equivalentList);
		}

		equivalentList.ensureCapacity(equivalentList.size() + 1);
		equivalentList.add(bipartition);

		Bipartition[] newBP = new Bipartition[bipartitions.length+1];
		System.arraycopy(bipartitions, 0, newBP, 0, bipartitions.length);
		newBP[bipartitions.length] = bipartition;
		bipartitions = newBP;
	}

	/**
	 * Sets the list of taxa for the Bipartitions in this set.
	 * The list must be in the correct order for proper 
	 * functioning.
	 * 
	 * @param taxa
	 */
	public void setTaxa(String[] taxa) {
		this.taxa = new String[taxa.length];
		System.arraycopy(taxa, 0, this.taxa, 0, taxa.length);
	}

	/**
	 * Returns all non-trivial bipartitions. A trivial bipartition is
	 * one where one side of the partition has zero or one member.
	 * @return
	 */
	private Bipartition[] getNonTrivialBipartitions() {
		Bipartition[] r = null;
		HashSet biparts = new HashSet();
		for(int i = 0; i < bipartitions.length; i++) {
			biparts.add(bipartitions[i]);
		}

		//remove trivial bipartitions
		for(int i = 0; i < bipartitions.length; i++) {
			int cardinality = bipartitions[i].getSmallerSide().cardinality();
			if(cardinality <= 1 || cardinality >= taxa.length-1) {
				biparts.remove(bipartitions[i]);
			}
		}
		r = new Bipartition[biparts.size()];
		biparts.toArray(r);
		return r;
	}

	/**
	 * Prints information about all non-trivial Bipartitions.
	 */
	public void printNonTrivialBipartitionsAndCounts() {
		printBipartitionsAndCounts(getNonTrivialBipartitions());
	}

	public SequenceAlignment getBipartitionsAsSequenceAlignment() {
		SequenceAlignment r = new SequenceAlignment();
		String present = "1";
		String absent = "0";
		String missing = "?";

		//determine total number of bipartitions (including duplicates)
		Integer[] counts = new Integer[bipartitionToCount.size()];
		bipartitionToCount.values().toArray(counts);
        int length = 0;
        for(int i = 0; i < counts.length; i++) {
        	length += counts[i].intValue();
        }
		
		int taxonCount = taxa.length;

		String title = "";
		StringBuffer sequenceBuffer = null; //for SequenceAlignment
		
		for(int i = 0; i < taxa.length; i++) {
			title = taxa[i];
			sequenceBuffer = new StringBuffer();

			for(int j = 0; j < bipartitions.length; j++) {

				if(bipartitions[j].getSmallerSide().get(i)) {
						sequenceBuffer.append(present);
				} else if(bipartitions[j].getParticipatingTaxonSet().get(i)) {
						sequenceBuffer.append(absent);
				} else {
						sequenceBuffer.append(missing);
				}
			}

			r.addSequence(sequenceBuffer.toString(), title);
		}
		r.setTaxa(taxa);
		return r;
	}

	/**
	 * Checks the given Bipartition for compatibility with all Bipartitions
	 * in this set. These checks take into account the taxa that are present
	 * in the Bipartitions being compared. True is returned if no conflicting

	 * @param bp
	 * @return
	 */
	boolean isCompatible(Bipartition bp) {
		boolean r = true;

		for(int i = 0; r && i < bipartitions.length; i++) {
			r = bipartitions[i].isSupertreeCompatible(bp);
		}
		return r;
	}
	/**
	 * print information about the given list of Bipartitions.
	 * 
	 * @param biparts
	 */
	private void printBipartitionsAndCounts(Bipartition[] biparts) {
		//print taxon list
		for(int i = 0; i < taxa.length; i++) {
			System.out.println(i + ": " + taxa[i]);
		}
		System.out.println();
		Bipartition[] bipartitions = new Bipartition[biparts.length];
		System.arraycopy(biparts, 0, bipartitions, 0, biparts.length);
		System.out.println("sort bipartitions...");
		Arrays.sort(bipartitions, new BipartitionSupportComparator());
		System.out.println("done sorting");
		System.out.println("i\tbipart\tcount\tcost\tsupport\tparticipating\tequivalent\tsmall set");

		for(int i = 0; i < bipartitions.length; i++) {
			StringBuffer sb = new StringBuffer();

			sb.append(i);
			sb.append("\t");
			sb.append(bipartitions[i].getString());
			sb.append("\t");
			sb.append(getCount(bipartitions[i]));
			sb.append("\t");
			sb.append(getCost(bipartitions[i]));
			sb.append("\t");
			sb.append(truncateDouble(getSupport(bipartitions[i]), 10000));
			if(bipartitions[i].getParticipatingTaxonSet() != null) {
				sb.append("\t");
				sb.append(bipartitions[i].getParticipatingTaxonSet().getString(taxa.length));
			}
			if(equivalentBipartitionToIndividualBipartitionList.containsKey(bipartitions[i])) {
				ArrayList equivalentList = 
					(ArrayList) equivalentBipartitionToIndividualBipartitionList.get(bipartitions[i]);
				sb.append("\t");
				sb.append(equivalentList.size());
			} else {
				sb.append("\t");
				sb.append("");
			}
			sb.append("\t");
			sb.append(bipartitions[i].getSmallerSide().toString());
			System.out.println(sb);
		}
	}

	public int getCount(Bipartition bipartition) {
		int r = 0;
		Integer count = (Integer) bipartitionToCount.get(bipartition);
		if(count != null) {
			r = count.intValue();
		}
		return r;
	}

	public Bipartition[] getBipartitions() {
		Bipartition[] r = new Bipartition[bipartitionToCount.size()];
		bipartitionToCount.keySet().toArray(r);
		return r;
	}

	/**
	 * Finds and returns a set of mutually compatible Bipartitions.
	 * Uses a greedy algorithm that removes the least supported Bipartition
	 * within the set until all remaining Bipartitions are compatible.
	 * 
	 * @return
	 */
	public BipartitionSet findCompatibleBipartitionSet() {
		BipartitionSet r = null;
		Bipartition[] biparts = getBipartitions();

		calculateBipartitionCosts();
		biparts = removeConflictorsForSupported(biparts, 0.5);
		int length;
		do {
			length = biparts.length;
			biparts = removeMostConflicted(biparts);
		} while(biparts.length < length);

		printBipartitionsAndCounts(biparts);

		r = new BipartitionSet();
		String[] taxaCopy = new String[taxa.length];
		System.arraycopy(taxa, 0, taxaCopy, 0, taxa.length);
		r.setTaxa(taxaCopy);
		r.bipartitionToSupport = new HashMap();
		r.bipartitionToCost = new HashMap();
		for(int i = 0; i < biparts.length; i++) {
			r.bipartitionToCount.put(biparts[i], 
					bipartitionToCount.get(biparts[i]));
			r.bipartitionToSupport.put(biparts[i],
					bipartitionToSupport.get(biparts[i]));
			r.bipartitionToCost.put(biparts[i], 
					bipartitionToCost.get(biparts[i]));
			r.equivalentBipartitionToIndividualBipartitionList.put(biparts[i], 
					equivalentBipartitionToIndividualBipartitionList.get(biparts[i]));
		}
		r.printNonTrivialBipartitionsAndCounts();
		return r;
	}

	public BipartitionSet findCompatibleBipartitionSet2() {
		System.out.println("BipartitionsSet.findCompatibleBipartitionSet2()");
		BipartitionSet r = null;
		//create a graph with Bipartitions as nodes, and edge between
		//Bipartitions that are compatible, and the edge weights being the sum
		//of the counts of the Bipartitions connected by the edge.

		ArrayList nodePairList = new ArrayList();
		for(int i = 0; i < bipartitions.length; i++) {
			for(int j = i+1; j < bipartitions.length; j++) {
				if(bipartitions[i].isSupertreeCompatible(bipartitions[j]));
				nodePairList.add(new Object[]{bipartitions[i], bipartitions[j]});
			}
		}

		Object[][] nodePairs = new Object[nodePairList.size()][];
		nodePairList.toArray(nodePairs);

		double[] edgeWeights = new double[nodePairs.length];
		for(int i = 0; i < edgeWeights.length; i++) {
			edgeWeights[i] = getCount((Bipartition) nodePairs[i][0]) + 
			getCount((Bipartition) nodePairs[i][1]);
		}
		return r;
	}

	/**
	 * Removes any Bipartitions that conflict with a Bipartition that has
	 * a support value > the given support. 
	 * A support value > 0.5 means that there are more columns directly
	 * supporting the bipartition than columns conflicting with the bipartition.
	 * @param biparts
	 * @return
	 */
	private Bipartition[] removeConflictorsForSupported(Bipartition[] biparts,
			double support) {
		Bipartition[] r = null;
		//load all Bipartitions into a HashSet, making them easier to remove
		//along the way
		HashSet retainedBiparts = new HashSet((int)(biparts.length*1.25));
		for(int i = 0; i < biparts.length; i++) {
			retainedBiparts.add(biparts[i]);
		}

		//loop through the biparts.
		for(int i = 0; i < biparts.length; i++) {
			double thisSupport = getSupport(biparts[i]);
			//if the bipart has support > the cutoff support parameter, look for 
			//conflicting biparts
			if(thisSupport > support) {
				for(int j = 0; j < biparts.length; j++) {
					//if a conflicting bipart is found, remove it from the HashSet
					if(!biparts[i].isCompatible(biparts[j]) && 
							!biparts[i].isSupertreeCompatible(biparts[j])) {
						retainedBiparts.remove(biparts[j]);
					}
				}
			}
		}

		//create a new array for the Bipartitions in the HashSet
		//sort the array and return it
		r = new Bipartition[retainedBiparts.size()];
		retainedBiparts.toArray(r);
		Arrays.sort(r, new BipartitionCountComparator());
		return r;
	}

	/**
	 * Removes the Bipartition with the lowest support value within
	 * the given set. If all of the bipartitions are mutually
	 * compatible nothing is removed and the input array is
	 * returned.
	 * 
	 * @param biparts
	 * @return
	 */
	private Bipartition[] removeMostConflicted(Bipartition[] biparts) {
		Bipartition[] r = null;
		//these are the number of columns supporting each bipartition
		int counts[] = new int[biparts.length];
		for(int i = 0; i < counts.length; i++) {
			counts[i] = getCount(biparts[i]);
		}
		double[] supports = new double[biparts.length];
		double minSupport = 1.0;
		int minSupportIndex = -1;
		for(int i = 0; i < biparts.length; i++) {
			int compatibleCount = 0; 
			int incompatibleCount = 0;
			int compatibleScore = 0;
			int incompatibleScore = 0;
			for(int j = 0; j < biparts.length; j++) {
				boolean compatible = biparts[i].isCompatible(biparts[j]);
				compatible |= biparts[i].isSupertreeCompatible(biparts[j]);
				if(compatible) {
					compatibleCount++;
					compatibleScore += counts[j];
				} else {
					incompatibleCount++;
					incompatibleScore += counts[j];
				}
			}
			supports[i] = 
				(double) counts[i] / (counts[i] + incompatibleScore);
			if(supports[i] < minSupport) {
				minSupport = supports[i];
				minSupportIndex = i;
			}
		}

		if(minSupportIndex >= 0) {
			//create an array with all biparts except the one at minIndex
			r = new Bipartition[biparts.length-1];
			int secondPartLength = biparts.length - minSupportIndex - 1;
			System.arraycopy(biparts, 0, r, 0, minSupportIndex);
			System.arraycopy(biparts, minSupportIndex+1, r, 
					minSupportIndex, secondPartLength);
		} else {
			r = biparts;
		}
		return r;
	}

	private Bipartition[] removeLeastSupported(Bipartition[] biparts) {
		Bipartition[] r = null;

		int minSupportIndex = 0;
		double minSupportValue = 0.0;

		return r;
	}
	/**
	 * Calculates the support values for all Bipartitions. The 
	 * support value is the number of columns in the alignment 
	 * directly supporting the Bipartition divided by the sum of
	 * the columns supporting the Bipartition and the columns 
	 * conflicting with the Bipartition. This does not involve
	 * columns that are compatible with, but do not directly 
	 * support, the Bipartition.
	 */
	private void calculateFullSupports() {
		Bipartition[] biparts = new Bipartition[bipartitionToCount.size()];
		bipartitionToCount.keySet().toArray(biparts);
		//these are the number of columns supporting each bipartition
		int counts[] = new int[biparts.length];
		for(int i = 0; i < counts.length; i++) {
			counts[i] = getCount(biparts[i]);
		}
		double[] directSupports = new double[biparts.length];
		for(int i = 0; i < biparts.length; i++) {
			if((i+1)%100 == 0) {
				System.out.print(".");
				if((i+1)%5000 == 0) {
					System.out.println("\t" + (i+1) + " of " + biparts.length);
				}
			}
			int compatibleCount = 0; 
			int incompatibleCount = 0;
			int compatibleScore = 0;
			int incompatibleScore = 0;
			for(int j = 0; j < biparts.length; j++) {
				boolean compatible = biparts[i].isCompatible(biparts[j]);
				compatible |= biparts[i].isSupertreeCompatible(biparts[j]);
				if(compatible) {
					compatibleCount++;
					compatibleScore += counts[j];
				} else {
					incompatibleCount++;
					incompatibleScore += counts[j];
				}
			}
			directSupports[i] = (double)counts[i] / (double)(counts[i] + incompatibleScore);
		}

		bipartitionToSupport = new HashMap();
		for(int i = 0; i < biparts.length; i++) {
			Double d = new Double(directSupports[i]);
			bipartitionToSupport.put(biparts[i], d);
		}
	}

	/**
	 * Calculates the cost for each Bipartition and stores the values
	 * in the bipartitionToCost HashMap for later retrieval.
	 * 
	 */
	private void calculateBipartitionCosts() {
		Bipartition[] biparts = getBipartitions();
		int[] costs = new int[biparts.length];

		for(int i = 0; i < biparts.length; i++) {
			for(int j = 0; j < biparts.length; j++) {
				boolean isCompatible = biparts[i].isCompatible(biparts[j]);
				if(!isCompatible) {
					costs[i] += getCount(biparts[j]);
				}
			}
		}

		bipartitionToCost = new HashMap((int) (biparts.length * 1.25));
		for(int i = 0; i < biparts.length; i++) {
			bipartitionToCost.put(biparts[i], new Integer(costs[i]));
		}
	}

	private void calculateCostForBipartition(Bipartition bipart) {
		int cost = 0;
		for(int j = 0; j < bipartitions.length; j++) {
			boolean isCompatible = bipart.isCompatible(bipartitions[j]);
			if(!isCompatible) {
				cost += getCount(bipartitions[j]);
			}
		}
		bipartitionToCost.put(bipart, new Integer(cost));
	}

	/**
	 * Returns the support value for the given Bipartition. The 
	 * support value is the number of columns in the alignment 
	 * directly supporting the Bipartition divided by the sum of
	 * the columns supporting the Bipartition and the columns 
	 * conflicting with the Bipartition. This does not involve
	 * columns that are compatible with, but do not directly 
	 * support, the Bipartition.
	 * 
	 * @param bipartition
	 * @return
	 */
	double getSupport(Bipartition bipartition) {
		double r = Double.NaN;
		if(bipartitionToSupport == null) {
			calculateFullSupports();
		}
		Double support = (Double)bipartitionToSupport.get(bipartition);
		if(support == null) {
			r = 0.0;
		} else {
			r = support.doubleValue();
		}
		return r;
	}

	/**
	 * Returns the cost for the given Bipartition. The cost is the number
	 * of Bipartitions that conflict with the given Bipartition. That is,
	 * the number of Bipartitions that have to be discarded in order 
	 * to retain this one.
	 * 
	 * @param bipartition
	 * @return
	 */
	public int getCost(Bipartition bipartition) {
		int r = -1;
		if(bipartitionToCost == null) {
			calculateBipartitionCosts();
		} if(!bipartitionToCost.containsKey(bipartition)) {
			calculateCostForBipartition(bipartition);
		}

		r = ((Integer)bipartitionToCost.get(bipartition)).intValue();
		return r;
	}

	/**
	 * Simple convenience method to remove unwanted digits from a double.
	 * 
	 * @param d
	 * @param factor
	 * @return
	 */
	private double truncateDouble(double d, int factor) {
		double r = 0;
		int trunc = (int) (d * factor);
		r = trunc / (double)factor;
		return r;
	}

	/**
	 * Used to sort Bipartitions based on their count.
	 * @author enordber
	 *
	 */
	private class BipartitionCountComparator implements Comparator {

		public int compare(Object o1, Object o2) {
			int count1 = getCount((Bipartition) o1);
			int count2 = getCount((Bipartition) o2);

			return count2 - count1;
		}
	}

	/**
	 * Used to sort Bipartitions based on support value.
	 * @author enordber
	 *
	 */
	private class BipartitionSupportComparator implements Comparator {

		public int compare(Object o1, Object o2) {
			int r = 0;
			double sup1 = getSupport((Bipartition) o1);
			double sup2 = getSupport((Bipartition) o2);
			double diff = sup2 - sup1;

			if(diff < 0) {
				r = -1;
			} else if(diff > 0) {
				r = 1;
			}
			return r;
		}
	}

	/**
	 * Used to sort Bipartitions based on their cost.
	 * @author enordber
	 *
	 */
	private class BipartitionCostComparator implements Comparator {

		public int compare(Object o1, Object o2) {
			int cost1 = getCost((Bipartition) o1);
			int cost2 = getCost((Bipartition) o2);

			return cost1 - cost2;
		}
	}

}
