package edu.vt.vbi.ci.pepr.tree;

import edu.vt.vbi.ci.util.ExtendedBitSet;

public class Bipartition implements Comparable{

	/*
	 * participatingTaxonSet supports supertree construction.
	 * A particular bipartition may be determined based on a
	 * subset of the taxa in the full analysis. This BitSet
	 * records which taxa were actually present in the data
	 * used to construct this Bipartition. This is used in
	 * compatibility tests with Bipartitions from other,
	 * overlapping, taxon sets. 
	 * 
	 * Compatibility tests change from:
	 * is A compatible with B 
	 * to:
	 *     is (A & B.participatingTaxonSet) compatible 
	 *   with (B & A.participatingTaxonSet) 
	 */
	private ExtendedBitSet participatingTaxonSet;
	private boolean useParticipatingTaxonSetForEqualityCheck = false;

	private ExtendedBitSet smallerSide;
	private ExtendedBitSet largerSide;
	private int size;
	private String bipartitionString;

	/**
	 * This constructor can be used when the two sides of the bipart 
	 * are known to be non-overlapping. That is, one side can be provided
	 * and the other side is the complement of the provided side. The size 
	 * must be provided to allow proper construction of the complement. When the
	 * two sides of the bipart (not technically correct) are not known to
	 * be non-overlapping, use the constructor that takes both sides of the 
	 * bipart as separate ExtendedBitSet parameters.
	 * @param bipart
	 * @param size
	 */
	public Bipartition(ExtendedBitSet bipart, int size) {
		this.size = size;
		ExtendedBitSet complement = bipart.getComplement(size);
		int bpCard = bipart.cardinality();
		int compCard = complement.cardinality();

		if(bpCard < compCard) {
			smallerSide = bipart;
			largerSide = complement;
		} else if(bpCard > compCard) {
			smallerSide = complement;
			largerSide = bipart;
		} else {
			int bpFirst = bipart.nextSetBit(0);
			int compFirst = complement.nextSetBit(0);
			if(bpFirst < compFirst) {
				smallerSide = bipart;
				largerSide = complement;
			} else {
				smallerSide = complement;
				largerSide = bipart;
			}
		}
	}

	public Bipartition(ExtendedBitSet sideA, ExtendedBitSet sideB) {
		int cardA = sideA.cardinality();
		int cardB = sideB.cardinality();
		if(cardA < cardB) {
			smallerSide = sideA;
			largerSide = sideB;
		} else if(cardB < cardA){
			largerSide = sideA;
			smallerSide = sideB;
		} else {
			int sideAFirst = sideA.nextSetBit(0);
			int sideBFirst = sideB.nextSetBit(0);
			if(sideAFirst < sideBFirst) {
				smallerSide = sideA;
				largerSide = sideB;
			} else {
				smallerSide = sideB;
				largerSide = sideA;
			}
		}

		size = Math.max(smallerSide.length(), largerSide.length());
	}
	/**
	 * This constructor allows a size to be specified that may be different
	 * from the size of the union of the two sides. This is useful when this 
	 * partition is part of a set of Bipartitons that may have different sizes.
	 * Its best if the size is not smaller than the union of the two sides.
	 * @param sideA
	 * @param sideB
	 */
	public Bipartition(ExtendedBitSet sideA, ExtendedBitSet sideB, int size) {
		this(sideA, sideB);
		this.size = size;
	}

	public ExtendedBitSet getSmallerSide() {
		return smallerSide;
	}

	public ExtendedBitSet getLargerSide() {
		return largerSide;
	}

	/**
	 * Bipartitions are equal if the both BitSets are equal.
	 */
	public boolean equals(Object other) {
		boolean r = false;
		if(other != null) {
			Bipartition otherBP = (Bipartition)other;
			r = this.getSmallerSide().equals(otherBP.getSmallerSide())
			&& this.getLargerSide().equals(otherBP.getLargerSide());
			if(r && useParticipatingTaxonSetForEqualityCheck) {
				r = isSupertreeEquivalent(otherBP);
			} 
		}
		return r;
	}

	public boolean isSupertreeEquivalent(Bipartition otherBP) {
		boolean r = false;
		r = this.getParticipatingTaxonSet().equals( 
				otherBP.getParticipatingTaxonSet());

		return r;
	}

	/**
	 * The hash code of a Bipartition is the hash code of the 
	 * smallerSide BitSet.
	 */
	public int hashCode() {
		int hashCode = getSmallerSide().hashCode();
		if(useParticipatingTaxonSetForEqualityCheck &&  participatingTaxonSet != null) {
			hashCode += participatingTaxonSet.hashCode();
		}
		return hashCode;
	}

	/**
	 * Returns true if the given Bipartition and this Bipartition
	 * may coexist in the same tree.
	 * 
	 * @param otherBP
	 * @return
	 */
	public boolean isCompatible(Bipartition otherBP) {
		boolean r = false;

		//Compatibility is tested by comparing the smallerSides from the
		//two Bipartitions
		//The two possibilities for these biparts to be compatible are:
		//1. one smallerSide is a subset of the other
		//2. the two smallerSides do not overlap
		//These are tested by getting the AND of the two smallerSides
		// If the AND is equal to one of the original smaller sides,
		//then one is a subset of the other. If the AND is zero (empty set),
		//then the two do not intersect. If neither of these is true, then
		//the bipartitions conflict with each other.

		ExtendedBitSet andSet = 
			this.getSmallerSide().getAnd(otherBP.getSmallerSide());
		r = andSet.cardinality() == 0;

		if(!r) {
			r = andSet.equals(this.getSmallerSide()) || 
			andSet.equals(otherBP.getSmallerSide());
		}

		return r;
	}

	/**
	 * Calculates and returns the number of changes required to make
	 * this Bipartition compatible with the given Bipartition.
	 * A cost of zero indicates the two Bipartitions are compatible.
	 * 
	 * @param otherBP
	 * @return
	 */
	public int getCompatibilityCost(Bipartition otherBP) {
		int r = 0;

		if(!this.isCompatible(otherBP)) {
			//			ExtendedBitSet andNotSetSmaller = 
			//				this.getSmallerSide().getAndNot(otherBP.getSmallerSide());
			//			ExtendedBitSet andNotSetLarger = 
			//				this.getSmallerSide().getAndNot(otherBP.getLargerSide());

			ExtendedBitSet andNotSetSmaller = null;
			ExtendedBitSet andNotSetLarger = null;
			if(useParticipatingTaxonSetForEqualityCheck) {
				andNotSetSmaller = 
					this.getSmallerSide().getAnd(otherBP.getParticipatingTaxonSet())
					.getAndNot(otherBP.getSmallerSide().getAnd(this.getParticipatingTaxonSet()));
				andNotSetLarger = 
					this.getSmallerSide().getAnd(otherBP.getParticipatingTaxonSet()).
					getAndNot(otherBP.getLargerSide().getAnd(this.getParticipatingTaxonSet()));
			} else {
				andNotSetSmaller = 
					this.getSmallerSide()
					.getAndNot(otherBP.getSmallerSide());
				andNotSetLarger = 
					this.getSmallerSide().
					getAndNot(otherBP.getLargerSide());

			}
			int smallCard = andNotSetSmaller.cardinality();
			int largeCard = andNotSetLarger.cardinality();

			r = Math.min(smallCard, largeCard);


			r = Math.min(r,
					this.getSmallerSide().
					getAnd(otherBP.getSmallerSide()).cardinality());

			r = Math.min(r,
					this.getSmallerSide().
					getAnd(otherBP.getLargerSide()).cardinality());

		}
		return r;
	}

	public boolean isSupertreeCompatible(Bipartition otherBP) {
		boolean r = false;

		//test if the taxon subsets intersect. if they do not, then
		//these bipartitions are compatible

		r = this.getParticipatingTaxonSet()
		.getAnd(otherBP.getParticipatingTaxonSet())
		.cardinality() == 0;
		if(!r) {
			ExtendedBitSet thisAdjustedSmallerSide = 
				this.getSmallerSide()
				.getAnd(otherBP.getParticipatingTaxonSet());
			ExtendedBitSet thisAdjustedLargerSide = 
				this.getLargerSide()
				.getAnd(otherBP.getParticipatingTaxonSet());
			ExtendedBitSet otherAdjustedSmallerSide = 
				otherBP.getSmallerSide()
				.getAnd(this.getParticipatingTaxonSet());
			ExtendedBitSet otherAdjustedLargerSide = 
				otherBP.getLargerSide()
				.getAnd(this.getParticipatingTaxonSet());

			ExtendedBitSet testSet = 
				thisAdjustedSmallerSide.getAnd(otherAdjustedSmallerSide);
			//		    System.out.println("a " + r);
			r = testSet.cardinality() == 0;
			//		    System.out.println("b " + r);
			r = r || testSet.cardinality()
			== thisAdjustedSmallerSide.cardinality();
			//		    System.out.println("c " + r);
			r = r || testSet.cardinality() 
			== otherAdjustedSmallerSide.cardinality();
			//		    System.out.println("d " + r);		    
		}
		//		System.out.println(this.getSmallerSide());
		//		System.out.println(otherBP.getSmallerSide());
		//		System.out.println("BPSet.isSupertreeCompatible(): " + r);
		return r;
	}

	public String getString() {
		if(bipartitionString == null) {
			ExtendedBitSet bs = getSmallerSide();
			char present = '*';
			char absent = '.';

			StringBuffer sb = new StringBuffer();
			for(int i = 0; i < size; i++) {
				if(bs.get(i)) {
					sb.append(present);
				} else {
					sb .append(absent);
				}
			}
			bipartitionString = sb.toString();
		}
		return bipartitionString;
	}

	ExtendedBitSet getParticipatingTaxonSet() {
		if(participatingTaxonSet == null) {
			//no participating taxon set has been provided, so assume all taxa
			//are participating
			participatingTaxonSet = new ExtendedBitSet();
			participatingTaxonSet.set(0, size);
			
		}
		return participatingTaxonSet;
	}

	public void setParticipatingTaxonSet(ExtendedBitSet participatingTaxonSet) {
		this.participatingTaxonSet = participatingTaxonSet;
	}
	void setUseParticipatingTaxonSetForEqualityCheck(boolean b) {
		useParticipatingTaxonSetForEqualityCheck = b;
	}

	public int compareTo(Object o) {
		int r = 0;
		Bipartition otherBP = (Bipartition) o;
		r = otherBP.hashCode() - this.hashCode();
		return r;
	}
	
	public int getSize() {
		return size;
	}

}
