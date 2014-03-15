package edu.vt.vbi.ci.util;

/**
 * @author ericnordberg
 */
public class IntPair implements Comparable {

	int a, b;
	int hashCode;
	double score = Double.NaN;
	
	public IntPair(int intA, int intB) {
		if(intA > intB) {
			a = intA;
			b = intB;
		} else {
			a = intB;
			b = intA;
		}
		hashCode = a+b;
	}
	
	public IntPair(int intA, int intB, double score) {
		this(intA, intB);
		this.score = score;
	}

	public void setScore(double s) {
		score = s;
	}
	
	public double getScore() {
		return score;
	}
	
	public int compareTo(Object arg0) {
		int r = 0;
		
		IntPair o = (IntPair) arg0;	
		
		r = a - o.a;
		if(r == 0) {
			r = b - o.b;
		}

		return r;
	}
	
	public void setPairValues(int intA, int intB) {
		if(intA > intB) {
			a = intA;
			b = intB;
		} else {
			a = intB;
			b = intA;
		}
		
		hashCode = a+b;
	}
	
	public int[] getPairValues() {
		int[] r = new int[]{a,b};
		return r;
	}
	
	public boolean equals(Object arg0) {
		boolean r = false;
		
		try {
			IntPair o = (IntPair) arg0;
			r = (a == o.a && b == o.b);
		} catch(Exception e) {}
		
		return r;
	}

	public int hashCode() {
		return hashCode;
	}
}
