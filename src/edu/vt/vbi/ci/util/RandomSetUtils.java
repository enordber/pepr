package edu.vt.vbi.ci.util;

import java.util.Arrays;
import java.util.Random;

import edu.vt.vbi.ci.pepr.stats.Hypergeometric;

public class RandomSetUtils {

	
	public static int[] getRandomSet(int length, int minValue, int maxValue,
			boolean allowReuse) {
		int[] r = null;
		int range = maxValue-minValue+1; //+1 so both extremes are included
		if(allowReuse || range >= length) {
			//it is ok to proceed. If these conditions are not met, then an
			//infinite loop will occur
			
			Random random = new Random();
			r = new int[length];
			
			boolean[] used = new boolean[range];
			int nextPosition; //value between 0 and range
			int nextValue; //value between minValue and maxValue
			
			for(int i = 0; i < r.length; i++) {
				nextPosition = random.nextInt(range);
				while(!allowReuse && used[nextPosition]) {
					nextPosition = random.nextInt(range);
				}
				used[nextPosition] = true;
				nextValue = minValue + nextPosition;
				r[i] = nextValue;
			}
		}
		return r;
	}
	
	/**
	 * Returns an array of random sets, with no set being exactly
	 * duplicated (although there is no guaranteed amount of difference
	 * between sets). All sets will be returned in sorted order, but will 
	 * contain a random composition of elements.
	 * When allowResuse is false, not only must (maxValue-minValue) 
	 * be >= length, but it must be sufficiently greater to allow the 
	 * requested number of different sets. That is: 
	 * (maxValue-minValue) choose length >= sets.
	 * 
	 * @param sets number of sets
	 * @param length length of each set
	 * @param minValue inclusive
	 * @param maxValue inclusive
	 * @param allowReuse indicates if a value can occur more than once
	 *                    within a set. This only applies within-set.
	 * @return
	 */
	public static int[][] getDifferentRandomSets(int sets, int length, 
			int minValue, int maxValue, boolean allowReuse) {
		int[][] r = null;
		
		int range = maxValue-minValue;
		int maxDiffSets = (int) getBinomial(range, length);
		if(allowReuse || maxDiffSets >= sets) {
			r = new int[sets][];
			for(int i = 0; i < r.length; i++) {
				int[] newSet = null;
				
				//see if this is a unique set
				boolean isUnique = false;
				while(!isUnique) {
					newSet = getRandomSet(length, minValue, 
							maxValue, allowReuse);
					Arrays.sort(newSet);
					isUnique = true;
					for(int j = 0; isUnique && j < i; j++) {
						isUnique = !Arrays.equals(r[j], newSet);
					}
				}
				r[i] = newSet;
			}
		}
		return r;
	}
	
    /**
     * Returns n choose k.
     * n! / (k! * (n-k)!)
     * @param n
     * @param k
     * @return
     */
    public static double getBinomial(int n, int k) {
    	    double r = 0;
    	    //validate n and k
    	    if(n > 0 && k > 0 && n >= k) {
    	    	    double currentValue = 1;
    	    	    
    	    	    while(k > 0) {
    	    	    	    double nextFactor = (double)n / (double)k;
    	    	    	    currentValue *= nextFactor;
    	    	    	    k--;
    	    	    	    n--;
    	    	    }
    	    	    r = currentValue;
    	    } else if(k == 0) {
    	    	    //n choose k is 1 if k = 0
    	    	    r = 1;
    	    }
    	
    	    return r;
    }

}
