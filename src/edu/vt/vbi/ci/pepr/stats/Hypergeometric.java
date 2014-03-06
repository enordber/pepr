package edu.vt.vbi.ci.pepr.stats;

/**
 *  Supports calculation of p values from a hypergeometric
 *  distribution.
 *  Information on hypergeometric distribution obtained from:
 *  Eric W. Weisstein. "Hypergeometric Distribution." From MathWorld--
 *      A Wolfram Web Resource.
 *  http://mathworld.wolfram.com/HypergeometricDistribution.html
 *
 *  Variable names used in this class are taken from the above source.
 *  i = number of successful selections
 *  n = number of selections made
 *  N = number of "good" selections possible
 *  m = number of "bad" selections possible
 *  r = m + N (size of the pool from which items are selected)
 */
public class Hypergeometric
{
     
    /**
     * ((n choose i) * (m choose (N-i))) / (r choose N)
     * 
     * @param i number of successful selections
     * @param n number of selections made
     * @param N number of "good" selections possible
     * @param r size of the pool from which items are selected
     * @return
     */
    public static double getHypergeometricLnPValue(int i, int n, int N, int r) {
    	double lnp = Double.NEGATIVE_INFINITY;
    	
    	//N is the number of 'good' choices, m is the number of
    	//'bad' choices. These numbers can be reversed and the result
    	//will be the same pValue. The calculation below only works
    	//properly when n < m, so make sure that is true.
    	
    	int rMinusN = r - N;
    	int m;
    	
    		m = rMinusN;
    	if(i >= n-m) { 
    		double nChoosei = getLnBinomial(n, i);
    		double mChooseNMinusi = getLnBinomial(m, (N-i));
    		double mPlusnChooseN = getLnBinomial((m+n), N);
    		
    		lnp = nChoosei + mChooseNMinusi - mPlusnChooseN;
    	}
    	return lnp;
    }

    public static double getPValueGreaterThanOrEqualTo(int i, int n, int N, int r) {
		double pSum = 0;
		for( ; (i <= n) && (i <= N); i++) {
		    double lnPValue = getHypergeometricLnPValue(i, n, N, r);
		    double pValue = Math.pow(Math.E, lnPValue);
		    pSum += pValue;
	    }

		return pSum;
    }
    
    public static double getPValueLessThanOrEqualTo(int i, int n, int N, int r) {
		double pSum = 0;
		for( ; i > -1; i--) {
		    double lnPValue = getHypergeometricLnPValue(i, n, N, r);
		    double pValue = Math.pow(Math.E, lnPValue);
		    pSum += pValue;
	    }

		return pSum;
    }

    /**
     * ((n choose i) * (m choose (N-i))) / (r choose N)
     * 
     * @param i number of successful selections
     * @param n number of selections made
     * @param N number of "good" selections possible
     * @param r size of the pool from which items are selected
     * @return
     */
    public static double getHypergeometricPValue(int i, int n, int N, int r) {
    	    double p = Double.NaN;
    	    
        //n is the number of 'good' choices, m is the number of
    	    //'bad' choices. These numbers can be reversed and the result
    	    //will be the same pValue. The calculation below only works
    	    //properly when n < m, so make sure that is true.

    	    int rMinusN = r - N;
    	    int m;
    	    
    	    if(rMinusN < N) {
    	    	    m = N;
    	    	    N = rMinusN;
    	    } else {
    	    	    m = rMinusN;
    	    }
    	        	    
     	double nChoosei = getBinomial(n, i);
     	double mChooseNMinusi = getBinomial(m, (N-i));
   	    double mPlusnChooseN = getBinomial((m+n), N);

     	p = nChoosei * mChooseNMinusi / mPlusnChooseN;
    	
     	return p;
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
    
    /**
     * Returns ln(n choose k).
     * ln(n!) - ln(k! * (n-k)!)
     * @param n
     * @param k
     * @return
     */
    public static double getLnBinomial(int n, int k) {
    	    double r = 0;
    	    double nd = n;
    	    double kd = k;
    	    //validate n and k
    	    if(nd > 0 && kd > 0 && nd >= kd) {
    	    	    double currentValue = 0;
    	    	    
    	    	    while(kd > 0) {
    	    	    	    double nextFactor = Math.log(nd) - Math.log(kd);
    	    	    	    currentValue += nextFactor;
    	    	    	    kd -= 1.0;
    	    	    	    nd -= 1.0;
    	    	    }
    	    	    r = currentValue;
    	    } 
    	
    	    return r;
    }


}
