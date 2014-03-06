/*
 * Created on Oct 14, 2006
 */
package edu.vt.vbi.ci.pepr.stats;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

/**
 * @author ericnordberg
 */
public class StatisticsUtilities {
	
	/**
	 * Returns the correlation coefficient for the two arrays.
	 * 
	 * @param x
	 * @param y
	 * @return
	 */
	public static double getRSquared(double[] x, double[] y) {
		double xMean = getMean(x);
		double yMean = getMean(y);
		double SSxy = sumOfSquares(x, y, xMean, yMean);
		double SSxx = sumOfSquares(x, xMean);
		double SSyy = sumOfSquares(y, yMean);
		
		double rSquared = (SSxy * SSxy) / (SSxx * SSyy);
//		System.out.println("xMean: " + xMean);
//		System.out.println("yMean: " + yMean);
//		System.out.println("SSxx: " + SSxx);
//		System.out.println("SSyy: " + SSyy);
//		System.out.println("SSxy: " + SSxy);
		if(SSxy < 0) {
			rSquared = -rSquared;
		}
		return rSquared;
	}
	
	/**
	 * Returns the correlation coefficient for the two 
	 * two-dimension arrays.
	 * 
	 * @param x
	 * @param y
	 * @return correlation coefficient
	 */
	public static double getRSquared(double[][] x, double[][] y) {
		double[] xMean = getCentroidVector(x);
		double[] yMean = getCentroidVector(y);
		double SSxy = sumOfSquaresManhattan(x, y, xMean, yMean);
		double SSxx = sumOfSquaresManhattan(x, xMean);
		double SSyy = sumOfSquaresManhattan(y, yMean);
		
		double rSquared = (SSxy * SSxy) / (SSxx * SSyy);
//		System.out.println("xMean[0]: " + xMean[0]);
//		System.out.println("yMean[0]: " + yMean[0]);
//		System.out.println("SSxx: " + SSxx);
//		System.out.println("SSyy: " + SSyy);
//		System.out.println("SSxy: " + SSxy);
		if(SSxy < 0) {
			rSquared = -rSquared;
		}
		return rSquared;
	}
	
	private static double sumOfSquaresManhattan(double[][] d, double[] mean) {
		double sumOfSquaresOfDistances = 0;
		for(int i = 0; i < d.length; i++) {
		    //make sure every element in this array is 
			//not NaN
			boolean noNaNFoundYet = true;
			for(int j = 0; noNaNFoundYet && j < d[i].length; j++) {
//				instead of using Double.isNaN(), avoid the method call
				//by taking advantage of the fact that NaN != NaN
				if(d[i][j] != d[i][j]) {
					noNaNFoundYet = false;
				}
			}
			
			if(noNaNFoundYet) {
				//determine distance between this element of d and
				//the mean. I'm using Manhattan distance, so this
				//will be the sum of the differences between the 
				//individual elements
				double sumOfDistances = 0;
				for(int j = 0; j < mean.length; j++) {
					double thisDistance = Math.abs(mean[j] - d[i][j]);
					double squareOfDistance = thisDistance * thisDistance;
					sumOfSquaresOfDistances += squareOfDistance;
				}
				
			}
		}
		return sumOfSquaresOfDistances;
	}

	private static double sumOfSquaresManhattan(double[][] x, double[][] y, double[] meanX, double[] meanY) {
		double sumOfProductsOfDistances = 0;
		for(int i = 0; i < x.length; i++) {
		    //make sure every element in this array is 
			//not NaN
			boolean noNaNFoundYet = true;
			for(int j = 0; noNaNFoundYet && j < x[i].length; j++) {
//				instead of using Double.isNaN(), avoid the method call
				//by taking advantage of the fact that NaN != NaN
				if(x[i][j] != x[i][j] || y[i][j] != y[i][j]) {
					noNaNFoundYet = false;
				}
			}
			
			if(noNaNFoundYet) {
				//determine distance between this element of d and
				//the mean. I'm using Manhattan distance, so this
				//will be the sum of the differences between the 
				//individual elements
				double sumOfDistances = 0;
				for(int j = 0; j < meanX.length; j++) {
					double thisXDistance = meanX[j] - x[i][j];
					double thisYDistance = meanY[j] - y[i][j];
					double productOfDistance = thisXDistance * thisYDistance;
					sumOfProductsOfDistances += productOfDistance;
				}
				
			}
		}
		return sumOfProductsOfDistances;
	}

	/**
	 * Returns the sum of squares of the values in d.
	 * This is the sum of the squares of the difference
	 * between each value and the mean of all values. IF any
	 * NaN values are in the vector d, they are ignored.
	 */
	public static double sumOfSquares(double[] d, double mean) {
		double r = 0;
		for(int i = 0; i < d.length; i++) {
			//instead of using Double.isNaN(), avoid the method call
			//by taking advantage of the fact that NaN != NaN

			if(d[i] == d[i]) {
				double diff = d[i] - mean;
				r += diff * diff;
			}
		}
		return r;
	}
	
	/**
	 * Calculates and returns the mean of the
	 * values in d. Any NaN values are ignored.
	 * 
	 * @param d
	 * @return
	 */
	public static double getMean(double[] d) {
		double r = 0.0;
		int denominator = 0;
		for(int i = 0; i < d.length; i++) {
			//instead of using Double.isNaN(), avoid the method call
			//by taking advantage of the fact that NaN != NaN

			if(d[i] == d[i]) {
				r += d[i];
				denominator++;
			}
		}
		
		r /= denominator;
		return r;
	}
	
	/**
	 * Returns the minimum value in the given array.
	 */
	public static double getMinimum(double[] d) {
		double r = Double.MAX_VALUE;
		
		for(int i = 0; i < d.length; i++) {
			double di = d[i];
			if(di == di) { //this is false for NaN
			    r = Math.min(r, di);
			}
		}
		
		return r;
	}
	
	/**
	 * Returns the maximum value in the given array.
	 */
	public static double getMaximum(double[] d) {
		double r = Double.NEGATIVE_INFINITY;
		
		for(int i = 0; i < d.length; i++) {
			double di = d[i];
			if(di == di) { //this is false for NaN
			    r = Math.max(r, di);
			}
		}
		
		return r;
	}
	
	/**
	 * Calculates and returns the mean of the
	 * values in f. Any NaN values are ignored.
	 * 
	 * @param f
	 * @return
	 */
	public static float getMean(float[] f) {
		float r = 0.0f;
		int denominator = 0;
		for(int i = 0; i < f.length; i++) {
			//instead of using Float.isNaN(), avoid the method call
			//by taking advantage of the fact that NaN != NaN

			if(f[i] == f[i]) {
				r += f[i];
				denominator++;
			}
		}
		
		r /= denominator;
		return r;
	}
	
	
	
	/**
	 * Returns an array of the means of the given array 
	 * of doubles. The given array must be balanced. The
	 * returned array will have the same length as the 
	 * second dimension.
	 * 
	 * @param d
	 * @return
	 */
	public static double[] getMean(double[][] d) {
		double[] r = null;
		int length = d[0].length;
		double[] sums = new double[length];
		r = new double[length];
		
		for(int i = 0; i < length; i++) {
			sums[i] = 0;
		}
		
		for(int i = 0; i < d.length; i++) {
			for(int j = 0; j < length; j++) {
				sums[j] += d[i][j];
			}
		}
		
		for(int i = 0; i < length; i++) {
			r[i] = sums[i]/d.length;
		}
		
		return r;
	}
	
	/**
	 * Returns the sum of squares of the two paired double
	 * array.
	 * 
	 * @param x
	 * @param y
	 * @param xMean
	 * @param yMean
	 * @return
	 */
	public static double sumOfSquares(double[] x, double[] y, double xMean, double yMean) {
		double r = 0.0;
		for(int i = 0; i < x.length; i++) {
			double diffX = x[i] - xMean;
			double diffY = y[i] - yMean;
			//instead of using Double.isNaN(), avoid the method call
			//by taking advantage of the fact that NaN != NaN
			if(diffX == diffX && diffY == diffY) {
			    r += diffX * diffY;
			}
		}
		return r;
	}
	
	/**
	 * Calculates and returns a centroid vector or the given set of 
	 * vectors. Each element in the centroid is the mean of that 
	 * element in the input vectors. All input vectors must be the
	 * same length. NaN values are ignored in calculating the
	 * centroid vector. The mean of any non-NaN values is used
	 * for each element in the centroid.
	 * 
	 * @param vectors
	 * @return
	 */
	public static double[] getCentroidVector(double[][] vectors) {
		double[] r = null;
		
		int vectorLength = vectors[0].length;
		int vectorCount = vectors.length;
		r = new double[vectorLength];
		
		for(int i = 0; i < vectorLength; i++) {
			//calculate the mean for the ith element
			int denominator = 0;
			
			double sum = 0.0;
			for(int j = 0; j < vectorCount; j++) {
				double value = vectors[j][i];
				//instead of using Double.isNaN(), avoid the method call
				//by taking advantage of the fact that NaN != NaN

				if(value == value) {
					sum += value;
					denominator++;
				}
				}
			
			double mean = sum/denominator;
			r[i] = mean;
		}
		
		return r;
	}
	
	/**
	 * Returns the sum of squares distance between the
	 * two input vectors. Input vectors must be the same
	 * length. If any position has NaN ineither vector, that
	 * position is ignored for distance calculation.
	 * 
	 * @return
	 */
	public static double getDistanceBetweenVectors(double[] x, double[] y) {
		double r = 0.0;
		
		for(int i = 0; i < x.length; i++) {
			//the product will be NaN if either factor is NaN,
			//and I want to ignore if either factor is NaN
			//instead of using Double.isNaN(), avoid the method call
			//by taking advantage of the fact that NaN != NaN
			double product = x[i] * y[i];

			if(product == product) {
				double diff = x[i] - y[i];
				r += (diff * diff);
			}
		}
		
		return r;
	}
	
	/**
	 * Performs k-means clustering in the input vectors.
	 * The returned array is the same length as the
	 * number of input vectors. The returned[i] is
	 * the number of the cluster to which the vector
	 * vectors[i] belongs.
	 * 
	 * All input vectors must be the same length.
	 * 
	 * The number of vectors must be >= k. If it's equal
	 * to k, that's pretty stupid; each vector is its own
	 * cluster.
	 *  
	 * @param k
	 * @param vectors
	 * @return
	 */
	public static int[] getKMeansClusterAssignments(int k, double[][] vectors) {
		System.out.println("StatisticsUtilities.getKMeansClusterAssignments()");
		int[] r = null;
		Random random = new Random();
		
		int vectorCount = vectors.length;
		int vectorLength = vectors[0].length;
		
		//need initial centroids
		System.out.println("get initial centoids...");
		double[][] centroids = //getRandomInitialCentroids(vectors, k);
			getDistantInitialCentroids(vectors, k);
		
		System.out.println("perform initial cluster assignments...");
		//do initial assignment of vectors to centroids
		int[] assignments = new int[vectorCount];
		for(int i = 0; i < vectorCount; i++) {
			//get the distance between this vector and each centroid
			//store the index of the closets centroid in assignments[i]
			double minDistanceSoFar = Double.POSITIVE_INFINITY;
			double[] thisVector = vectors[i];
			for(int j = 0; j < k; j++) {
				double thisDistance =
					getDistanceBetweenVectors(thisVector, centroids[j]);
				if(thisDistance < minDistanceSoFar) {
					minDistanceSoFar = thisDistance;
					assignments[i] = j;
				}
			}
		}
		
		//System.out.println("initial assignments: ");
		//for(int i = 0; i < assignments.length; i++) {
		//	System.out.println("" + i + ": " + assignments[i]);
		//}
		
		double[][][] clusters = getClusters(vectors, assignments);
		//System.out.println("initial:");
		//printClusters(clusters);
		
		int[] previousAssignments = new int[0]; 
		int round = 0;
		System.out.println("perform cluster refinement...");
		while(!Arrays.equals(previousAssignments, assignments)) {
			System.out.println("round " + round);
			previousAssignments = assignments;
			//until assignments no longer change
			//calculate new centroids
			for(int i = 0; i < k; i++) {
				if(clusters.length > 0 && clusters[i].length > 0) {
					centroids[i] = getCentroidVector(clusters[i]);
				}
			}
			//determine new assignments of vectors to centroids
			assignments = new int[vectorCount];
			for(int i = 0; i < vectorCount; i++) {
				//get the distance between this vector and each centroid
				//store the index of the closest centroid in assignments[i]
				double minDistanceSoFar = Double.POSITIVE_INFINITY;
				double[] thisVector = vectors[i];
				for(int j = 0; j < k; j++) {
					double thisDistance =
						getDistanceBetweenVectors(thisVector, centroids[j]);
					if(thisDistance < minDistanceSoFar) {
						minDistanceSoFar = thisDistance;
						assignments[i] = j;
					}
				}
			}
			
			clusters = getClusters(vectors, assignments);
			
			//printClustersAndCentroids(clusters, centroids);
			round++;
		}
		
		//print clusters
		//printClusters(clusters);
		
		r = assignments;
		return r;
	}

	/**
	 * Performs k-means clustering in the input vectors, based
	 * on correlation rather than vector distance.
	 * The returned array is the same length as the
	 * number of input vectors. The returned[i] is
	 * the number of the cluster to which the vector
	 * vectors[i] belongs.
	 * 
	 * All input vectors must be the same length.
	 * 
	 * The number of vectors must be >= k. If it's equal
	 * to k, that's pretty stupid; each vector is its own
	 * cluster.
	 *  
	 * @param k
	 * @param vectors
	 * @return
	 */
	public static int[] getKMeansCorrelationClusterAssignments(int k, double[][] vectors) {
		System.out.println("StatisticsUtilities.getKMeansCorrelationClusterAssignments() k = " + k);
		int[] r = null;
		//Random random = new Random();
		
		int vectorCount = vectors.length;
		//int vectorLength = vectors[0].length;
		
		//need initial centroids
		System.out.println("get initial centoids...");
		double[][] centroids = //getRandomInitialCentroids(vectors, k);
			getDistantInitialCentroids(vectors, k);
		
//		printCentroids(centroids);
		
		System.out.println("perform initial cluster assignments...");
		//do initial assignment of vectors to centroids
		int[] assignments = new int[vectorCount];
		for(int i = 0; i < vectorCount; i++) {
			//get the distance between this vector and each centroid
			//store the index of the closets centroid in assignments[i]
			double maxRSquaredSoFar = -2;
			double[] thisVector = vectors[i];
			for(int j = 0; j < k; j++) {
				double thisRSquared =
					getRSquared(thisVector, centroids[j]);
				
//				System.out.println("vector " + i + ", centoid " + j +
//						": r2 = " + thisRSquared);
				if(thisRSquared > maxRSquaredSoFar) {
					maxRSquaredSoFar = thisRSquared;
					assignments[i] = j;
				}
			}
		}
		
//		System.out.println("initial assignments: ");
//		for(int i = 0; i < assignments.length; i++) {
//			System.out.println("" + i + ": " + assignments[i]);
//		}
		
		double[][][] clusters = getClusters(vectors, assignments);
		//System.out.println("initial:");
		//printClusters(clusters);
		
		int[] previousAssignments = new int[0]; 
		int round = 0;
		System.out.println("perform cluster refinement...");
		while(!Arrays.equals(previousAssignments, assignments)) {
			System.out.println("round " + round);
			previousAssignments = assignments;
			//until assignments no longer change
			//calculate new centroids
			for(int i = 0; i < k; i++) {
				try {
				centroids[i] = getCentroidVector(clusters[i]);
				} catch(ArrayIndexOutOfBoundsException aioobe) {
					System.out.println("ArrayIndexOutOfBoundsException");
					System.out.println("i: " + i);
					System.out.println("centroids.length: " + centroids.length);
					System.out.println("clusters.length: " + clusters.length);
				    aioobe.printStackTrace();
				}
			}
			//determine new assignments of vectors to centroids
			assignments = new int[vectorCount];
			for(int i = 0; i < vectorCount; i++) {
				//get the distance between this vector and each centroid
				//store the index of the closets centroid in assignments[i]
				double maxRSquaredSoFar = -2;
				double[] thisVector = vectors[i];
				for(int j = 0; j < k; j++) {
					double thisRSquared =
						getRSquared(thisVector, centroids[j]);
					if(thisRSquared > maxRSquaredSoFar) {
						maxRSquaredSoFar = thisRSquared;
						assignments[i] = j;
					}
				}
			}
			
			clusters = getClusters(vectors, assignments);
			
			//printClustersAndCentroids(clusters, centroids);
			round++;
		}
		
		
		//print clusters
		//printClusters(clusters);
		
		r = assignments;
		return r;
	}

	
	/**
	 * The goal of this method is to choose initial centroids
	 * that are at the boundries of the vector space. I want the 
	 * starting centroids to be far away from each other. I think
	 * this will produce better clusters than selecting initial
	 * centroids at random.
	 * 
	 * This works by first determining the centroid of all vectors.
	 * This centroid is not one of the k initial centroids, it is
	 * just used to choose the k initial centroids. I'll call it the 
	 * 'seed' centroid, just to have a name for it.
	 * 
	 * Then, to choose the first initial centroid, find the vector 
	 * that has the greatest distance from the seed centroid. 
	 * The second of the k centroids will be the one farthest from
	 * the mean of the the seed and the first. This will be repeated
	 * until k centroids (in addition to the seed) have been choosen.
	 * 
	 * @param vectors
	 * @param k
	 * @return
	 */
	public static double[][] getDistantInitialCentroids(double[][] vectors, int k) {
		double[][] centroids = new double[k][];
		
		int vectorCount = vectors.length;
		int vectorLength = vectors[0].length;
		
		//use a boolean array to track which vectors have been choosen as centroids
		boolean[] assigned = new boolean[vectorCount];
		ArrayList unassignedVectors = new ArrayList(vectorCount);
		for(int i = 0; i < vectorCount; i++) {
			unassignedVectors.add(vectors[i]);
		}
		
		//get the seed vector - centroid of all vectors
		double[] seedCentroid = getCentroidVector(vectors);
		double[] referenceCentroid = seedCentroid;
		ArrayList centroidsWithSeed = new ArrayList(k+1);
		centroidsWithSeed.add(seedCentroid);
		
		//get the first k centroids, that are farthest from each other 
		//(or at least close to being farthest from each other)
		for(int i = 0; i < k; i++) {
			//get the reference centroid, which is the mean of
			//the seed vector and all previously selected centroids
			double[][] centroidsWithSeedVectors = 
				new double[centroidsWithSeed.size()][];
			centroidsWithSeed.toArray(centroidsWithSeedVectors);
			referenceCentroid = getCentroidVector(centroidsWithSeedVectors);
			double[] newCentroid = getVectorWithMaxDistance(unassignedVectors,
					referenceCentroid);
			
			//remove the new centroid from unassignedVectors list
			unassignedVectors.remove(newCentroid);
			
			//add the new centroid to centroidsWitSeed list
			centroidsWithSeed.add(newCentroid);
		}
		
		//remove the original seed centoid, and get the k centroids
		//as double arrays to return
		centroidsWithSeed.remove(seedCentroid);
		centroidsWithSeed.toArray(centroids);
		return centroids;
	}
	
	/**
	 * @param unassignedVectors
	 * @return
	 */
	private static double[] getVectorWithMaxDistance(ArrayList unassignedVectors,
			double[] referenceCentroid) {
		double[] r = null;
		double maxDistance = Double.MIN_VALUE;
		double[][] vectors = new double[unassignedVectors.size()][];
		unassignedVectors.toArray(vectors);
		for(int i = 0; i < vectors.length; i++) {
			double thisDistance = getDistanceBetweenVectors(referenceCentroid,
					vectors[i]);
			if(thisDistance > maxDistance) {
				maxDistance = thisDistance;
				r = vectors[i];
			}
		}
		return r;
	}
	
	private static double[][] getRandomInitialCentroids(double[][] vectors, int k) {
		double[][] centroids = new double[k][];
		
		int vectorCount = vectors.length;
		int vectorLength = vectors[0].length;
		
		Random random = new Random();
		//for now, select centroids randomly from the input vectors
		boolean[] assigned = new boolean[vectorCount];
		for(int i = 0; i < k;) {
			int index = random.nextInt(vectorCount);
			if(!assigned[index]) {
				centroids[i] = vectors[index];
				assigned[index] = true;
				i++;
			}
		}
		
		return centroids;
	}
	
	/**
	 * @param clusters
	 * @param centroids
	 */
	private static void printClustersAndCentroids(double[][][] clusters, double[][] centroids) {
		System.out.println("clusters:");
		for(int i = 0; i < clusters.length; i++) {
			System.out.println("cluster " + i + ", size " + clusters[i].length);
			String centroidString = "centroid: ";
			for(int j = 0; j < centroids[i].length; j++) {
				centroidString = centroidString + "\t" + centroids[i][j];
			}
			System.out.println(centroidString);
			for(int j = 0; j < clusters[i].length; j++) {
				String s = "";
				for(int l = 0; l < clusters[i][j].length; l++) {
					s = s + "\t" + clusters[i][j][l];
				}
				System.out.println(s);
			}
			System.out.println();
		}	
		}
	
	private static void printCentroids(double[][] centroids) {
		for(int i = 0; i < centroids.length; i++) {
			String centroidString = "centroid " + i + ": ";
			for(int j = 0; j < centroids[i].length; j++) {
				centroidString = centroidString + "\t" + centroids[i][j];
			}
			System.out.println(centroidString);
		}
	}
	/**
	 * @param clusters
	 */
	private static void printClusters(double[][][] clusters) {
		System.out.println("clusters:");
		for(int i = 0; i < clusters.length; i++) {
			System.out.println("cluster " + i + ", size " + clusters[i].length);
			for(int j = 0; j < clusters[i].length; j++) {
				String s = "";
				for(int l = 0; l < clusters[i][j].length; l++) {
					s = s + "\t" + clusters[i][j][l];
				}
				System.out.println(s);
			}
			System.out.println();
		}
	}
	
	/**
	 * Clusters the input vectors according to the given assignemnts.
	 * Returned value is a three-dimension array. First dimension is
	 * cluster. Each cluster has a group of vectors. Third dimension is
	 * vector.
	 * @param vectors
	 * @param assignments
	 * @return
	 */
	public static double[][][] getClusters(double[][] vectors, int[] assignments) {
		double[][][] r = null;
		
		//determine the number of clusters. This should be the maximum value
		//in assignments +1
		int maxAssignment = -1;
		for(int i = 0; i < assignments.length; i++) {
			if(assignments[i] > maxAssignment) {
				maxAssignment = assignments[i];
			}
		}
		
		//determine sizes of clusters
		int[] clusterSizes = new int[maxAssignment+1];
		for(int i = 0; i < assignments.length; i++) {
			clusterSizes[assignments[i]]++;
		}
		
		r = new double[maxAssignment+1][][];		
		
		//initialize clusters to correct sizes
		for(int i = 0; i < r.length; i++) {
			r[i] = new double[clusterSizes[i]][];
		}
		
		//assign vectors to the correct clusters
		for(int i = 0; i < assignments.length; i++) {
			r[assignments[i]][--clusterSizes[assignments[i]]] = vectors[i];
		}
		
		return r;
	}
	
	/**
	 * Returns an array with the z scores for entries
	 * in d. NaN values are ignored in the calculations.
	 * Any index in d with NaN will also have NaN in r.
	 * @param d
	 * @return
	 */
	public static double[] getZScores(double[] d) {
		double[] r = null;
		
		//count non-NaN entries, and track them in a boolean array
		//to allow easier checking of this in subsequent code
		int nonNaN = 0;
		boolean[] hasValue= new boolean[d.length];
		for(int i = 0; i < d.length; i++) {
			//instead of using Double.isNaN(), avoid the method call
			//by taking advantage of the fact that NaN != NaN
			if(d[i] == d[i]) {
				nonNaN++;
				hasValue[i] = true;
			}
		}
		
		//determine the mean of all non-NaN values
		double sum = 0.0;
		for(int i = 0; i < d.length; i++) {
			if(hasValue[i]) {
				sum += d[i];
			}
		}
		
		double mean = sum / nonNaN;
		//System.out.println("mean: " + mean);
		//determine the standard deviation from the mean
		//of all non-NaN values
		double sos = 0.0;
		for(int i = 0; i < d.length; i++) {
			if(hasValue[i]) {
				double diff = mean - d[i];  
				sos += diff * diff;
			}
		}
		
		double sd = Math.sqrt((sos/(nonNaN - 1)));
		//System.out.println("sd: " + sd);
		//determine the z-score of all non-NaN values
		r = new double[d.length];
		for(int i = 0; i < d.length; i++) {
			if(hasValue[i]) {
				r[i] = (d[i] - mean) / sd;
			} else {
				r[i] = Double.NaN;
			}
		}
		
		return r;
	}
	
	public static double getStandardDeviation(double[] d, double mean) {
		double sos = 0.0;
		for(int i = 0; i < d.length; i++) {
			double diff = mean - d[i];  
			sos += diff * diff;
		}
		
		double sd = Math.sqrt((sos/(d.length - 1)));
		return sd;
	}

	public static double getCovariance(double[] d1, double[] d2) {
		double r = 0;
		if(d1.length == d2.length) {
			int n = d1.length;
			int denominator = 0;
			double mean1 = getMean(d1);
			double mean2 = getMean(d2);

			//only calculare covariance if both means are numbers.
			//That is, if either one is NaN, return 0 as covariance.
			//0 may not be the best value for this (?)
			if(!Double.isNaN(mean1) && !Double.isNaN(mean2)) {
				double s = 0;
				double diff1, diff2;
				for(int i = 0; i < n; i++) {
					diff1 = d1[i] - mean1;
					diff2 = d2[i] - mean2;
					double product = diff1 * diff2;
					//only use this element if neither value is NaN
					//or Infinity
					if(!Double.isNaN(product) 
							&& !Double.isInfinite(product)) {
						denominator++;
						s += product;
					}
				}

				if(denominator > 1 
						&& !Double.isNaN(denominator)
						&& !Double.isInfinite(denominator)) {
					r = s/(denominator-1);
				}
			}
		}
		return r;
	}

	public static double getCovariance(double[][] d1, double[][] d2) {
		double r = Double.NaN;
		if(d1.length == d2.length) {
			int n = d1.length;
			double[] mean1 = getMean(d1);
			double[] mean2 = getMean(d2);
			
			double s = 0;
			double diff1, diff2;
			for(int i = 0; i < n; i++) {
				diff1 = (d1[i][0] - mean1[0]) + (d1[i][1] - mean1[1]);
				diff2 = (d2[i][0] - mean2[0]) + (d2[i][1] - mean2[1]);
				s += (diff1 * diff2);
			}
			
			r = s/(n-1);
		}
		return r;
	}
	
	/**
	 * Returns the value that is min + (max - min)/2.
	 * That is, the value half way between the minimum
	 * and the maximum value in the given array.
	 *
	 * @param a
	 * @return
	 */
	public static int getMidpoint(int[] a) {
		int r = 0;
		int min = Integer.MAX_VALUE;
		int max = Integer.MIN_VALUE;
		
		for(int i = 0 ; i < a.length; i++) {
			if(a[i] > max) {
				max = a[i];
			}
			if(a[i] < min) {
				min = a[i];
			}
		}
		r = min + ((max - min)/2);
		return r;
	}
	
	/**
	 * Returns the maximum value in the array.
	 */
	public static int getMax(int[] a) {
		int max = Integer.MIN_VALUE;
		for(int i = 0; i < a.length; i++) {
			if(a[i] > max) {
				max = a[i];
			}
		}
		
		return max;
	}
	
	public static double getMax(double[] a) {
		double max = Double.NaN;
		if(a != null && a.length > 0) {
			max = a[0];
		for(int i = 0; i < a.length; i++) {
			if(a[i] > max) {
				max = a[i];
			}
		}
		}
		return max;
	}
	
	public static double[] floatsToDoubles(float[] f) {
		double[] d = new double[f.length];
		
		for(int i = 0; i < f.length; i++) {
			d[i] = f[i];
		}
		return d;
	}
	
	public static double getTTestPValue(double[] d1, double[] d2) {
		double pValue = Double.NaN;
//		TTest tTest = new TTestImpl();
//		try {
//			pValue = tTest.tTest(d1, d2);
//		} catch (IllegalArgumentException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		} catch (MathException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}

		return pValue;
	}
	
	public static double getTTestPValue(double mu, double[] d) {
		double pValue = Double.NaN;
//		TTest tTest = new TTestImpl();
//		try {
//			pValue = tTest.tTest(mu, d);
//		} catch (IllegalArgumentException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		} catch (MathException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}

		return pValue;
	}
	
	/**
	 * Only works for non-negative values.
	 * 
	 * @param a
	 * @return
	 */
	public static int[] getDistribution(int[] a) {
		int[] r = null;
		int max = 0;
		for(int i = 0; i < a.length; i++) {
			max = Math.max(a[i], max);
		}
		
		r = new int[max+1];
		for(int i = 0; i < a.length; i++) {
			r[a[i]]++;
		}

		return r;
	}
	
	/**
	 * Only works for non-negative values.
	 * @param a
	 */
	public static void printDistribution(int[] a) {
		int max = 0;
		for(int i = 0; i < a.length; i++) {
			max = Math.max(a[i], max);
		}
		
		int[] counts = new int[max+1];
		for(int i = 0; i < a.length; i++) {
			counts[a[i]]++;
		}
		
		for(int i = 0; i < counts.length; i++) {
			if(counts[i] > 0) {
				System.out.println("\t" + i + ": " + counts[i]);
			}
		}
	}

}
