package edu.vt.vbi.ci.pepr.tree;

import java.util.Arrays;
import java.util.HashSet;

import edu.vt.vbi.ci.util.ExtendedBitSet;

/**
 * Construct phylogenetic trees from distance matrices
 * or similarity matrices.
 * 
 * First version will only construct topology, without
 * branch lengths, using neighbor-joining algorithm.
 * 
 * @author enordber
 *
 */
public class TreeBuilder {

	public static final int DISTANCE = 0;
	public static final int SIMILARITY = 1;
	private static final boolean debug = false;

	/**
	 * Returns a String representation of the neighbor-joining 
	 * tree of the given nodes, based on the given distance matrix.
	 * Tree is returned in newick format.
	 * @param nodes
	 * @param matrix
	 * @param matrixType Type of matrix being provided.
	 * 	                  DISTANCE or SIMILARITY
	 * @return
	 */
	public static String getNJTopology(String[] nodes, 
			double[][] matrix, int matrixType) {
		String r = null;

		//find pair of nodes with minimum distance

		/*
		 * indices for pair of nodes with minimum distance
		 */
		int mergeI = -1; 
		int mergeJ = -1; 

		int[] mergeCoordinates = null;
		if(matrixType == DISTANCE) {
			mergeCoordinates = getMinValueCoordinates(matrix);
		} else if(matrixType == SIMILARITY) {
			mergeCoordinates = getMaxValueCoordinates(matrix);
		} else {
			throw new IllegalArgumentException("matrixType must be " +
					"either TreeBuider.DISTANCE or " +
					"TreeBuilder.SIMILARITY in " +
			"TreeBuilder.getNJTopology()");
		}
		mergeI = mergeCoordinates[0];
		mergeJ = mergeCoordinates[1];

		if(debug) {
			System.out.println("merging nodes " + mergeI + " and " + mergeJ
					+ ": " + nodes[mergeI] + ", " + nodes[mergeJ] + 
					" distance: " + matrix[mergeI][mergeJ]);
		}
		//create new node, consisting of the two closest nodes.
		StringBuffer combinedNode = new StringBuffer();
		combinedNode.append("(");
		combinedNode.append(nodes[mergeI]);
		combinedNode.append(",");
		combinedNode.append(nodes[mergeJ]);
		combinedNode.append(")");

		//if more than one node is left, make a recursive call to 
		//this method
		if(nodes.length > 2) { //Actually, this check could be done
			//before searching for minimum distance.
			//If there are only two nodes remaining,
			//just combine them and return the result.

			//create new distance matrix, using mean distances for the new node
			//New distance matrix should be one column and one row smaller
			//than input matrix, because two nodes have been combined into
			//a single node.
			double[][] newDistances = 
				new double[nodes.length-1][nodes.length-1];
			//Determine new index for each original node.
			//The new node will be at index 0
			int[] newIndices = new int[nodes.length];
			int newIndexOffset = 1;
			for(int i = 0; i < nodes.length; i++) {
				if(i == mergeI || i == mergeJ) {
					newIndexOffset--;
					newIndices[i] = 0;
				} else {
					newIndices[i] = i+newIndexOffset;
				}
			}

			//copy unchanged values to new distances
			for(int i = 0; i < matrix.length; i++) {
				if(i != mergeI) {
					int newI = newIndices[i];
					for(int j = i; j < matrix.length; j++) {
						if(j != mergeJ) {
							int newJ = newIndices[j];
							newDistances[newI][newJ] = 
								matrix[i][j];
							newDistances[newJ][newI] = 
								matrix[i][j];
						}
					}
				}
			}

			//determine average distances for new combined node
			for(int i = 0; i < matrix.length; i++) {
				if(newIndices[i] != 0) {
					double valueI = matrix[mergeI][i];
					double valueJ = matrix[mergeJ][i];
					double mean = (valueI + valueJ)/2;
					newDistances[0][newIndices[i]] = mean;
					newDistances[newIndices[i]][0] = mean;
				}
			}

			String[] newNodes = new String[nodes.length-1];
			for(int i = 0; i < nodes.length; i++){
				newNodes[newIndices[i]] = nodes[i];
			}
			newNodes[0] = combinedNode.toString();

			r = getNJTopology(newNodes, newDistances, matrixType);
		}	else {	
			//if only one node remains, return the String 
			//representation of the node
			combinedNode.append(";");
			r = combinedNode.toString();
		}
		return r;
	}

	/**
	 * Returns a String representation of the neighbor-joining 
	 * tree of the given nodes, based on the given distance matrix.
	 * Tree is returned in newick format.
	 * @param nodes
	 * @param distanceMatrix
	 * @param matrixType Type of matrix being provided.
	 * 	                  DISTANCE or SIMILARITY
	 * @return
	 */
	public static String getNJTreeString(String[] nodes, 
			double[][] distanceMatrix, int matrixType) {
		String r = null;

		//if similarity matrix is passed in, convert it to a distance matrix
		if(matrixType == SIMILARITY) {
			distanceMatrix = convertSimilarityToDistanceMatrix(distanceMatrix);
			matrixType = DISTANCE;
		}
		//create Q matrix from distance Matrix
		//find pair of nodes with minimum distance/maximum similarity
		double[][] Q = new double[distanceMatrix.length][distanceMatrix.length];

		//calculate sums for each row
		double[] sums = new double[Q.length];
		for(int i = 0; i < sums.length; i++) {
			sums[i] = 0;
			for(int j = 0; j < Q.length; j++) {
				sums[i] += distanceMatrix[i][j];
			}
		}

		int rMinus2 = Q.length - 2;

		for(int i = 0; i < Q.length; i++) {
			//initialize diagonals to 0.0
			Q[i][i] = 0.0;
			for(int j = i+1; j < Q.length; j++) {
				double value = (rMinus2 * distanceMatrix[i][j]) -
				sums[i] - sums[j];
				Q[i][j] = value;
				Q[j][i] = value;
			}
		}

		if(debug) {
			//print distance matrix and Q matrix
			printMatrices(nodes, distanceMatrix, Q);
		}

		/*
		 * indices for pair of nodes with minimum distance
		 */
		int mergeI = -1; 
		int mergeJ = -1; 

		int[] mergeCoordinates = null;
		if(matrixType == DISTANCE) {
			mergeCoordinates = getMinValueCoordinates(Q);
		} else if(matrixType == SIMILARITY) {
			mergeCoordinates = getMaxValueCoordinates(Q);
		} else {
			throw new IllegalArgumentException("matrixType must be " +
					"either TreeBuider.DISTANCE or " +
					"TreeBuilder.SIMILARITY in " +
			"TreeBuilder.getNJTopology()");
		}
		mergeI = mergeCoordinates[0];
		mergeJ = mergeCoordinates[1];

		if(debug) {
			System.out.println("merging nodes " + mergeI + " and " + mergeJ
					+ ": " + nodes[mergeI] + ", " + nodes[mergeJ] + 
					" distance: " + distanceMatrix[mergeI][mergeJ] + " Q: " + Q[mergeI][mergeJ]);
		}

		//determine distances from the new node to each of the merged nodes

		//from mergeI to combined
		double mergeIToCombinedDistance = (distanceMatrix[mergeI][mergeJ] /2) + 
		( (1/(2d*rMinus2)) * (sums[mergeI] - sums[mergeJ]) );

		//from mergeJ to combined
		double mergeJToCombined = distanceMatrix[mergeI][mergeJ] - mergeIToCombinedDistance;

		//this is a kludge - prevent negative branch lengths
		mergeIToCombinedDistance = Math.max(mergeIToCombinedDistance, 0);
		mergeJToCombined = Math.max(mergeJToCombined, 0);
		
		//create new node, consisting of the two closest nodes.
		StringBuffer combinedNode = new StringBuffer();
		combinedNode.append("(");
		combinedNode.append(nodes[mergeI]);
		combinedNode.append(":" + mergeIToCombinedDistance);
		combinedNode.append(",");
		combinedNode.append(nodes[mergeJ]);
		combinedNode.append(":" + mergeJToCombined);
		combinedNode.append(")");

		//if more than one node is left, make a recursive call to 
		//this method
		if(nodes.length > 0) { 

			//create new distance matrix, using mean distances for the new node
			//New distance matrix should be one column and one row smaller
			//than input matrix, because two nodes have been combined into
			//a single node.
			double[][] newDistances = 
				new double[nodes.length-1][nodes.length-1];
			//Determine new index for each original node.
			//The new node will be at index 0
			int[] newIndices = new int[nodes.length];
			int newIndexOffset = 1;
			for(int i = 0; i < nodes.length; i++) {
				if(i == mergeI || i == mergeJ) {
					newIndexOffset--;
					newIndices[i] = 0;
				} else {
					newIndices[i] = i+newIndexOffset;
				}
			}

			//copy unchanged values to new distances
			for(int i = 0; i < distanceMatrix.length; i++) {
				if(i != mergeI) {
					int newI = newIndices[i];
					for(int j = i; j < distanceMatrix.length; j++) {
						if(j != mergeJ) {
							int newJ = newIndices[j];
							newDistances[newI][newJ] = 
								distanceMatrix[i][j];
							newDistances[newJ][newI] = 
								distanceMatrix[i][j];
						}
					}
				}
			}

			//determine average distances for new combined node
			for(int i = 0; i < distanceMatrix.length; i++) {
				if(newIndices[i] != 0) {
					double valueI = distanceMatrix[mergeI][i];
					double valueJ = distanceMatrix[mergeJ][i];
					double mean = (valueI + valueJ - distanceMatrix[mergeI][mergeJ])/2;
					newDistances[0][newIndices[i]] = mean;
					newDistances[newIndices[i]][0] = mean;
				}
			}

			String[] newNodes = new String[nodes.length-1];
			for(int i = 0; i < nodes.length; i++){
				newNodes[newIndices[i]] = nodes[i];
			}
			newNodes[0] = combinedNode.toString();


			if(newNodes.length > 3) {
			    //set the distanceMatrix to null. Otherwise, if the recursion goes too
				//deep, we'll use too much memory
				distanceMatrix = null;
				r = getNJTreeString(newNodes, newDistances, matrixType);
			} else if(newNodes.length == 3) {
				//there are three nodes remaining, which means there is no more
				//merging to be done. 
				//Need to determine the branch lengths for the three remaining 
				//branches (from each of the three remaining nodes to the 'root')

				if(debug) {
					System.out.println("done merging. these are the 3 nodes remaining: "); 
					printMatrices(newNodes, newDistances, null);
					for(int i = 0; i < newNodes.length; i++) {
						System.out.println(newNodes[i]);
					}
				}

				sums = new double[newDistances.length];
				for(int i = 0; i < newDistances.length; i++) {
					sums[i] = 0;
					for(int j = 0; j < newDistances[i].length; j++) {
						sums[i] += newDistances[i][j];
					}
				}

				//determine distances for three remaining nodes
				double distance0 = (newDistances[0][1] /2) + 
				(0.5 * (sums[0] - sums[1]));
				double distance1 = (newDistances[1][2] /2) + 
				(0.5 * (sums[1] - sums[2]));
				double distance2 = (newDistances[2][0] /2) + 
				(0.5 * (sums[2] - sums[0]));

				if(debug) {
					System.out.println("distance0: " + distance0);
					System.out.println("distance1: " + distance1);
					System.out.println("distance2: " + distance2);
				}

				r = "(" + newNodes[0] + ":" + distance0 + ","
				+ newNodes[1] + ":" + distance1 + "," + newNodes[2] + ":" + distance2 + ");";
			}
		}
		return r;
	}

	private static double[][] convertSimilarityToDistanceMatrix(double[][] similarityMatrix) {
		double[][] r = new double[similarityMatrix.length][similarityMatrix.length];
		//find maximum similarity and subtract everything from that to get distances
		//assume maximum similarity is on the diagional, because it will be a self-similarity score
		double maxSimilarity = 0;
		for(int i = 0; i < similarityMatrix.length; i++){
			if(similarityMatrix[i][i] > maxSimilarity) {
				maxSimilarity = similarityMatrix[i][i];
			}
		}
		for(int i = 0; i < r.length; i++){
			for(int j = 0; j < r.length; j++) {
				r[i][j] = maxSimilarity - similarityMatrix[i][j];
			}
		}
		return r;
	}

	/**
	 * Adds branch lengths to the given tree topology using neighbor joining
	 *  based on the given distance matrix.
	 * Tree is returned in newick format.
	 * @param nodes
	 * @param distanceMatrix
	 * @return
	 */
	public static String addNJBranchLengths(String topologyTreeString, String[] nodes, 
			double[][] distanceMatrix) {
		String r = null;
//		System.out.println("TreeBuilder.addNJBranchLengths() (public)");
//		System.out.println("\tnodes.length: " + nodes.length);
//		System.out.println("\tdistanceMatrix.length: " + distanceMatrix.length);

		AdvancedTree topologyTree = new AdvancedTree(topologyTreeString);
		
		String[] sortedNodes = new String[nodes.length];
		System.arraycopy(nodes, 0, sortedNodes, 0, nodes.length);
		Arrays.sort(sortedNodes);

		Bipartition[] topologyBiparts = topologyTree.getBipartitions(sortedNodes);
		HashSet<ExtendedBitSet> topologyBipartsSet = new HashSet();
		for(int i = 0; i < topologyBiparts.length; i++) {
			topologyBipartsSet.add(topologyBiparts[i].getSmallerSide());
		}
		
		ExtendedBitSet[] nodeBiparts = new ExtendedBitSet[nodes.length];
		
		for(int i = 0; i < nodeBiparts.length; i++) {
			int nodeIndex = Arrays.binarySearch(sortedNodes, nodes[i]);
			if(nodeIndex > -1) {
				ExtendedBitSet bitSet = new ExtendedBitSet();
				bitSet.set(nodeIndex);
				nodeBiparts[i] = bitSet;
			}
		}
		
		r = addNJBranchLengths(topologyBipartsSet, nodeBiparts, nodes, distanceMatrix);
		return r;
	}
	
	private static String addNJBranchLengths(HashSet topologyBiparts, ExtendedBitSet[] nodeBiparts, String[] nodes, double[][] distanceMatrix) {
		String r = null;
//		if(nodes.length != distanceMatrix.length) {
//			System.out.println("TreeBuilder.addNJBranchLengths. modes.length != distanceMatrix.length:");
//			System.out.println("\tnodes.length: " + nodes.length);
//			System.out.println("\tdistanceMatrix.length: " + distanceMatrix.length);
//		} else {
//			System.out.println("TreeBuilder.addNJBranchLengths. modes.length == distanceMatrix.length.  It's cool");
//		}
		//create Q matrix from distance Matrix
		//find pair of nodes with minimum distance
		double[][] Q = new double[distanceMatrix.length][distanceMatrix.length];

		//calculate sums for each row
		double[] sums = new double[Q.length];
		for(int i = 0; i < sums.length; i++) {
			sums[i] = 0;
			for(int j = 0; j < Q.length; j++) {
				sums[i] += distanceMatrix[i][j];
			}
		}

		int rMinus2 = Q.length - 2;

		for(int i = 0; i < Q.length; i++) {
			//initialize diagonals to 0.0
			Q[i][i] = 0.0;
			for(int j = i+1; j < Q.length; j++) {
				double value = (rMinus2 * distanceMatrix[i][j]) -
				sums[i] - sums[j];
				Q[i][j] = value;
				Q[j][i] = value;
			}
		}

		if(debug) {
			//print distance matrix and Q matrix
			printMatrices(nodes, distanceMatrix, Q);
		}

		/*
		 * indices for pair of nodes with minimum distance
		 */
		int mergeI = -1; 
		int mergeJ = -1; 

		int[] mergeCoordinates = getMergeIndices(topologyBiparts, nodeBiparts, Q); //getMinValueCoordinates(Q);
		mergeI = mergeCoordinates[0];
		mergeJ = mergeCoordinates[1];

		if(debug) {
			System.out.println("merging nodes " + mergeI + " and " + mergeJ
					+ ": " + nodes[mergeI] + ", " + nodes[mergeJ] + 
					" distance: " + distanceMatrix[mergeI][mergeJ] + " Q: " + Q[mergeI][mergeJ]);
		}

		//determine distances from the new node to each of the merged nodes

		//from mergeI to combined
		double mergeIToCombinedDistance = (distanceMatrix[mergeI][mergeJ] /2) + 
		( (1/(2d*rMinus2)) * (sums[mergeI] - sums[mergeJ]) );

		//from mergeJ to combined
		double mergeJToCombined = distanceMatrix[mergeI][mergeJ] - mergeIToCombinedDistance;

		//this is a kludge - prevent negative branch lengths
		mergeIToCombinedDistance = Math.max(mergeIToCombinedDistance, 0);
		mergeJToCombined = Math.max(mergeJToCombined, 0);
		
		//create new node, consisting of the two closest nodes.
		StringBuffer combinedNode = new StringBuffer();
		combinedNode.append("(");
		combinedNode.append(nodes[mergeI]);
		combinedNode.append(":" + mergeIToCombinedDistance);
		combinedNode.append(",");
		combinedNode.append(nodes[mergeJ]);
		combinedNode.append(":" + mergeJToCombined);
		combinedNode.append(")");

		//if more than one node is left, make a recursive call to 
		//this method
		if(nodes.length > 0) { 

			//create new distance matrix, using mean distances for the new node
			//New distance matrix should be one column and one row smaller
			//than input matrix, because two nodes have been combined into
			//a single node.
			double[][] newDistances = 
				new double[nodes.length-1][nodes.length-1];
			//Determine new index for each original node.
			//The new node will be at index 0
			int[] newIndices = new int[nodes.length];
			int newIndexOffset = 1;
			for(int i = 0; i < nodes.length; i++) {
				if(i == mergeI || i == mergeJ) {
					newIndexOffset--;
					newIndices[i] = 0;
				} else {
					newIndices[i] = i+newIndexOffset;
				}
			}

			//copy unchanged values to new distances
			for(int i = 0; i < distanceMatrix.length; i++) {
				if(i != mergeI) {
					int newI = newIndices[i];
					for(int j = i; j < distanceMatrix.length; j++) {
						if(j != mergeJ) {
							if(j >= newIndices.length) {
								System.out.println("\tj >= newIndices.length");
								System.out.println("\tj: " + j );
								System.out.println("\tnewIndices.length: " + newIndices.length);
								System.out.println("\tnodes.length: " + nodes.length);
								System.out.println("\tdistanceMatrix.length: " + distanceMatrix.length);
							}
							int newJ = newIndices[j]; //there's an "out of bounds" happening here sometimes
							newDistances[newI][newJ] = 
								distanceMatrix[i][j];
							newDistances[newJ][newI] = 
								distanceMatrix[i][j];
						}
					}
				}
			}

			//determine average distances for new combined node
			for(int i = 0; i < distanceMatrix.length; i++) {
				if(newIndices[i] != 0) {
					double valueI = distanceMatrix[mergeI][i];
					double valueJ = distanceMatrix[mergeJ][i];
					double mean = (valueI + valueJ - distanceMatrix[mergeI][mergeJ])/2;
					newDistances[0][newIndices[i]] = mean;
					newDistances[newIndices[i]][0] = mean;
				}
			}

			String[] newNodes = new String[nodes.length-1];
			for(int i = 0; i < nodes.length; i++){
				newNodes[newIndices[i]] = nodes[i];
			}
			newNodes[0] = combinedNode.toString();
			
			ExtendedBitSet[] newNodeBiparts = new ExtendedBitSet[nodeBiparts.length-1];
			for(int i = 0; i < nodeBiparts.length; i++) {
				newNodeBiparts[newIndices[i]] = nodeBiparts[i];
			}
			newNodeBiparts[0] = nodeBiparts[mergeI].getOr(nodeBiparts[mergeJ]);


			if(newNodes.length > 3) {
			    //set the distanceMatrix, and other variables, to null. Otherwise, if the recursion goes too
				//deep, we'll use too much memory
				distanceMatrix = null;
				nodeBiparts = null;
				nodes = null;
				r = addNJBranchLengths(topologyBiparts, newNodeBiparts, newNodes, newDistances);
			} else if(newNodes.length == 3) {
				//there are three nodes remaining, which means there is no more
				//merging to be done. 
				//Need to determine the branch lengths for the three remaining 
				//branches (from each of the three remaining nodes to the 'root')

				if(debug) {
					System.out.println("done merging. these are the 3 nodes remaining: "); 
					printMatrices(newNodes, newDistances, null);
					for(int i = 0; i < newNodes.length; i++) {
						System.out.println(newNodes[i]);
					}
				}

				sums = new double[newDistances.length];
				for(int i = 0; i < newDistances.length; i++) {
					sums[i] = 0;
					for(int j = 0; j < newDistances[i].length; j++) {
						sums[i] += newDistances[i][j];
					}
				}

				//determine distances for three remaining nodes
				double distance0 = (newDistances[0][1] /2) + 
				(0.5 * (sums[0] - sums[1]));
				double distance1 = (newDistances[1][2] /2) + 
				(0.5 * (sums[1] - sums[2]));
				double distance2 = (newDistances[2][0] /2) + 
				(0.5 * (sums[2] - sums[0]));

				if(debug) {
					System.out.println("distance0: " + distance0);
					System.out.println("distance1: " + distance1);
					System.out.println("distance2: " + distance2);
				}

				r = "(" + newNodes[0] + ":" + distance0 + ","
				+ newNodes[1] + ":" + distance1 + "," + newNodes[2] + ":" + distance2 + ");";
			}
		}
		return r;
	}

	private static int[] getMergeIndices(HashSet topologyBiparts,
			ExtendedBitSet[] nodeBiparts, double[][] Q) {
		int[] r = null;
		double minPairDistanceFound = Double.MAX_VALUE;
		
		//find two entries in nodeBiparts that can be combined (union) to 
		//produce one of the entries in topologyBiparts
		for(int i = 0; i < nodeBiparts.length; i++) {
			ExtendedBitSet iSet = nodeBiparts[i];
			for(int j = i+1; j < nodeBiparts.length; j++) {
				//create the union bipartion and see if it is in the set of
				//topology bipartitions
				ExtendedBitSet unionSet = iSet.getOr(nodeBiparts[j]);
				if(topologyBiparts.contains(unionSet)) {
					double distance = Q[i][j];
					if(distance < minPairDistanceFound) {
						r = new int[]{i,j};
						minPairDistanceFound = distance;
//						System.out.println("merge pair: " + i + ", " + j + ": " + distance);
					}
				} 
			}
		}
		return r;
	}

	private static void printMatrices(String[] nodes,
			double[][] distanceMatrix, double[][] q) {


		if(nodes != null) {
			System.out.println(nodes.length + " nodes: ");
			for(int i = 0; i < nodes.length; i++) {
				System.out.print(nodes[i] + "\t");
			}
			System.out.println();
		}

		if(distanceMatrix != null) {
			System.out.println("distance matrix:");
			for(int i = 0; i < distanceMatrix.length; i++) {
				for(int j = 0; j < distanceMatrix.length; j++) {
					System.out.print(distanceMatrix[i][j] + "\t");
				}
				System.out.println();
			}

			System.out.println();
		}

		if(q != null) {
			System.out.println("Q matrix:");
			for(int i = 0; i < q.length; i++) {
				for(int j = 0; j < q.length; j++) {
					System.out.print(q[i][j] + "\t");
				}
				System.out.println();
			}
		}
	}

	private static double[][] getQMatrix(double[][] distanceMatrix) {
		double[][] Q = new double[distanceMatrix.length][distanceMatrix.length];

		//calculate sums for each row
		double[] sums = new double[Q.length];
		for(int i = 0; i < sums.length; i++) {
			sums[i] = 0;
			for(int j = 0; j < Q.length; j++) {
				sums[i] += distanceMatrix[i][j];
			}
		}

		int rMinus2 = Q.length - 2;

		for(int i = 0; i < Q.length; i++) {
			//initialize diagonals to 0.0
			Q[i][i] = 0.0;
			for(int j = i+1; j < Q.length; j++) {
				double value = (rMinus2 * distanceMatrix[i][j]) -
				sums[i] - sums[j];
				Q[i][j] = value;
				Q[j][i] = value;
			}
		}

		return Q;
	}

	/**
	 * Internal utility method for finding the values of
	 * i and j that have the minimum value in the given
	 * matrix. This currently assumes a square symmetrical
	 * matrix.
	 * 
	 * @return
	 */
	private static int[] getMinValueCoordinates(double[][] matrix) {
		int[] r = null;
		//find pair of indices with minimum value
		int minI = -1, minJ = -1; //indices for pair with minimum value
		double minDistance = Double.MAX_VALUE; //minimum value so far

		for(int i = 0; i < matrix.length; i++) {
			for(int j = i+1; j < matrix[i].length; j++) {
				if(matrix[i][j] < minDistance) {
					minDistance = matrix[i][j];
					minI = i;
					minJ = j;
				}
			}
		}

		r = new int[]{minI, minJ};
		return r;
	}

	/**
	 * Internal utility method for finding the values of
	 * i and j that have the minimum value in the given
	 * matrix. This currently assumes a square symmetrical 
	 * matrix.
	 * 
	 * @return
	 */
	private static int[] getMaxValueCoordinates(double[][] matrix) {
		int[] r = null;
		//find pair of indices with maximum value
		int maxI = -1, maxJ = -1; //indices for pair with maximum value
		double maxDistance = Double.NEGATIVE_INFINITY; //maximum value so far

		for(int i = 0; i < matrix.length; i++) {
			for(int j = i+1; j < matrix[i].length; j++) {
				if(matrix[i][j] > maxDistance) {
					maxDistance = matrix[i][j];
					maxI = i;
					maxJ = j;
				}
			}
		}

		r = new int[]{maxI, maxJ};
		return r;
	}
}
