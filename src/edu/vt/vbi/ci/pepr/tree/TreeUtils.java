package edu.vt.vbi.ci.pepr.tree;

import java.util.HashMap;

import edu.vt.vbi.ci.util.StringPair;


public class TreeUtils {

	public static void main(String[] args) {
		/*
		CommandLineProperties clp = new CommandLineProperties(args);
		String treeFile = clp.getValues(HandyConstants.TREE_FILE)[0];
		String[] taxaToRemove = clp.getValues(HandyConstants.REMOVE);
		String treeString;
		try {
			treeString = new TextFile(treeFile).toString().trim();
			String trimmedTree = removeTaxaFromTree(treeString, taxaToRemove);
			System.out.println(trimmedTree);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
        */
		
		String in = "Yersinia pestis KIM 10 +";
		String out = compressTaxonNameForComparison(in);
		System.out.println(in + " --> " + out);
	}
	
	public static int getRobinsonFouldsDistance(BasicTree t1, BasicTree t2) {
		int r = -1;

		return r;
	}

	public static void loadGenomePairDistancesFromTree(String treeString,
			HashMap<StringPair, double[]> pairToDistances) {
		if(treeString != null && pairToDistances != null) {
			BasicTree tree = new BasicTree(treeString);
			double[][] distanceMatrix = tree.getDistanceMatrix();
			String[] genomeNames = tree.getLeaves();

			for(int i = 0; i < genomeNames.length; i++) {
				for(int j = i+1; j < genomeNames.length; j++) {
					StringPair genomePair = new StringPair(genomeNames[i], genomeNames[j]);
					double pairDistance = distanceMatrix[i][j];
 					double[] pairDistances = pairToDistances.get(genomePair);
					if(pairDistances == null) {
						pairDistances = new double[]{pairDistance};
					} else {
						double[] newDistances = new double[pairDistances.length+1];
						System.arraycopy(pairDistances, 0, newDistances, 0, pairDistances.length);
						newDistances[pairDistances.length] = pairDistance;
						pairDistances = newDistances;
					}
					pairToDistances.put(genomePair, pairDistances);
				}
			}
		}
	}
	
	public static String removeTaxaFromTree(String treeString, String[] taxaToRemove) {
		String r = null;
		//get BasicTree object for tree
		BasicTree tree = new BasicTree(treeString);
		
		for(int i = 0; i < taxaToRemove.length; i++) {
			tree.removeTaxon(taxaToRemove[i]);
		}
		
		r = tree.getTreeString(true, true);
		
		return r;
	}

	/**
	 * Returns the index of q in s, or -1 if q is not in s.
	 * 
	 * @param s
	 * @param q
	 * @return
	 */
	private static int indexOf(String[] s, String q) {
		int r = -1;
		for(int i = 0; i < s.length; i++) {
			if(q.equals(s[i])) {
				r = i;
				break;
			}
		}
		return r;
	}

	/**
	 * Remove the row and column for the specified doomedIndex form 
	 * the distance matrix. distances should be square and symmetrical.
	 *  
	 * @param distances
	 * @param doomedIndex
	 * @return
	 */
	private static double[][] removeIndexFromMatrix(double[][] distances,
			int doomedIndex) {
		double[][] r = null;
		if(distances != null && distances.length > 0 
				&& distances.length == distances[0].length) {
			int n = distances.length;
			r = new double[n-1][];
			int copyIndex = 0;
			for(int i = 0; i < n; i++) {
				if(i != doomedIndex) {
					r[copyIndex++] = remove(distances[i], doomedIndex);
				}
			}
		}
		
	    return r;
	}

	private static double[] remove(double[] d, int doomedIndex) {
		double[] r = null;
		if(d != null && doomedIndex >=0 && d.length > doomedIndex) {
			r = new double[d.length -1];
			System.arraycopy(d, 0, r, 0, doomedIndex);
			System.arraycopy(d, doomedIndex+1, r, doomedIndex, r.length-doomedIndex);
		}
		return r;
	}

	private static String[] remove(String[] s, int doomedIndex) {
		String[] r = null;
		if(s != null && doomedIndex >=0 && s.length > doomedIndex) {
			r = new String[s.length -1];
			System.arraycopy(s, 0, r, 0, doomedIndex);
			System.arraycopy(s, doomedIndex+1, r, doomedIndex, r.length-doomedIndex);
		}
		
		return r;
	}

	/**
	 * Takes a taxon name in a variety of forms and tries to 
	 * convert it to a compressed form for comparisons. Any non-alphanumeric
	 * characters are removed, and everything is converted to lowercase. This
	 * is so taxon names coming from different sources 
	 * (eg sequence file name, phylogeny taxon, fasta title line taxon) can
	 * still be compared. There is a risk that two different taxon names could
	 * compress to the same String. This seems unlikely, but should be watched
	 * for by users of this method
	 * 
	 * @return
	 */
	public static String compressTaxonNameForComparison(String fullName) {
		String r = fullName;
		//if fullName looks like it has an extension at the end (ie it's a file name)
		//the remove the extension
		if(r.matches(".*\\.f.+")) {
			r = r.substring(0, r.lastIndexOf('.'));
		}

		//if it's a PATRIC file name, it may have '.PATRIC' at the end
		if(r.endsWith(".PATRIC")) {
			r = r.substring(0, r.lastIndexOf('.'));			
		}
		
		//remove non-alphanumerics --oops, this is a problem. It turns out that
		//some PATRIC genome names differ by only spaces and non-alphanumeric 
		//characters
		//r = r.replaceAll("[^A-Za-z0-9]", "");
		
		//remove period, underscore, comma, and space
		r = r.replaceAll("[\\._,\\s]", "");
		
		//turn everything lowercase
		r = r.toLowerCase();
		
		return r;
	}
	
	public static String[] compressTaxonNamesForComparison(String[] fullNames) {
		String[] r= null;
		if(fullNames != null) {
			r = new String[fullNames.length];
			
			for(int i = 0; i < r.length; i++) {
				r[i] = compressTaxonNameForComparison(fullNames[i]);
			}
		}
		
		return r;
	}


}
