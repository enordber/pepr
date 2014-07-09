package edu.vt.vbi.ci.pepr.tree;

import java.util.HashMap;

import edu.vt.vbi.ci.util.StringPair;


public class TreeUtils {
	
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
		//then remove the extension
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
