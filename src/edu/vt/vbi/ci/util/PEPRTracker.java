package edu.vt.vbi.ci.util;

import java.util.ArrayList;

import edu.vt.vbi.ci.util.file.FastaSequenceFile;

/**
 * This class is for storing and printing out summary
 * information about a run of PEPR.
 * 
 * @author enordber
 *
 */
public class PEPRTracker {
	
	private static class TreeReport {
		private String runName;
		private String currentTreeString; //the tree built this round - either the initial full tree, or a refined subtree
		private String fullTreeString; //the current 'best', or most refined, tree. If this is a refining round, the refined subtree is included
		private long startTimeMillis;
		private long endTimeMillis;
		private int taxa;
		private int genes;
		private int alignedPositions;
		private double treeLength;
		private double fullTreeLength;
		
		private FastaSequenceFile[] inputSequenceFiles;
		private FastaSequenceFile[] outgroupSequenceFiles;

		private TreeReport() {
			startTimeMillis = System.currentTimeMillis();
		}
		
		public void setRunName(String runName) {
			this.runName = runName;
			endTimeMillis = System.currentTimeMillis();
		}
		
		public void setTree(String tree) {
			currentTreeString = tree;
			endTimeMillis = System.currentTimeMillis();
		}
		
		public void setTreeLength(double length) {
			treeLength = length;
		}
		
		public void setFullTree(String tree) {
			fullTreeString = tree;
			endTimeMillis = System.currentTimeMillis();
		}

		public void setFullTreeLength(double length) {
			fullTreeLength = length;
		}
		
		public void setInputSequenceFiles(FastaSequenceFile[] inputSequenceFiles) {
			this.inputSequenceFiles = inputSequenceFiles;
			endTimeMillis = System.currentTimeMillis();
		}

		public void setOutgroupPoolSequenceFiles(FastaSequenceFile[] outgroupSequenceFiles) {
			this.outgroupSequenceFiles = outgroupSequenceFiles;
		}

		public void setTaxonNumber(int taxa) {
			this.taxa = taxa;
		}

		public void setGeneNumber(int genes) {
			this.genes = genes;
		}

		public void setAlignedPositions(int alignedPositions) {
			this.alignedPositions = alignedPositions;
		}
	}
	
	
	private static long startTimeMillis = System.currentTimeMillis();
	private static long endTimeMillis;
	private static int totalGenomeCount;
	private static TreeReport openReport;
	private static ArrayList<TreeReport> treeReports;
	
	public static void newTree(String runName) {
		if(treeReports == null) {
			treeReports = new ArrayList<TreeReport>();
		}
		openReport = new TreeReport();
		treeReports.add(openReport);
	}


	public static void setRunName(String runName) {
		openReport.setRunName(runName);
	}
	
	public static void setTree(String tree) {
		openReport.setTree(tree);
	}

	public void setTreeLength(double length) {
		openReport.setTreeLength(length);
	}
	
	public void setFullTree(String tree) {
		openReport.setFullTree(tree);
	}

	public void setFullTreeLength(double length) {
		openReport.setFullTreeLength(length);
	}
	
	public static void setInputSequenceFiles(FastaSequenceFile[] inputSequenceFiles) {
		openReport.setInputSequenceFiles(inputSequenceFiles);
	}

	public static void setOutgroupPoolSequenceFiles(FastaSequenceFile[] outgroupSequenceFiles) {
		openReport.setOutgroupPoolSequenceFiles(outgroupSequenceFiles);
	}

	public static void setTaxonNumber(int taxa) {
		openReport.setTaxonNumber(taxa);
	}

	public static void setGeneNumber(int genes) {
		openReport.setGeneNumber(genes);
	}

	public static void setAlignedPositions(int alignedPositions) {
		openReport.setAlignedPositions(alignedPositions);
	}
	
	public String getReport() {
		String r = null;
		StringBuffer sb = new StringBuffer();
		r = sb.toString();
		return r;
	}

}
