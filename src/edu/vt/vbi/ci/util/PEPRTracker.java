package edu.vt.vbi.ci.util;

import java.util.ArrayList;
import java.util.Date;
import java.util.regex.Pattern;

import edu.vt.vbi.ci.pepr.tree.BasicTree;
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
		private String fullTreeMethod;
		private String supportTreeMethod;
		private String fullTreeOptions = "";
		private String supportTreeOptions = "";

		private FastaSequenceFile[] inputSequenceFiles;
		private FastaSequenceFile[] outgroupSequenceFiles;

		private TreeReport() {
			startTimeMillis = System.currentTimeMillis();
		}

		private void setRunName(String runName) {
			this.runName = runName;
			endTimeMillis = System.currentTimeMillis();
		}

		private void setTree(String tree) {
			currentTreeString = tree;
			BasicTree basicTree = new BasicTree(tree);
//			double length = 0;
//			for(double branchLength: basicTree.getBranchLengths()) {
//				length += branchLength;
//			}
//			setTreeLength(length);
			endTimeMillis = System.currentTimeMillis();
		}

		private void setTreeLength(double length) {
			treeLength = length;
		}

		private double getTreeLength() {
			return treeLength;
		}

		private void setFullTree(String tree) {
			fullTreeString = tree;
			BasicTree basicTree = new BasicTree(tree);
//			double length = 0;
//			for(double branchLength: basicTree.getBranchLengths()) {
//				length += branchLength;
//			}
//			setFullTreeLength(length);
			endTimeMillis = System.currentTimeMillis();
		}

		private void setFullTreeLength(double length) {
			fullTreeLength = length;
		}

		private double getFullTreeLength() {
			return fullTreeLength;
		}

		private void setInputSequenceFiles(FastaSequenceFile[] inputSequenceFiles) {
			this.inputSequenceFiles = inputSequenceFiles;
			endTimeMillis = System.currentTimeMillis();
		}

		private void setOutgroupPoolSequenceFiles(FastaSequenceFile[] outgroupSequenceFiles) {
			this.outgroupSequenceFiles = outgroupSequenceFiles;
		}

		private void setTaxonCount(int taxa) {
			this.taxa = taxa;
		}

		private int getTaxonCount() {
			return taxa;
		}

		private void setGeneNumber(int genes) {
			this.genes = genes;
		}

		private int getGeneNumber() {
			return genes;
		}

		private void setAlignedPositions(int alignedPositions) {
			this.alignedPositions = alignedPositions;
		}

		private int getAlignedPositions() {
			return alignedPositions;
		}

		private long getStartTimeMillis() {
			return startTimeMillis;
		}

		private void setStartTimeMillis(long startTimeMillis) {
			this.startTimeMillis = startTimeMillis;
		}

		private long getEndTimeMillis() {
			return endTimeMillis;
		}

		private void setEndTimeMillis(long endTimeMillis) {
			this.endTimeMillis = endTimeMillis;
		}

		private FastaSequenceFile[] getInputSequenceFiles() {
			return inputSequenceFiles;
		}

		private FastaSequenceFile[] getOutgroupSequenceFiles() {
			return outgroupSequenceFiles;
		}

		private String getCurrentTreeString() {
			return currentTreeString;
		}

		private String getFullTreeString() {
			return fullTreeString;
		}

		private String getRunName() {
			return runName;
		}

		private String getFullTreeMethod() {
			return fullTreeMethod;
		}

		private void setFullTreeMethod(String fullTreeMethod) {
			this.fullTreeMethod = fullTreeMethod;
		}

		private String getSupportTreeMethod() {
			return supportTreeMethod;
		}

		private void setSupportTreeMethod(String supportTreeMethod) {
			this.supportTreeMethod = supportTreeMethod;
		}

		private void setTreeOptions(String treeOptions) {
			if(this.fullTreeOptions == null || this.fullTreeOptions.length() == 0) {
				this.fullTreeOptions = treeOptions;
			} else {
				this.supportTreeOptions = treeOptions;
			}
		}

		private String getFullTreeOptions() {
			return fullTreeOptions;
		}

		private String getSupportTreeOptions() {
			return supportTreeOptions;
		}
	}


	private static long startTimeMillis = System.currentTimeMillis();
	private static long endTimeMillis;
	private static int totalGenomeCount;
	private static String runName = null;
	private static TreeReport openReport;
	private static ArrayList<TreeReport> treeReports;

	public static void newTree(String name) {
		if(treeReports == null) {
			treeReports = new ArrayList<TreeReport>();
		}
		openReport = new TreeReport();
		openReport.setRunName(name);
		if(runName == null) {
			runName = name;
		}
		treeReports.add(openReport);
	}


	public static void setRunName(String name) {
		openReport.setRunName(name);
		if(runName == null) {
			runName = name;
		}
	}

	public static void setFullTreeMethod(String fullTreeMethod) {
		openReport.setFullTreeMethod(fullTreeMethod);
	}

	public static void setSupportTreeMethod(String supportTreeMethod) {
		openReport.setSupportTreeMethod(supportTreeMethod);
	}

	public static void setTree(String tree) {
		openReport.setTree(tree);
		endTimeMillis = System.currentTimeMillis();
	}

	public static void setTreeLength(double length) {
		openReport.setTreeLength(length);
		endTimeMillis = System.currentTimeMillis();
	}

	public static void setFullTree(String tree) {
		openReport.setFullTree(tree);
		endTimeMillis = System.currentTimeMillis();
	}

	public static void setFullTreeLength(double length) {
		openReport.setFullTreeLength(length);
		endTimeMillis = System.currentTimeMillis();
	}

	public static void setInputSequenceFiles(FastaSequenceFile[] inputSequenceFiles) {
		openReport.setInputSequenceFiles(inputSequenceFiles);
	}

	public static void setOutgroupPoolSequenceFiles(FastaSequenceFile[] outgroupSequenceFiles) {
		openReport.setOutgroupPoolSequenceFiles(outgroupSequenceFiles);
	}

	public static void setTaxonNumber(int taxa) {
		openReport.setTaxonCount(taxa);
	}

	public static void setGeneNumber(int genes) {
		openReport.setGeneNumber(genes);
	}

	public static void setAlignedPositions(int alignedPositions) {
		openReport.setAlignedPositions(alignedPositions);
	}

	public static void setTreeOptions(String treeOptions) {
		openReport.setTreeOptions(treeOptions);
	}

	public static String getReport() {
		String r = null;
		Pattern genomeNameDelimiter = Pattern.compile("_@_");
		StringBuffer sb = new StringBuffer();
		sb.append("<pepr_report>");
		sb.append("\n");

		long millisPerSecond = 1000;
		long millisPerMinute = millisPerSecond*60;
		long millisPerHour = millisPerMinute*60;
		long millisPerDay = millisPerHour*24;
		sb.append("<run_name>");
		sb.append(runName);
		sb.append("</run_name>");
		sb.append("\n");

		sb.append("<start_time>");
		sb.append(new Date(startTimeMillis));
		sb.append("</start_time>");
		sb.append("\n");

		sb.append("<end_time>");
		sb.append(new Date(endTimeMillis));
		sb.append("</end_time>");
		sb.append("\n");

		long elapsedMillis = endTimeMillis - startTimeMillis;
		long elapsedHours = elapsedMillis / millisPerHour;
		long elapsedMinutes = (elapsedMillis / millisPerMinute) % 60;
		long elapsedSeconds = (elapsedMillis / millisPerSecond) % 60;

		sb.append("<total_elapsed_time>");
		sb.append(elapsedHours + "h " + elapsedMinutes + "m " + elapsedSeconds + "s");
		sb.append("</total_elapsed_time>");
		sb.append("\n");
		sb.append("<total_elapsed_seconds>");
		sb.append(elapsedMillis/millisPerSecond);
		sb.append("</total_elapsed_seconds>");
		sb.append("\n");

		sb.append("<trees>");
		sb.append("\n");
		for(TreeReport report: treeReports) {
			long startTime = report.getStartTimeMillis();
			long endTime = report.getEndTimeMillis();
			elapsedMillis = endTime - startTime;
			elapsedHours = elapsedMillis / millisPerHour;
			elapsedMinutes = (elapsedMillis / millisPerMinute) % 60;
			elapsedSeconds = (elapsedMillis / millisPerSecond) % 60;

			sb.append("<tree>");
			sb.append("\n");

			sb.append("<run_name>");
			sb.append(report.getRunName());
			sb.append("</run_name>");
			sb.append("\n");

			sb.append("<start_time>");
			sb.append(new Date(startTime));
			sb.append("</start_time>");
			sb.append("\n");

			sb.append("<end_time>");
			sb.append(new Date(endTime));
			sb.append("</end_time>");
			sb.append("\n");

			sb.append("<elapsed_time>");
			sb.append(elapsedHours + "h " + elapsedMinutes + "m " + elapsedSeconds + "s");
			sb.append("</elapsed_time>");
			sb.append("\n");
			sb.append("<elapsed_seconds>");
			sb.append(elapsedMillis/millisPerSecond);
			sb.append("</elapsed_seconds>");
			sb.append("\n");

			sb.append("<full_tree_method>");
			sb.append(report.getFullTreeMethod());
			sb.append("</full_tree_method>");
			sb.append("\n");

			sb.append("<full_tree_options>");
			sb.append(report.getFullTreeOptions());
			sb.append("</full_tree_options>");
			sb.append("\n");

			sb.append("<support_tree_method>");
			sb.append(report.getSupportTreeMethod());
			sb.append("</support_tree_method>");
			sb.append("\n");

			sb.append("<support_tree_options>");
			sb.append(report.getSupportTreeOptions());
			sb.append("</support_tree_options>");
			sb.append("\n");

			sb.append("<taxon_count>");
			sb.append(report.getTaxonCount());
			sb.append("</taxon_count>");
			sb.append("\n");

			sb.append("<gene_count>");
			sb.append(report.getGeneNumber());
			sb.append("</gene_count>");
			sb.append("\n");

			sb.append("<aligned_positions>");
			sb.append(report.getAlignedPositions());
			sb.append("</aligned_positions>");
			sb.append("\n");

			sb.append("<ingroup>");
			sb.append("\n");
			FastaSequenceFile[] ingroupFiles = report.getInputSequenceFiles();
			sb.append("<count>");
			sb.append(ingroupFiles.length);
			sb.append("</count>");
			sb.append("\n");
			sb.append("<ingroup_genomes>");
			sb.append("\n");
			for(FastaSequenceFile file: ingroupFiles) {
				if(file.getTaxa().length > 0) {
					sb.append("<genome>");						
					sb.append("\n");						
					String[] nameFields = genomeNameDelimiter.split(file.getTaxa()[0]);
					sb.append("<genome_name>");
					sb.append(nameFields[0]);
					sb.append("</genome_name>");
					sb.append("\n");
					if(nameFields.length > 1) {
						sb.append("<genome_id>");
						sb.append(nameFields[1]);
						sb.append("</genome_id>");
						sb.append("\n");						
					}
					sb.append("</genome>");
					sb.append("\n");						
				}
			}
			sb.append("</ingroup_genomes>");
			sb.append("\n");
			sb.append("</ingroup>");
			sb.append("\n");

			sb.append("<outgroup>");
			sb.append("\n");
			FastaSequenceFile[] outgroupFiles = report.getOutgroupSequenceFiles();
			sb.append("<count>");
			sb.append(outgroupFiles.length);
			sb.append("</count>");
			sb.append("\n");
			sb.append("<outgroup_genomes>");
			sb.append("\n");
			for(FastaSequenceFile file: outgroupFiles) {
				if(file.getTaxa().length > 0) {
					sb.append("<genome>");						
					sb.append("\n");						
					String[] nameFields = genomeNameDelimiter.split(file.getTaxa()[0]);
					sb.append("<genome_name>");
					sb.append(nameFields[0]);
					sb.append("</genome_name>");
					sb.append("\n");
					if(nameFields.length > 1) {
						sb.append("<genome_id>");
						sb.append(nameFields[1]);
						sb.append("</genome_id>");
						sb.append("\n");						
					}
					sb.append("</genome>");
					sb.append("\n");						
				}
			}
			sb.append("</outgroup_genomes>");
			sb.append("\n");
			sb.append("</outgroup>");
			sb.append("\n");

			sb.append("<current_round_tree>");
			sb.append(report.getCurrentTreeString());
			sb.append("</current_round_tree>");
			sb.append("\n");

//			sb.append("<current_round_tree_length>");
//			sb.append(report.getTreeLength());
//			sb.append("</current_round_tree_length>");
//			sb.append("\n");

			sb.append("<current_full_tree>");
			sb.append(report.getFullTreeString());
			sb.append("</current_full_tree>");
			sb.append("\n");

//			sb.append("<current_full_tree_length>");
//			sb.append(report.getFullTreeLength());
//			sb.append("</current_full_tree_length>");
//			sb.append("\n");

			sb.append("</tree>");
			sb.append("\n");
		}

		sb.append("<final_tree>");
		if(treeReports.size() > 0) {
			sb.append(treeReports.get(treeReports.size()-1).getFullTreeString());
		} else {
			sb.append("No tree was produced. Check log file for possible error messages.");
		}
		sb.append("</final_tree>");
		sb.append("\n");

		sb.append("</trees>");
		sb.append("\n");
		sb.append("</pepr_report>");
		sb.append("\n");
		r = sb.toString();
		return r;
	}

}
