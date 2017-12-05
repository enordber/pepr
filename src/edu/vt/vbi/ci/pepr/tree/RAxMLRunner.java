
package edu.vt.vbi.ci.pepr.tree;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashSet;
import java.util.Random;

import org.apache.log4j.Logger;

import edu.vt.vbi.ci.pepr.alignment.SequenceAlignment;
import edu.vt.vbi.ci.pepr.alignment.SequenceAlignmentParser;
import edu.vt.vbi.ci.util.CommandResults;
import edu.vt.vbi.ci.util.ExecUtilities;
import edu.vt.vbi.ci.util.HandyConstants;
import edu.vt.vbi.ci.util.PEPRTracker;
import edu.vt.vbi.ci.util.file.TextFile;

/**
 * this class runs the maximum likelihood phylogeny program RAxML.
 * 
 * @author enordber
 *
 */
public class RAxMLRunner {
	private static final int ML_ALGORITHM = 0;
	private static final int PARSIMONY_ALGORITHM = 1;
	private static final int PARSIMONY_WITH_BL_ALGORITHM = 2;
	private static final int PER_SITE_LL_ALGORITHM = 3;
	
	private static final String RAXML_HPC = "raxmlHPC";
	private static final String RAXML = "raxml";
	private static Random random = new Random();
	private static HashSet usedRunNames = new HashSet();
	private Logger logger = Logger.getLogger(getClass());
	private SequenceAlignment alignment;
	private int bootstrapReps = 10;
	private String runName;
	private int threadCount;
	private boolean useTaxonNames = true;
	private boolean doPerSiteLikelihoods = false;
	private boolean deleteGeneratedFiles = true;
	private static int lastFileNumberUsed = 0;
	private int algorithm = ML_ALGORITHM;
	private String matrix = "PROTGAMMALGF";
	private String[] perSiteLLTrees;
	private String perSiteLLResultFileName;

	/*
	 * if parsimonyOnly is true, raxml is run with the '-y' flag,
	 * which causes it to produce the parsimony starting tree then
	 * exit. The parsimony tree can be retrieved with getParsimonyTree()
	 */
	private boolean parsimonyOnly = false;

	/*
	 * if parsimonyWithBL is true, run raxml twice. The first run generates a 
	 * parsimony tree, with topology only. The second run determines ML branch 
	 * lengths for the parsimony tree.
	 */
	private boolean parsimonyWithBL = false;

	public RAxMLRunner(int threads) {
		setThreadCount(threads);

		String bootstrapProp = System.getProperty(HandyConstants.SUPPORT_REPS);
		if(bootstrapProp != null) {
			try {
				int bs = Integer.parseInt(bootstrapProp);
				setBootstrapReps(bs);
			} catch(NumberFormatException nfe) {
				logger.error("Value for " + HandyConstants.SUPPORT_REPS + 
						" must be an integer, not: " + bootstrapProp);
			}
		}
	}

	public void run() {
		runName = null;

		if(algorithm == PARSIMONY_WITH_BL_ALGORITHM) {
			runRaxmlParsimonyWithBranchLengths();
		} else if(algorithm == PER_SITE_LL_ALGORITHM) {
			runRaxmlPerSiteLL();
		} else {
			String raxmlPath = ExecUtilities.getCommandPath(RAXML_HPC);
			if(raxmlPath == null || raxmlPath.trim().startsWith("no ")) {
				raxmlPath = ExecUtilities.getCommandPath(RAXML);
			}

			if(getThreadCount() > 1) {
				raxmlPath = ExecUtilities.getCommandPath("raxmlHPC-PTHREADS");
			}
			String raxmlOptions = "";

			//determine the file name for the alignment file
			File workingDir = new File(System.getProperty("user.dir"));
			try {
				File tempFile = File.createTempFile(getAlignment().getName(), ".phy", workingDir);
				
				if(deleteGeneratedFiles()) {
					tempFile.deleteOnExit();
				}
				String infileName = tempFile.getName();
				//write the alignment to the file, in phylip format
				writeAlignmentToFile(getAlignment(), infileName);
				//determine the raxml run name, and the various result file names
				//that will be generated
				String runName = getRaxmlRunName();

				int randomSeed = getOddRandomSeed();
				int bootstrapReps = getBootstrapReps();
				//configure raxml command string
				String algorithmOption = " -f a ";
				if(bootstrapReps == 0) {
					algorithmOption = " -f d ";
				}
				
				String bootstrapOption = " -x " + randomSeed + 
				" -N " + bootstrapReps;
				if(bootstrapReps == 0 || algorithm == PER_SITE_LL_ALGORITHM) {
					bootstrapOption = "";
				}
				raxmlOptions = algorithmOption + 
				" -m " + 
				getMatrix() + 
				" -s " 
				+ infileName + " -n " + runName + bootstrapOption;
				if(getThreadCount() > 1) {
					raxmlOptions = raxmlOptions + " -T " + getThreadCount();
				}

				if(algorithm == PARSIMONY_ALGORITHM) {
					if(getBootstrapReps() > 0) {
						raxmlOptions = raxmlOptions + " -Y -N " + getBootstrapReps();
					} else {
						raxmlOptions = raxmlOptions + " -y ";
					}
				}
				
				String raxmlCommand = "";
				raxmlCommand = raxmlPath + raxmlOptions;
				String executableNameOnly = new File(raxmlPath).getName();
				PEPRTracker.setTreeOptions(executableNameOnly + " " + raxmlOptions);
				//run raxml, and wait for it to finish
				CommandResults results = ExecUtilities.exec(raxmlCommand);
			} catch (IOException ioe) {
				logger.error(ioe);
			}
		}
	}
	
	private String getMatrix() {
		return matrix;
	}
	
	public void setMatrix(String m) {
		matrix = m;
	}

	private void runRaxmlPerSiteLL() {

		String raxmlPath = ExecUtilities.getCommandPath(RAXML_HPC);
		if(raxmlPath == null || raxmlPath.trim().startsWith("no ")) {
			raxmlPath = ExecUtilities.getCommandPath(RAXML);
		}

		if(getThreadCount() > 1) {
			raxmlPath = ExecUtilities.getCommandPath("raxmlHPC-PTHREADS");
		}
		String raxmlOptions = "";

		//determine the file name for the alignment file
		File workingDir = new File(System.getProperty("user.dir"));
		try {
			File tempFile = File.createTempFile("raxml_", ".phy", workingDir);
			if(deleteGeneratedFiles()) {
				tempFile.deleteOnExit();
			}
			String infileName = tempFile.getName();
			//write the alignment to the file, in phylip format
			writeAlignmentToFile(getAlignment(), infileName);
			
			//determine file name for tree file
			tempFile = File.createTempFile("raxml_", ".nwk", workingDir);
			String treeFileName = tempFile.getName();
			//write the trees to the file
			writeTreeFile(getPerSiteLLTrees(), treeFileName);
			
			//determine the raxml run name, and the various result file names
			//that will be generated
			String runName = getRaxmlRunName();

			//configure raxml command string
			String algorithmOption = " -f g ";

			raxmlOptions = algorithmOption + 
			" -m " + getMatrix() + " -s " 
			+ infileName + " -n " + runName + " -z " + treeFileName;
			if(getThreadCount() > 1) {
				raxmlOptions = raxmlOptions + " -T " + getThreadCount();
			}
			perSiteLLResultFileName = 
				"RAxML_perSiteLLs." + runName;
			String raxmlCommand = raxmlPath + raxmlOptions;
			System.out.println("ML branch length raxml command: " + raxmlCommand);
			CommandResults secondResults = ExecUtilities.exec(raxmlCommand);

		} catch (IOException ioe) {
			logger.error(ioe);
		}
	}
	
	private void runRaxmlParsimonyWithBranchLengths() {
		String raxmlPath = ExecUtilities.getCommandPath(RAXML_HPC);
		if(raxmlPath == null || raxmlPath.trim().startsWith("no ")) {
			raxmlPath = ExecUtilities.getCommandPath(RAXML);
		}

		if(getThreadCount() > 1) {
			raxmlPath = ExecUtilities.getCommandPath("raxmlHPC-PTHREADS");
		}
		String raxmlOptions = "";

		//determine the file name for the alignment file
		File workingDir = new File(System.getProperty("user.dir"));
		try {
			File tempFile = File.createTempFile("raxml_", ".phy", workingDir);
			if(deleteGeneratedFiles()) {
				tempFile.deleteOnExit();
			}
			String infileName = tempFile.getName();
			//write the alignment to the file, in phylip format
			writeAlignmentToFile(getAlignment(), infileName);
			//determine the raxml run name, and the various result file names
			//that will be generated
			String runName = getRaxmlRunName();

			//configure raxml command string
			String algorithmOption = " -f d ";

			raxmlOptions = algorithmOption + 
			" -m PROTGAMMAWAGF -s " 
			+ infileName + " -n " + runName + " -Y";
			if(getThreadCount() > 1) {
				raxmlOptions = raxmlOptions + " -T " + getThreadCount();
			}
			//run raxml, and wait for it to finish
			String raxmlCommand = raxmlPath + raxmlOptions;
			CommandResults firstResults = ExecUtilities.exec(raxmlCommand);

			//do a second raxml run to determine the ML branch lengths 
			//given the parsimony tree that was just calculated
			String parsimonyTreeFileName = "RAxML_parsimonyTree."
				+ getRaxmlRunName();

			//set algorithmOption to calculate branch lengths 
			//for a given start tree
			algorithmOption = " -f e";

			raxmlOptions = algorithmOption + 
			" -m " + matrix + " -s " 
			+ infileName + " -n " + runName + "BL" + " -t " + parsimonyTreeFileName;
			if(getThreadCount() > 1) {
				raxmlCommand = raxmlCommand + " -T " + getThreadCount();
			}
			
			raxmlCommand = raxmlPath + raxmlOptions;
			String executableNameOnly = new File(raxmlPath).getName();
			PEPRTracker.setTreeOptions(executableNameOnly + " " + raxmlOptions);
			CommandResults secondResults = ExecUtilities.exec(raxmlCommand);
			if(deleteGeneratedFiles()) {
				new File(parsimonyTreeFileName).deleteOnExit();
			}

		} catch (IOException ioe) {
			logger.error(ioe);
		}
	}
	
	public void setPerSiteLLTrees(String[] trees) {
		perSiteLLTrees = trees;
	}
	
	public String[] getPerSiteLLTrees() {
		return perSiteLLTrees;
	}
	
	public String[] getPerSiteLLResultFile() {
		String[] r = null;
		try {
			TextFile tf = new TextFile(perSiteLLResultFileName);
			r = tf.getAllLines();
		} catch (IOException ioe) {
			logger.error(ioe);
		}
		return r;
	}

	
	public String getBestTreeWithSupports() {
		String bestTreeWithSupportsFileName = "RAxML_bipartitions." 
			+ getRaxmlRunName();
		TextFile treeFile;
		StringBuffer sb = new StringBuffer();
		try {
			treeFile = new TextFile(bestTreeWithSupportsFileName);
			String[] treeFileLines = treeFile.getAllLines();
			for(int i = 0; i < treeFileLines.length; i++) {
				sb.append(treeFileLines[i]);
			}
		} catch (IOException ioe) {
			logger.error(ioe);
		}

		return sb.toString();
	}

	public String getBestTree() {
		String bestTreeWithSupportsFileName = "RAxML_result." 
			+ getRaxmlRunName();
		TextFile treeFile;
		StringBuffer sb = new StringBuffer();
		try {
			treeFile = new TextFile(bestTreeWithSupportsFileName);
			String[] treeFileLines = treeFile.getAllLines();
			for(int i = 0; i < treeFileLines.length; i++) {
				sb.append(treeFileLines[i]);
			}
		} catch (IOException ioe) {
			logger.error(ioe);
		}

		return sb.toString();
	}

	public String getParsimonyTree() {
		String r = null;
		String bestTreeWithSupportsFileName = "RAxML_parsimonyTree." 
			+ getRaxmlRunName();
		TextFile treeFile;
		StringBuffer sb = new StringBuffer();
		try {
			treeFile = new TextFile(bestTreeWithSupportsFileName);
			String[] treeFileLines = treeFile.getAllLines();
			for(int i = 0; i < treeFileLines.length; i++) {
				sb.append(treeFileLines[i]);
			}
			r = sb.toString();
			if(deleteGeneratedFiles()) {
				//since the results have been retrieved, the file can be deleted
				treeFile.getFile().deleteOnExit();
			}
		} catch (IOException e) {
		}

		return r;
	}

	public String getParsimonyWithBLTree() {
		String parsimonyTreeWithBL = "RAxML_result." 
			+ getRaxmlRunName() + "BL";
		TextFile treeFile;
		StringBuffer sb = new StringBuffer();
		try {
			treeFile = new TextFile(parsimonyTreeWithBL);
			String[] treeFileLines = treeFile.getAllLines();
			for(int i = 0; i < treeFileLines.length; i++) {
				sb.append(treeFileLines[i]);
			}
			if(deleteGeneratedFiles()) {
				//since the results have been retrieved, the file can be deleted
				treeFile.getFile().deleteOnExit();
			}
		} catch (IOException ioe) {
			logger.error(ioe);
		}

		return sb.toString();

	}

	private int getOddRandomSeed() {
		int r = (random.nextInt(1000000) * 2) + 1;
		return r;
	}

	private int getBootstrapReps() {
		return bootstrapReps;
	}

	public void setBootstrapReps(int reps) {
		bootstrapReps = reps;
	}

	public void setUseTaxonNames(boolean b) {
		useTaxonNames = b;
	}

	private void writeAlignmentToFile(SequenceAlignment sa,
			String infileName) throws IOException {
		FileWriter fw = new FileWriter(infileName);
		String phylipAlignment;
		if(useTaxonNames) {
			phylipAlignment = sa.getAlignmentAsExtendedPhylipUsingTaxonNames();
		} else {
			phylipAlignment = sa.getAlignmentAsExtendedPhylipUsingSequenceNames();
		}
		fw.write(phylipAlignment);
		fw.flush();
		fw.close();
	}
	
	private void writeTreeFile(String[] trees, String fileName) throws IOException {
		FileWriter fw = new FileWriter(fileName);
		for(int i = 0; i < trees.length; i++) {
			fw.write(trees[i]);
			fw.write("\n");
		}
		fw.flush();
		fw.close();
	}

	public SequenceAlignment getAlignment() {
		return alignment;
	}

	public void setAlignment(SequenceAlignment alignment) {
		this.alignment = alignment;
	}

	private String getRaxmlRunName() {
		if(runName == null) {
			runName = getAvailableRaxmlRunName();
		}
		if(deleteGeneratedFiles()) {
			//target info file for deleting on exit
			new File("RAxML_info." + runName).deleteOnExit();
		}
		return runName;
	}

	/**
	 * Uses a set of trees to decorate a single tree with support values.
	 * This is currently only done on localhost, which shouldn't be a problem 
	 * because it's not very intensive.
	 * 
	 * @param mainTree
	 * @param supprtTrees
	 * @return
	 */
	public String getSupportDecoratedTree(SequenceAlignment alignment, 
			String mainTree, String[] supportTrees) {
		String r = null;
		
		String fullTreeFileName;
		String supportTreeFileName;
		File workingDir = new File(System.getProperty("user.dir"));
		try {
			//write full tree file
			File fullTreeFile = File.createTempFile("raxml_", ".nwk", workingDir);
			//fullTreeFile.deleteOnExit();
			fullTreeFileName = fullTreeFile.getName();
			logger.info("write full tree to " + fullTreeFileName);
			logger.info("mainTree: " + mainTree);
			FileWriter fullTreeWriter = new FileWriter(fullTreeFile);
			fullTreeWriter.write(mainTree + "\n");
			fullTreeWriter.flush();
			fullTreeWriter.close();
		
		    //write the support trees to a file
			File supportTreeFile = File.createTempFile("raxml_", ".sup", workingDir);
			//supportTreeFile.deleteOnExit();
			supportTreeFileName = supportTreeFile.getName();
			logger.info("write support trees to " + supportTreeFileName);
			FileWriter supportTreeWriter = new FileWriter(supportTreeFile);
			for(int i= 0; i < supportTrees.length; i++) {
				supportTreeWriter.write(supportTrees[i] + "\n");
			}
			supportTreeWriter.flush();
			supportTreeWriter.close();
			
			//write the sequence alignment file
			File alignmentFile = File.createTempFile("raxml_", ".phy", workingDir);
			String alignmentFileName = alignmentFile.getName();
			writeAlignmentToFile(alignment, alignmentFileName);

			//get the raxml command path
			String raxmlPath = ExecUtilities.getCommandPath("raxmlHPC-PTHREADS");
			if(getThreadCount() <= 1) {
				raxmlPath = ExecUtilities.getCommandPath("raxmlHPC");
			}
			//create raxml command 
			String runName = getAvailableRaxmlRunName();
			String raxmlDecorateCommand = raxmlPath + " -f b" +
			" -z " + supportTreeFileName + " -t " + fullTreeFileName + 
			" -T " + getThreadCount() + " -m " + matrix + " -n " + runName +
			" -s " + alignmentFileName;
			
			//run RAxML command
			CommandResults results = 
					ExecUtilities.exec(raxmlDecorateCommand);
			//read result from file
			String resultFileName = "RAxML_bipartitions." + runName;
			TextFile resultFile = new TextFile(resultFileName);
			r = resultFile.toString();
			if(deleteGeneratedFiles()) {
				resultFile.getFile().deleteOnExit();
			}
		} catch (IOException ioe) {
			logger.error(ioe);
		}

		return r;
	}
	
	private static synchronized String getAvailableRaxmlRunName() {
		String r = null;

		int i = lastFileNumberUsed++;
		r = "raxmlRun"+ i;

		while(usedRunNames.contains(r) || new File("RAxML_info." + r).exists()) {
			i++;
			r = "raxmlRun"+ i;			
		}

		lastFileNumberUsed =i;
		usedRunNames.add(r);
		return r;
	}

	private int getThreadCount() {
		return threadCount;
	}

	private void setThreadCount(int threadCount) {
		this.threadCount = threadCount;
	}

	public void setParsimonyOnly(boolean parsOnly) {
		parsimonyOnly = parsOnly;
		if(parsOnly) {
			algorithm = PARSIMONY_ALGORITHM;
			setParsimonyWithBL(false);
		}
	}

	public void setParsimonyWithBL(boolean parsBL) {
		parsimonyWithBL = parsBL;
		if(parsBL) {
			algorithm = PARSIMONY_WITH_BL_ALGORITHM;
			setParsimonyOnly(false);
		}
	}

	public void setPerSiteLogLikelihoods(boolean psll) {
		if(psll) {
			algorithm = PER_SITE_LL_ALGORITHM;
		}
	}
	
	public boolean deleteGeneratedFiles() {
		return deleteGeneratedFiles;
	}

	public void setDeleteGeneratedFiles(boolean deleteGeneratedFiles) {
		this.deleteGeneratedFiles = deleteGeneratedFiles;
	}
}
