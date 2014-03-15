
package edu.vt.vbi.ci.pepr.tree;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashSet;
import java.util.Random;

import edu.vt.vbi.ci.pepr.alignment.SequenceAlignment;
import edu.vt.vbi.ci.pepr.alignment.SequenceAlignmentParser;
import edu.vt.vbi.ci.util.CommandResults;
import edu.vt.vbi.ci.util.ExecUtilities;
import edu.vt.vbi.ci.util.HandyConstants;
import edu.vt.vbi.ci.util.RemoteHost;
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
	private SequenceAlignment alignment;
	private int bootstrapReps = 10;
	private String runName;
	private int threadCount;
	private RemoteHost host;
	private boolean useTaxonNames = true;
	private boolean doPerSiteLikelihoods = false;
	private boolean deleteGeneratedFiles = true;
	private static int lastFileNumberUsed = 0;
	private int algorithm = ML_ALGORITHM;
	private String matrix = "PROTGAMMALGF";
	private String[] perSiteLLTrees;
	private String perSiteLLResultFileName;

	/**
	 * These are the Remote Hosts for running bootstrap replicates.
	 * If a Host can run multiple simultaneous jobs it can be listed
	 * multiple times. For now, there is no sophisticated load balancing.
	 * The number of bootstrap replicates is divided by the number of
	 * bootstrapHosts entries and that number of replicated is sent to each
	 * host entry (which might mean sending multiple jobs to the same actual
	 * host).
	 */
	private RemoteHost[] bootstrapHosts;

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

	public static void main(String[] args) {
		String alignmentFileName = "/Users/enordber/vbi/rickettsia/set_2029_align";
		try {
			SequenceAlignment sa = 
				SequenceAlignmentParser.parseFastaAlignmentFile(alignmentFileName);
			RAxMLRunner raxmlRunner = new RAxMLRunner(1);
			raxmlRunner.setAlignment(sa);
			raxmlRunner.setBootstrapReps(10);
			raxmlRunner.run();
			System.out.println(raxmlRunner.getBestTreeWithSupports());
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public RAxMLRunner(int threads) {
		setThreadCount(threads);
		host = RemoteHost.getLocalHost();

		String bootstrapProp = System.getProperty(HandyConstants.SUPPORT_REPS);
		if(bootstrapProp != null) {
			try {
				int bs = Integer.parseInt(bootstrapProp);
				setBootstrapReps(bs);
			} catch(NumberFormatException nfe) {
				System.out.println("Value for " + HandyConstants.SUPPORT_REPS + 
						" must be an integer, not: " + bootstrapProp);
			}
		}
	}

	public void run() {
//		System.out.println("RAxMLRunner.run()");
//		System.out.println("algorithm: " + algorithm);
		runName = null;

		if(algorithm == PARSIMONY_WITH_BL_ALGORITHM) {
			runRaxmlParsimonyWithBranchLengths();
		} else if(algorithm == PER_SITE_LL_ALGORITHM) {
			runRaxmlPerSiteLL();
		} else {
			String raxmlPath = host.getCommandPath(RAXML_HPC);
			if(raxmlPath == null || raxmlPath.trim().startsWith("no ")) {
				raxmlPath = host.getCommandPath(RAXML);
			}

			if(getThreadCount() > 1) {
				raxmlPath = host.getCommandPath("raxmlHPC-PTHREADS");
			}
			//		System.out.println("raxml path: " + raxmlPath);
			String raxmlCommand = "";

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
				host.copyToRemoteHost(infileName, infileName);
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
				raxmlCommand = raxmlPath + algorithmOption + 
				" -m " + 
				getMatrix() + 
				" -s " 
				+ infileName + " -n " + runName + bootstrapOption;
				if(getThreadCount() > 1) {
					raxmlCommand = raxmlCommand + " -T " + getThreadCount();
				}

				if(algorithm == PARSIMONY_ALGORITHM) {
					if(getBootstrapReps() > 0) {
						raxmlCommand = raxmlCommand + " -Y -N " + getBootstrapReps();
					} else {
						raxmlCommand = raxmlCommand + " -y ";
					}
				}
				System.out.println("raxml command: " + raxmlCommand);
				//run raxml, and wait for it to finish
				CommandResults results = host.executeCommand(raxmlCommand);
				//grab the tree string from the result file

			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
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

		String raxmlPath = host.getCommandPath(RAXML_HPC);
		if(raxmlPath == null || raxmlPath.trim().startsWith("no ")) {
			raxmlPath = host.getCommandPath(RAXML);
		}

		if(getThreadCount() > 1) {
			raxmlPath = host.getCommandPath("raxmlHPC-PTHREADS");
		}
		//		System.out.println("raxml path: " + raxmlPath);
		String raxmlCommand = "";

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
			host.copyToRemoteHost(infileName, infileName);
			String runName = getRaxmlRunName();

			//configure raxml command string
			String algorithmOption = " -f g ";

			raxmlCommand = raxmlPath + algorithmOption + 
			" -m " + getMatrix() + " -s " 
			+ infileName + " -n " + runName + " -z " + treeFileName;
			if(getThreadCount() > 1) {
				raxmlCommand = raxmlCommand + " -T " + getThreadCount();
			}
			perSiteLLResultFileName = 
				"RAxML_perSiteLLs." + runName;
			System.out.println("ML branch length raxml command: " + raxmlCommand);
			CommandResults secondResults = host.executeCommand(raxmlCommand);

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	private void runRaxmlParsimonyWithBranchLengths() {
		System.out.println("RAxMLRunner.runRaxmlParsimonyWithBranchLengths()");

		String raxmlPath = host.getCommandPath(RAXML_HPC);
		if(raxmlPath == null || raxmlPath.trim().startsWith("no ")) {
			raxmlPath = host.getCommandPath(RAXML);
		}

		if(getThreadCount() > 1) {
			raxmlPath = host.getCommandPath("raxmlHPC-PTHREADS");
		}
		//		System.out.println("raxml path: " + raxmlPath);
		String raxmlCommand = "";

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
			host.copyToRemoteHost(infileName, infileName);
			String runName = getRaxmlRunName();

			//configure raxml command string
			String algorithmOption = " -f d ";

			raxmlCommand = raxmlPath + algorithmOption + 
			" -m PROTGAMMAWAGF -s " 
			+ infileName + " -n " + runName + " -Y";
			if(getThreadCount() > 1) {
				raxmlCommand = raxmlCommand + " -T " + getThreadCount();
			}
			System.out.println("parsimony tree raxml command: " + raxmlCommand);
			//run raxml, and wait for it to finish
			CommandResults firstResults = host.executeCommand(raxmlCommand);

			//do a second raxml run to determine the ML branch lengths 
			//given the parsimony tree that was just calculated
			String parsimonyTreeFileName = "RAxML_parsimonyTree."
				+ getRaxmlRunName();

			//set algorithmOption to calculate branch lengths 
			//for a given start tree
			algorithmOption = " -f e";

			raxmlCommand = raxmlPath + algorithmOption + 
			" -m " + matrix + " -s " 
			+ infileName + " -n " + runName + "BL" + " -t " + parsimonyTreeFileName;
			if(getThreadCount() > 1) {
				raxmlCommand = raxmlCommand + " -T " + getThreadCount();
			}
			System.out.println("ML branch length raxml command: " + raxmlCommand);
			CommandResults secondResults = host.executeCommand(raxmlCommand);
			if(deleteGeneratedFiles()) {
				new File(parsimonyTreeFileName).deleteOnExit();
			}

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
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
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return r;
	}

	
	public String getBestTreeWithSupports() {
		String bestTreeWithSupportsFileName = "RAxML_bipartitions." 
			+ getRaxmlRunName();
		host.copyFromRemoteHost(bestTreeWithSupportsFileName,
				bestTreeWithSupportsFileName);
		TextFile treeFile;
		StringBuffer sb = new StringBuffer();
		try {
			treeFile = new TextFile(bestTreeWithSupportsFileName);
			String[] treeFileLines = treeFile.getAllLines();
			for(int i = 0; i < treeFileLines.length; i++) {
				sb.append(treeFileLines[i]);
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		return sb.toString();
	}

	public String getBestTree() {
		String bestTreeWithSupportsFileName = "RAxML_result." 
			+ getRaxmlRunName();
		host.copyFromRemoteHost(bestTreeWithSupportsFileName,
				bestTreeWithSupportsFileName);
		TextFile treeFile;
		StringBuffer sb = new StringBuffer();
		try {
			treeFile = new TextFile(bestTreeWithSupportsFileName);
			String[] treeFileLines = treeFile.getAllLines();
			for(int i = 0; i < treeFileLines.length; i++) {
				sb.append(treeFileLines[i]);
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		return sb.toString();
	}

	public String getParsimonyTree() {
		String r = null;
		String bestTreeWithSupportsFileName = "RAxML_parsimonyTree." 
			+ getRaxmlRunName();
		host.copyFromRemoteHost(bestTreeWithSupportsFileName,
				bestTreeWithSupportsFileName);
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
			// TODO Auto-generated catch block
//			e.printStackTrace();
		}

		return r;
	}

	public String getParsimonyWithBLTree() {
		String parsimonyTreeWithBL = "RAxML_result." 
			+ getRaxmlRunName() + "BL";
		host.copyFromRemoteHost(parsimonyTreeWithBL,
				parsimonyTreeWithBL);
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
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
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
		//		if(reps > 1) {
		bootstrapReps = reps;
		//		}
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
			runName = getAvailableRaxmlRunName(getHost());
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
		System.out.println("RAxMLRunner.getSupportDecoratedTree()");
		String r = null;
		RemoteHost localHost = RemoteHost.getLocalHost();
		
		String fullTreeFileName;
		String supportTreeFileName;
		File workingDir = new File(System.getProperty("user.dir"));
		try {
			//write full tree file
			File fullTreeFile = File.createTempFile("raxml_", ".nwk", workingDir);
			//fullTreeFile.deleteOnExit();
			fullTreeFileName = fullTreeFile.getName();
			System.out.println("write full tree to " + fullTreeFileName);
			System.out.println("mainTree: " + mainTree);
			FileWriter fullTreeWriter = new FileWriter(fullTreeFile);
			fullTreeWriter.write(mainTree + "\n");
			fullTreeWriter.flush();
			fullTreeWriter.close();
		
		//write the support trees to a file
			File supportTreeFile = File.createTempFile("raxml_", ".sup", workingDir);
			//supportTreeFile.deleteOnExit();
			supportTreeFileName = supportTreeFile.getName();
			System.out.println("write support trees to " + supportTreeFileName);
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
			String raxmlPath = localHost.getCommandPath("raxmlHPC-PTHREADS");
			if(getThreadCount() <= 1) {
				raxmlPath = localHost.getCommandPath("raxmlHPC");
			}
			//create raxml command 
			String runName = getAvailableRaxmlRunName(localHost);
			String raxmlDecorateCommand = raxmlPath + " -f b" +
			" -z " + supportTreeFileName + " -t " + fullTreeFileName + 
			" -T " + getThreadCount() + " -m " + matrix + " -n " + runName +
			" -s " + alignmentFileName;
			
			//run RAxML command
			System.out.println("raxml decorate command: " + raxmlDecorateCommand);
			CommandResults results = 
				localHost.executeCommand(raxmlDecorateCommand);
			//read result from file
			String resultFileName = "RAxML_bipartitions." + runName;
			TextFile resultFile = new TextFile(resultFileName);
			r = resultFile.toString();
			if(deleteGeneratedFiles()) {
				resultFile.getFile().deleteOnExit();
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		return r;
	}
	
	private static synchronized String getAvailableRaxmlRunName(RemoteHost host) {
		String r = null;

		int i = lastFileNumberUsed++;
		r = "raxmlRun"+ i;

		while(usedRunNames.contains(r) || host.fileExists("RAxML_info." + r)) {
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

	RemoteHost getHost() {
		return host;
	}

	public void setHost(RemoteHost host) {
		this.host = host;
	}

	public void setParsimonyOnly(boolean parsOnly) {
		parsimonyOnly = parsOnly;
		if(parsOnly) {
			algorithm = PARSIMONY_ALGORITHM;
			setParsimonyWithBL(false);
		}
	}

	public void setParsimonyWithBL(boolean parsBL) {
//		System.out.println("RAxMLRunner.setParsimonyWithBL() " + parsBL);
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
	
	public void setBootstrapHosts(RemoteHost[] hosts) {
		bootstrapHosts = hosts;
	}

	public boolean deleteGeneratedFiles() {
		return deleteGeneratedFiles;
	}

	public void setDeleteGeneratedFiles(boolean deleteGeneratedFiles) {
		this.deleteGeneratedFiles = deleteGeneratedFiles;
	}
}
