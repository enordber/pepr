package edu.vt.vbi.ci.pepr.tree;

import java.io.IOException;

import edu.vt.vbi.ci.pepr.alignment.SequenceAlignment;
import edu.vt.vbi.ci.pepr.alignment.SequenceAlignmentParser;
import edu.vt.vbi.ci.util.CommandLineProperties;
import edu.vt.vbi.ci.util.HandyConstants;
import edu.vt.vbi.ci.util.file.TextFile;

public class TreeBuilderComparator {

	public static void main(String[] args) {
		TreeBuilderComparator tbc = new TreeBuilderComparator(args);
	}
	
	public TreeBuilderComparator(String[] args) {
		CommandLineProperties clp = new CommandLineProperties(args);
		
		//get names of alignment files
		String[] alignmentFileNames = clp.getValues(HandyConstants.FILE);
		//load alignments
		//assume input files are phylip-format
		SequenceAlignment[] alignments = null;
		if(clp.getValues(HandyConstants.FASTA, HandyConstants.FALSE)[0].equalsIgnoreCase(HandyConstants.TRUE)){
			alignments = loadFastaAlignments(alignmentFileNames);
		} else {
		    alignments = loadPhylipAlignments(alignmentFileNames);
		}
		
		String matrix = clp.getValues(HandyConstants.ML_MATRIX, "PROTGAMMAWAGF")[0];
		
		//get number of threads to use for raxml
		int threads = Integer.parseInt(clp. 
				getValues(HandyConstants.TREE_THREADS, "1")[0]);
		
		//do each comparison
		for(int i = 0; i < alignments.length; i++) {
			compareFastTreeAndRaxmlForAlignment(alignments[i], matrix, threads);
		}
	}
	
	public TreeBuilderComparator() {
		// TODO Auto-generated constructor stub
	}

	private SequenceAlignment[] loadPhylipAlignments(String[] fileNames) {
		SequenceAlignment[] r = new SequenceAlignment[fileNames.length];

		for(int i = 0; i < r.length; i++) {
			try{
				TextFile alignmentFile = new TextFile(fileNames[i]);
				String alignmentString = alignmentFile.toString();
				r[i] = 
					SequenceAlignmentParser.parsePhylipAlignment(alignmentString);
			} catch(IOException ioe) {
				ioe.printStackTrace();
			}
		}
		return r;
	}
	
	private SequenceAlignment[] loadFastaAlignments(String[] fileNames) {
		SequenceAlignment[] r = new SequenceAlignment[fileNames.length];

		for(int i = 0; i < r.length; i++) {
			try{
				r[i] = 
					SequenceAlignmentParser.parseFastaAlignmentFile(fileNames[i]);
			} catch(IOException ioe) {
				ioe.printStackTrace();
			}
		}
		return r;
	}
	
	public void compareFastTreeAndRaxmlForAlignment(SequenceAlignment alignment, 
			String matrix, int threads) {
		FastTreeRunner runner = new FastTreeRunner();
		runner.setAlignment(alignment);
		System.out.println("run fastTreeRunner");
		long beforeFastTree = System.currentTimeMillis();
		runner.run();
		String fastTreeString = runner.getResult();
		long afterFastTree = System.currentTimeMillis();
		System.out.println("fasttree : " + fastTreeString);
		long fastSeconds = (afterFastTree-beforeFastTree)/1000;
		System.out.println("FastTree s: " + fastSeconds);

		//run raxml on same alignment
		RAxMLRunner raxmlRunner = new RAxMLRunner(threads);
		if(!alignment.hasTaxonNames()) {
			raxmlRunner.setUseTaxonNames(false);
		}
		raxmlRunner.setMatrix(matrix);
		raxmlRunner.setAlignment(alignment);
		raxmlRunner.setBootstrapReps(0);
		long beforeRaxml = System.currentTimeMillis();
		raxmlRunner.run();
		String raxmlTreeString = raxmlRunner.getBestTree();
		long afterRaxml = System.currentTimeMillis();
		System.out.println("raxml tree: " + raxmlTreeString);
		long raxmlSeconds = (afterRaxml-beforeRaxml)/1000;
		System.out.println("raxml s: " + raxmlSeconds);

		TreeI fastTree = TreeParser.parseTreeString(fastTreeString);
		TreeI raxmlTree = TreeParser.parseTreeString(raxmlTreeString);

		int rfDistance = fastTree.getRobinsonFouldsDistance(raxmlTree);
		System.out.println("RF distance between trees: " + rfDistance);

		String compareOut = "COMPARE OUT " +
		"\tntax: " + alignment.getNTax() + "\tlength: " + alignment.getLength()
		+ "\t FastTree sec: " + fastSeconds + "\tRAxML sec: " + raxmlSeconds 
		+ "\tRF: " + rfDistance;
		
		System.out.println(compareOut);

	}

}
