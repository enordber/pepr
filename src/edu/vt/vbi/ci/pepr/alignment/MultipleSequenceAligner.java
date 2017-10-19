package edu.vt.vbi.ci.pepr.alignment;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import edu.vt.vbi.ci.util.CommandResults;
import edu.vt.vbi.ci.util.ExecUtilities;
import edu.vt.vbi.ci.util.HandyConstants;
import edu.vt.vbi.ci.util.file.FastaSequenceFile;
import edu.vt.vbi.ci.util.file.FastaSequenceSet;
import edu.vt.vbi.ci.util.file.FastaSequenceSetImpl;

/**
 * This class will serve as a wrapper for muscle,
 * so other java classes can easily do MSA with muscle
 * without all the details of dealing with muscle.
 * For now, I will hard code the info for muscle,
 * but that can be changed later
 * 
 * inputs:
 * 		single .faa or .ffn file
 * 		multiple .faa or .ffn files (must all be of same type)

 * 
 * outputs:
 * 		fasta-format alignment
 * 		nexus format alignment
 * 		output file for each input file
 * 		single output file with all alignments (concatenated)
 * 
 * @author enordber
 *
 */
public class MultipleSequenceAligner {

	private String inputFormat;
	private String outputFormat;
	private boolean concatenatedOutput;
	private boolean individualOutput;

	/*
	 * HashMap to store alignments when caching is enabled.
	 */
	private HashMap alignmentCache;

	private ArrayList fastaInputFiles;

	/**
	 * Store and command-line options to be used for muscle
	 */
	private ArrayList muscleOptions;


	/*
	 * The RemoteHost to be used for running muscle
	 */
//	private RemoteHost muscleHost;

	/*
	 * The directory to use on the muscle host
	 */
	private String muscleHostWorkingDir = "";

	/*
	 * Path to muscle executable on muscle host
	 */
	private String muscleHostMusclePath = "/usr/bin/muscle";

	/*
	 * working directory. The location for files related to tasks of this class.
	 */
	private File workingDir;

	/*
	 * The Alignment.
	 */
	private SequenceAlignment msa;

	public MultipleSequenceAligner() {
//		muscleHost = getLocalHost();
		determineMusclePath();
		addMuscleOption("-fasta");
		addMuscleOption("-stable");
		addMuscleOption("-quiet");
	}

	public SequenceAlignment getMSA(FastaSequenceSet in) {
		SequenceAlignment r = null;
		try {
			//create temporary sequence file for input
			File workingDir = getWorkingDir();
			File tempLocalIn = File.createTempFile("tmp_", ".seq", workingDir);
			FileWriter fw = new FileWriter(tempLocalIn);
			int sequenceCount = in.getSequenceCount();
			for(int i = 0; i < sequenceCount; i++) {
				String[] seq = in.getSequence(i);
				for(int j = 0; j < seq.length; j++) {
					fw.write(seq[j]);
					fw.write("\n");
				}
			}
			fw.close();

			String muscleInFileName = tempLocalIn.getPath();			

			String log = " -verbose -log " + muscleInFileName + ".log ";
			log = ""; //comment out to turn on muscle logging
			
			//run muscle, and capture output
			String muscleCommand = muscleHostMusclePath + "  " +
			getMuscleOptionString() + log +
			"-in " + muscleInFileName;	

			CommandResults muscleOut = 
				ExecUtilities.exec(muscleCommand);

			StringBuffer alignmentBuff = new StringBuffer();
			String[] stdout = muscleOut.getStdout();
			for(int i = 0; i < stdout.length; i++) {
				alignmentBuff.append(stdout[i]);
				alignmentBuff.append("\n");
			}
			String alignmentString = alignmentBuff.toString();

			//delete input file on local host
			tempLocalIn.delete();

			//create SequenceAlignment Object from alignment result
			FastaSequenceSet alignment =
				new FastaSequenceSetImpl(alignmentString);
			r = SequenceAlignmentParser.parseFastaAlignment(alignment);
			r.setName(in.getName() + HandyConstants.ALIGNMENT_FILE_SUFFIX);
		} catch(IOException ioe) {
			ioe.printStackTrace();
		}

		return r;
	}

	public SequenceAlignment getProfileMSA(FastaSequenceSet in1, FastaSequenceSet in2) {
		SequenceAlignment r = null;
		try {
			//create temporary sequence file for input
			File workingDir = getWorkingDir();
			File tempLocalIn1 = File.createTempFile("tmp_", ".seq", workingDir);
			FileWriter fw1 = new FileWriter(tempLocalIn1);
			int sequenceCount = in1.getSequenceCount();
			for(int i = 0; i < sequenceCount; i++) {
				String[] seq = in1.getSequence(i);
				for(int j = 0; j < seq.length; j++) {
					fw1.write(seq[j]);
					fw1.write("\n");
				}
			}
			fw1.close();

			File tempLocalIn2 = File.createTempFile("tmp_", ".seq", workingDir);
			FileWriter fw2 = new FileWriter(tempLocalIn2);
			sequenceCount = in2.getSequenceCount();
			for(int i = 0; i < sequenceCount; i++) {
				String[] seq = in2.getSequence(i);
				for(int j = 0; j < seq.length; j++) {
					fw2.write(seq[j]);
					fw2.write("\n");
				}
			}
			fw2.close();

			String muscleIn1FileName = tempLocalIn1.getPath();
			String muscleIn2FileName = tempLocalIn2.getPath();

			String log = " -verbose -log " + muscleIn1FileName + ".log ";
			log = ""; //comment out to turn on muscle logging
			
			//run muscle, and capture output
			String muscleCommand = muscleHostMusclePath + " -profile  " +
			getMuscleOptionString() + log +
			" -in1 " + muscleIn1FileName + " -in2 " + muscleIn2FileName;	

			System.out.println(muscleCommand);
			
			CommandResults muscleOut = 
					ExecUtilities.exec(muscleCommand);

			StringBuffer alignmentBuff = new StringBuffer();
			String[] stdout = muscleOut.getStdout();
			for(int i = 0; i < stdout.length; i++) {
				alignmentBuff.append(stdout[i]);
				alignmentBuff.append("\n");
			}
			String alignmentString = alignmentBuff.toString();

			//create SequenceAlignment Object from alignment result
			FastaSequenceSet alignment =
				new FastaSequenceSetImpl(alignmentString);
			r = SequenceAlignmentParser.parseFastaAlignment(alignment);
			r.setName(in1.getName() + HandyConstants.ALIGNMENT_FILE_SUFFIX);
		} catch(IOException ioe) {
			ioe.printStackTrace();
		}

		return r;
	}

	
	private File getWorkingDir() {
		if(workingDir == null) {
			workingDir =  new File(System.getProperty("user.dir"));
		}
		return workingDir;
	}

	public void setWorkingDir(File dir) {
		workingDir = dir;
	}

	private String getMuscleOptionString() {
		StringBuffer sb = new StringBuffer();
		if(muscleOptions != null) {
			String[] options = new String[muscleOptions.size()];
			muscleOptions.toArray(options);
			for(int i = 0; i < options.length; i++) {
				sb.append(options[i]);
				sb.append(" ");
			}
		}
		return sb.toString();
	}

	public static SequenceAlignment getMSA(String fastaSequence) {
		SequenceAlignment r = null;
		return r;
	}

	public String getInputFormat() {
		return inputFormat;
	}

	public void setInputFormat(String inputFormat) {
		this.inputFormat = inputFormat;
	}

	public String getOutputFormat() {
		return outputFormat;
	}

	public void setOutputFormat(String outputFormat) {
		this.outputFormat = outputFormat;
	}

	public boolean isConcatenatedOutput() {
		return concatenatedOutput;
	}

	public void setConcatenatedOutput(boolean concatenatedOutput) {
		this.concatenatedOutput = concatenatedOutput;
	}

	public boolean isIndividualOutput() {
		return individualOutput;
	}

	public void setIndividualOutput(boolean individualOutput) {
		this.individualOutput = individualOutput;
	}

	public void addInputFile(FastaSequenceFile infile) {
		if(fastaInputFiles == null) {
			fastaInputFiles = new ArrayList();
		}
		fastaInputFiles.ensureCapacity(fastaInputFiles.size() + 1);
		fastaInputFiles.add(infile);
	}

	public SequenceAlignment getMSA() {
		if(msa == null) {
			//do msa
			FastaSequenceSet infile = (FastaSequenceSet)fastaInputFiles.get(0);
			msa = getMSA(infile);
		}
		return msa;
	}

	public String getMuscleHostWorkingDir() {
		return muscleHostWorkingDir;
	}

	public void setMuscleHostWorkingDir(String muscleHostWorkingDir) {
		this.muscleHostWorkingDir = muscleHostWorkingDir;
	}

	private void determineMusclePath() {
		muscleHostMusclePath = ExecUtilities.getCommandPath("muscle-3.6");
	}

	public void addMuscleOption(String opt) {
		if(muscleOptions == null){
			muscleOptions = new ArrayList(1);
		}
		muscleOptions.ensureCapacity(muscleOptions.size()+1);
		muscleOptions.add(opt);
	}
}
