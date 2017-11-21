package edu.vt.vbi.ci.util;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;

import org.apache.log4j.Logger;

public class ExecUtilities {

	private static Logger logger = Logger.getLogger(ExecUtilities.class);
	/**
	 * Executes the given command via Runtime.exec()
	 * and returns the results of stdout and stderr
	 * as a CommandResults object.
	 * 
	 * @param command
	 * @return
	 */
	public static CommandResults exec(String command) {
		CommandResults r = null;
		try {
			logger.info(command);
			Process proc = Runtime.getRuntime().exec(command);
			r = getResultFromProcess(proc);
			logger.info(command + " has rc " + r.getRc() + "\n");
			if(r.getRc() != 0) {
				for (String s : r.getStderr())
				{
					logger.error("  " + s);
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
		return r;
	}

	/**
	 * Executes the given command via Runtime.exec()
	 * and returns the results of stdout and stderr
	 * as a CommandResults object.
	 * 
	 * @param command
	 * @return
	 */
	public static CommandResults exec(String[] command) {
		CommandResults r = null;
		try {
			Process proc = Runtime.getRuntime().exec(command);
			r = getResultFromProcess(proc);
		} catch (IOException e) {
			e.printStackTrace();
		} catch (InterruptedException e) {
			e.printStackTrace();
		}

		if (r.getRc() != 0)
		{
			logger.warn("Error return " + r.getRc() + " from command " + command);
		}


		return r;
	}

	/**
	 * Captures the output from the Process and wraps it in a 
	 * CommandResults Object.
	 * 
	 * @param proc
	 * @return
	 * @throws InterruptedException
	 */
	private static CommandResults getResultFromProcess(Process proc) throws InterruptedException {
		CommandResults r = null;
		InputStream is = proc.getInputStream();
		InputStream errStream = proc.getErrorStream();
		BufferedReader reader = 
				new BufferedReader(new InputStreamReader(is));
		BufferedReader errorReader = 
				new BufferedReader(new InputStreamReader(errStream));

		ReaderRunnable errorRead = new ReaderRunnable(errorReader);
		Thread errorThread = new Thread(errorRead);
		errorThread.start();

		ReaderRunnable outReader = new ReaderRunnable(reader);
		Thread outThread = new Thread(outReader);
		outThread.start();

		errorThread.join();
		outThread.join();

		//call destroy() on process to free resources. The streams
		//are not closed automatically, and may eventually 
		//accumulate and result in an "IOException: Too many open files"
		int rc = proc.waitFor();
		r = new CommandResults(outReader.getLinesRead(),
				errorRead.getLinesRead(), rc);

		try {
			proc.getErrorStream().close();
			proc.getInputStream().close();
			proc.getOutputStream().close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		proc.destroy();

		return r;

	}
	/**
	 * This is a Runnable to monitor and capture the
	 *  output from a command.
	 * 
	 * @author enordber
	 *
	 */
	private static class ReaderRunnable implements Runnable {
		private boolean debug = false;
		private BufferedReader reader;
		private ArrayList<String> linesRead;

		public ReaderRunnable(BufferedReader reader) {
			this.reader = reader;
			linesRead = new ArrayList<String>();
		}
		public void run() {
			String line;
			try {
				while((line = reader.readLine()) != null) {
					linesRead.add(line);
					if(debug) {
						System.out.println(line);
					}
				}
			} catch (IOException e) {
				e.printStackTrace();
			}
		}

		/**
		 * Returns an array of all lines read from the reader.
		 * @return
		 */
		public String[] getLinesRead(){
			String[] r = new String[linesRead.size()];
			linesRead.toArray(r);
			return r;
		}


	}

	/**
	 * Returns the absolute path to the given executable. This method
	 * relies on the 'which' utility. Therefore, it only works on
	 * Unix-type systems that have 'which' at /usr/bin/which. 
	 */
	public static String getCommandPath(String cmd) {
		String r = null;

		//check System properties for a path for this command
		r = System.getProperty(cmd);

		if(r == null ){
			//if no path is in System properties, use 'which' 
			//utility to look for command path
			String whichCmd = "/usr/bin/which " + cmd;
			CommandResults whichRes= ExecUtilities.exec(whichCmd);

			if(whichRes == null) {
				logger.info("unable to get command path for '" + cmd + "'");
			} else {
				String[] stdout = whichRes.getStdout();
				if(stdout.length > 0) {
					r = stdout[0];
				}
			}
		}
		return r;
	}

}
