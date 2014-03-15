package edu.vt.vbi.ci.util;

import java.net.InetAddress;
import java.net.UnknownHostException;
import java.util.HashMap;

/**
 * This class enables the execution of commands on 
 * a remote host. Commands are executed via ssh 
 * (or scp in the case of file transfers).
 * This requires that the remote host has a public key
 * corresponding to a private key available on 
 * the local host. This avoids having ssh prompt
 * for a password. 
 * 
 * This class requires ssh on the local machine. As far as 
 * I know, this means it must be some flavor of Unix.
 * 
 * This class primarily formats the command, which is
 * the executed using ExecUtilites.exec().
 * 
 * Despite the name, a RemoteHost Object may represent 
 * the local host. The class attempts to determine when this
 * is the case and bypass the use of ssh and scp, to avoid
 * unnecessary overhead.
 * 
 * This is being developed as part of the phylogenetic
 * tree building pipeline with Kelly Williams and
 * Allan Dickerman.
 * 
 * @author enordber
 *
 */
public class RemoteHost {

	private boolean debug = false;

	private String remoteHostName;
	private String remoteUserName;
	private String privateKeyFilePath;	
	private String sshCommandPrefix;
	private String copyFileCommand;
	private String copyRemoteFilePrefix;
	private String workingDirectory = "";

	/*
	 * The local host name is determined and compared
	 * to the remote host name. If they are the same, 
	 * then commandPrefix is constructed so commands are
	 * executed locally, instead of via ssh. Doing this
	 * allows the same code to be used on different machines
	 * without having to worry about the use of this class
	 * causing trouble.
	 */
	private String localHostName;

	/*
	 * A lot of my commands will be executing java programs.
	 * As a convenience this class has special executeJava
	 * methods. The javaPath is the absolute path on the 
	 * remote host to the java executable.
	 */
	private String javaPath = ExecUtilities.getCommandPath("java");
	private String javaCommandPrefix = javaPath;
	
	public static void main(String[] args) {
		String host = "nob.vbi.vt.edu";
		String user = "enordber";
		String keyPath = "/Users/enordber/.ssh/id_rsa";
		
		
		RemoteHost remoteHost = new RemoteHost(host, user, keyPath);
		remoteHost.disableStrictHostKeyChecking();
		String cmd = "ls RAxML_info.raxmlRun0";
		String[] stdout = remoteHost.executeCommand(cmd).getStdout();
		for(int i = 0; i < stdout.length; i++) {
			System.out.println(stdout[i]);
		}
		
		remoteHost.fileExists("RAxML_info.raxmlRun0");
	}
	
	/**
	 * Returns an instance of RemoteHost that represents the local host.
	 * This may be poorly-named, but it is convenient.
	 * @return
	 */
	public static RemoteHost getLocalHost() {
		RemoteHost r = null;
		String hostName = null;
		try {
			hostName =
				InetAddress.getLocalHost().getHostName();
		} catch (UnknownHostException e) {
			hostName = "unknown";
			e.printStackTrace();
		}	

		String username = System.getProperty("user.name");

		/*
		 * key path is not needed for localhost
		 */
		String keyPath = null;
		r = new RemoteHost(hostName, username, keyPath);
		r.setRemoteWorkingDirectory(System.getProperty("user.dir"));
		return r;
	}

	public RemoteHost(String remoteHost, String remoteUser,
			String privateKeyPath) {
		this.remoteHostName = remoteHost;
		this.remoteUserName = remoteUser;
		this.privateKeyFilePath = privateKeyPath;

		try {
			localHostName =
				InetAddress.getLocalHost().getHostName();
		} catch (UnknownHostException e) {
			localHostName = "unknown";
			e.printStackTrace();
		}	

		if(remoteHostName != null) {
			if(remoteHostName.equals(localHostName)) {
				this.sshCommandPrefix = "";
				this.copyFileCommand = "cp ";
				this.copyRemoteFilePrefix = "";
			} else {
				this.sshCommandPrefix = "ssh -l " + remoteUserName +
				" -i " + privateKeyFilePath + 
				" " + remoteHostName + " ";
				this.copyFileCommand = "scp -i " + privateKeyFilePath 
				+ " ";
				this.copyRemoteFilePrefix = remoteUserName + "@" +
				remoteHostName + ":";
			}
		}
	}
	
	void setRemoteHostName(String remoteHost) {
		this.remoteHostName = remoteHost;

		try {
			localHostName =
				InetAddress.getLocalHost().getHostName();
		} catch (UnknownHostException e) {
			localHostName = "unknown";
			e.printStackTrace();
		}	

		if(remoteHostName.equals(localHostName)) {
			this.sshCommandPrefix = "";
			this.copyFileCommand = "cp ";
			this.copyRemoteFilePrefix = "";
		} else {
			this.sshCommandPrefix = "ssh -l " + remoteUserName +
			" -i " + privateKeyFilePath + 
			" " + remoteHostName + " ";
			this.copyFileCommand = "scp -i " + privateKeyFilePath 
			+ " ";
			this.copyRemoteFilePrefix = remoteUserName + "@" +
			remoteHostName + ":";
		}
	}

	public boolean isLocalHost() {
		return remoteHostName.equals(localHostName);
	}

	public RemoteHost(String remoteHost, String remoteUser,
			String privateKeyPath, String remoteWorkingDirectory) {
		this(remoteHost, remoteUser, privateKeyPath);
		setRemoteWorkingDirectory(remoteWorkingDirectory);
	}

	/**
	 * Sets the path to the java executable on the 
	 * remote host. This needs to be set is the 
	 * executeJava methods will be used. Otherwise, 
	 * it is not needed.
	 * @param path
	 */
	public void setRemoteJavaPath(String path) {
		if(debug) {
			System.out.println("RemoteHost.setJavaPath(): " 
					+ path);
		}
		this.javaPath = path;
		this.javaCommandPrefix = sshCommandPrefix +
		javaPath + " ";
	}

	/**
	 * Set the working directory for the RemoteHost. This directory
	 * will be used for calls to executeCommandInWorkingDirectory(),
	 * executeJavaInWorkingDirectory(), copyToRemoteWorkingDirectory(), 
	 * and copyFromRemoteWorkingDirectory().
	 * 
	 * @param directory
	 */
	public void setRemoteWorkingDirectory(String directory) {
		workingDirectory = directory;
	}

	/**
	 * Returns the working directory for the remote host.
	 * 
	 * @return
	 */
	public String getRemoteWorkingDirectory() {
		return workingDirectory;
	}

	/**
	 * Executes the given command in the current working directory 
	 * (provided in the 4-argument constructor or set via 
	 * setWorkingDirectory()) on the remote host.
	 * 
	 * @param command
	 * @return
	 */
	public CommandResults executeCommandInWorkingDirectory(String command) {
		return executeCommand(workingDirectory, command);
	}

	/**
	 * Convenience method for executing java commands in the 
	 * current working directory (provided in the 4-argument 
	 * constructor or set via setWorkingDirectory()) on the 
	 * remote host.
	 * 
	 * @param command
	 * @return
	 */
	public CommandResults executeJavaInWorkingDirectory(String command) {
		return executeJava(workingDirectory, command);
	}

	/**
	 * Copies the specified local file to the working directory
	 * (provided in the 4-argument constructor or set via 
	 * setWorkingDirectory()) on the remote host.
	 * remote host.
	 * @param localFileName
	 * @param remoteFileName
	 * @return
	 */
	public CommandResults copyToRemoteWorkingDirectory(String localFileName,
			String remoteFileName) {
		remoteFileName = workingDirectory + "/" + remoteFileName;
		return copyToRemoteHost(localFileName, remoteFileName);
	}

	/**
	 * Copies a file in the workingDirectory of the remote host
	 * (provided in the 4-argument constructor or set via 
	 * setWorkingDirectory()) to the local host.
	 * 
	 * @param remoteFileName
	 * @param localFileName
	 * @return
	 */
	public CommandResults copyFromRemoteWorkingDirectory(String remoteFileName,
			String localFileName) {
		remoteFileName = workingDirectory + "/" + remoteFileName;
		return copyFromRemoteHost(remoteFileName, localFileName);

	}

	/**
	 * Executes the given command on the remote host.
	 * 
	 * @param command
	 * @return
	 */
	public CommandResults executeCommand(String command) {
		CommandResults r = null;
		command = sshCommandPrefix + command;
		r = ExecUtilities.exec(command);
		return r;
	}

	/**
	 * Changes to the specified directory on the remote
	 * host prior to executing the command.
	 * 
	 * @param directory
	 * @param command
	 * @return
	 */
	public CommandResults executeCommand(String directory, 
			String command) {
		CommandResults r = null;
		command = sshCommandPrefix + "cd " + 
		directory + " ; " +command;
		r = ExecUtilities.exec(command);
		return r;
	}

	/**
	 * Convenience for running java programs. Java
	 * programs can be run using the executeCommand 
	 * also, but this method doesn't require you to 
	 * provide the full java executable path each time.
	 * You must, however, set the java executable path
	 * with setRemoteJavaPath() prior to calling this 
	 * method.
	 * 
	 * @param command
	 * @return
	 */
	public CommandResults executeJava(String command) {
		CommandResults r = null;
		command = javaCommandPrefix + command;
		r = ExecUtilities.exec(command);
		return r;
	}

	/**
	 * Changes to the specified directory on the remote
	 * host prior to executing the command.
	 * @param directory
	 * @param command
	 * @return
	 */
	public CommandResults executeJava(String directory, String command) {
		CommandResults r = null;
		if(System.getProperty("user.dir").equals(directory)) {
			command = javaCommandPrefix + command;
		} else {
			command = sshCommandPrefix + "cd " + 
			directory + " ; " + javaPath + " " + command;
		}

		if(debug) {
			System.out.println("RemoteHost.executeJava() command:");
			System.out.println(command);
		}
		r = ExecUtilities.exec(command);
		return r;
	}

	/**
	 * Copies a file from the specified location on the 
	 * local machine to the specified location on the 
	 * remote host.
	 * 
	 * @param localFileName
	 * @param remoteFileName
	 * @return
	 */
	public CommandResults copyToRemoteHost(String localFileName,
			String remoteFileName) {
		CommandResults r = null;
		String command = copyFileCommand + localFileName + " " 
		+ copyRemoteFilePrefix + remoteFileName;
		r = ExecUtilities.exec(command);
		return r;
	}

	/**
	 * Copies a file from the specified location on the
	 * remote host to the specified location on the local
	 * machine.
	 * 
	 * @param remoteFileName
	 * @param localFileName
	 * @return
	 */
	public CommandResults copyFromRemoteHost(String remoteFileName,
			String localFileName) {
		CommandResults r = null;
		String command = copyFileCommand + copyRemoteFilePrefix 
		+ remoteFileName + " " + localFileName;
		r = ExecUtilities.exec(command);
		return r;
	}

	/**
	 * Returns the absolute path to the given executable. This method
	 * relies on the 'which' utility. Therefore, it only works on
	 * Unix-type systems that have 'which' at /usr/bin/which. 
	 */
	public String getCommandPath(String cmd) {
		String r = null;
		//look for command path in System Properties
		r = System.getProperty(cmd);

		if(r == null) {
			//use 'which' utility to look for command path
			String whichCmd = "/usr/bin/which " + cmd;
			CommandResults whichRes = executeCommand(whichCmd);

			String[] stdout = whichRes.getStdout();
			if(stdout.length > 0) {
				r = stdout[0];
			}
			try  {
			    System.setProperty(cmd, r);
			} catch(NullPointerException npe) {
				//this is throwing an exception sometimes when several threads
				//are trying to do this at the same time. Since it's it not 
				//critical, and the command is just being stored for future use
				//to avoid looking it up again, this exception can just be
				//ignored. It just means that the command may not be stored 
				//this time, so it will have to be looked up again.
			}
		}
		return r;
	}

	/**
	 * Return value of true indicates success in disabling 
	 * strict host key checking.
	 * @return
	 */
	public boolean disableStrictHostKeyChecking() {
		boolean r = false;
		System.out.println(">RemoteHost.disableStrictHostKeyChecking() ");
		//disable strict checking
		String disableDommand = sshCommandPrefix + 
		" -o StrictHostKeyChecking=no " + 
		"-i " + privateKeyFilePath + " " + remoteUserName 
		+ "@" +remoteHostName + " ; ls -la";
		String[] results = ExecUtilities.exec(disableDommand).getStdout();
		for(int i = 0; i < results.length; i++) {
			System.out.println(results[i]);
		}
		if(results != null && results.length > 0 
				&& results[0].startsWith("total")) {
			r = true;
		}
		
		System.out.println("<RemoteHost.disableStrictHostKeyChecking() ");
		return r;

	}
	
	public void enableStrictHostKeyChecking() {
		System.out.println(">RemoteHost.enableStrictHostKeyChecking() ");
		//disable strict checking
		String command = sshCommandPrefix + 
		" -o StrictHostKeyChecking=yes " + remoteHostName + " ; ls -la";
		ExecUtilities.exec(command);
		System.out.println("<RemoteHost.enableStrictHostKeyChecking() ");

		
	}

	public boolean fileExists(String fileName) {
		boolean r = false;
		String lsCommand = "ls " + fileName;
		String[] res = executeCommand(lsCommand).getStdout();
		if(debug) {
			System.out.println("results of " + lsCommand);
			for(int i = 0; i < res.length; i++) {
				System.out.println(res[i]);
			}
		}
		if(res.length > 0) {
			if(res[0].startsWith(fileName)) {
				r = true;
			}
		}
		if(debug) {
			System.out.println("RemoteHost.fileExists() "+ fileName 
					+ ": " + r);
		}
		return r;
	}

}
