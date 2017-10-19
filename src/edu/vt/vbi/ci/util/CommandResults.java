package edu.vt.vbi.ci.util;

/**
 * This class stores and provides access to
 * the results of an executed command. The
 * results consist of the STDOUT output and
 * the STDERR output. Both are stored and made
 * available as String arrays.
 * 
 * @author enordber
 *
 */
public class CommandResults {

	private String[] stdout;
	private String[] stderr;
    private int rc;
	
    public CommandResults(String[] stdout, String[] stderr, int rc) {
		this.stdout = stdout;
		this.stderr = stderr;
		this.rc = rc;
	}
    public int getRc() { return rc; }
	
	public String[] getStdout() {
		return stdout;
	}
	
	public String[] getStderr() {
		return stderr;
	}
}
