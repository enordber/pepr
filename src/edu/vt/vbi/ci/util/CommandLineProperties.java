package edu.vt.vbi.ci.util;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import edu.vt.vbi.ci.util.file.TextFile;

/**
 * Parses a String[] of arguments passed to a main method
 * as command-line arguments. The arguments are parsed to populate
 * this HashMap subclass.
 * 
 * A flag is preceded with a dash '-'.
 * If the following argument is also a flag, then the initial flag is
 * added, with HandyConstants.TRUE as the value.
 * 
 * If the following argument is not a flag, then it is added as the value.
 * 
 * If multiple non-flag arguments follow a flag, they are added as
 * a String[] as the value for that flag. If a flag appears multiple
 * times, the values are all added as a single String[] as the value
 * for that flag. The multiple values will be stored in the order they
 * occur in the argument list.
 * 
 * @author enordber
 *
 */
public class CommandLineProperties extends HashMap {

	private static final String dash = "-";
	private static final String TRUE = HandyConstants.TRUE;

	private HashMap commandToDescription;
	private String[] knownCommands;
	
	public CommandLineProperties() {
		commandToDescription = new HashMap();
		knownCommands = new String[0];
	}
	
	public CommandLineProperties(String[] args) {
		this();
		addArgs(args);
	}

	public void addArgs(String[] args) {
		String activeFlag = "";
		String value = "";

		if(args != null) {
			for(int i = 0; i < args.length; i++) {
				if(isFlag(args[i])) {
					activeFlag = args[i].substring(1);
					if(args.length == i+1 || isFlag(args[i+1])) {
						addValue(activeFlag, TRUE);
					}
				} else {
					value = args[i];
					addValue(activeFlag, value);
				}
			}
		}
	}

	public void addArgs(CommandLineProperties clp) {
		this.addArgs(clp.getArgs());
	}

	private boolean isFlag(String s) {
		boolean r = false;
		if(s != null) {
			r = s.startsWith(dash);
		}
		return r;
	}

	private void addValue(String key, String value) {
		//If this key does not exist, add value
		String[] values = getValues(key);
		if(values == null) {
			values = new String[]{value};
		} else {
			String[] newValues = new String[values.length+1];
			//add new value to element 0, so values added last will be first in
			//the list. This is to allow grabbing element 0 to always
			//give the last-added value
			System.arraycopy(values, 0, newValues, 1, values.length);
			newValues[0] = value;
			values = newValues;
		}
		put(key, values);
	}

	public String[] getValues(String arg) {
		String[] r = null;
		r = (String[])get(arg);
		return r;
	}

	public String[] getValues(String arg, String defaultValue) {
		String[] r = getValues(arg);
		if(r == null) {
			r = new String[]{defaultValue};
		}
		return r;
	}

	/**
	 * Returns a String[] that, if used as the starting value for
	 * a CommandLineProperties set, would produce a set identical 
	 * to this one. This is not necessarily identical to the one 
	 * that was used to initialize this CommandLineProperties.
	 * @return
	 */
	public String[] getArgs() {
		String[] r = null;
		String[] keys = new String[size()];
		keySet().toArray(keys);
		Arrays.sort(keys);
		ArrayList args = new ArrayList();
		for(int i = 0; i < keys.length; i++) {
			args.add(dash+keys[i]);
			String[] values = getValues(keys[i]);
			//values are returned in reverse order, because they are added
			//in reverse order. This allows the returned args array to be used
			//to recreate an equivalent CommandLineProperties
			for(int j = values.length-1; j >= 0; j--) {
				args.add(values[j]);
			}
		}
		r = new String[args.size()];
		args.toArray(r);
		return r;
	}
	
	public void setKnownCommand(String command, String description) {
		String[] newCommands = new String[knownCommands.length+1];
		System.arraycopy(knownCommands, 0, newCommands, 0, knownCommands.length);
		newCommands[knownCommands.length] = command;
		knownCommands = newCommands;
		
		commandToDescription.put(command, description);
	}
	
	public String getKnownCommandsAndDescriptions() {
		String r = null;
		StringBuffer sb = new StringBuffer();
		for(int i = 0; i < knownCommands.length; i++) {
			String command = knownCommands[i];
			String description = (String)commandToDescription.get(command);
			String line = "\t-" + command + "\t\t" + description;
			sb.append(line);
			sb.append("\n");
		}
		r = sb.toString();
		return r;
	}
	
	public static CommandLineProperties loadFromFile(String fileName) throws IOException {
		CommandLineProperties r = null;
		TextFile tf = new TextFile(fileName);
		String[] args = tf.getAllLines();
		r = new CommandLineProperties(args);
		return r;
	}

	public static void saveToFile(CommandLineProperties clp, String fileName) throws IOException {
		FileWriter fw = new FileWriter(fileName);
		String[] args = clp.getArgs();
		for(int i = 0; i < args.length; i++) {
			fw.write(args[i]+"\n");
		}
		fw.flush();
		fw.close();
	}
	
	
}