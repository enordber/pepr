package edu.vt.vbi.ci.util.file;

public class FastaUtilities {
	/**
	 * Searches through each title for something that
	 *  looks like a taxon name. Returns an array containing
	 *  these possible taxon names for each title, if any are
	 *  found.
	 *  
	 * @param titles
	 * @return
	 */
	public static String[] getTaxaFromTitles(String[] titles) {
		String[] taxa = new String[titles.length];
		String openSB = " [";
		String closeSB = "]";
		String[] forbiddenChars = new String[]{
				" ",
				"\\(",
				"\\)",
				":",
				",",
				"\\[",
				"\\]"
		};
		for(int i = 0; i < titles.length; i++) {
			taxa[i] = getTaxonFromTitle(titles[i]);

			//some of the taxon names contain extra stuff,
			//so try to remove that here
			int pipeIndex = taxa[i].indexOf("|");
			if(pipeIndex > -1) {
				taxa[i] = taxa[i].substring(0, pipeIndex).trim();
			}

			//replace forbidden characters with underscore
			for(int j = 0; j < forbiddenChars.length; j++) {
				taxa[i] = taxa[i].replaceAll(forbiddenChars[j], "_");						
				taxa[i] = taxa[i].replaceAll("__", "_");						
			}
		}
		
		return taxa;
	}
	
	public static String getTaxonFromTitle(String title) {
		String r = null;
		char openSB = '[';
		char closeSB = ']';
		//some titles have nested square brackets inside the taxon name. 
		//this is an alternate parsing method trying to deal with this.
		
		//find the last closeSB
		char[] titleChars = title.toCharArray();
		int lastCloseSB = -1;
		for(int i = titleChars.length-1; i > 0; i--) {
			if(titleChars[i] == closeSB) {
				lastCloseSB = i;
				break;
			}
		}

		//Start from the last closeSB then look for the matching open bracket.
		int openSBToIgnore = 0;
		int matchingOpenIndex = -1;
		for(int i = lastCloseSB-1; i >= 0; i--) {
			if(titleChars[i] == closeSB) {
				openSBToIgnore++;
			} else if(titleChars[i] == openSB) {
				if(openSBToIgnore == 0) {
					matchingOpenIndex = i;
					break;
				} else {
					openSBToIgnore--;
				}
			}			
		}
		
		if(matchingOpenIndex > -1 && lastCloseSB > -1) {
			r = title.substring(matchingOpenIndex+1, lastCloseSB);
		} else {
			//no matching square brackets found, so return full title, 
			//after removing ">" if it is present
			if(title.startsWith(">")) {
			    r = title.substring(1);
			} else {
				r = title;
			}
		}
		return r;
	}

}
