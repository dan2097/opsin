package uk.ac.cam.ch.wwmm.ptclib.string;

import java.io.UnsupportedEncodingException;
import java.net.URLEncoder;
import java.nio.ByteBuffer;
import java.nio.CharBuffer;
import java.nio.charset.CharacterCodingException;
import java.nio.charset.Charset;
import java.nio.charset.CharsetDecoder;
import java.nio.charset.CharsetEncoder;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**Static routines for string manipulation.
 * 
 * @author ptc24
 *
 */
public final class StringTools {
	
	//public static String base62 = "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
	/** Lowercase Greek Unicode characters */
	public static String lowerGreek = "\u03b1\u03b2\u03b3\u03b4\u03b5\u03b6\u03b7\u03b8\u03b9\u03ba\u03bb\u03bc\u03bd\u03be\u03bf\u03c0\u03c1\u03c2\u03c3\u03c4\u03c5\u03c6\u03c7\u03c8\u03c9";
	/** Quotation marks of various Unicode forms */
	public static String quoteMarks = "\"'\u2018\u2019\u201A\u201B\u201C\u201D\u201E\u201F";
	/** Hyphens, dashes and the like */
	public static String hyphens = "-\u2010\u2011\u2012\u2013\u2014\u2015";
	/** A regex fragment for any hyphen or other dash */
	public static String hyphensRe = "(?:-|\u2010|\u2011|\u2012|\u2013|\u2014|\u2015)";
	/** Apostrophes, backticks, primess etc */
	public static String primes = "'`\u2032\u2033\u2034";
	/** The en dash */
	public static String enDash = "\u2013";
	/** The em dash */
	public static String emDash = "\u2014";
	/** Three dots at mid level. Commmonly used for hydrogen bonds. */
	public static String midElipsis = "\u22ef";
	/** Less than, greater than, equals, and other related characters. */
	public static String relations = "=<>\u2260\u2261\u2262\u2263\u2264\u2265\u2266\u2267\u2268" +
			"\u2269\u226a\u226b";
	/** Tests for the presence of two adjacent lowercase letters. */
	public static Pattern twoLowerPattern = Pattern.compile(".*[a-z][a-z].*");
	/** Whitespace characters. */
	public static String whiteSpace = "\u0020\u0085\u00a0\u1680\u180e\u2000\u2001\u2002\u2003" +
			"\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u2028\u2029\u202f\u205f\u3000";
		
	private static Pattern optionalPattern = Pattern.compile("\\(([0-9 ]+)\\)\\?");
	
	/**Finds the first uppercase letter in a name that might sensibly be 
	 * converted to lowercase.
	 */
	public static Pattern firstLowerCaseable = Pattern.compile("(^|[^a-z])([A-Z][a-z][a-z])");
	
	/**Removes junk from a string. <br>
	 * Where junk = initial openbracket, opensquarebracket, <br>
	 * terminal perion, comma, semicolon, colon, pling, query, 
	 * closebracket, closesquarebracket
	 * 
	 * @param s The string to dejunk.
	 * @return The dejunked string.
	 */
	public static String scrubWord(String s) {
		Matcher m = Pattern.compile("^[\\(\\[]*(.*?)[\\.,;:!\\?\\)\\]]*$").matcher(s);
		if(m.find()) {
			String txt = m.group(1);
			txt =  txt.replaceAll("\u00AD", "");
			return txt;
		} else {
			return s;
		}
	}
	
	/**Removes the letter "s" from the end of a string, if present.
	 * 
	 * @param s The string.
	 * @return The potentially modified string.
	 */
	public static String removeTerminalS(String s) {
		if(s.endsWith("s")) return s.substring(0, s.length()-1);
		else return s;
	}

	/**Counts the number of open brackets in a string.
	 * 
	 * @param s The string.
	 * @return The number of open brackets.
	 */
	public static int countOpenBrackets(String s) {
		int c = 0;
		for(int i=0;i<s.length();i++) {
			if(s.charAt(i) == '(') c++;
		}
		return c;
	}
	
	/**Converts a list of characters into a string.
	 * 
	 * @param l A list of characters.
	 * @return The corresponding string.
	 */
	public static String charListToString(List<Character> l) {
		StringBuffer sb = new StringBuffer();
		for(char c : l) {
			sb.append(c);
		}
		return sb.toString();
	}

	/**Converts a list of strings into a string.
	 * 
	 * @param l A list of characters.
	 * @return The corresponding string.
	 */
	public static String stringListToString(List<String> l) {
		StringBuffer sb = new StringBuffer();
		for(String s : l) {
			sb.append(s);
		}
		return sb.toString();
	}

	/**Converts a list of objects into a string.
	 * 
	 * @param l A list of characters.
	 * @return The corresponding string.
	 */
	public static String objectListToString(List l, String separator) {
		StringBuffer sb = new StringBuffer();
		for(int i=0;i<l.size();i++) {
			sb.append(l.get(i).toString());
			if(separator != null && i < l.size()-1) sb.append(separator);
		}
		return sb.toString();
	}
	
	/**Converts a string to a list of characters.
	 * 
	 * @param s A string.
	 * @return The corresponding list of characters.
	 */
	public static List<Character> stringToList(String s) {
		List<Character> cl = new ArrayList<Character>();
		for(int i=0;i<s.length();i++) {
			cl.add(s.charAt(i));
		}
		return cl;
	}
	
	/**Produce repetitions of a string. Eg. HelloWorld * 2 = HelloWorldHelloWorld.
	 * 
	 * @param s The string to multiply.
	 * @param n The number of times to multiply it.
	 * @return The multiplied string.
	 */
	public static String multiplyString(String s, int n) {
		StringBuffer sb = new StringBuffer();
		for(int i=0;i<n;i++) {
			sb.append(s);
		}
		return sb.toString();
	}
	
	/**Checks to see if placing an open bracket (normal, square or curly) on the
	 * front of the string would cause the string to have balanced brackets.
	 * 
	 * @param s The string to test.
	 * @return The result of the test.
	 */
	public static boolean isLackingOpenBracket(String s) {
		int bracketLevel = 0;
		for(int i=0;i<s.length();i++) {
			char c = s.charAt(i);
			if(c == '(' || c == '[' || c == '{') bracketLevel++;
			if(c == ')' || c == ']' || c == '}') bracketLevel--;
			if(bracketLevel == -1) return true;
		}
		return false;
	}

	/**Checks to see if placing an close bracket (normal, square or curly) on 
	 * the end of the string would cause the string to have balanced brackets.
	 * 
	 * @param s The string to test.
	 * @return The result of the test.
	 */
	public static boolean isLackingCloseBracket(String s) {
		int bracketLevel = 0;
		for(int i=s.length()-1;i>=0;i--) {
			char c = s.charAt(i);
			if(c == '(' || c == '[' || c == '{') bracketLevel--;
			if(c == ')' || c == ']' || c == '}') bracketLevel++;
			if(bracketLevel == -1) return true;
		}
		return false;
	}
	
	/**Checks to see whether the brackets in the string are balanced. Note that
	 * this does not distinguish between normal, square and curly brackets.
	 * Furthermore, "a)(b" does not count as balanced.
	 * 
	 * @param s The string to test.
	 * @return The result of the test.
	 */
	public static boolean bracketsAreBalanced(String s) {
		int bracketLevel = 0;
		for(int i=0;i<s.length();i++) {
			char c = s.charAt(i);
			if(c == '(' || c == '[' || c == '{') bracketLevel++;
			if(c == ')' || c == ']' || c == '}') bracketLevel--;
			if(bracketLevel == -1) return false;
		}
		if(bracketLevel == 0) return true;
		return false;
	}
	
	/**Checks to see if the string is surrounded by a pair of (normal, square
	 * or curly) brackets, and that those brackets match each other.
	 * 
	 * @param s The string to test.
	 * @return The result of the test.
	 */
	public static boolean isBracketed(String s) {
		if(s == null || s.length() < 3) return false;
		char first = s.charAt(0);
		char last = s.charAt(s.length()-1);
		if(!((first == '(' && last == ')') || (first == '[' && last == ']') || (first == '{' && last == '}'))) return false;
		if(!bracketsAreBalanced(s.substring(1, s.length()-1))) return false;
		return true;
	}
	
	/**Joins a collection of strings into a single string.
	 * 
	 * @param strings The strings to join together.
	 * @param separator The separator to use.
	 * @return The resulting string.
	 */
	public static String collectionToString(Collection<String> strings, String separator) {
		Object [] objArray = strings.toArray();
		String [] stringArray = new String[objArray.length];
		for(int i=0;i<objArray.length;i++) {
			stringArray[i] = (String)objArray[i];
		}
		return arrayToString(stringArray, separator);
	}

	/**Joins an array of strings into a single string.
	 * 
	 * @param stringArray The strings to join together.
	 * @param separator The separator to use.
	 * @return The resulting string.
	 */
	public static String arrayToString(String [] stringArray, String separator) {
		StringBuffer sb = new StringBuffer();
		for(int i=0;i<stringArray.length-1;i++) {
			sb.append(stringArray[i]);
			sb.append(separator);
		}
		sb.append(stringArray[stringArray.length-1]);
		return sb.toString();
	}
	
	/**Finds the position of the last hyphen in the string.
	 * 
	 * @param s The string to test.
	 * @return The position of the last hyphen, or -1 if there isn't one.
	 */
	public static int lastIndexOfHyphen(String s) {
		int idx = -1;
		for(int i=0;i<hyphens.length();i++) {
			idx = Math.max(idx, s.lastIndexOf(hyphens.codePointAt(i)));
		}
		return idx;
	}
	
	/**URLEncodes a long string, adding newlines if necessary.
	 * 
	 * @param s The string to URLEncode.
	 * @return The URLEncoded string.
	 */
	public static String urlEncodeLongString(String s) {
		StringBuffer sb = new StringBuffer();
		int chunks = s.length() / 50;
		for(int i=0;i<chunks;i++) {
				sb.append(urlEncodeUTF8NoThrow(s.substring(i*50, (i+1)*50)));
				sb.append("\n");
		}
		sb.append(urlEncodeUTF8NoThrow(s.substring(chunks*50)));
		return sb.toString();
	}
	
	/**URLEncodes a string for UTF-8. This should not throw an exception as
	 * UTF-8 is unlikely to be an unsupported encoding.
	 * 
	 * @param s The string to encode.
	 * @return The encoded string.
	 */
	public static String urlEncodeUTF8NoThrow(String s) {
		try {
			return(URLEncoder.encode(s, "UTF-8"));
		} catch (UnsupportedEncodingException e) {
			throw new Error("Wot no UTF-8 for URLEncode?");
		}		
	}
	
	/*public static String intToBase62(int i) {
		String s = "";
		String start = "";
		if(i<0) {
			start = "-";
			i = -i;
		}
		do {
			int j = i % 62;
			s = base62.substring(j, j+1) + s;
			i = i / 62;
		} while (i > 0);
		return start+s;
	}*/
		
	/**Replace whitespace with a single space, remove soft hyphens, and convert
	 * (whitespace-delimited) tokens to lowercase if two adjacent lowercase
	 * characters are detected.
	 * 
	 * @param name The name to convert.
	 * @return The normalised name.
	 */
	public static String normaliseName(String name) {
		String [] subStrings = name.split("\\s+");
		for(int i=0;i<subStrings.length;i++) {
			if(twoLowerPattern.matcher(subStrings[i]).matches()) {
				subStrings[i] = subStrings[i].toLowerCase();
			}
			subStrings[i] = subStrings[i].replaceAll("\u00ad", "");
		}
		if(subStrings.length == 0) return "";
		if(subStrings.length == 1) return subStrings[0];
		return arrayToString(subStrings, " ");
	}
	
	/**As normalise name, but with a better heuristic for deciding when to
	 * convert to lowercase.
	 * 
	 * @param name The name to convert.
	 * @return The normalised name.
	 */
	public static String normaliseName2(String name) {
		String [] subStrings = name.split("\\s+");
		for(int i=0;i<subStrings.length;i++) {
			subStrings[i] = subStrings[i].replaceAll("\u00ad", "");
			Matcher m = firstLowerCaseable.matcher(subStrings[i]);
			if(m.find()) {
				subStrings[i] = subStrings[i].substring(0,m.start()) + 
					subStrings[i].substring(m.start(),m.end()).toLowerCase() + 
					subStrings[i].substring(m.end());
			}
		}
		if(subStrings.length == 0) return "";
		if(subStrings.length == 1) return subStrings[0];
		return arrayToString(subStrings, " ");
	}
	
	/**Converts a unicode string into ISO-8859-1, converting greek letters
	 * to their names, and difficult characters to underscore.
	 * 
	 * @param s The string to convert.
	 * @return The converted string.
	 */
	public static String unicodeToLatin(String s) {
		boolean hasUnicode = false;
		for(int i=0;i<s.length();i++) {
			char c = s.charAt(i);
			if(c >= 128) {
				hasUnicode = true;
				break;
			}
		}
		if(!hasUnicode) return s;
		s = s.replace("\u03b1", "alpha");
		s = s.replace("\u03b2", "beta");
		s = s.replace("\u03b3", "gamma");
		s = s.replace("\u03b4", "delta");
		s = s.replace("\u03b5", "epsilon");
		s = s.replace("\u03b6", "zeta");
		s = s.replace("\u03b7", "eta");
		s = s.replace("\u03b8", "theta");
		s = s.replace("\u03b9", "iota");
		s = s.replace("\u03ba", "kappa");
		s = s.replace("\u03bb", "lambda");
		s = s.replace("\u03bc", "mu");
		s = s.replace("\u03bd", "nu");
		s = s.replace("\u03be", "xi");
		s = s.replace("\u03bf", "omicron");
		s = s.replace("\u03c0", "pi");
		s = s.replace("\u03c1", "rho");
		s = s.replace("\u03c2", "stigma");
		s = s.replace("\u03c3", "sigma");
		s = s.replace("\u03c4", "tau");
		s = s.replace("\u03c5", "upsilon");
		s = s.replace("\u03c6", "phi");
		s = s.replace("\u03c7", "chi");
		s = s.replace("\u03c8", "psi");
		s = s.replace("\u03c9", "omega");
		
		Charset charset = Charset.forName("ISO-8859-1");
	    CharsetDecoder decoder = charset.newDecoder();
	    CharsetEncoder encoder = charset.newEncoder();	    
	    try {
	        ByteBuffer bbuf = encoder.encode(CharBuffer.wrap(s));
	        CharBuffer cbuf = decoder.decode(bbuf);
	        return cbuf.toString();
	    } catch (CharacterCodingException e) {
	    	s = s.replaceAll("[^A-Za-z0-9_+-]", "_");
		    try {
		    	ByteBuffer bbuf = encoder.encode(CharBuffer.wrap(s));
		    	CharBuffer cbuf = decoder.decode(bbuf);
		    	return cbuf.toString();
		    } catch (CharacterCodingException ee) {
		    	return null;
		    }
	    }
	}
	
	/**Tests to see whether one string could be an acronym of the other.
	 * 
	 * @param potentialAcroymn The potential acronym.
	 * @param reference The potential acronym expansion.
	 * @return Whether there is a match.
	 */
	public static boolean testForAcronym(String potentialAcroymn, String reference) {
		potentialAcroymn = potentialAcroymn.toLowerCase();
		reference = reference.toLowerCase();
		int refpoint = 0;
		for(int i=0;i<potentialAcroymn.length();i++) {
			refpoint = reference.indexOf(potentialAcroymn.charAt(i), refpoint);
			if(refpoint == -1) return false;
			refpoint++;
		}
		return true;
	}

	/**Sorts a list of strings, in reverse order of the values that they are 
	 * mapped to.
	 * 
	 * @param list The list to sort.
	 * @param map The mapping.
	 */
	public static void sortStringList(List<String> list, Map<String,? extends Comparable> map) {
		final Map<String,? extends Comparable> fmap = map;
		Collections.sort(list, Collections.reverseOrder(new Comparator<String>() {
			@SuppressWarnings("unchecked")
			public int compare(String o1, String o2) {
				// TODO Auto-generated method stub
				return fmap.get(o1).compareTo(fmap.get(o2));
			}
		}));
	}
	
	/**Extracts all of the keys from the map, and returns them sorted in 
	 * reverse order of the values.
	 * 
	 * @param map The map.
	 * @return The sorted list.
	 */
	public static List<String> getSortedList(Map<String,? extends Comparable> map) {
		List<String> list = new ArrayList<String>(map.keySet());
		sortStringList(list, map);
		return list;
	}
	
	/**Converts a string array to an ArrayList.
	 * 
	 * @param array The array.
	 * @return The ArrayList.
	 */
	public static List<String> arrayToList(String [] array) {
		List<String> list = new ArrayList<String>();
		for(int i=0;i<array.length;i++) {
			list.add(array[i]);
		}
		return list;
	}
	
	/**Merges two space-separated lists, removing duplicate items.
	 * 
	 * @param ssSet1 The first list.
	 * @param ssSet2 The second list.
	 * @return The combined list.
	 */
	public static String mergeSpaceSeparatedSets(String ssSet1, String ssSet2) {
		String [] array1 = ssSet1.split("\\s+");
		String [] array2 = ssSet2.split("\\s+");
		Set<String> outSet = new LinkedHashSet<String>();
		for(int i=0;i<array1.length;i++) {
			outSet.add(array1[i]);
		}
		for(int i=0;i<array2.length;i++) {
			outSet.add(array2[i]);
		}
		return collectionToString(outSet, " ");
	}
	
	/**Ensures that a string is no longer than the given length, if necessary
	 * by discarding the end.
	 * 
	 * @param s The string to shorten.
	 * @param maxlen The desired maximum length.
	 * @return The resulting string.
	 */
	public static String shorten(String s, int maxlen) {
		if(s.length() <= maxlen) return s;
		return s.substring(0, maxlen);
	}
	
	/**Turns a collection of strings into a single string, by sorting it
	 * and concatenating. Useful for hashing.
	 * 
	 * @param coll The collection of strings.
	 * @return The resulting concatenated string.
	 */
	public static String collectionToStableString(Collection<String> coll) {
		List<String> list = new ArrayList<String>(coll);
		Collections.sort(list);
		StringBuffer sb = new StringBuffer();
		for(String s : list) sb.append(s);
		return sb.toString();
	}

	/**Converts a string-to-string mapping into a string that is useful for
	 * hashing.
	 * 
	 * @param map The mapping.
	 * @return The resulting string.
	 */
	public static String mapToStableString(Map<String,String> map) {
		List<String> list = new ArrayList<String>(map.keySet());
		Collections.sort(list);
		StringBuffer sb = new StringBuffer();
		for(String s : list) sb.append(s + " -> " + map.get(s));
		return sb.toString();
	}
	
	/**Takes a space-separated list, and produces all of the possible strings
	 * that are subsets (including the whole set and the empty set) of that
	 * set.
	 * 
	 * @param ssList The space separated list.
	 * @return The possibilities.
	 */
	public static List<String> spaceSepListToSubLists(String ssList) {
		List<String> possibilities = new ArrayList<String>();
		possibilities.add("");
		String [] subStrs = ssList.split("\\s");
		for(int i=0;i<subStrs.length;i++) {
			for(String s : new ArrayList<String>(possibilities)) {
				if(s.length() == 0) {
					possibilities.add(subStrs[i]);
				} else {
					possibilities.add(s + " " + subStrs[i]);
				}
			}
		}
		return possibilities;
	}
	
	/*public static List<String> makeNGrams(List<List<String>> nGramParts) {
		List<String> nGrams = new ArrayList<String>();
		nGrams.add("");
		boolean addUnderscore = false;
		for(List<String> ngpl : nGramParts) {
			List<String> newNGrams = new ArrayList<String>();
			for(String nGram : nGrams) {
				if(addUnderscore) nGram = nGram+"_";
				for(String ngp : ngpl) {
					newNGrams.add(nGram + ngp);
				}
			}
			nGrams = newNGrams;
			addUnderscore = true;
		}
		//System.out.println(nGramParts);
		//System.out.println(nGrams);
		return nGrams;
	}*/
	
	
	/** Expands a string consiting of digits, whitespace and regex characters
	 * into a finite set of digits/whitespace only strings if possible.
	 * 
	 * @param regex The regex to expand.
	 * @return The strings that the regex can match.
	 */
	public static Set<String> expandRegex(String regex) {
		Set<String> results = new LinkedHashSet<String>();
		if(regex == null || regex.length() == 0) return results;
		if(regex.matches("[0-9 ]+")) {
			results.add(regex);
			return results;
		} else if(regex.matches("[0-9 ()?]+")) {
			Matcher m = optionalPattern.matcher(regex);
			if(m.find()) {
				String before = regex.substring(0, m.start());
				String middle = regex.substring(m.start(1), m.end(1));
				String after = regex.substring(m.end());
				results.addAll(expandRegex(before + after));
				results.addAll(expandRegex(before + middle + after));
			} else {
				results.add(regex);
			}
			return results;
		} else {
			results.add(regex);
			return results;			
		}
	}
}
