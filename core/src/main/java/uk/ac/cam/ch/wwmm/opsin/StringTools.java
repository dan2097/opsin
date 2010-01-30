package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.regex.Pattern;

/**Static routines for string manipulation.
 * This is a specially tailored version of StringTools as found in OSCAR for use in OPSIN
 *
 * @author ptc24/dl387
 *
 */
public final class StringTools {

	private static final Pattern MATCH_WHITESPACE = Pattern.compile("\\s+");

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

	/**
	 * Converts a list of strings into a single string delimited by the given separator
	 *
	 * @param l A list of strings.
	 * @return The corresponding string.
	 */
	public static String stringListToString(List<String> l, String separator) {
		StringBuffer sb = new StringBuffer();
		for(int i=0;i<l.size();i++) {
			sb.append(l.get(i));
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

	/**Checks to see if placing one or more close brackets (normal, square or curly) on 
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
			if(bracketLevel < 0) return true;
		}
		return false;
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

	/**Converts a unicode string into ISO-8859-1, converting greek letters
	 * to their names, and difficult characters to underscore.
	 *
	 * @param s The string to convert.
	 * @return The converted string.
	 * @throws PreProcessingException
	 */
	public static String convertNonAsciiAndNormaliseRepresentation(String s) throws PreProcessingException {
		s = MATCH_WHITESPACE.matcher(s).replaceAll(" ");//normalise white space
		StringBuilder sb = new StringBuilder(s);
		for(int i=0;i<sb.length();i++) {
			char c = sb.charAt(i);
			if(c >= 128) {
				sb.replace(i, i+1, getReplacementForNonASCIIChar(c));//replace non ascii characters with hard coded ascii strings
			}
			else if (c ==96){
				sb.replace(i, i+1, "'");//replace back ticks with apostrophe
			}
		}
		return sb.toString();
	}

    private static String getReplacementForNonASCIIChar(char c) throws PreProcessingException {
        switch (c) {
            case '\u03b1': return "alpha";//greeks
            case '\u03b2': return "beta";
            case '\u03b3': return "gamma";
            case '\u03b4': return "delta";
            case '\u03b5': return "epsilon";
            case '\u03b6': return "zeta";
            case '\u03b7': return "eta";
            case '\u03b8': return "theta";
            case '\u03b9': return "iota";
            case '\u03ba': return "kappa";
            case '\u03bb': return "lambda";
            case '\u03bc': return "mu";
            case '\u03bd': return "nu";
            case '\u03be': return "xi";
            case '\u03bf': return "omicron";
            case '\u03c0': return "pi";
            case '\u03c1': return "rho";
            case '\u03c2': return "stigma";
            case '\u03c3': return "sigma";
            case '\u03c4': return "tau";
            case '\u03c5': return "upsilon";
            case '\u03c6': return "phi";
            case '\u03c7': return "chi";
            case '\u03c8': return "psi";
            case '\u03c9': return "omega";

            case '\u2018': return "'";//quotation marks and primes (map to apostrophe/s)
            case '\u2019': return "'";
            case '\u201B': return "'";
            case '\u2032': return "'";//primes
            case '\u2033': return "''";
            case '\u2034': return "'''";
            case '\u2057': return "''''";
            case '\u2035': return "'";//back primes
            case '\u2036': return "''";
            case '\u2037': return "'''";

            case '\u2010': return "-";//dashes, hyphens and the minus sign
            case '\u2011': return "-";
            case '\u2012': return "-";
            case '\u2013': return "-";
            case '\u2014': return "-";
            case '\u2015': return "-";
            case '\u2212': return "-";
            
            case '\u00A0': return " ";//Non-breaking spaces
            case '\u2007': return " ";
            case '\u202F': return " ";
            
            case '\u200b': return "";//zero width space
            case '\u200d': return "";//zero width joiner

            case '\uFEFF': return "";//BOM-found at the start of some UTF files

            default: throw new PreProcessingException("Unrecognised unicode character: " + c);
        }
    }

	/**Converts a string array to an ArrayList.
	 *
	 * @param array The array.
	 * @return The ArrayList.
	 */
	public static List<String> arrayToList(String [] array) {
		List<String> list = new ArrayList<String>();
        list.addAll(Arrays.asList(array));
		return list;
	}

	/**
	 * If a dash is the last character it is removed
	 * @param locantText
	 * @return
	 */
	public static String removeDashIfPresent(String locantText){
		if(locantText.endsWith("-")) {
			locantText = locantText.substring(0, locantText.length()-1);
		}
		return locantText;
	}

	/**
	 * Any primes at the end of the string are removed
	 * @param locantText
	 * @return
	 */
	public static String removePrimesIfPresent(String locantText){
		while(locantText.endsWith("'")) {
			locantText = locantText.substring(0, locantText.length()-1);
		}
		return locantText;
	}
}
