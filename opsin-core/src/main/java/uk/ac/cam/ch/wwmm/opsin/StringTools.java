package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import static uk.ac.cam.ch.wwmm.opsin.OpsinTools.*;

/**Static routines for string manipulation.
 * This is a specially tailored version of StringTools as found in OSCAR for use in OPSIN
 *
 * @author ptc24
 * @author dl387
 *
 */
public final class StringTools {

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
		for(int i=n-1;i>=0;i--) {
			sb.append(s);
		}
		return sb.toString();
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
			else if (c ==34){
				sb.replace(i, i+1, "''");//replace quotation mark with two primes
			}
			else if (c<=31){
				sb.replace(i, i+1, "");//ignore control characters
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
            
            case '\u1D05': return "D";//small capitals
            case '\u029F': return "L";

            case '\u00B1': return "+-";//plus minus symbol
            case '\u2213': return "-+";
            
            case '\u00C6': return "AE";//common ligatures
            case '\u00E6': return "ae";
            case '\u0152': return "OE";
            case '\u0153': return "oe";
            case '\u0132': return "IJ";
            case '\u0133': return "ij";
            case '\u1D6B': return "ue";
            case '\uFB00': return "ff";
            case '\uFB01': return "fi";
            case '\uFB02': return "fl";
            case '\uFB03': return "ffi";
            case '\uFB04': return "ffl";
            case '\uFB06': return "st";
            
            case '\u00E0': return "a";//diacritics
            case '\u00C0': return "A";
            case '\u00E1': return "a";
            case '\u00C1': return "A";
            case '\u00E2': return "a";
            case '\u00C2': return "A";
            case '\u00E3': return "a";
            case '\u00C3': return "A";
            case '\u00E4': return "a";
            case '\u00C4': return "A";
            case '\u00E5': return "a";
            case '\u00C5': return "A";
            case '\u00E7': return "c";
            case '\u00C7': return "C";
            case '\u00E8': return "e";
            case '\u00C8': return "E";
            case '\u00E9': return "e";
            case '\u00C9': return "E";
            case '\u00EA': return "e";
            case '\u00CA': return "E";
            case '\u00EB': return "e";
            case '\u00CB': return "E";
            case '\u00EC': return "i";
            case '\u00CC': return "I";
            case '\u00ED': return "i";
            case '\u00CD': return "I";
            case '\u00EE': return "i";
            case '\u00CE': return "I";
            case '\u00EF': return "i";
            case '\u00CF': return "I";
            case '\u00F2': return "o";
            case '\u00D2': return "O";
            case '\u00F3': return "o";
            case '\u00D3': return "O";
            case '\u00F4': return "o";
            case '\u00D4': return "O";
            case '\u00F5': return "o";
            case '\u00D5': return "O";
            case '\u00F6': return "o";
            case '\u00D6': return "O";
            case '\u00F9': return "u";
            case '\u00D9': return "U";
            case '\u00FA': return "u";
            case '\u00DA': return "U";
            case '\u00FB': return "u";
            case '\u00DB': return "U";
            case '\u00FC': return "u";
            case '\u00DC': return "U";
            case '\u00FD': return "y";
            case '\u00DD': return "Y";

            case '\u0115': return "e";
            case '\u0114': return "E";
            case '\u0117': return "e";
            case '\u0116': return "E";
            
            case '\u2018': return "'";//quotation marks and primes (map to apostrophe/s)
            case '\u2019': return "'";
            case '\u201B': return "'";
            case '\u201C': return "''";
            case '\u201D': return "''";
            case '\u2032': return "'";//primes
            case '\u2033': return "''";
            case '\u2034': return "'''";
            case '\u2057': return "''''";
            case '\u2035': return "'";//back primes
            case '\u2036': return "''";
            case '\u2037': return "'''";
            case '\u00B4': return "'";//accents
            case '\u02CA': return "'";
            case '\u0301': return "'";
            case '\u02DD': return "''";
            case '\u030B': return "''";

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
	 * Counts the number of primes at the end of a locant
	 * @param locantText
	 * @return
	 */
	public static int countTerminalPrimes(String locantText){
		int numberOfPrimes = 0;
		for(int k = locantText.length() -1; k>0; k--){
			if (locantText.charAt(k)=='\''){
				numberOfPrimes++;
			}
			else{
				break;
			}
		}
		return numberOfPrimes;
	}
	
	/**
	 * Tests if this string ends with the specified suffix ignoring case.
	 * @param str
	 * @param suffix
	 * @return
	 */
	public static boolean endsWithCaseInsensitive(String str, String suffix) {
		if (suffix.length() > str.length()) {
			return false;
		}
		int strOffset = str.length() - suffix.length();
		return str.regionMatches(true, strOffset, suffix, 0, suffix.length());
	}
}
