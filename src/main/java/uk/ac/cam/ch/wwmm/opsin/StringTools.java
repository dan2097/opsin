package uk.ac.cam.ch.wwmm.opsin;

import java.nio.ByteBuffer;
import java.nio.CharBuffer;
import java.nio.charset.CharacterCodingException;
import java.nio.charset.Charset;
import java.nio.charset.CharsetDecoder;
import java.nio.charset.CharsetEncoder;
import java.util.ArrayList;
import java.util.List;

/**Static routines for string manipulation.
 * This is a specially tailored version of StringTools as found in OSCAR for use in OPSIN
 *
 * @author ptc24/dl387
 *
 */
public final class StringTools {

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
	/** Whitespace characters. */
	public static String whiteSpace = "\u0020\u0085\u00a0\u1680\u180e\u2000\u2001\u2002\u2003" +
			"\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u2028\u2029\u202f\u205f\u3000";


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
}
