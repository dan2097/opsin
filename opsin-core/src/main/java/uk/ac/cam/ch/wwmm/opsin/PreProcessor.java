package uk.ac.cam.ch.wwmm.opsin;

import java.util.HashMap;
import java.util.Locale;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Takes a name:
 * strips leading/trailing white space
 * Normalises representation of greeks and some other characters
 * @author dl387
 *
 */
class PreProcessor {
	private static final Map<String, String> DOTENCLOSED_TO_DESIRED = new HashMap<>();
	private static final Map<String, String> XMLENTITY_TO_DESIRED = new HashMap<>();
	private static final Pattern ESTER_CHALCOGEN_PAIR_LOCANT = Pattern.compile("(?<![A-Za-z0-9-])(O|S|Se|Te)(O|S|Se|Te)(?=[,-])(?!-[a-z]*perox)");
	//not preceded by a symbol pair, a multiplier (dithioperoxo has identical atoms) or a first prefix of selenothioperoxo
	private static final Pattern UNORIENTED_PEROXO_INFIX = Pattern.compile("(?<!di)(?<!tri)(?<!-)(?<!seleno)(?<!telluro)(?<!thio)((?:seleno|telluro|thio)(?:thio|seleno|telluro)?perox)", Pattern.CASE_INSENSITIVE);
	private static final String ACID_SYMBOLS = "(?-i:(?<pair>(?<first>O|S|Se|Te)(?<last>O|S|Se|Te))(?:,\\k<pair>)*)-";
	private static final Pattern PEROXO_ACID_IN_BRACKET = Pattern.compile("([(\\[{])([a-z]*perox[a-z]*(?:ic|ous)) ?([)\\]}]) ?(?:" + ACID_SYMBOLS + ")?acid", Pattern.CASE_INSENSITIVE);
	private static final Pattern PEROXO_INFIX_IN_BRACKET = Pattern.compile("([(\\[{])([a-z]*perox[a-z]*o)([)\\]}])([a-z]+(?:ic|ous)) " + ACID_SYMBOLS + "acid", Pattern.CASE_INSENSITIVE);
	private static final Pattern PEROXO_ACID_OO = Pattern.compile("(perox[a-z]*(?:ic|ous)) (?-i:OO)-acid", Pattern.CASE_INSENSITIVE);
	private static final Pattern DOUBLE_BRACKETED_PEROXO = Pattern.compile("\\(\\(([^()]*perox[^()]*)\\)\\)", Pattern.CASE_INSENSITIVE);

	static {
		DOTENCLOSED_TO_DESIRED.put("a", "alpha");
		DOTENCLOSED_TO_DESIRED.put("b", "beta");
		DOTENCLOSED_TO_DESIRED.put("g", "gamma");
		DOTENCLOSED_TO_DESIRED.put("d", "delta");
		DOTENCLOSED_TO_DESIRED.put("e", "epsilon");
		DOTENCLOSED_TO_DESIRED.put("l", "lambda");
		DOTENCLOSED_TO_DESIRED.put("x", "xi");
		DOTENCLOSED_TO_DESIRED.put("alpha", "alpha");
		DOTENCLOSED_TO_DESIRED.put("beta", "beta");
		DOTENCLOSED_TO_DESIRED.put("gamma", "gamma");
		DOTENCLOSED_TO_DESIRED.put("delta", "delta");
		DOTENCLOSED_TO_DESIRED.put("epsilon", "epsilon");
		DOTENCLOSED_TO_DESIRED.put("zeta", "zeta");
		DOTENCLOSED_TO_DESIRED.put("eta", "eta");
		DOTENCLOSED_TO_DESIRED.put("lambda", "lambda");
		DOTENCLOSED_TO_DESIRED.put("xi", "xi");
		DOTENCLOSED_TO_DESIRED.put("omega", "omega");
		DOTENCLOSED_TO_DESIRED.put("fwdarw", "->");

		XMLENTITY_TO_DESIRED.put("alpha", "alpha");
		XMLENTITY_TO_DESIRED.put("beta", "beta");
		XMLENTITY_TO_DESIRED.put("gamma", "gamma");
		XMLENTITY_TO_DESIRED.put("delta", "delta");
		XMLENTITY_TO_DESIRED.put("epsilon", "epsilon");
		XMLENTITY_TO_DESIRED.put("zeta", "zeta");
		XMLENTITY_TO_DESIRED.put("eta", "eta");
		XMLENTITY_TO_DESIRED.put("lambda", "lambda");
		XMLENTITY_TO_DESIRED.put("xi", "xi");
		XMLENTITY_TO_DESIRED.put("omega", "omega");
	}

	/**
	 * Master method for PreProcessing
	 * @param chemicalName
	 * @return
	 * @throws PreProcessingException 
	 */
	static String preProcess(String chemicalName) throws PreProcessingException {
		chemicalName = chemicalName.trim();//remove leading and trailing whitespace
		if (chemicalName.length() == 0){
			throw new PreProcessingException("Input chemical name was blank!");
		}
		
		chemicalName = performMultiCharacterReplacements(chemicalName);
		if (containsIgnoreCase(chemicalName, "perox")) {
			chemicalName = movePeroxoAcidSymbolsIntoInfix(chemicalName);
			chemicalName = copyEsterChalcogenPairIntoPeroxoInfix(chemicalName);
		}
		chemicalName = StringTools.convertNonAsciiAndNormaliseRepresentation(chemicalName);
		return chemicalName;
	}

	private static String performMultiCharacterReplacements(String chemicalName) {
		StringBuilder sb = new StringBuilder(chemicalName.length());
		for (int i = 0, nameLength = chemicalName.length(); i < nameLength; i++) {
			char ch = chemicalName.charAt(i);
			switch (ch) {
			case '$':
				if (i + 1 < nameLength){
					char letter = chemicalName.charAt(i + 1);
					String replacement = getReplacementForDollarGreek(letter);
					if (replacement != null){
						sb.append(replacement);
						i++;
						break;
					}
				}
				sb.append(ch);
				break;
			case '.':
				//e.g. .alpha.
				String dotEnclosedString = getLowerCasedDotEnclosedString(chemicalName, i);
				String dotEnclosedReplacement = DOTENCLOSED_TO_DESIRED.get(dotEnclosedString);
				if (dotEnclosedReplacement != null){
					sb.append(dotEnclosedReplacement);
					i = i + dotEnclosedString.length() + 1;
					break;
				}
				sb.append(ch);
				break;
			case '&':
				{
				//e.g. &alpha;
				String xmlEntityString = getLowerCasedXmlEntityString(chemicalName, i);
				String xmlEntityReplacement = XMLENTITY_TO_DESIRED.get(xmlEntityString);
				if (xmlEntityReplacement != null){
					sb.append(xmlEntityReplacement);
					i = i + xmlEntityReplacement.length() + 1;
					break;
				}
				sb.append(ch);
				break;
				}
			case '\u2070':
			case '\u00B9':
			case '\u00B2':
			case '\u00B3':
			case '\u2074':
			case '\u2075':
			case '\u2076':
			case '\u2077':
			case '\u2078':
			case '\u2079'://unicode superscript digits after a digit e.g. 1\u2074 --> 1^4 (otherwise isotope notation e.g. \u00B2H)
				if (i > 0 && isSuperscriptDigit(chemicalName.charAt(i - 1))) {
					sb.append(superscriptDigitToDigit(ch));
				}
				else if (i > 0 && Character.isDigit(chemicalName.charAt(i - 1))) {
					sb.append('^').append(superscriptDigitToDigit(ch));
				}
				else {
					sb.append(ch);
				}
				break;
			case '<':
				if (chemicalName.regionMatches(true, i, "<sup>", 0, 5)) {//e.g. 1<sup>4</sup> --> 1^4
					sb.append('^');
					i = i + 4;
					break;
				}
				if (chemicalName.regionMatches(true, i, "</sup>", 0, 6)) {
					i = i + 5;
					break;
				}
				sb.append(ch);
				break;
			case 's':
			case 'S'://correct British spelling to the IUPAC spelling
				if (chemicalName.regionMatches(true, i + 1, "ulph", 0, 4)){
					sb.append("sulf");
					i = i + 4;
					break;
				}
				sb.append(ch);
				break;
			default:
				sb.append(ch);
			}
		}
		return sb.toString();
	}

	/**
	 * e.g. ethane(thioperoxoic) OS-acid --> ethane(OS-thioperoxoic acid)
	 * The pair of symbols names the chain of the peroxo group from the acyl atom to the acidic atom
	 */
	private static String movePeroxoAcidSymbolsIntoInfix(String chemicalName) {
		Matcher m = PEROXO_ACID_IN_BRACKET.matcher(chemicalName);
		StringBuffer sb = new StringBuffer();
		while (m.find()) {
			String symbolPrefix = getPeroxoSymbolPrefix(m.group("first"), m.group("last"));
			if (symbolPrefix != null) {
				m.appendReplacement(sb, Matcher.quoteReplacement(m.group(1) + symbolPrefix + m.group(2) + " acid" + m.group(3)));
			}
		}
		m.appendTail(sb);
		chemicalName = sb.toString();
		m = PEROXO_INFIX_IN_BRACKET.matcher(chemicalName);
		sb = new StringBuffer();
		while (m.find()) {
			String symbolPrefix = getPeroxoSymbolPrefix(m.group("first"), m.group("last"));
			if (symbolPrefix != null) {
				m.appendReplacement(sb, Matcher.quoteReplacement(m.group(1) + symbolPrefix + m.group(2) + m.group(3) + m.group(4) + " acid"));
			}
		}
		m.appendTail(sb);
		chemicalName = DOUBLE_BRACKETED_PEROXO.matcher(sb.toString()).replaceAll("($1)");
		return PEROXO_ACID_OO.matcher(chemicalName).replaceAll("$1 acid");
	}

	private static String getPeroxoSymbolPrefix(String firstSymbol, String lastSymbol) {
		if (firstSymbol == null) {
			return "";
		}
		return firstSymbol.equals(lastSymbol) ? null : firstSymbol + lastSymbol + "-";
	}

	/**
	 * e.g. OS-methyl ethane(thioperoxoate) --> OS-methyl ethane(OS-thioperoxoate)
	 * The pair of symbols of an ester locant states which atom of the peroxy group is replaced and which carries the organyl
	 */
	private static String copyEsterChalcogenPairIntoPeroxoInfix(String chemicalName) {
		Matcher m = ESTER_CHALCOGEN_PAIR_LOCANT.matcher(chemicalName);
		String pair = null;
		while (m.find()) {
			if (m.group(1).equals(m.group(2))) {
				continue;
			}
			String foundPair = m.group(1) + m.group(2);
			if (pair != null && !pair.equals(foundPair)) {
				return chemicalName;
			}
			pair = foundPair;
		}
		if (pair == null) {
			return chemicalName;
		}
		return UNORIENTED_PEROXO_INFIX.matcher(chemicalName).replaceAll(Matcher.quoteReplacement(pair) + "-$1");
	}

	private static boolean containsIgnoreCase(String chemicalName, String text) {
		for (int i = 0, last = chemicalName.length() - text.length(); i <= last; i++) {
			if (chemicalName.regionMatches(true, i, text, 0, text.length())) {
				return true;
			}
		}
		return false;
	}

	private static boolean isSuperscriptDigit(char ch) {
		return ch == '\u2070' || ch == '\u00B9' || ch == '\u00B2' || ch == '\u00B3' || (ch >= '\u2074' && ch <= '\u2079');
	}

	private static char superscriptDigitToDigit(char ch) {
		switch (ch) {
		case '\u2070': return '0';
		case '\u00B9': return '1';
		case '\u00B2': return '2';
		case '\u00B3': return '3';
		default: return (char) ('4' + (ch - '\u2074'));
		}
	}

	private static String getLowerCasedDotEnclosedString(String chemicalName, int indexOfFirstDot) {
		int end = -1;
		int limit = Math.min(indexOfFirstDot + 9, chemicalName.length());
		for (int j = indexOfFirstDot + 1; j < limit; j++) {
			if (chemicalName.charAt(j) == '.'){
				end = j;
				break;
			}
		}
		if (end > 0){
			return chemicalName.substring(indexOfFirstDot + 1, end).toLowerCase(Locale.ROOT);
		}
		return null;
	}
	
	private static String getLowerCasedXmlEntityString(String chemicalName, int indexOfAmpersand) {
		int end = -1;
		int limit = Math.min(indexOfAmpersand + 9, chemicalName.length());
		for (int j = indexOfAmpersand + 1; j < limit; j++) {
			if (chemicalName.charAt(j) == ';'){
				end = j;
				break;
			}
		}
		if (end > 0){
			return chemicalName.substring(indexOfAmpersand + 1, end).toLowerCase(Locale.ROOT);
		}
		return null;
	}
	
	private static String getReplacementForDollarGreek(char ch) {
		switch (ch) {
		case 'a' :
			return "alpha";
		case 'b' :
			return "beta";
		case 'g' :
			return "gamma";
		case 'd' :
			return "delta";
		case 'e' :
			return "epsilon";
		case 'l' :
			return "lambda";
		default:
			return null;
		}
	}

}
