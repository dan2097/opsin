package uk.ac.cam.ch.wwmm.opsin;

import java.util.HashMap;
import java.util.Locale;
import java.util.Map;

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
