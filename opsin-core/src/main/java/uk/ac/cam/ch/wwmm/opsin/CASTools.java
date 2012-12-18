package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import static uk.ac.cam.ch.wwmm.opsin.OpsinTools.*;

/**
 * Tools for converting CAS nomenclature into IUPAC nomenclature.
 * @author dl387
 */
class CASTools {

	private static final Pattern matchCasCollectiveIndex = Pattern.compile("([\\[\\(\\{]([1-9][0-9]?[cC][iI][, ]?)+[\\]\\)\\}])+|[1-9][0-9]?[cC][iI]", Pattern.CASE_INSENSITIVE);
	private static final Pattern matchAcid = Pattern.compile("acid[\\]\\)\\}]*", Pattern.CASE_INSENSITIVE);
	private static final Pattern matchCommaSpace = Pattern.compile(", ");
	private static final Pattern matchCompoundWithPhrase = Pattern.compile("(compd\\. with|compound with|and) ", Pattern.CASE_INSENSITIVE);
	private static final Pattern matchFunctionalTermAllowingSubstituentPrefix = Pattern.compile("(amide|hydrazide|(thi|selen|tellur)?oxime|hydrazone|(iso)?(semicarbazone|thiosemicarbazone|selenosemicarbazone|tellurosemicarbazone)|imide|imine|semioxamazone)[\\]\\)\\}]*", Pattern.CASE_INSENSITIVE);
	
	/**
	 * Inverts a CAS name.
	 * Throws an exception is OPSIN is unable to determine whether something is a substituent or functional term
	 * or if something unexpected in a CAS name is encountered
	 * @param name
	 * @return
	 * @throws ParsingException
	 */
	static String uninvertCASName(String name, ParseRules parseRules) throws ParsingException {
		List<String> nameComponents = new ArrayList<String>(Arrays.asList(matchCommaSpace.split(name)));
		List<String> substituents = new ArrayList<String>();
		List<String> seperateWordSubstituents = new ArrayList<String>();
		List<String> functionalTerms = new ArrayList<String>();

		String parent = nameComponents.get(0);
		String[] parentNameParts = MATCH_SPACE.split(parent);
		if (parentNameParts.length != 1) {
			if (matchCasCollectiveIndex.matcher(parentNameParts[parentNameParts.length - 1]).matches()) {//CAS collective index description should be ignored
				StringBuilder parentSB = new StringBuilder();
				for (int i = 0; i < parentNameParts.length - 1; i++) {
					parentSB.append(parentNameParts[i]);
				}
				parent = parentSB.toString();
				parentNameParts = MATCH_SPACE.split(parent);
			}
			for (int i = 1; i < parentNameParts.length; i++) {
				if (!matchAcid.matcher(parentNameParts[i]).matches()) {
					ParseRulesResults results = parseRules.getParses(parentNameParts[i]);
					List<ParseTokens> parseTokens = results.getParseTokensList();
					if (parseTokens.isEmpty()) {
						throw new ParsingException("Invalid CAS name. Parent compound was followed by an unexpected term");
					}
				}
			}
		}
		boolean addedBracket = false;
		boolean esterEncountered = false;
		for (int i = 1; i < nameComponents.size(); i++) {
			String nameComponent = nameComponents.get(i);
			Matcher m = matchCompoundWithPhrase.matcher(nameComponent);
			boolean compoundWithcomponent = false;
			if (m.lookingAt()) {
				nameComponent = nameComponent.substring(m.group().length());
				compoundWithcomponent = true;
			}
			String[] components = MATCH_SPACE.split(nameComponents.get(i));
			for (String component : components) {
				if (compoundWithcomponent) {
					functionalTerms.add(component);
					continue;
				}
				if (component.endsWith("-")) {
					Character missingCloseBracket = missingCloseBracketCharIfApplicable(component);
					if (missingCloseBracket !=null) {
						if (addedBracket) {
							throw new ParsingException("Close bracket appears to be missing");
						}
						parent += missingCloseBracket;
						addedBracket = true;
					}
					substituents.add(component);
				} else {
					ParseRulesResults results = parseRules.getParses(component);
					List<ParseTokens> parseTokens = results.getParseTokensList();
					if (parseTokens.size() > 0) {
						List<ParseWord> parseWords = WordTools.splitIntoParseWords(parseTokens, component);

						List<ParseTokens> firstParseWordTokens = parseWords.get(0).getParseTokens();
						WordType firstWordType = OpsinTools.determineWordType(firstParseWordTokens.get(0).getAnnotations());
						for (int j = 1; j < firstParseWordTokens.size(); j++) {
							if (!firstWordType.equals(OpsinTools.determineWordType(firstParseWordTokens.get(j).getAnnotations()))) {
								throw new ParsingException(component + "can be interpreted in multiple ways. For the sake of precision OPSIN has decided not to process this as a CAS name");
							}
						}
						
						if (parseWords.size() ==1){
							if (firstWordType.equals(WordType.functionalTerm)) {
								if (component.equalsIgnoreCase("ester")) {
									if (esterEncountered) {
										throw new ParsingException("ester formation was mentioned more than once in CAS name!");
									}
									parent = uninvertEster(parent);
									esterEncountered = true;
								} else {
									functionalTerms.add(component);
								}
							} else if (firstWordType.equals(WordType.substituent)) {
								seperateWordSubstituents.add(component);
							} else if (firstWordType.equals(WordType.full)) {
								if (StringTools.endsWithCaseInsensitive(component, "ate") || StringTools.endsWithCaseInsensitive(component, "ite")//e.g. Piperazinium, 1,1-dimethyl-, 2,2,2-trifluoroacetate hydrochloride
												|| component.equalsIgnoreCase("hydrofluoride") || component.equalsIgnoreCase("hydrochloride") || component.equalsIgnoreCase("hydrobromide") || component.equalsIgnoreCase("hydroiodide")) {
									functionalTerms.add(component);
								} else {
									throw new ParsingException("Unable to interpret: " + component + " (as part of a CAS index name)- A full word was encountered where a substituent or functionalTerm was expected");
								}
							}
						}
						else if (parseWords.size() == 2 && firstWordType.equals(WordType.substituent)) {
							//could be something like O-methyloxime which is parsed as [O-methyl] [oxime]
							List<ParseTokens> secondParseWordTokens = parseWords.get(1).getParseTokens();
							WordType secondWordType = OpsinTools.determineWordType(secondParseWordTokens.get(0).getAnnotations());
							for (int j = 1; j < secondParseWordTokens.size(); j++) {
								if (!secondWordType.equals(OpsinTools.determineWordType(secondParseWordTokens.get(j).getAnnotations()))) {
									throw new ParsingException(component + "can be interpreted in multiple ways. For the sake of precision OPSIN has decided not to process this as a CAS name");
								}
							}
							if (secondWordType.equals(WordType.functionalTerm) && 
									matchFunctionalTermAllowingSubstituentPrefix.matcher(parseWords.get(1).getWord()).matches()){
								functionalTerms.add(component);
							}
							else{
								throw new ParsingException("Unrecognised CAS index name form, could have a missing space?");
							}
						}
						else {
							throw new ParsingException("Unrecognised CAS index name form");
						}
					} else {
						if (!matchCasCollectiveIndex.matcher(component).matches()) {//CAS collective index description should be ignored
							throw new ParsingException("Unable to interpret: " + component + " (as part of a CAS index name)");
						}
					}
				}
			}
		}
		StringBuilder casName = new StringBuilder();
		for (String prefixFunctionalTerm : seperateWordSubstituents) {
			casName.append(prefixFunctionalTerm);
			casName.append(" ");
		}
		for (String substituent : substituents) {
			casName.append(substituent);
		}
		casName.append(parent);
		for (String functionalTerm : functionalTerms) {
			casName.append(" ");
			casName.append(functionalTerm);
		}
		return casName.toString();
	}

	private static Character missingCloseBracketCharIfApplicable(String component) {
		char[] characters = component.toCharArray();
		int bracketLevel =0;
		Character missingCloseBracket =null;
		for (int i = 0; i < characters.length; i++) {
			char character = characters[i];
			if (character == '(' || character == '[' || character == '{') {
				bracketLevel++;
				if (bracketLevel ==1){
					missingCloseBracket = character;
				}
			}
			if (character == ')' || character == ']' || character == '}') {
				bracketLevel--;
				if (bracketLevel<0){
					return null;
				}
			}
		}
		if (bracketLevel == 1){
			if (missingCloseBracket == '('){
				return ')';
			}
			if (missingCloseBracket == '['){
				return ']';
			}
			if (missingCloseBracket == '{'){
				return '}';
			}
		}
		return null;
	}

	/**
	 * Modifies the name of the parent acid from ic to ate (or ous to ite)
	 * hence allowing the formation of the uninverted ester
	 * @param parent
	 * @return
	 * @throws ParsingException
	 */
	private static String uninvertEster(String parent) throws ParsingException {
		int len = parent.length();
		if (len < 9) {
			throw new ParsingException("Failed to uninvert CAS ester");
		}
		char lastChar = parent.charAt(len - 1);
		if (lastChar == ')' || lastChar == ']' || lastChar == '}') {
			if (parent.substring(parent.length() - 8).equalsIgnoreCase("ic acid)")) {
				parent = parent.substring(0, parent.length() - 8) + "ate)";
			} else if (parent.substring(parent.length() - 9).equalsIgnoreCase("ous acid)")) {
				parent = parent.substring(0, parent.length() - 9) + "ite)";
			} else {
				throw new ParsingException("Failed to uninvert CAS ester");
			}
		} else {
			if (parent.substring(parent.length() - 7).equalsIgnoreCase("ic acid")) {
				parent = parent.substring(0, parent.length() - 7) + "ate";
			} else if (parent.substring(parent.length() - 8).equalsIgnoreCase("ous acid")) {
				parent = parent.substring(0, parent.length() - 8) + "ite";
			} else {
				throw new ParsingException("Failed to uninvert CAS ester");
			}
		}
		return parent;
	}
}
