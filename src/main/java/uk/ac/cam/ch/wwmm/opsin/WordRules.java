package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import dk.brics.automaton.RunAutomaton;

import nu.xom.Element;
import nu.xom.Elements;

/**The rules by which names are divided into words.
 *
 * @author ptc24
 *
 */
class WordRules {

	private final static char NON_CHEM_CHAR = '\u00A3';

	/**A rule by which names are divided into words.*/
	class WordRule {

		/**The regular expression to recognise the applicability of a wordRule.*/
		Pattern regex;
		/**The string used by the automaton to build a regex */
		String regexForAutomaton;
		/**The type of each word in the rule.*/
		List<String> wordTypes;
		/**The name of the wordRule.*/
		String name;
		/**Does this wordRule employ techniques not implemented by automaton e.g. lookbehind.*/
		boolean hasFeaturesUnsupportedByAutomaton =false;

		/**Makes a wordRule based on an XML element.
		 *
		 * @param elem The XML element containing the information for the wordRule.
		 * @param id Used to identify this wordRule
		 * @throws ParsingException
		 */
		WordRule(Element elem, int id) throws ParsingException {
			name = elem.getAttributeValue("name");
			String rawRegex =elem.getAttributeValue("regex");
			Matcher m = matchStringReplacement.matcher(rawRegex);
			String newValue = "";
			int position = 0;
			while (m.find()){
				newValue += rawRegex.substring(position, m.start());
				if (stringReplacements.get(m.group())==null){
					throw new ParsingException("Regex entry for: " + m.group() + " missing! Check wordRules.xml");
				}
				newValue += stringReplacements.get(m.group());
				position = m.end();
			}
			newValue += rawRegex.substring(position);
			regex = Pattern.compile(newValue,Pattern.CASE_INSENSITIVE );

			/* The \u00A3 is a delimiter which is expected never to occur in chemical names
			 * The idea of this is that the states after the acceptState of the automaton can be examined and from the characters that
			 * were accepted to reach these states the id can be retrieved hence allowing the identification of which regex matched
			 * (the automaton has all word rules fed to it!)
			 */
			regexForAutomaton =newValue + "(" + NON_CHEM_CHAR + id +")?";
			m = matchLookBehind.matcher(regexForAutomaton);
			if (m.find()){
				hasFeaturesUnsupportedByAutomaton = true;
				regexForAutomaton =m.replaceAll("");
			}
			regexForAutomaton =matchNonCapturing.matcher(regexForAutomaton).replaceAll("");
			wordTypes = new ArrayList<String>();
			Elements wordElems = elem.getChildElements("word");
			for(int i=0;i<wordElems.size();i++) {
				wordTypes.add(wordElems.get(i).getAttributeValue("type"));
			}
		}
	}

	/**The wordRules themselves.*/
	private List<WordRule> wordRuleList;

	private HashMap<String, String> stringReplacements;
	private Pattern matchStringReplacement = Pattern.compile("%.*?%");
	private Pattern matchNonCapturing = Pattern.compile("\\?:");
	private Pattern matchLookBehind = Pattern.compile("\\(\\?<!(([^(]+?\\))|(\\([^(]+?\\)\\)))");
	private RunAutomaton allWordRulesRunAutomaton;
	private Pattern matchPoly = Pattern.compile("(?:poly|oligo)[\\[\\(\\{](.+)[\\]\\)\\}]", Pattern.CASE_INSENSITIVE );//poly or oligo followed by a bracket name

	/**Initialises the WordRules.
	 * @param resourceGetter
	 *
	 * @throws Exception If the data file can't be read properly.
	 */
	WordRules(ResourceGetter resourceGetter) throws Exception {
		Element wordRules =resourceGetter.getXMLDocument("wordRules.xml").getRootElement();
		Elements stringReplacementEls = wordRules.getChildElements("stringReplacement");
		Elements rules = wordRules.getChildElements("wordRule");
		stringReplacements = new HashMap<String, String>();
		for (int i = 0; i < stringReplacementEls.size(); i++) {
			String stringReplacement =stringReplacementEls.get(i).getAttributeValue("value");
			String name =stringReplacementEls.get(i).getAttributeValue("name");
			Matcher m = matchStringReplacement.matcher(stringReplacement);
			String newValue = "";
			int position = 0;
			while (m.find()){
				newValue += stringReplacement.substring(position, m.start());
				if (stringReplacements.get(m.group())==null){
					throw new ParsingException("Regex entry for: " + m.group() + " missing! Check wordRules.xml");
				}
				newValue += stringReplacements.get(m.group());
				position = m.end();
			}
			newValue += stringReplacement.substring(position);
			stringReplacements.put(name, newValue);
		}
		wordRuleList = new ArrayList<WordRule>();
		for(int i=0;i<rules.size();i++) {
			WordRule wr = new WordRule(rules.get(i), i);
			wordRuleList.add(wr);
		}
		String allWordRulesRegex ="";//given to automaton
		for (WordRule wr : wordRuleList) {
			allWordRulesRegex +="(";
			allWordRulesRegex += wr.regexForAutomaton;
			allWordRulesRegex +=")|";
		}
		allWordRulesRegex = allWordRulesRegex.substring(0, allWordRulesRegex.length()-1);
		allWordRulesRunAutomaton =AutomatonInitialiser.getAutomaton("wordRules", allWordRulesRegex, resourceGetter);
	}

	/**Takes a chemical name, breaks it up into words, and works out which
	 * wordRule should apply to it.
	 *
	 * @param p The Parse object containing the name, into which the results will be put.
	 * @throws ParsingException
	 */
	void parse(Parse p) throws ParsingException {
		String chemicalName =p.name;

		//TODO do this properly
		Matcher polyMatcher = matchPoly.matcher(chemicalName);//temporary kludge until wordRules system is rewritten
		if (polyMatcher.matches()){//Name is a polymer/oligomer name
			chemicalName = polyMatcher.group(1);//strip off starting poly( and trailing closing bracket
			String [] wordArray = chemicalName.split("\\s+");
			if (wordArray.length==1){
				p.wordRule="polymer";
				ParseWord pw = new ParseWord();
				pw.word = wordArray[0];
				pw.wordType = "substituent";
				p.words.add(pw);
				return;
			}
			else{
				throw new ParsingException("Unsupported Polymer name");
			}
		}
		String chemicalNameLowerCase =chemicalName.toLowerCase();
		char[] chemicalNameArray=chemicalName.toCharArray();
		List<Integer> parsingStartingPoints =new ArrayList<Integer>();
		/* Word rules may be applied from the start of the word, and from after any spaces.
		 * Only at maximum one word rule will be successfully applied
		 */
		parsingStartingPoints.add(0);
		for (int i = 0; i < chemicalNameArray.length -1; i++) {
			if (chemicalNameArray[i]==' '){
				parsingStartingPoints.add(i+1);
			}
		}

		for (Integer startingPoint : parsingStartingPoints) {
			int state =allWordRulesRunAutomaton.getInitialState();
			ArrayList<Integer> acceptStates =new ArrayList<Integer>();//typically will only be one or no accept states
			for (int i = startingPoint; i < chemicalNameLowerCase.length(); i++) {
				state = allWordRulesRunAutomaton.step(state, chemicalNameLowerCase.charAt(i));
				if (state == -1){
					break;
				}
				if (allWordRulesRunAutomaton.isAccept(state)){
					acceptStates.add(state);
				}
			}
			if (acceptStates.size()==0){
				continue;
			}
			Collections.reverse(acceptStates);//those which describe more of the chemical name will be tried first

			for (Integer acceptState : acceptStates) {
				state = allWordRulesRunAutomaton.step(acceptState, NON_CHEM_CHAR);

				/*
				 * Determine which regexes succeeded in matching.
				 * The regexes with lower ids are more specific and hence the one with the lowest id will be used in preference
				 */
				ArrayList<Integer> wordIDs = determineAllWordRulesThatMatched(state, "", new ArrayList<Integer>());
				Collections.sort(wordIDs);
				for (Integer idOfWordRule : wordIDs) {//typically there will only be one wordID
					int wordRule =Integer.valueOf(idOfWordRule);

					WordRule r =wordRuleList.get(wordRule);
					Matcher m = r.regex.matcher(chemicalNameLowerCase);
					if(m.find(startingPoint)) {
						int end =m.end();
						if (end <chemicalNameLowerCase.length()){
							if (chemicalNameLowerCase.charAt(end)!=' '){//regex has seen a functional class name in the middle of a word e.g. thiolate
								continue;
							}
						}
						//System.out.println(r.name);
						p.wordRule = r.name;
						for (int i = 1; i <= m.groupCount(); i++) {
							ParseWord pw = new ParseWord();
							pw.word = chemicalName.substring(m.start(i), m.end(i));
							if(i-1 < r.wordTypes.size()){
								pw.wordType = r.wordTypes.get(i-1);
							}
							else{
								pw.wordType = "full";//extra words are likely to be counter ions or it's a mixture
							}
							p.words.add(pw);
						}
						int start =m.start();
						if (start!=0){
							String [] wordArray = chemicalName.substring(0,start).trim().split("\\s+");
							for(int i=0;i<wordArray.length;i++) {
								ParseWord pw = new ParseWord();
								pw.word = wordArray[i];
								pw.wordType = "full";//extra words are likely to be counter ions or it's a mixture
								p.words.add(pw);
							}
						}

						if (end!=chemicalName.length()){
							String [] wordArray = chemicalName.substring(end, chemicalName.length()).trim().split("\\s+");
							for(int i=0;i<wordArray.length;i++) {
								ParseWord pw = new ParseWord();
								pw.word = wordArray[i];
								pw.wordType = "full";//extra words are likely to be counter ions or it's a mixture
								p.words.add(pw);
							}
						}
						return;
					}
					else if (!r.hasFeaturesUnsupportedByAutomaton){
						throw new ParsingException("Fault in parser, automaton found a match, but JAVA regex couldn't");
					}
				}
			}
		}
		//none of the recognised word rules
		String [] wordArray = chemicalName.split("\\s+");
		if (wordArray.length==1){
			p.wordRule="simple";
			ParseWord pw = new ParseWord();
			pw.word = wordArray[0];
			pw.wordType = "full";
			p.words.add(pw);
		}
		else{
			p.wordRule="binaryOrOther";
			for(int i=0;i<wordArray.length;i++) {
				ParseWord pw = new ParseWord();
				pw.word = wordArray[i];
				pw.wordType = "full";
				p.words.add(pw);
			}
		}
	}

	/**
	 * Called recursively to populate wordIDs array
	 * @param state
	 * @param wordIDSoFar
	 * @param wordIDs
	 * @return
	 */
	private ArrayList<Integer> determineAllWordRulesThatMatched(int state, String wordIDSoFar, ArrayList<Integer> wordIDs) {
		char[] stateSymbols = allWordRulesRunAutomaton.getCharIntervals();
        for (char stateSymbol : stateSymbols) {
            int potentialNextState = allWordRulesRunAutomaton.step(state, stateSymbol);
            if (potentialNextState != -1) {
                int currentState = potentialNextState;
                wordIDSoFar += stateSymbol;
                if (allWordRulesRunAutomaton.isAccept(currentState)) {
                    wordIDs.add(Integer.valueOf(wordIDSoFar));
                } else {
                    determineAllWordRulesThatMatched(currentState, wordIDSoFar, wordIDs);
                }
            }
        }
		return wordIDs;
	}
}
