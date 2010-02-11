package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;

import uk.ac.cam.ch.wwmm.opsin.ParseWord.WordType;

import nu.xom.Attribute;
import nu.xom.Element;
import nu.xom.Elements;

/**The rules by which words are grouped together (e.g. in functional class nomenclature)
 *
 * @author dl387
 *
 */
class WordRules {

	/**
	 * The currently supported wordRules
	 * @author dl387
	 *
	 */
	enum WordRule{
		acid,
		amide,
		anhydride,
		carbonylDerivative,
		divalentFunctionalGroup,
		ester,
		full,
		functionalClassEster,
		functionGroupAsGroup,
		glycol,
		monovalentFunctionalGroup,
		multiEster,
		oxide,
		polymer,
		simple
	}

	/**
	 * Describes a word that a wordRule is looking for
	 * @author dl387
	 *
	 */
	class WordDescription {
		/**Whether the word is a full word, substituent word or functionalTerm word*/
		private final WordType type;

		/**A case insensitive pattern which attempts to match the end of the String value of the word*/
		private Pattern endsWithPattern = null;

		/**The case insensitive String value of the word */
		private String value = null;

		/** Only applicable for functionalTerms. The string value of the functionalTerm's type attribute*/
		private String functionalGroupType = null;

		/** The value of the value attribute of the last group element in the word e.g. maybe a SMILES string*/
		private String endsWithGroupValueAtr = null;

		WordType getType() {
			return type;
		}

		Pattern getEndsWithPattern() {
			return endsWithPattern;
		}

		void setEndsWithPattern(Pattern endsWithPattern) {
			this.endsWithPattern = endsWithPattern;
		}

		String getValue() {
			return value;
		}

		void setValue(String value) {
			this.value = value.toLowerCase();
		}

		String getFunctionalGroupType() {
			return functionalGroupType;
		}

		void setFunctionalGroupType(String functionalGroupType) {
			this.functionalGroupType = functionalGroupType;
		}

		String getEndsWithElementValueAtr() {
			return endsWithGroupValueAtr;
		}

		void setEndsWithElementValueAtr(String endsWithElementValueAtr) {
			this.endsWithGroupValueAtr = endsWithElementValueAtr;
		}

		/**
		 * Makes a description of a word to looks for
		 * @param wordType
		 */
		WordDescription(WordType wordType){
			type =wordType;
		}
	}

	/**
	 * A representation of a wordRule element from wordRules.xml
	 * @author dl387
	 *
	 */
	class WordRuleDescription {
		private final List<WordDescription> wordDescriptions = new ArrayList<WordDescription>();
		private final String ruleName;
		private final String ruleType;

		String getRuleName() {
			return ruleName;
		}

		String getRuleType() {
			return ruleType;
		}

		/**
		 * Creates a wordRule from a wordRule element found in wordRules.xml
		 * @param wordRuleEl
		 */
		WordRuleDescription(Element wordRuleEl) {
			ruleName =wordRuleEl.getAttributeValue("name");
			ruleType =wordRuleEl.getAttributeValue("type");
			Elements words = wordRuleEl.getChildElements();
			for (int i = 0; i < words.size(); i++) {
				Element word = words.get(i);
				WordDescription wd = new WordDescription(WordType.valueOf(word.getAttributeValue("type")));
				if (word.getAttribute("value")!=null){
					wd.setValue(word.getAttributeValue("value"));
				}
				if (word.getAttribute("functionalGroupType")!=null){
					wd.setFunctionalGroupType(word.getAttributeValue("functionalGroupType"));
				}
				if (word.getAttribute("endsWithRegex")!=null){
					wd.setEndsWithPattern(Pattern.compile(word.getAttributeValue("endsWithRegex") +"$", Pattern.CASE_INSENSITIVE));
				}
				if (word.getAttribute("endsWithGroupValueAtr")!=null){
					wd.setEndsWithElementValueAtr(word.getAttributeValue("endsWithGroupValueAtr"));
				}
				wordDescriptions.add(wd);
			}
		}
	}


	/**Initialises the WordRules.
	 * @param resourceGetter
	 *
	 * @throws Exception If the data file can't be read properly.
	 */
	WordRules(ResourceGetter resourceGetter) throws Exception {
		Element wordRules =resourceGetter.getXMLDocument("wordRules.xml").getRootElement();
		Elements rules = wordRules.getChildElements("wordRule");
		for (int i = 0; i < rules.size(); i++) {
			wordRuleList.add (new WordRuleDescription(rules.get(i)));
		}
	}

	/**The wordRules themselves.*/
	private final List<WordRuleDescription> wordRuleList = new ArrayList<WordRuleDescription>();

	/**Takes a molecule element and places the word elements into wordRule elements
	 *
	 * @param moleculeEl A molecule element with word children
	 * @throws ParsingException
	 */
	void groupWordsIntoWordRules(Element moleculeEl) throws ParsingException {
		List<Element> wordEls = XOMTools.getChildElementsWithTagNames(moleculeEl, new String[]{"word"});
		//note that multiple words in wordEls may be later replaced by a wordRule element
		for (int i = 0; i <wordEls.size(); i++) {
			if (matchWordRule(wordEls, i)){
				i=-1;//if function did something
			}
		}
		Elements wordRuleEls = moleculeEl.getChildElements();
		for (int i = 0; i < wordRuleEls.size(); i++) {
			Element wordRuleEl = wordRuleEls.get(i);
			if (!wordRuleEl.getLocalName().equals("wordRule")){
				throw new ParsingException("Unable to assign wordRule to: " + wordRuleEl.getValue());
			}
		}
	}

	private boolean matchWordRule(List<Element> wordEls, int indexOfFirstWord) throws ParsingException {
		wordRuleLoop: for (WordRuleDescription wordRule : wordRuleList) {
			int i =indexOfFirstWord;
			int wordsInWordRule = wordRule.wordDescriptions.size();
			if (i + wordsInWordRule -1 < wordEls.size()){//need sufficient words to match the word rule
				for (int j = 0; j < wordsInWordRule; j++) {
					Element wordEl = wordEls.get(i+j);
					WordDescription wd = wordRule.wordDescriptions.get(j);
					if (!wd.type.toString().equals(wordEl.getAttributeValue("type"))){
						continue wordRuleLoop;//type mismatch;
					}
					if (wd.getValue() !=null && !wordEl.getAttributeValue("value").toLowerCase().equals(wd.getValue())){//word string contents mismatch
						continue wordRuleLoop;
					}
					if (wd.functionalGroupType !=null){
						if (WordType.functionalTerm.toString().equals(wordEl.getAttributeValue("type"))){
							Elements children = wordEl.getChildElements();
							Element lastChild = children.get(children.size()-1);
							while (lastChild.getChildElements().size()!=0){
								children = lastChild.getChildElements();
								lastChild = children.get(children.size()-1);
							}
							if (lastChild.getLocalName().equals("closebracket")){
								lastChild = (Element) XOMTools.getPreviousSibling(lastChild);
							}
							if (lastChild==null){
								throw new ParsingException("OPSIN Bug: Cannot find the functional element in a functionalTerm");
							}
							if (!wd.getFunctionalGroupType().equals(lastChild.getAttributeValue("type"))){
								continue wordRuleLoop;
							}
						}
					}
					if (wd.endsWithPattern !=null){
						if (!wd.endsWithPattern.matcher(wordEl.getAttributeValue("value")).find()){
							continue wordRuleLoop;
						}
					}
					if (wd.endsWithGroupValueAtr !=null){
						Elements children = wordEl.getChildElements();
						Element lastChild = children.get(children.size()-1);
						while (lastChild.getChildElements().size()!=0){
							children = lastChild.getChildElements();
							lastChild = children.get(children.size()-1);
						}
						Element groupToExamine;
						if (lastChild.getLocalName().equals("group")){
							groupToExamine = lastChild;
						}
						else{
							Elements groups = ((Element)lastChild.getParent()).getChildElements("group");
							if (groups.size()>0){
								groupToExamine = groups.get(groups.size()-1);
							}
							else{
								continue wordRuleLoop;
							}
						}
						if (groupToExamine.getAttribute("value")==null || !groupToExamine.getAttributeValue("value").equals(wd.endsWithGroupValueAtr)){
							continue wordRuleLoop;
						}
					}
				}
				//Word Rule matches!
				Element wordRuleEl = new Element("wordRule");
				wordRuleEl.addAttribute(new Attribute("type", wordRule.getRuleType()));
				wordRuleEl.addAttribute(new Attribute("wordRule", wordRule.getRuleName()));

				/*
				 * Some wordRules can not be entirely processed at the structure building stage
				 */
				if (wordRule.getRuleName().equals(WordRule.functionGroupAsGroup.toString())){//convert the functional term into a full term
					if (wordsInWordRule!=1){
						throw new ParsingException("OPSIN bug: Problem with functionGroupAsGroup wordRule");
					}
					convertFunctionalGroupIntoGroup(wordEls.get(i));
					wordRuleEl.getAttribute("wordRule").setValue(WordRule.simple.toString());
				}
				else if (wordRule.getRuleName().equals(WordRule.carbonylDerivative.toString())){//e.g. 4,4-diphenylsemicarbazone. This is better expressed as a full word as the substituent actually locants onto the functional term
					if (wordsInWordRule==3){//substituent present
						joinWords(wordEls, i+1, wordEls.get(i+1), wordEls.get(i+2));
						wordsInWordRule--;
						List<Element> functionalTerm = XOMTools.getDescendantElementsWithTagName(wordEls.get(i+1), "functionalTerm");//rename functionalTerm element to root
						if (functionalTerm.size()!=1){
							throw new ParsingException("OPSIN bug: Problem with carbonylDerivative wordRule");
						}
						functionalTerm.get(0).setLocalName("root");
						wordEls.get(i+1).getAttribute("type").setValue(WordType.full.toString());
					}
				}



				List<String> wordValues = new ArrayList<String>();
				Element parentEl = (Element) wordEls.get(i).getParent();
				int indexToInsertAt = parentEl.indexOf(wordEls.get(i));
				for (int j = 0; j < wordsInWordRule; j++) {
					Element wordEl = wordEls.remove(i);
					wordEl.detach();
					wordRuleEl.appendChild(wordEl);
					wordValues.add(wordEl.getAttributeValue("value"));
				}
				wordRuleEl.addAttribute(new Attribute("value", StringTools.stringListToString(wordValues, " ")));//The bare string of all the words under this wordRule
				parentEl.insertChild(wordRuleEl, indexToInsertAt);
				wordEls.add(i, wordRuleEl);
				return true;
			}
		}
		Element firstWord = wordEls.get(indexOfFirstWord);
		if (firstWord.getLocalName().equals("word") && WordType.full.toString().equals(firstWord.getAttributeValue("type"))){//No wordRule -->wordRule="simple"
			applySimpleWordRule(wordEls, indexOfFirstWord, firstWord);
			return false;
		}
		else if (WordType.substituent.toString().equals(firstWord.getAttributeValue("type"))){
			/*
			 * substituents may join together or to a full e.g. 2-ethyl toluene -->2-ethyltoluene
			 * 1-chloro 2-bromo ethane --> 1-chloro-2-bromo ethane then subsequently 1-chloro-2-bromo-ethane
			 */
			if (indexOfFirstWord +1 < wordEls.size()){
				Element wordToPotentiallyCombineWith = wordEls.get(indexOfFirstWord +1);
				if (WordType.full.toString().equals(wordToPotentiallyCombineWith.getAttributeValue("type")) ||
				WordType.substituent.toString().equals(wordToPotentiallyCombineWith.getAttributeValue("type"))){
					joinWords(wordEls, indexOfFirstWord, firstWord, wordToPotentiallyCombineWith);
					return true;
				}
			}
		}
		return false;
	}

	private void applySimpleWordRule(List<Element> wordEls, int indexOfFirstWord, Element firstWord) {
		Element parentEl = (Element) firstWord.getParent();
		int indexToInsertAt = firstWord.getParent().indexOf(firstWord);
		Element wordRuleEl = new Element("wordRule");
		wordRuleEl.addAttribute(new Attribute("wordRule", WordRule.simple.toString()));//No wordRule
		wordRuleEl.addAttribute(new Attribute("type", WordType.full.toString()));
		wordRuleEl.addAttribute(new Attribute("value", firstWord.getAttributeValue("value")));
		wordEls.remove(indexOfFirstWord);
		firstWord.detach();
		wordRuleEl.appendChild(firstWord);
		wordEls.add(indexOfFirstWord, wordRuleEl);
		parentEl.insertChild(wordRuleEl, indexToInsertAt);
	}

	/**
	 * Takes the list of wordEls, the indice of a word element, that element and the word element following it
	 * Merges the latter word element into the former element
	 * @param wordEls
	 * @param indexOfFirstWord
	 * @param firstWord
	 * @param wordToPotentiallyCombineWith
	 * @throws ParsingException
	 */
	private void joinWords(List<Element> wordEls, int indexOfFirstWord, Element firstWord, Element wordToPotentiallyCombineWith) throws ParsingException {
		wordEls.remove(indexOfFirstWord +1);
		wordToPotentiallyCombineWith.detach();
		Element assumedHyphen = new Element("hyphen");
		assumedHyphen.appendChild("-");
		Elements substituentEls = firstWord.getChildElements("substituent");
		if (substituentEls.size()==0){
			throw new ParsingException("OPSIN Bug: Substituent element not found where substituent element expected");
		}
		Element finalSubstituent = substituentEls.get(substituentEls.size()-1);
		finalSubstituent.appendChild(assumedHyphen);
		Elements elementsToMergeIntoSubstituent = wordToPotentiallyCombineWith.getChildElements();
		for (int j =  elementsToMergeIntoSubstituent.size() -1 ; j >=0; j--) {
			Element el = elementsToMergeIntoSubstituent.get(j);
			el.detach();
			XOMTools.insertAfter(finalSubstituent, el);
		}
		if (WordType.full.toString().equals(wordToPotentiallyCombineWith.getAttributeValue("type"))){
			firstWord.getAttribute("type").setValue(WordType.full.toString());
		}
		firstWord.getAttribute("value").setValue(firstWord.getAttributeValue("value") + wordToPotentiallyCombineWith.getAttributeValue("value"));
	}

	private void convertFunctionalGroupIntoGroup(Element word) throws ParsingException {
		word.getAttribute("type").setValue(WordType.full.toString());
		List<Element> functionalTerms = XOMTools.getDescendantElementsWithTagName(word, "functionalTerm");
		if (functionalTerms.size()!=1){
			throw new ParsingException("OPSIN Bug: Exactly 1 functionalTerm expected in functionalGroupAsGroup wordRule");
		}
		functionalTerms.get(0).setLocalName("root");
		Elements functionalGroups = functionalTerms.get(0).getChildElements("functionalGroup");
		if (functionalGroups.size()!=1){
			throw new ParsingException("OPSIN Bug: Exactly 1 functionalGroup expected in functionalGroupAsGroup wordRule");
		}
		functionalGroups.get(0).setLocalName("group");
		functionalGroups.get(0).getAttribute("type").setValue("simpleGroup");
		functionalGroups.get(0).addAttribute(new Attribute("subType", "functionalGroupAsGroup"));
	}

}
