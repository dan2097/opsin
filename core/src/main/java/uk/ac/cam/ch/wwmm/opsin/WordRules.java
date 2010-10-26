package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;

import uk.ac.cam.ch.wwmm.opsin.ParseWord.WordType;

import nu.xom.Attribute;
import nu.xom.Element;
import nu.xom.Elements;

import static uk.ac.cam.ch.wwmm.opsin.XmlDeclarations.*;

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
		acetal,
		additionCompound,
		acidHalideOrPseudoHalide,
		amide,
		anhydride,
		biochemicalEster,
		carbonylDerivative,
		divalentFunctionalGroup,
		ester,
		functionalClassEster,
		functionGroupAsGroup,
		glycol,
		glycolEther,
		hydrazide,
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
		
		/** The value of the type attribute of the last group element in the word e.g. maybe aminoAcid*/
		private String endsWithGroupType = null;
		
		/** The value of the subType attribute of the last group element in the word e.g. maybe elementaryAtom*/
		private String endsWithGroupSubType = null;

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

		String getEndsWithGroupValueAtr() {
			return endsWithGroupValueAtr;
		}

		void setEndsWithGroupValueAtr(String endsWithElementValueAtr) {
			this.endsWithGroupValueAtr = endsWithElementValueAtr;
		}

		String getEndsWithGroupType() {
			return endsWithGroupType;
		}

		void setEndsWithGroupType(String endsWithElementType) {
			this.endsWithGroupType = endsWithElementType;
		}
		
		String getEndsWithGroupSubType() {
			return endsWithGroupSubType;
		}

		void setEndsWithGroupSubType(String endsWithElementSubType) {
			this.endsWithGroupSubType = endsWithElementSubType;
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
					wd.setEndsWithGroupValueAtr(word.getAttributeValue("endsWithGroupValueAtr"));
				}
				if (word.getAttribute("endsWithGroupType")!=null){
					wd.setEndsWithGroupType(word.getAttributeValue("endsWithGroupType"));
				}
				if (word.getAttribute("endsWithGroupSubType")!=null){
					wd.setEndsWithGroupSubType(word.getAttributeValue("endsWithGroupSubType"));
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
	 * @param n2sConfig 
	 *
	 * @param moleculeEl A molecule element with word children
	 * @param allowSpaceRemoval 
	 * @throws ParsingException
	 */
	void groupWordsIntoWordRules(NameToStructureConfig n2sConfig, Element moleculeEl, boolean allowSpaceRemoval) throws ParsingException {
		List<Element> wordEls = XOMTools.getChildElementsWithTagName(moleculeEl, WORD_EL);
		//note that multiple words in wordEls may be later replaced by a wordRule element
		for (int i = 0; i <wordEls.size(); i++) {
			if (matchWordRule(n2sConfig, wordEls, i, allowSpaceRemoval)){
				i=-1;//if function did something
			}
		}
		Elements wordRuleEls = moleculeEl.getChildElements();
		for (int i = 0; i < wordRuleEls.size(); i++) {
			Element wordRuleEl = wordRuleEls.get(i);
			if (!wordRuleEl.getLocalName().equals(WORDRULE_EL)){
				throw new ParsingException("Unable to assign wordRule to: " + wordRuleEl.getValue());
			}
		}
	}

	private boolean matchWordRule(NameToStructureConfig n2sConfig, List<Element> wordEls, int indexOfFirstWord, boolean allowSpaceRemoval) throws ParsingException {
		wordRuleLoop: for (WordRuleDescription wordRule : wordRuleList) {
			int i =indexOfFirstWord;
			int wordsInWordRule = wordRule.wordDescriptions.size();
			if (i + wordsInWordRule -1 < wordEls.size()){//need sufficient words to match the word rule
				for (int j = 0; j < wordsInWordRule; j++) {
					Element wordEl = wordEls.get(i+j);
					WordDescription wd = wordRule.wordDescriptions.get(j);
					if (!wd.type.toString().equals(wordEl.getAttributeValue(TYPE_ATR))){
						continue wordRuleLoop;//type mismatch;
					}
					if (wd.getValue() !=null && !wordEl.getAttributeValue(VALUE_ATR).toLowerCase().equals(wd.getValue())){//word string contents mismatch
						continue wordRuleLoop;
					}
					if (wd.functionalGroupType !=null){
						if (WordType.functionalTerm.toString().equals(wordEl.getAttributeValue(TYPE_ATR))){
							Elements children = wordEl.getChildElements();
							Element lastChild = children.get(children.size()-1);
							while (lastChild.getChildElements().size()!=0){
								children = lastChild.getChildElements();
								lastChild = children.get(children.size()-1);
							}
							if (lastChild.getLocalName().equals(CLOSEBRACKET_EL)){
								lastChild = (Element) XOMTools.getPreviousSibling(lastChild);
							}
							if (lastChild==null){
								throw new ParsingException("OPSIN Bug: Cannot find the functional element in a functionalTerm");
							}
							if (!wd.getFunctionalGroupType().equals(lastChild.getAttributeValue(TYPE_ATR))){
								continue wordRuleLoop;
							}
						}
					}
					if (wd.endsWithPattern !=null){
						if (!wd.endsWithPattern.matcher(wordEl.getAttributeValue(VALUE_ATR)).find()){
							continue wordRuleLoop;
						}
					}
					if (wd.endsWithGroupValueAtr !=null){
						Element lastGroupInWordRule = getLastGroupInWordRule(wordEl);
						if (lastGroupInWordRule==null || !wd.endsWithGroupValueAtr.equals(lastGroupInWordRule.getAttributeValue(VALUE_ATR))){
							continue wordRuleLoop;
						}
					}
					if (wd.endsWithGroupType !=null){
						Element lastGroupInWordRule = getLastGroupInWordRule(wordEl);
						if (lastGroupInWordRule==null || !wd.endsWithGroupType.equals(lastGroupInWordRule.getAttributeValue(TYPE_ATR))){
							continue wordRuleLoop;
						}
					}
					if (wd.endsWithGroupSubType !=null){
						Element lastGroupInWordRule = getLastGroupInWordRule(wordEl);
						if (lastGroupInWordRule==null || !wd.endsWithGroupSubType.equals(lastGroupInWordRule.getAttributeValue(SUBTYPE_ATR))){
							continue wordRuleLoop;
						}
					}
				}
				//Word Rule matches!
				Element wordRuleEl = new Element(WORDRULE_EL);
				wordRuleEl.addAttribute(new Attribute(TYPE_ATR, wordRule.getRuleType()));
				wordRuleEl.addAttribute(new Attribute(WORDRULE_EL, wordRule.getRuleName()));

				/*
				 * Some wordRules can not be entirely processed at the structure building stage
				 */
				String wordRuleName = wordRule.getRuleName();
				if (wordRuleName.equals(WordRule.functionGroupAsGroup.toString())){//convert the functional term into a full term
					if (wordsInWordRule!=1){
						throw new ParsingException("OPSIN bug: Problem with functionGroupAsGroup wordRule");
					}
					convertFunctionalGroupIntoGroup(wordEls.get(i));
					wordRuleEl.getAttribute(WORDRULE_ATR).setValue(WordRule.simple.toString());
				}
				else if (wordRuleName.equals(WordRule.carbonylDerivative.toString())){//e.g. acetone 4,4-diphenylsemicarbazone. This is better expressed as a full word as the substituent actually locants onto the functional term
					if (wordsInWordRule==3){//substituent present
						joinWords(wordEls, i+1, wordEls.get(i+1), wordEls.get(i+2));
						wordsInWordRule--;
						List<Element> functionalTerm = XOMTools.getDescendantElementsWithTagName(wordEls.get(i+1), FUNCTIONALTERM_EL);//rename functionalTerm element to root
						if (functionalTerm.size()!=1){
							throw new ParsingException("OPSIN bug: Problem with carbonylDerivative wordRule");
						}
						functionalTerm.get(0).setLocalName(ROOT_EL);
						List<Element> functionalGroups = XOMTools.getDescendantElementsWithTagName(functionalTerm.get(0), FUNCTIONALGROUP_EL);//rename functionalGroup element to group
						if (functionalGroups.size()!=1){
							throw new ParsingException("OPSIN bug: Problem with carbonylDerivative wordRule");
						}
						functionalGroups.get(0).setLocalName(GROUP_EL);
						wordEls.get(i+1).getAttribute(TYPE_ATR).setValue(WordType.full.toString());
					}
				}
				else if (wordRuleName.equals(WordRule.additionCompound.toString()) || wordRuleName.equals(WordRule.oxide.toString())){//is the halide/pseudohalide/oxide actually a counterion rather than covalently bonded
					Element possibleElementaryAtom = wordEls.get(i);
					List<Element> elementaryAtoms = XOMTools.getDescendantElementsWithTagNameAndAttribute(possibleElementaryAtom, GROUP_EL, SUBTYPE_ATR, ELEMENTARYATOM_SUBTYPE_VAL);
					if (elementaryAtoms.size()==1){
						for (int j = 1; j < wordsInWordRule; j++) {
							if (bondWillBeIonic(elementaryAtoms.get(0), wordEls.get(i+j))){//use separate word rules for ionic components
								continue wordRuleLoop;
							}
						}
					}
				}



				List<String> wordValues = new ArrayList<String>();
				Element parentEl = (Element) wordEls.get(i).getParent();
				int indexToInsertAt = parentEl.indexOf(wordEls.get(i));
				for (int j = 0; j < wordsInWordRule; j++) {
					Element wordEl = wordEls.remove(i);
					wordEl.detach();
					wordRuleEl.appendChild(wordEl);
					wordValues.add(wordEl.getAttributeValue(VALUE_ATR));
				}
				wordRuleEl.addAttribute(new Attribute(VALUE_ATR, StringTools.stringListToString(wordValues, " ")));//The bare string of all the words under this wordRule
				parentEl.insertChild(wordRuleEl, indexToInsertAt);
				wordEls.add(i, wordRuleEl);
				return true;
			}
		}
		Element firstWord = wordEls.get(indexOfFirstWord);
		if (firstWord.getLocalName().equals(WORD_EL) && WordType.full.toString().equals(firstWord.getAttributeValue(TYPE_ATR))){//No wordRule -->wordRule="simple"
			applySimpleWordRule(wordEls, indexOfFirstWord, firstWord);
			return false;
		}
		else if (allowSpaceRemoval && WordType.substituent.toString().equals(firstWord.getAttributeValue(TYPE_ATR))){
			/*
			 * substituents may join together or to a full e.g. 2-ethyl toluene -->2-ethyltoluene
			 * 1-chloro 2-bromo ethane --> 1-chloro-2-bromo ethane then subsequently 1-chloro-2-bromo-ethane
			 */
			if (indexOfFirstWord +1 < wordEls.size()){
				Element wordToPotentiallyCombineWith = wordEls.get(indexOfFirstWord +1);
				if (WordType.full.toString().equals(wordToPotentiallyCombineWith.getAttributeValue(TYPE_ATR)) ||
				WordType.substituent.toString().equals(wordToPotentiallyCombineWith.getAttributeValue(TYPE_ATR))){
					joinWords(wordEls, indexOfFirstWord, firstWord, wordToPotentiallyCombineWith);
					return true;
				}
			}
		}
		if (n2sConfig.isAllowRadicals() && wordEls.size()==1 && indexOfFirstWord==0 && firstWord.getLocalName().equals(WORD_EL) && WordType.substituent.toString().equals(firstWord.getAttributeValue(TYPE_ATR))){
			applySimpleWordRule(wordEls, indexOfFirstWord, firstWord);
		}
		return false;
	}

	private Element getLastGroupInWordRule(Element wordEl) {
		Elements children = wordEl.getChildElements();
		Element lastChild = children.get(children.size()-1);
		while (lastChild.getChildElements().size()!=0){
			children = lastChild.getChildElements();
			lastChild = children.get(children.size()-1);
		}
		if (lastChild.getLocalName().equals(GROUP_EL)){
			return lastChild;
		}
		else{
			Elements groups = ((Element)lastChild.getParent()).getChildElements(GROUP_EL);
			if (groups.size()>0){
				return groups.get(groups.size()-1);
			}
		}
		return null;
	}

	private void applySimpleWordRule(List<Element> wordEls, int indexOfFirstWord, Element firstWord) {
		Element parentEl = (Element) firstWord.getParent();
		int indexToInsertAt = parentEl.indexOf(firstWord);
		Element wordRuleEl = new Element(WORDRULE_ATR);
		wordRuleEl.addAttribute(new Attribute(WORDRULE_ATR, WordRule.simple.toString()));//No wordRule
		wordRuleEl.addAttribute(new Attribute(TYPE_ATR, WordType.full.toString()));
		wordRuleEl.addAttribute(new Attribute(VALUE_ATR, firstWord.getAttributeValue(VALUE_ATR)));
		firstWord.detach();
		wordRuleEl.appendChild(firstWord);
		wordEls.set(indexOfFirstWord, wordRuleEl);
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
		Element assumedHyphen = new Element(HYPHEN_EL);
		assumedHyphen.appendChild("-");
		Elements substituentEls = firstWord.getChildElements(SUBSTITUENT_EL);
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
		if (WordType.full.toString().equals(wordToPotentiallyCombineWith.getAttributeValue(TYPE_ATR))){
			firstWord.getAttribute(TYPE_ATR).setValue(WordType.full.toString());
		}
		firstWord.getAttribute(VALUE_ATR).setValue(firstWord.getAttributeValue(VALUE_ATR) + wordToPotentiallyCombineWith.getAttributeValue(VALUE_ATR));
	}

	private void convertFunctionalGroupIntoGroup(Element word) throws ParsingException {
		word.getAttribute(TYPE_ATR).setValue(WordType.full.toString());
		List<Element> functionalTerms = XOMTools.getDescendantElementsWithTagName(word, FUNCTIONALTERM_EL);
		if (functionalTerms.size()!=1){
			throw new ParsingException("OPSIN Bug: Exactly 1 functionalTerm expected in functionalGroupAsGroup wordRule");
		}
		functionalTerms.get(0).setLocalName(ROOT_EL);
		Elements functionalGroups = functionalTerms.get(0).getChildElements(FUNCTIONALGROUP_EL);
		if (functionalGroups.size()!=1){
			throw new ParsingException("OPSIN Bug: Exactly 1 functionalGroup expected in functionalGroupAsGroup wordRule");
		}
		functionalGroups.get(0).setLocalName(GROUP_EL);
		functionalGroups.get(0).getAttribute(TYPE_ATR).setValue(SIMPLEGROUP_TYPE_VAL);
		functionalGroups.get(0).addAttribute(new Attribute(SUBTYPE_ATR, SIMPLEGROUP_SUBTYPE_VAL));
	}

	
	/**
	 * Checks whether the bond that will be formed will be ionic by inspection of the SMILES
	 * @param elementaryAtomEl
	 * @param functionalWord
	 * @return
	 * @throws ParsingException 
	 */
	private boolean bondWillBeIonic(Element elementaryAtomEl, Element functionalWord) throws ParsingException {
		String element1 = elementaryAtomEl.getAttributeValue(VALUE_ATR);
		if (element1.startsWith("[")){
			element1 = element1.substring(1, element1.length()-1);
		}
		List<Element> functionalGroups = XOMTools.getDescendantElementsWithTagName(functionalWord, FUNCTIONALGROUP_EL);
		if (functionalGroups.size()!=1){
			throw new ParsingException("OPSIN bug: Unable to find functional group in oxide or addition compound rule");
		}
		String smiles = functionalGroups.get(0).getAttributeValue(VALUE_ATR);
		String element2 ="";
		for (int i = 0; i <smiles.length(); i++) {
			if (Character.isUpperCase(smiles.charAt(i))){
				element2 += smiles.charAt(i);
				if (i+1 <smiles.length() && Character.isLowerCase(smiles.charAt(i +1))){
					element2 += smiles.charAt(i +1);
				}
				break;
			}
		}
		return !FragmentTools.isCovalent(element1, element2);
	}
}
