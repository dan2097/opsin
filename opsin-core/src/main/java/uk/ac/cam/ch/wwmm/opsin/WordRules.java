package uk.ac.cam.ch.wwmm.opsin;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.regex.Pattern;

import javax.xml.stream.XMLStreamConstants;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamReader;

import static uk.ac.cam.ch.wwmm.opsin.XmlDeclarations.*;

/**The rules by which words are grouped together (e.g. in functional class nomenclature)
 *
 * @author dl387
 *
 */
class WordRules {

	/**The wordRules themselves.*/
	private final List<WordRuleDescription> wordRuleList;

	/**
	 * Describes a word that a wordRule is looking for
	 * @author dl387
	 *
	 */
	private static class WordDescription {
		/**Whether the word is a full word, substituent word or functionalTerm word*/
		private final WordType type;

		/**A case insensitive pattern which attempts to match the end of the String value of the word*/
		private final Pattern endsWithPattern;

		/**The case insensitive String value of the word */
		private final String value;

		/** Only applicable for functionalTerms. The string value of the functionalTerm's type attribute*/
		private final String functionalGroupType;

		/** The value of the value attribute of the last group element in the word e.g. maybe a SMILES string*/
		private final String endsWithGroupValueAtr;
		
		/** The value of the type attribute of the last group element in the word e.g. maybe aminoAcid*/
		private final String endsWithGroupType;
		
		/** The value of the subType attribute of the last group element in the word e.g. maybe elementaryAtom*/
		private final String endsWithGroupSubType;
		
		/**
		 * Makes a description of a word to looks for
		 * @param reader
		 */
		WordDescription(XMLStreamReader reader){
			WordType type = null;
			Pattern endsWithPattern = null;
			String value = null;
			String functionalGroupType = null;
			String endsWithGroupValueAtr = null;
			String endsWithGroupType = null;
			String endsWithGroupSubType = null;
			for (int i = 0, l = reader.getAttributeCount(); i < l; i++) {
				String atrName = reader.getAttributeLocalName(i);
				String atrValue = reader.getAttributeValue(i);
				if (atrName.equals("type")){
					type = WordType.valueOf(atrValue);
				}
				else if (atrName.equals("value")){
					value = atrValue;
				}
				else if (atrName.equals("functionalGroupType")){
					functionalGroupType = atrValue;
				}
				else if (atrName.equals("endsWithRegex")){
					endsWithPattern = Pattern.compile(atrValue +"$", Pattern.CASE_INSENSITIVE);
				}
				else if (atrName.equals("endsWithGroupValueAtr")){
					endsWithGroupValueAtr = atrValue;
				}
				else if (atrName.equals("endsWithGroupType")){
					endsWithGroupType = atrValue;
				}
				else if (atrName.equals("endsWithGroupSubType")){
					endsWithGroupSubType = atrValue;
				}
			}
			if (type == null) {
				throw new RuntimeException("Malformed wordRules");
			}
			this.type = type;
			this.endsWithPattern = endsWithPattern;
			this.value = value;
			this.functionalGroupType = functionalGroupType;
			this.endsWithGroupValueAtr = endsWithGroupValueAtr;
			this.endsWithGroupType = endsWithGroupType;
			this.endsWithGroupSubType = endsWithGroupSubType;
		}

		WordType getType() {
			return type;
		}

		Pattern getEndsWithPattern() {
			return endsWithPattern;
		}

		String getValue() {
			return value;
		}

		String getFunctionalGroupType() {
			return functionalGroupType;
		}

		String getEndsWithGroupValueAtr() {
			return endsWithGroupValueAtr;
		}

		String getEndsWithGroupType() {
			return endsWithGroupType;
		}
		
		String getEndsWithGroupSubType() {
			return endsWithGroupSubType;
		}
	}

	/**
	 * A representation of a wordRule element from wordRules.xml
	 * @author dl387
	 *
	 */
	private static class WordRuleDescription {
		private final List<WordDescription> wordDescriptions;
		private final WordRule ruleName;
		private final WordType ruleType;
		
		List<WordDescription> getWordDescriptions() {
			return wordDescriptions;
		}

		WordRule getRuleName() {
			return ruleName;
		}

		WordType getRuleType() {
			return ruleType;
		}

		/**
		 * Creates a wordRule from a wordRule element found in wordRules.xml
		 * @param reader
		 * @throws XMLStreamException 
		 */
		WordRuleDescription(XMLStreamReader reader) throws XMLStreamException {
			List<WordDescription> wordDescriptions = new ArrayList<WordDescription>();
			ruleName = WordRule.valueOf(reader.getAttributeValue(null, "name"));
			ruleType = WordType.valueOf(reader.getAttributeValue(null,"type"));
			while (reader.hasNext()) {
				int event = reader.next();
				if (event == XMLStreamConstants.START_ELEMENT) {
					if (reader.getLocalName().equals("word")) {
						wordDescriptions.add(new WordDescription(reader));
					}
				}
				else if (event == XMLStreamConstants.END_ELEMENT) {
					if (reader.getLocalName().equals("wordRule")) {
						break;
					}
				}
			}
			this.wordDescriptions = Collections.unmodifiableList(wordDescriptions);
		}
	}


	/**Initialises the WordRules.
	 * @param resourceGetter
	 * @throws IOException 
	 */
	WordRules(ResourceGetter resourceGetter) throws IOException {
		List<WordRuleDescription> wordRuleList = new ArrayList<WordRuleDescription>();
		XMLStreamReader reader = resourceGetter.getXMLStreamReader("wordRules.xml");
		try {
			while (reader.hasNext()) {
				if (reader.next() == XMLStreamConstants.START_ELEMENT && 
						reader.getLocalName().equals("wordRule")) {
					wordRuleList.add(new WordRuleDescription(reader));
				}
			}
		}
		catch (XMLStreamException e) {
			throw new IOException("Parsing exception occurred while reading wordRules.xml", e);
		}
		finally {
			try {
				reader.close();
			} catch (XMLStreamException e) {
				throw new IOException("Parsing exception occurred while reading wordRules.xml", e);
			}
		}
		this.wordRuleList = Collections.unmodifiableList(wordRuleList);
	}

	/**Takes a molecule element and places the word elements into wordRule elements
	 * @param n2sConfig 
	 *
	 * @param moleculeEl A molecule element with word children
	 * @param allowSpaceRemoval 
	 * @throws ParsingException
	 */
	void groupWordsIntoWordRules(NameToStructureConfig n2sConfig, Element moleculeEl, boolean allowSpaceRemoval) throws ParsingException {
		List<Element> wordEls = OpsinTools.getChildElementsWithTagName(moleculeEl, WORD_EL);
		//note that multiple words in wordEls may be later replaced by a wordRule element
		for (int i = 0; i <wordEls.size(); i++) {
			if (matchWordRule(n2sConfig, wordEls, i, allowSpaceRemoval)){
				i=-1;//if function did something
			}
		}
		List<Element> wordRuleEls = moleculeEl.getChildElements();
		for (int i = 0; i < wordRuleEls.size(); i++) {
			Element wordRuleEl = wordRuleEls.get(i);
			if (!wordRuleEl.getName().equals(WORDRULE_EL)){
				throw new ParsingException("Unable to assign wordRule to: " + wordRuleEl.getValue());
			}
		}
	}

	private boolean matchWordRule(NameToStructureConfig n2sConfig, List<Element> wordEls, int indexOfFirstWord, boolean allowSpaceRemoval) throws ParsingException {
		wordRuleLoop: for (WordRuleDescription wordRuleDesc : wordRuleList) {
			int i =indexOfFirstWord;
			List<WordDescription> wordDescriptions = wordRuleDesc.getWordDescriptions();
			int wordsInWordRule = wordDescriptions.size();
			if (i + wordsInWordRule -1 < wordEls.size()){//need sufficient words to match the word rule
				for (int j = 0; j < wordsInWordRule; j++) {
					Element wordEl = wordEls.get(i+j);
					WordDescription wd = wordDescriptions.get(j);
					if (!wd.getType().toString().equals(wordEl.getAttributeValue(TYPE_ATR))){
						continue wordRuleLoop;//type mismatch;
					}
					if (wd.getValue() !=null && !wordEl.getAttributeValue(VALUE_ATR).toLowerCase().equals(wd.getValue())){//word string contents mismatch
						continue wordRuleLoop;
					}
					if (wd.getFunctionalGroupType() !=null){
						if (WordType.functionalTerm.toString().equals(wordEl.getAttributeValue(TYPE_ATR))){
							List<Element> children = wordEl.getChildElements();
							Element lastChild = children.get(children.size()-1);
							while (lastChild.getChildCount() != 0){
								children = lastChild.getChildElements();
								lastChild = children.get(children.size()-1);
							}
							if (lastChild.getName().equals(CLOSEBRACKET_EL)){
								lastChild = OpsinTools.getPreviousSibling(lastChild);
							}
							if (lastChild==null){
								throw new ParsingException("OPSIN Bug: Cannot find the functional element in a functionalTerm");
							}
							if (!wd.getFunctionalGroupType().equals(lastChild.getAttributeValue(TYPE_ATR))){
								continue wordRuleLoop;
							}
						}
					}
					if (wd.getEndsWithPattern() !=null){
						if (!wd.getEndsWithPattern().matcher(wordEl.getAttributeValue(VALUE_ATR)).find()){
							continue wordRuleLoop;
						}
					}
					if (wd.getEndsWithGroupValueAtr() !=null){
						Element lastGroupInWordRule = getLastGroupInWordRule(wordEl);
						if (lastGroupInWordRule==null || !wd.getEndsWithGroupValueAtr().equals(lastGroupInWordRule.getAttributeValue(VALUE_ATR))){
							continue wordRuleLoop;
						}
					}
					if (wd.getEndsWithGroupType() !=null){
						Element lastGroupInWordRule = getLastGroupInWordRule(wordEl);
						if (lastGroupInWordRule==null || !wd.getEndsWithGroupType().equals(lastGroupInWordRule.getAttributeValue(TYPE_ATR))){
							continue wordRuleLoop;
						}
					}
					if (wd.getEndsWithGroupSubType() !=null){
						Element lastGroupInWordRule = getLastGroupInWordRule(wordEl);
						if (lastGroupInWordRule==null || !wd.getEndsWithGroupSubType().equals(lastGroupInWordRule.getAttributeValue(SUBTYPE_ATR))){
							continue wordRuleLoop;
						}
					}
				}
				//Word Rule matches!
				Element wordRuleEl = new Element(WORDRULE_EL);
				wordRuleEl.addAttribute(new Attribute(TYPE_ATR, wordRuleDesc.getRuleType().toString()));
				wordRuleEl.addAttribute(new Attribute(WORDRULE_EL, wordRuleDesc.getRuleName().toString()));

				/*
				 * Some wordRules can not be entirely processed at the structure building stage
				 */
				WordRule wordRule = wordRuleDesc.getRuleName();
				if (wordRule == WordRule.functionGroupAsGroup){//convert the functional term into a full term
					Element functionalWord = wordEls.get(i + wordsInWordRule -1);
					if (!functionalWord.getAttributeValue(TYPE_ATR).equals(FUNCTIONALTERM_EL) || wordsInWordRule>2){
						throw new ParsingException("OPSIN bug: Problem with functionGroupAsGroup wordRule");
					}
					convertFunctionalGroupIntoGroup(functionalWord);
					if (wordsInWordRule==2){
						joinWords(wordEls, i, wordEls.get(i), functionalWord);
						wordsInWordRule =1;
					}
					wordRuleEl.getAttribute(WORDRULE_ATR).setValue(WordRule.simple.toString());
				}
				else if (wordRule == WordRule.carbonylDerivative || wordRule == WordRule.acidReplacingFunctionalGroup){//e.g. acetone 4,4-diphenylsemicarbazone. This is better expressed as a full word as the substituent actually locants onto the functional term
					if (wordsInWordRule==3){//substituent present
						joinWords(wordEls, i+1, wordEls.get(i+1), wordEls.get(i+2));
						wordsInWordRule--;
						List<Element> functionalTerm = OpsinTools.getDescendantElementsWithTagName(wordEls.get(i+1), FUNCTIONALTERM_EL);//rename functionalTerm element to root
						if (functionalTerm.size()!=1){
							throw new ParsingException("OPSIN bug: Problem with "+ wordRule +" wordRule");
						}
						functionalTerm.get(0).setName(ROOT_EL);
						List<Element> functionalGroups = OpsinTools.getDescendantElementsWithTagName(functionalTerm.get(0), FUNCTIONALGROUP_EL);//rename functionalGroup element to group
						if (functionalGroups.size()!=1){
							throw new ParsingException("OPSIN bug: Problem with "+ wordRule +" wordRule");
						}
						functionalGroups.get(0).setName(GROUP_EL);
						wordEls.get(i+1).getAttribute(TYPE_ATR).setValue(WordType.full.toString());
					}
				}
				else if (wordRule == WordRule.additionCompound || wordRule == WordRule.oxide){//is the halide/pseudohalide/oxide actually a counterion rather than covalently bonded
					Element possibleElementaryAtom = wordEls.get(i);
					List<Element> elementaryAtoms = OpsinTools.getDescendantElementsWithTagNameAndAttribute(possibleElementaryAtom, GROUP_EL, SUBTYPE_ATR, ELEMENTARYATOM_SUBTYPE_VAL);
					if (elementaryAtoms.size()==1){
						for (int j = 1; j < wordsInWordRule; j++) {
							if (bondWillBeIonic(elementaryAtoms.get(0), wordEls.get(i+j))){//use separate word rules for ionic components
								continue wordRuleLoop;
							}
						}
					}
				}



				List<String> wordValues = new ArrayList<String>();
				Element parentEl = wordEls.get(i).getParent();
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
		if (firstWord.getName().equals(WORD_EL) && WordType.full.toString().equals(firstWord.getAttributeValue(TYPE_ATR))){//No wordRule -->wordRule="simple"
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
		if (n2sConfig.isAllowRadicals() && wordEls.size()==1 && indexOfFirstWord==0 && firstWord.getName().equals(WORD_EL) && WordType.substituent.toString().equals(firstWord.getAttributeValue(TYPE_ATR))){
			applySubstituentWordRule(wordEls, indexOfFirstWord, firstWord);
		}
		return false;
	}

	private Element getLastGroupInWordRule(Element wordEl) {
		List<Element> children = wordEl.getChildElements();
		Element lastChild = children.get(children.size()-1);
		while (lastChild.getChildCount() != 0){
			children = lastChild.getChildElements();
			lastChild = children.get(children.size()-1);
		}
		if (lastChild.getName().equals(GROUP_EL)){
			return lastChild;
		}
		else{
			List<Element> groups = lastChild.getParent().getChildElements(GROUP_EL);
			if (groups.size() > 0){
				return groups.get(groups.size() - 1);
			}
		}
		return null;
	}

	private void applySimpleWordRule(List<Element> wordEls, int indexOfFirstWord, Element firstWord) {
		Element parentEl = firstWord.getParent();
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
	

	private void applySubstituentWordRule(List<Element> wordEls, int indexOfFirstWord, Element firstWord) {
		Element parentEl = firstWord.getParent();
		int indexToInsertAt = parentEl.indexOf(firstWord);
		Element wordRuleEl = new Element(WORDRULE_ATR);
		wordRuleEl.addAttribute(new Attribute(WORDRULE_ATR, WordRule.substituent.toString()));
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
		List<Element> substituentEls = firstWord.getChildElements(SUBSTITUENT_EL);
		if (substituentEls.size()==0){
			throw new ParsingException("OPSIN Bug: Substituent element not found where substituent element expected");
		}
		Element finalSubstituent = substituentEls.get(substituentEls.size()-1);
		List<Element> finalSubstituentChildren = finalSubstituent.getChildElements();
		if (!finalSubstituentChildren.get(finalSubstituentChildren.size()-1).getName().equals(HYPHEN_EL)){//add an implicit hyphen if one is not already present
			Element implicitHyphen = new Element(HYPHEN_EL);
			implicitHyphen.setValue("-");
			finalSubstituent.appendChild(implicitHyphen);
		}
		List<Element> elementsToMergeIntoSubstituent = wordToPotentiallyCombineWith.getChildElements();
		for (int j =  elementsToMergeIntoSubstituent.size() -1 ; j >=0; j--) {
			Element el = elementsToMergeIntoSubstituent.get(j);
			el.detach();
			OpsinTools.insertAfter(finalSubstituent, el);
		}
		if (WordType.full.toString().equals(wordToPotentiallyCombineWith.getAttributeValue(TYPE_ATR))){
			firstWord.getAttribute(TYPE_ATR).setValue(WordType.full.toString());
		}
		firstWord.getAttribute(VALUE_ATR).setValue(firstWord.getAttributeValue(VALUE_ATR) + wordToPotentiallyCombineWith.getAttributeValue(VALUE_ATR));
	}

	private void convertFunctionalGroupIntoGroup(Element word) throws ParsingException {
		word.getAttribute(TYPE_ATR).setValue(WordType.full.toString());
		List<Element> functionalTerms = OpsinTools.getDescendantElementsWithTagName(word, FUNCTIONALTERM_EL);
		if (functionalTerms.size()!=1){
			throw new ParsingException("OPSIN Bug: Exactly 1 functionalTerm expected in functionalGroupAsGroup wordRule");
		}
		functionalTerms.get(0).setName(ROOT_EL);
		List<Element> functionalGroups = functionalTerms.get(0).getChildElements(FUNCTIONALGROUP_EL);
		if (functionalGroups.size()!=1){
			throw new ParsingException("OPSIN Bug: Exactly 1 functionalGroup expected in functionalGroupAsGroup wordRule");
		}
		functionalGroups.get(0).setName(GROUP_EL);
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
		List<Element> functionalGroups = OpsinTools.getDescendantElementsWithTagName(functionalWord, FUNCTIONALGROUP_EL);
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
