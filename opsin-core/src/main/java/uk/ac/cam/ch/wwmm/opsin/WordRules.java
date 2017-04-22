package uk.ac.cam.ch.wwmm.opsin;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Locale;
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
	
	enum EndsWithGroup {
		acid,
		ateGroup;
	}
	
	private static final Pattern icOrOusAcid = Pattern.compile("(ic|ous)([ ]?acid)?$");
	private static final Pattern ateOrIteOrAmide = Pattern.compile("(at|it|amid)e?$");

	/**
	 * Describes a word that a wordRule is looking for
	 * @author dl387
	 *
	 */
	private static class WordDescription {
		/**Whether the word is a full word, substituent word or functionalTerm word*/
		private final WordType type;
		
		/**A group with a hardcoded method for efficient detection */
		private final EndsWithGroup endsWithGroup;
		
		/**A case insensitive pattern which attempts to match the end of the String value of the word*/
		private final Pattern endsWithPattern;

		/**The case insensitive String value of the word */
		private final String value;

		/** Only applicable for functionalTerms. The string value of the functionalTerm's type attribute*/
		private final String functionalGroupType;
		
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
			String value = null;
			EndsWithGroup endsWithGroup = null;
			Pattern endsWithPattern = null;
			String functionalGroupType = null;
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
				else if (atrName.equals("endsWith")){
					endsWithGroup = EndsWithGroup.valueOf(atrValue);
				}
				else if (atrName.equals("endsWithRegex")){
					endsWithPattern = Pattern.compile(atrValue +"$", Pattern.CASE_INSENSITIVE);
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
			this.endsWithGroup = endsWithGroup;
			this.endsWithPattern = endsWithPattern;
			this.value = value;
			this.functionalGroupType = functionalGroupType;
			this.endsWithGroupType = endsWithGroupType;
			this.endsWithGroupSubType = endsWithGroupSubType;
		}

		WordType getType() {
			return type;
		}
		
		EndsWithGroup getEndsWithGroup() {
			return endsWithGroup;
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
	 * @param moleculeEl A molecule element with word children
	 * @param n2sConfig 
	 * @param allowSpaceRemoval 
	 * @param componentRatios 
	 * @throws ParsingException
	 */
	void groupWordsIntoWordRules(Element moleculeEl, NameToStructureConfig n2sConfig, boolean allowSpaceRemoval, Integer[] componentRatios) throws ParsingException {
		WordRulesInstance instance = new WordRulesInstance(moleculeEl, n2sConfig, allowSpaceRemoval, componentRatios);
		List<Element> wordEls = moleculeEl.getChildElements(WORD_EL);
		//note that multiple words in wordEls may be later replaced by a wordRule element
		for (int i = 0; i <wordEls.size(); i++) {
			if (instance.matchWordRule(wordEls, i)) {
				i=-1;//if function did something
			}
		}
		List<Element> wordRuleEls = moleculeEl.getChildElements();
		for (Element wordRuleEl : wordRuleEls) {
			if (!wordRuleEl.getName().equals(WORDRULE_EL)){
				throw new ParsingException("Unable to assign wordRule to: " + wordRuleEl.getAttributeValue(VALUE_ATR));
			}
		}
	}
	
	private class WordRulesInstance {
		private final Element moleculeEl;
		private final boolean allowRadicals;
		private final boolean allowSpaceRemoval;
		private final Integer expectedNumOfComponents;
		
		WordRulesInstance(Element moleculeEl, NameToStructureConfig n2sConfig, boolean allowSpaceRemoval, Integer[] componentRatios) {
			this.moleculeEl = moleculeEl;
			this.allowRadicals = n2sConfig.isAllowRadicals();
			this.allowSpaceRemoval = allowSpaceRemoval;
			this.expectedNumOfComponents = componentRatios != null ? componentRatios.length : null;
		}
		
		private boolean matchWordRule(List<Element> wordEls, int indexOfFirstWord) throws ParsingException {
			wordRuleLoop: for (WordRuleDescription wordRuleDesc : wordRuleList) {
				int i = indexOfFirstWord;
				List<WordDescription> wordDescriptions = wordRuleDesc.getWordDescriptions();
				int wordsInWordRule = wordDescriptions.size();
				if (i + wordsInWordRule <= wordEls.size()) {//need sufficient words to match the word rule
					for (int j = 0; j < wordsInWordRule; j++) {
						Element wordEl = wordEls.get(i + j);
						WordDescription wd = wordDescriptions.get(j);
						if (!wd.getType().toString().equals(wordEl.getAttributeValue(TYPE_ATR))){
							continue wordRuleLoop;//type mismatch;
						}

						String functionalGroupTypePredicate = wd.getFunctionalGroupType();
						if (functionalGroupTypePredicate != null) {
							if (!WordType.functionalTerm.toString().equals(wordEl.getAttributeValue(TYPE_ATR))){
								continue wordRuleLoop;
							}
							Element lastEl = getLastElementInWord(wordEl);
							if (lastEl == null) {
								throw new ParsingException("OPSIN Bug: Cannot find the functional element in a functionalTerm");
							}
							while (lastEl.getName().equals(CLOSEBRACKET_EL) || lastEl.getName().equals(STRUCTURALCLOSEBRACKET_EL)) {
								lastEl = OpsinTools.getPreviousSibling(lastEl);
								if (lastEl == null) {
									throw new ParsingException("OPSIN Bug: Cannot find the functional element in a functionalTerm");
								}
							}
							if (!functionalGroupTypePredicate.equals(lastEl.getAttributeValue(TYPE_ATR))) {
								continue wordRuleLoop;
							}
						}
						
						EndsWithGroup endsWithGroupPredicate = wd.getEndsWithGroup();
						if (endsWithGroupPredicate != null && !endsWithGroupPredicateSatisfied(wordEl, endsWithGroupPredicate)) {
							continue wordRuleLoop;
						}

						String valuePredicate = wd.getValue();
						if (valuePredicate != null && !wordEl.getAttributeValue(VALUE_ATR).toLowerCase(Locale.ROOT).equals(valuePredicate)){//word string contents mismatch
							continue wordRuleLoop;
						}

						Pattern endsWithPatternPredicate = wd.getEndsWithPattern();
						if (endsWithPatternPredicate != null) {
							if (!endsWithPatternPredicate.matcher(wordEl.getAttributeValue(VALUE_ATR)).find()){
								continue wordRuleLoop;
							}
						}

						String endsWithGroupTypePredicate = wd.getEndsWithGroupType();
						if (endsWithGroupTypePredicate != null) {
							Element lastGroupInWordRule = getLastGroupInWordRule(wordEl);
							if (lastGroupInWordRule == null || !endsWithGroupTypePredicate.equals(lastGroupInWordRule.getAttributeValue(TYPE_ATR))){
								continue wordRuleLoop;
							}
						}

						String endsWithSubGroupTypePredicate = wd.getEndsWithGroupSubType();
						if (endsWithSubGroupTypePredicate != null) {
							Element lastGroupInWordRule = getLastGroupInWordRule(wordEl);
							if (lastGroupInWordRule == null || !endsWithSubGroupTypePredicate.equals(lastGroupInWordRule.getAttributeValue(SUBTYPE_ATR))){
								continue wordRuleLoop;
							}
						}
					}
					//Word Rule matches!
					Element wordRuleEl = new GroupingEl(WORDRULE_EL);
					WordRule wordRule = wordRuleDesc.getRuleName();
					wordRuleEl.addAttribute(new Attribute(TYPE_ATR, wordRuleDesc.getRuleType().toString()));
					wordRuleEl.addAttribute(new Attribute(WORDRULE_EL, wordRule.toString()));

					/*
					 * Some wordRules can not be entirely processed at the structure building stage
					 */
					switch (wordRule) {
					case functionGroupAsGroup:
						//convert the functional term into a full term
						Element functionalWord = wordEls.get(i + wordsInWordRule -1);
						if (!functionalWord.getAttributeValue(TYPE_ATR).equals(FUNCTIONALTERM_EL) || wordsInWordRule>2){
							throw new ParsingException("OPSIN bug: Problem with functionGroupAsGroup wordRule");
						}
						convertFunctionalGroupIntoGroup(functionalWord);
						if (wordsInWordRule==2){
							joinWords(wordEls, wordEls.get(i), functionalWord);
							wordsInWordRule =1;
						}
						wordRuleEl.getAttribute(WORDRULE_ATR).setValue(WordRule.simple.toString());
						break;
					case carbonylDerivative:
					case acidReplacingFunctionalGroup:
						//e.g. acetone 4,4-diphenylsemicarbazone. This is better expressed as a full word as the substituent actually locants onto the functional term
						for (int j = 1; j < (wordsInWordRule - 1); j++) {
							Element wordEl = wordEls.get(i + j);
							if (WordType.substituent.toString().equals(wordEl.getAttributeValue(TYPE_ATR))) {
								joinWords(wordEls, wordEls.get(i + j), wordEls.get(i + j + 1));
								wordsInWordRule--;
								List<Element> functionalTerm = OpsinTools.getDescendantElementsWithTagName(wordEls.get(i + j), FUNCTIONALTERM_EL);//rename functionalTerm element to root
								if (functionalTerm.size() != 1){
									throw new ParsingException("OPSIN bug: Problem with "+ wordRule +" wordRule");
								}
								functionalTerm.get(0).setName(ROOT_EL);
								List<Element> functionalGroups = OpsinTools.getDescendantElementsWithTagName(functionalTerm.get(0), FUNCTIONALGROUP_EL);//rename functionalGroup element to group
								if (functionalGroups.size() != 1){
									throw new ParsingException("OPSIN bug: Problem with "+ wordRule +" wordRule");
								}
								functionalGroups.get(0).setName(GROUP_EL);
								wordEls.get(i + j).getAttribute(TYPE_ATR).setValue(WordType.full.toString());
							}
						}
						break;
					case additionCompound:
					case oxide:
						//is the halide/pseudohalide/oxide actually a counterion rather than covalently bonded
						Element possibleElementaryAtomContainingWord = wordEls.get(i);
						List<Element> elementaryAtoms = OpsinTools.getDescendantElementsWithTagNameAndAttribute(possibleElementaryAtomContainingWord, GROUP_EL, SUBTYPE_ATR, ELEMENTARYATOM_SUBTYPE_VAL);
						if (elementaryAtoms.size() == 1) {
							Element elementaryAtom = elementaryAtoms.get(0);
							ChemEl chemEl1 = getChemElFromElementaryAtomEl(elementaryAtom);
							if (wordRule == WordRule.oxide) {
								if (wordsInWordRule != 2){
									throw new ParsingException("OPSIN bug: Problem with "+ wordRule +" wordRule");
								}
								Element oxideWord = wordEls.get(i + 1);
								ChemEl chemEl2 = getChemElFromWordWithFunctionalGroup(oxideWord);
								if (!FragmentTools.isCovalent(chemEl1, chemEl2)){
									Element oxideGroup = convertFunctionalGroupIntoGroup(oxideWord);
									setOxideStructureAppropriately(oxideGroup, elementaryAtom);
									applySimpleWordRule(wordEls, indexOfFirstWord, possibleElementaryAtomContainingWord);
									continue wordRuleLoop;
								}
							}
							else {
								for (int j = 1; j < wordsInWordRule; j++) {
									Element functionalGroup = wordEls.get(i + j);
									ChemEl chemEl2 = getChemElFromWordWithFunctionalGroup(functionalGroup);
									if (!FragmentTools.isCovalent(chemEl1, chemEl2)) {//use separate word rules for ionic components
										boolean specialCaseCovalency = false;
										if (chemEl2.isHalogen() && wordsInWordRule == 2) {
											switch (chemEl1) {
											case Mg:
												if (possibleElementaryAtomContainingWord.getChildCount() > 1) {
													//treat grignards (i.e. substitutedmagnesium halides) as covalent
													specialCaseCovalency = true;
												}
												break;
											case Al:
												if (chemEl2 == ChemEl.Cl || chemEl2 == ChemEl.Br || chemEl2 == ChemEl.I) {
													specialCaseCovalency = true;
												}
												break;
											case Ti:
												if (oxidationNumberOrMultiplierIs(elementaryAtom, functionalGroup, 4) && 
														(chemEl2 == ChemEl.Cl || chemEl2 == ChemEl.Br || chemEl2 == ChemEl.I)) {
													specialCaseCovalency = true;
												}
												break;
											case V:
												if (oxidationNumberOrMultiplierIs(elementaryAtom, functionalGroup, 4) && 
														chemEl2 == ChemEl.Cl) {
													specialCaseCovalency = true;
												}
												break;
											case Zr:
											case Hf:
												if (oxidationNumberOrMultiplierIs(elementaryAtom, functionalGroup, 4) && 
														chemEl2 == ChemEl.Br) {
													specialCaseCovalency = true;
												}
												break;
											case U:
												if (oxidationNumberOrMultiplierIs(elementaryAtom, functionalGroup, 6) && 
														(chemEl2 == ChemEl.F || chemEl2 == ChemEl.Cl)) {
													specialCaseCovalency = true;
												}
												break;
											case Np:
											case Pu:
												if (oxidationNumberOrMultiplierIs(elementaryAtom, functionalGroup, 6) && 
														chemEl2 == ChemEl.F) {
													specialCaseCovalency = true;
												}
												break;
											default:
												break;
											}
										}
										else if (chemEl2 == ChemEl.H && wordsInWordRule == 2) {
											if (chemEl1 == ChemEl.Al) {
												//aluminium hydrides are covalent
												specialCaseCovalency = true;
											}
										}
										if (!specialCaseCovalency) {
											continue wordRuleLoop;
										}
									}
								}
							}
						}
						break;
					case potentialAlcoholEster:
						if (expectedNumOfComponents != null && expectedNumOfComponents == moleculeEl.getChildCount()) {
							//don't apply this wordRule if doing so makes the number of components incorrect
							continue wordRuleLoop;
						}
						break;
					default:
						break;
					}

					List<String> wordValues = new ArrayList<String>();
					Element parentEl = wordEls.get(i).getParent();
					int indexToInsertAt = parentEl.indexOf(wordEls.get(i));
					for (int j = 0; j < wordsInWordRule; j++) {
						Element wordEl = wordEls.remove(i);
						wordEl.detach();
						wordRuleEl.addChild(wordEl);
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
						joinWords(wordEls, firstWord, wordToPotentiallyCombineWith);
						return true;
					}
				}
			}
			else if (WordType.functionalTerm.toString().equals(firstWord.getAttributeValue(TYPE_ATR)) &&  firstWord.getAttributeValue(VALUE_ATR).toLowerCase(Locale.ROOT).equals("salt")) {
				wordEls.remove(indexOfFirstWord);
				firstWord.detach();
				if (moleculeEl.getAttribute(ISSALT_ATR) == null) {
					moleculeEl.addAttribute(ISSALT_ATR, "yes");
				}
				return true;
			}
			if (allowRadicals && wordEls.size() == 1 && indexOfFirstWord == 0 && firstWord.getName().equals(WORD_EL) && WordType.substituent.toString().equals(firstWord.getAttributeValue(TYPE_ATR))){
				//name is all one substituent, make this a substituent and finish
				applySubstituentWordRule(wordEls, indexOfFirstWord, firstWord);
			}
			return false;
		}

		private boolean endsWithGroupPredicateSatisfied(Element wordEl, EndsWithGroup endsWithGroupPredicate) throws ParsingException {
			Element lastEl = getLastElementInWord(wordEl);
			if (lastEl == null) {
				return false;
			}
			String elName = lastEl.getName();
			while (elName.equals(CLOSEBRACKET_EL) || 
					elName.equals(STRUCTURALCLOSEBRACKET_EL) ||
					elName.equals(ISOTOPESPECIFICATION_EL)) {
				lastEl = OpsinTools.getPreviousSibling(lastEl);
				if (lastEl == null) {
					return false;
				}
				elName = lastEl.getName();
			}

			if (endsWithGroupPredicate == EndsWithGroup.acid) {
				if (elName.equals(SUFFIX_EL)) {
					if (icOrOusAcid.matcher(lastEl.getAttributeValue(VALUE_ATR)).find()) {
						return true;
					}
				}
				else if (elName.equals(GROUP_EL)) {
					if (lastEl.getAttribute(FUNCTIONALIDS_ATR) != null && icOrOusAcid.matcher(lastEl.getValue()).find()) {
						return true;
					}
				}
			}
			else if (endsWithGroupPredicate == EndsWithGroup.ateGroup) {
				if (elName.equals(GROUP_EL)) {
					if (lastEl.getAttribute(FUNCTIONALIDS_ATR) != null && ateOrIteOrAmide.matcher(lastEl.getValue()).find()) {
						return true;
					}
				}
				else {
					while (lastEl != null && elName.equals(SUFFIX_EL)) {
						String suffixValAtr = lastEl.getAttributeValue(VALUE_ATR);
						if (ateOrIteOrAmide.matcher(suffixValAtr).find() || suffixValAtr.equals("glycoside")) {
							return true;
						}
						//glycoside is not always the last suffix
						lastEl = OpsinTools.getPreviousSibling(lastEl, SUFFIX_EL);
					}
				}
			}
			return false;
		}

		private boolean oxidationNumberOrMultiplierIs(Element elementaryAtomEl, Element functionalGroupWord, int expectedVal) throws ParsingException {
			List<Element> functionalGroups = OpsinTools.getDescendantElementsWithTagName(functionalGroupWord, FUNCTIONALGROUP_EL);
			if (functionalGroups.size() != 1) {
				throw new ParsingException("OPSIN bug: Unable to find functional group in oxide or addition compound rule");
			}
			Element possibleMultiplier = OpsinTools.getPreviousSibling(functionalGroups.get(0));
			if (possibleMultiplier != null && possibleMultiplier.getName().equals(MULTIPLIER_EL)) {
				return Integer.parseInt(possibleMultiplier.getAttributeValue(VALUE_ATR)) == expectedVal;
					
			}
			else {
				Element possibleOxidationNumber = OpsinTools.getNextSibling(elementaryAtomEl);
				if(possibleOxidationNumber != null && possibleOxidationNumber.getName().equals(OXIDATIONNUMBERSPECIFIER_EL)) {
					return Integer.parseInt(possibleOxidationNumber.getAttributeValue(VALUE_ATR)) == expectedVal;
				}
			}
			return false;
		}

		private Element getLastGroupInWordRule(Element wordEl) {
			Element lastEl = getLastElementInWord(wordEl);
			if (lastEl.getName().equals(GROUP_EL)) {
				return lastEl;
			}
			else{
				List<Element> groups = lastEl.getParent().getChildElements(GROUP_EL);
				if (groups.size() > 0) {
					return groups.get(groups.size() - 1);
				}
			}
			return null;
		}
		
		private Element getLastElementInWord(Element wordEl) {
			List<Element> children = wordEl.getChildElements();
			Element lastChild = children.get(children.size() - 1);
			while (lastChild.getChildCount() != 0) {
				children = lastChild.getChildElements();
				lastChild = children.get(children.size() - 1);
			}
			return lastChild;
		}

		private void applySimpleWordRule(List<Element> wordEls, int indexOfFirstWord, Element firstWord) {
			Element parentEl = firstWord.getParent();
			int indexToInsertAt = parentEl.indexOf(firstWord);
			Element wordRuleEl = new GroupingEl(WORDRULE_EL);
			wordRuleEl.addAttribute(new Attribute(WORDRULE_ATR, WordRule.simple.toString()));//No wordRule
			wordRuleEl.addAttribute(new Attribute(TYPE_ATR, WordType.full.toString()));
			wordRuleEl.addAttribute(new Attribute(VALUE_ATR, firstWord.getAttributeValue(VALUE_ATR)));
			firstWord.detach();
			wordRuleEl.addChild(firstWord);
			wordEls.set(indexOfFirstWord, wordRuleEl);
			parentEl.insertChild(wordRuleEl, indexToInsertAt);
		}
		

		private void applySubstituentWordRule(List<Element> wordEls, int indexOfFirstWord, Element firstWord) {
			Element parentEl = firstWord.getParent();
			int indexToInsertAt = parentEl.indexOf(firstWord);
			Element wordRuleEl = new GroupingEl(WORDRULE_EL);
			wordRuleEl.addAttribute(new Attribute(WORDRULE_ATR, WordRule.substituent.toString()));
			wordRuleEl.addAttribute(new Attribute(TYPE_ATR, WordType.full.toString()));
			wordRuleEl.addAttribute(new Attribute(VALUE_ATR, firstWord.getAttributeValue(VALUE_ATR)));
			firstWord.detach();
			wordRuleEl.addChild(firstWord);
			wordEls.set(indexOfFirstWord, wordRuleEl);
			parentEl.insertChild(wordRuleEl, indexToInsertAt);
		}

		/**
		 * Merges two adjacent words
		 * The latter word (wordToPotentiallyCombineWith) is merged into the former and removed from wordEls
		 * @param wordEls
		 * @param firstWord
		 * @param wordToPotentiallyCombineWith
		 * @throws ParsingException
		 */
		private void joinWords(List<Element> wordEls, Element firstWord, Element wordToPotentiallyCombineWith) throws ParsingException {
			wordEls.remove(wordToPotentiallyCombineWith);
			wordToPotentiallyCombineWith.detach();
			List<Element> substituentEls = firstWord.getChildElements(SUBSTITUENT_EL);
			if (substituentEls.size()==0){
				throw new ParsingException("OPSIN Bug: Substituent element not found where substituent element expected");
			}
			Element finalSubstituent = substituentEls.get(substituentEls.size() - 1);
			List<Element> finalSubstituentChildren = finalSubstituent.getChildElements();
			if (!finalSubstituentChildren.get(finalSubstituentChildren.size() - 1).getName().equals(HYPHEN_EL)){//add an implicit hyphen if one is not already present
				Element implicitHyphen = new TokenEl(HYPHEN_EL, "-");
				finalSubstituent.addChild(implicitHyphen);
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

		private Element convertFunctionalGroupIntoGroup(Element word) throws ParsingException {
			word.getAttribute(TYPE_ATR).setValue(WordType.full.toString());
			List<Element> functionalTerms = OpsinTools.getDescendantElementsWithTagName(word, FUNCTIONALTERM_EL);
			if (functionalTerms.size() != 1){
				throw new ParsingException("OPSIN Bug: Exactly 1 functionalTerm expected in functionalGroupAsGroup wordRule");
			}
			Element functionalTerm = functionalTerms.get(0);
			functionalTerm.setName(ROOT_EL);
			List<Element> functionalGroups = functionalTerm.getChildElements(FUNCTIONALGROUP_EL);
			if (functionalGroups.size() != 1){
				throw new ParsingException("OPSIN Bug: Exactly 1 functionalGroup expected in functionalGroupAsGroup wordRule");
			}
			Element functionalGroup = functionalGroups.get(0);
			functionalGroup.setName(GROUP_EL);
			functionalGroup.getAttribute(TYPE_ATR).setValue(SIMPLEGROUP_TYPE_VAL);
			functionalGroup.addAttribute(new Attribute(SUBTYPE_ATR, SIMPLEGROUP_SUBTYPE_VAL));
			return functionalGroup;
		}
		

		/**
		 * Sets the SMILES of the oxide group to be something like [O-2]
		 * ... unless the oxide group is multiplied and the elementaryAtom has no oxidation states greater 2
		 * in which case [O-][O-] would be assumed
		 * @param oxideGroup
		 * @param elementaryAtom
		 */
		private void setOxideStructureAppropriately(Element oxideGroup, Element elementaryAtom) {
			boolean chainInterpretation = false;
			Integer multiplierVal = null;
			Element possibleMultiplier = OpsinTools.getPreviousSibling(oxideGroup);
			if (possibleMultiplier != null && 
					possibleMultiplier.getName().equals(MULTIPLIER_EL)){
				multiplierVal = Integer.parseInt(possibleMultiplier.getAttributeValue(VALUE_ATR));
				if (multiplierVal > 1) {
					String commonOxidationStatesAndMax = elementaryAtom.getAttributeValue(COMMONOXIDATIONSTATESANDMAX_ATR);
					if (commonOxidationStatesAndMax == null ||
							Integer.parseInt(commonOxidationStatesAndMax.split(":")[1]) <= 2){
						chainInterpretation = true;
					}
				}
			}

			Attribute value = oxideGroup.getAttribute(VALUE_ATR);
			String smiles = value.getValue();
			String element;
			if (smiles.equals("O")){
				element = "O";
			}
			else if (smiles.equals("S")){
				element = "S";
			}
			else if (smiles.startsWith("[Se")){
				element = "Se";
			}
			else if (smiles.startsWith("[Te")){
				element = "Te";
			}
			else{
				throw new RuntimeException("OPSIN Bug: Unexpected smiles for oxideGroup: " + smiles);
			}
			if (chainInterpretation){
				StringBuilder sb = new StringBuilder();
				sb.append('[');
				sb.append(element);
				sb.append("-]");
				for (int i = 2; i < multiplierVal; i++) {
					sb.append('[');
					sb.append(element);
					sb.append(']');
				}
				sb.append('[');
				sb.append(element);
				sb.append("-]");
				value.setValue(sb.toString());
				possibleMultiplier.detach();
			}
			else{
				value.setValue("[" + element + "-2]");
			}
		}
		

		private ChemEl getChemElFromElementaryAtomEl(Element elementaryAtomEl)  {
			String elementStr = elementaryAtomEl.getAttributeValue(VALUE_ATR);
			if (elementStr.startsWith("[")) {
				int len = elementStr.length() - 1;
				for (int i = 1; i < len; i++) {
					char ch = elementStr.charAt(i);
					if ((ch >= 'A' && ch <='Z') || (ch >= 'a' && ch <='z')) {
						if (i + 1 < len) {
							char ch2 = elementStr.charAt(i + 1);
							if ((ch2 >= 'A' && ch2 <='Z') || (ch2 >= 'a' && ch2 <='z')) {
								//two letter element
								elementStr = elementStr.substring(i, i + 2);
								break;
							}
						}
						//one letter element
						elementStr = elementStr.substring(i, i + 1);
						break;
					}
				}
			}
			return ChemEl.valueOf(elementStr);
		}
		
		private ChemEl getChemElFromWordWithFunctionalGroup(Element functionalWord) throws ParsingException  {
			List<Element> functionalGroups = OpsinTools.getDescendantElementsWithTagName(functionalWord, FUNCTIONALGROUP_EL);
			if (functionalGroups.size() != 1){
				throw new ParsingException("OPSIN bug: Unable to find functional group in oxide or addition compound rule");
			}
			String smiles = functionalGroups.get(0).getAttributeValue(VALUE_ATR);
			String elementStr = "";
			for (int i = 0; i < smiles.length(); i++) {
				if (Character.isUpperCase(smiles.charAt(i))){
					elementStr += smiles.charAt(i);
					if (i + 1 <smiles.length() && Character.isLowerCase(smiles.charAt(i + 1))){
						elementStr += smiles.charAt(i + 1);
					}
					break;
				}
			}
			return ChemEl.valueOf(elementStr);
		}
	}

}
