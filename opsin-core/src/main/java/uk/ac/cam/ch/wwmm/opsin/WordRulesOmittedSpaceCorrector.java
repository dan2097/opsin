package uk.ac.cam.ch.wwmm.opsin;

import static uk.ac.cam.ch.wwmm.opsin.XmlDeclarations.*;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.regex.Pattern;

class WordRulesOmittedSpaceCorrector {
	private final static Pattern matchAteOrIteEnding = Pattern.compile("[ai]t[e]?$", Pattern.CASE_INSENSITIVE);
	
	private final BuildState state;
	private final Element parse;

	WordRulesOmittedSpaceCorrector(BuildState state, Element parse) {
		this.state =state;
		this.parse = parse;
	}

	void correctOmittedSpaces() throws StructureBuildingException {
		List<Element> wordRules = OpsinTools.getDescendantElementsWithTagName(parse, WORDRULE_EL);
		for (Element wordRule : wordRules) {
			WordRule wordRuleVal = WordRule.valueOf(wordRule.getAttributeValue(WORDRULE_ATR));
			if (wordRuleVal == WordRule.divalentFunctionalGroup){
				checkAndCorrectOmittedSpacesInDivalentFunctionalGroupRule(wordRule);
			}
			else if (wordRuleVal == WordRule.simple){
				//note that this function may change the word rule to ester
				checkAndCorrectOmittedSpaceEster(wordRule);
			}
		}
	}

	/**
	 * Corrects cases like "methylethyl ether" to "methyl ethyl ether"
	 * @param wordRule
	 */
	private void checkAndCorrectOmittedSpacesInDivalentFunctionalGroupRule(Element divalentFunctionalGroupWordRule)  {
		List<Element> substituentWords = OpsinTools.getChildElementsWithTagNameAndAttribute(divalentFunctionalGroupWordRule, WORD_EL, TYPE_ATR, SUBSTITUENT_TYPE_VAL);
		if (substituentWords.size() == 1){//potentially has been "wrongly" interpreted e.g. ethylmethyl ketone is more likely to mean ethyl methyl ketone
			List<Element> children = substituentWords.get(0).getChildElements();
			if (children.size() == 2){
				Element firstSubstituent = children.get(0);
				//rule out correct usage e.g. diethyl ether and locanted substituents e.g. 2-methylpropyl ether
				if (firstSubstituent.getAttribute(LOCANT_ATR)==null && firstSubstituent.getAttribute(MULTIPLIER_ATR)==null){
					Element subToMove =children.get(1);
					subToMove.detach();
					Element newWord =new GroupingEl(WORD_EL);
					newWord.addAttribute(new Attribute(TYPE_ATR, SUBSTITUENT_TYPE_VAL));
					newWord.addChild(subToMove);
					OpsinTools.insertAfter(substituentWords.get(0), newWord);
				}
			}
		}
	}
	
	
	/**
	 * Corrects cases like methyl-2-ethylacetate --> methyl 2-ethylacetate
	 * @param wordRule
	 * @throws StructureBuildingException 
	 */
	private void checkAndCorrectOmittedSpaceEster(Element wordRule) throws StructureBuildingException {
		List<Element> words = wordRule.getChildElements(WORD_EL);
		if (words.size()!=1){
			return;
		}
		Element word =words.get(0);
		String wordRuleContents = wordRule.getAttributeValue(VALUE_ATR);
		if (matchAteOrIteEnding.matcher(wordRuleContents).find()){
			List<Element> roots = word.getChildElements(ROOT_EL);
			if (roots.size()==1){
				Element rootGroup = roots.get(0).getFirstChildElement(GROUP_EL);
				if (AMINOACID_TYPE_VAL.equals(rootGroup.getAttributeValue(TYPE_ATR))){
					return;//amino acids are implicitly N locanted
				}
				Fragment rootFrag = state.xmlFragmentMap.get(rootGroup);
				int functionalAtomsCount = rootFrag.getFunctionalAtomCount();
				if (functionalAtomsCount >0){
					List<Element> substituentsAndBrackets = OpsinTools.getChildElementsWithTagNames(word, new String[]{SUBSTITUENT_EL, BRACKET_EL});
					if (substituentsAndBrackets.size()==0){
						return;
					}
					Element firstChild = substituentsAndBrackets.get(0);
					if (!checkSuitabilityOfSubstituentForEsterFormation(firstChild, functionalAtomsCount)){
						return;
					}
					if (substituentsAndBrackets.size()>1 && (allBarFirstSubstituentHaveLocants(substituentsAndBrackets) || insufficientSubstitutableHydrogenForSubstition(substituentsAndBrackets, rootFrag))){
						transformToEster(wordRule, firstChild);
					}
					else if (substituentsAndBrackets.size()==1){
						String multiplierValue = firstChild.getAttributeValue(MULTIPLIER_ATR);
						int multiplier = 1;
						if (multiplierValue!=null){
							multiplier= Integer.parseInt(multiplierValue);
						}
						if (specialCaseWhereEsterPreferred(getRightMostGroup(firstChild), multiplierValue, wordRuleContents) ||
								substitutionWouldBeAmbiguous(rootGroup, rootFrag, multiplier)){
							transformToEster(wordRule, firstChild);
						}
					}
				}
			}
		}
	}

	private boolean allBarFirstSubstituentHaveLocants(List<Element> substituentsAndBrackets) {
		if (substituentsAndBrackets.size() <=1){
			return false;
		}
		for (int i = 1; i < substituentsAndBrackets.size(); i++) {
			if (substituentsAndBrackets.get(i).getAttribute(LOCANT_ATR)==null){
				return false;
			}
		}
		return true;
	}

	private boolean insufficientSubstitutableHydrogenForSubstition(List<Element> substituentsAndBrackets, Fragment frag) {
		int substitutableHydrogens = getAtomForEachSubstitutableHydrogen(frag).size();
		for (int i = 1; i < substituentsAndBrackets.size(); i++) {
			Element subOrBracket = substituentsAndBrackets.get(i);
			Fragment f = state.xmlFragmentMap.get(getRightMostGroup(subOrBracket));
			String multiplierValue = subOrBracket.getAttributeValue(MULTIPLIER_ATR);
			int multiplier = 1;
			if (multiplierValue!=null){
				multiplier= Integer.parseInt(multiplierValue);
			}
			substitutableHydrogens -= (getTotalOutAtomValency(f) * multiplier);
		}
		int firstFragSubstitutableHydrogenRequired = getTotalOutAtomValency(state.xmlFragmentMap.get(getRightMostGroup(substituentsAndBrackets.get(0))));
		String multiplierValue = substituentsAndBrackets.get(0).getAttributeValue(MULTIPLIER_ATR);
		int multiplier = 1;
		if (multiplierValue!=null){
			multiplier= Integer.parseInt(multiplierValue);
		}
		if (substitutableHydrogens >=0 && (substitutableHydrogens - (firstFragSubstitutableHydrogenRequired * multiplier)) <0){
			return true;
		}
		return false;
	}

	private int getTotalOutAtomValency(Fragment f) {
		int outAtomValency = 0;
		for (int i = 0, l = f.getOutAtomCount(); i < l; i++) {
			outAtomValency += f.getOutAtom(i).getValency();
		}
		return outAtomValency;
	}

	/**
	 * Ester form preferred when mono is used and when an alkyl chain is used on formate/acetate
	 * e.g. ethylacetate
	 * @param substituentGroupEl
	 * @param multiplierValue
	 * @param wordRuleContents 
	 * @return
	 */
	private boolean specialCaseWhereEsterPreferred(Element substituentGroupEl, String multiplierValue, String wordRuleContents) {
		if (multiplierValue!=null && Integer.parseInt(multiplierValue)==1){
			return true;
		}
		if (substituentGroupEl.getAttributeValue(TYPE_ATR).equals(CHAIN_TYPE_VAL) && ALKANESTEM_SUBTYPE_VAL.equals(substituentGroupEl.getAttributeValue(SUBTYPE_ATR))){
			String potentialString = "(?i)" + substituentGroupEl.getValue() + "yl[\\-]?(form|methan|acet|ethan)[o]?ate";
			if (wordRuleContents.matches(potentialString)){
				return true;
			}
		}
		return false;
	}

	private boolean substitutionWouldBeAmbiguous(Element rootGroup, Fragment frag, int multiplier) {
		if (multiplier ==1 && (rootGroup.getAttribute(DEFAULTINID_ATR)!=null || rootGroup.getAttribute(DEFAULTINLOCANT_ATR)!=null)){
			return false;
		}
		List<Atom> atomForEachSubstitutableHydrogen = getAtomForEachSubstitutableHydrogen(frag);
		StereoAnalyser analyzer = new StereoAnalyser(frag);
		Set<Integer> uniqueEnvironments = new HashSet<Integer>();
		for (Atom a : atomForEachSubstitutableHydrogen) {
			uniqueEnvironments.add(analyzer.getAtomEnvironmentNumber(a));
		}
		if (atomForEachSubstitutableHydrogen.size()==multiplier){
			return false;
		}
		if (uniqueEnvironments.size()==1 && (multiplier==1 || multiplier == atomForEachSubstitutableHydrogen.size()-1)){
			return false;
		}
		return true;
	}

	private boolean checkSuitabilityOfSubstituentForEsterFormation(Element subOrBracket, int rootFunctionalAtomsCount) {
		if (subOrBracket.getAttribute(LOCANT_ATR)!=null){
			return false;
		}
		Fragment rightMostGroup = state.xmlFragmentMap.get(getRightMostGroup(subOrBracket));
		if (rightMostGroup.getOutAtomCount() != 1 || rightMostGroup.getOutAtom(0).getValency()!=1){
			return false;
		}
		String multiplierStr = subOrBracket.getAttributeValue(MULTIPLIER_ATR);
		if (multiplierStr!=null){
			int multiplier = Integer.parseInt(multiplierStr);
			if (multiplier > rootFunctionalAtomsCount){
				return false;
			}
		}
		return true;
	}

	/**
	 * Returns the right most group
	 * @param subOrBracket
	 * @return
	 */
	private Element getRightMostGroup (Element subOrBracket) {
		Element group;
		if (subOrBracket.getName().equals(BRACKET_EL)){
			group = StructureBuildingMethods.findRightMostGroupInBracket(subOrBracket);
		}
		else{
			group = subOrBracket.getFirstChildElement(GROUP_EL);
		}
		return group;
	}

	private List<Atom> getAtomForEachSubstitutableHydrogen(Fragment frag) {
		List<Atom> substitutableAtoms = new ArrayList<Atom>();
		List<Atom> atomList = frag.getAtomList();
		for (Atom atom : atomList) {
			if (FragmentTools.isCharacteristicAtom(atom)){
				continue;
			}
			int currentExpectedValency = atom.determineValency(true);
			int currentValency =  (atom.getIncomingValency() + (atom.hasSpareValency() ? 1 : 0) + atom.getOutValency());
			for (int i = currentValency; i < currentExpectedValency; i++) {
				substitutableAtoms.add(atom);
			}
		}
		return substitutableAtoms;
	}

	private void transformToEster(Element parentSimpleWordRule, Element substituentOrBracket) throws StructureBuildingException {
		parentSimpleWordRule.getAttribute(WORDRULE_ATR).setValue(WordRule.ester.toString());
		List<Element> childElsOfSub = substituentOrBracket.getChildElements();
		Element lastChildElOfSub =childElsOfSub.get(childElsOfSub.size()-1);
		if (lastChildElOfSub.getName().equals(HYPHEN_EL)){
			lastChildElOfSub.detach();
		}
		substituentOrBracket.detach();
		Element newSubstituentWord = new GroupingEl(WORD_EL);
		newSubstituentWord.addAttribute(new Attribute(TYPE_ATR, SUBSTITUENT_TYPE_VAL));
		newSubstituentWord.addChild(substituentOrBracket);
		parentSimpleWordRule.insertChild(newSubstituentWord, 0);
		String multiplierStr = substituentOrBracket.getAttributeValue(MULTIPLIER_ATR);
		if (multiplierStr!=null){
			substituentOrBracket.removeAttribute(substituentOrBracket.getAttribute(MULTIPLIER_ATR));
			int multiplier = Integer.parseInt(multiplierStr);
			for (int i = 1; i < multiplier; i++) {
				Element clone = state.fragManager.cloneElement(state, newSubstituentWord);
				OpsinTools.insertAfter(newSubstituentWord, clone);
			}
		}
	}
}
