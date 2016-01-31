package uk.ac.cam.ch.wwmm.opsin;

import static uk.ac.cam.ch.wwmm.opsin.XmlDeclarations.*;
import static uk.ac.cam.ch.wwmm.opsin.StructureBuildingMethods.findRightMostGroupInSubBracketOrRoot;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.regex.Pattern;

class WordRulesOmittedSpaceCorrector {
	private static final Pattern matchAteOrIteEnding = Pattern.compile("[ai]t[e]?[\\])}]*$", Pattern.CASE_INSENSITIVE);
	
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
	 * @param divalentFunctionalGroupWordRule
	 */
	private void checkAndCorrectOmittedSpacesInDivalentFunctionalGroupRule(Element divalentFunctionalGroupWordRule)  {
		List<Element> substituentWords = OpsinTools.getChildElementsWithTagNameAndAttribute(divalentFunctionalGroupWordRule, WORD_EL, TYPE_ATR, SUBSTITUENT_TYPE_VAL);
		if (substituentWords.size() == 1){//potentially has been "wrongly" interpreted e.g. ethylmethyl ketone is more likely to mean ethyl methyl ketone
			List<Element> children = OpsinTools.getChildElementsWithTagNames(substituentWords.get(0), new String[]{SUBSTITUENT_EL, BRACKET_EL});
			if (children.size() == 2) {
				Element firstSubOrbracket = children.get(0);
				//rule out correct usage e.g. diethyl ether and locanted substituents e.g. 2-methylpropyl ether
				if (firstSubOrbracket.getAttribute(LOCANT_ATR) == null && firstSubOrbracket.getAttribute(MULTIPLIER_ATR) == null) {
					Element firstGroup = findRightMostGroupInSubBracketOrRoot(firstSubOrbracket);
					Fragment firstFrag = firstGroup.getFrag();
					if (hasSingleMonoValentCarbonOrSiliconRadical(firstFrag)) {
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
	}

	/**
	 * Corrects cases like methyl-2-ethylacetate --> methyl 2-ethylacetate
	 * @param wordRule
	 * @throws StructureBuildingException 
	 */
	private void checkAndCorrectOmittedSpaceEster(Element wordRule) throws StructureBuildingException {
		List<Element> words = wordRule.getChildElements(WORD_EL);
		if (words.size() != 1) {
			return;
		}
		Element word = words.get(0);
		String wordRuleContents = wordRule.getAttributeValue(VALUE_ATR);
		if (matchAteOrIteEnding.matcher(wordRuleContents).find()) {
			List<Element> children = OpsinTools.getChildElementsWithTagNames(word, new String[]{SUBSTITUENT_EL, BRACKET_EL, ROOT_EL});
			if (children.size() >= 2) {
				Element rootEl = children.get(children.size() - 1);
				Element rootGroup = findRightMostGroupInSubBracketOrRoot(rootEl);
				Fragment rootFrag = rootGroup.getFrag();
				int functionalAtomsCount = rootFrag.getFunctionalAtomCount();
				int rootMultiplier = 1;
				{
					String rootElMultiplierAtrVal = rootEl.getAttributeValue(MULTIPLIER_ATR);
					if (rootElMultiplierAtrVal != null) {
						rootMultiplier = Integer.parseInt(rootElMultiplierAtrVal);
						functionalAtomsCount *= rootMultiplier;
					}
				}
				if (functionalAtomsCount > 0){
					List<Element> substituents = children.subList(0, children.size() - 1);
					int substituentCount = substituents.size();
					if (substituentCount == 1 && rootMultiplier > 1) {
						return;
					}
					Element firstChild = substituents.get(0);
					if (!checkSuitabilityOfSubstituentForEsterFormation(firstChild, functionalAtomsCount)){
						return;
					}
					String multiplierValue = firstChild.getAttributeValue(MULTIPLIER_ATR);
					if (specialCaseWhereEsterPreferred(findRightMostGroupInSubBracketOrRoot(firstChild), multiplierValue, rootGroup, substituentCount)) {
						transformToEster(wordRule, firstChild);
					}
					else if (substituentCount > 1 && 
							(allBarFirstSubstituentHaveLocants(substituents) || insufficientSubstitutableHydrogenForSubstitution(substituents, rootFrag, rootMultiplier))){
						transformToEster(wordRule, firstChild);
					}
					else if ((substituentCount == 1 || rootMultiplier > 1) && substitutionWouldBeAmbiguous(rootFrag, multiplierValue)) {
						//either 1 substituent or multiplicative nomenclature (in the multiplicative nomenclature case many substituents will not have locants)
						transformToEster(wordRule, firstChild);
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

	private boolean insufficientSubstitutableHydrogenForSubstitution(List<Element> substituentsAndBrackets, Fragment frag, int rootMultiplier) {
		int substitutableHydrogens = getAtomForEachSubstitutableHydrogen(frag).size() * rootMultiplier;
		for (int i = 1; i < substituentsAndBrackets.size(); i++) {
			Element subOrBracket = substituentsAndBrackets.get(i);
			Fragment f = findRightMostGroupInSubBracketOrRoot(subOrBracket).getFrag();
			String multiplierValue = subOrBracket.getAttributeValue(MULTIPLIER_ATR);
			int multiplier = 1;
			if (multiplierValue != null){
				multiplier = Integer.parseInt(multiplierValue);
			}
			substitutableHydrogens -= (getTotalOutAtomValency(f) * multiplier);
		}
		Element potentialEsterSub = substituentsAndBrackets.get(0);
		int firstFragSubstitutableHydrogenRequired = getTotalOutAtomValency(findRightMostGroupInSubBracketOrRoot(potentialEsterSub).getFrag());
		String multiplierValue = potentialEsterSub.getAttributeValue(MULTIPLIER_ATR);
		int multiplier = 1;
		if (multiplierValue != null){
			multiplier = Integer.parseInt(multiplierValue);
		}
		if (substitutableHydrogens >=0 && (substitutableHydrogens - (firstFragSubstitutableHydrogenRequired * multiplier)) < 0){
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
	 * Ester form preferred when:
	 * mono is used on substituent
	 * alkyl chain is used on formate/acetate e.g. ethylacetate
	 * Root is carbamate, >=2 substituents, and this is the only word rule 
	 * (ester and non-ester carbamates differ only by whether or not there is a space, heuristically the ester is almost always intended under these conditions)
	 * @param substituentGroupEl
	 * @param multiplierValue
	 * @param rootGroup 
	 * @param numOfSubstituents 
	 * @return
	 */
	private boolean specialCaseWhereEsterPreferred(Element substituentGroupEl, String multiplierValue, Element rootGroup, int numOfSubstituents) {
		if (multiplierValue != null && Integer.parseInt(multiplierValue) == 1){
			return true;
		}
		String rootGroupName = rootGroup.getParent().getValue();
		if (substituentGroupEl.getAttributeValue(TYPE_ATR).equals(CHAIN_TYPE_VAL) &&
				ALKANESTEM_SUBTYPE_VAL.equals(substituentGroupEl.getAttributeValue(SUBTYPE_ATR))) {
			if (substituentGroupEl.getParent().getValue().matches(substituentGroupEl.getValue() + "yl-?") &&
					rootGroupName.matches(".*(form|methan|acet|ethan)[o]?ate?")) {
				return true;
			}
		}
		if ((rootGroupName.endsWith("carbamate") || rootGroupName.endsWith("carbamat")) && numOfSubstituents >= 2) {
			Element temp = substituentGroupEl.getParent();
			while (temp.getParent() != null) {
				temp = temp.getParent();
			}
			if (temp.getChildElements(WORDRULE_EL).size() == 1) {
				return true;
			}
		}
		return false;
	}

	private boolean substitutionWouldBeAmbiguous(Fragment frag, String multiplierValue) {
		int multiplier = 1;
		if (multiplierValue != null){
			multiplier = Integer.parseInt(multiplierValue);
		}
		if (multiplier == 1 && frag.getDefaultInAtom() != null) {
			return false;
		}
		List<Atom> atomForEachSubstitutableHydrogen = getAtomForEachSubstitutableHydrogen(frag);
		if (atomForEachSubstitutableHydrogen.size() == multiplier){
			return false;
		}
		StereoAnalyser analyser = new StereoAnalyser(frag);
		Set<String> uniqueEnvironments = new HashSet<String>();
		for (Atom a : atomForEachSubstitutableHydrogen) {
			uniqueEnvironments.add(AmbiguityChecker.getAtomEnviron(analyser, a));
		}
		if (uniqueEnvironments.size() == 1 && (multiplier == 1 || multiplier == atomForEachSubstitutableHydrogen.size() - 1)){
			return false;
		}
		return true;
	}

	private boolean checkSuitabilityOfSubstituentForEsterFormation(Element subOrBracket, int rootFunctionalAtomsCount) {
		if (subOrBracket.getAttribute(LOCANT_ATR) != null){
			return false;
		}
		Fragment rightMostGroup = findRightMostGroupInSubBracketOrRoot(subOrBracket).getFrag();
		if (!hasSingleMonoValentCarbonOrSiliconRadical(rightMostGroup)) {
			return false;
		}
		String multiplierStr = subOrBracket.getAttributeValue(MULTIPLIER_ATR);
		if (multiplierStr != null) {
			int multiplier = Integer.parseInt(multiplierStr);
			if (multiplier > rootFunctionalAtomsCount) {
				return false;
			}
		}
		return true;
	}
	
	private boolean hasSingleMonoValentCarbonOrSiliconRadical(Fragment frag) {
		if (frag.getOutAtomCount() == 1) {
			OutAtom outAtom = frag.getOutAtom(0);
			if (outAtom.getValency() == 1 && 
					(outAtom.getAtom().getElement() == ChemEl.C || outAtom.getAtom().getElement() == ChemEl.Si)) {
				return true;
			}
		}
		return false;
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
