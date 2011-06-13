package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;
import java.util.Set;
import java.util.Stack;

import uk.ac.cam.ch.wwmm.opsin.ParseWord.WordType;
import uk.ac.cam.ch.wwmm.opsin.WordRules.WordRule;

import nu.xom.Element;
import nu.xom.Elements;
import static uk.ac.cam.ch.wwmm.opsin.XmlDeclarations.*;
import static uk.ac.cam.ch.wwmm.opsin.OpsinTools.*;
import static uk.ac.cam.ch.wwmm.opsin.StructureBuildingMethods.*;

/**Constructs a single OPSIN fragment which describes the molecule from the ComponentGenerator/ComponentProcessor results.
 *
 * @author ptc24
 * @author dl387
 *
 */
class StructureBuilder {
	/**	Builds a molecule as a Fragment based on ComponentProcessor output.
	 * @param state
	 * @param molecule The ComponentProcessor output.
	 * @return A single Fragment - the built molecule.
	 * @throws StructureBuildingException If the molecule won't build - there may be many reasons.
	 */
	Fragment buildFragment(BuildState state, Element molecule) throws StructureBuildingException {
		Elements wordRules = molecule.getChildElements(WORDRULE_EL);
		if (wordRules.size()==0){
			throw new StructureBuildingException("Molecule contains no words!?");
		}
		Stack<Element> wordRuleStack = new Stack<Element>();
		for (int i = wordRules.size() -1; i >=0; i--) {
			wordRuleStack.add(wordRules.get(i));
		}
		
		List<Fragment> rGroups = new ArrayList<Fragment>();//rGroups need to represented as normal atoms for the purpose of working out stereochemistry. They will be converted to a suitable representation later
		List<Element> wordRulesVisited = new ArrayList<Element>();
		while (wordRuleStack.size()>0) {
			Element nextWordRuleEl = wordRuleStack.peek();//just has a look what's next
			if(!wordRulesVisited.contains(nextWordRuleEl)){
				wordRulesVisited.add(nextWordRuleEl);
				Elements wordRuleChildren = nextWordRuleEl.getChildElements(WORDRULE_EL);
				if (wordRuleChildren.size()!=0){//nested word rules
					for (int i = wordRuleChildren.size() -1; i >=0; i--) {
						wordRuleStack.add(wordRuleChildren.get(i));
					}
					continue;
				}
			}
			Element currentWordRuleEl = wordRuleStack.pop();
			WordRule wordRule = WordRule.valueOf(currentWordRuleEl.getAttributeValue(WORDRULE_ATR));
			List<Element> words = XOMTools.getChildElementsWithTagNames(currentWordRuleEl, new String[]{WORD_EL, WORDRULE_EL});
			state.currentWordRule =wordRule;
			if(wordRule == WordRule.simple || wordRule == WordRule.substituent) {
				for (Element word : words) {
					if (!word.getLocalName().equals(WORD_EL) || !word.getAttributeValue(TYPE_ATR).equals(WordType.full.toString()) && (!state.n2sConfig.isAllowRadicals() || !word.getAttributeValue(TYPE_ATR).equals(WordType.substituent.toString()))){
						throw new StructureBuildingException("OPSIN bug: Unexpected contents of 'simple' wordRule");
					}
					resolveWordOrBracket(state, word);
				}
			}
			else if(wordRule == WordRule.ester || wordRule == WordRule.multiEster) {
				buildEster(state, words);//e.g. ethyl ethanoate, dimethyl terephthalate,  methyl propanamide
			}
			else if (wordRule == WordRule.divalentFunctionalGroup){
				buildDiValentFunctionalGroup(state, words);// diethyl ether or methyl propyl ketone
			}
			else if (wordRule == WordRule.monovalentFunctionalGroup){
				buildMonovalentFunctionalGroup(state, words);// ethyl chloride, isophthaloyl dichloride, diethyl ether, ethyl alcohol
			}
			else if(wordRule == WordRule.functionalClassEster) {
				buildFunctionalClassEster(state, words);//e.g. ethanoic acid ethyl ester, tetrathioterephthalic acid dimethyl ester
			}
			else if (wordRule == WordRule.amide){
				//e.g. ethanoic acid ethyl amide, terephthalic acid dimethyl amide, ethanoic acid amide
				//already processed by the ComponentProcessor
				for (Element word : words) {
					resolveWordOrBracket(state, word);
				}
			}
			else if(wordRule == WordRule.oxide) {
				buildOxide(state, words);//e.g. styrene oxide, triphenylphosphane oxide, thianthrene 5,5-dioxide, propan-2-one oxide
			}
			else if(wordRule == WordRule.carbonylDerivative) {
				buildCarbonylDerivative(state, words);//e.g. Imidazole-2-carboxamide O-ethyloxime, pentan-3-one oxime
			}
			else if(wordRule == WordRule.anhydride) {//e.g. acetic anhydride
				buildAnhydride(state, words);
			}
			else if(wordRule == WordRule.acidHalideOrPseudoHalide) {//e.g. phosphinimidic chloride
				buildAcidHalideOrPseudoHalide(state, words);
			}
			else if(wordRule == WordRule.additionCompound) {//e.g. carbon tetrachloride
				buildAdditionCompound(state, words);
			}
			else if (wordRule == WordRule.glycol){
				buildGlycol(state, words);//e.g. ethylene glycol
			}
			else if (wordRule == WordRule.glycolEther){
				buildGlycolEther(state, words);//e.g. octaethyleneglycol monododecyl ether
			}
			else if(wordRule == WordRule.acetal) {
				buildAcetal(state, words);//e.g. propanal diethyl acetal
			}
			else if(wordRule == WordRule.biochemicalEster) {
				buildBiochemicalEster(state, words, wordRules.size());//e.g. uridine 5'-(tetrahydrogen triphosphate)
			}
			else if(wordRule == WordRule.polymer) {
				rGroups.addAll(buildPolymer(state, words));
			}
			else if(wordRule == WordRule.hydrazide) {
				//e.g. carbonic dihydrazide
				//already processed by the ComponentProcessor
				for (Element word : words) {
					resolveWordOrBracket(state, word);
				}
			}
			else{
				throw new StructureBuildingException("Unknown Word Rule");
			}
		}
		
		List<Element> groupElements = XOMTools.getDescendantElementsWithTagName(molecule, GROUP_EL);
		processOxidoSpecialCase(state, groupElements);
		processOxidationNumbers(state, groupElements);
		state.fragManager.convertSpareValenciesToDoubleBonds();
		state.fragManager.checkValencies();
		
		boolean explicitStoichometryPresent = applyExplicitStoichometryIfProvided(state, wordRules);
		int overallCharge = state.fragManager.getOverallCharge();
		if (overallCharge!=0 && wordRules.size() >1){//a net charge is present! Could just mean the counterion has not been specified though
			balanceChargeIfPossible(state, molecule, overallCharge, explicitStoichometryPresent);
		}
		makeHydrogensExplicit(state);

		Fragment uniFrag = state.fragManager.getUnifiedFragment();
		processStereochemistry(state, molecule, uniFrag);

		for (Fragment rGroup : rGroups) {
			Atom rAtom = rGroup.getFirstAtom();
			rAtom.setElement("R");
		}

		if (uniFrag.getInAtoms().size()>0){
			throw new StructureBuildingException("In atoms as used in multiplicative nomenclature were never used, maybe the name is malformed");
		}
		if (uniFrag.getOutAtoms().size()>0 && !state.n2sConfig.isAllowRadicals()){
			throw new StructureBuildingException("Radicals are currently set to not convert to structures");
		}
		return uniFrag;
	}

	private void buildEster(BuildState state, List<Element> words) throws StructureBuildingException {
		int wordIndice = 0;
		Element currentWord=words.get(wordIndice);
		BuildResults substituentsBr = new BuildResults();
		while (currentWord.getAttributeValue(TYPE_ATR).equals(WordType.substituent.toString())){
			resolveWordOrBracket(state, currentWord);
			BuildResults substituentBr = new BuildResults(state, currentWord);
			int outAtomCount = substituentBr.getOutAtomCount();
			boolean traditionalEster =false;
			for (int i = 0; i < outAtomCount; i++) {
				OutAtom out = substituentBr.getOutAtom(i);
				if (out.getValency()>1){
					FragmentTools.splitOutAtomIntoValency1OutAtoms(out);
					traditionalEster =true;
				}
			}
			if (traditionalEster){//e.g. ethylidene dipropanoate
				substituentBr = new BuildResults(state, currentWord);
				outAtomCount = substituentBr.getOutAtomCount();
			}
			if (outAtomCount ==1){//TODO add support for locanted terepthaloyl
				String locantForSubstituent = currentWord.getAttributeValue(LOCANT_ATR);
				if (locantForSubstituent!=null){
					substituentBr.getFirstOutAtom().setLocant(locantForSubstituent);//indexes which functional atom to connect to when there is a choice. Also can disambiguate which atom is a S in things like thioates
				}
			}
			else if (outAtomCount ==0){
				throw new StructureBuildingException("Substituent was expected to have at least one outAtom");
			}
			substituentsBr.mergeBuildResults(substituentBr);
			currentWord=words.get(++wordIndice);
		}

		if (wordIndice == words.size() || !words.get(wordIndice).getAttributeValue(TYPE_ATR).equals(WordType.full.toString())){
			throw new StructureBuildingException("Full word not found where full word expected: missing ate group in ester");
		}

		List<BuildResults> ateGroups = new ArrayList<BuildResults>();
		Map<BuildResults, String> buildResultsToLocant = new HashMap<BuildResults, String>();//typically locant will be null
		for (; wordIndice < words.size(); wordIndice++) {
			Element word =words.get(wordIndice);
			if (word.getAttributeValue(TYPE_ATR).equals(WordType.full.toString())){
				String locant = word.getAttributeValue(LOCANT_ATR);//specifying a locant for an ateWord is very unusual as this information is typically redundant c.f. dodecamethylene 1,12-bis(chloroformate)
				resolveWordOrBracket(state, word);
				BuildResults ateBR = new BuildResults(state, word);
				if (ateBR.getFunctionalAtomCount()<1){
					throw new StructureBuildingException("bug? ate group did not have any functional atoms!");
				}
				ateGroups.add(ateBR);
				buildResultsToLocant.put(ateBR, locant);
			}
			else{
				throw new StructureBuildingException("Non full word found where only full words were expected");
			}
		}
		int outAtomCount =substituentsBr.getOutAtomCount();
		int esterIdCount = 0;
		for (BuildResults br : ateGroups) {
			esterIdCount += br.getFunctionalAtomCount();
		}
		if (outAtomCount > esterIdCount){
			throw new StructureBuildingException("There are more radicals in the substituents(" + outAtomCount +") than there are places to form esters("+esterIdCount+")");
		}
		for(int i=0; i< outAtomCount; i++) {
			BuildResults ateBr = ateGroups.get(i % ateGroups.size());
			Atom ateAtom;
			if (substituentsBr.getFirstOutAtom().getLocant()!=null){
				ateAtom =determineFunctionalAtomToUse(substituentsBr.getFirstOutAtom().getLocant(), ateBr);
			}
			else{
				ateAtom =ateBr.getFunctionalAtom(0);
				ateBr.removeFunctionalAtom(0);
			}
			String locant = buildResultsToLocant.get(ateBr);
			if (locant ==null){//typical case
				Atom atomOnSubstituentToUse =substituentsBr.getOutAtomTakingIntoAccountWhetherSetExplicitly(0);
				state.fragManager.createBond(ateAtom, atomOnSubstituentToUse, 1);
				substituentsBr.removeOutAtom(0);
			}
			else{
				Integer outAtomPosition =null;
				for (int j = 0; j < substituentsBr.getOutAtomCount(); j++) {
					if (substituentsBr.getOutAtom(j).getAtom().hasLocant(locant)){
						outAtomPosition = j;
						break;
					}
				}
				if (outAtomPosition ==null){
					throw new StructureBuildingException("Unable to find substituent with locant: " + locant + " to form ester!");
				}
				Atom atomOnSubstituentToUse = substituentsBr.getOutAtom(outAtomPosition).getAtom();
				state.fragManager.createBond(ateAtom, atomOnSubstituentToUse, 1);
				substituentsBr.removeOutAtom(outAtomPosition);
			}
			ateAtom.neutraliseCharge();
		}
	}



	private void buildDiValentFunctionalGroup(BuildState state, List<Element> words) throws StructureBuildingException {
		int wordIndice = 0;
		if (!words.get(wordIndice).getAttributeValue(TYPE_ATR).equals(WordType.substituent.toString())) {
			throw new StructureBuildingException("word: " +wordIndice +" was expected to be a substituent");
		}
		resolveWordOrBracket(state, words.get(wordIndice));
		BuildResults substituent1 =new BuildResults(state, words.get(wordIndice));
		if (substituent1.getOutAtom(0).getValency() !=1){
			throw new StructureBuildingException("OutAtom has unexpected valency. Expected 1. Actual: " + substituent1.getOutAtom(0).getValency());
		}
		BuildResults substituent2;
		if (substituent1.getOutAtomCount()==2){// e.g. tetramethylene sulfone
			if (substituent1.getOutAtom(1).getValency() !=1){
				throw new StructureBuildingException("OutAtom has unexpected valency. Expected 1. Actual: " + substituent1.getOutAtom(1).getValency());
			}
			substituent2 = substituent1;
		}
		else{
			if (substituent1.getOutAtomCount()!=1){
				throw new StructureBuildingException("Expected one outAtom. Found " + substituent1.getOutAtomCount() );
			}
			wordIndice++;
			if (words.get(wordIndice).getAttributeValue(TYPE_ATR).equals(WordType.functionalTerm.toString())) {//e.g. methyl sulfoxide rather than dimethyl sulfoxide
				Element clone = state.fragManager.cloneElement(state, words.get(0));
				XOMTools.insertAfter(words.get(0), clone);
				words = OpsinTools.elementsToElementArrayList(((Element)words.get(0).getParent()).getChildElements());
			}
			else{
				resolveWordOrBracket(state, words.get(wordIndice));
			}
			substituent2 =new BuildResults(state, words.get(wordIndice));
			if (substituent2.getOutAtomCount()!=1){
				throw new StructureBuildingException("Expected one outAtom. Found " + substituent2.getOutAtomCount() );
			}
			if (substituent2.getOutAtom(0).getValency() !=1){
				throw new StructureBuildingException("OutAtom has unexpected valency. Expected 1. Actual: " + substituent2.getOutAtom(0).getValency());
			}
		}
		wordIndice++;
		if (words.get(wordIndice) ==null || !words.get(wordIndice).getAttributeValue(TYPE_ATR).equals(WordType.functionalTerm.toString())) {
			throw new StructureBuildingException(words.get(wordIndice).getValue()+" was expected to be a functionalTerm");
		}
		List<Element> functionalGroup = XOMTools.getDescendantElementsWithTagName(words.get(wordIndice), FUNCTIONALGROUP_EL);
		if (functionalGroup.size()!=1){
			throw new StructureBuildingException("Unexpected number of functionalGroups found, could be a bug in OPSIN's grammar");
		}
		String smilesOfGroup = functionalGroup.get(0).getAttributeValue(VALUE_ATR);
		Fragment diValentGroup =state.fragManager.buildSMILES(smilesOfGroup, FUNCTIONALCLASS_TYPE_VAL, NONE_LABELS_VAL);

		Atom outAtom1 =substituent1.getOutAtomTakingIntoAccountWhetherSetExplicitly(0);
		substituent1.removeOutAtom(0);
		Atom outAtom2 = substituent2.getOutAtomTakingIntoAccountWhetherSetExplicitly(0);
		substituent2.removeOutAtom(0);
		if (diValentGroup.getOutAtoms().size()==1){//c.f. peroxide where it is a linker
			state.fragManager.createBond(outAtom1, diValentGroup.getOutAtom(0).getAtom(), 1);
			diValentGroup.removeOutAtom(0);
			state.fragManager.createBond(outAtom2, diValentGroup.getFirstAtom(), 1);
		}
		else{
			if (outAtom1 != outAtom2){//general case
				state.fragManager.createBond(outAtom1, diValentGroup.getFirstAtom(), 1);
				state.fragManager.createBond(outAtom2, diValentGroup.getFirstAtom(), 1);
			}
			else{//e.g. carbonyl sulfide
				state.fragManager.createBond(outAtom1, diValentGroup.getFirstAtom(), 2);
			}
		}
		state.fragManager.incorporateFragment(diValentGroup, outAtom1.getFrag());
	}

	private void buildMonovalentFunctionalGroup(BuildState state, List<Element> words) throws StructureBuildingException {
		resolveWordOrBracket(state, words.get(0));
		List<Element> groups = XOMTools.getDescendantElementsWithTagName(words.get(0), GROUP_EL);
		for (Element group : groups) {//replaces outAtoms with valency greater than 1 with multiple outAtoms; e.g. ylidene -->diyl
			Fragment frag = state.xmlFragmentMap.get(group);
			for (int i = frag.getOutAtoms().size()-1; i>=0; i--) {
				OutAtom outAtom =frag.getOutAtom(i);
				if (outAtom.getValency()>1){
					FragmentTools.splitOutAtomIntoValency1OutAtoms(outAtom);
				}
			}
		}
		BuildResults substituentBR = new BuildResults(state, words.get(0));

		List<Fragment> functionalGroupFragments = new ArrayList<Fragment>();
		for (int i=1; i<words.size(); i++ ) {
			Element functionalGroupWord =words.get(i);
			List<Element> functionalGroups = XOMTools.getDescendantElementsWithTagName(functionalGroupWord, FUNCTIONALGROUP_EL);
			if (functionalGroups.size()!=1){
				throw new StructureBuildingException("Expected exactly 1 functionalGroup. Found " + functionalGroups.size());
			}
			
			Fragment monoValentFunctionGroup =state.fragManager.buildSMILES(functionalGroups.get(0).getAttributeValue(VALUE_ATR), FUNCTIONALCLASS_TYPE_VAL, NONE_LABELS_VAL);
			if (functionalGroups.get(0).getAttributeValue(TYPE_ATR).equals(MONOVALENTSTANDALONEGROUP_TYPE_VAL)){
				Atom ideAtom = monoValentFunctionGroup.getDefaultInAtom();
				ideAtom.addChargeAndProtons(1, 1);//e.g. make cyanide charge netural
			}
			Element possibleMultiplier = (Element) XOMTools.getPreviousSibling(functionalGroups.get(0));
			functionalGroupFragments.add(monoValentFunctionGroup);
			if (possibleMultiplier!=null){
				int multiplierValue = Integer.parseInt(possibleMultiplier.getAttributeValue(VALUE_ATR));
				for (int j = 1; j < multiplierValue; j++) {
					functionalGroupFragments.add(state.fragManager.copyFragment(monoValentFunctionGroup));
				}
				possibleMultiplier.detach();
			}
		}
		
		int outAtomCount =substituentBR.getOutAtomCount();
		if (outAtomCount > functionalGroupFragments.size()){//something like isophthaloyl chloride (more precisely written isophthaloyl dichloride)
			if (functionalGroupFragments.size()!=1){
				throw new StructureBuildingException("Incorrect number of functional groups found to balance outAtoms");
			}
			Fragment monoValentFunctionGroup = functionalGroupFragments.get(0);
			for (int j = 1; j < outAtomCount; j++) {
				functionalGroupFragments.add(state.fragManager.copyFragment(monoValentFunctionGroup));
			}
		}
		else if (functionalGroupFragments.size() > outAtomCount){
			throw new StructureBuildingException("There are more function groups to attach than there are positions to attach them to!");
		}
		for (int i = 0; i < outAtomCount; i++) {
			Fragment ideFrag =functionalGroupFragments.get(i);
			Atom ideAtom = ideFrag.getDefaultInAtom();
			Atom subAtom=substituentBR.getOutAtomTakingIntoAccountWhetherSetExplicitly(0);
			state.fragManager.createBond(ideAtom, subAtom, 1);
			substituentBR.removeOutAtom(0);
			state.fragManager.incorporateFragment(ideFrag, subAtom.getFrag());
		}
	}

	private void buildFunctionalClassEster(BuildState state, List<Element> words) throws StructureBuildingException {
		if (!words.get(0).getAttributeValue(TYPE_ATR).equals(WordType.full.toString())){
			throw new StructureBuildingException("Don't alter wordRules.xml without checking the consequences!");
		}
		resolveWordOrBracket(state, words.get(0));//the group
		BuildResults acidBr = new BuildResults(state, words.get(0));

		if (acidBr.getFunctionalAtomCount()==0){
			throw new StructureBuildingException("No functionalAtoms detected!");
		}

		int i=1;
		Element currentWord = words.get(i);
		while (currentWord.getAttributeValue(TYPE_ATR).equals(WordType.substituent.toString())){
			if (acidBr.getFunctionalAtomCount()==0){
				throw new StructureBuildingException("Insufficient functionalAtoms on acid");
			}
			resolveWordOrBracket(state, currentWord);
			BuildResults substituentBr = new BuildResults(state, currentWord);
			if (substituentBr.getOutAtomCount() ==1){
				String locantForSubstituent = currentWord.getAttributeValue(LOCANT_ATR);
				Atom functionalAtom;
				if (locantForSubstituent!=null){
					functionalAtom =determineFunctionalAtomToUse(locantForSubstituent, acidBr);
				}
				else{
					functionalAtom =acidBr.getFunctionalAtom(0);
					acidBr.removeFunctionalAtom(0);
				}
				if (substituentBr.getOutAtom(0).getValency()!=1){
					throw new StructureBuildingException("Substituent was expected to have only have an outgoing valency of 1");
				}
				state.fragManager.createBond(functionalAtom,substituentBr.getOutAtomTakingIntoAccountWhetherSetExplicitly(0), 1);
				if (functionalAtom.getCharge()==-1){
					functionalAtom.neutraliseCharge();
				}
				substituentBr.removeOutAtom(0);
			}
			else {
				throw new StructureBuildingException("Substituent was expected to have one outAtom");
			}
			currentWord=words.get(++i);
		}
		if (!words.get(i++).getAttributeValue(TYPE_ATR).equals(WordType.functionalTerm.toString())){//ester
			throw new StructureBuildingException("Number of words different to expectations; did not find ester");
		}
	}
	
	/**
	 * Handles names like thiophene 1,1-dioxide; carbon dioxide; benzene oxide
	 * Does the same for sulfide/selenide/telluride
	 * @param state
	 * @param words
	 * @throws StructureBuildingException
	 */
	private void buildOxide(BuildState state, List<Element> words) throws StructureBuildingException {
		resolveWordOrBracket(state, words.get(0));//the group
		List<Fragment> oxideFragments = new ArrayList<Fragment>();
		List<String> locantsForOxide =new ArrayList<String>();//often not specified
		if (!words.get(1).getAttributeValue(TYPE_ATR).equals(WordType.functionalTerm.toString())){
			throw new StructureBuildingException("Oxide functional term not found where expected!");
		}
		Element rightMostGroup;
		if (words.get(0).getLocalName().equals(WORDRULE_EL)){//e.g. Nicotinic acid N-oxide
			List<Element> fullWords = XOMTools.getDescendantElementsWithTagNameAndAttribute(words.get(0), WORD_EL, TYPE_ATR, WordType.full.toString());
			if (fullWords.size()==0){
				throw new StructureBuildingException("OPSIN is entirely unsure where the oxide goes so has decided not to guess");
			}
			rightMostGroup = findRightMostGroupInBracket(fullWords.get(fullWords.size()-1));
		}
		else{
			rightMostGroup = findRightMostGroupInBracket(words.get(0));
		}
		
		int numberOfOxygenToAdd =1;
		List<Element> multipliers =XOMTools.getDescendantElementsWithTagName(words.get(1), MULTIPLIER_EL);
		if (multipliers.size() >1){
			throw new StructureBuildingException("Expected 0 or 1 multiplier found: " + multipliers.size());
		}
		if (multipliers.size()==1){
			numberOfOxygenToAdd = Integer.parseInt(multipliers.get(0).getAttributeValue(VALUE_ATR));
			multipliers.get(0).detach();
		}
		else{
			if (ELEMENTARYATOM_SUBTYPE_VAL.equals(rightMostGroup.getAttributeValue(SUBTYPE_ATR))){
				Atom elementaryAtom = state.xmlFragmentMap.get(rightMostGroup).getFirstAtom();
				int charge = elementaryAtom.getCharge();
				if (charge >0 && charge %2 ==0){
					numberOfOxygenToAdd = charge/2;
				}
				else if (elementaryAtom.getProperty(Atom.OXIDATION_NUMBER)!=null){
					int valency = elementaryAtom.getProperty(Atom.OXIDATION_NUMBER) - elementaryAtom.getIncomingValency();
					if (valency >0 && valency %2 ==0){
						numberOfOxygenToAdd = valency/2;
					}
				}
			}
		}
		List<Element> functionalGroup =XOMTools.getDescendantElementsWithTagName(words.get(1), FUNCTIONALGROUP_EL);
		if (functionalGroup.size()!=1){
			throw new StructureBuildingException("Expected 1 group element found: " + functionalGroup.size());
		}
		String smilesReplacement = functionalGroup.get(0).getAttributeValue(VALUE_ATR);
		String labels =  functionalGroup.get(0).getAttributeValue(LABELS_ATR);
		for (int i = 0; i < numberOfOxygenToAdd; i++) {
			oxideFragments.add(state.fragManager.buildSMILES(smilesReplacement, FUNCTIONALCLASS_TYPE_VAL, labels));
		}
		List<Element> locantEls =XOMTools.getDescendantElementsWithTagName(words.get(1), LOCANT_EL);
		if (locantEls.size() >1){
			throw new StructureBuildingException("Expected 0 or 1 locant elements found: " + locantEls.size());
		}
		if (locantEls.size()==1){
			String[] locants = MATCH_COMMA.split(StringTools.removeDashIfPresent(locantEls.get(0).getValue()));
            locantsForOxide.addAll(Arrays.asList(locants));
			locantEls.get(0).detach();
		}
		if (!locantsForOxide.isEmpty() && locantsForOxide.size()!=oxideFragments.size()){
			throw new StructureBuildingException("Mismatch between number of locants and number of oxides specified");
		}
		List<Fragment> orderedPossibleFragments = new ArrayList<Fragment>();//In preference suffixes are substituted onto e.g. acetonitrile oxide
		Elements suffixEls = ((Element)rightMostGroup.getParent()).getChildElements(SUFFIX_EL);
		for (int i = suffixEls.size()-1; i >=0; i--) {//suffixes (if any) from right to left
			Element suffixEl = suffixEls.get(i);
			Fragment suffixFrag =state.xmlFragmentMap.get(suffixEl);
			if (suffixFrag!=null){
				orderedPossibleFragments.add(suffixFrag);
			}
		}
		Fragment groupToModify = state.xmlFragmentMap.get(rightMostGroup);//all the suffixes are actually part of this fragment already
		orderedPossibleFragments.add(groupToModify);
		mainLoop: for (int i = 0; i < oxideFragments.size(); i++) {
			Atom oxideAtom = oxideFragments.get(i).getFirstAtom();
			if (!locantsForOxide.isEmpty()){
				Atom atomToAddOxideTo =groupToModify.getAtomByLocantOrThrow(locantsForOxide.get(i));
				formAppropriateBondToOxideAndAdjustCharges(state, atomToAddOxideTo, oxideAtom);
			}
			else{
				for (Fragment frag : orderedPossibleFragments) {
					String subTypeVal = state.xmlFragmentMap.getElement(frag).getAttributeValue(SUBTYPE_ATR);
					if (ELEMENTARYATOM_SUBTYPE_VAL.equals(subTypeVal)){
						Atom elementaryAtom= frag.getFirstAtom();
						formAppropriateBondToOxideAndAdjustCharges(state, elementaryAtom, oxideAtom);//e.g. carbon dioxide
						int chargeOnAtom =elementaryAtom.getCharge();
						if (chargeOnAtom>=2){
							elementaryAtom.setCharge(chargeOnAtom-2);
						}
						continue mainLoop;
					}
					else{
						List<Atom> atomList = frag.getAtomList();
						for (Atom atom : atomList) {
							if (!atom.getElement().equals("C") && !atom.getElement().equals("O")){
								formAppropriateBondToOxideAndAdjustCharges(state, atom, oxideAtom);
								continue mainLoop;
							}
						}
					}
				}
				//No heteroatoms could be found. Perhaps it's supposed to be something like styrene oxide
				Set<Bond> bondSet = groupToModify.getBondSet();//looking for double bond
				for (Bond bond : bondSet) {
					if (bond.getOrder()==2 && bond.getFromAtom().getElement().equals("C") && bond.getToAtom().getElement().equals("C")){
						bond.setOrder(1);
						state.fragManager.createBond(bond.getFromAtom(), oxideAtom, 1);
						state.fragManager.createBond(bond.getToAtom(), oxideAtom, 1);
						continue mainLoop;
					}
				}
				
				//...or maybe something a bit iffy nomenclature wise like benzene oxide :-S
				for (Bond bond : bondSet) {
					Atom fromAtom =bond.getFromAtom();
					Atom toAtom = bond.getToAtom();
					if (fromAtom.hasSpareValency() && toAtom.hasSpareValency() &&fromAtom.getElement().equals("C") && toAtom.getElement().equals("C")){
						fromAtom.setSpareValency(false);
						toAtom.setSpareValency(false);
						state.fragManager.createBond(fromAtom, oxideAtom, 1);
						state.fragManager.createBond(toAtom, oxideAtom, 1);
						continue mainLoop;
					}
				}
				for (Fragment frag : orderedPossibleFragments) {//something like where oxide goes on an oxygen propan-2-one oxide
					List<Atom> atomList = frag.getAtomList();
					for (Atom atom : atomList) {
						if (!atom.getElement().equals("C")){
							formAppropriateBondToOxideAndAdjustCharges(state, atom, oxideAtom);
							continue mainLoop;
						}
					}
				}
				throw new StructureBuildingException("Unable to find suitable atom or a double bond to add oxide to");
			}
		}
		for (Fragment oxide : oxideFragments) {
			state.fragManager.incorporateFragment(oxide, groupToModify);
		}
	}

	/**
	 * Decides whether an oxide should double bond e.g. P=O or single bond as a zwitterionic form e.g. [N+]-[O-]
	 * Corrects the charges if necessary and forms the bond
	 * @param state
	 * @param atomToAddOxideTo
	 * @param oxideAtom
	 * @throws StructureBuildingException 
	 */
	private void formAppropriateBondToOxideAndAdjustCharges(BuildState state, Atom atomToAddOxideTo, Atom oxideAtom) throws StructureBuildingException {
		Integer maxVal = ValencyChecker.getMaximumValency(atomToAddOxideTo.getElement(), atomToAddOxideTo.getCharge());
		if (maxVal ==null || (atomToAddOxideTo.getIncomingValency() + atomToAddOxideTo.getOutValency() +2) <= maxVal){
			if (atomToAddOxideTo.getLambdaConventionValency()==null || !ValencyChecker.checkValencyAvailableForBond(atomToAddOxideTo, 2)){//probably in well formed names 2 protons should always be added but some names use the lambdaConvention to specify the valency after oxide has been applied
				atomToAddOxideTo.addChargeAndProtons(0, 2);//this is an additive operation, up the proton count by 2
			}
			state.fragManager.createBond(atomToAddOxideTo, oxideAtom, 2);
		}
		else{
			if (atomToAddOxideTo.getCharge()!=0 || oxideAtom.getCharge()!=0){
				throw new StructureBuildingException("Oxide appeared to refer to an atom that has insufficent valency to accept the addition of oxygen");
			}
			atomToAddOxideTo.addChargeAndProtons(1, 1);
			oxideAtom.addChargeAndProtons(-1, -1);
			maxVal = ValencyChecker.getMaximumValency(atomToAddOxideTo.getElement(), atomToAddOxideTo.getCharge());
			if (maxVal !=null && (atomToAddOxideTo.getIncomingValency() + atomToAddOxideTo.getOutValency() +1) > maxVal){
				throw new StructureBuildingException("Oxide appeared to refer to an atom that has insufficent valency to accept the addition of oxygen");
			}
			state.fragManager.createBond(atomToAddOxideTo, oxideAtom, 1);
		}
	}

	private void buildCarbonylDerivative(BuildState state, List<Element> words) throws StructureBuildingException {
		if (!WordType.full.toString().equals(words.get(0).getAttributeValue(TYPE_ATR))){
			throw new StructureBuildingException("OPSIN bug: Wrong word type encountered when applying carbonylDerivative wordRule");
		}
		List<Fragment> replacementFragments = new ArrayList<Fragment>();
		List<String> locantForFunctionalTerm =new ArrayList<String>();//usually not specified
		if (!words.get(1).getAttributeValue(TYPE_ATR).equals(WordType.functionalTerm.toString())){//e.g. acetone O-ethyloxime or acetone 1-chloro-1-methylhydrazone
			for (int i = 1; i < words.size(); i++) {
				Fragment frag = state.xmlFragmentMap.get(findRightMostGroupInWordOrWordRule(words.get(i)));
				replacementFragments.add(frag);
				Elements children =words.get(i).getChildElements();
				if (children.size()==1 && children.get(0).getLocalName().equals(BRACKET_EL) && children.get(0).getAttribute(LOCANT_ATR)!=null){
					locantForFunctionalTerm.add(children.get(0).getAttributeValue(LOCANT_ATR));
				}
				else if (children.size()==2 && children.get(0).getAttribute(LOCANT_ATR)!=null ){
					String locant =children.get(0).getAttributeValue(LOCANT_ATR);
					if (children.get(1).getLocalName().equals(ROOT_EL) && !frag.hasLocant(locant) && MATCH_NUMERIC_LOCANT.matcher(locant).matches()){ //e.g. 1,3-benzothiazole-2-carbaldehyde 2-phenylhydrazone
						locantForFunctionalTerm.add(children.get(0).getAttributeValue(LOCANT_ATR));
						children.get(0).removeAttribute(children.get(0).getAttribute(LOCANT_ATR));
					}
				}
			}
		}
		else{//e.g. butan-2,3-dione dioxime or hexan2,3-dione 2-oxime
			int numberOfCarbonylReplacements =1;
			List<Element> multipliers =XOMTools.getDescendantElementsWithTagName(words.get(1), MULTIPLIER_EL);
			if (multipliers.size() >1){
				throw new StructureBuildingException("Expected 0 or 1 multiplier found: " + multipliers.size());
			}
			if (multipliers.size()==1){
				numberOfCarbonylReplacements = Integer.parseInt(multipliers.get(0).getAttributeValue(VALUE_ATR));
				multipliers.get(0).detach();
			}
			List<Element> functionalGroup =XOMTools.getDescendantElementsWithTagName(words.get(1), FUNCTIONALGROUP_EL);
			if (functionalGroup.size()!=1){
				throw new StructureBuildingException("Expected 1 functionalGroup element found: " + functionalGroup.size());
			}
			String smilesReplacement = functionalGroup.get(0).getAttributeValue(VALUE_ATR);
			String labels =  functionalGroup.get(0).getAttributeValue(LABELS_ATR);
			for (int i = 0; i < numberOfCarbonylReplacements; i++) {
				Fragment replacementFragment = state.fragManager.buildSMILES(smilesReplacement, FUNCTIONALCLASS_TYPE_VAL, labels);
				if (i >0){
					FragmentTools.relabelLocants(replacementFragment.getAtomList(), StringTools.multiplyString("'", i));
				}
				List<Atom> atomList = replacementFragment.getAtomList();
				for (Atom atom : atomList) {
					atom.removeLocantsOtherThanElementSymbolLocants();//prevents numeric locant locanted substitution from outside the functional word
				}
				replacementFragments.add(replacementFragment);
			}
			List<Element> locantEls =XOMTools.getDescendantElementsWithTagName(words.get(1), LOCANT_EL);
			if (locantEls.size() >1){
				throw new StructureBuildingException("Expected 0 or 1 locant elements found: " + locantEls.size());
			}
			if (locantEls.size()==1){
				String[] locants = MATCH_COMMA.split(StringTools.removeDashIfPresent(locantEls.get(0).getValue()));
                locantForFunctionalTerm.addAll(Arrays.asList(locants));
				locantEls.get(0).detach();
			}
		}
		if (!locantForFunctionalTerm.isEmpty() && locantForFunctionalTerm.size()!=replacementFragments.size()){
			throw new StructureBuildingException("Mismatch between number of locants and number of carbonyl replacements");
		}

		Element rightMostGroup = findRightMostGroupInWordOrWordRule(words.get(0));
		Element parent = (Element) rightMostGroup.getParent();
		boolean multiplied =false;
		while (!parent.equals(words.get(0))){
			if (parent.getAttribute(MULTIPLIER_ATR)!=null){
				multiplied =true;
			}
			parent =(Element) parent.getParent();
		}
		if (!multiplied){
			List<Atom> carbonylOxygens = findCarbonylOxygens(state.xmlFragmentMap.get(rightMostGroup), locantForFunctionalTerm);
			int replacementsToPerform = Math.min(replacementFragments.size(), carbonylOxygens.size());
			replaceCarbonylOxygenWithReplacementFragments(state, words, replacementFragments, carbonylOxygens, replacementsToPerform);
		}

		resolveWordOrBracket(state, words.get(0));//the component
		if (replacementFragments.size() >0){
			//Note that the right most group may be multiplied e.g. 3,3'-methylenebis(2,4,6-trimethylbenzaldehyde) disemicarbazone
			//or the carbonyl may not even be on the right most group e.g.  4-oxocyclohexa-2,5-diene-1-carboxylic acid 4-oxime
			BuildResults br = new BuildResults(state, words.get(0));
			List<Atom> carbonylOxygens = new ArrayList<Atom>();
			List<Fragment> fragments = new ArrayList<Fragment>(br.getFragments());
			for (ListIterator<Fragment> iterator = fragments.listIterator(fragments.size()); iterator.hasPrevious();) {//iterate in reverse order - right most groups preferred
				carbonylOxygens.addAll(findCarbonylOxygens(iterator.previous(), locantForFunctionalTerm));
			}
			replaceCarbonylOxygenWithReplacementFragments(state, words, replacementFragments, carbonylOxygens, replacementFragments.size());
		}
	}

	private void replaceCarbonylOxygenWithReplacementFragments(BuildState state, List<Element> words, List<Fragment> replacementFragments, List<Atom> carbonylOxygens, int functionalReplacementsToPerform) throws StructureBuildingException {
		if (functionalReplacementsToPerform > carbonylOxygens.size()){
			throw new StructureBuildingException("Insufficient carbonyl groups found!");
		}
		for (int i = 0; i < functionalReplacementsToPerform; i++) {
			Atom carbonylOxygen =carbonylOxygens.remove(0);//the oxygen of the carbonyl
			Fragment carbonylFrag = carbonylOxygen.getFrag();
			Fragment replacementFrag = replacementFragments.remove(0);
			List<Atom> atomList = replacementFrag.getAtomList();
			Atom atomToReplaceCarbonylOxygen = atomList.get(atomList.size()-1);
			Atom numericLocantAtomConnectedToCarbonyl = OpsinTools.depthFirstSearchForAtomWithNumericLocant(carbonylOxygen);
			if (numericLocantAtomConnectedToCarbonyl!=null){
				atomList.get(0).addLocant(atomList.get(0).getElement() + numericLocantAtomConnectedToCarbonyl.getFirstLocant());//adds a locant like O1 giving another way of referencing this atom
			}
			if (!words.get(1).getAttributeValue(TYPE_ATR).equals(WordType.functionalTerm.toString())){
				resolveWordOrBracket(state, words.get(1 +i));
			}
			for (Atom atom : atomList) {
				atom.removeLocantsOtherThanElementSymbolLocants();//prevents numeric locant locanted substitution from outside the functional word
				List<String> locants =atom.getLocants();
				for (int j = locants.size() -1; j >=0; j--) {
					String locant = locants.get(j);
					if (carbonylFrag.hasLocant(locant)){
						atom.removeLocant(locant);
					}
				}
			}
			state.fragManager.replaceAtomWithAnotherAtomPreservingConnectivity(carbonylOxygen, atomToReplaceCarbonylOxygen);
			atomToReplaceCarbonylOxygen.setType(carbonylOxygen.getType());//copy the type e.g. if the carbonyl was a suffix this should appear as a suffix
			if (state.xmlFragmentMap.getElement(replacementFrag)==null){//incorporate only for the case that replacementFrag came from a functional class element
				state.fragManager.incorporateFragment(replacementFrag, carbonylFrag);
			}
		}
	}

	/**
	 * Given a fragment and optionally a list of locants finds carbonyl atoms
	 * If locants are given the carbonyl must be assoicated with one of the given locants
	 * @param fragment
	 * @param locantForCarbonylAtom
	 * @return
	 * @throws StructureBuildingException
	 */
	private List<Atom> findCarbonylOxygens(Fragment fragment, List<String> locantForCarbonylAtom) throws StructureBuildingException {
		List<Atom> matches = new ArrayList<Atom>();
		List<Atom> rootFragAtomList = fragment.getAtomList();
		for (Atom atom : rootFragAtomList) {//find all carbonyl oxygen
			if (atom.getElement().equals("O") && atom.getCharge()==0){
				List<Atom> neighbours =atom.getAtomNeighbours();
				if (neighbours.size()==1){
					if (neighbours.get(0).getElement().equals("C")){
						if (!locantForCarbonylAtom.isEmpty()){
							Atom numericLocantAtomConnectedToCarbonyl = OpsinTools.depthFirstSearchForAtomWithNumericLocant(atom);
							if (numericLocantAtomConnectedToCarbonyl!=null){//could be the carbon of the carbonyl or the ring the carbonyl connects to in say a carbaldehyde
								boolean matchesLocant = false;
								for (String locant : locantForCarbonylAtom) {
									if (numericLocantAtomConnectedToCarbonyl.hasLocant(locant)){
										matchesLocant =true;
									}
								}
								if (!matchesLocant){
									continue;
								}
							}
							else{
								continue;
							}
						}
						Bond b = atom.getBondToAtomOrThrow(neighbours.get(0));
						if (b.getOrder()==2){
							matches.add(atom);
						}
					}
				}
			}
		}
		return matches;
	}

	private void buildAnhydride(BuildState state, List<Element> words) throws StructureBuildingException {
		if (words.size()!=2 && words.size()!=3){
			throw new StructureBuildingException("Unexpected number of words in anhydride. Check wordRules.xml, this is probably a bug");
		}
		Element anhydrideWord = words.get(words.size()-1);
		List<Element> functionalClass =XOMTools.getDescendantElementsWithTagName(anhydrideWord, FUNCTIONALGROUP_EL);
		if (functionalClass.size()!=1){
			throw new StructureBuildingException("Expected 1 group element found: " + functionalClass.size());
		}
		String anhydrideSmiles = functionalClass.get(0).getAttributeValue(VALUE_ATR);
		int numberOfAnhydrideLinkages =1;
		List<Element> multipliers =XOMTools.getDescendantElementsWithTagName(anhydrideWord, MULTIPLIER_EL);
		if (multipliers.size() >1){
			throw new StructureBuildingException("Expected 0 or 1 multiplier found: " + multipliers.size());
		}
		if (multipliers.size()==1){
			numberOfAnhydrideLinkages = Integer.parseInt(multipliers.get(0).getAttributeValue(VALUE_ATR));
			multipliers.get(0).detach();
		}
		String anhydrideLocant = null;
		List<Element> anhydrideLocants =XOMTools.getDescendantElementsWithTagNames(anhydrideWord, new String[]{LOCANT_EL, COLONSEPERATEDLOCANT_EL});
		if (anhydrideLocants.size() >1){
			throw new StructureBuildingException("Expected 0 or 1 anhydrideLocants found: " + anhydrideLocants.size());
		}
		if (anhydrideLocants.size()==1){
			anhydrideLocant = anhydrideLocants.get(0).getValue();
			anhydrideLocants.get(0).detach();
		}
		resolveWordOrBracket(state, words.get(0));
		BuildResults br1 = new BuildResults(state, words.get(0));
		if (br1.getFunctionalAtomCount() ==0){
			throw new StructureBuildingException("Cannot find functionalAtom to form anhydride");
		}
		if (words.size()==3){//asymmetric anhydride
			if (anhydrideLocant!=null){
				throw new StructureBuildingException("Unsupported or invalid anhydride");
			}
			resolveWordOrBracket(state, words.get(1));
			BuildResults br2 = new BuildResults(state, words.get(1));
			if (br2.getFunctionalAtomCount() ==0){
				throw new StructureBuildingException("Cannot find functionalAtom to form anhydride");
			}
			if (numberOfAnhydrideLinkages>1){
				for (int i = numberOfAnhydrideLinkages-1; i >=0 ; i--) {
					if (br2.getFunctionalAtomCount()==0){
						throw new StructureBuildingException("Cannot find functionalAtom to form anhydride");
					}
					BuildResults newAcidBr;
					if (i!=0){
						Element newAcid = state.fragManager.cloneElement(state, words.get(0));
						XOMTools.insertAfter(words.get(0), newAcid);
						newAcidBr = new BuildResults(state, newAcid);
					}
					else{
						newAcidBr =br1;
					}
					formAnhydrideLink(state, anhydrideSmiles, newAcidBr, br2);
				}
				
			}
			else{
				if (br1.getFunctionalAtomCount()!=1 && br2.getFunctionalAtomCount()!=1 ) {
					throw new StructureBuildingException("Invalid anhydride description");
				}
				formAnhydrideLink(state, anhydrideSmiles, br1, br2);
			}
		}
		else{//symmetric anhydride
			if (br1.getFunctionalAtomCount()>1){//cyclic anhydride
				if (br1.getFunctionalAtomCount()==2){
					if (numberOfAnhydrideLinkages!=1 || anhydrideLocant !=null ){
						throw new StructureBuildingException("Unsupported or invalid anhydride");
					}
					formAnhydrideLink(state, anhydrideSmiles, br1, br1);
				}
				else{//cyclic anhydride where group has more than 2 acids
					if (anhydrideLocant ==null){
						throw new StructureBuildingException("Anhydride formation appears to be ambiguous; More than 2 acids, no locants");
					}
					String[] acidLocants =MATCH_COLON.split(StringTools.removeDashIfPresent(anhydrideLocant));
					if (acidLocants.length != numberOfAnhydrideLinkages){
						throw new StructureBuildingException("Mismatch between number of locants and number of anhydride linkages to form");
					}
					if (br1.getFunctionalAtomCount() < (numberOfAnhydrideLinkages *2)){
						throw new StructureBuildingException("Mismatch between number of acid atoms and number of anhydride linkages to form");
					}
					List<Atom> functionalAtoms = new ArrayList<Atom>();
					for (int i = 0; i < br1.getFunctionalAtomCount(); i++) {
						functionalAtoms.add(br1.getFunctionalAtom(i));
					}
			
					for (int i = 0; i < numberOfAnhydrideLinkages; i++) {
						String[] locants = MATCH_COMMA.split(acidLocants[i]);
						Atom oxygen1 =null;
						for (int j = functionalAtoms.size() -1; j >=0; j--) {
							Atom functionalAtom = functionalAtoms.get(j);
							Atom numericLocantAtomConnectedToFunctionalAtom = OpsinTools.depthFirstSearchForAtomWithNumericLocant(functionalAtom);
							if (numericLocantAtomConnectedToFunctionalAtom.hasLocant(locants[0])){
								oxygen1=functionalAtom;
								functionalAtoms.remove(j);
								break;
							}
						}
						Atom oxygen2 =null;
						for (int j = functionalAtoms.size() -1; j >=0; j--) {
							Atom functionalAtom = functionalAtoms.get(j);
							Atom numericLocantAtomConnectedToFunctionalAtom = OpsinTools.depthFirstSearchForAtomWithNumericLocant(functionalAtom);
							if (numericLocantAtomConnectedToFunctionalAtom.hasLocant(locants[1])){
								oxygen2=functionalAtom;
								functionalAtoms.remove(j);
								break;
							}
						}
						if (oxygen1 ==null || oxygen2==null){
							throw new StructureBuildingException("Unable to find locanted atom for anhydride formation");
						}
						formAnhydrideLink(state, anhydrideSmiles, oxygen1, oxygen2);
					}
				}
			}
			else{
				if (numberOfAnhydrideLinkages!=1 || anhydrideLocant !=null ){
					throw new StructureBuildingException("Unsupported or invalid anhydride");
				}
				Element newAcid = state.fragManager.cloneElement(state, words.get(0));
				XOMTools.insertAfter(words.get(0), newAcid);
				BuildResults br2 = new BuildResults(state, newAcid);
				formAnhydrideLink(state, anhydrideSmiles, br1, br2);
			}
		}
	}

	/**
	 * Given buildResults for both the acids and the SMILES of the anhydride forms the anhydride bond using the first functionalAtom on each BuildResults
	 * @param state
	 * @param anhydrideSmiles
	 * @param acidBr1
	 * @param acidBr2
	 * @throws StructureBuildingException
	 */
	private void formAnhydrideLink(BuildState state, String anhydrideSmiles, BuildResults acidBr1, BuildResults acidBr2)throws StructureBuildingException {
		Atom oxygen1 = acidBr1.getFunctionalAtom(0);
		acidBr1.removeFunctionalAtom(0);
		Atom oxygen2 = acidBr2.getFunctionalAtom(0);
		acidBr2.removeFunctionalAtom(0);
		formAnhydrideLink(state, anhydrideSmiles, oxygen1, oxygen2);
	}
	
	/**
	 * Given two atoms and the SMILES of the anhydride forms the anhydride bond
	 * @param state
	 * @param anhydrideSmiles
	 * @param oxygen1
	 * @param oxygen2
	 * @throws StructureBuildingException
	 */
	private void formAnhydrideLink(BuildState state, String anhydrideSmiles, Atom oxygen1, Atom oxygen2)throws StructureBuildingException {
		if (!oxygen1.getElement().equals("O")||!oxygen2.getElement().equals("O") || oxygen1.getBonds().size()!=1 ||oxygen2.getBonds().size()!=1) {
			throw new StructureBuildingException("Problem building anhydride");
		}
		Atom atomOnSecondAcidToConnectTo = oxygen2.getAtomNeighbours().get(0);
		state.fragManager.removeAtomAndAssociatedBonds(oxygen2);
		Fragment anhydride = state.fragManager.buildSMILES(anhydrideSmiles, FUNCTIONALCLASS_TYPE_VAL, NONE_LABELS_VAL);
		Fragment acidFragment1 = oxygen1.getFrag();
		state.fragManager.replaceAtomWithAnotherAtomPreservingConnectivity(oxygen1, anhydride.getFirstAtom());
		List<Atom> atomsInAnhydrideLinkage = anhydride.getAtomList();
		state.fragManager.createBond(atomsInAnhydrideLinkage.get(atomsInAnhydrideLinkage.size()-1), atomOnSecondAcidToConnectTo, 1);
		state.fragManager.incorporateFragment(anhydride, acidFragment1);
	}
	
	private void buildAcidHalideOrPseudoHalide(BuildState state, List<Element> words) throws StructureBuildingException {
		if (!words.get(0).getAttributeValue(TYPE_ATR).equals(WordType.full.toString())){
			throw new StructureBuildingException("Don't alter wordRules.xml without checking the consequences!");
		}
		resolveWordOrBracket(state, words.get(0));
		BuildResults acidBr = new BuildResults(state, words.get(0));
		int functionalAtomCount =acidBr.getFunctionalAtomCount();
		if (functionalAtomCount==0){
			throw new StructureBuildingException("No functionalAtoms detected!");
		}

		boolean monoMultiplierDetected =false;
		List<Fragment> functionalGroupFragments = new ArrayList<Fragment>();
		for (int i=1; i<words.size(); i++ ) {
			Element functionalGroupWord =words.get(i);
			List<Element> functionalGroups = XOMTools.getDescendantElementsWithTagName(functionalGroupWord, FUNCTIONALGROUP_EL);
			if (functionalGroups.size()!=1){
				throw new StructureBuildingException("Expected exactly 1 functionalGroup. Found " + functionalGroups.size());
			}
			
			Fragment monoValentFunctionGroup =state.fragManager.buildSMILES(functionalGroups.get(0).getAttributeValue(VALUE_ATR), FUNCTIONALCLASS_TYPE_VAL, NONE_LABELS_VAL);
			if (functionalGroups.get(0).getAttributeValue(TYPE_ATR).equals(MONOVALENTSTANDALONEGROUP_TYPE_VAL)){
				Atom ideAtom = monoValentFunctionGroup.getDefaultInAtom();
				ideAtom.addChargeAndProtons(1, 1);//e.g. make cyanide charge netural
			}
			Element possibleMultiplier = (Element) XOMTools.getPreviousSibling(functionalGroups.get(0));
			functionalGroupFragments.add(monoValentFunctionGroup);
			if (possibleMultiplier!=null){
				int multiplierValue = Integer.parseInt(possibleMultiplier.getAttributeValue(VALUE_ATR));
				if (multiplierValue==1){
					monoMultiplierDetected = true;
				}
				for (int j = 1; j < multiplierValue; j++) {
					functionalGroupFragments.add(state.fragManager.copyFragment(monoValentFunctionGroup));
				}
				possibleMultiplier.detach();
			}
		}
		int halideCount = functionalGroupFragments.size();
		if (halideCount > functionalAtomCount || (!monoMultiplierDetected && halideCount <functionalAtomCount)){
			throw new StructureBuildingException("Mismatch between number of halide/pseudo halide fragments and acidic oxygens");
		}
		for (int i = halideCount - 1; i>=0; i--) {
			Fragment ideFrag =functionalGroupFragments.get(i);
			Atom ideAtom = ideFrag.getDefaultInAtom();
			Atom acidAtom = acidBr.getFunctionalAtom(i);
			if (!acidAtom.getElement().equals("O")){
				throw new StructureBuildingException("Atom type expected to be oxygen but was: " +acidAtom.getElement());
			}
			acidBr.removeFunctionalAtom(i);
			Fragment acidFragment =acidAtom.getFrag();
			state.fragManager.replaceAtomWithAnotherAtomPreservingConnectivity(acidAtom, ideAtom);
			state.fragManager.incorporateFragment(ideFrag, acidFragment);
		}
	}
	
	private void buildAdditionCompound(BuildState state, List<Element> words) throws StructureBuildingException {
		if (!words.get(0).getAttributeValue(TYPE_ATR).equals(WordType.full.toString())){
			throw new StructureBuildingException("Don't alter wordRules.xml without checking the consequences!");
		}
		resolveWordOrBracket(state, words.get(0));
		Element elementaryAtomEl = StructureBuildingMethods.findRightMostGroupInBracket(words.get(0));
		Fragment elementaryAtomFrag = state.xmlFragmentMap.get(elementaryAtomEl);
		Atom elementaryAtom = elementaryAtomFrag.getFirstAtom();
		int charge = elementaryAtom.getCharge();
		List<Fragment> functionalGroupFragments = new ArrayList<Fragment>();
		for (int i=1; i<words.size(); i++ ) {
			Element functionalGroupWord =words.get(i);
			List<Element> functionalGroups = XOMTools.getDescendantElementsWithTagName(functionalGroupWord, FUNCTIONALGROUP_EL);
			if (functionalGroups.size()!=1){
				throw new StructureBuildingException("Expected exactly 1 functionalGroup. Found " + functionalGroups.size());
			}
			
			Fragment monoValentFunctionGroup =state.fragManager.buildSMILES(functionalGroups.get(0).getAttributeValue(VALUE_ATR), FUNCTIONALCLASS_TYPE_VAL, NONE_LABELS_VAL);
			if (functionalGroups.get(0).getAttributeValue(TYPE_ATR).equals(MONOVALENTSTANDALONEGROUP_TYPE_VAL)){
				Atom ideAtom = monoValentFunctionGroup.getDefaultInAtom();
				ideAtom.addChargeAndProtons(1, 1);//e.g. make cyanide charge netural
			}
			Element possibleMultiplier = (Element) XOMTools.getPreviousSibling(functionalGroups.get(0));
			functionalGroupFragments.add(monoValentFunctionGroup);
			if (possibleMultiplier!=null){
				int multiplierValue = Integer.parseInt(possibleMultiplier.getAttributeValue(VALUE_ATR));
				for (int j = 1; j < multiplierValue; j++) {
					functionalGroupFragments.add(state.fragManager.copyFragment(monoValentFunctionGroup));
				}
				possibleMultiplier.detach();
			}
			else if (words.size()==2){//silicon chloride -->silicon tetrachloride
				int incomingBondOrder =elementaryAtom.getIncomingValency();
				int expectedValency;
				if (charge > 0){
					expectedValency = incomingBondOrder + charge;
				}
				else{
					if (elementaryAtom.getProperty(Atom.OXIDATION_NUMBER)!=null){
						expectedValency = elementaryAtom.getProperty(Atom.OXIDATION_NUMBER);
					}
					else{
						if (elementaryAtomEl.getAttribute(COMMONOXIDATIONSTATESANDMAX_ATR)!=null){
							String[] typicalOxidationStates = MATCH_COMMA.split(MATCH_COLON.split(elementaryAtomEl.getAttributeValue(COMMONOXIDATIONSTATESANDMAX_ATR))[0]);
							expectedValency = Integer.parseInt(typicalOxidationStates[0]);
						}
						else{
							expectedValency = ValencyChecker.getPossibleValencies(elementaryAtom.getElement(), charge)[0];
						}
					}
				}
				int implicitMultiplier = expectedValency -incomingBondOrder >1 ? expectedValency -incomingBondOrder : 1;
				for (int j = 1; j < implicitMultiplier; j++) {
					functionalGroupFragments.add(state.fragManager.copyFragment(monoValentFunctionGroup));
				}
			}
		}
		int halideCount = functionalGroupFragments.size();
		if (charge>0){
			elementaryAtom.setCharge(charge - halideCount);
		}
		Integer maximumVal = ValencyChecker.getMaximumValency(elementaryAtom.getElement(), elementaryAtom.getCharge());
		if (maximumVal!=null && halideCount > maximumVal){
			throw new StructureBuildingException("Too many halides/psuedo halides addded to " +elementaryAtom.getElement());
		}
		for (int i = halideCount - 1; i>=0; i--) {
			Fragment ideFrag =functionalGroupFragments.get(i);
			Atom ideAtom = ideFrag.getDefaultInAtom();
			state.fragManager.incorporateFragment(ideFrag, ideAtom, elementaryAtomFrag, elementaryAtom, 1);
		}
	}
	

	private void buildGlycol(BuildState state, List<Element> words) throws StructureBuildingException {
		int wordIndice  = 0;
		resolveWordOrBracket(state, words.get(wordIndice));//the group
		Element finalGroup = findRightMostGroupInWordOrWordRule(words.get(wordIndice));
		Fragment theDiRadical = state.xmlFragmentMap.get(finalGroup);
		if (theDiRadical.getOutAtoms().size()!=2){
			throw new StructureBuildingException("Glycol class names (e.g. ethylene glycol) expect two outAtoms. Found: " + theDiRadical.getOutAtoms() );
		}
		wordIndice++;
		if (wordIndice >= words.size() || !words.get(wordIndice).getAttributeValue(TYPE_ATR).equals(WordType.functionalTerm.toString())){
			throw new StructureBuildingException("Glycol functionalTerm word expected");
		}
		List<Element> functionalClassEls = XOMTools.getDescendantElementsWithTagName(words.get(wordIndice), FUNCTIONALCLASS_EL);
		if (functionalClassEls.size()!=1){
			throw new StructureBuildingException("Glycol functional class not found where expected");
		}
		
		Atom outAtom1 = theDiRadical.getAtomOrNextSuitableAtomOrThrow(theDiRadical.getOutAtom(0).getAtom(), theDiRadical.getOutAtom(0).getValency(), false);
		Fragment functionalFrag =state.fragManager.buildSMILES(functionalClassEls.get(0).getAttributeValue(VALUE_ATR), FUNCTIONALCLASS_TYPE_VAL, NONE_LABELS_VAL);
		if (theDiRadical.getOutAtom(0).getValency() !=1){
			throw new StructureBuildingException("OutAtom has unexpected valency. Expected 1. Actual: " + theDiRadical.getOutAtom(0).getValency());
		}
		state.fragManager.createBond(outAtom1, functionalFrag.getFirstAtom(), 1);
		state.fragManager.incorporateFragment(functionalFrag, theDiRadical);
		
		Atom outAtom2 = theDiRadical.getAtomOrNextSuitableAtomOrThrow(theDiRadical.getOutAtom(1).getAtom(), theDiRadical.getOutAtom(1).getValency(), false);
		Fragment hydroxy =state.fragManager.buildSMILES("O", FUNCTIONALCLASS_TYPE_VAL, NONE_LABELS_VAL);
		if (theDiRadical.getOutAtom(1).getValency() !=1){
			throw new StructureBuildingException("OutAtom has unexpected valency. Expected 1. Actual: " + theDiRadical.getOutAtom(1).getValency());
		}
		state.fragManager.createBond(outAtom2, hydroxy.getFirstAtom(), 1);
		state.fragManager.incorporateFragment(hydroxy, theDiRadical);
		theDiRadical.removeOutAtom(1);
		theDiRadical.removeOutAtom(0);
	}
	

	/**
	 * Handles Glcyol ethers nomenclature e.g.
	 * triethylene glycol n-butyl ether
	 * tripropylene glycol methyl ether
	 * dipropylene glycol methyl ether acetate
	 * @param state
	 * @param words
	 * @throws StructureBuildingException
	 */
	private void buildGlycolEther(BuildState state, List<Element> words) throws StructureBuildingException {
		List<Element> wordsToAttachToGlcyol = new ArrayList<Element>();
		Element glycol =words.get(0);
		if (!glycol.getAttributeValue(TYPE_ATR).equals(WordType.full.toString())){
			throw new StructureBuildingException("OPSIN Bug: Cannot find glycol word!");
		}
		for (int i = 1; i < words.size(); i++) {
			Element wordOrWordRule =words.get(i);
			//ether ignored
			if (!wordOrWordRule.getAttributeValue(TYPE_ATR).equals(WordType.functionalTerm.toString())){
				wordsToAttachToGlcyol.add(wordOrWordRule);
			}
			else if (!wordOrWordRule.getAttributeValue(VALUE_ATR).equals("ether")){
				throw new StructureBuildingException("Unexpected word encountered when applying glycol ether word rule " + wordOrWordRule.getAttributeValue(VALUE_ATR));
			}
		}
		if (wordsToAttachToGlcyol.size() !=1 && wordsToAttachToGlcyol.size() !=2 ){
			throw new StructureBuildingException("Unexpected number of substituents for glycol ether. Expected 1 or 2 found: " +wordsToAttachToGlcyol.size());
		}
		Element finalGroup = findRightMostGroupInWordOrWordRule(glycol);
		Fragment theDiRadical = state.xmlFragmentMap.get(finalGroup);
		List<Atom> atomList = theDiRadical.getAtomList();
		List<Atom> glycolAtoms = new ArrayList<Atom>();
		for (Atom atom : atomList) {
			if (atom.getElement().equals("O")&& atom.getType().equals(FUNCTIONALCLASS_TYPE_VAL)){
				glycolAtoms.add(atom);
			}
		}
		if (glycolAtoms.size()!=2){
			throw new StructureBuildingException("OPSIN bug: unable to find the two glycol oxygens");
		}
		BuildResults br1 = new BuildResults(state, wordsToAttachToGlcyol.get(0));
		if (br1.getOutAtomCount() ==0){
			throw new StructureBuildingException("Substituent had no outAtom to form glycol ether");
		}
		state.fragManager.createBond(glycolAtoms.get(0), br1.getOutAtom(0).getAtom(), 1);
		br1.removeOutAtom(0);
		if (wordsToAttachToGlcyol.size()==2){
			BuildResults br2 = new BuildResults(state, wordsToAttachToGlcyol.get(1));
			if (br2.getOutAtomCount() >0){//form ether
				state.fragManager.createBond(glycolAtoms.get(1), br2.getOutAtom(0).getAtom(), 1);
				br2.removeOutAtom(0);
			}
			else if (br2.getFunctionalAtomCount() >0){//form ester
				Atom ateAtom = br2.getFunctionalAtom(0);
				ateAtom.neutraliseCharge();
				state.fragManager.replaceAtomWithAnotherAtomPreservingConnectivity(glycolAtoms.get(1), br2.getFunctionalAtom(0));
				br2.removeFunctionalAtom(0);
			}
			else{
				throw new StructureBuildingException("Word had neither an outAtom or a functionalAtom! hence neither and ether or ester could be formed : " + wordsToAttachToGlcyol.get(1).getAttributeValue(VALUE_ATR));
			}
		}
	}

	/**
	 * Builds acetals/ketals/hemiacetals/hemiketals and chalcogen analogues
	 * The distinction between acetals and ketals is not enforced (ketals are a subset of acetals)
	 * @param state
	 * @param words
	 * @throws StructureBuildingException
	 */
	private void buildAcetal(BuildState state, List<Element> words) throws StructureBuildingException {
		for (int i = 0; i < words.size()-1; i++) {
			resolveWordOrBracket(state, words.get(i));
		}
		BuildResults substituentsBr = new BuildResults();
		for (int i = 1; i < words.size()-1; i++) {
			Element currentWord = words.get(i);
			BuildResults substituentBr = new BuildResults(state, currentWord);
			int outAtomCount = substituentBr.getOutAtomCount();
			if (outAtomCount ==1){
				String locantForSubstituent = currentWord.getAttributeValue(LOCANT_ATR);
				if (locantForSubstituent!=null){
					substituentBr.getFirstOutAtom().setLocant(locantForSubstituent);
				}
			}
			else if (outAtomCount ==0){
				throw new StructureBuildingException("Substituent was expected to have at least one outAtom");
			}
			substituentsBr.mergeBuildResults(substituentBr);
		}	
		Element rightMostGroup = findRightMostGroupInWordOrWordRule(words.get(0));
		Fragment rootFragment = state.xmlFragmentMap.get(rightMostGroup);//the group which will be modified
		List<Atom> carbonylOxygen= findCarbonylOxygens(rootFragment, new ArrayList<String>());
		Element functionalWord = words.get(words.size()-1);
		List<Element> functionalClasses = XOMTools.getDescendantElementsWithTagName(functionalWord, FUNCTIONALCLASS_EL);
		if (functionalClasses.size()!=1){
			throw new StructureBuildingException("OPSIN bug: unable to find acetal functionalClass");
		}
		Element functionalClassEl = functionalClasses.get(0);
		String functionalClass = functionalClassEl.getValue();
		Element beforeAcetal = (Element) XOMTools.getPreviousSibling(functionalClassEl);
		int numberOfAcetals =1;
		List<String> elements = null;
		if (beforeAcetal!=null){
			if (beforeAcetal.getLocalName().equals(MULTIPLIER_EL)){
				numberOfAcetals = Integer.parseInt(beforeAcetal.getAttributeValue(VALUE_ATR));
			}
			else{
				elements = determineChalcogenReplacementOfAcetal(functionalClassEl);
				if (elements.size()>2){
					throw new StructureBuildingException(functionalClass + " only has two oxygen");
				}
				if (elements.size()==1){
					elements.add("O");
				}
			}
		}
		if (elements==null){
			elements = new ArrayList<String>();
			elements.add("O");
			elements.add("O");
		}

		if (carbonylOxygen.size() < numberOfAcetals){
			throw new StructureBuildingException("Insufficient carbonyls to form " + numberOfAcetals +" " + functionalClass );
		}
		boolean hemiacetal = functionalClass.contains("hemi");
		List<Fragment> acetalFrags = new ArrayList<Fragment>();
		for (int i = 0; i < numberOfAcetals; i++) {
			acetalFrags.add(formAcetal(state, carbonylOxygen, elements));
		}
		int bondsToForm = hemiacetal ? numberOfAcetals : 2*numberOfAcetals;
		if (substituentsBr.getOutAtomCount()!=bondsToForm){
			throw new StructureBuildingException("incorrect number of susbtituents when forming " + functionalClass);
		}
		connectSubstituentsToAcetal(state, acetalFrags, substituentsBr, hemiacetal);
	}

	private List<String> determineChalcogenReplacementOfAcetal(Element functionalClassEl) throws StructureBuildingException {
		Element currentEl = (Element) functionalClassEl.getParent().getChild(0);
		int multiplier =1;
		List<String> elements = new ArrayList<String>();
		while(currentEl !=functionalClassEl){
			if (currentEl.getLocalName().equals(MULTIPLIER_EL)){
				multiplier = Integer.parseInt(currentEl.getAttributeValue(VALUE_ATR));
			}
			else if (currentEl.getLocalName().equals(GROUP_EL)){
				for (int i = 0; i < multiplier; i++) {
					elements.add(currentEl.getAttributeValue(VALUE_ATR));
				}
			}
			else{
				throw new StructureBuildingException("Unexpected element before acetal");
			}
			currentEl =(Element) XOMTools.getNextSibling(currentEl);
		}
		return elements;
	}

	private Fragment formAcetal(BuildState state, List<Atom> carbonylOxygen, List<String> elements) throws StructureBuildingException {
		Atom neighbouringCarbon = carbonylOxygen.get(0).getAtomNeighbours().get(0);
		state.fragManager.removeAtomAndAssociatedBonds(carbonylOxygen.get(0));
		carbonylOxygen.remove(0);
		Fragment acetalFrag = state.fragManager.buildSMILES(StringTools.stringListToString(elements, "."));
		FragmentTools.assignElementLocants(acetalFrag, new ArrayList<Fragment>());
		List<Atom> acetalAtomList = acetalFrag.getAtomList();
		Atom atom1 = acetalAtomList.get(0);
		state.fragManager.createBond(neighbouringCarbon, atom1, 1);
		Atom atom2 = acetalAtomList.get(1);
		state.fragManager.createBond(neighbouringCarbon, atom2, 1);
		state.fragManager.incorporateFragment(acetalFrag, neighbouringCarbon.getFrag());
		return acetalFrag;
	}
	
	private void buildBiochemicalEster(BuildState state, List<Element> words, int numberOfWordRules) throws StructureBuildingException {
		for (Element word : words) {
			if (!WordType.full.toString().equals(word.getAttributeValue(TYPE_ATR))){
				throw new StructureBuildingException("Bug in word rule for biochemicalEster");
			}
			resolveWordOrBracket(state, word);
		}
		for (int i = 1; i < words.size(); i++) {
			Element ateWord = words.get(i);
			BuildResults br = new BuildResults(state, ateWord);
			String locant = ateWord.getAttributeValue(LOCANT_ATR);
			if (br.getFunctionalAtomCount()==0){
				throw new StructureBuildingException("Unable to find functional atom to form biochemical ester");
			}
			Atom functionalAtom =br.getFunctionalAtom(0);
			br.removeFunctionalAtom(0);
			Atom atomOnBiochemicalFragment;
			functionalAtom.neutraliseCharge();
			Fragment biochemicalFragment = state.xmlFragmentMap.get(findRightMostGroupInWordOrWordRule(words.get(0)));
			if (locant!=null){
				atomOnBiochemicalFragment = biochemicalFragment.getAtomByLocantOrThrow(locant);
				if (atomOnBiochemicalFragment.getBonds().size()!=1){
					atomOnBiochemicalFragment = biochemicalFragment.getAtomByLocantOrThrow("O" + locant);
				}
			}
			else{
				atomOnBiochemicalFragment = biochemicalFragment.getAtomByLocant("O5'");//take a guess at it being 5' ;-)
				if (atomOnBiochemicalFragment==null){
					List<Atom> atoms = biochemicalFragment.getAtomList();
					for (Atom atom : atoms) {
						if (atom.getElement().equals("O") && atom.getBonds().size()==1  && atom.getFirstBond().getOrder()==1){
							Atom adjacentAtom = atom.getAtomNeighbours().get(0);
							List<Atom> neighbours = adjacentAtom.getAtomNeighbours();
							if (adjacentAtom.getElement().equals("C") && neighbours.size()==3){
								neighbours.remove(atom);
								if (neighbours.get(0).getElement().equals("O") && adjacentAtom.getBondToAtomOrThrow(neighbours.get(0)).getOrder()==2){
									continue;
								}
								if (neighbours.get(1).getElement().equals("O") && adjacentAtom.getBondToAtomOrThrow(neighbours.get(1)).getOrder()==2){
									continue;
								}
							}
							atomOnBiochemicalFragment= atom;//find a hydroxy - not a carboxylic acid
						}
					}
				}
			}
			String element = atomOnBiochemicalFragment !=null ? atomOnBiochemicalFragment.getElement() : null;
			if (atomOnBiochemicalFragment ==null || 
					(atomOnBiochemicalFragment.getBonds().size()!=1 && !element.equals("O") && !element.equals("S") && !element.equals("Se") && !element.equals("Te"))){
				throw new StructureBuildingException("Failed to find hydroxy group on biochemical fragment");
			}
			state.fragManager.replaceAtomWithAnotherAtomPreservingConnectivity(functionalAtom, atomOnBiochemicalFragment);
			Element ateGroup = findRightMostGroupInWordOrWordRule(ateWord);
			if (ateGroup.getAttribute(NUMBEROFFUNCTIONALATOMSTOREMOVE_ATR)==null && numberOfWordRules==1){
				//by convention [O-] are implicitly converted to [OH] to balance charge
				for (int j = br.getFunctionalAtomCount() -1; j>=0; j--) {
					Atom atomToDefunctionalise =br.getFunctionalAtom(j);
					br.removeFunctionalAtom(j);
					atomToDefunctionalise.neutraliseCharge();
				}
			}
		}
	}

	private void connectSubstituentsToAcetal(BuildState state, List<Fragment> acetalFrags, BuildResults subBr, boolean hemiacetal) throws StructureBuildingException {
		Map<Fragment,Integer> usageMap= new HashMap<Fragment, Integer>();
		for (int i = subBr.getOutAtomCount() -1; i>=0; i--) {
			OutAtom out = subBr.getOutAtom(i);
			subBr.removeOutAtom(i);
			Atom atomToUse = null;
			if (out.getLocant()!=null){
				boolean numericLocant = MATCH_NUMERIC_LOCANT.matcher(out.getLocant()).matches();
				for (Fragment possibleAcetalFrag : acetalFrags) {
					if (numericLocant){
						Atom a  =OpsinTools.depthFirstSearchForNonSuffixAtomWithLocant(possibleAcetalFrag.getFirstAtom(), out.getLocant());
						if (a!=null){
							List<Atom> atomList =  possibleAcetalFrag.getAtomList();
							if (atomList.get(0).getBonds().size()==1){
								atomToUse = atomList.get(0);
								break;
							}
							else if (atomList.get(1).getBonds().size()==1){
								atomToUse = atomList.get(1);
								break;
							}
						}
					}
					else if (possibleAcetalFrag.hasLocant(out.getLocant())){
						atomToUse = possibleAcetalFrag.getAtomByLocantOrThrow(out.getLocant());
						break;
					}
				}
				if (atomToUse==null){
					throw new StructureBuildingException("Unable to find suitable acetalFrag");
				}
			}
			else{
				List<Atom> atomList =  acetalFrags.get(0).getAtomList();
				if (atomList.get(0).getBonds().size()==1){
					atomToUse = atomList.get(0);
				}
				else if (atomList.get(1).getBonds().size()==1){
					atomToUse = atomList.get(1);
				}
				else{
					throw new StructureBuildingException("OPSIN bug: unable to find acetal atom");
				}
			}
			Fragment acetalFrag = atomToUse.getFrag();
			int usage = usageMap.get(acetalFrag) !=null ? usageMap.get(acetalFrag) : 0;
			state.fragManager.createBond(out.getAtom(), atomToUse, out.getValency());
			usage++;
			if (usage >=2 || hemiacetal){
				acetalFrags.remove(acetalFrag);
			}
			usageMap.put(acetalFrag, usage);
		}
	}

	private List<Fragment> buildPolymer(BuildState state, List<Element> words) throws StructureBuildingException {
		if (words.size()!=2){
			throw new StructureBuildingException("Currently unsupported polymer name type");
		}
		Element polymer = words.get(1);
		resolveWordOrBracket(state, polymer);
		BuildResults polymerBr = new BuildResults(state, polymer);
		List<Fragment> rGroups = new ArrayList<Fragment>();
		if (polymerBr.getOutAtomCount() ==2 && polymerBr.getInAtomCount()==0){
			Atom inAtom =polymerBr.getOutAtomTakingIntoAccountWhetherSetExplicitly(0);
			Atom outAtom =polymerBr.getOutAtomTakingIntoAccountWhetherSetExplicitly(1);
			/*
			 * We assume the polymer repeats so as an approximation we create an R group with the same element as the group at the other end of polymer (with valency equal to the bondorder of the Rgroup so no H added)
			 */
			Fragment rGroup1 =state.fragManager.buildSMILES("[" + outAtom.getElement() + "|" + polymerBr.getOutAtom(0).getValency() + "]");
			state.fragManager.createBond(inAtom, rGroup1.getFirstAtom(), polymerBr.getOutAtom(0).getValency());

			Fragment rGroup2 =state.fragManager.buildSMILES("[" + inAtom.getElement() + "|" + polymerBr.getOutAtom(1).getValency() + "]");
			state.fragManager.createBond(outAtom, rGroup2.getFirstAtom(), polymerBr.getOutAtom(1).getValency());
			rGroups.add(rGroup1);
			rGroups.add(rGroup2);
			polymerBr.removeAllOutAtoms();
		}
		else{
			throw new StructureBuildingException("Polymer building failed: Two termini were not found; Expected 2 outAtoms, found: " +polymerBr.getOutAtomCount() +" ,expected 0 inAtoms, found: " +polymerBr.getInAtomCount());
		}
		return rGroups;
	}

	/**
	 * Finds a suitable functional atom corresponding to the given locant
	 * Takes into account situations where function replacement may have resulted in the wrong atoms being functional atoms
	 * @param locant
	 * @param mainGroupBR
	 * @return functionalAtomToUse
	 * @throws StructureBuildingException
	 */
	private Atom determineFunctionalAtomToUse(String locant, BuildResults mainGroupBR) throws StructureBuildingException {
		for (int i = 0; i < mainGroupBR.getFunctionalAtomCount(); i++) {
			//look for exact locant match
			Atom possibleAtom = mainGroupBR.getFunctionalAtom(i);
			if (possibleAtom.hasLocant(locant)){
				mainGroupBR.removeFunctionalAtom(i);
				if (possibleAtom.getProperty(Atom.AMBIGUOUS_ELEMENT_ASSIGNMENT)!=null){
					possibleAtom.getProperty(Atom.AMBIGUOUS_ELEMENT_ASSIGNMENT).remove(possibleAtom);
				}
				return possibleAtom;
			}
		}
		if (MATCH_NUMERIC_LOCANT.matcher(locant).matches()){
			//None of the functional atoms had an appropriate locant. Look for the case whether the locant refers to the backbone. e.g. 5-methyl 2-aminopentanedioate
			for (int i = 0; i < mainGroupBR.getFunctionalAtomCount(); i++) {
				Atom possibleAtom = mainGroupBR.getFunctionalAtom(i);
				if (OpsinTools.depthFirstSearchForNonSuffixAtomWithLocant(possibleAtom, locant)!=null){
					mainGroupBR.removeFunctionalAtom(i);
					if (possibleAtom.getProperty(Atom.AMBIGUOUS_ELEMENT_ASSIGNMENT)!=null){
						possibleAtom.getProperty(Atom.AMBIGUOUS_ELEMENT_ASSIGNMENT).remove(possibleAtom);
					}
					return possibleAtom;
				}
			}
		}
		else if (MATCH_ELEMENT_SYMBOL_LOCANT.matcher(locant).matches()){
			//None of the functional atoms had an appropriate locant. Look for the special cases:
			//	Where the lack of primes on an element symbol locant should be ignored e.g. O,O-diethyl carbonate
			//	Where the locant is used to decide on the ester configuration c.f. O-methyl ..thioate and S-methyl ..thioate
			boolean isElementSymbol = MATCH_ELEMENT_SYMBOL.matcher(locant).matches();
			for (int i = 0; i < mainGroupBR.getFunctionalAtomCount(); i++) {
				Atom possibleAtom = mainGroupBR.getFunctionalAtom(i);
				if (possibleAtom.getProperty(Atom.AMBIGUOUS_ELEMENT_ASSIGNMENT)!=null){
					Set<Atom> atoms =possibleAtom.getProperty(Atom.AMBIGUOUS_ELEMENT_ASSIGNMENT);
					boolean foundAtom = false;
                    for (Atom a : atoms) {
                        if (a.hasLocant(locant) || (isElementSymbol && a.getElement().equals(locant))){
                            //swap locants and element type
                            List<String> tempLocants = new ArrayList<String>(a.getLocants());
                            List<String> tempLocants2 = new ArrayList<String>(possibleAtom.getLocants());
                            a.clearLocants();
                            possibleAtom.clearLocants();
                            for (String l : tempLocants) {
                                possibleAtom.addLocant(l);
                            }
                            for (String l : tempLocants2) {
                                a.addLocant(l);
                            }
                            String originalElement = possibleAtom.getElement();
                            possibleAtom.setElement(a.getElement());
                            a.setElement(originalElement);
                            mainGroupBR.removeFunctionalAtom(i);
                            foundAtom =true;
                            break;
                        }
                    }
                    if (foundAtom){
                    	atoms.remove(possibleAtom);
                        return possibleAtom;
                    }
				}
				if (isElementSymbol && possibleAtom.getElement().equals(locant)){
				    return possibleAtom;
				}
			}
		}

		throw new StructureBuildingException("Cannot find functional atom with locant: " +locant + " to form an ester with");
	}

	/**
	 * Valency is used to determine the expected number of hydrogen
	 * Hydrogens are then added to bring the number of connections up to the minimum required to satisfy the atom's valency
	 * This allows the valency of the atom to be encoded e.g. phopshane-3 hydrogen, phosphorane-5 hydrogen.
	 * It is also neccesary when considering stereochemistry as a hydrogen beats nothing in the CIP rules
	 * @param state
	 * @throws StructureBuildingException
	 */
	static void makeHydrogensExplicit(BuildState state) throws StructureBuildingException {
		Set<Fragment> fragments = state.fragManager.getFragPile();
		for (Fragment fragment : fragments) {
			if (fragment.getSubType().equals(ELEMENTARYATOM_SUBTYPE_VAL)){//these do not have implicit hydrogen e.g. phosphorus is literally just a phosphorus atom
				continue;
			}
			List<Atom> atomList =fragment.getAtomList();
			for (Atom parentAtom : atomList) {
				int explicitHydrogensToAdd = StructureBuildingMethods.calculateSubstitutableHydrogenAtoms(parentAtom);
				for (int i = 0; i < explicitHydrogensToAdd; i++) {
					Atom hydrogen = state.fragManager.createAtom("H", fragment);
					state.fragManager.createBond(parentAtom, hydrogen, 1);
				}
				if (parentAtom.getAtomParity()!=null){
					if (explicitHydrogensToAdd >1){
						//Cannot have tetrahedral chirality and more than 2 hydrogens
						parentAtom.setAtomParity(null);//probably caused by deoxy
					}
					else{
						modifyAtomParityToTakeIntoAccountExplicitHydrogen(parentAtom);
					}
				}
			}
		}
	}

	private static void modifyAtomParityToTakeIntoAccountExplicitHydrogen(Atom atom) throws StructureBuildingException {
		AtomParity atomParity = atom.getAtomParity();
		if (!StereoAnalyser.isPossiblyStereogenic(atom)){
			//no longer a stereoCentre e.g. due to unsaturation
			atom.setAtomParity(null);
		}
		else{
			Atom[] atomRefs4 = atomParity.getAtomRefs4();
			Integer positionOfImplicitHydrogen = null;
			Integer positionOfDeoxyHydrogen = null;
			for (int i = 0; i < atomRefs4.length; i++) {
				if (atomRefs4[i].equals(AtomParity.hydrogen)){
					positionOfImplicitHydrogen = i;
				}
				else if (atomRefs4[i].equals(AtomParity.deoxyHydrogen)){
					positionOfDeoxyHydrogen = i;
				}
			}
			if (positionOfImplicitHydrogen !=null || positionOfDeoxyHydrogen !=null){
				//atom parity was set in SMILES, the dummy hydrogen atom has now been substituted
				List<Atom> neighbours = atom.getAtomNeighbours();
				for (Atom atomRef : atomRefs4) {
					neighbours.remove(atomRef);
				}
				if (neighbours.size()==0){
					throw new StructureBuildingException("OPSIN Bug: Unable to determine which atom has substitued a hydrogen at stereocentre");
				}
				else if (neighbours.size()==1 && positionOfDeoxyHydrogen!=null){
					atomRefs4[positionOfDeoxyHydrogen] = neighbours.get(0);
					if (positionOfImplicitHydrogen != null){
						throw new StructureBuildingException("OPSIN Bug: Unable to determine which atom has substitued a hydrogen at stereocentre");
					}
				}
				else if (neighbours.size()==1 && positionOfImplicitHydrogen!=null){
					atomRefs4[positionOfImplicitHydrogen] = neighbours.get(0);
				}
				else if (neighbours.size()==2 && positionOfDeoxyHydrogen!=null && positionOfImplicitHydrogen!=null){
					//TODO get CIP priority of neighbours. Refactor stereoanalyzer into two classes so that get CIP priority is decoupled from symmetry analysis
					atomRefs4[positionOfDeoxyHydrogen] = neighbours.get(0);
					atomRefs4[positionOfImplicitHydrogen] = neighbours.get(1);
				}
				else{
					throw new StructureBuildingException("OPSIN Bug: Unable to determine which atom has substitued a hydrogen at stereocentre");
				}
			}
		}
	}

	private boolean applyExplicitStoichometryIfProvided(BuildState state, Elements wordRules) throws StructureBuildingException {
		boolean explicitStoichometryPresent =false;
		for (int i = 0; i < wordRules.size(); i++) {
			Element wordRule = wordRules.get(i);
			if (wordRule.getAttribute(STOICHOMETRY_ATR)!=null){
				int stoichometry = Integer.parseInt(wordRule.getAttributeValue(STOICHOMETRY_ATR));
				wordRule.removeAttribute(wordRule.getAttribute(STOICHOMETRY_ATR));
				for (int j = 1; j < stoichometry; j++) {
					Element clone = state.fragManager.cloneElement(state, wordRule);
					XOMTools.insertAfter(wordRule, clone);
				}
				explicitStoichometryPresent =true;
			}
		}
		return explicitStoichometryPresent;
	}

	/**
	 * A net charge is present; Given the molecule element the overallCharge is there an unambiguous way of 
	 * multiplying fragments to make the net charge 0
	 * metals without specified charge may be given an implicit positive charge
	 * 
	 * If this fails look for the case where there are multiple molecules and the mixture is only negative due to negatively charged functional Atoms e.g. pyridine acetate and remove the negative charge
	 * @param state
	 * @param molecule
	 * @param explicitStoichometryPresent 
	 * @param overallCharge
	 * @throws StructureBuildingException
	 */
	private void balanceChargeIfPossible(BuildState state, Element molecule, int overallCharge, boolean explicitStoichometryPresent) throws StructureBuildingException {
		List<Element> wordRules = XOMTools.getChildElementsWithTagName(molecule, WORDRULE_ATR);

		List<Element> positivelyChargedComponents = new ArrayList<Element>();
		List<Element> negativelyChargedComponents = new ArrayList<Element>();
		HashMap<Element, Integer> componentToChargeMapping = new HashMap<Element, Integer>();
		HashMap<Element, BuildResults> componentToBR = new HashMap<Element, BuildResults>();
		
		List<Element> cationicElements = new ArrayList<Element>();
		List<Element> elementaryAtoms = XOMTools.getDescendantElementsWithTagNameAndAttribute(molecule, GROUP_EL, SUBTYPE_ATR, ELEMENTARYATOM_SUBTYPE_VAL);
		for (Element elementaryAtom : elementaryAtoms) {
			if (elementaryAtom.getAttribute(COMMONOXIDATIONSTATESANDMAX_ATR)!=null){
				Fragment cationicFrag =state.xmlFragmentMap.get(elementaryAtom);
				if (cationicFrag.getFirstAtom().getCharge()==0){//if not 0 charge cannot be implicitly modified
					String[] typicalOxidationStates = MATCH_COMMA.split(MATCH_COLON.split(elementaryAtom.getAttributeValue(COMMONOXIDATIONSTATESANDMAX_ATR))[0]);
					int typicalCharge = Integer.parseInt(typicalOxidationStates[typicalOxidationStates.length-1]);
					if (typicalCharge > cationicFrag.getFirstAtom().getAtomNeighbours().size()){
						cationicElements.add(elementaryAtom);
					}
				}
			}
		}
		overallCharge = setCationicElementsToTypicalCharge(state, cationicElements, overallCharge);
		if (overallCharge==0){
			return;
		}
		if (cationicElements.size() ==1 && overallCharge <0){//e.g. nickel tetrachloride [Ni2+]-->[Ni4+]
			boolean success = setChargeOnCationicElementAppropriately(state, overallCharge, cationicElements.get(0));
			if (success){
				return;
			}
		}
		for (Element wordRule : wordRules) {
			BuildResults br = new BuildResults(state, wordRule);
			componentToBR.put(wordRule, br);
			int charge = br.getCharge();
			if (charge>0){
				positivelyChargedComponents.add(wordRule);
			}
			else if (charge <0){
				negativelyChargedComponents.add(wordRule);
			}
			componentToChargeMapping.put(wordRule, charge);
		}
		if (!explicitStoichometryPresent &&
				(positivelyChargedComponents.size()==1 && cationicElements.size() ==0 && negativelyChargedComponents.size() >=1 || positivelyChargedComponents.size()>=1 && negativelyChargedComponents.size() ==1 )){
			boolean success = multiplyChargedComponents(state, negativelyChargedComponents, positivelyChargedComponents, componentToChargeMapping, overallCharge);
			if (success){
				return;
			}
		}
		if (cationicElements.size() ==1){//e.g. magnesium monochloride [Mg2+]-->[Mg+]
			boolean success = setChargeOnCationicElementAppropriately(state, overallCharge, cationicElements.get(0));
			if (success){
				return;
			}
		}
		if (overallCharge <0){//neutralise functionalAtoms if they are the sole cause of the negative charge and multiple molecules are present
			int chargeOnFunctionalAtoms = 0;
			for (Element wordRule : wordRules) {
				BuildResults br = componentToBR.get(wordRule);
				int functionalAtomCount = br.getFunctionalAtomCount();
				for (int i = functionalAtomCount -1; i >=0; i--) {
					chargeOnFunctionalAtoms += br.getFunctionalAtom(i).getCharge();
				}
			}
			if (chargeOnFunctionalAtoms <= overallCharge){
				for (Element wordRule : wordRules) {
					BuildResults br = componentToBR.get(wordRule);
					int functionalAtomCount = br.getFunctionalAtomCount();
					for (int i = functionalAtomCount -1; i >=0; i--) {
						if (overallCharge==0){
							return;
						}
						overallCharge-=br.getFunctionalAtom(i).getCharge();
						br.getFunctionalAtom(i).neutraliseCharge();
						br.removeFunctionalAtom(i);
					}
				}
			}
		}
	}


	/**
	 * Sets the cationicElements to the lowest typical charge as specified by the COMMONOXIDATIONSTATESANDMAX_ATR that is >= incoming valency
	 * The valency incoming to the cationicElement is taken into account e.g. phenylmagnesium chloride is [Mg+]
	 * @param state
	 * @param cationicElements
	 * @param overallCharge
	 * @return
	 */
	private int setCationicElementsToTypicalCharge(BuildState state, List<Element> cationicElements, int overallCharge)  {
		for (Element cationicElement : cationicElements) {
			Fragment cationicFrag = state.xmlFragmentMap.get(cationicElement);
			String[] typicalOxidationStates = MATCH_COMMA.split(MATCH_COLON.split(cationicElement.getAttributeValue(COMMONOXIDATIONSTATESANDMAX_ATR))[0]);
			int incomingValency = cationicFrag.getFirstAtom().getIncomingValency();
			for (String typicalOxidationState : typicalOxidationStates) {
				int charge = Integer.parseInt(typicalOxidationState);
				if (charge>= incomingValency){
					charge -= incomingValency;
					overallCharge += charge;
					cationicFrag.getFirstAtom().setCharge(charge);
					break;
				}
			}
		}
		return overallCharge;
	}

	/**
	 * Multiplies out charged word rules to balance charge
	 * Return true if balancing was possible else false
	 * @param state
	 * @param negativelyChargedComponents
	 * @param positivelyChargedComponents
	 * @param componentToChargeMapping
	 * @param overallCharge
	 * @return
	 * @throws StructureBuildingException
	 */
	private boolean multiplyChargedComponents(BuildState state, List<Element>negativelyChargedComponents,List<Element> positivelyChargedComponents,HashMap<Element, Integer> componentToChargeMapping, int overallCharge) throws StructureBuildingException {
		Element componentToMultiply;
		if (overallCharge >0){
			if (negativelyChargedComponents.size() >1){
				return false;//ambiguous as to which to multiply
			}
			componentToMultiply = negativelyChargedComponents.get(0);
		}
		else{
			if (positivelyChargedComponents.size() >1){
				return false;//ambiguous as to which to multiply
			}
			componentToMultiply = positivelyChargedComponents.get(0);
		}

		int charge = componentToChargeMapping.get(componentToMultiply);
		if (overallCharge % charge ==0){//e.g. magnesium chloride
			if (!componentCanBeMultiplied(componentToMultiply)){
				return false;
			}
			int timesToDuplicate = Math.abs(overallCharge/charge);
			for (int i = 0; i < timesToDuplicate; i++) {
				XOMTools.insertAfter(componentToMultiply, state.fragManager.cloneElement(state, componentToMultiply));
			}
		}
		else{//e.g. iron(3+) sulfate -->2:3 mixture
			if (positivelyChargedComponents.size() >1 || !componentCanBeMultiplied(positivelyChargedComponents.get(0))){
				return false;
			}
			if (negativelyChargedComponents.size() >1 || !componentCanBeMultiplied(negativelyChargedComponents.get(0))){
				return false;
			}
			int positiveCharge = componentToChargeMapping.get(positivelyChargedComponents.get(0));
			int negativeCharge = Math.abs(componentToChargeMapping.get(negativelyChargedComponents.get(0)));
			int targetTotalAbsoluteCharge = positiveCharge * negativeCharge;
			for (int i = (targetTotalAbsoluteCharge/negativeCharge); i >1; i--) {
				XOMTools.insertAfter(negativelyChargedComponents.get(0), state.fragManager.cloneElement(state, negativelyChargedComponents.get(0)));
			}
			for (int i = (targetTotalAbsoluteCharge/positiveCharge); i >1; i--) {
				XOMTools.insertAfter(positivelyChargedComponents.get(0), state.fragManager.cloneElement(state, positivelyChargedComponents.get(0)));
			}
		}
		return true;
	}

	private boolean componentCanBeMultiplied(Element componentToMultiply) {
		if (componentToMultiply.getAttributeValue(WORDRULE_ATR).equals(WordRule.simple.toString()) && XOMTools.getChildElementsWithTagNameAndAttribute(componentToMultiply, WORD_EL, TYPE_ATR, WordType.full.toString()).size()>1){
			return false;//already has been multiplied e.g. dichloride
		}
		Element firstChild = (Element) componentToMultiply.getChild(0);
		while (firstChild.getChildElements().size() !=0){
			firstChild = (Element) firstChild.getChild(0);
		}
		if (firstChild.getLocalName().equals(MULTIPLIER_EL)){//e.g. monochloride. Allows specification of explicit stoichiometry
			return false;
		}
		return true;
	}
	
	private boolean setChargeOnCationicElementAppropriately(BuildState state, int overallCharge, Element cationicElement)  {
		Atom cation = state.xmlFragmentMap.get(cationicElement).getFirstAtom();
		int chargeOnCationNeeded = -(overallCharge -cation.getCharge());
		int maximumCharge = Integer.parseInt(MATCH_COLON.split(cationicElement.getAttributeValue(COMMONOXIDATIONSTATESANDMAX_ATR))[1]);
		if (chargeOnCationNeeded >=0 && chargeOnCationNeeded <= maximumCharge){
			cation.setCharge(chargeOnCationNeeded);
			return true;
		}
		return false;
	}

	private Element findRightMostGroupInWordOrWordRule(Element wordOrWordRule) throws StructureBuildingException {
		if (wordOrWordRule.getLocalName().equals(WORDRULE_EL)){
			List<Element> words = XOMTools.getDescendantElementsWithTagName(wordOrWordRule, WORD_EL);
			for (int i = words.size() -1 ; i >=0; i--) {//ignore functionalTerm Words
				if (words.get(i).getAttributeValue(TYPE_ATR).equals(WordType.functionalTerm.toString())){
					words.remove(words.get(i));
				}
			}
			if (words.size()==0){
				throw new StructureBuildingException("OPSIN bug: word element not found where expected");
			}
			return StructureBuildingMethods.findRightMostGroupInBracket(words.get(words.size()-1));
		}
		else if (wordOrWordRule.getLocalName().equals(WORD_EL)){//word element can be treated just like a bracket
			return StructureBuildingMethods.findRightMostGroupInBracket(wordOrWordRule);
		}
		else{
			throw new StructureBuildingException("OPSIN bug: expected word or wordRule");
		}
	}

	/**
	 * Nasty special case to cope with oxido and related groups acting as O= or even [O-][N+]
	 * This nasty behaviour is in generated ChemDraw names and is supported by most nameToStructure tools so it is supported here
	 * Acting as O= notably is often correct behaviour for inorganics
	 * @param state
	 * @param groups
	 */
	private void processOxidoSpecialCase(BuildState state, List<Element> groups)  {
		for (Element group : groups) {
			if (OXIDOLIKE_SUBTYPE_VAL.equals(group.getAttributeValue(SUBTYPE_ATR))){
				Atom oxidoAtom = state.xmlFragmentMap.get(group).getFirstAtom();
				Atom connectedAtom = oxidoAtom.getAtomNeighbours().get(0);
				String element = connectedAtom.getElement();
				if (checkForConnectedOxo(state, connectedAtom)){//e.g. not oxido(trioxo)ruthenium
					continue;
				}
				if (ELEMENTARYATOM_SUBTYPE_VAL.equals(connectedAtom.getFrag().getSubType()) ||
						((element.equals("S") || element.equals("P")) && connectedAtom.getCharge() ==0 && ValencyChecker.checkValencyAvailableForBond(connectedAtom, 1))){
					oxidoAtom.neutraliseCharge();
					oxidoAtom.getFirstBond().setOrder(2);
				}
				else if (element.equals("N") && connectedAtom.getCharge()==0){
					int incomingValency = connectedAtom.getIncomingValency();
					if ((incomingValency + connectedAtom.getOutValency()) ==3 && connectedAtom.hasSpareValency()){
						connectedAtom.addChargeAndProtons(1, 1);//e.g. N-oxidopyridine
					}
					else if ((incomingValency + connectedAtom.getOutValency()) ==4){
						if (connectedAtom.getLambdaConventionValency()!=null && connectedAtom.getLambdaConventionValency()==5){
							oxidoAtom.setCharge(0);
							oxidoAtom.setProtonsExplicitlyAddedOrRemoved(0);
							oxidoAtom.getFirstBond().setOrder(2);
						}
						else{
							connectedAtom.addChargeAndProtons(1, 1);
						}
					}
				}
			}
		}
	}

	/**
	 * Is the atom connected to an atom whose fragment has an xml entry called "oxo"
	 * @param atom
	 * @return
	 */
	private boolean checkForConnectedOxo(BuildState state, Atom atom) {
		Set<Bond> bonds = atom.getBonds();
		for (Bond bond : bonds) {
			Atom connectedAtom;
			if (bond.getFromAtom() == atom){
				connectedAtom = bond.getToAtom();
			}
			else{
				connectedAtom = bond.getFromAtom();
			}
			Element correspondingEl = state.xmlFragmentMap.getElement(connectedAtom.getFrag());
			if (correspondingEl.getValue().equals("oxo")){
				return true;
			}
		}
		return false;
	}
	

	/**
	 * Sets the charge according to the oxidation number if the oxidation number atom property has been set
	 * @param state
	 * @param groups
	 * @throws StructureBuildingException 
	 */
	private void processOxidationNumbers(BuildState state, List<Element> groups) throws StructureBuildingException  {
		for (Element group : groups) {
			if (ELEMENTARYATOM_SUBTYPE_VAL.equals(group.getAttributeValue(SUBTYPE_ATR))){
				Atom atom = state.xmlFragmentMap.get(group).getFirstAtom();
				if (atom.getProperty(Atom.OXIDATION_NUMBER)!=null){
					List<Atom> neighbours = atom.getAtomNeighbours();
					int chargeThatWouldFormIfLigandsWereRemoved =0;
					for (Atom neighbour : neighbours) {
						Element neighbourEl = state.xmlFragmentMap.getElement(neighbour.getFrag());
						Bond b = atom.getBondToAtomOrThrow(neighbour);
						//carbonyl and nitrosyl are neutral ligands
						if (!((neighbourEl.getValue().equals("carbon") && NONCARBOXYLICACID_TYPE_VAL.equals(neighbourEl.getAttributeValue(TYPE_ATR)))
								|| neighbourEl.getValue().equals("nitrosyl"))){
							chargeThatWouldFormIfLigandsWereRemoved+=b.getOrder();
						}
					}
					
					atom.setCharge(atom.getProperty(Atom.OXIDATION_NUMBER)-chargeThatWouldFormIfLigandsWereRemoved);
				}
			}
		}
	}

	/**
	 * Handles the application of stereochemistry and checking
	 * existing stereochemical specification is still relevant.
	 * @param state
	 * @param molecule
	 * @param uniFrag
	 * @throws StructureBuildingException
	 */
	private void processStereochemistry(BuildState state, Element molecule, Fragment uniFrag) throws StructureBuildingException {
		List<Element> stereoChemistryEls = findStereochemistryElsInProcessingOrder(molecule);
		List<Atom> atomList = uniFrag.getAtomList();
		List<Atom> atomsWithPreDefinedAtomParity = new ArrayList<Atom>();
		for (Atom atom : atomList) {
			if (atom.getAtomParity()!=null){
				atomsWithPreDefinedAtomParity.add(atom);
			}
		}
		Set<Bond> bonds = uniFrag.getBondSet();
		List<Bond> bondsWithPreDefinedBondStereo = new ArrayList<Bond>();
		for (Bond bond : bonds) {
			if (bond.getBondStereo()!=null){
				bondsWithPreDefinedBondStereo.add(bond);
			}
		}
		if (stereoChemistryEls.size() >0 || atomsWithPreDefinedAtomParity.size() >0 || bondsWithPreDefinedBondStereo.size() >0){
			StereochemistryHandler.processStereochemicalElements(state, uniFrag, stereoChemistryEls, atomsWithPreDefinedAtomParity, bondsWithPreDefinedBondStereo);
		}
	}

	/**
	 * Finds stereochemistry els in a recursive right to left manner.
	 * Within the same scope though stereochemistry els are found left to right
	 * @param molecule
	 * @return
	 */
	private List<Element> findStereochemistryElsInProcessingOrder(Element parentEl) {
		List<Element> matchingElements = new ArrayList<Element>();
		Elements children =parentEl.getChildElements();
		List<Element> stereochemistryElsAtThisLevel = new ArrayList<Element>();
		for (int i = children.size()-1; i >=0; i--) {
			Element child = children.get(i);
			if (child.getLocalName().equals(STEREOCHEMISTRY_EL)){
				stereochemistryElsAtThisLevel.add(child);
			}
			else{
				matchingElements.addAll(findStereochemistryElsInProcessingOrder(child));
			}
		}
		Collections.reverse(stereochemistryElsAtThisLevel);
		matchingElements.addAll(stereochemistryElsAtThisLevel);
		return matchingElements;
	}
}
