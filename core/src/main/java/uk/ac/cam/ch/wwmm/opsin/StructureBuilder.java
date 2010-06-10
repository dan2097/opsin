package uk.ac.cam.ch.wwmm.opsin;

import java.util.*;
import java.util.regex.Pattern;

import uk.ac.cam.ch.wwmm.opsin.ParseWord.WordType;
import uk.ac.cam.ch.wwmm.opsin.WordRules.WordRule;

import nu.xom.Element;
import nu.xom.Elements;
import static uk.ac.cam.ch.wwmm.opsin.XmlDeclarations.*;
import static uk.ac.cam.ch.wwmm.opsin.StructureBuildingMethods.*;

/**Constructs a single OPSIN fragment which describes the molecule from the postprocessor results.
 *
 * @author ptc24/dl387
 *
 */
class StructureBuilder {

	private final Pattern matchDigits = Pattern.compile("\\d+");
	private final Pattern matchComma =Pattern.compile(",");
	private final Pattern matchColon =Pattern.compile(":");
	private final Pattern matchNumericLocant =Pattern.compile("\\d+[a-z]?'*");
	private final Pattern matchElementSymbolLocant =Pattern.compile("[A-Z][a-z]?'*");
	private final Pattern matchElementSymbol =Pattern.compile("[A-Z][a-z]?");

	/**	Builds a molecule as a Fragment based on preStructurebuilder output.
	 * @param state
	 * @param molecule The preStructurebuilderr output.
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
			if(wordRule == WordRule.simple) {
				for (Element word : words) {
					if (!word.getLocalName().equals(WORD_EL) || !word.getAttributeValue(TYPE_ATR).equals(WordType.full.toString())){
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
				buildAmide(state, words);//e.g. ethanoic acid ethyl amide, terephthalic acid dimethyl amide, ethanoic acid amide
			}
			else if (wordRule == WordRule.glycol){
				buildGlycol(state, words);//e.g. ethylene glycol
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
			else if(wordRule == WordRule.acetal) {
				buildAcetal(state, words);//e.g. propanal diethyl acetal
			}
			else if(wordRule == WordRule.polymer) {
				rGroups.addAll(buildPolymer(state, words));
			}
			else{
				throw new StructureBuildingException("Unknown Word Rule");
			}
		}

		state.fragManager.convertSpareValenciesToDoubleBonds();
		state.fragManager.checkValencies();
		int overallCharge = state.fragManager.getOverallCharge();
		if (overallCharge!=0 && wordRules.size() >1){//a net charge is present! Could just mean the counterion has not been specified though
			balanceChargeIfPossible(state, molecule, overallCharge);
		}
		makeHydrogensExplicit(state);

		Fragment uniFrag = state.fragManager.getUnifiedFragment();
		List<Element> stereoChemistryEls = XOMTools.getDescendantElementsWithTagName(molecule, STEREOCHEMISTRY_EL);
		if (stereoChemistryEls.size() >0){
			StereochemistryHandler.processStereochemicalElements(state, uniFrag, stereoChemistryEls);
		}

		for (Fragment rGroup : rGroups) {
			Atom rAtom = rGroup.getFirstAtom();
			rAtom.setElement("R");
		}

		if (uniFrag.getOutAtoms().size()>0 || uniFrag.getInAtoms().size()>0){
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
		for (; wordIndice < words.size(); wordIndice++) {
			Element word =words.get(wordIndice);
			if (word.getAttributeValue(TYPE_ATR).equals(WordType.full.toString())){
				resolveWordOrBracket(state, word);
				BuildResults ateBR = new BuildResults(state, word);
				if (ateBR.getFunctionalAtomCount()<1){
					throw new StructureBuildingException("bug? ate group did not have any functional atoms!");
				}
				ateGroups.add(ateBR);
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
			state.fragManager.createBond(ateAtom,substituentsBr.getOutAtomTakingIntoAccountWhetherSetExplicitly(0), 1);
			substituentsBr.removeOutAtom(0);
			ateAtom.setCharge(0);
		}
	}



	private void buildDiValentFunctionalGroup(BuildState state, List<Element> words) throws StructureBuildingException {
		int wordIndice = 0;
		if (!words.get(wordIndice).getAttributeValue(TYPE_ATR).equals(WordType.substituent.toString())) {
			throw new StructureBuildingException("word: " +wordIndice +" was expected to be a substituent");
		}
		if (words.get(wordIndice +1).getAttributeValue(TYPE_ATR).equals(WordType.functionalTerm.toString())) {//e.g. methyl sulfoxide rather than dimethyl sulfoxide
			Element clone = state.fragManager.cloneElement(state, words.get(0));
			XOMTools.insertAfter(words.get(0), clone);
			words = OpsinTools.elementsToElementArrayList(((Element)words.get(0).getParent()).getChildElements());
		}
		resolveWordOrBracket(state, words.get(wordIndice));
		BuildResults substituent1 =new BuildResults(state, words.get(wordIndice));
		if (substituent1.getOutAtomCount()!=1){
			throw new StructureBuildingException("Expected one outAtom. Found " + substituent1.getOutAtomCount() );
		}
		if (substituent1.getOutAtom(0).getValency() !=1){
			throw new StructureBuildingException("OutAtom has unexpected valency. Expected 1. Actual: " + substituent1.getOutAtom(0).getValency());
		}
		wordIndice++;
		resolveWordOrBracket(state, words.get(wordIndice));
		BuildResults substituent2 =new BuildResults(state, words.get(wordIndice));
		if (substituent2.getOutAtomCount()!=1){
			throw new StructureBuildingException("Expected one outAtom. Found " + substituent2.getOutAtomCount() );
		}
		if (substituent2.getOutAtom(0).getValency() !=1){
			throw new StructureBuildingException("OutAtom has unexpected valency. Expected 1. Actual: " + substituent2.getOutAtom(0).getValency());
		}
		wordIndice++;
		if (words.get(wordIndice) ==null || !words.get(wordIndice).getAttributeValue(TYPE_ATR).equals(WordType.functionalTerm.toString())) {
			throw new StructureBuildingException("word: " +wordIndice +" was expected to be a functionalTerm");
		}
		List<Element> functionalGroup = XOMTools.getDescendantElementsWithTagName(words.get(wordIndice), "functionalGroup");
		if (functionalGroup.size()!=1){
			throw new StructureBuildingException("Unexpected number of functionalGroups found, could be a bug in OPSIN's grammar");
		}
		String smilesOfGroup = functionalGroup.get(0).getAttributeValue(VALUE_ATR);
		Fragment diValentGroup =state.fragManager.buildSMILES(smilesOfGroup, FUNCTIONALCLASS_TYPE_VAL, NONE_LABELS_VAL);

		Atom outAtom =substituent1.getOutAtomTakingIntoAccountWhetherSetExplicitly(0);
		substituent1.removeOutAtom(0);
		if (diValentGroup.getOutAtoms().size()==1){//c.f. peroxide where it is a linker
			state.fragManager.createBond(outAtom, diValentGroup.getOutAtom(0).getAtom(), 1);
			diValentGroup.removeOutAtom(0);
		}
		else{
			state.fragManager.createBond(outAtom, diValentGroup.getAtomByIDOrThrow(diValentGroup.getIdOfFirstAtom()), 1);
		}
		outAtom = substituent2.getOutAtomTakingIntoAccountWhetherSetExplicitly(0);
		substituent2.removeOutAtom(0);
		state.fragManager.createBond(outAtom, diValentGroup.getAtomByIDOrThrow(diValentGroup.getIdOfFirstAtom()), 1);
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
				ideAtom.setCharge(ideAtom.getCharge()+1);//e.g. make cyanide charge netural
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

	private void buildAmide(BuildState state, List<Element> words) throws StructureBuildingException {
		if (!words.get(0).getAttributeValue(TYPE_ATR).equals(WordType.full.toString())){
			throw new StructureBuildingException("Don't alter wordRules.xml without checking the consequences!");
		}
		BuildResults acidBr = new BuildResults(state, words.get(0));
		if (acidBr.getFunctionalAtomCount()==0){
			throw new StructureBuildingException("No functionalAtoms detected!");
		}
		int wordIndice =1;
		
		if (words.get(wordIndice).getAttributeValue(TYPE_ATR).equals(WordType.functionalTerm.toString())){//"amide"
			for (int i =0; i < acidBr.getFunctionalAtomCount(); i++) {
				Atom functionalAtom = acidBr.getFunctionalAtom(i);
				if (!functionalAtom.getElement().equals("O")){
					throw new StructureBuildingException("Expected oxygen functional atom found:" + functionalAtom.getElement());
				}
				functionalAtom.setElement("N");
				functionalAtom.replaceLocant("N" +StringTools.multiplyString("'", i));
			}
		}
		else if (words.get(wordIndice).getAttributeValue(TYPE_ATR).equals(WordType.full.toString())){//substituentamide
			for (; wordIndice < words.size(); wordIndice++) {
				resolveWordOrBracket(state, words.get(wordIndice));//the substituted amide ([N-]) group
				List<Element> root = XOMTools.getDescendantElementsWithTagName(words.get(wordIndice), ROOT_EL);
				if (root.size()!=1){
					throw new StructureBuildingException("Cannot find root element");
				}
				Fragment amide = state.xmlFragmentMap.get(root.get(0).getFirstChildElement(GROUP_EL));
				if (acidBr.getFunctionalAtomCount()==0){
					throw new StructureBuildingException("Cannot find oxygen to replace with nitrogen for amide formation!");
				}
				Atom functionalAtom = acidBr.getFunctionalAtom(0);
				if (!functionalAtom.getElement().equals("O")){
					throw new StructureBuildingException("Expected oxygen functional atom found:" + functionalAtom.getElement());
				}
				Atom amideNitrogen = null;
				for (Atom a : amide.getAtomList()) {
					if (a.getElement().equals("N")){
						amideNitrogen = a;
						break;
					}
				}
				if (amideNitrogen ==null){
					throw new StructureBuildingException("Cannot find nitrogen atom in amide!");
				}
				acidBr.removeFunctionalAtom(0);
				amideNitrogen.setCharge(0);
				state.fragManager.replaceTerminalAtomWithFragment(functionalAtom, amideNitrogen);
			}
		}
		if (wordIndice +1 < words.size()){
			throw new StructureBuildingException("More words than expected when applying amide word rule!");
		}
		resolveWordOrBracket(state, words.get(0));//the group that had amide functionality added
	}

	private void buildGlycol(BuildState state, List<Element> words) throws StructureBuildingException {
		int wordIndice  = 0;
		resolveWordOrBracket(state, words.get(wordIndice));//the group
		BuildResults theDiRadical = new BuildResults(state, words.get(wordIndice));
		if (theDiRadical.getOutAtomCount()!=2){
			throw new StructureBuildingException("Glycol class names (e.g. ethylene glycol) expect two outAtoms. Found: " + theDiRadical.getOutAtomCount() );
		}
		if (wordIndice +1 >= words.size() || !words.get(wordIndice+1).getAttributeValue(TYPE_ATR).equals(WordType.functionalTerm.toString())){
			throw new StructureBuildingException("Glycol functionalTerm word expected");
		}
		for (int i = theDiRadical.getOutAtomCount() -1; i >=0 ; i--) {
			Atom outAtom =theDiRadical.getOutAtomTakingIntoAccountWhetherSetExplicitly(0);
			Fragment glycol =state.fragManager.buildSMILES("O", FUNCTIONALCLASS_TYPE_VAL, NONE_LABELS_VAL);
			if (theDiRadical.getOutAtom(0).getValency() !=1){
				throw new StructureBuildingException("OutAtom has unexpected valency. Expected 1. Actual: " + theDiRadical.getOutAtom(0).getValency());
			}
			theDiRadical.removeOutAtom(0);
			state.fragManager.createBond(outAtom, glycol.getAtomByIDOrThrow(glycol.getIdOfFirstAtom()), 1);
		}
	}
	
	private void buildOxide(BuildState state, List<Element> words) throws StructureBuildingException {
		resolveWordOrBracket(state, words.get(0));//the group
		List<Fragment> oxideFragments = new ArrayList<Fragment>();
		List<String> locantsForOxide =new ArrayList<String>();//often not specified
		if (!words.get(1).getAttributeValue(TYPE_ATR).equals(WordType.functionalTerm.toString())){
			throw new StructureBuildingException("Oxide functional term not found where expected!");
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
		List<Element> functionalClass =XOMTools.getDescendantElementsWithTagName(words.get(1), GROUP_EL);
		if (functionalClass.size()!=1){
			throw new StructureBuildingException("Expected 1 group element found: " + functionalClass.size());
		}
		String smilesReplacement = functionalClass.get(0).getAttributeValue(VALUE_ATR);
		String labels =  functionalClass.get(0).getAttributeValue(LABELS_ATR);
		for (int i = 0; i < numberOfOxygenToAdd; i++) {
			oxideFragments.add(state.fragManager.buildSMILES(smilesReplacement, "", labels));
		}
		List<Element> locantEls =XOMTools.getDescendantElementsWithTagName(words.get(1), LOCANT_EL);
		if (locantEls.size() >1){
			throw new StructureBuildingException("Expected 0 or 1 locant elements found: " + locantEls.size());
		}
		if (locantEls.size()==1){
			String[] locants = matchComma.split(StringTools.removeDashIfPresent(locantEls.get(0).getValue()));
            locantsForOxide.addAll(Arrays.asList(locants));
			locantEls.get(0).detach();
		}
		if (!locantsForOxide.isEmpty() && locantsForOxide.size()!=oxideFragments.size()){
			throw new StructureBuildingException("Mismatch between number of locants and number of oxides specified");
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
					List<Atom> atomList = frag.getAtomList();
					for (Atom atom : atomList) {
						if (!atom.getElement().equals("C") && !atom.getElement().equals("O")){
							formAppropriateBondToOxideAndAdjustCharges(state, atom, oxideAtom);
							continue mainLoop;
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
				for (Fragment frag : orderedPossibleFragments) {//something like carbon dioxide
					String subTypeVal = state.xmlFragmentMap.getElement(frag).getAttributeValue(SUBTYPE_ATR);
					if (ELEMENTARYATOMINORGANIC_SUBTYPE_VAL.equals(subTypeVal) || ELEMENTARYATOMORGANIC_SUBTYPE_VAL.equals(subTypeVal)){
						formAppropriateBondToOxideAndAdjustCharges(state, frag.getFirstAtom(), oxideAtom);
						continue mainLoop;
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
		if (ValencyChecker.checkValencyAvailableForBond(atomToAddOxideTo, 2)){
			state.fragManager.createBond(atomToAddOxideTo, oxideAtom, 2);
		}
		else{
			if (atomToAddOxideTo.getCharge()!=0 || oxideAtom.getCharge()!=0){
				throw new StructureBuildingException("Oxide appeared to refer to an atom that has insufficent valency to accept the addition of oxygen");
			}
			atomToAddOxideTo.addChargeAndProtons(1, 1);
			oxideAtom.setCharge(-1);
			if (!ValencyChecker.checkValencyAvailableForBond(atomToAddOxideTo, 1)){
				throw new StructureBuildingException("Oxide appeared to refer to an atom that has insufficent valency to accept the addition of oxygen");
			}
			state.fragManager.createBond(atomToAddOxideTo, oxideAtom, 1);
		}
	}

	private void buildCarbonylDerivative(BuildState state, List<Element> words) throws StructureBuildingException {
		if (!WordType.full.toString().equals(words.get(0).getAttributeValue(TYPE_ATR))){
			throw new StructureBuildingException("OPSIN bug: Wrong word type encountered when applying carbonylDerivative wordRule");
		}
		Element rightmostGroup = findRightMostGroupInWordOrWordRule(words.get(0));
		Fragment rootFragment = state.xmlFragmentMap.get(rightmostGroup);//the group which will be modified
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
					if (children.get(1).getLocalName().equals(ROOT_EL) && !frag.hasLocant(locant) && matchNumericLocant.matcher(locant).matches()){ //e.g. 1,3-benzothiazole-2-carbaldehyde 2-phenylhydrazone
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
			List<Element> functionalClass =XOMTools.getDescendantElementsWithTagName(words.get(1), GROUP_EL);
			if (functionalClass.size()!=1){
				throw new StructureBuildingException("Expected 1 group element found: " + functionalClass.size());
			}
			String smilesReplacement = functionalClass.get(0).getAttributeValue(VALUE_ATR);
			String labels =  functionalClass.get(0).getAttributeValue(LABELS_ATR);
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
				String[] locants = matchComma.split(StringTools.removeDashIfPresent(locantEls.get(0).getValue()));
                locantForFunctionalTerm.addAll(Arrays.asList(locants));
				locantEls.get(0).detach();
			}
		}
		if (!locantForFunctionalTerm.isEmpty() && locantForFunctionalTerm.size()!=replacementFragments.size()){
			throw new StructureBuildingException("Mismatch between number of locants and number of carbonyl replacements");
		}
		List<Atom> carbonylOxygens = findCarbonylOxygens(rootFragment, locantForFunctionalTerm);
		int functionalReplacementsToPerform = Math.min(replacementFragments.size(), carbonylOxygens.size());
		replaceCarbonylOxygenWithReplacementFragments(state, words, replacementFragments, carbonylOxygens, functionalReplacementsToPerform);
		resolveWordOrBracket(state, words.get(0));//the group
		
		if (!replacementFragments.isEmpty()){//Look for any more carbonyls that have appeared due to substitution e.g. 4-oxocyclohexa-2,5-diene-1-carboxylic acid 4-oxime
			BuildResults br = new BuildResults(state, words.get(0));
			for (Fragment f : br.getFragments()) {
				carbonylOxygens.addAll(findCarbonylOxygens(f, locantForFunctionalTerm));
			}
			if (carbonylOxygens.size() < replacementFragments.size()){
				throw new StructureBuildingException("Insufficient carbonyl groups found!");
			}
			functionalReplacementsToPerform = replacementFragments.size();
			replaceCarbonylOxygenWithReplacementFragments(state, words, replacementFragments, carbonylOxygens, functionalReplacementsToPerform);
		}
	}

	private void replaceCarbonylOxygenWithReplacementFragments(BuildState state, List<Element> words, List<Fragment> replacementFragments, List<Atom> carbonylOxygens, int functionalReplacementsToPerform) throws StructureBuildingException {
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
			state.fragManager.incorporateFragment(replacementFrag, carbonylFrag);
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
						Bond b = fragment.findBondOrThrow(atom, neighbours.get(0));
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
		List<Element> functionalClass =XOMTools.getDescendantElementsWithTagName(anhydrideWord, GROUP_EL);
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
		List<Element> anhydrideLocants =XOMTools.getDescendantElementsWithTagName(anhydrideWord, ANHYDRIDELOCANT_EL);
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
					String[] acidLocants =matchColon.split(StringTools.removeDashIfPresent(anhydrideLocant));
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
						String[] locants = matchComma.split(acidLocants[i]);
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
		state.fragManager.replaceTerminalAtomWithFragment(oxygen1, anhydride.getFirstAtom());
		List<Atom> atomsInAnhydrideLinkage = anhydride.getAtomList();
		if (atomsInAnhydrideLinkage.size()==0){//e.g. anhydride
			state.fragManager.createBond(oxygen1, atomOnSecondAcidToConnectTo, 1);
		}
		else{//e.g. peroxyanhydride
			state.fragManager.createBond(atomsInAnhydrideLinkage.get(atomsInAnhydrideLinkage.size()-1), atomOnSecondAcidToConnectTo, 1);
		}
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
				ideAtom.setCharge(ideAtom.getCharge()+1);//e.g. make cyanide charge netural
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
			state.fragManager.replaceTerminalAtomWithFragment(acidAtom, ideAtom);
		}
	}
	
	private void buildAdditionCompound(BuildState state, List<Element> words) throws StructureBuildingException {
		if (!words.get(0).getAttributeValue(TYPE_ATR).equals(WordType.full.toString())){
			throw new StructureBuildingException("Don't alter wordRules.xml without checking the consequences!");
		}
		resolveWordOrBracket(state, words.get(0));
		Atom elementaryAtom = state.xmlFragmentMap.get(StructureBuildingMethods.findRightMostGroupInBracket(words.get(0))).getFirstAtom();

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
				ideAtom.setCharge(ideAtom.getCharge()+1);//e.g. make cyanide charge netural
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
		int halideCount = functionalGroupFragments.size();
		Integer maximumVal = ValencyChecker.getMaximumValency(elementaryAtom.getElement(), elementaryAtom.getCharge());
		if (maximumVal!=null && halideCount > maximumVal){
			throw new StructureBuildingException("Too many halides/psuedo halides addded to " +elementaryAtom.getElement());
		}
		for (int i = halideCount - 1; i>=0; i--) {
			Fragment ideFrag =functionalGroupFragments.get(i);
			Atom ideAtom = ideFrag.getDefaultInAtom();
			state.fragManager.createBond(elementaryAtom, ideAtom, 1);
		}
	}

	private void buildAcetal(BuildState state, List<Element> words) throws StructureBuildingException {
		throw new StructureBuildingException("Not yet Supported");
		//TODO Add support for acetal/ketal/hemiketal/hemiacetal
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
			//look fo exact locant match
			Atom possibleAtom = mainGroupBR.getFunctionalAtom(i);
			if (possibleAtom.hasLocant(locant)){
				mainGroupBR.removeFunctionalAtom(i);
				if (possibleAtom.getProperty(Atom.AMBIGUOUS_ELEMENT_ASSIGNMENT)!=null){
					possibleAtom.getProperty(Atom.AMBIGUOUS_ELEMENT_ASSIGNMENT).remove(possibleAtom);
				}
				return possibleAtom;
			}
		}
		if (matchDigits.matcher(locant).matches()){
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
		else if (matchElementSymbolLocant.matcher(locant).matches()){
			//None of the functional atoms had an appropriate locant. Look for the special cases:
			//	Where the lack of primes on an element symbol locant should be ignored e.g. O,O-diethyl carbonate
			//	Where the locant is used to decide on the ester configuration c.f. O-methyl ..thioate and S-methyl ..thioate
			boolean isElementSymbol = matchElementSymbol.matcher(locant).matches();
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
	private void makeHydrogensExplicit(BuildState state) throws StructureBuildingException {
		Set<Fragment> fragments = state.fragManager.getFragPile();
		for (Fragment fragment : fragments) {
			if (fragment.getSubType().equals(ELEMENTARYATOMINORGANIC_SUBTYPE_VAL) || fragment.getSubType().equals(ELEMENTARYATOMORGANIC_SUBTYPE_VAL)){//these do not have implicit hydrogen e.g. phosphorus is literally just a phosphorus atom
				continue;
			}
			List<Atom> atomList =fragment.getAtomList();
			for (Atom parentAtom : atomList) {
				int explicitHydrogensToAdd = StructureBuildingMethods.calculateSubstitutableHydrogenAtoms(parentAtom);
				Atom hydrogen = null;
				for (int i = 0; i < explicitHydrogensToAdd; i++) {
					hydrogen = state.fragManager.createAtom("H", fragment);
					state.fragManager.createBond(parentAtom, hydrogen, 1);
				}
				parentAtom.setHydrogenCount(explicitHydrogensToAdd);
				if (parentAtom.getAtomParity()!=null){
					if (explicitHydrogensToAdd >1){
						throw new StructureBuildingException("Cannot have tetrahedral chirality and more than 2 hydrogens");
					}
					AtomParity atomParity = parentAtom.getAtomParity();
					Atom[] atomRefs4 = atomParity.getAtomRefs4();
					if (explicitHydrogensToAdd ==1){
						//atom parity was set in SMILES but at this stage the id of the hydrogen was not known, now it is so replace the dummy atom
						boolean setAtom =false;
						for (int i = 0; i < atomRefs4.length; i++) {
							if (atomRefs4[i].equals(AtomParity.hydrogen)){
								atomRefs4[i] = hydrogen;
								setAtom =true;
								break;
							}
						}
						if (!setAtom){
							throw new StructureBuildingException("OPSIN Bug: Unable to find implicit hydrogen to replace with explicit hydrogen in atomRefs4");
						}
					}
					else{
						//atom parity was set in SMILES, the dummy hydrogen atom has now been substituted
						List<Atom> neighbours = parentAtom.getAtomNeighbours();
						for (Atom atom : atomRefs4) {
							neighbours.remove(atom);
						}
						if (neighbours.size()!=1){
							throw new StructureBuildingException("OPSIN Bug: Unable to determine which atom has substitued a hydrogen at stereocentre");
						}
						boolean setAtom =false;
						for (int i = 0; i < atomRefs4.length; i++) {
							if (atomRefs4[i].equals(AtomParity.hydrogen)){
								atomRefs4[i] = neighbours.get(0);
								setAtom =true;
								break;
							}
						}
						if (!setAtom){
							throw new StructureBuildingException("OPSIN Bug: Unable to find implicit hydrogen to replace with explicit hydrogen in atomRefs4");
						}
					}
				}
			}
		}
	}

	/**
	 * A net charge is present; Given the molecule element the overallCharge is there an unambiguous way of 
	 * multiplying fragments to make the net charge 0
	 * metals without specified charge may be given an implicit positive charge
	 * 
	 * If this fails look for the case where there are multiple molecules and the mixture is only negative due to negatively charged functional Atoms e.g. pyridine acetate and remove the negative charge
	 * @param state
	 * @param molecule
	 * @param overallCharge
	 * @throws StructureBuildingException
	 */
	private void balanceChargeIfPossible(BuildState state, Element molecule, int overallCharge) throws StructureBuildingException {
		List<Element> words = XOMTools.getDescendantElementsWithTagNameAndAttribute(molecule, WORD_EL, TYPE_ATR, WordType.full.toString());
		if (words.size() <2){
			return;
		}
		List<Element> positivelyChargedWords = new ArrayList<Element>();
		List<Element> negativelyChargedWords = new ArrayList<Element>();
		HashMap<Element, Integer> wordToChargeMapping = new HashMap<Element, Integer>();
		HashMap<Element, BuildResults> wordToBR = new HashMap<Element, BuildResults>();
		
		List<Element> cationicElements = new ArrayList<Element>();
		List<Element> inorganicElementaryAtoms = XOMTools.getDescendantElementsWithTagNameAndAttribute(molecule, GROUP_EL, SUBTYPE_ATR, ELEMENTARYATOMINORGANIC_SUBTYPE_VAL);
		for (Element elementaryAtom : inorganicElementaryAtoms) {
			if (elementaryAtom.getAttribute(TYPICALANDMAXIMUMCHARGE_ATR)!=null){
				Fragment cationicFrag =state.xmlFragmentMap.get(elementaryAtom);
				int typicalCharge = Integer.parseInt(matchComma.split(elementaryAtom.getAttributeValue(TYPICALANDMAXIMUMCHARGE_ATR))[0]);
				if (typicalCharge > cationicFrag.getFirstAtom().getAtomNeighbours().size()){
					cationicElements.add(elementaryAtom);
				}
			}
		}
		overallCharge = setCationicElementsToTypicalCharge(state, cationicElements, overallCharge);
		if (overallCharge==0){
			return;
		}
		for (Element word : words) {
			BuildResults br = new BuildResults(state, word);
			wordToBR.put(word, br);
			int charge = br.getCharge();
			if (charge>0){
				positivelyChargedWords.add(word);
			}
			else if (charge <0){
				negativelyChargedWords.add(word);
			}
			wordToChargeMapping.put(word, charge);
		}
		if (positivelyChargedWords.size()==1 && cationicElements.size() ==0 && negativelyChargedWords.size() >=1 || positivelyChargedWords.size()>=1 && negativelyChargedWords.size() ==1 ){
			boolean success = multiplyChargedWords(state, negativelyChargedWords, positivelyChargedWords, wordToChargeMapping, overallCharge);
			if (success){
				return;
			}
		}
		if (cationicElements.size() ==1){
			boolean success = setChargeOnCationicElementAppropriately(state, overallCharge, cationicElements.get(0));
			if (success){
				return;
			}
		}
		if (overallCharge <0){//neutralise functionalAtoms if they are the sole cause of the negative charge and multiple molecules are present
			int chargeOnFunctionalAtoms = 0;
			for (Element word : words) {
				BuildResults br = wordToBR.get(word);
				int functionalAtomCount = br.getFunctionalAtomCount();
				for (int i = functionalAtomCount -1; i >=0; i--) {
					chargeOnFunctionalAtoms += br.getFunctionalAtom(i).getCharge();
				}
			}
			if (chargeOnFunctionalAtoms == overallCharge){
				for (Element word : words) {
					BuildResults br = wordToBR.get(word);
					int functionalAtomCount = br.getFunctionalAtomCount();
					for (int i = functionalAtomCount -1; i >=0; i--) {
						br.getFunctionalAtom(i).setCharge(0);
						br.removeFunctionalAtom(i);
					}
				}
			}
		}
	}


	/**
	 * Sets the cationicElements to their typical charge as specified by the TYPICALANDMAXIMUMCHARGE_ATR
	 * The number of atoms the cationicElement is connected to is taken into account e.g. phenylmagnesium chloride is [Mg+]
	 * @param state
	 * @param cationicElements
	 * @param overallCharge
	 * @return
	 * @throws StructureBuildingException
	 */
	private int setCationicElementsToTypicalCharge(BuildState state, List<Element> cationicElements, int overallCharge) throws StructureBuildingException {
		for (Element cationicElement : cationicElements) {
			Fragment cationicFrag = state.xmlFragmentMap.get(cationicElement);
			int charge = Integer.parseInt(matchComma.split(cationicElement.getAttributeValue(TYPICALANDMAXIMUMCHARGE_ATR))[0]);
			charge -= cationicFrag.getFirstAtom().getAtomNeighbours().size();
			overallCharge += charge;
			cationicFrag.getFirstAtom().setCharge(charge);
		}
		return overallCharge;
	}

	/**
	 * Multiplies out charged Words to balance charge
	 * Return true if balancing was possible else false
	 * @param state
	 * @param negativelyChargedWords
	 * @param positivelyChargedWords
	 * @param wordToChargeMapping
	 * @param overallCharge
	 * @return
	 * @throws StructureBuildingException
	 */
	private boolean multiplyChargedWords(BuildState state, List<Element>negativelyChargedWords,List<Element> positivelyChargedWords,HashMap<Element, Integer> wordToChargeMapping, int overallCharge) throws StructureBuildingException {
		Element wordToMultiply;
		if (overallCharge >0){
			if (negativelyChargedWords.size() >1){
				return false;//ambiguous as to which to multiply
			}
			wordToMultiply = negativelyChargedWords.get(0);
		}
		else{
			if (positivelyChargedWords.size() >1){
				return false;//ambiguous as to which to multiply
			}
			wordToMultiply = positivelyChargedWords.get(0);
		}
		Element firstChild = (Element) wordToMultiply.getChild(0);
		while (firstChild.getChildElements().size() !=0){
			firstChild = (Element) firstChild.getChild(0);
		}
		if (firstChild.getLocalName().equals(MULTIPLIER_EL)){//e.g. monochloride. Allows specification of explicit stoichiometry
			return false;
		}
		int charge = wordToChargeMapping.get(wordToMultiply);
		if (overallCharge % charge ==0){
			int timesToDuplicate = Math.abs(overallCharge/charge);
			for (int i = 0; i < timesToDuplicate; i++) {
				XOMTools.insertAfter(wordToMultiply, state.fragManager.cloneElement(state, wordToMultiply));
			}
		}
		return true;
	}
	
	private boolean setChargeOnCationicElementAppropriately(BuildState state, int overallCharge, Element cationicElement) throws StructureBuildingException {
		Atom cation = state.xmlFragmentMap.get(cationicElement).getFirstAtom();
		int chargeOnCationNeeded = -(overallCharge -cation.getCharge());
		int maximumCharge = Integer.parseInt(matchComma.split(cationicElement.getAttributeValue(TYPICALANDMAXIMUMCHARGE_ATR))[1]);
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
}
