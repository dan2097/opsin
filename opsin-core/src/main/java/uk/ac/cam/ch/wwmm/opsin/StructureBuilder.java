package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;

import uk.ac.cam.ch.wwmm.opsin.StereoAnalyser.StereoBond;
import uk.ac.cam.ch.wwmm.opsin.StereoAnalyser.StereoCentre;
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
	private final BuildState state;
	private final List<Atom> polymerAttachmentPoints = new ArrayList<Atom>();//rGroups need to be represented as normal atoms for the purpose of working out stereochemistry. They will be converted to a suitable representation later
	
	private int currentTopLevelWordRuleCount;
	
	StructureBuilder(BuildState state) {
		this.state = state;
	}

	/**	Builds a molecule as a Fragment based on ComponentProcessor output.
	 * @param molecule The ComponentProcessor output.
	 * @return A single Fragment - the built molecule.
	 * @throws StructureBuildingException If the molecule won't build - there may be many reasons.
	 */
	Fragment buildFragment(Element molecule) throws StructureBuildingException {
		List<Element> wordRules = molecule.getChildElements(WORDRULE_EL);

		currentTopLevelWordRuleCount = wordRules.size();
		if (currentTopLevelWordRuleCount == 0) {
			throw new StructureBuildingException("Molecule contains no word rules!?");
		}
		
		for (Element wordRule : wordRules) {
			processWordRuleChildrenThenRule(wordRule);
		}
		
		if (currentTopLevelWordRuleCount != wordRules.size()) {
			wordRules = molecule.getChildElements(WORDRULE_EL);//very rarely a word rule adds a top level word rule
		}

		List<Element> groupElements = OpsinTools.getDescendantElementsWithTagName(molecule, GROUP_EL);
		processSpecialCases(groupElements);
		processOxidationNumbers(groupElements);
		state.fragManager.convertSpareValenciesToDoubleBonds();
		state.fragManager.checkValencies();
		
		manipulateStoichiometry(molecule, wordRules);
		
		state.fragManager.makeHydrogensExplicit();

		Fragment uniFrag = state.fragManager.getUnifiedFragment();
		processStereochemistry(molecule, uniFrag);

		if (uniFrag.getOutAtomCount() > 0) {
			if (!state.n2sConfig.isAllowRadicals()) {
				throw new StructureBuildingException("Radicals are currently set to not convert to structures");
			}
			if (state.n2sConfig.isOutputRadicalsAsWildCardAtoms()) {
				convertOutAtomsToAttachmentAtoms(uniFrag);
			}
		}
		
		if (polymerAttachmentPoints.size() > 0) {
			for (Atom rAtom : polymerAttachmentPoints) {
				rAtom.setElement(ChemEl.R);
			}
			uniFrag.setPolymerAttachmentPoints(polymerAttachmentPoints);
		}
		return uniFrag;
	}


	private void processWordRuleChildrenThenRule(Element wordRule) throws StructureBuildingException {
		List<Element> wordRuleChildren = wordRule.getChildElements(WORDRULE_EL);
		for (Element wordRuleChild : wordRuleChildren) {
			processWordRuleChildrenThenRule(wordRuleChild);
		}
		processWordRule(wordRule);
	}
	
	private void processWordRule(Element wordRuleEl) throws StructureBuildingException {
		WordRule wordRule = WordRule.valueOf(wordRuleEl.getAttributeValue(WORDRULE_ATR));
		List<Element> words = OpsinTools.getChildElementsWithTagNames(wordRuleEl, new String[]{WORD_EL, WORDRULE_EL});
		state.currentWordRule = wordRule;
		switch (wordRule) {
		case simple:
			for (Element word : words) {
				if (!word.getName().equals(WORD_EL) || !word.getAttributeValue(TYPE_ATR).equals(WordType.full.toString())){
					throw new StructureBuildingException("OPSIN bug: Unexpected contents of 'simple' wordRule");
				}
				resolveWordOrBracket(state, word);
			}
			break;
		case substituent:
			for (Element word : words) {
				if (!word.getName().equals(WORD_EL) || !word.getAttributeValue(TYPE_ATR).equals(WordType.substituent.toString()) || !state.n2sConfig.isAllowRadicals()){
					throw new StructureBuildingException("OPSIN bug: Unexpected contents of 'substituent' wordRule");
				}
				resolveWordOrBracket(state, word);
			}
			break;
		case ester:
		case multiEster:
			buildEster(words);//e.g. ethyl ethanoate, dimethyl terephthalate,  methyl propanamide
			break;
		case divalentFunctionalGroup:
			buildDiValentFunctionalGroup(words);// diethyl ether or methyl propyl ketone
			break;
		case monovalentFunctionalGroup:
			buildMonovalentFunctionalGroup(words);// ethyl chloride, isophthaloyl dichloride, diethyl ether, ethyl alcohol
			break;
		case functionalClassEster:
			buildFunctionalClassEster(words);//e.g. ethanoic acid ethyl ester, tetrathioterephthalic acid dimethyl ester
			break;
		case acidReplacingFunctionalGroup:
			//e.g. ethanoic acid ethyl amide, terephthalic acid dimethyl amide,
			//ethanoic acid amide, carbonic dihydrazide
			//already processed by the ComponentProcessor
			for (Element word : words) {
				resolveWordOrBracket(state, word);
			}
			break;
		case oxide:
			buildOxide(words);//e.g. styrene oxide, triphenylphosphane oxide, thianthrene 5,5-dioxide, propan-2-one oxide
			break;
		case carbonylDerivative:
			buildCarbonylDerivative(words);//e.g. Imidazole-2-carboxamide O-ethyloxime, pentan-3-one oxime
			break;
		case anhydride:
			buildAnhydride(words);//e.g. acetic anhydride
			break;
		case acidHalideOrPseudoHalide:
			buildAcidHalideOrPseudoHalide(words);//e.g. phosphinimidic chloride
			break;
		case additionCompound:
			buildAdditionCompound(words);//e.g. carbon tetrachloride
			break;
		case glycol:
			buildGlycol(words);//e.g. ethylene glycol
			break;
		case glycolEther:
			buildGlycolEther(words);//e.g. octaethyleneglycol monododecyl ether
			break;
		case acetal:
			buildAcetal(words);//e.g. propanal diethyl acetal
			break;
		case potentialAlcoholEster:
			//e.g. uridine 5'-(tetrahydrogen triphosphate)
			if (!buildAlcoholEster(words, currentTopLevelWordRuleCount)){
				//should be processed as two "simple" wordrules if no hydroxy found, hence number of top level word rules may change
				//These simple word rules have already been processed
				splitAlcoholEsterRuleIntoTwoSimpleWordRules(words);
				currentTopLevelWordRuleCount++;
			}
			break;
		case cyclicPeptide:
			buildCyclicPeptide(words);
			break;
		case amineDiConjunctiveSuffix:
			//e.g. glycine N,N-diacetic acid
			buildAmineDiConjunctiveSuffix(words);
			break;
		case polymer:
			buildPolymer(words);
			break;
		default:
			throw new StructureBuildingException("Unexpected Word Rule");
		}
	}


	private void buildEster(List<Element> words) throws StructureBuildingException {
		boolean inSubstituents = true;
		BuildResults substituentsBr = new BuildResults();
		List<BuildResults> ateGroups = new ArrayList<BuildResults>();
		Map<BuildResults, String> buildResultsToLocant = new HashMap<BuildResults, String>();//typically locant will be null
		
		for (Element word : words) {
			resolveWordOrBracket(state, word);
			BuildResults br = new BuildResults(word);
			if (inSubstituents && br.getFunctionalAtomCount() > 0){
				inSubstituents = false;
			}
			if (inSubstituents){
				if (!word.getAttributeValue(TYPE_ATR).equals(WordType.substituent.toString())){
					if (word.getAttributeValue(TYPE_ATR).equals(WordType.full.toString())){
						throw new StructureBuildingException("bug? ate group did not have any functional atoms!");
					}
					else{
						throw new StructureBuildingException("OPSIN bug: Non substituent word found where substituent expected in ester");
					}
				}
				int outAtomCount = br.getOutAtomCount();
				boolean traditionalEster =false;
				for (int i = 0; i < outAtomCount; i++) {
					OutAtom out = br.getOutAtom(i);
					if (out.getValency()>1){
						FragmentTools.splitOutAtomIntoValency1OutAtoms(out);
						traditionalEster =true;
					}
				}
				if (traditionalEster){//e.g. ethylidene dipropanoate
					br = new BuildResults(word);
					outAtomCount = br.getOutAtomCount();
				}
				if (outAtomCount ==1){//TODO add support for locanted terepthaloyl
					String locantForSubstituent = word.getAttributeValue(LOCANT_ATR);
					if (locantForSubstituent!=null){
						br.getFirstOutAtom().setLocant(locantForSubstituent);//indexes which functional atom to connect to when there is a choice. Also can disambiguate which atom is a S in things like thioates
					}
				}
				else if (outAtomCount ==0){
					throw new StructureBuildingException("Substituent was expected to have at least one outAtom");
				}
				substituentsBr.mergeBuildResults(br);
			}
			else{
				String locant = word.getAttributeValue(LOCANT_ATR);//specifying a locant for an ateWord is very unusual as this information is typically redundant c.f. dodecamethylene 1,12-bis(chloroformate)
				if (br.getFunctionalAtomCount()<1){
					throw new StructureBuildingException("bug? ate group did not have any functional atoms!");
				}
				ateGroups.add(br);
				buildResultsToLocant.put(br, locant);
			}
		}
		if (ateGroups.size() ==0){
			throw new StructureBuildingException("OPSIN bug: Missing ate group in ester");
		}

		int outAtomCount =substituentsBr.getOutAtomCount();
		if (outAtomCount ==0){
			throw new StructureBuildingException("OPSIN bug: Missing outatom on ester substituents");
		}
		int esterIdCount = 0;
		for (BuildResults br : ateGroups) {
			esterIdCount += br.getFunctionalAtomCount();
		}
		if (outAtomCount > esterIdCount){
			throw new StructureBuildingException("There are more radicals in the substituents(" + outAtomCount +") than there are places to form esters("+esterIdCount+")");
		}
		if (esterIdCount > outAtomCount && outAtomCount % ateGroups.size() !=0) {
			//actually checks if the same number of ester forming points would be used in each ate group e.g. ethyl diacetate is wrong
			throw new StructureBuildingException("There are less radicals in the substituents(" + outAtomCount +") than there are places to form esters("+esterIdCount+")");
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
				Atom atomOnSubstituentToUse = getOutAtomTakingIntoAccountWhetherSetExplicitly(substituentsBr, 0);
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



	private void buildDiValentFunctionalGroup(List<Element> words) throws StructureBuildingException {
		int wordIndice = 0;
		if (!words.get(wordIndice).getAttributeValue(TYPE_ATR).equals(WordType.substituent.toString())) {
			throw new StructureBuildingException("word: " +wordIndice +" was expected to be a substituent");
		}
		resolveWordOrBracket(state, words.get(wordIndice));
		BuildResults substituent1 =new BuildResults(words.get(wordIndice));
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
				OpsinTools.insertAfter(words.get(0), clone);
				words = words.get(0).getParent().getChildElements();
			}
			else{
				resolveWordOrBracket(state, words.get(wordIndice));
			}
			substituent2 =new BuildResults(words.get(wordIndice));
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
		List<Element> functionalGroup = OpsinTools.getDescendantElementsWithTagName(words.get(wordIndice), FUNCTIONALGROUP_EL);
		if (functionalGroup.size()!=1){
			throw new StructureBuildingException("Unexpected number of functionalGroups found, could be a bug in OPSIN's grammar");
		}
		String smilesOfGroup = functionalGroup.get(0).getAttributeValue(VALUE_ATR);
		Fragment diValentGroup =state.fragManager.buildSMILES(smilesOfGroup, FUNCTIONALCLASS_TYPE_VAL, NONE_LABELS_VAL);

		Atom outAtom1 = getOutAtomTakingIntoAccountWhetherSetExplicitly(substituent1, 0);
		substituent1.removeOutAtom(0);
		Atom outAtom2 = getOutAtomTakingIntoAccountWhetherSetExplicitly(substituent2, 0);
		substituent2.removeOutAtom(0);
		if (diValentGroup.getOutAtomCount()==1){//c.f. peroxide where it is a linker
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

	private void buildMonovalentFunctionalGroup(List<Element> words) throws StructureBuildingException {
		resolveWordOrBracket(state, words.get(0));
		List<Element> groups = OpsinTools.getDescendantElementsWithTagName(words.get(0), GROUP_EL);
		for (Element group : groups) {//replaces outAtoms with valency greater than 1 with multiple outAtoms; e.g. ylidene -->diyl
			Fragment frag = group.getFrag();
			for (int i = frag.getOutAtomCount()-1; i>=0; i--) {
				OutAtom outAtom =frag.getOutAtom(i);
				if (outAtom.getValency()>1){
					FragmentTools.splitOutAtomIntoValency1OutAtoms(outAtom);
				}
			}
		}
		BuildResults substituentBR = new BuildResults(words.get(0));

		List<Fragment> functionalGroupFragments = new ArrayList<Fragment>();
		for (int i=1; i<words.size(); i++ ) {
			Element functionalGroupWord =words.get(i);
			List<Element> functionalGroups = OpsinTools.getDescendantElementsWithTagName(functionalGroupWord, FUNCTIONALGROUP_EL);
			if (functionalGroups.size()!=1){
				throw new StructureBuildingException("Expected exactly 1 functionalGroup. Found " + functionalGroups.size());
			}
			
			Fragment monoValentFunctionGroup =state.fragManager.buildSMILES(functionalGroups.get(0).getAttributeValue(VALUE_ATR), FUNCTIONALCLASS_TYPE_VAL, NONE_LABELS_VAL);
			if (functionalGroups.get(0).getAttributeValue(TYPE_ATR).equals(MONOVALENTSTANDALONEGROUP_TYPE_VAL)){
				Atom ideAtom = monoValentFunctionGroup.getDefaultInAtomOrFirstAtom();
				ideAtom.addChargeAndProtons(1, 1);//e.g. make cyanide charge netural
			}
			Element possibleMultiplier = OpsinTools.getPreviousSibling(functionalGroups.get(0));
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
			Atom ideAtom = ideFrag.getDefaultInAtomOrFirstAtom();
			Atom subAtom = getOutAtomTakingIntoAccountWhetherSetExplicitly(substituentBR, 0);
			state.fragManager.createBond(ideAtom, subAtom, 1);
			substituentBR.removeOutAtom(0);
			state.fragManager.incorporateFragment(ideFrag, subAtom.getFrag());
		}
	}

	private void buildFunctionalClassEster(List<Element> words) throws StructureBuildingException {
		Element firstWord = words.get(0);
		if (!firstWord.getAttributeValue(TYPE_ATR).equals(WordType.full.toString())) {
			throw new StructureBuildingException("Don't alter wordRules.xml without checking the consequences!");
		}
		resolveWordOrBracket(state, firstWord);//the group
		BuildResults acidBr = new BuildResults(firstWord);

		if (acidBr.getFunctionalAtomCount()==0) {
			throw new StructureBuildingException("No functionalAtoms detected!");
		}

		int wordCountMinus1 = words.size() - 1;
		if (wordCountMinus1 < 2 || !words.get(wordCountMinus1).getAttributeValue(TYPE_ATR).equals(WordType.functionalTerm.toString())) {
			throw new StructureBuildingException("OPSIN Bug: Bug in functionalClassEster rule; 'ester' not found where it was expected");
		}
		
		for (int i = 1; i < wordCountMinus1; i++) {
			Element currentWord = words.get(i);
			String wordType  = currentWord.getAttributeValue(TYPE_ATR);
			if (!wordType.equals(WordType.substituent.toString())) {
				if (wordType.equals(WordType.functionalTerm.toString()) && currentWord.getAttributeValue(VALUE_ATR).equalsIgnoreCase("ester")) {
					//superfluous ester word
					continue;
				}
				throw new StructureBuildingException("OPSIN Bug: Bug in functionalClassEster rule; Encountered: " + currentWord.getAttributeValue(VALUE_ATR));
			}
			resolveWordOrBracket(state, currentWord);
			BuildResults substituentBr = new BuildResults(currentWord);
			int outAtomCount = substituentBr.getOutAtomCount();
			if (acidBr.getFunctionalAtomCount() < outAtomCount) {
				throw new StructureBuildingException("Insufficient functionalAtoms on acid");
			}
			for (int j = 0; j < outAtomCount; j++) {
				String locantForSubstituent = currentWord.getAttributeValue(LOCANT_ATR);
				Atom functionalAtom;
				if (locantForSubstituent != null) {
					functionalAtom = determineFunctionalAtomToUse(locantForSubstituent, acidBr);
				}
				else{
					functionalAtom = acidBr.getFunctionalAtom(0);
					acidBr.removeFunctionalAtom(0);
				}
				if (substituentBr.getOutAtom(j).getValency() != 1) {
					throw new StructureBuildingException("Substituent was expected to have only have an outgoing valency of 1");
				}
				state.fragManager.createBond(functionalAtom, getOutAtomTakingIntoAccountWhetherSetExplicitly(substituentBr, j), 1);
				if (functionalAtom.getCharge() == -1) {
					functionalAtom.neutraliseCharge();
				}
			}
			substituentBr.removeAllOutAtoms();
		}
	}
	
	/**
	 * Handles names like thiophene 1,1-dioxide; carbon dioxide; benzene oxide
	 * Does the same for sulfide/selenide/telluride
	 * @param words
	 * @throws StructureBuildingException
	 */
	private void buildOxide(List<Element> words) throws StructureBuildingException {
		resolveWordOrBracket(state, words.get(0));//the group
		List<Fragment> oxideFragments = new ArrayList<Fragment>();
		List<String> locantsForOxide =new ArrayList<String>();//often not specified
		if (!words.get(1).getAttributeValue(TYPE_ATR).equals(WordType.functionalTerm.toString())){
			throw new StructureBuildingException("Oxide functional term not found where expected!");
		}
		Element rightMostGroup;
		if (words.get(0).getName().equals(WORDRULE_EL)){//e.g. Nicotinic acid N-oxide
			List<Element> fullWords = OpsinTools.getDescendantElementsWithTagNameAndAttribute(words.get(0), WORD_EL, TYPE_ATR, WordType.full.toString());
			if (fullWords.size()==0){
				throw new StructureBuildingException("OPSIN is entirely unsure where the oxide goes so has decided not to guess");
			}
			rightMostGroup = findRightMostGroupInBracket(fullWords.get(fullWords.size()-1));
		}
		else{
			rightMostGroup = findRightMostGroupInBracket(words.get(0));
		}
		
		int numberOfOxygenToAdd =1;
		List<Element> multipliers =OpsinTools.getDescendantElementsWithTagName(words.get(1), MULTIPLIER_EL);
		if (multipliers.size() >1){
			throw new StructureBuildingException("Expected 0 or 1 multiplier found: " + multipliers.size());
		}
		if (multipliers.size()==1){
			numberOfOxygenToAdd = Integer.parseInt(multipliers.get(0).getAttributeValue(VALUE_ATR));
			multipliers.get(0).detach();
		}
		else{
			if (ELEMENTARYATOM_TYPE_VAL.equals(rightMostGroup.getAttributeValue(TYPE_ATR))){
				Atom elementaryAtom = rightMostGroup.getFrag().getFirstAtom();
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
		List<Element> functionalGroup =OpsinTools.getDescendantElementsWithTagName(words.get(1), FUNCTIONALGROUP_EL);
		if (functionalGroup.size()!=1){
			throw new StructureBuildingException("Expected 1 group element found: " + functionalGroup.size());
		}
		String smilesReplacement = functionalGroup.get(0).getAttributeValue(VALUE_ATR);
		String labels =  functionalGroup.get(0).getAttributeValue(LABELS_ATR);
		for (int i = 0; i < numberOfOxygenToAdd; i++) {
			oxideFragments.add(state.fragManager.buildSMILES(smilesReplacement, FUNCTIONALCLASS_TYPE_VAL, labels != null ? labels : NONE_LABELS_VAL));
		}
		List<Element> locantEls =OpsinTools.getDescendantElementsWithTagName(words.get(1), LOCANT_EL);
		if (locantEls.size() >1){
			throw new StructureBuildingException("Expected 0 or 1 locant elements found: " + locantEls.size());
		}
		if (locantEls.size()==1){
			String[] locants = StringTools.removeDashIfPresent(locantEls.get(0).getValue()).split(",");
			locantsForOxide.addAll(Arrays.asList(locants));
			locantEls.get(0).detach();
		}
		if (!locantsForOxide.isEmpty() && locantsForOxide.size()!=oxideFragments.size()){
			throw new StructureBuildingException("Mismatch between number of locants and number of oxides specified");
		}	
		Fragment groupToModify = rightMostGroup.getFrag();//all the suffixes are part of this fragment at this point
		mainLoop: for (int i = 0; i < oxideFragments.size(); i++) {
			Atom oxideAtom = oxideFragments.get(i).getFirstAtom();
			if (!locantsForOxide.isEmpty()){
				Atom atomToAddOxideTo =groupToModify.getAtomByLocantOrThrow(locantsForOxide.get(i));
				if (atomToAddOxideTo.getElement() == ChemEl.C && !ELEMENTARYATOM_TYPE_VAL.equals(groupToModify.getType())) {
					throw new StructureBuildingException("Locant " + locantsForOxide.get(i) + " indicated oxide applied to carbon, but this would lead to hypervalency!");
				}
				formAppropriateBondToOxideAndAdjustCharges(atomToAddOxideTo, oxideAtom);
			}
			else{
				if (ELEMENTARYATOM_TYPE_VAL.equals(groupToModify.getType())){
					Atom elementaryAtom= groupToModify.getFirstAtom();
					formAppropriateBondToOxideAndAdjustCharges(elementaryAtom, oxideAtom);//e.g. carbon dioxide
					int chargeOnAtom =elementaryAtom.getCharge();
					if (chargeOnAtom>=2){
						elementaryAtom.setCharge(chargeOnAtom-2);
					}
					continue mainLoop;
				}
				else{
					List<Atom> atomList = groupToModify.getAtomList();
					//In preference suffixes are substituted onto e.g. acetonitrile oxide
					for (Atom atom : atomList) {
						if (!atom.getType().equals(SUFFIX_TYPE_VAL)) {
							continue;
						}
						if (atom.getElement() != ChemEl.C && atom.getElement() != ChemEl.O) {
							formAppropriateBondToOxideAndAdjustCharges(atom, oxideAtom);
							continue mainLoop;
						}
					}
					for (Atom atom : atomList) {
						if (atom.getElement() != ChemEl.C && atom.getElement() != ChemEl.O) {
							formAppropriateBondToOxideAndAdjustCharges(atom, oxideAtom);
							continue mainLoop;
						}
					}
				}
				//No heteroatoms could be found. Perhaps it's supposed to be something like styrene oxide
				Set<Bond> bondSet = groupToModify.getBondSet();//looking for double bond
				for (Bond bond : bondSet) {
					if (bond.getOrder()==2 && bond.getFromAtom().getElement() == ChemEl.C && bond.getToAtom().getElement() == ChemEl.C){
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
					if (fromAtom.hasSpareValency() && toAtom.hasSpareValency() &&fromAtom.getElement() == ChemEl.C && toAtom.getElement() == ChemEl.C){
						fromAtom.setSpareValency(false);
						toAtom.setSpareValency(false);
						state.fragManager.createBond(fromAtom, oxideAtom, 1);
						state.fragManager.createBond(toAtom, oxideAtom, 1);
						continue mainLoop;
					}
				}
				//something like where oxide goes on an oxygen propan-2-one oxide
				List<Atom> atomList = groupToModify.getAtomList();
				for (Atom atom : atomList) {
					if (!atom.getType().equals(SUFFIX_TYPE_VAL)) {
						continue;
					}
					if (atom.getElement() != ChemEl.C) {
						formAppropriateBondToOxideAndAdjustCharges(atom, oxideAtom);
						continue mainLoop;
					}
				}
				for (Atom atom : atomList) {
					if (atom.getElement() != ChemEl.C) {
						formAppropriateBondToOxideAndAdjustCharges(atom, oxideAtom);
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
	 * @param atomToAddOxideTo
	 * @param oxideAtom
	 * @throws StructureBuildingException 
	 */
	private void formAppropriateBondToOxideAndAdjustCharges(Atom atomToAddOxideTo, Atom oxideAtom) throws StructureBuildingException {
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

	private void buildCarbonylDerivative(List<Element> words) throws StructureBuildingException {
		if (!WordType.full.toString().equals(words.get(0).getAttributeValue(TYPE_ATR))){
			throw new StructureBuildingException("OPSIN bug: Wrong word type encountered when applying carbonylDerivative wordRule");
		}
		List<Fragment> replacementFragments = new ArrayList<Fragment>();
		List<String> locantForFunctionalTerm =new ArrayList<String>();//usually not specified
		if (!words.get(1).getAttributeValue(TYPE_ATR).equals(WordType.functionalTerm.toString())){//e.g. acetone O-ethyloxime or acetone 1-chloro-1-methylhydrazone
			for (int i = 1; i < words.size(); i++) {
				Fragment frag = findRightMostGroupInWordOrWordRule(words.get(i)).getFrag();
				replacementFragments.add(frag);
				List<Element> children =words.get(i).getChildElements();
				if (children.size()==1 && children.get(0).getName().equals(BRACKET_EL) && children.get(0).getAttribute(LOCANT_ATR)!=null){
					locantForFunctionalTerm.add(children.get(0).getAttributeValue(LOCANT_ATR));
				}
				else if (children.size()==2 && children.get(0).getAttribute(LOCANT_ATR)!=null ){
					String locant =children.get(0).getAttributeValue(LOCANT_ATR);
					if (children.get(1).getName().equals(ROOT_EL) && !frag.hasLocant(locant) && MATCH_NUMERIC_LOCANT.matcher(locant).matches()){ //e.g. 1,3-benzothiazole-2-carbaldehyde 2-phenylhydrazone
						locantForFunctionalTerm.add(children.get(0).getAttributeValue(LOCANT_ATR));
						children.get(0).removeAttribute(children.get(0).getAttribute(LOCANT_ATR));
					}
				}
			}
		}
		else{//e.g. butan-2,3-dione dioxime or hexan2,3-dione 2-oxime
			int numberOfCarbonylReplacements =1;
			List<Element> multipliers =OpsinTools.getDescendantElementsWithTagName(words.get(1), MULTIPLIER_EL);
			if (multipliers.size() >1){
				throw new StructureBuildingException("Expected 0 or 1 multiplier found: " + multipliers.size());
			}
			if (multipliers.size()==1){
				numberOfCarbonylReplacements = Integer.parseInt(multipliers.get(0).getAttributeValue(VALUE_ATR));
				multipliers.get(0).detach();
			}
			List<Element> functionalGroup =OpsinTools.getDescendantElementsWithTagName(words.get(1), FUNCTIONALGROUP_EL);
			if (functionalGroup.size()!=1){
				throw new StructureBuildingException("Expected 1 functionalGroup element found: " + functionalGroup.size());
			}
			String smilesReplacement = functionalGroup.get(0).getAttributeValue(VALUE_ATR);
			String labels =  functionalGroup.get(0).getAttributeValue(LABELS_ATR);
			for (int i = 0; i < numberOfCarbonylReplacements; i++) {
				Fragment replacementFragment = state.fragManager.buildSMILES(smilesReplacement, FUNCTIONALCLASS_TYPE_VAL, labels != null ? labels : NONE_LABELS_VAL);
				if (i >0){
					FragmentTools.relabelLocants(replacementFragment.getAtomList(), StringTools.multiplyString("'", i));
				}
				List<Atom> atomList = replacementFragment.getAtomList();
				for (Atom atom : atomList) {
					atom.removeLocantsOtherThanElementSymbolLocants();//prevents numeric locant locanted substitution from outside the functional word
				}
				replacementFragments.add(replacementFragment);
			}
			List<Element> locantEls =OpsinTools.getDescendantElementsWithTagName(words.get(1), LOCANT_EL);
			if (locantEls.size() >1){
				throw new StructureBuildingException("Expected 0 or 1 locant elements found: " + locantEls.size());
			}
			if (locantEls.size() == 1) {
				String[] locants = StringTools.removeDashIfPresent(locantEls.get(0).getValue()).split(",");
				locantForFunctionalTerm.addAll(Arrays.asList(locants));
				locantEls.get(0).detach();
			}
		}
		if (!locantForFunctionalTerm.isEmpty() && locantForFunctionalTerm.size()!=replacementFragments.size()){
			throw new StructureBuildingException("Mismatch between number of locants and number of carbonyl replacements");
		}

		Element rightMostGroup = findRightMostGroupInWordOrWordRule(words.get(0));
		Element parent = rightMostGroup.getParent();
		boolean multiplied =false;
		while (!parent.equals(words.get(0))){
			if (parent.getAttribute(MULTIPLIER_ATR)!=null){
				multiplied =true;
			}
			parent = parent.getParent();
		}
		if (!multiplied){
			List<Atom> carbonylOxygens = findCarbonylOxygens(rightMostGroup.getFrag(), locantForFunctionalTerm);
			int replacementsToPerform = Math.min(replacementFragments.size(), carbonylOxygens.size());
			replaceCarbonylOxygenWithReplacementFragments(words, replacementFragments, carbonylOxygens, replacementsToPerform);
		}

		resolveWordOrBracket(state, words.get(0));//the component
		if (replacementFragments.size() >0){
			//Note that the right most group may be multiplied e.g. 3,3'-methylenebis(2,4,6-trimethylbenzaldehyde) disemicarbazone
			//or the carbonyl may not even be on the right most group e.g.  4-oxocyclohexa-2,5-diene-1-carboxylic acid 4-oxime
			BuildResults br = new BuildResults(words.get(0));
			List<Atom> carbonylOxygens = new ArrayList<Atom>();
			List<Fragment> fragments = new ArrayList<Fragment>(br.getFragments());
			for (ListIterator<Fragment> iterator = fragments.listIterator(fragments.size()); iterator.hasPrevious();) {//iterate in reverse order - right most groups preferred
				carbonylOxygens.addAll(findCarbonylOxygens(iterator.previous(), locantForFunctionalTerm));
			}
			replaceCarbonylOxygenWithReplacementFragments(words, replacementFragments, carbonylOxygens, replacementFragments.size());
		}
	}

	private void replaceCarbonylOxygenWithReplacementFragments(List<Element> words, List<Fragment> replacementFragments, List<Atom> carbonylOxygens, int functionalReplacementsToPerform) throws StructureBuildingException {
		if (functionalReplacementsToPerform > carbonylOxygens.size()){
			throw new StructureBuildingException("Insufficient carbonyl groups found!");
		}
		for (int i = 0; i < functionalReplacementsToPerform; i++) {
			Atom carbonylOxygen =carbonylOxygens.remove(0);//the oxygen of the carbonyl
			Fragment carbonylFrag = carbonylOxygen.getFrag();
			Fragment replacementFrag = replacementFragments.remove(0);
			List<Atom> atomList = replacementFrag.getAtomList();
			if (atomList.size() == 2){
				//special case for oxime
				//adds a locant like O1 giving another way of referencing this atom
				Atom numericLocantAtomConnectedToCarbonyl = OpsinTools.depthFirstSearchForAtomWithNumericLocant(carbonylOxygen);
				if (numericLocantAtomConnectedToCarbonyl != null) {
					Atom lastatom = atomList.get(1);
					lastatom.addLocant(lastatom.getElement().toString() + numericLocantAtomConnectedToCarbonyl.getFirstLocant());	
				}
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
			if (replacementFrag.getOutAtomCount() !=1) {
				throw new RuntimeException("OPSIN Bug: Carbonyl replacement fragment expected to have one outatom");
			}
			Atom atomToReplaceCarbonylOxygen = replacementFrag.getOutAtom(0).getAtom();
			replacementFrag.removeOutAtom(0);
			state.fragManager.replaceAtomWithAnotherAtomPreservingConnectivity(carbonylOxygen, atomToReplaceCarbonylOxygen);
			atomToReplaceCarbonylOxygen.setType(carbonylOxygen.getType());//copy the type e.g. if the carbonyl was a suffix this should appear as a suffix
			if (replacementFrag.getTokenEl().getParent() == null) {//incorporate only for the case that replacementFrag came from a functional class element
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
			if (atom.getElement() == ChemEl.O && atom.getCharge()==0){
				List<Atom> neighbours =atom.getAtomNeighbours();
				if (neighbours.size()==1){
					if (neighbours.get(0).getElement() == ChemEl.C){
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

	private void buildAnhydride(List<Element> words) throws StructureBuildingException {
		if (words.size()!=2 && words.size()!=3){
			throw new StructureBuildingException("Unexpected number of words in anhydride. Check wordRules.xml, this is probably a bug");
		}
		Element anhydrideWord = words.get(words.size()-1);
		List<Element> functionalClass =OpsinTools.getDescendantElementsWithTagName(anhydrideWord, FUNCTIONALGROUP_EL);
		if (functionalClass.size()!=1){
			throw new StructureBuildingException("Expected 1 group element found: " + functionalClass.size());
		}
		String anhydrideSmiles = functionalClass.get(0).getAttributeValue(VALUE_ATR);
		int numberOfAnhydrideLinkages =1;
		List<Element> multipliers =OpsinTools.getDescendantElementsWithTagName(anhydrideWord, MULTIPLIER_EL);
		if (multipliers.size() >1){
			throw new StructureBuildingException("Expected 0 or 1 multiplier found: " + multipliers.size());
		}
		if (multipliers.size()==1){
			numberOfAnhydrideLinkages = Integer.parseInt(multipliers.get(0).getAttributeValue(VALUE_ATR));
			multipliers.get(0).detach();
		}
		String anhydrideLocant = null;
		List<Element> anhydrideLocants =OpsinTools.getDescendantElementsWithTagNames(anhydrideWord, new String[]{LOCANT_EL, COLONORSEMICOLONDELIMITEDLOCANT_EL});
		if (anhydrideLocants.size() >1){
			throw new StructureBuildingException("Expected 0 or 1 anhydrideLocants found: " + anhydrideLocants.size());
		}
		if (anhydrideLocants.size()==1){
			anhydrideLocant = anhydrideLocants.get(0).getValue();
			anhydrideLocants.get(0).detach();
		}
		resolveWordOrBracket(state, words.get(0));
		BuildResults br1 = new BuildResults(words.get(0));
		if (br1.getFunctionalAtomCount() ==0){
			throw new StructureBuildingException("Cannot find functionalAtom to form anhydride");
		}
		if (words.size()==3){//asymmetric anhydride
			if (anhydrideLocant!=null){
				throw new StructureBuildingException("Unsupported or invalid anhydride");
			}
			resolveWordOrBracket(state, words.get(1));
			BuildResults br2 = new BuildResults(words.get(1));
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
						OpsinTools.insertAfter(words.get(0), newAcid);
						newAcidBr = new BuildResults(newAcid);
					}
					else{
						newAcidBr =br1;
					}
					formAnhydrideLink(anhydrideSmiles, newAcidBr, br2);
				}
				
			}
			else{
				if (br1.getFunctionalAtomCount()!=1 && br2.getFunctionalAtomCount()!=1 ) {
					throw new StructureBuildingException("Invalid anhydride description");
				}
				formAnhydrideLink(anhydrideSmiles, br1, br2);
			}
		}
		else{//symmetric anhydride
			if (br1.getFunctionalAtomCount()>1){//cyclic anhydride
				if (br1.getFunctionalAtomCount()==2){
					if (numberOfAnhydrideLinkages!=1 || anhydrideLocant !=null ){
						throw new StructureBuildingException("Unsupported or invalid anhydride");
					}
					formAnhydrideLink(anhydrideSmiles, br1, br1);
				}
				else{//cyclic anhydride where group has more than 2 acids
					if (anhydrideLocant ==null){
						throw new StructureBuildingException("Anhydride formation appears to be ambiguous; More than 2 acids, no locants");
					}
					String[] acidLocants =MATCH_COLONORSEMICOLON.split(StringTools.removeDashIfPresent(anhydrideLocant));
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
						String[] locants = acidLocants[i].split(",");
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
						formAnhydrideLink(anhydrideSmiles, oxygen1, oxygen2);
					}
				}
			}
			else{
				if (numberOfAnhydrideLinkages!=1 || anhydrideLocant !=null ){
					throw new StructureBuildingException("Unsupported or invalid anhydride");
				}
				Element newAcid = state.fragManager.cloneElement(state, words.get(0));
				OpsinTools.insertAfter(words.get(0), newAcid);
				BuildResults br2 = new BuildResults(newAcid);
				formAnhydrideLink(anhydrideSmiles, br1, br2);
			}
		}
	}

	/**
	 * Given buildResults for both the acids and the SMILES of the anhydride forms the anhydride bond using the first functionalAtom on each BuildResults
	 * @param anhydrideSmiles
	 * @param acidBr1
	 * @param acidBr2
	 * @throws StructureBuildingException
	 */
	private void formAnhydrideLink(String anhydrideSmiles, BuildResults acidBr1, BuildResults acidBr2)throws StructureBuildingException {
		Atom oxygen1 = acidBr1.getFunctionalAtom(0);
		acidBr1.removeFunctionalAtom(0);
		Atom oxygen2 = acidBr2.getFunctionalAtom(0);
		acidBr2.removeFunctionalAtom(0);
		formAnhydrideLink(anhydrideSmiles, oxygen1, oxygen2);
	}
	
	/**
	 * Given two atoms and the SMILES of the anhydride forms the anhydride bond
	 * @param anhydrideSmiles
	 * @param oxygen1
	 * @param oxygen2
	 * @throws StructureBuildingException
	 */
	private void formAnhydrideLink(String anhydrideSmiles, Atom oxygen1, Atom oxygen2)throws StructureBuildingException {
		if (oxygen1.getElement() != ChemEl.O || oxygen2.getElement() != ChemEl.O || oxygen1.getBondCount()!=1 ||oxygen2.getBondCount()!=1) {
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
	
	private void buildAcidHalideOrPseudoHalide(List<Element> words) throws StructureBuildingException {
		if (!words.get(0).getAttributeValue(TYPE_ATR).equals(WordType.full.toString())){
			throw new StructureBuildingException("Don't alter wordRules.xml without checking the consequences!");
		}
		resolveWordOrBracket(state, words.get(0));
		BuildResults acidBr = new BuildResults(words.get(0));
		int functionalAtomCount =acidBr.getFunctionalAtomCount();
		if (functionalAtomCount==0){
			throw new StructureBuildingException("No functionalAtoms detected!");
		}

		boolean monoMultiplierDetected =false;
		List<Fragment> functionalGroupFragments = new ArrayList<Fragment>();
		for (int i = 1; i < words.size(); i++ ) {
			Element functionalGroupWord =words.get(i);
			List<Element> functionalGroups = OpsinTools.getDescendantElementsWithTagName(functionalGroupWord, FUNCTIONALGROUP_EL);
			if (functionalGroups.size()!=1){
				throw new StructureBuildingException("Expected exactly 1 functionalGroup. Found " + functionalGroups.size());
			}
			
			Fragment monoValentFunctionGroup =state.fragManager.buildSMILES(functionalGroups.get(0).getAttributeValue(VALUE_ATR), FUNCTIONALCLASS_TYPE_VAL, NONE_LABELS_VAL);
			if (functionalGroups.get(0).getAttributeValue(TYPE_ATR).equals(MONOVALENTSTANDALONEGROUP_TYPE_VAL)){
				Atom ideAtom = monoValentFunctionGroup.getDefaultInAtomOrFirstAtom();
				ideAtom.addChargeAndProtons(1, 1);//e.g. make cyanide charge netural
			}
			Element possibleMultiplier = OpsinTools.getPreviousSibling(functionalGroups.get(0));
			functionalGroupFragments.add(monoValentFunctionGroup);
			if (possibleMultiplier != null){
				int multiplierValue = Integer.parseInt(possibleMultiplier.getAttributeValue(VALUE_ATR));
				if (multiplierValue == 1) {
					monoMultiplierDetected = true;
				}
				for (int j = 1; j < multiplierValue; j++) {
					functionalGroupFragments.add(state.fragManager.copyFragment(monoValentFunctionGroup));
				}
				possibleMultiplier.detach();
			}
		}
		int halideCount = functionalGroupFragments.size();
		if (halideCount < functionalAtomCount && halideCount == 1 && !monoMultiplierDetected) {
			//e.g. phosphoric chloride, chloride is implicitly multiplied
			Fragment ideFrag = functionalGroupFragments.get(0);
			for (int i = halideCount; i < functionalAtomCount; i++) {
				functionalGroupFragments.add(state.fragManager.copyFragment(ideFrag));
			}
			halideCount = functionalAtomCount;
		}
		else if (halideCount > functionalAtomCount || (!monoMultiplierDetected && halideCount <functionalAtomCount)){
			throw new StructureBuildingException("Mismatch between number of halide/pseudo halide fragments and acidic oxygens");
		}
		for (int i = halideCount - 1; i>=0; i--) {
			Fragment ideFrag =functionalGroupFragments.get(i);
			Atom ideAtom = ideFrag.getDefaultInAtomOrFirstAtom();
			Atom acidAtom = acidBr.getFunctionalAtom(i);
			if (acidAtom.getElement() != ChemEl.O){
				throw new StructureBuildingException("Atom type expected to be oxygen but was: " +acidAtom.getElement());
			}
			acidBr.removeFunctionalAtom(i);
			Fragment acidFragment =acidAtom.getFrag();
			state.fragManager.replaceAtomWithAnotherAtomPreservingConnectivity(acidAtom, ideAtom);
			state.fragManager.incorporateFragment(ideFrag, acidFragment);
		}
	}
	
	private void buildAdditionCompound(List<Element> words) throws StructureBuildingException {
		Element firstWord = words.get(0);
		if (!firstWord.getAttributeValue(TYPE_ATR).equals(WordType.full.toString())) {
			throw new StructureBuildingException("Don't alter wordRules.xml without checking the consequences!");
		}
		resolveWordOrBracket(state, firstWord);
		Element elementaryAtomEl = StructureBuildingMethods.findRightMostGroupInBracket(firstWord);
		Fragment elementaryAtomFrag = elementaryAtomEl.getFrag();
		Atom elementaryAtom = elementaryAtomFrag.getFirstAtom();
		int charge = elementaryAtom.getCharge();
		List<Fragment> functionalGroupFragments = new ArrayList<Fragment>();
		for (int i = 1; i < words.size(); i++ ) {
			Element functionalGroupWord = words.get(i);
			List<Element> functionalGroups = OpsinTools.getDescendantElementsWithTagName(functionalGroupWord, FUNCTIONALGROUP_EL);
			if (functionalGroups.size() != 1){
				throw new StructureBuildingException("Expected exactly 1 functionalGroup. Found " + functionalGroups.size());
			}
			Element functionGroup = functionalGroups.get(0);
			
			Fragment monoValentFunctionGroup = state.fragManager.buildSMILES(functionGroup.getAttributeValue(VALUE_ATR), FUNCTIONALCLASS_TYPE_VAL, NONE_LABELS_VAL);
			if (functionGroup.getAttributeValue(TYPE_ATR).equals(MONOVALENTSTANDALONEGROUP_TYPE_VAL)){
				Atom ideAtom = monoValentFunctionGroup.getDefaultInAtomOrFirstAtom();
				ideAtom.addChargeAndProtons(1, 1);//e.g. make cyanide and the like charge neutral
			}
			Element possibleMultiplier = OpsinTools.getPreviousSibling(functionGroup);
			functionalGroupFragments.add(monoValentFunctionGroup);
			if (possibleMultiplier != null) {
				int multiplierValue = Integer.parseInt(possibleMultiplier.getAttributeValue(VALUE_ATR));
				for (int j = 1; j < multiplierValue; j++) {
					functionalGroupFragments.add(state.fragManager.copyFragment(monoValentFunctionGroup));
				}
				possibleMultiplier.detach();
			}
			else if (words.size() == 2) {//silicon chloride -->silicon tetrachloride
				int incomingBondOrder = elementaryAtom.getIncomingValency();
				int expectedValency;
				if (charge > 0) {
					expectedValency = incomingBondOrder + charge;
				}
				else{
					if (elementaryAtom.getProperty(Atom.OXIDATION_NUMBER) != null) {
						expectedValency = elementaryAtom.getProperty(Atom.OXIDATION_NUMBER);
					}
					else{
						if (elementaryAtomEl.getAttribute(COMMONOXIDATIONSTATESANDMAX_ATR) != null) {
							String[] typicalOxidationStates = elementaryAtomEl.getAttributeValue(COMMONOXIDATIONSTATESANDMAX_ATR).split(":")[0].split(",");
							expectedValency = Integer.parseInt(typicalOxidationStates[0]);
						}
						else{
							expectedValency = ValencyChecker.getPossibleValencies(elementaryAtom.getElement(), charge)[0];
						}
					}
				}
				int implicitMultiplier = expectedValency - incomingBondOrder > 1 ? expectedValency - incomingBondOrder : 1;
				for (int j = 1; j < implicitMultiplier; j++) {
					functionalGroupFragments.add(state.fragManager.copyFragment(monoValentFunctionGroup));
				}
			}
		}
		if (charge > 0) {
			elementaryAtom.setCharge(charge - functionalGroupFragments.size());
		}
		
		//[AlH3] --> [AlH4-] , [AlH4] --> [AlH4-]
		applyAluminiumHydrideSpecialCase(firstWord, elementaryAtom, functionalGroupFragments);

		int halideCount = functionalGroupFragments.size();
		Integer maximumVal = ValencyChecker.getMaximumValency(elementaryAtom.getElement(), elementaryAtom.getCharge());
		if (maximumVal != null && halideCount > maximumVal) {
			throw new StructureBuildingException("Too many halides/psuedo halides addded to " +elementaryAtom.getElement());
		}
		for (int i = halideCount - 1; i >= 0; i--) {
			Fragment ideFrag = functionalGroupFragments.get(i);
			Atom ideAtom = ideFrag.getDefaultInAtomOrFirstAtom();
			state.fragManager.incorporateFragment(ideFrag, ideAtom, elementaryAtomFrag, elementaryAtom, 1);
		}
	}

	private void applyAluminiumHydrideSpecialCase(Element firstWord, Atom elementaryAtom,
			List<Fragment> functionalGroupFragments) throws StructureBuildingException {
		if ((elementaryAtom.getElement() == ChemEl.Al || elementaryAtom.getElement() == ChemEl.B)
				&& elementaryAtom.getCharge() == 0) {
			if (functionalGroupFragments.size() == 3) {
				if (functionalGroupFragments.get(0).getDefaultInAtomOrFirstAtom().getElement() == ChemEl.H
					&& functionalGroupFragments.get(1).getDefaultInAtomOrFirstAtom().getElement() == ChemEl.H
					&& functionalGroupFragments.get(2).getDefaultInAtomOrFirstAtom().getElement() == ChemEl.H) {
					Element counterCationWordRule = OpsinTools.getPreviousSibling(firstWord.getParent());
					if (counterCationWordRule != null && counterCationWordRule.getChildCount() == 1) {
						Element word =counterCationWordRule.getFirstChildElement(WORD_EL);
						if (word != null && word.getChildCount() ==1) {
							Element root = word.getFirstChildElement(ROOT_EL);
							if (root != null && root.getChildCount() ==1) {
								Element group = root.getFirstChildElement(GROUP_EL);
								if (group != null && ELEMENTARYATOM_TYPE_VAL.equals(group.getAttributeValue(TYPE_ATR))) {
									ChemEl chemEl = group.getFrag().getFirstAtom().getElement();
									if (chemEl == ChemEl.Li || chemEl == ChemEl.Na || chemEl == ChemEl.K || chemEl == ChemEl.Rb || chemEl == ChemEl.Cs) {
										functionalGroupFragments.add(state.fragManager.copyFragment(functionalGroupFragments.get(0)));
										elementaryAtom.setCharge(-1);
									}
								}
							}
						}
					}
				
				}
			}
			else if (functionalGroupFragments.size() == 4) {
				if (functionalGroupFragments.get(0).getDefaultInAtomOrFirstAtom().getElement() == ChemEl.H
					&& functionalGroupFragments.get(1).getDefaultInAtomOrFirstAtom().getElement() == ChemEl.H
					&& functionalGroupFragments.get(2).getDefaultInAtomOrFirstAtom().getElement() == ChemEl.H
					&& functionalGroupFragments.get(3).getDefaultInAtomOrFirstAtom().getElement() == ChemEl.H) {
					elementaryAtom.setCharge(-1);
				}
			}
		}
	}
	

	private void buildGlycol(List<Element> words) throws StructureBuildingException {
		int wordIndice  = 0;
		resolveWordOrBracket(state, words.get(wordIndice));//the group
		Element finalGroup = findRightMostGroupInWordOrWordRule(words.get(wordIndice));
		Fragment theDiRadical = finalGroup.getFrag();
		if (theDiRadical.getOutAtomCount()!=2){
			throw new StructureBuildingException("Glycol class names (e.g. ethylene glycol) expect two outAtoms. Found: " + theDiRadical.getOutAtomCount() );
		}
		wordIndice++;
		if (wordIndice >= words.size() || !words.get(wordIndice).getAttributeValue(TYPE_ATR).equals(WordType.functionalTerm.toString())){
			throw new StructureBuildingException("Glycol functionalTerm word expected");
		}
		List<Element> functionalClassEls = OpsinTools.getDescendantElementsWithTagName(words.get(wordIndice), FUNCTIONALCLASS_EL);
		if (functionalClassEls.size()!=1){
			throw new StructureBuildingException("Glycol functional class not found where expected");
		}
		
		OutAtom outAtom1 = theDiRadical.getOutAtom(0);
		Atom chosenAtom1 = outAtom1.isSetExplicitly() ? outAtom1.getAtom() : findAtomForUnlocantedRadical(state, theDiRadical, outAtom1);
		Fragment functionalFrag =state.fragManager.buildSMILES(functionalClassEls.get(0).getAttributeValue(VALUE_ATR), FUNCTIONALCLASS_TYPE_VAL, NONE_LABELS_VAL);
		if (outAtom1.getValency() != 1){
			throw new StructureBuildingException("OutAtom has unexpected valency. Expected 1. Actual: " + outAtom1.getValency());
		}
		state.fragManager.createBond(chosenAtom1, functionalFrag.getFirstAtom(), 1);
		state.fragManager.incorporateFragment(functionalFrag, theDiRadical);
		
		OutAtom outAtom2 = theDiRadical.getOutAtom(1);
		Atom chosenAtom2 = outAtom2.isSetExplicitly() ? outAtom2.getAtom() : findAtomForUnlocantedRadical(state, theDiRadical, outAtom2);
		Fragment hydroxy =state.fragManager.buildSMILES("O", FUNCTIONALCLASS_TYPE_VAL, NONE_LABELS_VAL);
		if (outAtom2.getValency() != 1){
			throw new StructureBuildingException("OutAtom has unexpected valency. Expected 1. Actual: " + outAtom2.getValency());
		}
		state.fragManager.createBond(chosenAtom2, hydroxy.getFirstAtom(), 1);
		state.fragManager.incorporateFragment(hydroxy, theDiRadical);
		theDiRadical.removeOutAtom(1);
		theDiRadical.removeOutAtom(0);
	}
	

	/**
	 * Handles Glcyol ethers nomenclature e.g.
	 * triethylene glycol n-butyl ether
	 * tripropylene glycol methyl ether
	 * dipropylene glycol methyl ether acetate
	 * @param words
	 * @throws StructureBuildingException
	 */
	private void buildGlycolEther(List<Element> words) throws StructureBuildingException {
		List<Element> wordsToAttachToGlycol = new ArrayList<Element>();
		Element glycol =words.get(0);
		resolveWordOrBracket(state, glycol);//if this actually is something like ethylene glycol this is a no-op as it will already have been resolved
		if (!glycol.getAttributeValue(TYPE_ATR).equals(WordType.full.toString())){
			throw new StructureBuildingException("OPSIN Bug: Cannot find glycol word!");
		}
		for (int i = 1; i < words.size(); i++) {
			Element wordOrWordRule =words.get(i);
			//ether ignored
			if (!wordOrWordRule.getAttributeValue(TYPE_ATR).equals(WordType.functionalTerm.toString())){
				resolveWordOrBracket(state, wordOrWordRule);//the substituent to attach
				wordsToAttachToGlycol.add(wordOrWordRule);
			}
			else if (!wordOrWordRule.getAttributeValue(VALUE_ATR).equalsIgnoreCase("ether")){
				throw new StructureBuildingException("Unexpected word encountered when applying glycol ether word rule " + wordOrWordRule.getAttributeValue(VALUE_ATR));
			}
		}
		int numOfEthers = wordsToAttachToGlycol.size();
		if (numOfEthers == 0) {
			throw new StructureBuildingException("OPSIN Bug: Unexpected number of substituents for glycol ether");
		}
		Element finalGroup = findRightMostGroupInWordOrWordRule(glycol);
		List<Atom> hydroxyAtoms = FragmentTools.findHydroxyGroups(finalGroup.getFrag());
		if (hydroxyAtoms.size() == 0) {
			throw new StructureBuildingException("No hydroxy groups found in: " + finalGroup.getValue() + " to form ether");
		}
		if (hydroxyAtoms.size() < numOfEthers) {
			throw new StructureBuildingException("Insufficient hydroxy groups found in: " + finalGroup.getValue() + " to form required number of ethers");
		}
		for (int i = 0; i < numOfEthers; i++) {
			BuildResults br = new BuildResults(wordsToAttachToGlycol.get(i));
			if (br.getOutAtomCount() >0){//form ether
				state.fragManager.createBond(hydroxyAtoms.get(i), br.getOutAtom(0).getAtom(), 1);
				br.removeOutAtom(0);
			}
			else if (br.getFunctionalAtomCount() >0){//form ester
				Atom ateAtom = br.getFunctionalAtom(0);
				ateAtom.neutraliseCharge();
				state.fragManager.replaceAtomWithAnotherAtomPreservingConnectivity(hydroxyAtoms.get(i), br.getFunctionalAtom(0));
				br.removeFunctionalAtom(0);
			}
			else{
				throw new StructureBuildingException("Word had neither an outAtom or a functionalAtom! hence neither and ether or ester could be formed : " + wordsToAttachToGlycol.get(i).getAttributeValue(VALUE_ATR));
			}
		}
	}

	/**
	 * Builds acetals/ketals/hemiacetals/hemiketals and chalcogen analogues
	 * The distinction between acetals and ketals is not enforced (ketals are a subset of acetals)
	 * @param words
	 * @throws StructureBuildingException
	 */
	private void buildAcetal(List<Element> words) throws StructureBuildingException {
		for (int i = 0; i < words.size()-1; i++) {
			resolveWordOrBracket(state, words.get(i));
		}
		BuildResults substituentsBr = new BuildResults();
		for (int i = 1; i < words.size()-1; i++) {
			Element currentWord = words.get(i);
			BuildResults substituentBr = new BuildResults(currentWord);
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
		Fragment rootFragment = rightMostGroup.getFrag();//the group which will be modified
		List<Atom> carbonylOxygen= findCarbonylOxygens(rootFragment, new ArrayList<String>());
		Element functionalWord = words.get(words.size()-1);
		List<Element> functionalClasses = OpsinTools.getDescendantElementsWithTagName(functionalWord, FUNCTIONALCLASS_EL);
		if (functionalClasses.size()!=1){
			throw new StructureBuildingException("OPSIN bug: unable to find acetal functionalClass");
		}
		Element functionalClassEl = functionalClasses.get(0);
		String functionalClass = functionalClassEl.getValue();
		Element beforeAcetal = OpsinTools.getPreviousSibling(functionalClassEl);
		int numberOfAcetals =1;
		String[] elements = functionalClassEl.getAttributeValue(VALUE_ATR).split(",");
		if (beforeAcetal != null){
			if (beforeAcetal.getName().equals(MULTIPLIER_EL)){
				numberOfAcetals = Integer.parseInt(beforeAcetal.getAttributeValue(VALUE_ATR));
			}
			else{
				replaceChalcogensInAcetal(functionalClassEl, elements);
			}
		}

		if (carbonylOxygen.size() < numberOfAcetals){
			throw new StructureBuildingException("Insufficient carbonyls to form " + numberOfAcetals +" " + functionalClass );
		}
		boolean hemiacetal = functionalClass.contains("hemi");
		List<Fragment> acetalFrags = new ArrayList<Fragment>();
		for (int i = 0; i < numberOfAcetals; i++) {
			acetalFrags.add(formAcetal(carbonylOxygen, elements));
		}
		int bondsToForm = hemiacetal ? numberOfAcetals : 2*numberOfAcetals;
		if (substituentsBr.getOutAtomCount()!=bondsToForm){
			throw new StructureBuildingException("incorrect number of susbtituents when forming " + functionalClass);
		}
		connectSubstituentsToAcetal(acetalFrags, substituentsBr, hemiacetal);
	}

	private void replaceChalcogensInAcetal(Element functionalClassEl, String[] elements) throws StructureBuildingException {
		Element currentEl = functionalClassEl.getParent().getChild(0);
		int multiplier = 1;
		if (currentEl.getName().equals(MULTIPLIER_EL)){
			multiplier = Integer.parseInt(currentEl.getAttributeValue(VALUE_ATR));
			if (multiplier > 2){
				throw new StructureBuildingException(functionalClassEl.getValue() + " only has two oxygen!");
			}
			currentEl = OpsinTools.getNextSibling(currentEl);
		}
		int i = 0;
		while(currentEl != functionalClassEl) {
			if (currentEl.getName().equals(GROUP_EL)) {
				for (int j = 0; j < multiplier; j++) {
					if (i == 2) {
						throw new StructureBuildingException(functionalClassEl.getValue() + " only has two oxygen!");
					}
					if (!elements[i].equals("O")){
						throw new StructureBuildingException("Replacement on " + functionalClassEl.getValue() + " can only be used to replace oxygen!");
					}
					elements[i++] = currentEl.getAttributeValue(VALUE_ATR);
				}
			}
			else {
				throw new StructureBuildingException("Unexpected element before acetal");
			}
			currentEl = OpsinTools.getNextSibling(currentEl);
		}
	}

	private Fragment formAcetal(List<Atom> carbonylOxygen, String[] elements) throws StructureBuildingException {
		Atom neighbouringCarbon = carbonylOxygen.get(0).getAtomNeighbours().get(0);
		state.fragManager.removeAtomAndAssociatedBonds(carbonylOxygen.get(0));
		carbonylOxygen.remove(0);
		Fragment acetalFrag = state.fragManager.buildSMILES(StringTools.arrayToString(elements, "."),"",NONE_LABELS_VAL);
		FragmentTools.assignElementLocants(acetalFrag, new ArrayList<Fragment>());
		List<Atom> acetalAtomList = acetalFrag.getAtomList();
		Atom atom1 = acetalAtomList.get(0);
		state.fragManager.createBond(neighbouringCarbon, atom1, 1);
		Atom atom2 = acetalAtomList.get(1);
		state.fragManager.createBond(neighbouringCarbon, atom2, 1);
		state.fragManager.incorporateFragment(acetalFrag, neighbouringCarbon.getFrag());
		return acetalFrag;
	}
	
	private boolean buildAlcoholEster(List<Element> words, int numberOfWordRules) throws StructureBuildingException {
		for (Element word : words) {
			if (!WordType.full.toString().equals(word.getAttributeValue(TYPE_ATR))){
				throw new StructureBuildingException("Bug in word rule for potentialAlcoholEster");
			}
			resolveWordOrBracket(state, word);
		}
		int ateWords = words.size() -1;
		if (ateWords < 1){
			throw new StructureBuildingException("Bug in word rule for potentialAlcoholEster");
		}
		
		Fragment potentialAlcoholFragment = findRightMostGroupInWordOrWordRule(words.get(0)).getFrag();
		List<Atom> hydroxyAtoms = FragmentTools.findHydroxyGroups(potentialAlcoholFragment);
		
		List<Atom> chosenHydroxyAtoms = new ArrayList<Atom>();
		List<BuildResults> ateBuildResults = new ArrayList<BuildResults>();
		for (int i = 1; i < words.size(); i++) {
			Element ateWord = words.get(i);
			BuildResults wordBr = new BuildResults(ateWord);
			if (isAppropriateAteGroupForAlcoholEster(ateWord, wordBr)) {
				String locant = ateWord.getAttributeValue(LOCANT_ATR);
				if (locant != null) {
					Atom atomOnAlcoholFragment = potentialAlcoholFragment.getAtomByLocantOrThrow(locant);
					if (!hydroxyAtoms.contains(atomOnAlcoholFragment) || chosenHydroxyAtoms.contains(atomOnAlcoholFragment)) {
						atomOnAlcoholFragment = potentialAlcoholFragment.getAtomByLocantOrThrow("O" + locant);
					}
					if (!hydroxyAtoms.contains(atomOnAlcoholFragment) || chosenHydroxyAtoms.contains(atomOnAlcoholFragment)) {
						throw new StructureBuildingException(locant + " did not point to a hydroxy group to be used for ester formation");
					}
					chosenHydroxyAtoms.add(atomOnAlcoholFragment);
				}
				else if (words.size() == 2) {
					//special case for adenosine triphosphate and the like
					//also true for pyridoxine derivatives
					//guess that locant might be 5'
					Atom atomOnAlcoholFragment = potentialAlcoholFragment.getAtomByLocant("O5'");
					if (hydroxyAtoms.contains(atomOnAlcoholFragment)) {
						chosenHydroxyAtoms.add(atomOnAlcoholFragment);
					}
				}
				ateBuildResults.add(wordBr);
			}
			else {
				return false;
			}
		}

		if (chosenHydroxyAtoms.size() < ateWords) {
			if (!chosenHydroxyAtoms.isEmpty()) {
				throw new RuntimeException("OPSIN Bug: Either all or none of the esters should be locanted in alcohol ester rule");
			}
			if (hydroxyAtoms.size() == ateWords  || hydroxyAtoms.size() > ateWords && (AmbiguityChecker.allAtomsEquivalent(hydroxyAtoms) || potentialAlcoholFragment.getTokenEl().getValue().equals("glycerol") )) {
				for (int i = 0; i < ateWords; i++) {
					chosenHydroxyAtoms.add(hydroxyAtoms.get(i));
				}
			}
			else {
				return false;
			}
		}

		for (int i = 0; i < ateWords; i++) {
			BuildResults br = ateBuildResults.get(i);
			Element ateWord = words.get(i + 1);
			Element ateGroup = findRightMostGroupInWordOrWordRule(ateWord);
			if (ateGroup.getAttribute(NUMBEROFFUNCTIONALATOMSTOREMOVE_ATR) == null && numberOfWordRules == 1) {
				//by convention [O-] are implicitly converted to [OH] when phosphates/sulfates are attached
				//If word rules is > 1 this will be done or not done as part of charge balancing
				for (int j = br.getFunctionalAtomCount() -1; j >= 1; j--) {
					Atom atomToDefunctionalise = br.getFunctionalAtom(j);
					br.removeFunctionalAtom(j);
					atomToDefunctionalise.neutraliseCharge();
				}
			}
			Atom functionalAtom = br.getFunctionalAtom(0);
			br.removeFunctionalAtom(0);
			functionalAtom.neutraliseCharge();
			state.fragManager.replaceAtomWithAnotherAtomPreservingConnectivity(functionalAtom, chosenHydroxyAtoms.get(i));
		}
		return true;
	}
	
	private void buildAmineDiConjunctiveSuffix(List<Element> words) throws StructureBuildingException {
		for (Element word : words) {
			if (!WordType.full.toString().equals(word.getAttributeValue(TYPE_ATR))){
				throw new StructureBuildingException("Bug in word rule for amineDiConjunctiveSuffix");
			}
			resolveWordOrBracket(state, word);
		}
		if (words.size() != 3) {
			throw new StructureBuildingException("Unexpected number of words encountered when processing name of type amineDiConjunctiveSuffix, expected 3 but found: " + words.size());
		}
		Element aminoAcid = findRightMostGroupInWordOrWordRule(words.get(0));
		if (aminoAcid == null) {
			throw new RuntimeException("OPSIN Bug: failed to find amino acid");
		}
		Atom amineAtom = aminoAcid.getFrag().getDefaultInAtom();
		if (amineAtom == null) {
			throw new StructureBuildingException("OPSIN did not know where the amino acid amine was located");
		}
		
		for (int i = 1; i < words.size(); i++) {
			Element word = words.get(i);
			Fragment suffixLikeGroup = findRightMostGroupInWordOrWordRule(word).getFrag();
			String locant = word.getAttributeValue(LOCANT_ATR);
			if (locant != null){
				if (!locant.equals("N")) {
					throw new RuntimeException("OPSIN Bug: locant expected to be N but was: " + locant);
				}
			}
			Atom atomToConnectToOnConjunctiveFrag = FragmentTools.lastNonSuffixCarbonWithSufficientValency(suffixLikeGroup);
			if (atomToConnectToOnConjunctiveFrag == null) {
				throw new StructureBuildingException("OPSIN Bug: Unable to find non suffix carbon with sufficient valency");
			}
			state.fragManager.createBond(atomToConnectToOnConjunctiveFrag, amineAtom, 1);
		}
	}

	private static final Pattern matchCommonCarboxylicSalt = Pattern.compile("tri-?fluoro-?acetate?$", Pattern.CASE_INSENSITIVE);
	private static final Pattern matchCommonEsterFormingInorganicSalt = Pattern.compile("(ortho-?)?(bor|phosphor|phosphate?|phosphite?)|carbam|carbon|sulfur|sulfate?|sulfite?|diphosphate?|triphosphate?", Pattern.CASE_INSENSITIVE);

	/**
	 * CAS endorses the use of ...ol ...ate names means esters
	 * but only for cases involving "common acids":
	 * Acetic acid; Benzenesulfonic acid; Benzenesulfonic acid, 4-methyl-; Benzoic acid and its monoamino, mononitro, and dinitro derivatives;
	 * Boric acid (H3BO3); Carbamic acid; Carbamic acid, N-methyl-; Carbamic acid, N-phenyl-; Carbonic acid; Formic acid; Methanesulfonic acid;
	 * Nitric acid; Phosphoric acid; Phosphorodithioic acid; Phosphorothioic acid; Phosphorous acid; Propanoic acid; Sulfuric acid; and Sulfurous acid.
	 * ...unless the alcohol component is also common.
	 * 
	 * As in practice a lot of use won't be from CAS names we use the following heuristic:
	 * Is locanted OR
	 * Has 1 functional atom  (And not common salt e.g. Trifluoroacetate) OR
	 * common phosphorus/sulfur ate including di/tri phosphate
	 * @param ateWord
	 * @param wordBr
	 * @return
	 * @throws StructureBuildingException 
	 */
	private boolean isAppropriateAteGroupForAlcoholEster(Element ateWord, BuildResults wordBr) throws StructureBuildingException {
		if (wordBr.getFunctionalAtomCount() > 0) {
			if (ateWord.getAttributeValue(LOCANT_ATR) != null) {
				//locanted, so locant must be used for this purpose
				return true;
			}
			if (wordBr.getFunctionalAtomCount() == 1) {
				if (matchCommonCarboxylicSalt.matcher(ateWord.getAttributeValue(VALUE_ATR)).find()) {
					return false;
				}
				return true;
			}
			String ateGroupText = findRightMostGroupInWordOrWordRule(ateWord).getValue();
			//e.g. triphosphate
			if (matchCommonEsterFormingInorganicSalt.matcher(ateGroupText).matches()) {
				return true;
			}
			
		}
		return false;
	}

	private void splitAlcoholEsterRuleIntoTwoSimpleWordRules(List<Element> words) {
		Element firstGroup = words.get(0);
		Element wordRule = firstGroup.getParent();
		wordRule.getAttribute(WORDRULE_ATR).setValue(WordRule.simple.toString());
		wordRule.getAttribute(VALUE_ATR).setValue(firstGroup.getAttributeValue(VALUE_ATR));

		Element newWordRule = new GroupingEl(WORDRULE_EL);
		newWordRule.addAttribute(TYPE_ATR, WordType.full.toString());
		newWordRule.addAttribute(WORDRULE_ATR, WordRule.simple.toString());
		newWordRule.addAttribute(VALUE_ATR, words.get(1).getAttributeValue(VALUE_ATR));
		OpsinTools.insertAfter(wordRule, newWordRule);
		for (int i = 1; i < words.size(); i++) {
			Element word = words.get(i);
			word.detach();
			newWordRule.addChild(word);
		}
	}

	private void connectSubstituentsToAcetal(List<Fragment> acetalFrags, BuildResults subBr, boolean hemiacetal) throws StructureBuildingException {
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
							if (atomList.get(0).getBondCount()==1){
								atomToUse = atomList.get(0);
								break;
							}
							else if (atomList.get(1).getBondCount()==1){
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
				if (atomList.get(0).getBondCount()==1){
					atomToUse = atomList.get(0);
				}
				else if (atomList.get(1).getBondCount()==1){
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

	private void buildCyclicPeptide(List<Element> words) throws StructureBuildingException {
		if (words.size() != 2){
			throw new StructureBuildingException("OPSIN Bug: Expected 2 words in cyclic peptide name, found: " + words.size());
		}
		Element peptide = words.get(1);
		resolveWordOrBracket(state, peptide);
		BuildResults peptideBr = new BuildResults(peptide);
		if (peptideBr.getOutAtomCount() ==1){
			Atom outAtom = getOutAtomTakingIntoAccountWhetherSetExplicitly(peptideBr, 0);
			List<Element> aminoAcids = OpsinTools.getDescendantElementsWithTagNameAndAttribute(peptide, GROUP_EL, TYPE_ATR, AMINOACID_TYPE_VAL);
			if (aminoAcids.size() < 2){
				throw new StructureBuildingException("Cyclic peptide building failed: Requires at least two amino acids!");
			}
			Atom inAtom = aminoAcids.get(0).getFrag().getDefaultInAtomOrFirstAtom();

			state.fragManager.createBond(outAtom, inAtom, peptideBr.getOutAtom(0).getValency());
			peptideBr.removeAllOutAtoms();
		}
		else{
			throw new StructureBuildingException("Cyclic peptide building failed: Expected 1 outAtoms, found: " +peptideBr.getOutAtomCount());
		}
	}

	private void buildPolymer(List<Element> words) throws StructureBuildingException {
		if (words.size()!=2){
			throw new StructureBuildingException("Currently unsupported polymer name type");
		}
		Element polymer = words.get(1);
		resolveWordOrBracket(state, polymer);
		BuildResults polymerBr = new BuildResults(polymer);
		if (polymerBr.getOutAtomCount() ==2){
			Atom inAtom = getOutAtomTakingIntoAccountWhetherSetExplicitly(polymerBr, 0);
			Atom outAtom = getOutAtomTakingIntoAccountWhetherSetExplicitly(polymerBr, 1);
			/*
			 * We assume the polymer repeats so as an approximation we create an R group with the same element as the group at the other end of polymer (with valency equal to the bondorder of the Rgroup so no H added)
			 */
			Atom rGroup1 =state.fragManager.buildSMILES("[" + outAtom.getElement().toString() + "|" + polymerBr.getOutAtom(0).getValency() + "]", "", "alpha").getFirstAtom();
			rGroup1.setProperty(Atom.ATOM_CLASS, 1);
			state.fragManager.createBond(inAtom, rGroup1, polymerBr.getOutAtom(0).getValency());

			Atom rGroup2 =state.fragManager.buildSMILES("[" + inAtom.getElement().toString() + "|" + polymerBr.getOutAtom(1).getValency() + "]", "", "omega").getFirstAtom();
			rGroup2.setProperty(Atom.ATOM_CLASS, 2);
			state.fragManager.createBond(outAtom, rGroup2, polymerBr.getOutAtom(1).getValency());
			polymerAttachmentPoints.add(rGroup1);
			polymerAttachmentPoints.add(rGroup2);
			polymerBr.removeAllOutAtoms();
		}
		else{
			throw new StructureBuildingException("Polymer building failed: Two termini were not found; Expected 2 outAtoms, found: " +polymerBr.getOutAtomCount());
		}
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
				Set<Atom> degenerateAtoms = possibleAtom.getProperty(Atom.AMBIGUOUS_ELEMENT_ASSIGNMENT);
				if (degenerateAtoms != null){
					degenerateAtoms.remove(possibleAtom);
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
					Set<Atom> degenerateAtoms = possibleAtom.getProperty(Atom.AMBIGUOUS_ELEMENT_ASSIGNMENT);
					if (degenerateAtoms != null){
						degenerateAtoms.remove(possibleAtom);
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
				if (isElementSymbol && possibleAtom.getElement().toString().equals(locant)){
					mainGroupBR.removeFunctionalAtom(i);
					return possibleAtom;
				}
				Set<Atom> degenerateAtoms = possibleAtom.getProperty(Atom.AMBIGUOUS_ELEMENT_ASSIGNMENT);
				if (degenerateAtoms != null){
					boolean foundAtom = false;
					for (Atom a : degenerateAtoms) {
						if (a.hasLocant(locant) || (isElementSymbol && a.getElement().toString().equals(locant))){
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
							ChemEl originalChemEl = possibleAtom.getElement();
							possibleAtom.setElement(a.getElement());
							a.setElement(originalChemEl);
							mainGroupBR.removeFunctionalAtom(i);
							foundAtom =true;
							break;
						}
					}
					if (foundAtom){
						degenerateAtoms.remove(possibleAtom);
						return possibleAtom;
					}
				}
			}
		}

		throw new StructureBuildingException("Cannot find functional atom with locant: " +locant + " to form an ester with");
	}

	/**
	 * Applies explicit stoichiometry, charge balancing and fractional multipliers
	 * @param molecule
	 * @param wordRules
	 * @throws StructureBuildingException
	 */
	private void manipulateStoichiometry(Element molecule, List<Element> wordRules) throws StructureBuildingException {
		boolean explicitStoichiometryPresent = applyExplicitStoichiometryIfProvided(wordRules);
		boolean chargedFractionalGroup = false;
		List<Element> wordRulesWithFractionalMultipliers = new ArrayList<Element>(0);
		for (Element wordRule : wordRules) {
			Element fractionalMultiplier = wordRule.getChild(0);
			while (fractionalMultiplier.getChildCount() != 0){
				fractionalMultiplier = fractionalMultiplier.getChild(0);
			}
			if (fractionalMultiplier.getName().equals(FRACTIONALMULTIPLIER_EL)) {
				if (explicitStoichiometryPresent) {
					throw new StructureBuildingException("Fractional multipliers should not be used in conjunction with explicit stoichiometry");
				}
				String[] value = fractionalMultiplier.getAttributeValue(VALUE_ATR).split("/");
				if (value.length != 2) {
					throw new RuntimeException("OPSIN Bug: malformed fractional multiplier: " + fractionalMultiplier.getAttributeValue(VALUE_ATR));
				}
				try {
					int numerator = Integer.parseInt(value[0]);
					int denominator = Integer.parseInt(value[1]);
					if (denominator != 2) {
						throw new RuntimeException("Only fractions of a 1/2 currently supported");
					}
					for (int j = 1; j < numerator; j++) {
						Element clone = state.fragManager.cloneElement(state, wordRule);
						OpsinTools.insertAfter(wordRule, clone);
						wordRulesWithFractionalMultipliers.add(clone);
					}
				}
				catch (NumberFormatException e) {
					throw new RuntimeException("OPSIN Bug: malformed fractional multiplier: " + fractionalMultiplier.getAttributeValue(VALUE_ATR));
				}
				//don't detach the fractional multiplier to avoid charge balancing multiplication (cf. handling of mono)
				wordRulesWithFractionalMultipliers.add(wordRule);
				if (new BuildResults(wordRule).getCharge() !=0){
					chargedFractionalGroup = true;
				}
			}
		}
		if (wordRulesWithFractionalMultipliers.size() > 0) {
			if (wordRules.size() == 1) {
				throw new StructureBuildingException("Unexpected fractional multiplier found at start of word");
			}
			if (chargedFractionalGroup) {
				for (Element wordRule : wordRules) {
					if (wordRulesWithFractionalMultipliers.contains(wordRule)) {
						continue;
					}
					Element clone = state.fragManager.cloneElement(state, wordRule);
					OpsinTools.insertAfter(wordRule, clone);
				}
			}
		}
		boolean saltExpected = molecule.getAttribute(ISSALT_ATR) != null;
		if (saltExpected) {
			deprotonateAcidIfSaltWithMetal(molecule);
		}
		int overallCharge = state.fragManager.getOverallCharge();
		if (overallCharge!=0 && wordRules.size() >1){//a net charge is present! Could just mean the counterion has not been specified though
			balanceChargeIfPossible(molecule, overallCharge, explicitStoichiometryPresent);
		}
		if (wordRulesWithFractionalMultipliers.size() > 0 && !chargedFractionalGroup) {
			for (Element wordRule : molecule.getChildElements(WORDRULE_EL)) {
				if (wordRulesWithFractionalMultipliers.contains(wordRule)) {
					continue;
				}
				Element clone = state.fragManager.cloneElement(state, wordRule);
				OpsinTools.insertAfter(wordRule, clone);
			}
		}

	}

	private boolean applyExplicitStoichiometryIfProvided(List<Element> wordRules) throws StructureBuildingException {
		boolean explicitStoichiometryPresent =false;
		for (Element wordRule : wordRules) {
			if (wordRule.getAttribute(STOICHIOMETRY_ATR)!=null){
				int stoichiometry = Integer.parseInt(wordRule.getAttributeValue(STOICHIOMETRY_ATR));
				wordRule.removeAttribute(wordRule.getAttribute(STOICHIOMETRY_ATR));
				for (int j = 1; j < stoichiometry; j++) {
					Element clone = state.fragManager.cloneElement(state, wordRule);
					OpsinTools.insertAfter(wordRule, clone);
				}
				explicitStoichiometryPresent =true;
			}
		}
		return explicitStoichiometryPresent;
	}
	

	private void deprotonateAcidIfSaltWithMetal(Element molecule) {
		List<BuildResults> positivelyChargedComponents = new ArrayList<BuildResults>();
		List<BuildResults> negativelyChargedComponents = new ArrayList<BuildResults>();
		List<BuildResults> neutralComponents = new ArrayList<BuildResults>();
		List<Element> wordRules = molecule.getChildElements(WORDRULE_ATR);
		for (Element wordRule : wordRules) {
			BuildResults br = new BuildResults(wordRule);
			int charge = br.getCharge();
			if (charge > 0) {
				positivelyChargedComponents.add(br);
			}
			else if (charge < 0) {
				negativelyChargedComponents.add(br);
			}
			else {
				neutralComponents.add(br);
			}
		}
		if (negativelyChargedComponents.size() == 0 && (positivelyChargedComponents.size() > 0 || getMetalsThatCanBeImplicitlyCations(molecule).size() > 0)) {
			for (int i = neutralComponents.size() - 1; i>=0; i--) {
				List<Atom> functionalAtoms = new ArrayList<Atom>();
				for (Fragment f : neutralComponents.get(i).getFragments()) {
					for (int j = 0; j < f.getFunctionalAtomCount(); j++) {
						functionalAtoms.add(f.getFunctionalAtom(j).getAtom());
					}
				}
				for (Atom functionalAtom : functionalAtoms) {
					if (functionalAtom.getCharge() == 0 && functionalAtom.getIncomingValency() == 1){
						functionalAtom.addChargeAndProtons(-1, -1);
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
	 * @param molecule
	 * @param explicitStoichiometryPresent 
	 * @param overallCharge
	 * @throws StructureBuildingException
	 */
	private void balanceChargeIfPossible(Element molecule, int overallCharge, boolean explicitStoichiometryPresent) throws StructureBuildingException {
		List<Element> wordRules = molecule.getChildElements(WORDRULE_ATR);

		List<Element> positivelyChargedComponents = new ArrayList<Element>();
		List<Element> negativelyChargedComponents = new ArrayList<Element>();
		Map<Element, Integer> componentToChargeMapping = new HashMap<Element, Integer>();
		Map<Element, BuildResults> componentToBR = new HashMap<Element, BuildResults>();
		
		List<Element> cationicElements = getMetalsThatCanBeImplicitlyCations(molecule);
		overallCharge = setCationicElementsToTypicalCharge(cationicElements, overallCharge);
		if (overallCharge == 0) {
			return;
		}
		if (overallCharge == -2 && triHalideSpecialCase(wordRules)) {
			//e.g. three iodides --> triiodide ion
			return;
		}
		
		for (Element wordRule : wordRules) {
			BuildResults br = new BuildResults(wordRule);
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
		if (cationicElements.size() ==1 && overallCharge < 0) {//e.g. manganese tetrachloride [Mn+2]-->[Mn+4]
			boolean mustBeCommonOxidationState = negativelyChargedComponents.size() == 1 && negativelyChargedComponents.get(0).getChildElements(WORD_EL).size() == 1;
			//For simple case e.g. silver oxide, constrain the metal to common oxidation states, otherwise allow any plausible oxidation state of metal
			if (setChargeOnCationicElementAppropriately(overallCharge, cationicElements.get(0), mustBeCommonOxidationState)) {
				return;
			}
		}
		
		if (!explicitStoichiometryPresent &&
				(positivelyChargedComponents.size()==1 && cationicElements.size() ==0 && negativelyChargedComponents.size() >=1 || positivelyChargedComponents.size()>=1 && negativelyChargedComponents.size() ==1 )){
			boolean success = multiplyChargedComponents(negativelyChargedComponents, positivelyChargedComponents, componentToChargeMapping, overallCharge);
			if (success){
				return;
			}
		}
		if (cationicElements.size() ==1){//e.g. magnesium monochloride [Mg2+]-->[Mg+]
			boolean success = setChargeOnCationicElementAppropriately(overallCharge, cationicElements.get(0), false);
			if (success){
				return;
			}
		}
		if (overallCharge <0){
			if (overallCharge == -1 && acetylideSpecialCase(wordRules)) {
				//e.g. acetylide dianion --> acetylide monoanion
				return;
			}
			//neutralise functionalAtoms if they are the sole cause of the negative charge and multiple molecules are present
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

	private List<Element> getMetalsThatCanBeImplicitlyCations(Element molecule) {
		List<Element> cationicElements = new ArrayList<Element>();
		List<Element> elementaryAtoms = OpsinTools.getDescendantElementsWithTagNameAndAttribute(molecule, GROUP_EL, TYPE_ATR, ELEMENTARYATOM_TYPE_VAL);
		for (Element elementaryAtom : elementaryAtoms) {
			if (elementaryAtom.getAttribute(COMMONOXIDATIONSTATESANDMAX_ATR)!=null){
				Atom metalAtom = elementaryAtom.getFrag().getFirstAtom();
				if (metalAtom.getCharge() == 0 && metalAtom.getProperty(Atom.OXIDATION_NUMBER) == null) {//if not 0 charge cannot be implicitly modified
					String[] typicalOxidationStates = elementaryAtom.getAttributeValue(COMMONOXIDATIONSTATESANDMAX_ATR).split(":")[0].split(",");
					int typicalCharge = Integer.parseInt(typicalOxidationStates[typicalOxidationStates.length-1]);
					if (typicalCharge > metalAtom.getBondCount()){
						cationicElements.add(elementaryAtom);
					}
				}
			}
		}
		return cationicElements;
	}


	/**
	 * Sets the cationicElements to the lowest typical charge as specified by the COMMONOXIDATIONSTATESANDMAX_ATR that is >= incoming valency
	 * The valency incoming to the cationicElement is taken into account e.g. phenylmagnesium chloride is [Mg+]
	 * @param cationicElements
	 * @param overallCharge
	 * @return
	 */
	private int setCationicElementsToTypicalCharge(List<Element> cationicElements, int overallCharge)  {
		for (Element cationicElement : cationicElements) {
			Fragment cationicFrag = cationicElement.getFrag();
			String[] typicalOxidationStates = cationicElement.getAttributeValue(COMMONOXIDATIONSTATESANDMAX_ATR).split(":")[0].split(",");
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
	 * Checks for tribromide/triodide and joins the ions if found
	 * @param wordRules
	 * @return
	 */
	private boolean triHalideSpecialCase(List<Element> wordRules) {
		for (Element wordRule : wordRules) {
			if (wordRule.getChildCount() == 3) {
				String value = wordRule.getAttributeValue(VALUE_ATR);
				if ("tribromide".equals(value) || "tribromid".equals(value) || "triiodide".equals(value) || "triiodid".equals(value)) {
					List<Element> groups1 = OpsinTools.getDescendantElementsWithTagName(wordRule.getChild(0), GROUP_EL);
					List<Element> groups2 = OpsinTools.getDescendantElementsWithTagName(wordRule.getChild(1), GROUP_EL);
					List<Element> groups3 = OpsinTools.getDescendantElementsWithTagName(wordRule.getChild(2), GROUP_EL);
					if (groups1.size() != 1 || groups2.size() != 1 || groups3.size() != 1) {
						throw new RuntimeException("OPSIN Bug: Unexpected trihalide representation");
					}
					Atom centralAtom = groups1.get(0).getFrag().getFirstAtom();
					Atom otherAtom1 = groups2.get(0).getFrag().getFirstAtom();
					otherAtom1.setCharge(0);
					Atom otherAtom2 = groups3.get(0).getFrag().getFirstAtom();
					otherAtom2.setCharge(0);
					state.fragManager.createBond(centralAtom, otherAtom1, 1);
					state.fragManager.createBond(centralAtom, otherAtom2, 1);
					return true;
				}
			}
		}
		return false;
	}
	
	private boolean acetylideSpecialCase(List<Element> wordRules) {
		for (Element wordRule : wordRules) {
			String value = wordRule.getAttributeValue(VALUE_ATR);
			if ("acetylide".equals(value) || "acetylid".equals(value)) {
				List<Element> groups = OpsinTools.getDescendantElementsWithTagName(wordRule, GROUP_EL);
				if (groups.size() != 1) {
					throw new RuntimeException("OPSIN Bug: Unexpected acetylide representation");
				}
				Fragment frag = groups.get(0).getFrag();
				Atom firstAtom = frag.getFirstAtom();
				if (frag.getCharge() == -2 && firstAtom.getCharge() == -1)  {
					firstAtom.addChargeAndProtons(1, 1);
					return true;
				}
			}
		}
		return false;
	}

	/**
	 * Multiplies out charged word rules to balance charge
	 * Return true if balancing was possible else false
	 * @param negativelyChargedComponents
	 * @param positivelyChargedComponents
	 * @param componentToChargeMapping
	 * @param overallCharge
	 * @return
	 * @throws StructureBuildingException
	 */
	private boolean multiplyChargedComponents(List<Element>negativelyChargedComponents, List<Element> positivelyChargedComponents, Map<Element, Integer> componentToChargeMapping, int overallCharge) throws StructureBuildingException {
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
				OpsinTools.insertAfter(componentToMultiply, state.fragManager.cloneElement(state, componentToMultiply));
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
				OpsinTools.insertAfter(negativelyChargedComponents.get(0), state.fragManager.cloneElement(state, negativelyChargedComponents.get(0)));
			}
			for (int i = (targetTotalAbsoluteCharge/positiveCharge); i >1; i--) {
				OpsinTools.insertAfter(positivelyChargedComponents.get(0), state.fragManager.cloneElement(state, positivelyChargedComponents.get(0)));
			}
		}
		return true;
	}

	private boolean componentCanBeMultiplied(Element componentToMultiply) {
		if (componentToMultiply.getAttributeValue(WORDRULE_ATR).equals(WordRule.simple.toString()) && OpsinTools.getChildElementsWithTagNameAndAttribute(componentToMultiply, WORD_EL, TYPE_ATR, WordType.full.toString()).size()>1){
			return false;//already has been multiplied e.g. dichloride
		}
		Element firstChild = componentToMultiply.getChild(0);
		while (firstChild.getChildCount() != 0){
			firstChild = firstChild.getChild(0);
		}
		if (firstChild.getName().equals(MULTIPLIER_EL) || firstChild.getName().equals(FRACTIONALMULTIPLIER_EL) ){//e.g. monochloride. Allows specification of explicit stoichiometry
			return false;
		}
		return true;
	}
	
	private boolean setChargeOnCationicElementAppropriately(int overallCharge, Element cationicElement, boolean mustBeCommonOxidationState)  {
		Atom cation = cationicElement.getFrag().getFirstAtom();
		int chargeOnCationNeeded = -(overallCharge -cation.getCharge());
		if (mustBeCommonOxidationState) {
			String[] typicalOxidationStates = cationicElement.getAttributeValue(COMMONOXIDATIONSTATESANDMAX_ATR).split(":")[0].split(",");
			for (String typicalOxidationState : typicalOxidationStates) {
				int charge = Integer.parseInt(typicalOxidationState);
				if (charge == chargeOnCationNeeded) {
					cation.setCharge(chargeOnCationNeeded);
					return true;
				}
			}
		}
		else {
			int maximumCharge = Integer.parseInt(cationicElement.getAttributeValue(COMMONOXIDATIONSTATESANDMAX_ATR).split(":")[1]);
			if (chargeOnCationNeeded >=0 && chargeOnCationNeeded <= maximumCharge){
				cation.setCharge(chargeOnCationNeeded);
				return true;
			}
		}
		return false;
	}

	private Element findRightMostGroupInWordOrWordRule(Element wordOrWordRule) throws StructureBuildingException {
		if (wordOrWordRule.getName().equals(WORDRULE_EL)){
			List<Element> words = OpsinTools.getDescendantElementsWithTagName(wordOrWordRule, WORD_EL);
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
		else if (wordOrWordRule.getName().equals(WORD_EL)){//word element can be treated just like a bracket
			return StructureBuildingMethods.findRightMostGroupInBracket(wordOrWordRule);
		}
		else{
			throw new StructureBuildingException("OPSIN bug: expected word or wordRule");
		}
	}

	/**
	 * Nasty special cases:
	 * Oxido and related groups can act as O= or even [O-][N+]
	 * This nasty behaviour is in generated ChemDraw names and is supported by most nameToStructure tools so it is supported here
	 * Acting as O= notably is often correct behaviour for inorganics
	 * 
	 * Methionine (and the like) when substituted at the sulfur/selenium/tellurium are implicitly positively charged
	 * 
	 * purine nucleosides/nucleotides are implicitly positively charged when 7-substituted
	 * @param groups
	 */
	private void processSpecialCases(List<Element> groups)  {
		for (Element group : groups) {
			String subType = group.getAttributeValue(SUBTYPE_ATR);
			if (OXIDOLIKE_SUBTYPE_VAL.equals(subType)){
				Atom oxidoAtom = group.getFrag().getFirstAtom();
				if (oxidoAtom.getBondCount() != 1) {
					continue;
				}
				Atom connectedAtom = oxidoAtom.getFirstBond().getOtherAtom(oxidoAtom);
				ChemEl chemEl = connectedAtom.getElement();
				if (checkForConnectedOxo(connectedAtom)){//e.g. not oxido(trioxo)ruthenium
					continue;
				}
				if (ELEMENTARYATOM_TYPE_VAL.equals(connectedAtom.getFrag().getType()) ||
						((chemEl == ChemEl.S || chemEl == ChemEl.P) && connectedAtom.getCharge() ==0 && ValencyChecker.checkValencyAvailableForBond(connectedAtom, 1))){
					oxidoAtom.neutraliseCharge();
					oxidoAtom.getFirstBond().setOrder(2);
				}
				else if (chemEl == ChemEl.N && connectedAtom.getCharge()==0){
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
			else if (AMINOACID_TYPE_VAL.equals(group.getAttributeValue(TYPE_ATR))) {
				for (Atom atom : group.getFrag().getAtomList()) {
					if (atom.getElement().isChalcogen() && atom.getElement() != ChemEl.O &&
							atom.getBondCount() == 3 && atom.getIncomingValency() == 3 && atom.getCharge() == 0) {
						atom.addChargeAndProtons(1, 1);
					}
				}
			}
			else if (BIOCHEMICAL_SUBTYPE_VAL.equals(subType)) {
				Fragment frag = group.getFrag();
				Atom atom = frag.getAtomByLocant("7");
				if (atom != null) {
					String groupName = group.getValue();
					if (groupName.equals("adenosin") || groupName.equals("guanosin") || groupName.equals("inosin") || 
							groupName.equals("thioinosin") || groupName.equals("xanthosin") || groupName.equals("nucleocidin") ||
							groupName.contains("adenylic") || groupName.contains("guanylic") || groupName.contains("inosinic") ||
							groupName.contains("xanthylic") || groupName.endsWith("adenylyl") || groupName.endsWith("adenosyl") ||
							groupName.endsWith("guanylyl") || groupName.endsWith("guanosyl") || groupName.endsWith("inosinylyl") ||
							groupName.endsWith("inosyl") || groupName.endsWith("xanthylyl") || groupName.endsWith("xanthosyl")) {
						if (atom.getElement() == ChemEl.N && atom.hasSpareValency() && atom.getBondCount() == 3 && atom.getIncomingValency() == 3 && atom.getCharge() == 0) {
							atom.addChargeAndProtons(1, 1);
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
	private boolean checkForConnectedOxo(Atom atom) {
		List<Bond> bonds = atom.getBonds();
		for (Bond bond : bonds) {
			Atom connectedAtom;
			if (bond.getFromAtom() == atom){
				connectedAtom = bond.getToAtom();
			}
			else{
				connectedAtom = bond.getFromAtom();
			}
			Element correspondingEl = connectedAtom.getFrag().getTokenEl();
			if (correspondingEl.getValue().equals("oxo")){
				return true;
			}
		}
		return false;
	}
	

	/**
	 * Sets the charge according to the oxidation number if the oxidation number atom property has been set
	 * @param groups
	 * @throws StructureBuildingException 
	 */
	private void processOxidationNumbers(List<Element> groups) throws StructureBuildingException  {
		for (Element group : groups) {
			if (ELEMENTARYATOM_TYPE_VAL.equals(group.getAttributeValue(TYPE_ATR))){
				Atom atom = group.getFrag().getFirstAtom();
				if (atom.getProperty(Atom.OXIDATION_NUMBER)!=null){
					List<Atom> neighbours = atom.getAtomNeighbours();
					int chargeThatWouldFormIfLigandsWereRemoved =0;
					for (Atom neighbour : neighbours) {
						Element neighbourEl = neighbour.getFrag().getTokenEl();
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
	 * @param molecule
	 * @param uniFrag
	 * @throws StructureBuildingException
	 */
	private void processStereochemistry(Element molecule, Fragment uniFrag) throws StructureBuildingException {
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
			StereoAnalyser stereoAnalyser = new StereoAnalyser(uniFrag);
			Map<Atom, StereoCentre> atomStereoCentreMap = new HashMap<Atom, StereoCentre>();//contains all atoms that are stereo centres with a mapping to the corresponding StereoCentre object
			List<StereoCentre> stereoCentres = stereoAnalyser.findStereoCentres();
			for (StereoCentre stereoCentre : stereoCentres) {
				atomStereoCentreMap.put(stereoCentre.getStereoAtom(),stereoCentre);
			}
			Map<Bond, StereoBond> bondStereoBondMap = new HashMap<Bond, StereoBond>();
			List<StereoBond> stereoBonds = stereoAnalyser.findStereoBonds();
			for (StereoBond stereoBond : stereoBonds) {
				Bond b = stereoBond.getBond();
				if (FragmentTools.notIn6MemberOrSmallerRing(b)){
					bondStereoBondMap.put(b, stereoBond);
				}
			}
			StereochemistryHandler stereoChemistryHandler = new StereochemistryHandler(state, atomStereoCentreMap, bondStereoBondMap);
			stereoChemistryHandler.applyStereochemicalElements(stereoChemistryEls);
			stereoChemistryHandler.removeRedundantStereoCentres(atomsWithPreDefinedAtomParity, bondsWithPreDefinedBondStereo);
		}
	}

	/**
	 * Finds stereochemistry els in a recursive right to left manner.
	 * Within the same scope though stereochemistry els are found left to right
	 * @param parentEl
	 * @return
	 */
	private List<Element> findStereochemistryElsInProcessingOrder(Element parentEl) {
		List<Element> matchingElements = new ArrayList<Element>();
		List<Element> children =parentEl.getChildElements();
		List<Element> stereochemistryElsAtThisLevel = new ArrayList<Element>();
		for (int i = children.size()-1; i >=0; i--) {
			Element child = children.get(i);
			if (child.getName().equals(STEREOCHEMISTRY_EL)){
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
	
	private void convertOutAtomsToAttachmentAtoms(Fragment uniFrag) throws StructureBuildingException {
		int outAtomCount = uniFrag.getOutAtomCount();
		for (int i = outAtomCount -1; i >=0; i--) {
			OutAtom outAtom = uniFrag.getOutAtom(i);
			uniFrag.removeOutAtom(i);
			Atom rGroup = state.fragManager.createAtom(ChemEl.R, uniFrag);
			state.fragManager.createBond(outAtom.getAtom(), rGroup, outAtom.getValency());
		}
	}

	/**
	 * Returns the atom corresponding to position i in the outAtoms list
	 * If the outAtom is not set explicitly a suitable atom will be found
	 * @param buildResults 
	 * @param i index
	 * @return atom
	 * @throws StructureBuildingException
	 */
	private Atom getOutAtomTakingIntoAccountWhetherSetExplicitly(BuildResults buildResults, int i) throws StructureBuildingException {
		OutAtom outAtom = buildResults.getOutAtom(i);
		if (outAtom.isSetExplicitly()){
			return outAtom.getAtom();
		}
		else{
			return findAtomForUnlocantedRadical(state, outAtom.getAtom().getFrag(), outAtom);
		}
	}
}
