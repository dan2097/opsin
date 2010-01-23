package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.Stack;
import java.util.regex.Pattern;

import uk.ac.cam.ch.wwmm.opsin.ParseWord.WordType;
import uk.ac.cam.ch.wwmm.opsin.WordRules.WordRule;

import nu.xom.Attribute;
import nu.xom.Element;
import nu.xom.Elements;
import static  uk.ac.cam.ch.wwmm.opsin.StructureBuildingMethods.*;

/**Constructs a single OPSIN fragment which describes the molecule from the postprocessor results.
 *
 * @author ptc24/dl387
 *
 */
class StructureBuilder {

	private Pattern matchDigits = Pattern.compile("\\d+");
	private Pattern matchComma =Pattern.compile(",");
	private Pattern matchNumericLocant =Pattern.compile("\\d+[a-z]?'*");

	/**	Builds a molecule as a Fragment based on preStructurebuilder output.
	 * @param state
	 * @param molecule The preStructurebuilderr output.
	 * @return A single Fragment - the built molecule.
	 * @throws StructureBuildingException If the molecule won't build - there may be many reasons.
	 */
	Fragment buildFragment(BuildState state, Element molecule) throws StructureBuildingException {
		Elements wordRules = molecule.getChildElements("wordRule");
		if (wordRules.size()==0){
			throw new StructureBuildingException("Molecule contains no words!?");
		}
		Stack<Element> wordRuleStack = new Stack<Element>();
		for (int i = wordRules.size() -1; i >=0; i--) {
			wordRuleStack.add(wordRules.get(i));
		}
		
		List<Element> wordRulesVisited = new ArrayList<Element>();
		while (wordRuleStack.size()>0) {
			Element nextWordRuleEl = wordRuleStack.peek();//just has a look what's next
			if(!wordRulesVisited.contains(nextWordRuleEl)){
				wordRulesVisited.add(nextWordRuleEl);
				Elements wordRuleChildren = nextWordRuleEl.getChildElements("wordRule");
				if (wordRuleChildren.size()!=0){//nested word rules
					for (int i = wordRuleChildren.size() -1; i >=0; i--) {
						wordRuleStack.add(wordRuleChildren.get(i));
					}
					continue;
				}
			}
			Element currentWordRuleEl = wordRuleStack.pop();
			WordRule wordRule = WordRule.valueOf(currentWordRuleEl.getAttributeValue("wordRule"));
			List<Element> words = XOMTools.getChildElementsWithTagNames(currentWordRuleEl, new String[]{"word","wordRule"});
			state.currentWordRule =wordRule;
			if(wordRule == WordRule.simple) {
				for (Element word : words) {
					if (!word.getLocalName().equals("word") || !word.getAttributeValue("type").equals(WordType.full.toString())){
						throw new StructureBuildingException("OPSIN bug: Unexpected contents of 'simple' wordRule");
					}
					resolveWordOrBracket(state, word);
				}
			}
			else if (wordRule == WordRule.acid){
				buildAcid(state,words);//ethanoic acid
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
			else if(wordRule == WordRule.polymer) {
				buildPolymer(state, words);
			}
			else{
				throw new StructureBuildingException("Unknown Word Rule");
			}
		}

		state.fragManager.convertSpareValenciesToDoubleBonds();
		state.fragManager.checkValencies();
		makeHydrogensExplicit(state);

		Fragment uniFrag = state.fragManager.getUnifiedFragment();
		List<Element> stereoChemistryEls = XOMTools.getDescendantElementsWithTagName(molecule, "stereoChemistry");
		if (stereoChemistryEls.size() >0){
			StereochemistryHandler.processStereochemicalElements(state, molecule, uniFrag, stereoChemistryEls);
		}

		//adds Xe group at all atoms which have unused outIDs
		//note that SMILES generated by CDK is not always correct
//		BuildResults moleculeBuildResults =  new BuildResults(state, molecule);
//		for (int i = 0; i < moleculeBuildResults.getOutIDCount(); i++) {
//			Atom outAtom =moleculeBuildResults.getOutAtomTakingIntoAccountWhetherSetExplicitly(i);
//			Fragment rGroup =state.fragManager.buildSMILES("[Xe]");
//			state.fragManager.attachFragments(outAtom, rGroup.getAtomByIDOrThrow(rGroup.getIdOfFirstAtom()), moleculeBuildResults.getOutID(i).valency);
//		}
//		if(wordRule.equals("polymer")){
//			List<Atom> atomList = uniFrag.getAtomList();
//			for (Atom atom : atomList) {
//				if  (atom.getElement().equals("No") || atom.getElement().equals("Lr")){
//					atom.setElement("R");
//				}
//			}
//		}

		if (uniFrag.getOutIDs().size()>0 || uniFrag.getInIDs().size()>0){
			throw new StructureBuildingException("Radicals are currently set to not convert to structures");
		}
		return uniFrag;
	}

	private void buildAcid(BuildState state, List<Element> words) throws StructureBuildingException {
		if (words.size()!=2 || !words.get(1).getAttributeValue("type").equals(WordType.functionalTerm.toString())){
			throw new StructureBuildingException("functionalTerm word acid missing");
		}
		resolveWordOrBracket(state, words.get(0));
	}

	private void buildEster(BuildState state, List<Element> words) throws StructureBuildingException {
		int wordIndice = 0;
		Element currentWord=words.get(wordIndice);
		BuildResults substituentsBr = new BuildResults();
		while (currentWord.getAttributeValue("type").equals(WordType.substituent.toString())){
			resolveWordOrBracket(state, currentWord);
			BuildResults substituentBr = new BuildResults(state, currentWord);
			if (substituentBr.getOutIDCount() ==1){//TODO add support for locanted terepthaloyl
				String locantForSubstituent = getLocantForSubstituent(currentWord);
				if (locantForSubstituent!=null){
					substituentBr.getFirstOutID().locant=locantForSubstituent;//indexes which functional atom to connect to when there is a choice. Also can disambiguate which atom is a S in things like thioates
				}
			}
			else if (substituentBr.getOutIDCount() ==0){
				throw new StructureBuildingException("Substituent was expected to have at least one outID");
			}
			substituentsBr.mergeBuildResults(substituentBr);
			currentWord=words.get(++wordIndice);
		}

		if (wordIndice == words.size() || !words.get(wordIndice).getAttributeValue("type").equals(WordType.full.toString())){
			throw new StructureBuildingException("Full word not found where full word expected: missing ate group in ester");
		}

		BuildResults ateGroups = new BuildResults();
		for (; wordIndice < words.size(); wordIndice++) {
			Element word =words.get(wordIndice);
			if (word.getAttributeValue("type").equals(WordType.full.toString())){
				resolveWordOrBracket(state, word);
				BuildResults ateBR = new BuildResults(state, word);
				ateGroups.mergeBuildResults(ateBR);
			}
			else{
				throw new StructureBuildingException("Non full word found where only full words were expected");
			}
		}
		int esterIdCount = ateGroups.getFunctionalIDCount();
		int outIDCount =substituentsBr.getOutIDCount();
		if (outIDCount > esterIdCount){
			throw new StructureBuildingException("There are more radicals in the substituents(" + outIDCount +") than there are places to form esters("+esterIdCount+")");
		}
		for(int i=0; i< outIDCount; i++) {
			Atom ateAtom;
			if (substituentsBr.getFirstOutID().locant!=null){
				ateAtom =determineFunctionalAtomToUse(substituentsBr.getFirstOutID().locant, ateGroups);
			}
			else{
				ateAtom =ateGroups.getFunctionalAtom(0);
				ateGroups.removeFunctionalID(0);
			}
			state.fragManager.createBond(ateAtom,substituentsBr.getOutAtomTakingIntoAccountWhetherSetExplicitly(0), 1);
			substituentsBr.removeOutID(0);
			ateAtom.setCharge(0);
		}
	}



	private void buildDiValentFunctionalGroup(BuildState state, List<Element> words) throws StructureBuildingException {
		int wordIndice = 0;
		if (!words.get(wordIndice).getAttributeValue("type").equals(WordType.substituent.toString())) {
			throw new StructureBuildingException("word: " +wordIndice +" was expected to be a substituent");
		}
		if (words.get(wordIndice +1).getAttributeValue("type").equals(WordType.functionalTerm.toString())) {//e.g. methyl sulfoxide rather than dimethyl sulfoxide
			Element clone = state.fragManager.cloneElement(state, words.get(0));
			XOMTools.insertAfter(words.get(0), clone);
			words = OpsinTools.elementsToElementArrayList(((Element)words.get(0).getParent()).getChildElements());
		}
		resolveWordOrBracket(state, words.get(wordIndice));
		BuildResults substituent1 =new BuildResults(state, words.get(wordIndice));
		if (substituent1.getOutIDCount()!=1){
			throw new StructureBuildingException("Expected one outID. Found " + substituent1.getOutIDCount() );
		}
		if (substituent1.getOutID(0).valency !=1){
			throw new StructureBuildingException("OutID has unexpected valency. Expected 1. Actual: " + substituent1.getOutID(0).valency);
		}
		wordIndice++;
		resolveWordOrBracket(state, words.get(wordIndice));
		BuildResults substituent2 =new BuildResults(state, words.get(wordIndice));
		if (substituent2.getOutIDCount()!=1){
			throw new StructureBuildingException("Expected one outID. Found " + substituent2.getOutIDCount() );
		}
		if (substituent2.getOutID(0).valency !=1){
			throw new StructureBuildingException("OutID has unexpected valency. Expected 1. Actual: " + substituent2.getOutID(0).valency);
		}
		wordIndice++;
		if (words.get(wordIndice) ==null || !words.get(wordIndice).getAttributeValue("type").equals(WordType.functionalTerm.toString())) {
			throw new StructureBuildingException("word: " +wordIndice +" was expected to be a functionalTerm");
		}
		List<Element> functionalGroup = XOMTools.getDescendantElementsWithTagName(words.get(wordIndice), "functionalGroup");
		if (functionalGroup.size()!=1){
			throw new StructureBuildingException("Unexpected number of functionalGroups found, could be a bug in OPSIN's grammar");
		}
		String smilesOfGroup = functionalGroup.get(0).getAttributeValue("value");
		Fragment diValentGroup =state.fragManager.buildSMILES(smilesOfGroup, "simpleGroup", "functionalGroup", "none");

		Atom outAtom =substituent1.getOutAtomTakingIntoAccountWhetherSetExplicitly(0);
		substituent1.removeOutID(0);
		if (diValentGroup.getOutIDs().size()==1){//c.f. peroxide where it is a linker
			state.fragManager.createBond(outAtom, diValentGroup.getAtomByIDOrThrow(diValentGroup.getOutID(0).id), 1);
			diValentGroup.removeOutID(0);
		}
		else{
			state.fragManager.createBond(outAtom, diValentGroup.getAtomByIDOrThrow(diValentGroup.getIdOfFirstAtom()), 1);
		}
		outAtom = substituent2.getOutAtomTakingIntoAccountWhetherSetExplicitly(0);
		substituent2.removeOutID(0);
		state.fragManager.createBond(outAtom, diValentGroup.getAtomByIDOrThrow(diValentGroup.getIdOfFirstAtom()), 1);
	}

	private void buildMonovalentFunctionalGroup(BuildState state, List<Element> words) throws StructureBuildingException {
		resolveWordOrBracket(state, words.get(0));
		List<Element> groups = XOMTools.getDescendantElementsWithTagName(words.get(0), "group");
		for (Element group : groups) {//replaces outIDs with valency greater than 1 with multiple outIDs; e.g. ylidene -->diyl
			Fragment frag = state.xmlFragmentMap.get(group);
			for (int i = frag.getOutIDs().size()-1; i>=0; i--) {
				OutID outID =frag.getOutID(i);
				if (outID.valency>1){
					for (int j = 2; j <= outID.valency; j++) {
						frag.addOutID(outID.id, 1, outID.setExplicitly);
					}
					outID.valency=1;
				}
			}
		}
		BuildResults substituentBR = new BuildResults(state, words.get(0));

		List<Fragment> functionalGroupFragments = new ArrayList<Fragment>();
		for (int i=1; i<words.size(); i++ ) {
			Element functionalGroupWord =words.get(i);
			List<Element> functionalGroups = XOMTools.getDescendantElementsWithTagName(functionalGroupWord,"functionalGroup");
			if (functionalGroups.size()!=1){
				throw new StructureBuildingException("Expected exactly 1 functionalGroup. Found " + functionalGroups.size());
			}
			
			Fragment monoValentFunctionGroup =state.fragManager.buildSMILES(functionalGroups.get(0).getAttributeValue("value"), "simpleGroup", "functionalGroup", "none");
			if (functionalGroups.get(0).getAttributeValue("type").equals("monoValentStandaloneGroup")){
				Atom ideAtom = monoValentFunctionGroup.getAtomByIDOrThrow(monoValentFunctionGroup.getDefaultInID());
				ideAtom.setCharge(ideAtom.getCharge()+1);//e.g. make cyanide charge netural
			}
			Element possibleMultiplier = (Element) XOMTools.getPreviousSibling(functionalGroups.get(0));
			functionalGroupFragments.add(monoValentFunctionGroup);
			if (possibleMultiplier!=null){
				int multiplierValue = Integer.parseInt(possibleMultiplier.getAttributeValue("value"));
				for (int j = 1; j < multiplierValue; j++) {
					functionalGroupFragments.add(state.fragManager.copyAndRelabel(monoValentFunctionGroup));
				}
				possibleMultiplier.detach();
			}
		}
		
		int numberOfOutIDs =substituentBR.getOutIDCount();
		if (numberOfOutIDs > functionalGroupFragments.size()){//something like isophthaloyl chloride (more precisely written isophthaloyl dichloride)
			if (functionalGroupFragments.size()!=1){
				throw new StructureBuildingException("Incorrect number of functional groups found to balance outIDs");
			}
			Fragment monoValentFunctionGroup = functionalGroupFragments.get(0);
			for (int j = 1; j < numberOfOutIDs; j++) {
				functionalGroupFragments.add(state.fragManager.copyAndRelabel(monoValentFunctionGroup));
			}
		}
		else if (functionalGroupFragments.size() > numberOfOutIDs){
			throw new StructureBuildingException("There are more function groups to attach than there are positions to attach them to!");
		}
		for (int i = 0; i < numberOfOutIDs; i++) {
			Fragment ideFrag =functionalGroupFragments.get(i);
			Atom ideAtom = ideFrag.getAtomByIDOrThrow(ideFrag.getDefaultInID());
			Atom subAtom=substituentBR.getOutAtomTakingIntoAccountWhetherSetExplicitly(0);
			state.fragManager.createBond(ideAtom, subAtom, 1);
			substituentBR.removeOutID(0);
		}
	}

	private void buildFunctionalClassEster(BuildState state, List<Element> words) throws StructureBuildingException {
		if (!words.get(0).getAttributeValue("type").equals(WordType.full.toString())){
			throw new StructureBuildingException("Don't alter wordRules.xml without checking the consequences!");
		}
		resolveWordOrBracket(state, words.get(0));//the group
		BuildResults acidBr = new BuildResults(state, words.get(0));
		if (!words.get(1).getAttributeValue("type").equals(WordType.functionalTerm.toString())){//acid
			throw new StructureBuildingException("Don't alter wordRules.xml without checking the consequences!");
		}
		if (acidBr.getFunctionalIDCount()==0){
			throw new StructureBuildingException("No functionalIds detected!");
		}

		int i=2;
		Element currentWord = words.get(i);
		while (currentWord.getAttributeValue("type").equals(WordType.substituent.toString())){
			if (acidBr.getFunctionalIDCount()==0){
				throw new StructureBuildingException("Insufficient functionalIDs on acid");
			}
			resolveWordOrBracket(state, currentWord);
			BuildResults substituentBr = new BuildResults(state, currentWord);
			if (substituentBr.getOutIDCount() ==1){
				String locantForSubstituent = getLocantForSubstituent(currentWord);
				Atom functionalAtom;
				if (locantForSubstituent!=null){
					functionalAtom =determineFunctionalAtomToUse(locantForSubstituent, acidBr);
				}
				else{
					functionalAtom =acidBr.getFunctionalAtom(0);
					acidBr.removeFunctionalID(0);
				}
				if (substituentBr.getOutID(0).valency!=1){
					throw new StructureBuildingException("Substituent was expected to have only have an outgoing valency of 1");
				}
				state.fragManager.createBond(functionalAtom,substituentBr.getOutAtomTakingIntoAccountWhetherSetExplicitly(0), 1);
				substituentBr.removeOutID(0);
			}
			else {
				throw new StructureBuildingException("Substituent was expected to have one outID");
			}
			currentWord=words.get(++i);
		}
		if (!words.get(i++).getAttributeValue("type").equals(WordType.functionalTerm.toString())){//ester
			throw new StructureBuildingException("Number of words different to expectations; did not find ester");
		}
	}

	private void buildAmide(BuildState state, List<Element> words) throws StructureBuildingException {
		if (!words.get(0).getAttributeValue("type").equals(WordType.full.toString())){
			throw new StructureBuildingException("Don't alter wordRules.xml without checking the consequences!");
		}
		resolveWordOrBracket(state, words.get(0));//the group
		BuildResults acidBr = new BuildResults(state, words.get(0));
		if (acidBr.getFunctionalIDCount()==0){
			throw new StructureBuildingException("No functionalIds detected!");
		}
		int wordIndice =1;
		if (words.get(wordIndice).getAttributeValue("type").equals(WordType.functionalTerm.toString())){//"acid"
			wordIndice++;
		}
		
		if (words.get(wordIndice).getAttributeValue("type").equals(WordType.functionalTerm.toString())){//"amide"
			for (int i =0; i < acidBr.getFunctionalIDCount(); i++) {
				Atom functionalAtom = acidBr.getFunctionalAtom(i);
				if (!functionalAtom.getElement().equals("O")){
					throw new StructureBuildingException("Expected oxygen functional atom found:" + functionalAtom.getElement());
				}
				functionalAtom.setElement("N");
				functionalAtom.replaceLocant("N" +StringTools.multiplyString("'", i));
			}
		}
		else if (words.get(wordIndice).getAttributeValue("type").equals(WordType.full.toString())){//substituentamide
			resolveWordOrBracket(state, words.get(wordIndice));//the substituted amide ([N-]) group
			List<Element> root = XOMTools.getDescendantElementsWithTagName(words.get(wordIndice), "root");
			if (root.size()!=1){
				throw new StructureBuildingException("Cannot find root element");
			}
			Fragment amide = state.xmlFragmentMap.get(root.get(0).getFirstChildElement("group"));
			if (acidBr.getFunctionalIDCount()!=1){
				throw new StructureBuildingException("Formation of amide is unexpectedly ambiguous!");
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
			acidBr.removeFunctionalID(0);
			amideNitrogen.setCharge(0);
			state.fragManager.replaceTerminalAtomWithFragment(functionalAtom, amideNitrogen);
		}
		if (wordIndice +1 < words.size()){
			throw new StructureBuildingException("More words than expected when applying amide word rule!");
		}
	}

	private void buildGlycol(BuildState state, List<Element> words) throws StructureBuildingException {
		int wordIndice  = 0;
		resolveWordOrBracket(state, words.get(wordIndice));//the group
		BuildResults theDiRadical = new BuildResults(state, words.get(wordIndice));
		if (theDiRadical.getOutIDCount()!=2){
			throw new StructureBuildingException("Glycol class names (e.g. ethylene glycol) expect two outIDs. Found: " + theDiRadical.getOutIDCount() );
		}
		if (wordIndice +1 >= words.size() || !words.get(wordIndice+1).getAttributeValue("type").equals(WordType.functionalTerm.toString())){
			throw new StructureBuildingException("Glycol functionalTerm word expected");
		}
		for (int i = theDiRadical.getOutIDCount() -1; i >=0 ; i--) {
			Atom outAtom =theDiRadical.getOutAtomTakingIntoAccountWhetherSetExplicitly(0);
			Fragment glycol =state.fragManager.buildSMILES("O", "glycol", "none");
			if (theDiRadical.getOutID(0).valency !=1){
				throw new StructureBuildingException("OutID has unexpected valency. Expected 1. Actual: " + theDiRadical.getOutID(0).valency);
			}
			theDiRadical.removeOutID(0);
			state.fragManager.createBond(outAtom, glycol.getAtomByIDOrThrow(glycol.getIdOfFirstAtom()), 1);
		}
	}
	
	private void buildOxide(BuildState state, List<Element> words) throws StructureBuildingException {
		resolveWordOrBracket(state, words.get(0));//the group
		List<Fragment> oxideFragments = new ArrayList<Fragment>();
		List<String> locantsForOxide =new ArrayList<String>();//often not specified
		int numberOfOxygenToAdd =1;
		List<Element> multipliers =XOMTools.getDescendantElementsWithTagName(words.get(1),"multiplier");
		if (multipliers.size() >1){
			throw new StructureBuildingException("Expected 0 or 1 multiplier found: " + multipliers.size());
		}
		if (multipliers.size()==1){
			numberOfOxygenToAdd = Integer.parseInt(multipliers.get(0).getAttributeValue("value"));
			multipliers.get(0).detach();
		}
		List<Element> functionalClass =XOMTools.getDescendantElementsWithTagName(words.get(1), "group");
		if (functionalClass.size()!=1){
			throw new StructureBuildingException("Expected 1 group element found: " + functionalClass.size());
		}
		String smilesReplacement = functionalClass.get(0).getAttributeValue("value");
		String labels =  functionalClass.get(0).getAttributeValue("labels");
		for (int i = 0; i < numberOfOxygenToAdd; i++) {
			oxideFragments.add(state.fragManager.buildSMILES(smilesReplacement, "", labels));
		}
		List<Element> locantEls =XOMTools.getDescendantElementsWithTagName(words.get(1), "locant");
		if (locantEls.size() >1){
			throw new StructureBuildingException("Expected 0 or 1 locant elements found: " + locantEls.size());
		}
		if (locantEls.size()==1){
			String[] locants = matchComma.split(StringTools.removeDashIfPresent(locantEls.get(0).getValue()));
			for (String locant : locants) {
				locantsForOxide.add(locant);
			}
			locantEls.get(0).detach();
		}
		if (!locantsForOxide.isEmpty() && locantsForOxide.size()!=oxideFragments.size()){
			throw new StructureBuildingException("Mismatch between number of locants and number of oxides specified");
		}
		Element rightMostGroup =findRightMostGroupInBracket(words.get(0));
		List<Fragment> orderedPossibleFragments = new ArrayList<Fragment>();//In preference suffixes are substituted onto e.g. acetonitrile oxide
		Elements suffixEls = ((Element)rightMostGroup.getParent()).getChildElements("suffix");
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
						if (!atom.getElement().equals("C")){
							formAppropriateBondToOxideAndAdjustCharges(state, atom, oxideAtom);
							continue mainLoop;
						}
					}
				}
				//No heteroatoms could be found. Perhaps it's supposed to be something like styrene oxide
				Set<Bond> bondSet = groupToModify.getBondSet();//looking for double bond
				for (Bond bond : bondSet) {
					if (bond.getOrder()==2){
						bond.setOrder(1);
						state.fragManager.createBond(bond.getFromAtom(), oxideAtom, 1);
						state.fragManager.createBond(bond.getToAtom(), oxideAtom, 1);
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
			atomToAddOxideTo.addChargeAndProtons(1, 1);
			oxideAtom.addChargeAndProtons(-1, -1);
			if (!ValencyChecker.checkValencyAvailableForBond(atomToAddOxideTo, 1)){
				throw new StructureBuildingException("Oxide appeared to refer to an atom that has insufficent valency to accept the addition of oxygen");
			}
			state.fragManager.createBond(atomToAddOxideTo, oxideAtom, 1);
		}
	}

	private void buildCarbonylDerivative(BuildState state, List<Element> words) throws StructureBuildingException {
		resolveWordOrBracket(state, words.get(0));//the group
		BuildResults moleculeToModify = new BuildResults(state, words.get(0));//the group which will be modified
		List<Fragment> replacementFragments = new ArrayList<Fragment>();
		List<String> locantForFunctionalTerm =new ArrayList<String>();//usually not specified
		if (!words.get(1).getAttributeValue("type").equals(WordType.functionalTerm.toString())){//e.g. acetone O-ethyloxime or acetone 1-chloro-1-methylhydrazone
			for (int i = 1; i < words.size(); i++) {
				Fragment frag = state.xmlFragmentMap.get(findRightMostGroupInBracket(words.get(i)));
				replacementFragments.add(frag);
				Elements children =words.get(i).getChildElements();
				if (children.size()==1 && children.get(0).getLocalName().equals("bracket") && children.get(0).getAttribute("locant")!=null){
					locantForFunctionalTerm.add(children.get(0).getAttributeValue("locant"));
				}
				else if (children.size()==2 && children.get(0).getAttribute("locant")!=null ){
					String locant =children.get(0).getAttributeValue("locant");
					if (children.get(1).getLocalName().equals("root") && !frag.hasLocant(locant) && matchNumericLocant.matcher(locant).matches()){
						locantForFunctionalTerm.add(children.get(0).getAttributeValue("locant"));
						children.get(0).removeAttribute(children.get(0).getAttribute("locant"));
					}
				}
			}
		}
		else{//e.g. butan-2,3-dione dioxime or hexan2,3-dione 2-oxime
			int numberOfCarbonylReplacements =1;
			List<Element> multipliers =XOMTools.getDescendantElementsWithTagName(words.get(1),"multiplier");
			if (multipliers.size() >1){
				throw new StructureBuildingException("Expected 0 or 1 multiplier found: " + multipliers.size());
			}
			if (multipliers.size()==1){
				numberOfCarbonylReplacements = Integer.parseInt(multipliers.get(0).getAttributeValue("value"));
				multipliers.get(0).detach();
			}
			List<Element> functionalClass =XOMTools.getDescendantElementsWithTagName(words.get(1), "group");
			if (functionalClass.size()!=1){
				throw new StructureBuildingException("Expected 1 group element found: " + functionalClass.size());
			}
			String smilesReplacement = functionalClass.get(0).getAttributeValue("value");
			String labels =  functionalClass.get(0).getAttributeValue("labels");
			for (int i = 0; i < numberOfCarbonylReplacements; i++) {
				replacementFragments.add(state.fragManager.buildSMILES(smilesReplacement, "", labels));
			}
			List<Element> locantEls =XOMTools.getDescendantElementsWithTagName(words.get(1), "locant");
			if (locantEls.size() >1){
				throw new StructureBuildingException("Expected 0 or 1 locant elements found: " + locantEls.size());
			}
			if (locantEls.size()==1){
				String[] locants = matchComma.split(StringTools.removeDashIfPresent(locantEls.get(0).getValue()));
				for (String locant : locants) {
					locantForFunctionalTerm.add(locant);
				}
				locantEls.get(0).detach();
			}
		}
		if (!locantForFunctionalTerm.isEmpty() && locantForFunctionalTerm.size()!=replacementFragments.size()){
			throw new StructureBuildingException("Mismatch between number of locants and number of carbonyl replacements");
		}
		List<Atom> matches = new ArrayList<Atom>();
		List<Fragment> fragmentsInGroupToModify = new ArrayList<Fragment>(moleculeToModify.getFragments());
		for (int i= fragmentsInGroupToModify.size()-1; i>=0; i-- ){//find all carbonyl oxygen
			Fragment frag = fragmentsInGroupToModify.get(i);
			List<Atom> atomList = frag.getAtomList();
			for (Atom atom : atomList) {
				if (atom.getElement().equals("O") && atom.getCharge()==0){
					List<Atom> neighbours =atom.getAtomNeighbours();
					if (neighbours.size()==1){
						if (neighbours.get(0).getElement().equals("C")){
							if (!locantForFunctionalTerm.isEmpty()){
								Atom numericLocantAtomConnectedToCarbonyl = OpsinTools.depthFirstSearchForAtomWithNumericLocant(atom);
								if (numericLocantAtomConnectedToCarbonyl!=null){//could be the carbon of the carbonyl or the ring the carbonyl connects to in say a carbaldehyde
									boolean matchesLocant = false;
									for (String locant : locantForFunctionalTerm) {
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
							Bond b = frag.findBondOrThrow(atom, neighbours.get(0));
							if (b.getOrder()==2){
								matches.add(atom);
							}
						}
					}
				}
			}
		}
		if (matches.size() < replacementFragments.size()){
			throw new StructureBuildingException("Insufficient carbonyl groups found!");
		}
		for (int i = 0; i < replacementFragments.size(); i++) {
			Atom atomToBeReplaced =matches.get(i);//the oxygen of the carbonyl
			Fragment replacementFrag = replacementFragments.get(i);
			List<Atom> atomList = replacementFrag.getAtomList();
			Atom atomToReplaceCarbonylOxygen = atomList.get(atomList.size()-1);
			Atom numericLocantAtomConnectedToCarbonyl = OpsinTools.depthFirstSearchForAtomWithNumericLocant(atomToBeReplaced);
			if (numericLocantAtomConnectedToCarbonyl!=null){
				atomList.get(0).addLocant(atomList.get(0).getElement() + numericLocantAtomConnectedToCarbonyl.getFirstLocant());
			}
			if (!words.get(1).getAttributeValue("type").equals(WordType.functionalTerm.toString())){
				resolveWordOrBracket(state, words.get(1 +i));
			}
			state.fragManager.replaceTerminalAtomWithFragment(atomToBeReplaced, atomToReplaceCarbonylOxygen);
		}
	}

	private void buildPolymer(BuildState state, List<Element> words) throws StructureBuildingException {
		if (words.size()!=2){
			throw new StructureBuildingException("Currently unsupported polymer name type");
		}
		Element polymer = words.get(1);
		resolveWordOrBracket(state, polymer);
		BuildResults polymerBr = new BuildResults(state, polymer);
		if (polymerBr.getOutIDCount() ==2 && polymerBr.getInIDCount()==0){
			Atom inAtom =polymerBr.getOutAtomTakingIntoAccountWhetherSetExplicitly(0);
			Fragment rGroup =state.fragManager.buildSMILES("[No]");//TODO stop using actual atoms (confuses E/Z!)
			state.fragManager.createBond(inAtom, rGroup.getAtomByIDOrThrow(rGroup.getIdOfFirstAtom()), polymerBr.getOutID(0).valency);
			Atom outAtom =polymerBr.getOutAtomTakingIntoAccountWhetherSetExplicitly(1);
			rGroup =state.fragManager.buildSMILES("[Lr]");
			state.fragManager.createBond(outAtom, rGroup.getAtomByIDOrThrow(rGroup.getIdOfFirstAtom()), polymerBr.getOutID(1).valency);
			polymerBr.removeAllOutIDs();
		}
		else{
			throw new StructureBuildingException("Polymer building failed: Two termini were not found; Expected 2 outIDs, found: " +polymerBr.getOutIDCount() +" ,expected 0 inIDs, found: " +polymerBr.getInIDCount());
		}
	}

	/**
	 * Return the locant associated with the first child of a word if there is only child
	 * else returns null
	 * @param currentWord
	 * @return locant or null
	 */
	private String getLocantForSubstituent(Element currentWord) {
		Elements children = currentWord.getChildElements();
		if (children.size()==1){
			Element child =children.get(0);
			if (child.getAttribute("locant")!=null){
				return child.getAttributeValue("locant");
			}
			else{
				return null;
			}
		}
		else{
			return null;
		}
	}

	/**
	 * Finds a suitable functional atom corresponding to the given locant
	 * @param locant
	 * @param mainGroupBR
	 * @return functionalAtomToUse
	 * @throws StructureBuildingException
	 */
	private Atom determineFunctionalAtomToUse(String locant, BuildResults mainGroupBR) throws StructureBuildingException {
		for (int i = 0; i < mainGroupBR.getFunctionalIDCount(); i++) {
			Atom possibleAtom = mainGroupBR.getFunctionalAtom(i);
			if (possibleAtom.hasLocant(locant)){
				mainGroupBR.removeFunctionalID(i);
				return possibleAtom;
			}
		}
		if (matchDigits.matcher(locant).matches()){
			//None of the functional atoms had an appropriate locant. Look for the case whether the locant refers to the backbone. e.g. 5-methyl 2-aminopentanedioate
			for (int i = 0; i < mainGroupBR.getFunctionalIDCount(); i++) {
				Atom possibleAtom = mainGroupBR.getFunctionalAtom(i);
				if (OpsinTools.depthFirstSearchForNonSuffixAtomWithLocant(possibleAtom, locant)!=null){
					mainGroupBR.removeFunctionalID(i);
					return possibleAtom;
				}
			}
		}
		else{
			//None of the functional atoms had an appropriate locant. Look for the special case where the locant is used to decide on the ester configuration c.f. O-methyl ..thioate and S-methyl ..thioate
			for (int i = 0; i < mainGroupBR.getFunctionalIDCount(); i++) {
				Atom possibleAtom = mainGroupBR.getFunctionalAtom(i);
				if (possibleAtom.getNote("ambiguousElementAssignment")!=null){
					String[] atomIDs =possibleAtom.getNote("ambiguousElementAssignment").split(",");
                    for (String atomID : atomIDs) {
                        Atom a = mainGroupBR.getAtomByIdOrThrow(Integer.parseInt(atomID));
                        if (a.hasLocant(locant)) {
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
                            mainGroupBR.removeFunctionalID(i);
                            return possibleAtom;
                        }
                    }
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
			List<Atom> atomList =fragment.getAtomList();
			for (Atom parentAtom : atomList) {
				Integer valency =parentAtom.determineValency(false);
				int explicitHydrogensToAdd=0;
				if (valency !=null){
					explicitHydrogensToAdd=valency-parentAtom.getIncomingValency();
					parentAtom.setExplicitHydrogens(explicitHydrogensToAdd);
				}
				int currentId = 0;
				for (int i = 1; i <= explicitHydrogensToAdd; i++) {
					Atom a = state.fragManager.createAtom("H", fragment);
					state.fragManager.createBond(parentAtom, a, 1);
				}
				if (parentAtom.getAtomParityElement()!=null){
					if (explicitHydrogensToAdd >1){
						throw new StructureBuildingException("Cannot have tetrahedral chirality and more than 2 hydrogens");
					}
					if (explicitHydrogensToAdd ==1){
						Element atomParityEl = parentAtom.getAtomParityElement();
						Attribute atomRefs4Atr = atomParityEl.getAttribute("atomRefs4");
						String atomRefs4 = atomRefs4Atr.getValue();
						atomRefs4 = atomRefs4.replaceFirst("a" + parentAtom.getID() +"_H", "a" +currentId);//atom parity was set in SMILES but at this stage the id of the hydrogen was not known, now it is so replace the dummy ID
						atomRefs4Atr.setValue(atomRefs4);
					}
				}
			}
		}
	}
}
