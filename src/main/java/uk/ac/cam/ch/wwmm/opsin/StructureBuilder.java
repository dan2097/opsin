package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import uk.ac.cam.ch.wwmm.ptclib.string.StringTools;
import uk.ac.cam.ch.wwmm.ptclib.xml.XOMTools;

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

	/**	Builds a molecule as a Fragment based on preStructurebuilder output.
	 * @param state 
	 * @param molecule The preStructurebuilderr output.
	 * @return A single Fragment - the built molecule.
	 * @throws StructureBuildingException If the molecule won't build - there may be many reasons.
	 */
	Fragment buildFragment(BuildState state, Element molecule) throws StructureBuildingException {
		String wordRule = state.wordRule;
		Elements words = molecule.getChildElements();
		if (words.size()==0){
			throw new StructureBuildingException("Molecule contains no words!?");
		}

		if(wordRule.equals("simple")) {
			Element word = molecule.getFirstChildElement("word");
			resolveWordOrBracket(state, word);
		}
		else if (wordRule.equals("acid")){
			buildAcid(state,words);//ethanoic acid
		}
		else if(wordRule.equals("ester") || wordRule.equals("diester")) {
			buildEster(state, words);//e.g. ethyl ethanoate, dimethyl terephthalate,  methyl propanamide
		}
		else if (wordRule.equals("divalentLiteralFunctionalGroup")){
			buildDiValentLiteralFunctionalGroup(state, words);// diethyl ether or methyl propyl ketone
		}
		else if (wordRule.equals("monovalentFunctionalGroup")){
			buildMonovalentFunctionalGroup(state, words);// ethyl chloride or isophthaloyl dichloride
		}
		else if (wordRule.equals("monovalentLiteralFunctionalGroup")){
			buildMonovalentLiteralFunctionalGroup(state, words);//e.g. diethyl ether, ethyl alcohol
		}
		else if(wordRule.equals("functionalClassEster")) {
			buildFunctionalClassEster(state, words);//e.g. ethanoic acid ethyl ester, tetrathioterephthalic acid dimethyl ester
		}
		else if (wordRule.equals("amide")){
			buildAmide(state, words);//e.g. ethanoic acid ethyl amide, terephthalic acid dimethyl amide, ethanoic acid amide
		}
		else if (wordRule.equals("glycol")){
			buildGlycol(state, words);//e.g. ethylene glycol
		}
		else if(wordRule.equals("binaryOrOther")) {
			for (int i = 0; i < words.size(); i++) {
				resolveWordOrBracket(state, words.get(i));
			}
		}
		else if(wordRule.equals("polymer")) {
			buildPolymer(state, words);
		}
		else{
			throw new StructureBuildingException("Unknown Word Rule");
		}
//		else if (wordRule.equals("oxime")){//e.g. Imidazole-2-carboxamide O-ethyloxime, pentan-3-one oxime
//			int substituentPresent;
//			if (words.size()==2){
//				substituentPresent=0;
//			}
//			else if (words.size()==3){
//				substituentPresent=1;
//			}
//			else{
//				throw new StructureBuildingException("Don't alter wordRules.xml without checking the consequences!");
//			}
//			
//			int numberOfOximes =1;
//			moleculeBuildResults =resolveWordOrBracket(words.get(0), null, new LinkedHashSet<Fragment>());//the group
//			List<List<Atom>> matches = moleculeBuildResults.subStructureSearch("O=C", state.fragManager);
//			System.out.println(matches.size());
//			if (matches.size() < numberOfOximes){
//				throw new StructureBuildingException("Insufficient carbonyl groups found!");
//			}
//
//			for (int i = 0; i < numberOfOximes; i++) {
//				List<Atom> atomList = matches.get(i);
//				Atom atomToBeReplaced =atomList.get(0);//the oxygen of the carbonyl
//				Fragment parentFrag =atomToBeReplaced.getFrag();
//				atomToBeReplaced.setElement("N");
//				int ID = state.idManager.getNextID();
//				List<Atom> parentFragAtomList =parentFrag.getAtomList();
//				for (Atom atom :parentFragAtomList) {
//					if (atom.hasLocant("O")){
//						atom.removeLocant("O");
//					}
//				}
//				Atom addedHydroxy = new Atom(ID, "O", "O", parentFrag);
//				parentFrag.addAtom(addedHydroxy);
//				Bond newBond =new Bond(atomToBeReplaced.getID(), addedHydroxy.getID(), 1);
//				parentFrag.addBond(newBond);
//				if (i== 0 && substituentPresent==1){
//					moleculeBuildResults.mergeBuildResults(resolveWordOrBracket(words.get(1), null, new LinkedHashSet<Fragment>()));
//					if (moleculeBuildResults.getOutIDCount() !=1){
//						throw new StructureBuildingException("Expected outID on substituent before oxime");
//					}
//					state.fragManager.attachFragments(addedHydroxy, moleculeBuildResults.getOutAtomTakingIntoAccountWhetherSetExplicitly(0), moleculeBuildResults.getOutID(0).valency);
//					moleculeBuildResults.removeOutID(0);
//				}
//			}
//		}

		state.fragManager.tidyUpFragments();//FIXME remove this?
		state.fragManager.convertSpareValenciesToDoubleBonds();
		state.fragManager.checkValencies();
		makeHydrogensExplicit(state);
		matchStereochemistryToAtomsAndBonds(state, molecule);
		Fragment uniFrag = state.fragManager.getUnifiedFragment();
		applyStereochemistry(uniFrag);
		
		//adds Xe group at all atoms which have unused outIDs
		//note that SMILES generated by CDK is not always correct
//		for (int i = 0; i < moleculeBuildResults.getOutIDCount(); i++) {
//			Atom outAtom =moleculeBuildResults.getOutAtomTakingIntoAccountWhetherSetExplicitly(i);
//			Fragment rGroup =state.fragManager.buildSMILES("[Xe]");
//			state.fragManager.attachFragments(outAtom, rGroup.getAtomByIDOrThrow(rGroup.getIdOfFirstAtom()), moleculeBuildResults.getOutID(i).valency);
//		}

		if (uniFrag.getOutIDs().size()>0 || uniFrag.getInIDs().size()>0){
			//throw new StructureBuildingException("Radicals are currently set to not convert to structures");
		}
		return uniFrag;
	}

	private void buildAcid(BuildState state, Elements words) throws StructureBuildingException {
		resolveWordOrBracket(state, words.get(0));
		if (words.size()<2 || !words.get(1).getAttributeValue("type").equals("literal")){
			throw new StructureBuildingException("literal word acid missing");
		}
		resolveTrailingFullWords(state, words, 2);
	}

	private void buildEster(BuildState state, Elements words) throws StructureBuildingException {
		int wordIndice = resolvePreceedingFullWordsAndReturnNonFullIndice(state, words);
		Element currentWord=words.get(wordIndice);
		BuildResults substituentsBr = new BuildResults();
		while (currentWord.getAttributeValue("type").equals("substituent")){
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
		
		if (wordIndice == words.size() || !words.get(wordIndice).getAttributeValue("type").equals("full")){
			throw new StructureBuildingException("Full word not found where full word expected: missing ate group in ester");
		}
		
		BuildResults ateGroups = new BuildResults();
		for (; wordIndice < words.size(); wordIndice++) {
			Element word =words.get(wordIndice);
			if (word.getAttributeValue("type").equals("full")){
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
			state.fragManager.attachFragments(ateAtom,substituentsBr.getOutAtomTakingIntoAccountWhetherSetExplicitly(0), 1);
			substituentsBr.removeOutID(0);
			ateAtom.setCharge(0);
		}
	}



	private void buildDiValentLiteralFunctionalGroup(BuildState state, Elements words) throws StructureBuildingException {
		int wordIndice = resolvePreceedingFullWordsAndReturnNonFullIndice(state, words);
		if (!words.get(wordIndice).getAttributeValue("type").equals("substituent")) {
			throw new StructureBuildingException("word: " +wordIndice +" was expected to be a substituent");
		}
		if (words.get(wordIndice +1).getAttributeValue("type").equals("literal")) {//e.g. methyl sulfoxide rather than dimethyl sulfoxide
			Element clone = state.fragManager.cloneElement(state, words.get(0));
			XOMTools.insertAfter(words.get(0), clone);
			words = ((Element)words.get(0).getParent()).getChildElements();
		}
		resolveWordOrBracket(state, words.get(wordIndice));
		BuildResults substituent1 =new BuildResults(state, words.get(wordIndice));
		if (substituent1.getOutIDCount()!=1){
			throw new StructureBuildingException("Expected one out outID. Found " + substituent1.getOutIDCount() );
		}
		if (substituent1.getOutID(0).valency !=1){
			throw new StructureBuildingException("OutID has unexpected valency. Expected 1. Actual: " + substituent1.getOutID(0).valency);
		}
		wordIndice++;
		resolveWordOrBracket(state, words.get(wordIndice));
		BuildResults substituent2 =new BuildResults(state, words.get(wordIndice));
		if (substituent2.getOutIDCount()!=1){
			throw new StructureBuildingException("Expected one out outID. Found " + substituent2.getOutIDCount() );
		}
		if (substituent2.getOutID(0).valency !=1){
			throw new StructureBuildingException("OutID has unexpected valency. Expected 1. Actual: " + substituent2.getOutID(0).valency);
		}
		wordIndice++;
		if (words.get(wordIndice) ==null || !words.get(wordIndice).getAttributeValue("type").equals("literal")) {
			throw new StructureBuildingException("word: " +wordIndice +" was expected to be a literal");
		}
		String smilesOfGroup = null;
		String functionalGroupName =words.get(wordIndice).getValue();
		if (functionalGroupName.equalsIgnoreCase("ether")){
			smilesOfGroup="O";
		}
		else if (functionalGroupName.equalsIgnoreCase("ketone")){
			smilesOfGroup="C=O";
		}
		else if (functionalGroupName.equalsIgnoreCase("peroxide")){
			smilesOfGroup="OO-";
		}
		else if (functionalGroupName.equalsIgnoreCase("selenide")){
			smilesOfGroup="[Se]";
		}
		else if (functionalGroupName.equalsIgnoreCase("selenone")){
			smilesOfGroup="[Se](=O)=O";
		}
		else if (functionalGroupName.equalsIgnoreCase("selenoxide")){
			smilesOfGroup="[Se]=O";
		}
		else if (functionalGroupName.equalsIgnoreCase("selone")){
			smilesOfGroup="C=[Se]";
		}
		else if (functionalGroupName.equalsIgnoreCase("selenoketone")){
			smilesOfGroup="[Se]=O";
		}
		else if (functionalGroupName.equalsIgnoreCase("sulfide")){
			smilesOfGroup="S";
		}
		else if (functionalGroupName.equalsIgnoreCase("sulfone")){
			smilesOfGroup="S(=O)=O";
		}
		else if (functionalGroupName.equalsIgnoreCase("sulfoxide")){
			smilesOfGroup="S=O";
		}
		else if (functionalGroupName.equalsIgnoreCase("telluride")){
			smilesOfGroup="[Te]";
		}
		else if (functionalGroupName.equalsIgnoreCase("telluroketone")){
			smilesOfGroup="[Te]=O";
		}
		else if (functionalGroupName.equalsIgnoreCase("tellurone")){
			smilesOfGroup="[Te](=O)=O";
		}
		else if (functionalGroupName.equalsIgnoreCase("telluroxide")){
			smilesOfGroup="[Te]=O";
		}
		else if (functionalGroupName.equalsIgnoreCase("thioketone")){
			smilesOfGroup="S=O";
		}
		else{
			throw new StructureBuildingException("Unknown functionalGroup: " + functionalGroupName);
		}
		Fragment diValentGroup =state.fragManager.buildSMILES(smilesOfGroup, "simpleGroup", "functionalGroup", "none");
		
		Atom outAtom =substituent1.getOutAtomTakingIntoAccountWhetherSetExplicitly(0);
		substituent1.removeOutID(0);
		if (diValentGroup.getOutIDs().size()==1){//c.f. peroxide where it is a linker
			state.fragManager.attachFragments(outAtom, diValentGroup.getAtomByIDOrThrow(diValentGroup.getOutID(0).id), 1);
			diValentGroup.removeOutID(0);
		}
		else{
			state.fragManager.attachFragments(outAtom, diValentGroup.getAtomByIDOrThrow(diValentGroup.getIdOfFirstAtom()), 1);
		}
		outAtom = substituent2.getOutAtomTakingIntoAccountWhetherSetExplicitly(0);
		substituent2.removeOutID(0);
		state.fragManager.attachFragments(outAtom, diValentGroup.getAtomByIDOrThrow(diValentGroup.getIdOfFirstAtom()), 1);
		resolveTrailingFullWords(state, words, wordIndice + 1);
	}

	private void buildMonovalentFunctionalGroup(BuildState state, Elements words) throws StructureBuildingException {
		int wordIndice = resolvePreceedingFullWordsAndReturnNonFullIndice(state, words);
		resolveWordOrBracket(state, words.get(wordIndice));
		List<Element> groups = OpsinTools.findDescendantElementsWithTagName(words.get(0), "group");
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
	
		int numberOfOutIDs =substituentBR.getOutIDCount();
		if (numberOfOutIDs > words.size()-1){//something like isophthaloyl chloride (more precisely written isophthaloyl dichloride)
			if (words.size()-wordIndice != 2){
				throw new StructureBuildingException("Incorrect number of functional class groups found to balance outIDs");
			}
			int diff =numberOfOutIDs - (words.size()-wordIndice -1);
			Element elementToBeCloned =words.get(words.size()-1);
			for (int i = 0; i < diff; i++) {
				Element clone =state.fragManager.cloneElement(state, elementToBeCloned);
				XOMTools.insertAfter(elementToBeCloned, clone);
				elementToBeCloned=clone;
			}
			words = ((Element)words.get(0).getParent()).getChildElements();//update words
		}
		for (int i = 0; i < numberOfOutIDs; i++) {
			wordIndice++;
			resolveWordOrBracket(state, words.get(wordIndice));
			Element ideGroup = StructureBuildingMethods.findRightMostGroupInBracket(words.get(wordIndice));
			Fragment ideFrag = state.xmlFragmentMap.get(ideGroup);
			Atom ideAtom = ideFrag.getAtomByIDOrThrow(ideFrag.getDefaultInID());
			Atom subAtom=substituentBR.getOutAtomTakingIntoAccountWhetherSetExplicitly(0);
			state.fragManager.attachFragments(ideAtom, subAtom, 1);
			substituentBR.removeOutID(0);
			ideAtom.setCharge(0);
		}
		resolveTrailingFullWords(state, words, ++wordIndice);//e.g. ethyl chloride hydrochloride
	}

	private void buildMonovalentLiteralFunctionalGroup(BuildState state,Elements words) throws StructureBuildingException {
		int wordIndice = resolvePreceedingFullWordsAndReturnNonFullIndice(state, words);
		resolveWordOrBracket(state, words.get(wordIndice));//the substituent
		BuildResults substituentBr = new BuildResults(state, words.get(wordIndice));
		wordIndice++;
		String smilesOfGroup = null;
		String functionalGroupName =words.get(wordIndice++).getValue();
		if (functionalGroupName.equalsIgnoreCase("alcohol")){
			smilesOfGroup="O";
		}
		else if (functionalGroupName.equalsIgnoreCase("selenol")){
			smilesOfGroup="[Se]";
		}
		else if (functionalGroupName.equalsIgnoreCase("thiol")){
			smilesOfGroup="S";
		}
		else{
			throw new StructureBuildingException("Unknown functionalGroup: " + functionalGroupName);
		}
		for (int i = 0; i < substituentBr.getOutIDCount(); i++) {
			Atom outAtom = substituentBr.getOutAtomTakingIntoAccountWhetherSetExplicitly(i);
			Fragment rGroup =state.fragManager.buildSMILES(smilesOfGroup, "monovalentgroup", "none");
			if (substituentBr.getOutID(i).valency !=1){
				throw new StructureBuildingException("OutID has unexpected valency. Expected 1. Actual: " + substituentBr.getOutID(i).valency);
			}
			state.fragManager.attachFragments(outAtom, rGroup.getAtomByIDOrThrow(rGroup.getIdOfFirstAtom()), 1);
			substituentBr.removeOutID(i);
		}
		resolveTrailingFullWords(state, words, wordIndice);
	}

	private void buildFunctionalClassEster(BuildState state, Elements words) throws StructureBuildingException {
		if (!words.get(0).getAttributeValue("type").equals("full")){
			throw new StructureBuildingException("Don't alter wordRules.xml without checking the consequences!");
		}
		resolveWordOrBracket(state, words.get(0));//the group
		BuildResults acidBr = new BuildResults(state, words.get(0));
		if (!words.get(1).getAttributeValue("type").equals("literal")){//acid
			throw new StructureBuildingException("Don't alter wordRules.xml without checking the consequences!");
		}
		if (acidBr.getFunctionalIDCount()==0){
			throw new StructureBuildingException("No functionalIds detected!");
		}
	
		int i=2;
		Element currentWord = words.get(i);
		while (currentWord.getAttributeValue("type").equals("substituent")){
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
				state.fragManager.attachFragments(functionalAtom,substituentBr.getOutAtomTakingIntoAccountWhetherSetExplicitly(0), 1);
				substituentBr.removeOutID(0);
			}
			else {
				throw new StructureBuildingException("Substituent was expected to have one outID");
			}
			currentWord=words.get(++i);
		}
		if (!words.get(i++).getAttributeValue("type").equals("literal")){//ester
			throw new StructureBuildingException("Number of words different to expectations; did not find ester");
		}
		resolveTrailingFullWords(state, words, i);
	}

	private void buildAmide(BuildState state, Elements words) throws StructureBuildingException {
		if (!words.get(0).getAttributeValue("type").equals("full")){
			throw new StructureBuildingException("Don't alter wordRules.xml without checking the consequences!");
		}
		resolveWordOrBracket(state, words.get(0));//the group
		BuildResults acidBr = new BuildResults(state, words.get(0));
		if (acidBr.getFunctionalIDCount()==0){
			throw new StructureBuildingException("No functionalIds detected!");
		}
		if (!words.get(1).getAttributeValue("type").equals("literal")){//acid
			throw new StructureBuildingException("Don't alter wordRules.xml without checking the consequences!");
		}
		int substituentCount = OpsinTools.findChildElementsWithTagNameAndAttribute((Element) words.get(0).getParent(), "word", "type", "substituent").size();
		if (acidBr.getFunctionalIDCount() < substituentCount){
			throw new StructureBuildingException("More substituents than acid functionalIDs detcted during amide construction!");
		}
		if (substituentCount==0){//e.g. ethanoic acid amide
			for (int i =0; i < acidBr.getFunctionalIDCount(); i++) {
				Atom functionalAtom = acidBr.getFunctionalAtom(i);
				if (!functionalAtom.getElement().equals("O")){
					throw new StructureBuildingException("Expected oxygen functional atom found:" + functionalAtom.getElement());
				}
				functionalAtom.setElement("N");
				functionalAtom.replaceLocant("N" +StringTools.multiplyString("'", i));
			}
		}
		int wordIndice;
		for (wordIndice = 2 ; wordIndice < 2 + substituentCount; wordIndice++) {
			if (!words.get(wordIndice).getAttributeValue("type").equals("substituent")){
				throw new StructureBuildingException("Word: " + wordIndice + " was expected to be a substituent");
			}
			//a fragment is constructed which is just the Nitrogen atom of the amide
			Element dummyRoot = new Element("root");
			Element dummGroup = new Element("group");
			dummyRoot.appendChild(dummGroup);
			Fragment dummyFrag = state.fragManager.buildSMILES("N", "suffix", "N");
			state.xmlFragmentMap.put(dummGroup, dummyFrag);
			words.get(wordIndice).appendChild(dummyRoot);
			resolveWordOrBracket(state, words.get(wordIndice));
			BuildResults substituentBr = new BuildResults(state, words.get(wordIndice));
			if (substituentBr.getOutIDCount()!=0){
				throw new StructureBuildingException("Substituent should have not OutIDs after forming amide");
			}
			
			Atom functionalAtom = acidBr.getFunctionalAtom(0);
			if (!functionalAtom.getElement().equals("O")){
				throw new StructureBuildingException("Expected oxygen functional atom found:" + functionalAtom.getElement());
			}
			acidBr.removeFunctionalID(0);
			state.fragManager.replaceTerminalAtomWithFragment(functionalAtom, dummyFrag.getFirstAtom());//the first atom will be the amide N
		}
		if (!words.get(wordIndice++).getAttributeValue("type").equals("literal")){//amide
			throw new StructureBuildingException("Number of words different to expectations; did not find amide");
		}
		resolveTrailingFullWords(state, words, wordIndice);
	}

	private void buildGlycol(BuildState state, Elements words) throws StructureBuildingException {
		int wordIndice  = resolvePreceedingFullWordsAndReturnNonFullIndice(state, words);
		resolveWordOrBracket(state, words.get(wordIndice));//the group
		BuildResults theDiRadical = new BuildResults(state, words.get(wordIndice));
		if (theDiRadical.getOutIDCount()!=2){
			throw new StructureBuildingException("Glycol class names (e.g. ethylene glycol) expect two outIDs. Found: " + theDiRadical.getOutIDCount() );
		}
		if (wordIndice +1 >= words.size() || words.get(wordIndice+1).getAttributeValue("type").equals("glycol")){
			throw new StructureBuildingException("Glycol literal word expected");
		}
		for (int i = 0; i < theDiRadical.getOutIDCount(); i++) {
			Atom outAtom =theDiRadical.getOutAtomTakingIntoAccountWhetherSetExplicitly(i);
			Fragment glycol =state.fragManager.buildSMILES("O", "glycol", "none");
			if (theDiRadical.getOutID(i).valency !=1){
				throw new StructureBuildingException("OutID has unexpected valency. Expected 1. Actual: " + theDiRadical.getOutID(i).valency);
			}
			theDiRadical.removeOutID(0);
			state.fragManager.attachFragments(outAtom, glycol.getAtomByIDOrThrow(glycol.getIdOfFirstAtom()), 1);
		}
		resolveTrailingFullWords(state, words, wordIndice + 2);
	}

	private void buildPolymer(BuildState state, Elements words) throws StructureBuildingException {
		if (words.size()>1){
			throw new StructureBuildingException("Currently unsupported polymer name type");
		}
		Element polymer = words.get(0);
		resolveWordOrBracket(state, polymer);
		BuildResults polymerBr = new BuildResults(state, polymer);
		if (polymerBr.getOutIDCount() ==2 && polymerBr.getInIDCount()==0){
			Atom inAtom =polymerBr.getOutAtomTakingIntoAccountWhetherSetExplicitly(0);
			Fragment rGroup =state.fragManager.buildSMILES("[Xe]");
			state.fragManager.attachFragments(inAtom, rGroup.getAtomByIDOrThrow(rGroup.getIdOfFirstAtom()), polymerBr.getOutID(0).valency);
			Atom outAtom =polymerBr.getOutAtomTakingIntoAccountWhetherSetExplicitly(1);
			rGroup =state.fragManager.buildSMILES("[Rn]");
			state.fragManager.attachFragments(outAtom, rGroup.getAtomByIDOrThrow(rGroup.getIdOfFirstAtom()), polymerBr.getOutID(1).valency);
			polymerBr.removeAllOutIDs();
		}
		else{
			throw new StructureBuildingException("Polymer building failed: Two termini were not found; Expected 2 outIDs, found: " +polymerBr.getOutIDCount() +" ,expected 0 inIDs, found: " +polymerBr.getInIDCount());
		}
	}

	private int resolvePreceedingFullWordsAndReturnNonFullIndice(BuildState state, Elements words) throws StructureBuildingException {
		for (int i = 0; i < words.size(); i++) {
			Element word =words.get(i);
			if (word.getAttributeValue("type").equals("full")){
				resolveWordOrBracket(state, word);
			}
			else{
				return i;
			}
		}
		throw new StructureBuildingException("No non full words were encountered!");
	}
	
	private void resolveTrailingFullWords(BuildState state, Elements words, int indice) throws StructureBuildingException {
		for (int i = indice; i < words.size(); i++) {
			Element word =words.get(i);
			if (word.getAttributeValue("type").equals("full")){
				resolveWordOrBracket(state, word);
			}
			else{
				throw new StructureBuildingException("Non full word found where only full words were expected");
			}
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
				String locant =child.getAttributeValue("locant");
				if (!locant.matches("\\d+")){//FIXME accept numeric locants
					return locant;
				}
				else{
					return null;
				}
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
		//None of the functional atoms had an appropriate locant. Look for the special case where the locant is used to decide on the ester configuration c.f. O-methyl ..thioate and S-methyl ..thioate
		for (int i = 0; i < mainGroupBR.getFunctionalIDCount(); i++) {
			Atom possibleAtom = mainGroupBR.getFunctionalAtom(i);
			if (possibleAtom.getNote("ambiguousElementAssignment")!=null){
				String[] atomIDs =possibleAtom.getNote("ambiguousElementAssignment").split(",");
				for (int j = 0; j < atomIDs.length; j++) {
					Atom a =mainGroupBR.getAtomByIdOrThrow(Integer.parseInt(atomIDs[j]));
					if (a.hasLocant(locant)){
						//swap locants and element type
						List<String> tempLocants =new ArrayList<String>(a.getLocants());
						List<String> tempLocants2 =new ArrayList<String>(possibleAtom.getLocants());
						a.clearLocants();
						possibleAtom.clearLocants();
						for (String l : tempLocants) {
							possibleAtom.addLocant(l);
						}
						for (String l : tempLocants2) {
							a.addLocant(l);
						}
						String originalElement =possibleAtom.getElement();
						possibleAtom.setElement(a.getElement());
						a.setElement(originalElement);
						return possibleAtom;
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
				Integer valency =parentAtom.determineValency();
				int explicitHydrogensToAdd=0;
				if (valency !=null){
					explicitHydrogensToAdd=valency-parentAtom.getIncomingValency();
					parentAtom.setExplicitHydrogens(explicitHydrogensToAdd);
				}
				int currentId = 0;
				for (int i = 1; i <= explicitHydrogensToAdd; i++) {
					currentId = state.idManager.getNextID();
					Atom a = new Atom(currentId, "H", fragment);
					fragment.addAtom(a);
					Bond b =new Bond(parentAtom, a, 1);
					fragment.addBond(b);
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

	/**
	 * Finds all stereochemistry elements and attempts to locate a suitable atom/bond for them
	 * Currently only locanted stereochemistry elements are supported
	 * @param state 
	 * @param elem
	 * @throws StructureBuildingException 
	 */
	private void matchStereochemistryToAtomsAndBonds(BuildState state, Element molecule) throws StructureBuildingException {
		List<Element> stereoChemistryEls = OpsinTools.findDescendantElementsWithTagName(molecule, "stereoChemistry");
		for (Element stereoChemistryEl : stereoChemistryEls) {
			String stereoChemistryType =stereoChemistryEl.getAttributeValue("type");
			Element parent = (Element) stereoChemistryEl.getParent().getParent();//want to iterate at the level above the containing substituent or bracket
			//generally the LAST group in this list will be the appropriate groups e.g. (5S)-5-ethyl-6-methylheptane where the heptane is the appropriate group
			List<Element> possibleGroups = OpsinTools.findDescendantElementsWithTagName(parent, "group");
			if (stereoChemistryType.equals("RorS")){
				String locant = getLocant(stereoChemistryEl);
				String rOrS = stereoChemistryEl.getAttributeValue("value");
				Atom a = null;
				for (int i = possibleGroups.size()-1; i >=0; i--) {
					a = state.xmlFragmentMap.get(possibleGroups.get(i)).getAtomByLocant(locant);
					if (a !=null && a.getBonds().size()==4 && (a.getExplicitHydrogens()==null || a.getExplicitHydrogens()<=1)){
						a.setStereochemistry(rOrS);//TODO check that atom is in fact chiral properly!
						break;
					}
				}
				if (a ==null && parent.getLocalName().equals("word") && parent.getAttributeValue("type").equals("substituent")){
					//something like (3R,4R,5R)-ethyl 4-acetamido-5-amino-3-(pentan-3-yloxy)cyclohex-1-enecarboxylate
					//I think this is a violation of the IUPAC rules...but anyway...
					List<Element> words = OpsinTools.findChildElementsWithTagNameAndAttribute(((Element)parent.getParent()), "word", "type", "full");
					wordLoop: for (Element word : words) {
						possibleGroups = OpsinTools.findDescendantElementsWithTagName(word, "group");
						for (int i = possibleGroups.size()-1; i >=0; i--) {
							a = state.xmlFragmentMap.get(possibleGroups.get(i)).getAtomByLocant(locant);
							if (a !=null && a.getBonds().size()==4 && (a.getExplicitHydrogens()==null || a.getExplicitHydrogens()<=1)){
								a.setStereochemistry(rOrS);//TODO check that atom is in fact chiral properly!
								break wordLoop;
							}
						}
					}
				}
				if (a==null){
					throw new StructureBuildingException("Could not find atom that: " + stereoChemistryEl.getValue() + " appeared to be referring to");
				}
			}
			else if (stereoChemistryType.equals("EorZ")){
				String locant = getLocant(stereoChemistryEl);
				String eOrZ = stereoChemistryEl.getAttributeValue("value");
				Atom a = null;
				Bond b = null;
				for (int i = possibleGroups.size()-1; i >=0; i--) {
					a = state.xmlFragmentMap.get(possibleGroups.get(i)).getAtomByLocant(locant);
					if (a !=null){
						List<Bond> bonds = a.getBonds();
						if (!checkAtomIsAppropriateForEZStereochemistry(a)){
							continue;
						}
						for (Bond potentialBond : bonds) {
							if (potentialBond.getOrder()==2){
								if (b == null){
									b = potentialBond;
								}
								else{//atom seems to connect to multiple double bonds!!!
									b = null;
								}
							}
						}
						if (b != null){
							Atom atomAtOtherEndOfBond;
							if (b.getFromAtom()==a){
								atomAtOtherEndOfBond = b.getToAtom();
							}
							else{
								atomAtOtherEndOfBond = b.getFromAtom();
							}
							if (!checkAtomIsAppropriateForEZStereochemistry(atomAtOtherEndOfBond)){//|| (a.getAtomIsInACycle()&& atomAtOtherEndOfBond.getAtomIsInACycle()
								b = null;
								continue;
							}
							b.setStereochemistry(eOrZ);
							break;
						}
					}
				}
				if ((a==null || b==null) && parent.getLocalName().equals("word") && parent.getAttributeValue("type").equals("substituent") ){
					//the element is in front of a substituent and may refer to the full group
					List<Element> words = OpsinTools.findChildElementsWithTagNameAndAttribute(((Element)parent.getParent()), "word", "type", "full");
					wordLoop: for (Element word : words) {
						possibleGroups = OpsinTools.findDescendantElementsWithTagName(word, "group");
						for (int i = possibleGroups.size()-1; i >=0; i--) {
							a = state.xmlFragmentMap.get(possibleGroups.get(i)).getAtomByLocant(locant);
							if (a !=null){
								List<Bond> bonds = a.getBonds();
								if (!checkAtomIsAppropriateForEZStereochemistry(a)){
									continue;
								}
								for (Bond potentialBond : bonds) {
									if (potentialBond.getOrder()==2){
										if (b == null){
											b = potentialBond;
										}
										else{//atom seems to connect to multiple double bonds!!!
											b = null;
										}
									}
								}
								if (b != null){
									Atom atomAtOtherEndOfBond;
									if (b.getFromAtom()==a){
										atomAtOtherEndOfBond = b.getToAtom();
									}
									else{
										atomAtOtherEndOfBond = b.getFromAtom();
									}
									if (!checkAtomIsAppropriateForEZStereochemistry(atomAtOtherEndOfBond)){
										b = null;
										continue;
									}
									b.setStereochemistry(eOrZ);
									break wordLoop;
								}
							}
						}
					}
				}
				if (a ==null){
					throw new StructureBuildingException("Could not find bond that: " + stereoChemistryEl.getValue() + " appeared to be referring to");
				}
				else if (b ==null){
					throw new StructureBuildingException("Bond that: " + stereoChemistryEl.getValue() + " appeared to be referring to was not a double bond or was not capable of having E/Z stereochemistry");
				}
			}
			else{
				//currently unsupported
			}
			stereoChemistryEl.detach();
		}
		
	}

	/**
	 * Confirms that EZ can be applied sensibly to this atom
	 * @param a
	 * @return
	 */
	private boolean checkAtomIsAppropriateForEZStereochemistry(Atom a) {
		List<Bond> bonds =a.getBonds();
		if (bonds.size()==3 && (a.getExplicitHydrogens()==null || a.getExplicitHydrogens()<=1)){
			return true;
		}
		if (bonds.size()==2 && (a.getElement().equals("N") || a.getElement().equals("P") || a.getElement().equals("As") || a.getElement().equals("Sb"))){
			//nitrogen has a lone pair!
			return true;
		}
		return false;
	}

	/**
	 * Searches for atoms that have been designating as having stereochemistry and applys the CIP rules
	 * to assign atom parity to the given atom
	 * @param uniFrag
	 * @throws StructureBuildingException 
	 */
	private void applyStereochemistry(Fragment uniFrag) throws StructureBuildingException {
		List<Atom> atomList =uniFrag.getAtomList();
		if (atomList.size()>300){return;}//FIXME temporary workaround for CIP code running out of memory
		for (Atom atom : atomList) {
			if (atom.getStereochemistry()!=null){
				List<Atom> cipOrderedAtoms = CIPresolver.getCIPOrderedAtoms(uniFrag, atom, atom.getAtomNeighbours());
				if (cipOrderedAtoms.size()!=4){
					throw new StructureBuildingException("Only tetrahedral chirality is currently supported");
				}
				String atomRefs4= "a"+cipOrderedAtoms.get(cipOrderedAtoms.size()-1).getID();
				for (int i = 0; i < cipOrderedAtoms.size() -1; i++) {//from highest to lowest (true for S) hence atomParity 1 for S
					atomRefs4 +=" a"+cipOrderedAtoms.get(i).getID();
				}
				if (atom.getStereochemistry().equals("R")){
					atom.setAtomParityElement(atomRefs4, -1);
				}
				else if (atom.getStereochemistry().equals("S")){
					atom.setAtomParityElement(atomRefs4, 1);
				}
				else{
					throw new StructureBuildingException("Unexpected stereochemistry type: " + atom.getStereochemistry());
				}
			}
		}
		List<Bond> bondList =uniFrag.getBondList();
		for (Bond bond : bondList) {
			if (bond.getStereochemistry()!=null){
				List<Atom> atomsAtFromEnd = new ArrayList<Atom>();
				atomsAtFromEnd.addAll(bond.getFromAtom().getAtomNeighbours());
				atomsAtFromEnd.remove(bond.getToAtom());
				List<Atom> orderedFromEndAtoms = CIPresolver.getCIPOrderedAtoms(uniFrag, bond.getFromAtom(), atomsAtFromEnd);
				List<Atom> atomsAtToEnd = new ArrayList<Atom>();
				atomsAtToEnd.addAll(bond.getToAtom().getAtomNeighbours());
				atomsAtToEnd.remove(bond.getFromAtom());
				List<Atom> orderedToEndAtoms = CIPresolver.getCIPOrderedAtoms(uniFrag, bond.getToAtom(), atomsAtToEnd);
				String atomRefs4= "a" + orderedFromEndAtoms.get(0).getID();//orderedNodes1.get(1) will be the higher priority group
				atomRefs4 +=" a" + bond.getFromAtom().getID();
				atomRefs4 +=" a" + bond.getToAtom().getID();
				atomRefs4 +=" a" + orderedToEndAtoms.get(0).getID();
				if (bond.getStereochemistry().equals("E")){
					bond.setBondStereoElement(atomRefs4, "T");
				}
				else if (bond.getStereochemistry().equals("Z")){
					bond.setBondStereoElement(atomRefs4, "C");
				}
				else{
					throw new StructureBuildingException("Unexpected stereochemistry type: " + bond.getStereochemistry());
				}
			}
		}
	}
}
