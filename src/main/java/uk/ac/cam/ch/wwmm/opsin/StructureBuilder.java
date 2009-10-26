package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.Stack;

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
		if(wordRule.equals("ester") || wordRule.equals("diester")) {
			buildEster(state, words);//e.g. ethyl ethanoate, dimethyl terephthalate,  methyl propanamide
		}
		else if(wordRule.equals("polymer")) {
			buildPolymer(state, words);
		}else if(wordRule.equals("simple") || wordRule.equals("acid")) {
			Element word = molecule.getFirstChildElement("word");
			resolveWordOrBracket(state, word);
		}
		else{
			throw new StructureBuildingException("Unknown Word Rule");
		}
		/*else if (wordRule.equals("functionalClassEster")){//e.g. ethanoic acid ethyl ester, tetrathioterephthalic acid dimethyl ester
			if (!words.get(0).getAttributeValue("type").equals("full")){
				throw new StructureBuildingException("Don't alter wordRules.xml without checking the consequences!");
			}
			moleculeBuildResults =resolveWordOrBracket(state, words.get(0), null, new LinkedHashSet<Fragment>());//the group
			if (!words.get(1).getAttributeValue("type").equals("literal")){//acid
				throw new StructureBuildingException("Don't alter wordRules.xml without checking the consequences!");
			}
			if (moleculeBuildResults.functionalIDs.size()==0){
				throw new StructureBuildingException("No functionalIds detected!");
			}
			int i;
			int subsToCheck = Math.min(moleculeBuildResults.functionalIDs.size(),substituentCount);
			for (i = 2 ; i < 2 + subsToCheck; i++) {
				if (!words.get(i).getAttributeValue("type").equals("substituent")){
					throw new StructureBuildingException("Word: " + i + " was expected to be a substituent to satisfy the number of functionalIds");
				}
				moleculeBuildResults.mergeBuildResults(resolveWordOrBracket(state, words.get(i), null, new LinkedHashSet<Fragment>()));
				
				String locantForSubstituent = getLocantForSubstituent(words.get(i));
				Atom functionalAtom;
				if (locantForSubstituent!=null){
					functionalAtom =determineFunctionalAtomToUse(locantForSubstituent, moleculeBuildResults);
				}
				else{
					functionalAtom =moleculeBuildResults.getFunctionalOutAtom(0);
					moleculeBuildResults.functionalIDs.remove(0);
				}
				state.fragManager.attachFragments(functionalAtom,moleculeBuildResults.getOutAtomTakingIntoAccountWhetherSetExplicitly(0), 1);
				moleculeBuildResults.removeOutID(0);
			}
			if (!words.get(i).getAttributeValue("type").equals("literal")){//ester
				throw new StructureBuildingException("Number of words different to expectations; did not find ester");
			}
			if (moleculeBuildResults.getOutIDCount()!=0){
				throw new StructureBuildingException("Could not find corresponding functionalId for all outIds!");
			}
		}
		else if (wordRule.equals("amide")){//e.g. ethanoic acid ethyl amide, terephthalic acid dimethyl amide, ethanoic acid amide
			if (!words.get(0).getAttributeValue("type").equals("full")){
				throw new StructureBuildingException("Don't alter wordRules.xml without checking the consequences!");
			}
			moleculeBuildResults =resolveWordOrBracket(state, words.get(0), null, new LinkedHashSet<Fragment>());//the group
			if (!words.get(1).getAttributeValue("type").equals("literal")){//acid
				throw new StructureBuildingException("Don't alter wordRules.xml without checking the consequences!");
			}
			if (moleculeBuildResults.functionalIDs.size()==0){
				throw new StructureBuildingException("No functionalIds detected!");
			}
			int subsToCheck = Math.min(moleculeBuildResults.functionalIDs.size(),substituentCount);
			if (subsToCheck==0){//e.g. ethanoic acid amide
				for (int i =0; i < moleculeBuildResults.functionalIDs.size(); i++) {
					Atom functionalAtom = moleculeBuildResults.getFunctionalOutAtom(i);
					if (!functionalAtom.getElement().equals("O")){
						throw new StructureBuildingException("Expected oxygen functional atom found:" + functionalAtom.getElement());
					}
					functionalAtom.setElement("N");
					functionalAtom.replaceLocant("N" +StringTools.multiplyString("'", i));
				}
			}
			int wordIndice;
			for (wordIndice = 2 ; wordIndice < 2 + subsToCheck; wordIndice++) {
				if (!words.get(wordIndice).getAttributeValue("type").equals("substituent")){
					throw new StructureBuildingException("Word: " + wordIndice + " was expected to be a substituent to satisfy the number of functionalIds");
				}
				//a fragment is constructed which is just the Nitrogen atom of the amide
				Element dummyRoot = new Element("root");
				Element dummGroup = new Element("group");
				dummyRoot.appendChild(dummGroup);
				Fragment dummyFrag = state.fragManager.buildSMILES("N", "suffix", "N");
				state.xmlFragmentMap.put(dummGroup, dummyFrag);
				words.get(wordIndice).appendChild(dummyRoot);
				moleculeBuildResults.mergeBuildResults(resolveWordOrBracket(state, words.get(wordIndice), null, new LinkedHashSet<Fragment>()));
				
				Atom functionalAtom =moleculeBuildResults.getFunctionalOutAtom(0);
				if (!functionalAtom.getElement().equals("O")){
					throw new StructureBuildingException("Expected oxygen functional atom found:" + functionalAtom.getElement());
				}
				moleculeBuildResults.functionalIDs.remove(0);
				state.fragManager.replaceTerminalAtomWithFragment(functionalAtom, dummyFrag.getFirstAtom());//the first atom will be the amide N
			}
			if (!words.get(wordIndice).getAttributeValue("type").equals("literal")){//amide
				throw new StructureBuildingException("Number of words different to expectations; did not find amide");
			}
			if (moleculeBuildResults.getOutIDCount()!=0){
				throw new StructureBuildingException("Could not find corresponding functionalId for all outIds!");
			}
		}
		else if (wordRule.equals("monovalentFunctionalGroup")){// ethyl chloride or isophthaloyl dichloride
			BuildResults substituentBR =resolveWordOrBracket(state, words.get(0), null, new LinkedHashSet<Fragment>());
			moleculeBuildResults =substituentBR;
			for (int i = 0; i < substituentBR.getOutIDCount(); i++) {//replaces outIDs with valency greater than 1 with multiple outIDs; e.g. ylidene -->diyl
				OutID outID =substituentBR.getOutID(i);
				if (outID.valency>1){
					for (int j = 2; j <= outID.valency; j++) {
						substituentBR.addOutID(outID.id, 1, outID.setExplicitly);
					}
					outID.valency=1;
				}
			}
			int numberOfOutIDs =substituentBR.getOutIDCount();
			if (numberOfOutIDs > words.size()-1){//something like isophthaloyl chloride (more precisely written isophthaloyl dichloride)
				if (words.size()-1 != 1){
					throw new StructureBuildingException("Incorrect number of functional class groups found to balance outIDs");
				}
				int diff =numberOfOutIDs - (words.size()-1);
				Element elementToBeCloned =words.get(words.size()-1);
				for (int i = 0; i < diff; i++) {
					Element clone =state.fragManager.cloneElement(elementToBeCloned, state);
					XOMTools.insertAfter(elementToBeCloned, clone);
					elementToBeCloned=clone;
				}
				words = molecule.getChildElements();
			}
			for (int i = 1; i <= numberOfOutIDs; i++) {
				BuildResults ideBR =resolveWordOrBracket(state, words.get(i), null, new LinkedHashSet<Fragment>());
				moleculeBuildResults.mergeBuildResults(ideBR);
				Atom ideAtom =ideBR.getMainFragment().getAtomByIDOrThrow(ideBR.getMainFragment().getDefaultInID());
				Atom subAtom=substituentBR.getOutAtomTakingIntoAccountWhetherSetExplicitly(0);
				state.fragManager.attachFragments(ideAtom, subAtom, 1);
				substituentBR.removeOutID(0);
				ideAtom.setCharge(0);
			}

			for (int i = numberOfOutIDs +1; i <words.size(); i++) {//resolve anything else
				moleculeBuildResults.mergeBuildResults(resolveWordOrBracket(state, words.get(i), null, new LinkedHashSet<Fragment>()));
			}
		}
		else if(wordRule.equals("monovalentLiteralFunctionalGroup")) {
			if (words.size()!=2){
				throw new StructureBuildingException("Don't alter wordRules.xml without checking the consequences!");
			}
			String smilesOfGroup = null;
			String functionalGroupName =words.get(1).getValue();
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
			moleculeBuildResults =resolveWordOrBracket(state, words.get(0), null, new LinkedHashSet<Fragment>());//the group
			for (int i = 0; i < moleculeBuildResults.getOutIDCount(); i++) {
				Atom outAtom =moleculeBuildResults.getOutAtomTakingIntoAccountWhetherSetExplicitly(i);
				Fragment rGroup =state.fragManager.buildSMILES(smilesOfGroup, "divalentgroup", "none");
				if (moleculeBuildResults.getOutID(i).valency !=1){
					throw new StructureBuildingException("OutID has unexpected valency. Expected 1. Actual: " + moleculeBuildResults.getOutID(i).valency);
				}
				state.fragManager.attachFragments(outAtom, rGroup.getAtomByIDOrThrow(rGroup.getIdOfFirstAtom()), 1);
			}
			moleculeBuildResults.removeAllOutIDs();
		}
		else if(wordRule.equals("divalentLiteralFunctionalGroup")) {
			String smilesOfGroup = null;
			String functionalGroupName =words.get(words.size()-1).getValue();
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

			Fragment rGroup =state.fragManager.buildSMILES(smilesOfGroup, "simpleGroup", "functionalGroup", "none");
			if (words.size()!=2 && words.size()!=3){
				throw new StructureBuildingException("Unexpected number of words. Expected 2 or 3, found: " + words.size());
			}

			if (words.size()==2){//e.g. methyl sulfoxide rather than dimethyl sulfoxide
				Element clone = state.fragManager.cloneElement(words.get(0), state);
				XOMTools.insertAfter(words.get(0), clone);
				words = molecule.getChildElements();
			}

			moleculeBuildResults =resolveWordOrBracket(state, words.get(0), null, new LinkedHashSet<Fragment>());
			if (moleculeBuildResults.getOutIDCount()!=1){
				throw new StructureBuildingException("Excpected one out outID. Found " + moleculeBuildResults.getOutIDCount() );
			}
			if (moleculeBuildResults.getOutID(0).valency !=1){
				throw new StructureBuildingException("OutID has unexpected valency. Expected 1. Actual: " + moleculeBuildResults.getOutID(0).valency);
			}
			Atom outAtom =moleculeBuildResults.getOutAtomTakingIntoAccountWhetherSetExplicitly(0);
			moleculeBuildResults.removeOutID(0);
			if (rGroup.getOutIDs().size()==1){//c.f. peroxide where it is a linker
				state.fragManager.attachFragments(outAtom, rGroup.getAtomByIDOrThrow(rGroup.getOutID(0).id), 1);
				rGroup.removeOutID(0);
			}
			else{
				state.fragManager.attachFragments(outAtom, rGroup.getAtomByIDOrThrow(rGroup.getIdOfFirstAtom()), 1);
			}

			moleculeBuildResults.mergeBuildResults(resolveWordOrBracket(state, words.get(1), null, new LinkedHashSet<Fragment>()));
			if (moleculeBuildResults.getOutIDCount()!=1){
				throw new StructureBuildingException("Excpected one out outID. Found " + moleculeBuildResults.getOutIDCount() );
			}
			if (moleculeBuildResults.getOutID(0).valency !=1){
				throw new StructureBuildingException("OutID has unexpected valency. Expected 1. Actual: " + moleculeBuildResults.getOutID(0).valency);
			}
			outAtom =moleculeBuildResults.getOutAtomTakingIntoAccountWhetherSetExplicitly(0);
			moleculeBuildResults.removeOutID(0);
			state.fragManager.attachFragments(outAtom, rGroup.getAtomByIDOrThrow(rGroup.getIdOfFirstAtom()), 1);
		}
		else if(wordRule.equals("glycol")) {
			if (words.size()!=2){
				throw new StructureBuildingException("Don't alter wordRules.xml without checking the consequences!");
			}
			moleculeBuildResults =resolveWordOrBracket(state, words.get(0), null, new LinkedHashSet<Fragment>());//the group
			if (moleculeBuildResults.getOutIDCount()!=2){
				throw new StructureBuildingException("Glycol class names (e.g. ethylene glycol) expect two outIDs. Found: " + moleculeBuildResults.getOutIDCount() );
			}
			for (int i = 0; i < moleculeBuildResults.getOutIDCount(); i++) {
				Atom outAtom =moleculeBuildResults.getOutAtomTakingIntoAccountWhetherSetExplicitly(i);
				Fragment glycol =state.fragManager.buildSMILES("O", "glycol", "none");
				if (moleculeBuildResults.getOutID(i).valency !=1){
					throw new StructureBuildingException("OutID has unexpected valency. Expected 1. Actual: " + moleculeBuildResults.getOutID(i).valency);
				}
				state.fragManager.attachFragments(outAtom, glycol.getAtomByIDOrThrow(glycol.getIdOfFirstAtom()), 1);
			}
			moleculeBuildResults.removeAllOutIDs();
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
		else if(wordRule.equals("binaryOrOther")) {
			moleculeBuildResults=resolveWordOrBracket(state, words.get(words.size()-1), null, new LinkedHashSet<Fragment>());
			for (int i = words.size()-2; i >=0; i--) {//allows resolving of counter ions first
				moleculeBuildResults.mergeBuildResults(resolveWordOrBracket(state, words.get(i), null, new LinkedHashSet<Fragment>()));
			}*/

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

	private void buildEster(BuildState state, Elements words) throws StructureBuildingException {
		int i = resolvePreceedingFullWordsAndReturnNonFullIndice(state, words);
		Element currentWord=words.get(i);
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
			currentWord=words.get(++i);
		}
		
		if (i == words.size() || !words.get(i).getAttributeValue("type").equals("full")){
			throw new StructureBuildingException("Full word not found where full word expected: missing ate group in ester");
		}
		resolveWordOrBracket(state, words.get(i));
		BuildResults mainGroupBr = new BuildResults(state, words.get(i));
		i++;
		

		int esterIdCount = mainGroupBr.getFunctionalIDCount();
		int outIDCount =substituentsBr.getOutIDCount();
		if (outIDCount > esterIdCount){
			throw new StructureBuildingException("There are more radicals in the substituents than there are places to form esters");
		}
		for(i=0; i< outIDCount; i++) {
			Atom ateAtom;
			if (substituentsBr.getFirstOutID().locant!=null){
				ateAtom =determineFunctionalAtomToUse(substituentsBr.getFirstOutID().locant, mainGroupBr);
			}
			else{
				ateAtom =mainGroupBr.getFunctionalAtom(0);
				mainGroupBr.removeFunctionalID(0);
			}
			state.fragManager.attachFragments(ateAtom,substituentsBr.getOutAtomTakingIntoAccountWhetherSetExplicitly(0), 1);
			substituentsBr.removeOutID(0);
			ateAtom.setCharge(0);
		}
		resolveTrailingFullWords(state, words, i);
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
	 * Finds the first acceptable sub or root child of the element given (which is a bracket). Any brackets that are encountered are recursively enumerated
	 * @param state 
	 * @param currentElement The bracket from which to search the children of
	 * @param roots Used to check whether the element is unacceptable
	 * @param elementsUsed Used to check whether the element is unacceptable
	 * @return element The element found or null
	 * @throws StructureBuildingException 
	 */
	private Element getFirstUnusedSubOrRootFromBracket(BuildState state, Element currentElement, ArrayList<Element> roots, ArrayList<Element> elementsUsed) throws StructureBuildingException {
		Elements bracketChildren =currentElement.getChildElements();
		if (bracketChildren.size() > 0){
			for (int i = 0; i < bracketChildren.size(); i++) {
				Element bracketChild = bracketChildren.get(i);
				if (bracketChild.getLocalName().equals("substituent") || bracketChild.getLocalName().equals("bracket") || bracketChild.getLocalName().equals("root")){
					if (bracketChild.getLocalName().equals("bracket") ){
						if (bracketChild.getAttribute("type")==null){
							Element foundEl = getFirstUnusedSubOrRootFromBracket(state, bracketChild, roots, elementsUsed);
							if (foundEl!=null){
								return foundEl;
							}
						}
					}
					else{
						if (!roots.contains(currentElement) && !elementsUsed.contains(currentElement)){
							Element group =bracketChild.getFirstChildElement("group");
							if (group.getAttribute("isAMultiRadical")!=null || state.xmlFragmentMap.get(group).getOutIDs().size()==0){//you want either a linker or a terminal. e.g. oxy, nitrilo etc. or a ender e.g. phenol, benzene, ethanol.
								return bracketChild;
							}
						}
					}
				}
			}
			return null;
		}
		else{
			throw new StructureBuildingException("Empty bracket!");
		}
	}

	/**
	 * Checks through the groups accessible from the currentElement taking into account brackets
	 * i.e. those that it is feasible that the group of the currentElement could substitute onto
	 * @param state 
	 * @param currentElement
	 * @param locant: the locant string to check for the presence of
	 * @return The fragment with the locant, or null
	 * @throws StructureBuildingException 
	 */
	private Fragment findAlternativeFragmentWithLocant(BuildState state, Element startingElement, String locant) throws StructureBuildingException {
		Stack<Element> s = new Stack<Element>();
		s.add(startingElement);
		
		boolean doneFirstIteration =false;//check on index only done on first iteration to only get elements with an index greater than the starting element
		while (s.size()>0){
			Element currentElement =s.pop();
			if (currentElement.getLocalName().equals("group")){
				Fragment groupFrag =state.xmlFragmentMap.get(currentElement);
				if (groupFrag.hasLocant(locant)){
					return groupFrag;
				}
				continue;
			}
			Element parent = (Element)currentElement.getParent();
			List<Element> siblings = OpsinTools.findChildElementsWithTagNames(parent, new String[]{"bracket", "substituent", "root"});

			for (Element bracketOrSub : siblings) {
				if (!doneFirstIteration && parent.indexOf(bracketOrSub )<=parent.indexOf(currentElement)){
					continue;
				}
				if (bracketOrSub.getLocalName().equals("bracket")){
					s.push((Element)bracketOrSub.getChild(0));
				}
				else{
					Element group = bracketOrSub.getFirstChildElement("group");
					s.push(group);
				}
			}
			doneFirstIteration =true;
		}
		return null;
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
