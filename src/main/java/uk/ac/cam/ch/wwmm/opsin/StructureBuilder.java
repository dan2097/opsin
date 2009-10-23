package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Stack;

import uk.ac.cam.ch.wwmm.opsin.PreProcessor.OpsinMode;
import uk.ac.cam.ch.wwmm.ptclib.string.StringTools;
import uk.ac.cam.ch.wwmm.ptclib.xml.XOMTools;

import nu.xom.Attribute;
import nu.xom.Element;
import nu.xom.Elements;

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
		BuildResults moleculeBuildResults =null;//contains every atom in the molecule. If a mixture was the input then both "molecules" will be in this
		int substituentCount = OpsinTools.findChildElementsWithTagNameAndAttribute(molecule, "word", "type", "substituent").size();
		//TODO factor these into separate methods and maybe the rest of the class into a separate class
		if(wordRule.equals("ester") || wordRule.equals("diester")) { //e.g. ethyl ethanoate, dimethyl terephthalate,  methyl propanamide
			int i=0;
			Element currentWord=words.get(i);
			do{
				BuildResults substituentBuildResults =resolveWordOrBracket(state, currentWord, null, new LinkedHashSet<Fragment>());
				if (moleculeBuildResults==null){
					moleculeBuildResults =substituentBuildResults;
				}
				else{
					moleculeBuildResults.mergeBuildResults(substituentBuildResults);
				}
				if (substituentBuildResults.getOutIDCount() ==1){
					String locantForSubstituent = getLocantForSubstituent(currentWord);
					if (locantForSubstituent!=null){
						substituentBuildResults.getFirstOutID().locant=locantForSubstituent;//indexes which functional atom to connect to when there is a choice. Also can disambiguate which atom is a S in things like thioates
					}
				}
				i++;
				currentWord=words.get(i);
			}
			while (currentWord.getAttributeValue("type").equals("substituent"));

			while(i <words.size() && words.get(i).getAttributeValue("type").equals("full")) {
				moleculeBuildResults.mergeBuildResults(resolveWordOrBracket(state, words.get(i), null, new LinkedHashSet<Fragment>()));
				i++;
			}

			int esterIdCount = moleculeBuildResults.functionalIDs.size();
			int esterBondsToForm =Math.min(moleculeBuildResults.functionalIDs.size(), moleculeBuildResults.getOutIDCount());

			if (substituentCount > esterIdCount){
				throw new StructureBuildingException("Name appears to have erroneous spaces interfering with ester interpretation");
			}

			for(i=0;i<esterBondsToForm;i++) {
				Atom ateAtom;
				if (moleculeBuildResults.getFirstOutID().locant!=null){
					ateAtom =determineFunctionalAtomToUse(moleculeBuildResults.getFirstOutID().locant, moleculeBuildResults);
				}
				else{
					ateAtom =moleculeBuildResults.getFunctionalOutAtom(0);
					moleculeBuildResults.functionalIDs.remove(0);
				}
				state.fragManager.attachFragments(ateAtom,moleculeBuildResults.getOutAtomTakingIntoAccountWhetherSetExplicitly(0), 1);
				moleculeBuildResults.removeOutID(0);
				ateAtom.setCharge(0);
			}
		}
		else if (wordRule.equals("functionalClassEster")){//e.g. ethanoic acid ethyl ester, tetrathioterephthalic acid dimethyl ester
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
		else if(wordRule.equals("polymer")) {
			for (Element word : state.firstMultiRadical.keySet()) {//don't use multiplicative nomenclature
				state.firstMultiRadical.put(word, null);
			}
			Element polymer = molecule.getFirstChildElement("word");
			moleculeBuildResults =resolveWordOrBracket(state, polymer, null, new LinkedHashSet<Fragment>());
			int inIdCount=0;
			Fragment firstFragment =null;
			for (Fragment frag : moleculeBuildResults.fragments) {
				inIdCount +=frag.getInIDs().size();
				if (frag.getInIDs().size()!=0){
					firstFragment =frag;
				}
			}
			if (moleculeBuildResults.getOutIDCount() ==2  && inIdCount ==0){//e.g. poly(ethylene)
				Atom inAtom =moleculeBuildResults.getOutAtomTakingIntoAccountWhetherSetExplicitly(0);
				Fragment rGroup =state.fragManager.buildSMILES("[Xe]");
				state.fragManager.attachFragments(inAtom, rGroup.getAtomByIDOrThrow(rGroup.getIdOfFirstAtom()), moleculeBuildResults.getOutID(0).valency);
				Atom outAtom =moleculeBuildResults.getOutAtomTakingIntoAccountWhetherSetExplicitly(1);
				rGroup =state.fragManager.buildSMILES("[Rn]");
				state.fragManager.attachFragments(outAtom, rGroup.getAtomByIDOrThrow(rGroup.getIdOfFirstAtom()), moleculeBuildResults.getOutID(1).valency);
				moleculeBuildResults.removeAllOutIDs();
			}
			else if (moleculeBuildResults.getOutIDCount() ==1 && inIdCount==1){//other cases
				OutID inId = firstFragment.getInID(0);
				Atom inAtom = firstFragment.getAtomByIdOrNextSuitableAtomOrThrow(inId.id, inId.valency);
				Fragment rGroup =state.fragManager.buildSMILES("[Xe]");
				state.fragManager.attachFragments(inAtom, rGroup.getAtomByIDOrThrow(rGroup.getIdOfFirstAtom()), inId.valency);
				firstFragment.removeInID(0);

				Atom outAtom =moleculeBuildResults.getOutAtomTakingIntoAccountWhetherSetExplicitly(0);
				rGroup =state.fragManager.buildSMILES("[Rn]");
				state.fragManager.attachFragments(outAtom, rGroup.getAtomByIDOrThrow(rGroup.getIdOfFirstAtom()), moleculeBuildResults.getOutID(0).valency);
				moleculeBuildResults.removeOutID(0);
			}
			else{
				throw new StructureBuildingException("Polymer building failed: Two termini were not found; Expected 1 outIDs, found: " +moleculeBuildResults.getOutIDCount() +" ,expected 1 inIDs, found: " +inIdCount);
			}
		}
		else if(wordRule.equals("binaryOrOther")) {
			moleculeBuildResults=resolveWordOrBracket(state, words.get(words.size()-1), null, new LinkedHashSet<Fragment>());
			for (int i = words.size()-2; i >=0; i--) {//allows resolving of counter ions first
				moleculeBuildResults.mergeBuildResults(resolveWordOrBracket(state, words.get(i), null, new LinkedHashSet<Fragment>()));
			}
		}else if(wordRule.equals("simple") || wordRule.equals("acid")) {
			Element word = molecule.getFirstChildElement("word");
			moleculeBuildResults =resolveWordOrBracket(state, word, null, new LinkedHashSet<Fragment>());
		}
		else{
			throw new StructureBuildingException("Unknown Word Rule");
		}
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
		if (moleculeBuildResults.getOutIDCount()>0){
			throw new StructureBuildingException("Radicals are currently set to not convert to structures");
		}
		return uniFrag;
	}

	/**Resolves the contents of a &lt;word&gt; or &lt;bracket&gt; tag, recursively.
	 * @param state 
	 *
	 * @param wob The &lt;word&gt; or &lt;bracket&gt; tag.
	 * @param parentFrag The fragment that the results are to be attached to. May be null.
	 * @param fragmentsToAvoid Usually empty list of fragments that should not be joined to their parentFrag
	 * @return A buildResults for the results of building.
	 * @throws StructureBuildingException If the contents won't build properly.
	 */
	BuildResults resolveWordOrBracket(BuildState state, Element wob, Fragment parentFrag, LinkedHashSet<Fragment> fragmentsToAvoid) throws StructureBuildingException {
		if (parentFrag==null && wob.getLocalName().equals("word")){
			if (state.firstMultiRadical.get(wob) !=null){
				ArrayList<Element> elements =new ArrayList<Element>();
				Element root =state.firstMultiRadical.get(wob);
				state.firstMultiRadical.remove(wob);
				elements.add(root);
				BuildResults buildResults = resolveRootOrSubstituent(state, root, null, new ArrayList<Element>(), fragmentsToAvoid);
				BuildResults buildResultsMultiplicative = resolveWordOrBracketLeftToRight(state, buildResults, elements);
				buildResultsMultiplicative.mergeBuildResults(resolveWordOrBracket(state, wob, null, buildResultsMultiplicative.fragments));

				List<Element> groups =OpsinTools.findDescendantElementsWithTagName(wob, "group");
				for (Element group : groups) {
					if (!buildResults.fragments.contains(state.xmlFragmentMap.get(group))){
						throw new StructureBuildingException("OPSIN's multiplicative nomenclature implementation has failed to resolve all groups!");
					}
				}
				return buildResultsMultiplicative;
			}
		}
		String locantStr = "0";

		if(wob.getAttribute("locant")!= null) {
			locantStr = wob.getAttributeValue("locant");
		}
		Element root = wob.getFirstChildElement("root");
		List<Element> bracketsOrSubs = OpsinTools.findChildElementsWithTagNames(wob, new String[]{"bracket", "substituent"});
		BuildResults buildResults = null;
		if(root == null) {
			Element lastBracketOrSub = bracketsOrSubs.get(bracketsOrSubs.size()-1);
			if (lastBracketOrSub.getLocalName().equals("bracket")){
				//if this is a bracket things get complicated as we need to have resolved everything in this recursively,
				//starting from the last deepest bracket. The criteria of being the last element takes priority over depth
				while (lastBracketOrSub.getLocalName().equals("bracket")){
					bracketsOrSubs = OpsinTools.findChildElementsWithTagNames(lastBracketOrSub, new String[]{"bracket", "substituent", "root"});
					lastBracketOrSub = bracketsOrSubs.get(bracketsOrSubs.size()-1);
				}
				root =lastBracketOrSub;
				buildResults = resolveRootOrSubstituent(state, lastBracketOrSub, null, new ArrayList<Element>(), fragmentsToAvoid );
				do {
					lastBracketOrSub = (Element) lastBracketOrSub.getParent();
					bracketsOrSubs = OpsinTools.findChildElementsWithTagNames(lastBracketOrSub, new String[]{"bracket", "substituent", "root"});
					for (int i = bracketsOrSubs.size() -2; i >=0; i--) {//don't use last element as that's the root
						Element currentChild= bracketsOrSubs.get(i);
						if (currentChild.getLocalName().equals("bracket")){
							buildResults.mergeBuildResults(resolveWordOrBracket(state, currentChild, buildResults.getMainFragment(), fragmentsToAvoid));
						}
						else{
							buildResults.mergeBuildResults(resolveRootOrSubstituent(state, currentChild , buildResults.getMainFragment(), new ArrayList<Element>(), fragmentsToAvoid));
						}
					}
				} while(!lastBracketOrSub.equals(wob));
			}
			else{
				root =lastBracketOrSub;
			}
		}
		if (buildResults==null){
			buildResults = resolveRootOrSubstituent(state, root, null, bracketsOrSubs, fragmentsToAvoid);
			Elements roots= wob.getChildElements("root");
			for (int i = 1; i < roots.size(); i++) {//special case where there are multiple roots e.g. naphthalene-1,5-diyl 4,4'-bis(2,5-dioxo-1H-pyrrol-1-yl)dibenzoate (this is not technically multiplicative)
				buildResults.mergeBuildResults(resolveRootOrSubstituent(state, roots.get(i), null, new ArrayList<Element>(), fragmentsToAvoid));
			}
		}

		if(parentFrag != null && !fragmentsToAvoid.contains(buildResults.getMainFragment())) {
			joinFragments(state, buildResults, (Element) root.getParent(), parentFrag, locantStr);
		}
		return buildResults;
	}


	/**Resolves the contents of a &lt;word&gt; or &lt;bracket&gt; tag, recursively BUT from left to right
	 * This is utilised by multiplicative nomenclature
	 * @param state
	 * @param buildResults The current building progress (this function is called recursively)
	 * @param roots The elements upon which building is occuring. Say for ethylenediaminetetraacetic acid this would initially contain ethylene, then when it's run next would contain both amines
	 * @return A buildResults for the results of building
	 * @throws StructureBuildingException If the contents won't build properly.
	 */
	private BuildResults resolveWordOrBracketLeftToRight(BuildState state, BuildResults buildResults, ArrayList<Element> roots) throws StructureBuildingException {
		ArrayList<Element> elementsUsed =new ArrayList<Element>();
		int numberOfOutIDs =buildResults.getOutIDCount();
		if (numberOfOutIDs <=1){//no longer multiplicative nomenclature (should be the end of the name...)
			if (numberOfOutIDs==1){
				for (Element root : roots) {//I think if this code is run something has gone wrong!
					Element group=root.getFirstChildElement("group");
					if (group!=null){
						buildResults.fragments.remove(state.xmlFragmentMap.get(group));
					}
				}
			}
			return buildResults;
		}

		int outIDsPerRoot =numberOfOutIDs/roots.size();
//		System.out.println("number of roots: " + roots.size());
//		System.out.println("outIdsPerRoot " + outIDsPerRoot);
		int extraOutIDs = numberOfOutIDs % roots.size();//if this is not 0 then the resulting structure will be a radical
		if (extraOutIDs!=0){
			throw new StructureBuildingException("Unhandled use of multiplicative nomenclature");
		}

		//resolve left to right
		mainLoop: for (int i = 0; i < roots.size(); i++) {
			Element root =roots.get(i);
			int outIDsForThisRoot =outIDsPerRoot;
			//System.out.println("Root " + i);
			Element parent =(Element) root.getParent();

			int indexOfRoot =parent.indexOf(root);
			Elements children =parent.getChildElements();
			ArrayList<Element> potentialSubsAndRoots =new ArrayList<Element>();
			for (int j = 0; j < children.size(); j++) {
				Element currentElement = children.get(j);
				if (j > indexOfRoot && (currentElement.getLocalName().equals("substituent") || currentElement.getLocalName().equals("bracket") || currentElement.getLocalName().equals("root"))
						&& !roots.contains(currentElement) && !elementsUsed.contains(currentElement)){
					
					if (currentElement.getLocalName().equals("bracket")){
						if (currentElement.getAttribute("type")==null){//don't want to explore implicit brackets
							Element foundEl =getFirstUnusedSubOrRootFromBracket(state, currentElement, roots, elementsUsed);
							if (foundEl!=null){
								potentialSubsAndRoots.add(foundEl);
							}
						}
					}
					else{
						Element group =currentElement.getFirstChildElement("group");
						if (group.getAttribute("isAMultiRadical")!=null || state.xmlFragmentMap.get(group).getOutIDs().size()==0){//you want either a linker or a terminal. e.g. oxy, nitrilo etc. or a ender e.g. phenol, benzene, ethanol.
							potentialSubsAndRoots.add(currentElement);
						}
					}
				}
			}
			if (potentialSubsAndRoots.size() ==0){break;}

//			System.out.println(root.toXML());
//			System.out.println("number of available options: " + potentialSubsAndRoots.size());


			//TODO Do this properly e.g. only do it in specific cases e.g. methylenecyclohexane, not methandiylcyclohexane
//			int bondOrder =outIDsPerRoot/potentialBracketsOrSubs.size();
//			if (potentialBracketsOrSubs.size() < outIDsForThisRoot){//some of the outID need to be removed and others have their valency increased
//				if (outIDsPerRoot % potentialBracketsOrSubs.size() !=0){
//					throw new StructureBuildingException("Unhandled use of multiplicative nomenclature");
//				}
//				int valencyCount=0;//how much valency has been changed
//				ArrayList<OutID> outIDsToRemove =new ArrayList<OutID>();
//				for (OutID outId : buildResults.outIDs) {
//					if (valencyCount>0 && outId.valency==1){
//						outIDsToRemove.add(outId);
//						valencyCount--;
//					}
//					if (outId.valency<bondOrder){
//						valencyCount= valencyCount+(bondOrder-outId.valency);
//						outId.valency=bondOrder;
//					}
//				}
//				buildResults.outIDs.removeAll(outIDsToRemove);
//			}

			for (int j = 0; j < potentialSubsAndRoots.size(); j++) {
				Element currentSubOrRoot = potentialSubsAndRoots.get(j);
				if (outIDsForThisRoot >=1){
					//System.out.println("current: " + currentSubOrRoot.toXML());
					if(!buildResults.fragments.contains(state.xmlFragmentMap.get(currentSubOrRoot.getFirstChildElement("group")))){
						BuildResults tempChildBuildResults = resolveRootOrSubstituent(state, currentSubOrRoot , null, new ArrayList<Element>(), new LinkedHashSet<Fragment>());
						joinFragmentsMultiplicative(state, tempChildBuildResults, buildResults);
						buildResults.mergeBuildResults(tempChildBuildResults);
						elementsUsed.add(currentSubOrRoot);
						outIDsForThisRoot--;
					}
				}
				else{
					continue mainLoop;
				}
			}
		}
		if (elementsUsed.size()!=0){
			return resolveWordOrBracketLeftToRight(state, buildResults, elementsUsed);
		}

		if (numberOfOutIDs !=0 && roots.size() >0){
			for (int i = 0; i < roots.size(); i++) {
				Element parent =(Element) roots.get(i).getParent();
				if (parent.getLocalName().equals("bracket")){
					roots.set(i, parent);
				}
				else{
					return buildResults;
					//as it isn't in a bracket all of the molecule should have been evaluated, so the only conclusion is that the molecule is a radical
				}
			}
			return resolveWordOrBracketLeftToRight(state, buildResults, roots);
		}
		return buildResults;
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
	

	/**Resolves the contents of a &lt;root&gt; or &lt;substituent&gt; tag.
	 * @param state 
	 * @param rootOrSub The &lt;root&gt; or &lt;substituent&gt; tag.
	 * @param parentFrag The fragment to which the built fragment will be attached. May be null.
	 * @param bracketsOrSubsNodes The elements which should be evaluated and connected to the fragment built from rootOrSub. Can be empty
	 * @param fragmentsToAvoid Usually empty list of fragments that should not be joined to their parentFrag
	 * @return A buildResults for the results of building.
	 * @throws StructureBuildingException If the contents won't build properly.
	 */
	BuildResults resolveRootOrSubstituent(BuildState state, Element rootOrSub, Fragment parentFrag, List<Element> bracketsOrSubs, LinkedHashSet<Fragment> fragmentsToAvoid) throws StructureBuildingException {
		Element group = rootOrSub.getFirstChildElement("group");
		Fragment thisFrag = state.xmlFragmentMap.get(group);

		ArrayList<Element> unsaturators = new ArrayList<Element>();
		ArrayList<Element> heteroatoms = new ArrayList<Element>();
		ArrayList<Element> hydrogenElements = new ArrayList<Element>();

		Elements children =rootOrSub.getChildElements();
		for (int i = 0; i < children.size(); i++) {
			Element currentEl =children.get(i);
			String elName =currentEl.getLocalName();
			if (elName.equals("unsaturator")){
				unsaturators.add(currentEl);
			}
			else if (elName.equals("heteroatom")){
				heteroatoms.add(currentEl);
			}
			else if (elName.equals("hydro")){
				hydrogenElements.add(currentEl);
			}
			else if (elName.equals("hydrogen")){
				hydrogenElements.add(currentEl);
			}
			else if (elName.equals("indicatedHydrogen")){
				hydrogenElements.add(currentEl);
			}
		}

		/*
		 *OutID queries from this point onwards should be directed to the buildResults
		 */
		BuildResults buildResults = new BuildResults(thisFrag);
		if (fragmentsToAvoid.contains(buildResults.getMainFragment())){
			buildResults.functionalIDs.clear();
			buildResults.removeAllOutIDs();
		}

		int idOfFirstAtomInFragment= thisFrag.getIdOfFirstAtom();
		//used if locants are not specified; whenever used it is incremented
		//all non-suffix atoms in the fragment are eventually checked if the defaultId does not correspond to a suitable atom
		//e.g. if the id was 3 in a 5 atom fragment which had ids 1-5, atoms would be checked for suitability in the order 3,4,5,1,2
		int defaultId = idOfFirstAtomInFragment;

		/*
		 * Add locanted functionality
		 */
		for(int i=hydrogenElements.size() -1;i >= 0;i--) {
			Element hydrogen = hydrogenElements.get(i);
			String locant = getLocant(hydrogen);
			if(!locant.equals("0")) {
				thisFrag.getAtomByLocantOrThrow(locant).subtractSpareValency(1);
				hydrogenElements.remove(hydrogen);
			}
			hydrogen.detach();
		}

		for(int i=unsaturators.size() -1;i >= 0;i--) {
			Element unsaturator = unsaturators.get(i);
			String locant = getLocant(unsaturator);
			int bondOrder = Integer.parseInt(unsaturator.getAttributeValue("value"));
			if(bondOrder <= 1) {
				continue;
			}
			if(!locant.equals("0")){
				unsaturators.remove(unsaturator);
				Integer idOfFirstAtomInMultipleBond=thisFrag.getIDFromLocantOrThrow(locant);
				if (unsaturator.getAttribute("compoundLocant")!=null){
					state.fragManager.unsaturate(idOfFirstAtomInMultipleBond, unsaturator.getAttributeValue("compoundLocant"), bondOrder, thisFrag);
				}
				else{
					state.fragManager.unsaturate(idOfFirstAtomInMultipleBond, bondOrder, thisFrag);
				}
			}
			unsaturator.detach();
		}

		for(int i=heteroatoms.size() -1;i >= 0;i--) {
			Element heteroatom = heteroatoms.get(i);
			String locant = getLocant(heteroatom);
			String atomSMILES = heteroatom.getAttributeValue("value");
			if(!locant.equals("0")) {
				state.fragManager.makeHeteroatom(thisFrag.getAtomByLocantOrThrow(locant), atomSMILES, true);
				heteroatoms.remove(heteroatom);
			}
			heteroatom.detach();
		}

		for(int i=bracketsOrSubs.size() -1;i >= 0;i--) {
			Element bracketOrSub = bracketsOrSubs.get(i);
			if (bracketOrSub.equals(rootOrSub)){continue;}
			if(bracketOrSub.getAttribute("locant") != null) {
				//slim possibility does exist that the locant is not referring to this fragment
				if (thisFrag.getAtomByLocant(bracketOrSub.getAttributeValue("locant"))!=null){
					if (bracketOrSub.getLocalName().equals("bracket")){
						buildResults.mergeBuildResults(resolveWordOrBracket(state, bracketOrSub, thisFrag, fragmentsToAvoid));
					}
					else{
						buildResults.mergeBuildResults(resolveRootOrSubstituent(state, bracketOrSub , thisFrag, new ArrayList<Element>(), fragmentsToAvoid));
					}
					bracketsOrSubs.remove(i);
				}
			}
		}

		if(parentFrag != null) {
			/* attach the Fragment that has just been built to the parentFrag */
			if(!fragmentsToAvoid.contains(buildResults.getMainFragment()) && rootOrSub.getAttribute("locant") != null) {
				String locantStr =rootOrSub.getAttributeValue("locant");
				joinFragments(state, buildResults, rootOrSub, parentFrag, locantStr);
			}
		}

		/*
		 * Add unlocanted functionality
		 */	

		if (hydrogenElements.size()>0){
			/*
			 * This function is not entirely straightforward as certain atoms definitely should have their spare valency reduced
			 * However names are not consistent as to whether they bother having the hydro tags do this!
			 * The atoms in atomsWithSV are in atom order those that can take a hydro element and then those that shouldn't really take a hydro element as it's absence is unambiguous
			 */
			LinkedList<Atom> atomsWithSV = new LinkedList<Atom>();
			LinkedList<Atom> atomsWhichImplicitlyWillHaveTheirSVRemoved = new LinkedList<Atom>();
			List<Atom> atomList =thisFrag.getAtomList();
			for (Atom atom : atomList) {
				if (atom.getType().equals("suffix")){
					break;
				}
				atom.ensureSVIsConsistantWithValency(false);//doesn't take into account suffixes
				if (atom.getSpareValency() >=1){
					if (atom.getNote("OneSuffixAttached")!=null){
						atomsWhichImplicitlyWillHaveTheirSVRemoved.add(atom);
					}
					else{
						atomsWithSV.add(atom);
					}
				}
			}
			atomsWithSV.addAll(atomsWhichImplicitlyWillHaveTheirSVRemoved);//these end up at the end of the list
			if (hydrogenElements.size()> atomsWithSV.size()){
				throw new StructureBuildingException("Cannot find atom to add hydrogen to (" +
						hydrogenElements.size() + " hydrogen adding tags but only " +  atomsWithSV.size() +" positions that can be hydrogenated)" );
			}
			for(int j=0;j<hydrogenElements.size();j++) {
				Atom atomToReduceSpareValencyOn=atomsWithSV.removeFirst();
				atomToReduceSpareValencyOn.subtractSpareValency(1);
				hydrogenElements.get(j).detach();
			}
		}

		for(int j=0;j<unsaturators.size();j++) {
			Element unsaturator = unsaturators.get(j);
			int bondOrder = Integer.parseInt(unsaturator.getAttributeValue("value"));
			if(bondOrder <= 1) {
				continue;
			}
			//checks if both atoms can accept an extra bond (if double bond) or two extra bonds (if triple bond)

			Atom currentAtom =thisFrag.getAtomByIDOrThrow(defaultId);
			Atom nextAtom =thisFrag.getAtomByIDOrThrow(defaultId +1);
			while (currentAtom.getSpareValency() != 0 || ValencyChecker.checkValencyAvailableForBond(currentAtom, bondOrder-1 + currentAtom.getOutValency()) != true ||
					nextAtom.getSpareValency() != 0 || ValencyChecker.checkValencyAvailableForBond(nextAtom, bondOrder-1 + nextAtom.getOutValency()) != true){
				defaultId++;
				currentAtom =thisFrag.getAtomByIDOrThrow(defaultId);
				nextAtom =thisFrag.getAtomByIDOrThrow(defaultId +1);
				if (currentAtom.getType().equals("suffix") || nextAtom.getType().equals("suffix")){
					throw new StructureBuildingException("No suitable atom found");
				}
			}
			Integer idOfFirstAtomInMultipleBond=currentAtom.getID();
			if (unsaturator.getAttribute("compoundLocant")!=null){
				state.fragManager.unsaturate(idOfFirstAtomInMultipleBond, unsaturator.getAttributeValue("compoundLocant"), bondOrder, thisFrag);
			}
			else{
				state.fragManager.unsaturate(idOfFirstAtomInMultipleBond, bondOrder, thisFrag);
			}
			defaultId=idOfFirstAtomInMultipleBond +2;
			unsaturator.detach();
		}
		defaultId = idOfFirstAtomInFragment;

		for(int j=0;j<heteroatoms.size();j++) {
			Element heteroatom = heteroatoms.get(j);
			String atomSMILES = heteroatom.getAttributeValue("value");
			//finds an atom for which changing it to the specified heteroatom will not cause valency to be violated
			Atom atomToReplaceWithHeteroAtom=thisFrag.getAtomByIDOrThrow(defaultId);
			while (ValencyChecker.checkValencyAvailableForReplacementByHeteroatom(atomToReplaceWithHeteroAtom, atomSMILES) != true){
				defaultId++;
				atomToReplaceWithHeteroAtom=thisFrag.getAtomByIDOrThrow(defaultId);
				if (atomToReplaceWithHeteroAtom.getType().equals("suffix")){
					throw new StructureBuildingException("No suitable atom found");
				}
			}
			state.fragManager.makeHeteroatom(atomToReplaceWithHeteroAtom, atomSMILES, true);
			defaultId++;
			heteroatom.detach();
		}
		defaultId = idOfFirstAtomInFragment;

		if (thisFrag.getOutIDs().size()>0){//assign any outIDs that have not been set to a specific atom to a specific atom
			for (OutID outID : thisFrag.getOutIDs()) {
				if (!outID.setExplicitly){
					defaultId=outID.id;
					Atom atomToAssociateOutIDWith=thisFrag.getAtomByIDOrThrow(defaultId);
					while (ValencyChecker.checkValencyAvailableForBond(atomToAssociateOutIDWith, atomToAssociateOutIDWith.getSpareValency() + atomToAssociateOutIDWith.getOutValency() +outID.valency) != true){
						defaultId++;
						atomToAssociateOutIDWith=thisFrag.getAtomByIDOrThrow(defaultId);
						if (atomToAssociateOutIDWith.getType().equals("suffix")){
							throw new StructureBuildingException("No suitable atom found");
						}
					}
					outID.id=defaultId;
					outID.setExplicitly=true;
					atomToAssociateOutIDWith.addOutValency(outID.valency);
					defaultId = idOfFirstAtomInFragment;
				}
			}
		}

		for(int i=bracketsOrSubs.size() -1;i >= 0;i--) {
			Element bracketOrSub = bracketsOrSubs.get(i);
			if (bracketOrSub.equals(rootOrSub)){continue;}
			if (bracketOrSub.getLocalName().equals("bracket")){
				buildResults.mergeBuildResults(resolveWordOrBracket(state, bracketOrSub, thisFrag, fragmentsToAvoid));
			}
			else{
				buildResults.mergeBuildResults(resolveRootOrSubstituent(state, bracketOrSub , thisFrag, new ArrayList<Element>(), fragmentsToAvoid));
			}
		}

		if(parentFrag != null) {
			/* attach the Fragment that has just been built to the parentFrag */
			if(!fragmentsToAvoid.contains(buildResults.getMainFragment()) && rootOrSub.getAttribute("locant") == null) {
				joinFragments(state, buildResults, rootOrSub, parentFrag, "0");
			}
		}

		return buildResults;
	}

	/**
	 * Used to join two fragments. If the locant is not present on the parent alternative feasible fragments will be accessed to
	 * try to find a fragment with a suitable locant which assumedly was the fragment to which bonding was intended
	 * @param state 
	 * @param fragToBeJoinedBuildResults: BuildResults for the fragment to be joined
	 * @param parentSubstituentOrBracketOfFrag: The XML corresponding to the substituent/bracket this fragment is in
	 * This is a bracket if called from resolveWordOrBracket or a substituent if called from resolveRootOrSubstituent
	 * @param parentFrag: The parent fragment to attach to
	 * @param locantStr: A locant on this fragment
	 * @throws StructureBuildingException
	 */
	private void joinFragments(BuildState state, BuildResults fragToBeJoinedBuildResults, Element parentSubstituentOrBracketOfFrag,
			Fragment parentFrag, String locantStr) throws StructureBuildingException {

		Fragment fragToBeJoined =fragToBeJoinedBuildResults.getMainFragment();
		boolean polymerMode =false;
		if (state.mode.equals(OpsinMode.poly) && fragToBeJoinedBuildResults.getOutIDCount() ==2){
			polymerMode =true;
		}
		if (fragToBeJoinedBuildResults.getOutIDCount() <=0){
			throw new StructureBuildingException("Fragment does not have any specified outIDs!");
		}
		Atom from;
		int bondOrder;
		if (polymerMode){
			OutID firstOutID = fragToBeJoined.getOutID(0);
			fragToBeJoined.addInID(firstOutID.id, firstOutID.valency, firstOutID.setExplicitly);
			fragToBeJoined.removeOutID(firstOutID);
			fragToBeJoinedBuildResults.removeOutID(0);
			from = fragToBeJoinedBuildResults.getOutAtom(0);
			bondOrder = fragToBeJoinedBuildResults.getOutID(0).valency;
			fragToBeJoined.removeOutID(fragToBeJoined.getOutID(0));
			fragToBeJoinedBuildResults.removeOutID(0);
		}
		else{
			from = fragToBeJoinedBuildResults.getOutAtom(0);
			bondOrder =fragToBeJoinedBuildResults.getFirstOutID().valency;
			if (fragToBeJoinedBuildResults.getFirstOutID().setExplicitly != true){//not set explicitly so may be an inappropriate atom
				from=fragToBeJoined.getAtomByIdOrNextSuitableAtomOrThrow(from.getID(), bondOrder);
			}
			fragToBeJoinedBuildResults.removeOutID(0);
		}


		Atom to =null;
		if (polymerMode){//TODO allow names like poly[oxy(dichloromethyl)methylene] e.g. search for next with two outIDs or an inId
			parentFrag = getNextInScopeMultiValentFragment(parentSubstituentOrBracketOfFrag, state);
			if (parentFrag ==null){
				throw new StructureBuildingException("Polymer building failed: No suitable fragment found to connect to");
			}
			int parentExpectedBondOrder;
			if (parentFrag.getOutIDs().size()==2){
				OutID outId =parentFrag.getOutID(0);
				to =parentFrag.getAtomByIDOrThrow(outId.id);
				parentExpectedBondOrder =outId.valency;
				outId.buildResults.outIDs.remove(outId);
				parentFrag.removeOutID(0);
			}
			else if (parentFrag.getInIDs().size()==1){
				to =parentFrag.getAtomByIDOrThrow(parentFrag.getInID(0).id);
				parentExpectedBondOrder =parentFrag.getInID(0).valency;
				parentFrag.removeInID(0);
			}
			else{
				throw new StructureBuildingException("Polymer building failed: Fragment was expected to be have two outIDs or an inId but in fact had: " +parentFrag.getOutIDs().size() +" outIDs and : " +parentFrag.getInIDs().size() +" inIDs");
			}
			if (parentExpectedBondOrder != bondOrder){//wrong outId on the frag to be joined has been utilised, expected to happen in ambiguous cases such as nitrilo
				OutID theRealOutId = fragToBeJoined.getInID(0);
				fragToBeJoined.removeInID(0);
				fragToBeJoined.addInID(from.getID(), bondOrder, true);
				from = fragToBeJoined.getAtomByIDOrThrow(theRealOutId.id);//probably the same atom as before
				bondOrder = theRealOutId.valency;
				if (parentExpectedBondOrder != bondOrder){
					throw new StructureBuildingException("Polymer building failed: bond order disagreement");
				}
			}
		}
		else{
			if(locantStr.equals("0")){
				to = parentFrag.getAtomByIdOrNextSuitableAtom(parentFrag.getDefaultInID(), bondOrder, true);
			}
			else{
				to =parentFrag.getAtomByLocant(locantStr);
			}
		}

		//case where you should actually be substituting onto the previous element e.g. 5-(4-methylphenylcarbonyl)pentane
		if (to==null){
			if(locantStr.equals("0")){
				ArrayList<Fragment> possibleParents =findAlternativeFragments(state, parentSubstituentOrBracketOfFrag);
				for (Fragment fragment : possibleParents) {
					to = fragment.getAtomByIdOrNextSuitableAtom(fragment.getDefaultInID(), bondOrder, true);
					if (to !=null){
						break;
					}
				}
				if (to  ==null){
					throw new StructureBuildingException("Cannot find fragment to reasonably attach atom with id: " +from.getID() +" to using unspecified locant");
				}
			}
			else{
				parentFrag =findAlternativeFragmentWithLocant(state, parentSubstituentOrBracketOfFrag, locantStr);
				if (parentFrag  ==null){
					throw new StructureBuildingException("Cannot find fragment to reasonably attach atom with id: " +from.getID() +" to using locant: " + locantStr);
				}
				to =parentFrag.getAtomByLocant(locantStr);
			}
		}
		state.fragManager.attachFragments(from, to, bondOrder);
	}

	/**
	 * Given a subsituent/bracket finds the next multi valent substituent/root that is in scope and hence its group
	 * e.g. for oxy(dichloromethyl)methylene given oxy substituent the methylene group would be found
	 * for oxy(dichloroethylene) given oxy substituent the ethylene group would be found
	 * for oxy(carbonylimino) given oxy carbonyl would be found
	 * @param substituentOrBracket
	 * @param state 
	 * @return frag
	 * @throws StructureBuildingException 
	 */
	private Fragment getNextInScopeMultiValentFragment(Element substituentOrBracket, BuildState state) throws StructureBuildingException {
		if (!substituentOrBracket.getLocalName().equals("substituent") && !substituentOrBracket.getLocalName().equals("bracket")){
			throw new StructureBuildingException("Input to this function should be a substituent or bracket");
		}
		if (substituentOrBracket.getParent()==null){
			throw new StructureBuildingException("substituent did not have a parent!");
		}
		Element parent =(Element) substituentOrBracket.getParent();
		
		List<Element> children = OpsinTools.findChildElementsWithTagNames(parent, new String[]{"substituent", "bracket", "root"});//will be returned in index order
		int indexOfSubstituent =parent.indexOf(substituentOrBracket);
		for (Element child : children) {
			if (parent.indexOf(child) <=indexOfSubstituent){//only want things after the input
				continue;
			}
			List<Element> childDescendants;
			if (child.getLocalName().equals("bracket")){
				childDescendants = OpsinTools.findDescendantElementsWithTagNames(child, new String[]{"substituent", "root"});//will be returned in depth-first order
			}
			else{
				childDescendants =new ArrayList<Element>();
				childDescendants.add(child);
			}
			for (Element descendantChild : childDescendants) {
				Element group = descendantChild.getFirstChildElement("group");
				if (group == null){
					throw new StructureBuildingException("substituent/root is missing its group");
				}
				Fragment possibleFrag = state.xmlFragmentMap.get(group);
				if ((possibleFrag.getOutIDs().size()==2 && possibleFrag.getInIDs().size()==0)||(possibleFrag.getOutIDs().size()==0 && possibleFrag.getInIDs().size()==1)){
					return possibleFrag;
				}
			}
		}
		return null;
	}

	/**
	 * Joins two fragments together. The main difference between this and joinFragments is that this function expects no specified locants
	 * and that in this on the parentBuildResult loses an OutID, not the fragment to be joined.
	 * @param state 
	 * @param fragToBeJoinedBuildResults
	 * @param parentFragBuildResults
	 * @throws StructureBuildingException
	 */
	private void joinFragmentsMultiplicative(BuildState state, BuildResults fragToBeJoinedBuildResults, BuildResults parentFragBuildResults) throws StructureBuildingException {
		Fragment fragToBeJoined =fragToBeJoinedBuildResults.getMainFragment();
		if (parentFragBuildResults.getOutIDCount() <=0){
			throw new StructureBuildingException("Fragment does not have any specified outIDs!");
		}
		Atom from =parentFragBuildResults.getOutAtom(0);
		int bondOrder =parentFragBuildResults.getFirstOutID().valency;
		if (parentFragBuildResults.getFirstOutID().setExplicitly != true){//not set explicitly so may be an inappropriate atom
			from=from.getFrag().getAtomByIdOrNextSuitableAtomOrThrow(from.getID(), bondOrder);
		}
		parentFragBuildResults.removeOutID(0);

		Atom to;
		if (fragToBeJoined.getInIDs().size() > 0){//inID may have been explicitly specified e.g. N,N'-(butane-2,2-diyl)bis(acrylamide) specifies the N position on the acrylamide
			to = fragToBeJoined.getAtomByIDOrThrow(fragToBeJoined.getInID(0).id);
			fragToBeJoined.removeInID(0);
		}
		else{
			 to = fragToBeJoined.getAtomByIdOrNextSuitableAtomOrThrow(fragToBeJoined.getDefaultInID(), bondOrder);
		}
		state.fragManager.attachFragments(from, to, bondOrder);
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
	 * Finds all the groups accessible from the currentElement taking into account brackets
	 * i.e. those that it is feasible that the group of the currentElement could substitute onto
	 * @param state 
	 * @param currentElement
	 * @return A list of fragments in the order to try them as possible parent fragments
	 */
	private ArrayList<Fragment> findAlternativeFragments(BuildState state, Element startingElement) {
		Stack<Element> s = new Stack<Element>();
		s.add(startingElement);
		ArrayList<Fragment> foundFragments =new ArrayList<Fragment>();
		boolean doneFirstIteration =false;//check on index only done on first iteration to only get elements with an index greater than the starting element
		while (s.size()>0){
			Element currentElement =s.pop();
			if (currentElement.getLocalName().equals("group")){
				Fragment groupFrag =state.xmlFragmentMap.get(currentElement);
				foundFragments.add(groupFrag);
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
		return foundFragments;
	}


	/**Gets the locant from a group/suffix tag, defaulting to "0"
	 *
	 * @param element
	 * @return The locant on the group/suffix tag.
	 */
	static String getLocant(Element element) {
		String locantStr = element.getAttributeValue("locant");
		if(locantStr == null) return "0";
		return locantStr;
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
	 * @param moleculeBuildResults
	 * @return functionalAtomToUse
	 * @throws StructureBuildingException
	 */
	private Atom determineFunctionalAtomToUse(String locant, BuildResults moleculeBuildResults) throws StructureBuildingException {
		for (int i = 0; i < moleculeBuildResults.functionalIDs.size(); i++) {
			Atom possibleAtom = moleculeBuildResults.getAtomById(moleculeBuildResults.functionalIDs.get(i));
			if (possibleAtom.hasLocant(locant)){
				moleculeBuildResults.functionalIDs.remove(i);
				return possibleAtom;
			}
		}
		//None of the functional atoms had an appropriate locant. Look for the special case where the locant is used to decide on the ester configuration c.f. O-methyl ..thioate and S-methyl ..thioate
		for (int i = 0; i < moleculeBuildResults.functionalIDs.size(); i++) {
			Atom possibleAtom = moleculeBuildResults.getAtomById(moleculeBuildResults.functionalIDs.get(i));
			if (possibleAtom.getNote("ambiguousElementAssignment")!=null){
				String[] atomIDs =possibleAtom.getNote("ambiguousElementAssignment").split(",");
				for (int j = 0; j < atomIDs.length; j++) {
					Atom a =moleculeBuildResults.getAtomById(Integer.parseInt(atomIDs[j]));
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
		List<Fragment> fragments = state.fragManager.getFragPile();
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
