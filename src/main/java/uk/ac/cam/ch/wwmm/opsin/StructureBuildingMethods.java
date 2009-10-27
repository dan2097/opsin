package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Stack;
import java.util.regex.Pattern;

import uk.ac.cam.ch.wwmm.ptclib.string.StringTools;
import uk.ac.cam.ch.wwmm.ptclib.xml.XOMTools;

import nu.xom.Element;
import nu.xom.Elements;

class StructureBuildingMethods {
	private static Pattern matchComma =Pattern.compile(",");
	private StructureBuildingMethods() {}

	/**
	 * Resolves a word/bracket:
	 * Locanted attributes of words are resolved onto their group
	 * Locanted substitution is performed
	 * Connections involving multi radicals are processed
	 * Unlocanted attributes of words are resolved onto their group
	 * 
	 * @param state
	 * @param word
	 * @throws StructureBuildingException
	 */
	static void resolveWordOrBracket(BuildState state, Element word) throws StructureBuildingException {
		if (!word.getLocalName().equals("word") && !word.getLocalName().equals("bracket")){
			throw new StructureBuildingException("A word or bracket is the expected input");
		}
		recursivelyResolveLocantedFeatures(state, word);
		recursivelyResolveUnLocantedFeatures(state, word);
	}
	
	/**
	 * Performs locanted attribute resolution
	 * then additive joining of fragments
	 * then locanted substitutive joining of fragments
	 * 
	 * @param state
	 * @param word
	 * @throws StructureBuildingException
	 */
	static void recursivelyResolveLocantedFeatures(BuildState state, Element word) throws StructureBuildingException {
		if (!word.getLocalName().equals("word") && !word.getLocalName().equals("bracket")){
			throw new StructureBuildingException("A word or bracket is the expected input");
		}
		List<Element> subsBracketsAndRoots = OpsinTools.findChildElementsWithTagNames(word, new String[]{"bracket","substituent","root"});
		//substitution occurs left to right so by doing this right to left you ensure that any groups that will come into existence
		//due to multipliers being expanded will be in existence
		for (int i =subsBracketsAndRoots.size()-1; i>=0; i--) {
			Element subBracketOrRoot = subsBracketsAndRoots.get(i);
			if (subBracketOrRoot.getLocalName().equals("bracket")){
				recursivelyResolveLocantedFeatures(state,subBracketOrRoot);
				if (potentiallyCanSubstitute(subBracketOrRoot)){
					performAdditiveOperations(state, subBracketOrRoot);
					performLocantedSubstitutiveOperations(state, subBracketOrRoot);
				}
			}
			else{
				resolveRootOrSubstituentLocanted(state, subBracketOrRoot);
			}
		}
	}
	
	/**
	 * Performs locanted attribute resolution
	 * then additive joining of fragments
	 * then locanted substitutive joining of fragments
	 * 
	 * @param state
	 * @param word
	 * @throws StructureBuildingException
	 */
	static void recursivelyResolveUnLocantedFeatures(BuildState state, Element word) throws StructureBuildingException {
		if (!word.getLocalName().equals("word") && !word.getLocalName().equals("bracket")){
			throw new StructureBuildingException("A word or bracket is the expected input");
		}
		List<Element> subsBracketsAndRoots = OpsinTools.findChildElementsWithTagNames(word, new String[]{"bracket","substituent","root"});
		//substitution occurs left to right so by doing this right to left you ensure that any groups that will come into existence
		//due to multipliers being expanded will be in existence
		for (int i =subsBracketsAndRoots.size()-1; i>=0; i--) {
			Element subBracketOrRoot = subsBracketsAndRoots.get(i);
			if (subBracketOrRoot.getLocalName().equals("bracket")){
				recursivelyResolveUnLocantedFeatures(state,subBracketOrRoot);
				if (potentiallyCanSubstitute(subBracketOrRoot)){
					performUnLocantedSubstitutiveOperations(state, subBracketOrRoot);
				}
			}
			else{
				resolveRootOrSubstituentUnLocanted(state, subBracketOrRoot);
			}
		}
	}
	
	static void resolveRootOrSubstituentLocanted(BuildState state, Element subOrRoot) throws StructureBuildingException {
		
		resolveLocantedFeatures(state, subOrRoot);//e.g. unsaturators, hydro groups and heteroatom replacement

		boolean foundSomethingToSubstitute = potentiallyCanSubstitute(subOrRoot);

		if (foundSomethingToSubstitute){
			performAdditiveOperations(state, subOrRoot);//e.g. ethylenediimino, oxyethylene (operations where two outIDs are used to produce the bond and no locant is required as groups)
			performLocantedSubstitutiveOperations(state, subOrRoot);//e.g. 2-methyltoluene
		}
	}
	
	static void resolveRootOrSubstituentUnLocanted(BuildState state, Element subOrRoot) throws StructureBuildingException {
		
		boolean foundSomethingToSubstitute = potentiallyCanSubstitute(subOrRoot);
		
		resolveUnLocantedFeatures(state, subOrRoot);//e.g. unsaturators, hydro groups and heteroatom replacement
		
		if (foundSomethingToSubstitute){
			performUnLocantedSubstitutiveOperations(state, subOrRoot);//e.g. tetramethylfuran
		}
	}
	
	
	private static void multiplyOutAndSubstitute(BuildState state, Element subOrBracket, String locantString) throws StructureBuildingException {
		int multiplier = Integer.parseInt(subOrBracket.getAttributeValue("multiplier"));
		String[] locants =null;
		if (locantString !=null){
			locants = matchComma.split(subOrBracket.getAttributeValue("locant"));
		}
		List<Element> clonedElements = new ArrayList<Element>();
		for (int i = multiplier -1; i >=0; i--) {
			Element currentElement;
			if (i!=0){
				currentElement = state.fragManager.cloneElement(subOrBracket, state, StringTools.multiplyString("'", i));
				clonedElements.add(currentElement);
			}
			else{
				currentElement=subOrBracket;
			}
			Element rootSubOrRoot;
			if (currentElement.getLocalName().equals("bracket")){
				List<Element> subsBracketsAndRoots = OpsinTools.findChildElementsWithTagNames(currentElement, new String[]{"bracket","substituent","root"});
				while (subsBracketsAndRoots.get(subsBracketsAndRoots.size()-1).getLocalName().equals("bracket")){
					subsBracketsAndRoots = OpsinTools.findChildElementsWithTagNames(subsBracketsAndRoots.get(subsBracketsAndRoots.size()-1), new String[]{"bracket","substituent","root"});
				}
				rootSubOrRoot =subsBracketsAndRoots.get(subsBracketsAndRoots.size()-1);
			}
			else{
				rootSubOrRoot = currentElement;
			}
			Fragment frag = state.xmlFragmentMap.get(rootSubOrRoot.getFirstChildElement("group"));
			if (frag.getOutIDs().size() !=1 ){
				throw new StructureBuildingException("Substitutive bond formation failure: Fragment expected to have one OutId but had: "+ frag.getOutIDs().size());
			}
			Atom parentAtom =null;
			if (locants !=null){
				String locant =locants[i];
				Fragment parentFrag = findFragmentWithLocant(state, subOrBracket, locant);
				if (parentFrag==null){
					throw new StructureBuildingException("Cannot find in scope atom with locant: " +locant + " to substitute onto");
				}
				parentAtom = parentFrag.getAtomByLocantOrThrow(locant);
			}
			else{
				parentAtom = findAtomForSubstitution(state, subOrBracket, frag.getOutID(0).valency);
				if (parentAtom==null){
					throw new StructureBuildingException("Cannot find in scope atom to substitute onto");
				}
			}
			
			joinFragmentsSubstitutively(state, frag, parentAtom);

		}
		for (Element clone : clonedElements) {//make sure cloned substituents don't substitute onto each other!
			XOMTools.insertAfter(subOrBracket, clone);
		}
	}

	private static void performLocantedSubstitutiveOperations(BuildState state, Element subBracketOrRoot) throws StructureBuildingException {
		Element group;
		if (subBracketOrRoot.getLocalName().equals("bracket")){
			group =findRightMostGroupInBracket(subBracketOrRoot);
		}
		else{
			group =subBracketOrRoot.getFirstChildElement("group");
		}
		Fragment frag = state.xmlFragmentMap.get(group);
		if (frag.getOutIDs().size() >=1 && subBracketOrRoot.getAttribute("locant")!=null){
			String locantString = subBracketOrRoot.getAttributeValue("locant");
			if (subBracketOrRoot.getAttribute("multiplier")!=null){//e.g. 1,2-diethyl
				multiplyOutAndSubstitute(state, subBracketOrRoot, locantString);
			}
			else{
				Fragment parentFrag = findFragmentWithLocant(state, subBracketOrRoot, locantString);
				joinFragmentsSubstitutively(state, frag, parentFrag.getAtomByLocantOrThrow(locantString));
			}
		}
	}

	private static void performUnLocantedSubstitutiveOperations(BuildState state, Element subBracketOrRoot) throws StructureBuildingException {
		Element group;
		if (subBracketOrRoot.getLocalName().equals("bracket")){
			group =findRightMostGroupInBracket(subBracketOrRoot);
		}
		else{
			group =subBracketOrRoot.getFirstChildElement("group");
		}
		Fragment frag = state.xmlFragmentMap.get(group);
		if (frag.getOutIDs().size() >=1){
			if (subBracketOrRoot.getAttribute("locant")!=null){
				throw new StructureBuildingException("Substituent has an unused outId and has a locant but locanted susbtitution should already been been performed!");
			}
			if (group.getAttribute("isAMultiRadical")!=null && frag.getOutIDs().size() ==1 ){//already been used for bonding once
				return;
			}
			if (subBracketOrRoot.getAttribute("multiplier")!=null){//e.g. 1,2-diethyl
				multiplyOutAndSubstitute(state, subBracketOrRoot, null);
			}
			else{
				Atom atomToJoinTo = findAtomForSubstitution(state, subBracketOrRoot, frag.getOutID(0).valency);
				if (atomToJoinTo ==null){
					throw new StructureBuildingException("Unlocanted substitution failed: unable to find suitable atom to bond atom with id: " + frag.getOutID(0).id + " to!");
				}
				joinFragmentsSubstitutively(state, frag, atomToJoinTo);
			}
		}
	}

	/**
	 * Adds locanted unsaturators, heteroatoms and hydrogen elements to the group within the sub or root
	 * @param state
	 * @param subOrRoot
	 * @throws StructureBuildingException 
	 */
	static void resolveLocantedFeatures(BuildState state, Element subOrRoot) throws StructureBuildingException {
		Elements groups = subOrRoot.getChildElements("group");
		if (groups.size()!=1){
			throw new StructureBuildingException("Each sub or root should only have one group element. This indicates a bug in OPSIN");
		}
		Element group = subOrRoot.getFirstChildElement("group");
		Fragment thisFrag = state.xmlFragmentMap.get(group);
	
		ArrayList<Element> unsaturators = new ArrayList<Element>();
		ArrayList<Element> heteroatoms = new ArrayList<Element>();
		ArrayList<Element> hydrogenElements = new ArrayList<Element>();
	
		Elements children =subOrRoot.getChildElements();
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
		 * Add locanted functionality
		 */
		for(int i=hydrogenElements.size() -1;i >= 0;i--) {
			Element hydrogen = hydrogenElements.get(i);
			String locant = getLocant(hydrogen);
			if(!locant.equals("0")) {
				thisFrag.getAtomByLocantOrThrow(locant).subtractSpareValency(1);
				hydrogenElements.remove(hydrogen);
				hydrogen.detach();
			}
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
				unsaturator.detach();
			}
		}
	
		for(int i=heteroatoms.size() -1;i >= 0;i--) {
			Element heteroatom = heteroatoms.get(i);
			String locant = getLocant(heteroatom);
			String atomSMILES = heteroatom.getAttributeValue("value");
			if(!locant.equals("0")) {
				state.fragManager.makeHeteroatom(thisFrag.getAtomByLocantOrThrow(locant), atomSMILES, true);
				heteroatoms.remove(heteroatom);
				heteroatom.detach();
			}
		}
	}

	/**
	 * Adds locanted unsaturators, heteroatoms and hydrogen elements to the group within the sub or root
	 * @param state
	 * @param subOrRoot
	 * @throws StructureBuildingException 
	 */
	static void resolveUnLocantedFeatures(BuildState state, Element subOrRoot) throws StructureBuildingException {
		Elements groups = subOrRoot.getChildElements("group");
		if (groups.size()!=1){
			throw new StructureBuildingException("Each sub or root should only have one group element. This indicates a bug in OPSIN");
		}
		Element group = subOrRoot.getFirstChildElement("group");
		Fragment thisFrag = state.xmlFragmentMap.get(group);
	
		ArrayList<Element> unsaturators = new ArrayList<Element>();
		ArrayList<Element> heteroatoms = new ArrayList<Element>();
		ArrayList<Element> hydrogenElements = new ArrayList<Element>();
	
		Elements children =subOrRoot.getChildElements();
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
		
		if (hydrogenElements.size()>0){
			/*
			 * This function is not entirely straightforward as certain atoms definitely should have their spare valency reduced
			 * However names are not consistent as to whether they bother having the hydro tags do this!
			 * The atoms in atomsWithSV are in atom order those that can take a hydro element and then those that shouldn't really take a hydro element as its absence is unambiguous
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
	
		int idOfFirstAtomInFragment= thisFrag.getIdOfFirstAtom();
		//defaultId corresponds to an id on thisFrag; whenever it is used it is incremented
		//all non-suffix atoms in the fragment are eventually checked if the defaultId does not correspond to a suitable atom
		//e.g. if the id was 3 in a 5 atom fragment which had ids 1-5, atoms would be checked for suitability in the order 3,4,5,1,2
		int defaultId = idOfFirstAtomInFragment;
	
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
	}

	private static void performAdditiveOperations(BuildState state, Element subBracketOrRoot) throws StructureBuildingException {
		Element group;
		if (subBracketOrRoot.getLocalName().equals("bracket")){
			group =findRightMostGroupInBracket(subBracketOrRoot);
		}
		else{
			group =subBracketOrRoot.getFirstChildElement("group");
		}
		Fragment frag = state.xmlFragmentMap.get(group);
		if (frag.getOutIDs().size() >=1){
			if (group.getAttribute("isAMultiRadical")!=null){
				//this should be either additive nomenclature or multiplicative nomenclature
				Fragment nextFrag = getNextInScopeMultiValentFragment(subBracketOrRoot, state);
				if (nextFrag!=null){
					Element parentEl = (Element) state.xmlFragmentMap.getElement(nextFrag).getParent();
					if (parentEl.getFirstChildElement("multiplier")!=null){//multiplicative nomenclature
						//TODO resupport multiplicative nomenclature
					}
					else{
						joinFragmentsAdditively(state, frag, nextFrag);
					}
				}
			}
			else if (subBracketOrRoot.getAttribute("locant")==null){//e.g. chlorocarbonyl
				Fragment nextFrag = getNextInScopeFragment(state, subBracketOrRoot);
				if (nextFrag!=null && nextFrag.getOutIDs().size()>=1 && state.xmlFragmentMap.getElement(nextFrag).getAttribute("isAMultiRadical")!=null){
					if (state.xmlFragmentMap.getElement(nextFrag).getAttribute("substitutionRemovesOutIDs")!=null){//e.g. aminocarbonyl
						joinFragmentsAdditively(state, frag, nextFrag);
					}
					//something like aminomethylene is substitutive
				}
			}
		}
		else if (!subBracketOrRoot.getLocalName().equals("root")){//FIXME move this
			throw new StructureBuildingException("Everything other than the root of the compound is expected to form bonds of some sort");
		}
	}

	private static Fragment getNextInScopeFragment(BuildState state, Element subBracketOrRoot) throws StructureBuildingException {
		if (subBracketOrRoot.getParent()==null){
			throw new StructureBuildingException("substituent/bracket/root did not have a parent!");
		}
		Element parent =(Element) subBracketOrRoot.getParent();
		
		List<Element> children = OpsinTools.findChildElementsWithTagNames(parent, new String[]{"substituent", "bracket", "root"});//will be returned in index order
		int indexOfSubstituent =parent.indexOf(subBracketOrRoot);
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
			Element descendantChild =childDescendants.get(0);
			Element group = descendantChild.getFirstChildElement("group");
			if (group == null){
				throw new StructureBuildingException("substituent/root is missing its group!");
			}
			Fragment possibleFrag = state.xmlFragmentMap.get(group);
			return possibleFrag;
		}
		return null;
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
	private static Fragment getNextInScopeMultiValentFragment(Element substituentOrBracket, BuildState state) throws StructureBuildingException {
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
				if (possibleFrag.getOutIDs().size() >=1 && group.getAttribute("isAMultiRadical")!=null ||(possibleFrag.getOutIDs().size()==0 && possibleFrag.getInIDs().size()==1)){
					return possibleFrag;
				}
			}
		}
		return null;
	}

	private static void joinFragmentsAdditively(BuildState state, Fragment fragToBeJoined, Fragment parentFrag) throws StructureBuildingException {
		int outIdCountOnParent = parentFrag.getOutIDs().size();
		if (outIdCountOnParent ==0){
			throw new StructureBuildingException("Additive bond formation failure: Fragment expected to have at least one OutId but had none");
		}
		OutID in = parentFrag.getOutID(0);
		Atom to = parentFrag.getAtomByIDOrThrow(in.id);
		int bondOrder = in.valency;
		if (in.setExplicitly != true){//not set explicitly so may be an inappropriate atom
			to = to.getFrag().getAtomByIdOrNextSuitableAtomOrThrow(to.getID(), bondOrder);
		}
		parentFrag.removeOutID(in);
		
		
		int outIdCount = fragToBeJoined.getOutIDs().size();
		if (outIdCount ==0){
			throw new StructureBuildingException("Additive bond formation failure: Fragment expected to have at least one OutId but had none");
		}
		OutID out =null;

		for (int i =outIdCount -1; i>=0; i--) {
			if (fragToBeJoined.getOutID(i).valency == bondOrder){
				out = fragToBeJoined.getOutID(i);
				break;
			}
		}

		if (out ==null){
			if (outIdCount >=bondOrder){//handles cases like nitrilo needing to be -N= (remove later outIds first as per usual)
				int valency =0;
				for (int i =outIdCount -1; i>=0; i--) {
					OutID nextOutID = fragToBeJoined.getOutID(i);
					valency += nextOutID.valency;
					if (valency==bondOrder){
						nextOutID.valency=valency;
						out = nextOutID;
						break;
					}
					fragToBeJoined.removeOutID(nextOutID);
				}
				if (out==null){
					throw new StructureBuildingException("Additive bond formation failure: bond order disagreement");
				}
			}
			else{
				throw new StructureBuildingException("Additive bond formation failure: bond order disagreement");
			}
		}
		
		Atom from = fragToBeJoined.getAtomByIDOrThrow(out.id);
		if (out.setExplicitly != true){//not set explicitly so may be an inappropriate atom
			from=from.getFrag().getAtomByIdOrNextSuitableAtomOrThrow(from.getID(), bondOrder);
		}
		fragToBeJoined.removeOutID(out);
	
		state.fragManager.attachFragments(from, to, bondOrder);
		System.out.println("Additively bonded to: " + from.getID() + " " + to.getID());
	}

	private static void joinFragmentsSubstitutively(BuildState state, Fragment fragToBeJoined, Atom atomToJoinTo) throws StructureBuildingException {
		int outIdCount = fragToBeJoined.getOutIDs().size();
		if (outIdCount >1 && fragToBeJoined.getAtomList().size()==1 && state.xmlFragmentMap.getElement(fragToBeJoined).getAttribute("substitutionRemovesOutIDs")==null){
			//special case e.g. methylenecyclohexane
			int valencyOfOutId =0;
			for (int i = outIdCount -1; i >=0 ; i--) {//remove all outIDs and add one with the total valency of all those that have been removed
				OutID out = fragToBeJoined.getOutID(i);
				valencyOfOutId +=out.valency;
				fragToBeJoined.removeOutID(i);
			}
			fragToBeJoined.addOutID(fragToBeJoined.getIdOfFirstAtom(), valencyOfOutId, true);
		}
		if (outIdCount ==0 ){
			throw new StructureBuildingException("Substitutive bond formation failure: Fragment expected to have one OutId but had none");
		}
		OutID out = fragToBeJoined.getOutID(0);
		Atom from = fragToBeJoined.getAtomByIDOrThrow(out.id);
		int bondOrder = out.valency;
		if (out.setExplicitly != true){//not set explicitly so may be an inappropriate atom
			from=from.getFrag().getAtomByIdOrNextSuitableAtomOrThrow(from.getID(), bondOrder);
		}
		fragToBeJoined.removeOutID(out);
		state.fragManager.attachFragments(from, atomToJoinTo, bondOrder);
		System.out.println("Substitutively bonded to: " + from.getID() + " " + atomToJoinTo.getID());
	}

	private static Atom findAtomForSubstitution(BuildState state, Element subOrBracket, int bondOrder) throws StructureBuildingException {
		//case where you should actually be substituting onto the previous element e.g. 5-(4-methylphenylcarbonyl)pentane
		Atom to =null;
		ArrayList<Fragment> possibleParents =findAlternativeFragments(state, subOrBracket);
		for (Fragment fragment : possibleParents) {
			to = fragment.getAtomByIdOrNextSuitableAtom(fragment.getDefaultInID(), bondOrder, true);
			if (to !=null){
				break;
			}
		}
		return to;
	}

	/**
	 * Finds all the groups accessible from the currentElement taking into account brackets
	 * i.e. those that it is feasible that the group of the currentElement could substitute onto
	 * @param state 
	 * @param currentElement
	 * @return A list of fragments in the order to try them as possible parent fragments
	 */
	private static ArrayList<Fragment> findAlternativeFragments(BuildState state, Element startingElement) {
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

	/**
	 * Checks through the groups accessible from the currentElement taking into account brackets
	 * i.e. those that it is feasible that the group of the currentElement could substitute onto
	 * @param state 
	 * @param currentElement
	 * @param locant: the locant string to check for the presence of
	 * @return The fragment with the locant, or null
	 * @throws StructureBuildingException 
	 */
	private static Fragment findFragmentWithLocant(BuildState state, Element startingElement, String locant) throws StructureBuildingException {
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

	private static Element findRightMostGroupInBracket(Element bracket) {
		List<Element> subsBracketsAndRoots = OpsinTools.findChildElementsWithTagNames(bracket, new String[]{"bracket","substituent","root"});
		while (subsBracketsAndRoots.get(subsBracketsAndRoots.size()-1).getLocalName().equals("bracket")){
			subsBracketsAndRoots = OpsinTools.findChildElementsWithTagNames(subsBracketsAndRoots.get(subsBracketsAndRoots.size()-1), new String[]{"bracket","substituent","root"});
		}
		return subsBracketsAndRoots.get(subsBracketsAndRoots.size()-1).getFirstChildElement("group");
	}
	
	private static boolean potentiallyCanSubstitute(Element subBracketOrRoot) {
		Element parent =(Element) subBracketOrRoot.getParent();
		Elements children =parent.getChildElements();
		for (int i = parent.indexOf(subBracketOrRoot) +1 ; i < children.size(); i++) {
			if (!children.get(i).getLocalName().equals("hyphen")){
				return true;
			}
		}
		return false;
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
}
