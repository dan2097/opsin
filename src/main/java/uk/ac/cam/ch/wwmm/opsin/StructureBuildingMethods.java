package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Stack;
import java.util.regex.Pattern;


import nu.xom.Attribute;
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
		//TODO check all things that can substitute have outIDs
		//TOOD think whether you can avoid the need to have a cansubstitute function by only using appropriate group
		List<Element> subsBracketsAndRoots = XOMTools.getDescendantElementsWithTagNames(word, new String[]{"bracket","substituent","root"});
		for (int i = 0; i < subsBracketsAndRoots.size(); i++) {
			Element subsBracketsAndRoot = subsBracketsAndRoots.get(i);
			if (subsBracketsAndRoot.getAttribute("multiplier")!=null){
				throw new StructureBuildingException("Structure building problem: multiplier on :" +subsBracketsAndRoot.getLocalName() + " was never used");
			}
		}
		List<Element> groups = XOMTools.getDescendantElementsWithTagName(word, "group");
		for (int i = 0; i < groups.size(); i++) {
			Element group = groups.get(i);
			if (group.getAttribute("resolved")==null && i!=groups.size()-1){
				throw new StructureBuildingException("Structure building problem: Bond was not made from :" +group.getValue() + " but one should of been");
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
	static void recursivelyResolveLocantedFeatures(BuildState state, Element word) throws StructureBuildingException {
		if (!word.getLocalName().equals("word") && !word.getLocalName().equals("bracket")){
			throw new StructureBuildingException("A word or bracket is the expected input");
		}
		List<Element> subsBracketsAndRoots = XOMTools.getChildElementsWithTagNames(word, new String[]{"bracket","substituent","root"});
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
		List<Element> subsBracketsAndRoots = XOMTools.getChildElementsWithTagNames(word, new String[]{"bracket","substituent","root"});
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
	
	
	private static void performLocantedSubstitutiveOperations(BuildState state, Element subBracketOrRoot) throws StructureBuildingException {
		Element group;
		if (subBracketOrRoot.getLocalName().equals("bracket")){
			group =findRightMostGroupInBracket(subBracketOrRoot);
		}
		else{
			group =subBracketOrRoot.getFirstChildElement("group");
		}
		if (group.getAttribute("resolved")!=null){
			return;
		}
		Fragment frag = state.xmlFragmentMap.get(group);
		if (frag.getOutIDs().size() >=1 && subBracketOrRoot.getAttribute("locant")!=null){
			String locantString = subBracketOrRoot.getAttributeValue("locant");
			if (frag.getOutIDs().size() >1){
				checkAndApplySpecialCaseWhereOutIdsCanBeCombinedOrThrow(state, frag);
			}
			if (subBracketOrRoot.getAttribute("multiplier")!=null){//e.g. 1,2-diethyl
				multiplyOutAndSubstitute(state, subBracketOrRoot);
			}
			else{
				Fragment parentFrag = findFragmentWithLocant(state, subBracketOrRoot, locantString);
				if (parentFrag==null){
					throw new StructureBuildingException("Cannot find in scope fragment with atom with locant " + locantString + ".");
				}
				group.addAttribute(new Attribute("resolved","yes"));
				Element groupToAttachTo = state.xmlFragmentMap.getElement(parentFrag);
				if (groupToAttachTo.getAttribute("acceptsAdditiveBonds")!=null && parentFrag.getOutIDs().size()>0 && groupToAttachTo.getAttribute("isAMultiRadical")!=null 
						&& parentFrag.getAtomByLocantOrThrow(locantString).getOutValency()>0 && frag.getOutID(0).valency==1){
					//horrible special case, probably not even technically allowed by the IUPAC rules e.g. C-methylcarbonimidoyl
					joinFragmentsAdditively(state, frag, parentFrag);
				}
				else{
					joinFragmentsSubstitutively(state, frag, parentFrag.getAtomByLocantOrThrow(locantString));
				}
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
		if (group.getAttribute("resolved")!=null){
			return;
		}
		Fragment frag = state.xmlFragmentMap.get(group);
		if (frag.getOutIDs().size() >=1){
			if (subBracketOrRoot.getAttribute("locant")!=null){
				throw new StructureBuildingException("Substituent has an unused outId and has a locant but locanted susbtitution should already been been performed!");
			}
			if (frag.getOutIDs().size() > 1){
				checkAndApplySpecialCaseWhereOutIdsCanBeCombinedOrThrow(state, frag);
			}
			if (subBracketOrRoot.getAttribute("multiplier")!=null){//e.g. diethyl
				multiplyOutAndSubstitute(state, subBracketOrRoot);
			}
			else{
				Atom atomToJoinTo = findAtomForSubstitution(state, subBracketOrRoot, frag.getOutID(0).valency);
				if (atomToJoinTo ==null){
					throw new StructureBuildingException("Unlocanted substitution failed: unable to find suitable atom to bond atom with id:" + frag.getOutID(0).id + " to!");
				}
				group.addAttribute(new Attribute("resolved","yes"));
				joinFragmentsSubstitutively(state, frag, atomToJoinTo);
			}
		}
	}

	/**
	 * Multiplies out groups/brakets and substitutes them. The attribute "locant" is checked for locants
	 * If it is present it should contain a comma seperated list of locants
	 * The strategy employed is to clone subOrBracket and its associated fragments as many times as the multiplier attribute
	 * perform(Un)LocantedSubstitutiveOperations is then called with on each call a different clone (or the original element) being in position
	 * Hence bonding between the clones is impossible
	 * @param state
	 * @param subOrBracket
	 * @throws StructureBuildingException
	 */
	private static void multiplyOutAndSubstitute(BuildState state, Element subOrBracket) throws StructureBuildingException {
		int multiplier = Integer.parseInt(subOrBracket.getAttributeValue("multiplier"));
		subOrBracket.removeAttribute(subOrBracket.getAttribute("multiplier"));
		String[] locants =null;
		if (subOrBracket.getAttribute("locant") !=null){
			locants = matchComma.split(subOrBracket.getAttributeValue("locant"));
		}
		Element parentWordOrBracket =(Element) subOrBracket.getParent();
		int indexOfSubOrBracket = parentWordOrBracket.indexOf(subOrBracket);
		subOrBracket.detach();
		List<Element> multipliedElements = new ArrayList<Element>();
		for (int i = multiplier -1; i >=0; i--) {
			Element currentElement;
			if (i!=0){
				currentElement = state.fragManager.cloneElement(state, subOrBracket, StringTools.multiplyString("'", i));
				if (currentElement.getLocalName().equals("bracket")){
					addPrimesToLocantedStereochemistryElements(currentElement, StringTools.multiplyString("'", i));
				}
			}
			else{
				currentElement = subOrBracket;
			}
			multipliedElements.add(currentElement);
			parentWordOrBracket.insertChild(currentElement, indexOfSubOrBracket);
			if (locants !=null){
				currentElement.getAttribute("locant").setValue(locants[i]);
				performLocantedSubstitutiveOperations(state, currentElement);
			}
			else{
				performUnLocantedSubstitutiveOperations(state, currentElement);
			}
			currentElement.detach();
		}
		for (Element multipliedElement : multipliedElements) {//attach all the multiplied subs/brackets
			parentWordOrBracket.insertChild(multipliedElement, indexOfSubOrBracket);
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
			if(!locant.equals("0")) {
				String atomSMILES = heteroatom.getAttributeValue("value");
				state.fragManager.makeHeteroatom(thisFrag.getAtomByLocantOrThrow(locant), atomSMILES, true);
				if (heteroatom.getAttribute("lambda")!=null){
					thisFrag.getAtomByLocantOrThrow(locant).setValency(Integer.parseInt(heteroatom.getAttributeValue("lambda")));
				}
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
			if (heteroatom.getAttribute("lambda")!=null){
				atomToReplaceWithHeteroAtom.setValency(Integer.parseInt(heteroatom.getAttributeValue("lambda")));
			}
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
		if (subBracketOrRoot.getAttribute("locant")!=null){//additive nomenclature does not employ locants
			return;
		}
		Element group;
		if (subBracketOrRoot.getLocalName().equals("bracket")){
			group =findRightMostGroupInBracket(subBracketOrRoot);
		}
		else{
			group =subBracketOrRoot.getFirstChildElement("group");
		}
		if (group.getAttribute("resolved")!=null){
			return;
		}
		Fragment frag = state.xmlFragmentMap.get(group);
		int outIDCount = frag.getOutIDs().size();
		if (outIDCount >=1){
			if (subBracketOrRoot.getAttribute("multiplier") ==null){
				Element nextSiblingEl = (Element) XOMTools.getNextSibling(subBracketOrRoot);
				if (nextSiblingEl.getAttribute("multiplier")!=null && 
						(outIDCount >= Integer.parseInt(nextSiblingEl.getAttributeValue("multiplier")) || //probably multiplicative nomenclature, should be as many outIds as the multiplier
						outIDCount==1 && frag.getOutID(0).valency==2 && Integer.parseInt(nextSiblingEl.getAttributeValue("multiplier"))==2) &&
						hasRootLikeOrMultiRadicalGroup(state, nextSiblingEl)){
					if (outIDCount==1){//special case e.g. 4,4'-(benzylidene)dianiline
						OutID out = frag.getOutID(0);
						out.valency = 1;
						frag.addOutID(out.id, 1, out.setExplicitly);
						//special case where something like benzylidene is being used as if it meant benzdiyl for multiplicative nomenclature
						//this is allowed in the IUPAC 79 recommendations but not recommended in the current recommendations
					}
					performMultiplicativeOperations(state, group, nextSiblingEl);
				}
				else if (group.getAttribute("isAMultiRadical")!=null){//additive nomenclature e.g. ethyleneoxy
					Fragment nextFrag = getNextInScopeMultiValentFragment(state, subBracketOrRoot);
					if (nextFrag!=null){
						Element nextMultiRadicalGroup = state.xmlFragmentMap.getElement(nextFrag);
						Element parentSubOrRoot = (Element) nextMultiRadicalGroup.getParent();
						if (!state.wordRule.equals("polymer")){//imino does not behave like a substituent in polymers only as a linker
							if (nextMultiRadicalGroup.getAttribute("iminoLike")!=null){//imino/methylene can just act as normal substituents, should an additive bond really be made???
								List<Fragment> alternativeFragments = findAlternativeFragments(state, subBracketOrRoot);
								if (nextFrag !=alternativeFragments.get(alternativeFragments.size()-1)){//imino is not the absolute next frag
									return;
								}
							}
							if (group.getAttribute("iminoLike")!=null &&  (XOMTools.getNextSibling(group.getParent())==null || XOMTools.getNextSibling(subBracketOrRoot)!=parentSubOrRoot)){
								return;//they are at different levels or again not adjacent e.g. ((chloroimino)ethylene)dibenzene
							}
						}
						if (parentSubOrRoot.getAttribute("multiplier")!=null){
							throw new StructureBuildingException("Attempted to form additive bond to a multiplied component");
						}
						group.addAttribute(new Attribute("resolved","yes"));
						joinFragmentsAdditively(state, frag, nextFrag);
					}
				}
				else {//e.g. chlorocarbonyl or hydroxy(sulfanyl)phosphoryl
					List<Fragment> siblingFragments = findAlternativeFragments(state, subBracketOrRoot);
					if (siblingFragments.size()>0){
						Fragment nextFrag = siblingFragments.get(siblingFragments.size()-1);
						Element nextGroup = state.xmlFragmentMap.getElement(nextFrag);
						if (nextGroup.getAttribute("acceptsAdditiveBonds")!=null && nextFrag!=null && nextGroup.getAttribute("isAMultiRadical")!=null  && (nextFrag.getOutIDs().size()>1|| nextGroup.getAttribute("resolved")!=null && nextFrag.getOutIDs().size()>=1 )){
							Atom toAtom = nextFrag.getAtomByIDOrThrow(nextFrag.getOutID(0).id);
							if (calculateSubstitutableHydrogenAtoms(toAtom) ==0){
								group.addAttribute(new Attribute("resolved","yes"));
								joinFragmentsAdditively(state, frag, nextFrag);//e.g. aminocarbonyl or aminothio
							}
						}
						if (group.getAttribute("resolved")==null && siblingFragments.size()>1){
							for (int i = 0; i< siblingFragments.size()-1; i++) {
								Fragment lastFrag = siblingFragments.get(i);
								Element lastGroup = state.xmlFragmentMap.getElement(lastFrag);
								if (lastGroup.getAttribute("acceptsAdditiveBonds")!=null && lastFrag!=null && lastGroup.getAttribute("isAMultiRadical")!=null  && (lastFrag.getOutIDs().size()>1|| lastGroup.getAttribute("resolved")!=null && lastFrag.getOutIDs().size()>=1 )){
									Atom toAtom = lastFrag.getAtomByIDOrThrow(lastFrag.getOutID(0).id);
									if (calculateSubstitutableHydrogenAtoms(toAtom) ==0){
										group.addAttribute(new Attribute("resolved","yes"));
										joinFragmentsAdditively(state, frag, lastFrag);//e.g. hydroxy(sulfanyl)phosphoryl
									}
									break;
								}
								
								//loop may continue if lastFrag was in fact completely unsubtitutable e.g. hydroxy...phosphoryloxy. The oxy is unsubstituable as the phosphoryl will already have bonded to it
								if (lastFrag.getAtomByIdOrNextSuitableAtom(lastFrag.getDefaultInID(), frag.getOutID(outIDCount-1).valency, true)!=null){
									break;
								}
							}
						}
					}
				}
			}
			else{// e.g. dimethoxyphosphoryl or bis(methylamino)phosphoryl
				List<Fragment> siblingFragments = findAlternativeFragments(state, subBracketOrRoot);	
				if (siblingFragments.size()>0){
					int multiplier = Integer.parseInt(subBracketOrRoot.getAttributeValue("multiplier"));
					Fragment nextFrag = siblingFragments.get(siblingFragments.size()-1);
					Element nextGroup = state.xmlFragmentMap.getElement(nextFrag);
					if (nextGroup.getAttribute("acceptsAdditiveBonds")!=null && nextFrag!=null && nextGroup.getAttribute("isAMultiRadical")!=null  && (nextFrag.getOutIDs().size()>=multiplier|| nextGroup.getAttribute("resolved")!=null && nextFrag.getOutIDs().size()>=multiplier +1 )){
						Atom toAtom = nextFrag.getAtomByIDOrThrow(nextFrag.getOutID(0).id);
						if (calculateSubstitutableHydrogenAtoms(toAtom) ==0){
							group.addAttribute(new Attribute("resolved","yes"));
							multiplyOutAndAdditivelyBond(state, subBracketOrRoot, nextFrag);//e.g.dihydroxyphosphoryl
						}
					}
					if (group.getAttribute("resolved")==null && siblingFragments.size()>1){
						for (int i = 0; i< siblingFragments.size()-1; i++) {
							Fragment lastFrag = siblingFragments.get(i);
							Element lastGroup = state.xmlFragmentMap.getElement(lastFrag);
							if (lastGroup.getAttribute("acceptsAdditiveBonds")!=null && lastFrag!=null && lastGroup.getAttribute("isAMultiRadical")!=null  &&  (lastFrag.getOutIDs().size()>=multiplier|| lastGroup.getAttribute("resolved")!=null && lastFrag.getOutIDs().size()>=multiplier +1 )){
								Atom toAtom = lastFrag.getAtomByIDOrThrow(lastFrag.getOutID(0).id);
								if (calculateSubstitutableHydrogenAtoms(toAtom) ==0){
									group.addAttribute(new Attribute("resolved","yes"));
									multiplyOutAndAdditivelyBond(state, subBracketOrRoot, lastFrag);//e.g. dihydroxyphosphoryloxy
								}
								break;
							}
							
							//loop may continue if lastFrag was in fact completely unsubtitutable e.g. hydroxy...phosphoryloxy. The oxy is unsubstituable as the phosphoryl will already have bonded to it
							if (lastFrag.getAtomByIdOrNextSuitableAtom(lastFrag.getDefaultInID(), frag.getOutID(outIDCount-1).valency, true)!=null){
								break;
							}
						}
					}
				}
			}
		}
	}

	/**
	 * Searches the input for something that either is a multiRadical or has no outIDs i.e. not dimethyl
	 * @param state 
	 * @param subBracketOrRoot
	 * @return
	 */
	private static boolean hasRootLikeOrMultiRadicalGroup(BuildState state, Element subBracketOrRoot) {
		List<Element> groups = XOMTools.getDescendantElementsWithTagName(subBracketOrRoot, "group");
		for (Element group : groups) {
			Fragment frag =state.xmlFragmentMap.get(group);
			int outIdCount =frag.getOutIDs().size();
			if (group.getAttribute("isAMultiRadical")!=null){
				if (outIdCount >=1 ){
					return true;//a multi radical
				}
			}
			else if (outIdCount ==0 && group.getAttribute("resolved")==null){
				return true;// a terminus
			}
		}
		return false;
	}

	/**
	 * Multiply out subOrBracket and additively bond all substituents to the specified fragment
	 * @param state
	 * @param subOrBracket
	 * @param fragToAdditivelyBondTo
	 * @throws StructureBuildingException
	 */
	private static void multiplyOutAndAdditivelyBond(BuildState state, Element subOrBracket, Fragment fragToAdditivelyBondTo) throws StructureBuildingException {
		int multiplier = Integer.parseInt(subOrBracket.getAttributeValue("multiplier"));
		subOrBracket.removeAttribute(subOrBracket.getAttribute("multiplier"));
		List<Element> clonedElements = new ArrayList<Element>();
		for (int i = multiplier -1; i >=0; i--) {
			Element currentElement;
			if (i!=0){
				currentElement = state.fragManager.cloneElement(state, subOrBracket, StringTools.multiplyString("'", i));
				if (currentElement.getLocalName().equals("bracket")){
					addPrimesToLocantedStereochemistryElements(currentElement, StringTools.multiplyString("'", i));
				}
				clonedElements.add(currentElement);
			}
			else{
				currentElement = subOrBracket;
			}
			Element group;
			if (currentElement.getLocalName().equals("bracket")){
				group = findRightMostGroupInBracket(currentElement);
			}
			else{
				group = currentElement.getFirstChildElement("group");
			}
			Fragment frag = state.xmlFragmentMap.get(group);
			if (frag.getOutIDs().size() !=1 ){
				throw new StructureBuildingException("Additive bond formation failure: Fragment expected to have one OutId in this case but had: "+ frag.getOutIDs().size());
			}
			joinFragmentsAdditively(state, frag, fragToAdditivelyBondTo);
		}
		for (Element clone : clonedElements) {//make sure cloned substituents don't substitute onto each other!
			XOMTools.insertAfter(subOrBracket, clone);
		}
	}

	/**
	 * Creates a build results from the input group for use as the input to the real performMultiplicativeOperations function
	 * @param state
	 * @param group
	 * @param multipliedParent
	 * @throws StructureBuildingException 
	 */
	private static void performMultiplicativeOperations(BuildState state, Element group, Element multipliedParent) throws StructureBuildingException{
		BuildResults multiRadicalBR = new BuildResults(state, (Element) group.getParent());
		performMultiplicativeOperations(state, multiRadicalBR, multipliedParent);
	}

	private static void performMultiplicativeOperations(BuildState state, BuildResults multiRadicalBR, Element multipliedParent) throws StructureBuildingException {
		int multiplier = Integer.parseInt(multipliedParent.getAttributeValue("multiplier"));
		if (multiplier > multiRadicalBR.getOutIDCount()){
			throw new StructureBuildingException("Multiplication bond formation failure: number of outIDs disagree with multiplier(multiplier: " + multiplier + ", outIDcount: " + multiRadicalBR.getOutIDCount()+ ") , this is an OPSIN bug");
		}
		if (state.debug){System.out.println(multiplier +" multiplicative bonds to be formed");}
		multipliedParent.removeAttribute(multipliedParent.getAttribute("multiplier"));
		List<String> inLocants = null;
		if (multipliedParent.getAttribute("inLocants")!=null){//true for the root of a multiplicative name
			String inLocantsString = multipliedParent.getAttributeValue("inLocants");
			if (inLocantsString.equals("default")){
				inLocants = new ArrayList<String>(multiplier);
				for (int i = 0; i < multiplier; i++) {
					inLocants.add("default");
				}
			}
			else{
				inLocants = StringTools.arrayToList(matchComma.split(inLocantsString));
				if (inLocants.size() != multiplier){
					throw new StructureBuildingException("Mismatch between multiplier and number inLocants in multiplicative nomenclature");
				}
			}
		}
		List<Element> clonedElements = new ArrayList<Element>();
		BuildResults newBr = new BuildResults();
		for (int i = multiplier -1; i >=0; i--) {
			Element currentElement;
			if (i!=0){
				currentElement = state.fragManager.cloneElement(state, multipliedParent, StringTools.multiplyString("'", i));
				clonedElements.add(currentElement);
			}
			else{
				currentElement=multipliedParent;
			}
			Element group;
			if (currentElement.getLocalName().equals("bracket")){
				group =getFirstMultiValentGroup(state, currentElement);
				if (group == null){//root will not have a multivalent group
					group = findRightMostGroupInBracket(currentElement);
				}
			}
			else{
				group = currentElement.getFirstChildElement("group");
			}
			Fragment frag = state.xmlFragmentMap.get(group);
			if (inLocants !=null){
				if (group.getAttribute("isAMultiRadical")!=null){//e.g. methylenedisulfonyl dichloride
					if (!multipliedParent.getAttributeValue("inLocants").equals("default")){
						throw new StructureBuildingException("inLocants should not be specified for a multiradical parent in multiplicative nomenclature");
					}
					Element rightMostGroup;
					if (currentElement.getLocalName().equals("bracket")){
						rightMostGroup = findRightMostGroupInBracket(currentElement);//this is the only one that can accept inIDs
					}
					else{
						rightMostGroup = currentElement.getFirstChildElement("group");
					}
					rightMostGroup.addAttribute(new Attribute("resolved","yes"));
				}
				else{
					group.addAttribute(new Attribute("resolved","yes"));
					boolean inIdAdded =false;
					for (int j = inLocants.size() -1; j >=0; j--) {
						String locant = inLocants.get(j);
						if (locant.equals("default")){//note that if one entry in inLocantArray is default then they all are "default"
							frag.addInID(frag.getDefaultInID(), 1);
							inIdAdded=true;
							inLocants.remove(j);
							break;
						}
						else{
							Atom inAtom = frag.getAtomByLocant(locant);
							if (inAtom!=null){
								frag.addInID(inAtom.getID(), 1);
								inIdAdded=true;
								inLocants.remove(j);
								break;
							}
						}
					}
					if (!inIdAdded){
						throw new StructureBuildingException("Locants for inIDs on the root were either misassigned to the root or were invalid: " + inLocants.toString() +" could not be assigned!");
					}
				}
			}
			if (frag.getInIDs().size()!=1 && frag.getOutIDs().size() ==0 ){
				throw new StructureBuildingException("Multiplication bond formation failure: OPSIN bug, input to joinFragmentsMultiplicatively was unexpected");
			}
			
			Element multiRadicalGroup =state.xmlFragmentMap.getElement(multiRadicalBR.getOutID(i).frag);
			if (multiRadicalGroup.getAttribute("resolved")==null){
				resolveUnLocantedFeatures(state, (Element) multiRadicalGroup.getParent());//the addition of unlocanted unsaturators can effect the position of radicals e.g. diazenyl
				multiRadicalGroup.addAttribute(new Attribute("resolved","yes"));
			}
			joinFragmentsAdditively(state, multiRadicalBR.getOutID(i).frag, frag);
			if (currentElement.getLocalName().equals("bracket")){
				recursivelyResolveUnLocantedFeatures(state, currentElement);//there may be outIDs that are involved in unlocanted substitution, these can be safely used now e.g. ...bis((3-hydroxy-4-methoxyphenyl)methylene) where (3-hydroxy-4-methoxyphenyl)methylene is the currentElement
			}
			
			if (inLocants ==null){
				//currentElement is not a root element. Need to build up a new BuildResults so as to call performMultiplicativeOperations again
				//at this stage an outID has been removed from the fragment within currentElement through an additive bond
				newBr.mergeBuildResults(new BuildResults(state, currentElement));
			}
		}

		if (newBr.getFragmentCount()==1){
			throw new StructureBuildingException("Multiplicative nomenclarture cannot yield only one temporary terminal fragment");
		}
		if (newBr.getFragmentCount()>=2){
			List<Element> siblings = XOMTools.getNextSiblingsOfTypes(multipliedParent, new String[]{"substituent", "bracket", "root"});
			if (siblings.size()==0){
				Element parentOfMultipliedEl = (Element) multipliedParent.getParent();
				if (parentOfMultipliedEl.getLocalName().equals("bracket")){//brackets are allowed
					siblings = XOMTools.getNextSiblingsOfTypes(parentOfMultipliedEl, new String[]{"substituent", "bracket", "root"});
					if (siblings.get(0).getAttribute("multiplier")==null){
						throw new StructureBuildingException("Multiplier not found where multiplier was expected for succesful multiplicative nomenclature");
					}
					performMultiplicativeOperations(state, newBr, siblings.get(0));
				}
				else{
					throw new StructureBuildingException("Could not find suitable element to continue multiplicative nomenclature");
				}
			}
			else{
				if (siblings.get(0).getAttribute("multiplier")==null){
					throw new StructureBuildingException("Multiplier not found where multiplier was expected for succesful multiplicative nomenclature");
				}
				performMultiplicativeOperations(state, newBr, siblings.get(0));
			}
		}
		
		for (Element clone : clonedElements) {//make sure cloned substituents don't substitute onto each other!
			XOMTools.insertAfter(multipliedParent, clone);
		}
	}

	/**
	 * Given a subsituent/bracket finds the next multi valent substituent/root that is in scope and hence its group
	 * e.g. for oxy(dichloromethyl)methylene given oxy substituent the methylene group would be found
	 * for oxy(dichloroethylene) given oxy substituent the ethylene group would be found
	 * for oxy(carbonylimino) given oxy carbonyl would be found
	 * @param state 
	 * @param substituentOrBracket
	 * @return frag
	 * @throws StructureBuildingException 
	 */
	private static Fragment getNextInScopeMultiValentFragment(BuildState state, Element substituentOrBracket) throws StructureBuildingException {
		if (!substituentOrBracket.getLocalName().equals("substituent") && !substituentOrBracket.getLocalName().equals("bracket")){
			throw new StructureBuildingException("Input to this function should be a substituent or bracket");
		}
		if (substituentOrBracket.getParent()==null){
			throw new StructureBuildingException("substituent did not have a parent!");
		}
		Element parent =(Element) substituentOrBracket.getParent();
		
		List<Element> children = XOMTools.getChildElementsWithTagNames(parent, new String[]{"substituent", "bracket", "root"});//will be returned in index order
		int indexOfSubstituent =parent.indexOf(substituentOrBracket);
		for (Element child : children) {
			if (parent.indexOf(child) <=indexOfSubstituent){//only want things after the input
				continue;
			}
			if (child.getAttribute("multiplier")!=null){
				continue;
			}
			List<Element> childDescendants;
			if (child.getLocalName().equals("bracket")){
				childDescendants = XOMTools.getDescendantElementsWithTagNames(child, new String[]{"substituent", "root"});//will be returned in depth-first order
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
				if (group.getAttribute("isAMultiRadical")!=null &&
						(possibleFrag.getOutIDs().size() >=2 || (possibleFrag.getOutIDs().size() >=1 && group.getAttribute("resolved")!=null ))){
					return possibleFrag;
				}
			}
		}
		return null;
	}
	
	/**
	 * Given a bracket searches in a depth first manner for the first multi valent group
	 * @param state 
	 * @param bracket
	 * @return group
	 * @throws StructureBuildingException 
	 */
	private static Element getFirstMultiValentGroup(BuildState state, Element bracket) throws StructureBuildingException {
		if (!bracket.getLocalName().equals("bracket")){
			throw new StructureBuildingException("Input to this function should be a bracket");
		}
		
		List<Element> groups = XOMTools.getDescendantElementsWithTagName(bracket, "group");//will be returned in index order
		for (Element group : groups) {
			Fragment possibleFrag = state.xmlFragmentMap.get(group);
			if (group.getAttribute("isAMultiRadical")!=null &&
					(possibleFrag.getOutIDs().size() >=2 || (possibleFrag.getOutIDs().size() >=1 && group.getAttribute("resolved")!=null ))){
				return group;
			}
		}
		return null;
	}

	private static void joinFragmentsAdditively(BuildState state, Fragment fragToBeJoined, Fragment parentFrag) throws StructureBuildingException {
		int outIdCount = fragToBeJoined.getOutIDs().size();
		if (outIdCount ==0){
			throw new StructureBuildingException("Additive bond formation failure: Fragment expected to have at least one OutId but had none");
		}
		
		Atom to;
		int bondOrder;
		if (parentFrag.getInIDs().size()==1){//special case for the parent of multiplicative nomenclature. This is really substitutive nomenclature
			InID in = parentFrag.getInID(0);
			to = parentFrag.getAtomByIDOrThrow(in.id);
			parentFrag.removeInID(in);
			bondOrder = fragToBeJoined.getOutID(outIdCount-1).valency;
		}
		else{
			int outIdCountOnParent = parentFrag.getOutIDs().size();
			if (outIdCountOnParent ==0){
				throw new StructureBuildingException("Additive bond formation failure: Fragment expected to have at least one OutId but had none");
			}
			OutID in = parentFrag.getOutID(0);
			to = parentFrag.getAtomByIDOrThrow(in.id);
			bondOrder = in.valency;
			if (in.setExplicitly != true){//not set explicitly so may be an inappropriate atom
				to = to.getFrag().getAtomByIdOrNextSuitableAtomOrThrow(to.getID(), bondOrder);
			}
			parentFrag.removeOutID(in);
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
				for (int i =outIdCount -1; i >= 0; i--) {
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
		if (state.debug){System.out.println("Additively bonded " + from.getID() + " (" +state.xmlFragmentMap.getElement(from.getFrag()).getValue()+") " + to.getID() + " (" +state.xmlFragmentMap.getElement(to.getFrag()).getValue()+")" );}
	}

	private static void joinFragmentsSubstitutively(BuildState state, Fragment fragToBeJoined, Atom atomToJoinTo) throws StructureBuildingException {
		int outIdCount = fragToBeJoined.getOutIDs().size();
		if (outIdCount >1){
			throw new StructureBuildingException("Substitutive bond formation failure: Fragment expected to have one OutId but had: "+ outIdCount);
		}
		if (outIdCount ==0 ){
			throw new StructureBuildingException("Substitutive bond formation failure: Fragment expected to have one OutId but had none");
		}
		if (state.xmlFragmentMap.getElement(fragToBeJoined).getAttribute("iminoLike")!=null){//special case for methylene/imino
			if (fragToBeJoined.getOutIDs().size()==1 && fragToBeJoined.getOutID(0).valency==1 ){
				fragToBeJoined.getOutID(0).valency=2;
			}
		}
		OutID out = fragToBeJoined.getOutID(0);
		Atom from = fragToBeJoined.getAtomByIDOrThrow(out.id);
		int bondOrder = out.valency;
		if (out.setExplicitly != true){//not set explicitly so may be an inappropriate atom
			from=from.getFrag().getAtomByIdOrNextSuitableAtomOrThrow(from.getID(), bondOrder);
		}
		fragToBeJoined.removeOutID(out);
		state.fragManager.attachFragments(from, atomToJoinTo, bondOrder);
		if (state.debug){System.out.println("Substitutively bonded " + from.getID() + " (" +state.xmlFragmentMap.getElement(from.getFrag()).getValue()+") " + atomToJoinTo.getID() + " (" +state.xmlFragmentMap.getElement(atomToJoinTo.getFrag()).getValue()+")");}
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
	 * @return A list of fragments in the order to try them as possible parent fragments (for substitutive operations)
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
			List<Element> siblings = XOMTools.getChildElementsWithTagNames(parent, new String[]{"bracket", "substituent", "root"});
	
			for (Element bracketOrSub : siblings) {
				if (bracketOrSub.getAttribute("multiplier")!=null){
					continue;
				}
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
			List<Element> siblings = XOMTools.getChildElementsWithTagNames(parent, new String[]{"bracket", "substituent", "root"});
	
			for (Element bracketOrSub : siblings) {
				if (bracketOrSub.getAttribute("multiplier")!=null){
					continue;
				}
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

	static Element findRightMostGroupInBracket(Element bracket) {
		List<Element> subsBracketsAndRoots = XOMTools.getChildElementsWithTagNames(bracket, new String[]{"bracket","substituent","root"});
		while (subsBracketsAndRoots.get(subsBracketsAndRoots.size()-1).getLocalName().equals("bracket")){
			subsBracketsAndRoots = XOMTools.getChildElementsWithTagNames(subsBracketsAndRoots.get(subsBracketsAndRoots.size()-1), new String[]{"bracket","substituent","root"});
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

	/**
	 * In cases such as methylenecyclohexane two outIDs are combined to form a single outID with valency
	 * equal to sum of the valency of the other outIDs.
	 * This is only allowed on substituents where all the outIDs are on the same atom
	 * @param state
	 * @param frag
	 * @throws StructureBuildingException 
	 */
	private static void checkAndApplySpecialCaseWhereOutIdsCanBeCombinedOrThrow(BuildState state, Fragment frag) throws StructureBuildingException {
		int outIdCount = frag.getOutIDs().size();
		if (outIdCount<=1){
			return;
		}
		//special case- all outIDs on same atom e.g. methylenecyclohexane
		int idOfFirstOutId = frag.getOutID(0).id;
		int valencyOfOutId =0;
		for (int i = outIdCount -1; i >=0 ; i--) {//remove all outIDs and add one with the total valency of all those that have been removed
			OutID out = frag.getOutID(i);
			if (out.id !=idOfFirstOutId){
				throw new StructureBuildingException("Substitutive bond formation failure: Fragment expected to have one OutId but had: "+ outIdCount);
			}
			valencyOfOutId +=out.valency;
			frag.removeOutID(i);
		}
		frag.addOutID(frag.getIdOfFirstAtom(), valencyOfOutId, true);
	}
	
	/**
	 * Calculates the number of substitutable hydrogen by taking into account:
	 * Specified valency if applicable, outIDs and the lowest valency state that will satisfy these 
	 * e.g. thio has 2 outIds and no bonds hence -->2 outgoing, lowest stable valency = 2 hence no substitutable hydrogen
	 * e.g. phosphonyl has 2 outIds and one double bond -->4 outgoing, lowest stable valency =5 hence 1 substitutable hydrogen
	 * @param atom
	 * @return
	 */
	private static int calculateSubstitutableHydrogenAtoms(Atom atom) {
		Integer valency = null;
		int currentValency =atom.getIncomingValency() + atom.getOutValency();
		if (atom.getValency() != null){
			valency = atom.getValency();
		}
		if (valency ==null){
			Integer defaultValency =ValencyChecker.getDefaultValency(atom.getElement(), atom.getCharge());
			if (defaultValency !=null && currentValency <= defaultValency){
				valency = defaultValency;
			}
			if (valency ==null){
				Integer[] possibleValencies =ValencyChecker.getPossibleValencies(atom.getElement(), atom.getCharge());
				if (possibleValencies!=null) {
					for (Integer possibleValency : possibleValencies) {
						if (currentValency <= possibleValency){
							valency = possibleValency;
							break;
						}
					}
				}
			}
		}
		if (valency ==null){
			return 0;
		}
		else{
			return valency-currentValency;
		}
	}
	
	/**
	 * Stereochemistry terms are assigned right at the end so that checks can be done on whether the indicated atom is in fact chiral.
	 * In the process of multiplication locants are primed. This function adds the appropriate number of primes to any locanted stereochemistry locants
	 * The primesString is the string containing the primes to add to each locant
	 * @param backet
	 * @param primesString
	 */
	private static void addPrimesToLocantedStereochemistryElements(Element bracket, String primesString) {
		List<Element> stereoChemistryElements =XOMTools.getDescendantElementsWithTagName(bracket, "stereoChemistry");
		for (Element stereoChemistryElement : stereoChemistryElements) {
			if (stereoChemistryElement.getParent().equals(bracket)){continue;}
			if (!getLocant(stereoChemistryElement).equals("0")){
				stereoChemistryElement.getAttribute("locant").setValue(getLocant(stereoChemistryElement) + primesString);
			}
		}
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
