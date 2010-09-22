package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.Stack;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import uk.ac.cam.ch.wwmm.opsin.WordRules.WordRule;
import static uk.ac.cam.ch.wwmm.opsin.XmlDeclarations.*;


import nu.xom.Attribute;
import nu.xom.Element;
import nu.xom.Elements;

class StructureBuildingMethods {
	private final static Pattern matchComma =Pattern.compile(",");
	private final static Pattern matchCompoundLocant =Pattern.compile("[\\[\\(\\{](\\d+[a-z]?'*)[\\]\\)\\}]");
	private StructureBuildingMethods() {}

	/**
	 * Resolves a word/bracket:
	 * Locanted attributes of words are resolved onto their group
	 * Locanted substitution is performed
	 * Connections involving multi radicals are processed
	 * Unlocanted attributes of words are resolved onto their group
	 *
	 * If word is a wordRule the function will instantly return
	 *
	 * @param state
	 * @param word
	 * @throws StructureBuildingException
	 */
	static void resolveWordOrBracket(BuildState state, Element word) throws StructureBuildingException {
		if (word.getLocalName().equals(WORDRULE_EL)){//already been resolved
			return;
		}
		if (!word.getLocalName().equals(WORD_EL) && !word.getLocalName().equals(BRACKET_EL)){
			throw new StructureBuildingException("A word or bracket is the expected input");
		}
		recursivelyResolveLocantedFeatures(state, word);
		recursivelyResolveUnLocantedFeatures(state, word);
		//TODO check all things that can substitute have outAtoms
		//TOOD think whether you can avoid the need to have a cansubstitute function by only using appropriate group
		List<Element> subsBracketsAndRoots = XOMTools.getDescendantElementsWithTagNames(word, new String[]{BRACKET_EL, SUBSTITUENT_EL, ROOT_EL});
        for (Element subsBracketsAndRoot : subsBracketsAndRoots) {
            if (subsBracketsAndRoot.getAttribute(MULTIPLIER_ATR) != null) {
                throw new StructureBuildingException("Structure building problem: multiplier on :" + subsBracketsAndRoot.getLocalName() + " was never used");
            }
        }
		List<Element> groups = XOMTools.getDescendantElementsWithTagName(word, GROUP_EL);
		for (int i = 0; i < groups.size(); i++) {
			Element group = groups.get(i);
			if (group.getAttribute(RESOLVED_ATR)==null && i!=groups.size()-1){
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
		if (!word.getLocalName().equals(WORD_EL) && !word.getLocalName().equals(BRACKET_EL)){
			throw new StructureBuildingException("A word or bracket is the expected input");
		}
		List<Element> subsBracketsAndRoots = XOMTools.getChildElementsWithTagNames(word, new String[]{BRACKET_EL, SUBSTITUENT_EL, ROOT_EL});
		//substitution occurs left to right so by doing this right to left you ensure that any groups that will come into existence
		//due to multipliers being expanded will be in existence
		for (int i =subsBracketsAndRoots.size()-1; i>=0; i--) {
			Element subBracketOrRoot = subsBracketsAndRoots.get(i);
			if (subBracketOrRoot.getLocalName().equals(BRACKET_EL)){
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
		if (!word.getLocalName().equals(WORD_EL) && !word.getLocalName().equals(BRACKET_EL)){
			throw new StructureBuildingException("A word or bracket is the expected input");
		}
		List<Element> subsBracketsAndRoots = XOMTools.getChildElementsWithTagNames(word, new String[]{BRACKET_EL, SUBSTITUENT_EL, ROOT_EL});
		//substitution occurs left to right so by doing this right to left you ensure that any groups that will come into existence
		//due to multipliers being expanded will be in existence
		for (int i =subsBracketsAndRoots.size()-1; i>=0; i--) {
			Element subBracketOrRoot = subsBracketsAndRoots.get(i);
			if (subBracketOrRoot.getLocalName().equals(BRACKET_EL)){
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
			performAdditiveOperations(state, subOrRoot);//e.g. ethylenediimino, oxyethylene (operations where two outAtoms are used to produce the bond and no locant is required as groups)
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
		if (subBracketOrRoot.getLocalName().equals(BRACKET_EL)){
			group =findRightMostGroupInBracket(subBracketOrRoot);
		}
		else{
			group =subBracketOrRoot.getFirstChildElement(GROUP_EL);
		}
		if (group.getAttribute(RESOLVED_ATR)!=null){
			return;
		}
		Fragment frag = state.xmlFragmentMap.get(group);
		if (frag.getOutAtoms().size() >=1 && subBracketOrRoot.getAttribute(LOCANT_ATR)!=null){
			String locantString = subBracketOrRoot.getAttributeValue(LOCANT_ATR);
			if (frag.getOutAtoms().size() >1){
				checkAndApplySpecialCaseWhereOutAtomsCanBeCombinedOrThrow(frag, group.getAttributeValue(SUBTYPE_ATR));
			}
			if (subBracketOrRoot.getAttribute(MULTIPLIER_ATR)!=null){//e.g. 1,2-diethyl
				multiplyOutAndSubstitute(state, subBracketOrRoot);
			}
			else{
				Fragment parentFrag = findFragmentWithLocant(state, subBracketOrRoot, locantString);
				if (parentFrag==null){
					throw new StructureBuildingException("Cannot find in scope fragment with atom with locant " + locantString + ".");
				}
				group.addAttribute(new Attribute(RESOLVED_ATR, "yes"));
				Element groupToAttachTo = state.xmlFragmentMap.getElement(parentFrag);
				if (groupToAttachTo.getAttribute(ACCEPTSADDITIVEBONDS_ATR)!=null && parentFrag.getOutAtoms().size()>0 && groupToAttachTo.getAttribute(ISAMULTIRADICAL_ATR)!=null
						&& parentFrag.getAtomByLocantOrThrow(locantString).getOutValency()>0 && frag.getOutAtom(0).getValency()==1){
					//horrible special case to allow C-methylcarbonimidoyl
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
		if (subBracketOrRoot.getLocalName().equals(BRACKET_EL)){
			group =findRightMostGroupInBracket(subBracketOrRoot);
		}
		else{
			group =subBracketOrRoot.getFirstChildElement(GROUP_EL);
		}
		if (group.getAttribute(RESOLVED_ATR)!=null){
			return;
		}
		Fragment frag = state.xmlFragmentMap.get(group);
		if (frag.getOutAtoms().size() >=1){
			if (subBracketOrRoot.getAttribute(LOCANT_ATR)!=null){
				throw new StructureBuildingException("Substituent has an unused outAtom and has a locant but locanted susbtitution should already been been performed!");
			}
			if (frag.getOutAtoms().size() > 1){
				checkAndApplySpecialCaseWhereOutAtomsCanBeCombinedOrThrow(frag, group.getAttributeValue(SUBTYPE_ATR));
			}
			if (subBracketOrRoot.getAttribute(MULTIPLIER_ATR)!=null){//e.g. diethyl
				multiplyOutAndSubstitute(state, subBracketOrRoot);
			}
			else{
				Atom atomToJoinTo = findAtomForSubstitution(state, subBracketOrRoot, frag.getOutAtom(0).getValency());
				if (atomToJoinTo ==null){
					throw new StructureBuildingException("Unlocanted substitution failed: unable to find suitable atom to bond atom with id:" + frag.getOutAtom(0).getAtom().getID() + " to!");
				}
				group.addAttribute(new Attribute(RESOLVED_ATR, "yes"));
				joinFragmentsSubstitutively(state, frag, atomToJoinTo);
			}
		}
	}

	/**
	 * Multiplies out groups/brakets and substitutes them. The attribute "locant" is checked for locants
	 * If it is present it should contain a comma separated list of locants
	 * The strategy employed is to clone subOrBracket and its associated fragments as many times as the multiplier attribute
	 * perform(Un)LocantedSubstitutiveOperations is then called with on each call a different clone (or the original element) being in position
	 * Hence bonding between the clones is impossible
	 * @param state
	 * @param subOrBracket
	 * @throws StructureBuildingException
	 */
	private static void multiplyOutAndSubstitute(BuildState state, Element subOrBracket) throws StructureBuildingException {
		int multiplier = Integer.parseInt(subOrBracket.getAttributeValue(MULTIPLIER_ATR));
		subOrBracket.removeAttribute(subOrBracket.getAttribute(MULTIPLIER_ATR));
		String[] locants =null;
		if (subOrBracket.getAttribute(LOCANT_ATR) !=null){
			locants = matchComma.split(subOrBracket.getAttributeValue(LOCANT_ATR));
		}
		Element parentWordOrBracket =(Element) subOrBracket.getParent();
		int indexOfSubOrBracket = parentWordOrBracket.indexOf(subOrBracket);
		subOrBracket.detach();

		List<Element> elementsNotToBeMultiplied = new ArrayList<Element>();//anything before the multiplier in the sub/bracket
		Element multiplierEl = subOrBracket.getFirstChildElement(MULTIPLIER_EL);
		if (multiplierEl ==null){
			throw new StructureBuildingException("Multiplier not found where multiplier expected");
		}
		for (int j = subOrBracket.indexOf(multiplierEl) -1 ; j >=0 ; j--) {
			Element el = (Element) subOrBracket.getChild(j);
			el.detach();
			elementsNotToBeMultiplied.add(el);
		}
		multiplierEl.detach();

		List<Element> multipliedElements = new ArrayList<Element>();
		for (int i = multiplier -1; i >=0; i--) {
			Element currentElement;
			if (i!=0){
				currentElement = state.fragManager.cloneElement(state, subOrBracket, StringTools.multiplyString("'", i));
				addPrimesToLocantedStereochemistryElements(currentElement, StringTools.multiplyString("'", i));//Stereochemistry elements with locants will need to have their locants primed (stereochemistry is only processed after structure building)
			}
			else{
				currentElement = subOrBracket;
			}
			multipliedElements.add(currentElement);
			parentWordOrBracket.insertChild(currentElement, indexOfSubOrBracket);
			if (locants !=null){
				currentElement.getAttribute(LOCANT_ATR).setValue(locants[i]);
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
		for (Element el : elementsNotToBeMultiplied) {//re-add anything before multiplier to original subOrBracket
			subOrBracket.insertChild(el, 0);
		}
	}

	/**
	 * Adds locanted unsaturators, heteroatoms and hydrogen elements to the group within the sub or root
	 * @param state
	 * @param subOrRoot
	 * @throws StructureBuildingException
	 */
	static void resolveLocantedFeatures(BuildState state, Element subOrRoot) throws StructureBuildingException {
		Elements groups = subOrRoot.getChildElements(GROUP_EL);
		if (groups.size()!=1){
			throw new StructureBuildingException("Each sub or root should only have one group element. This indicates a bug in OPSIN");
		}
		Element group = groups.get(0);
		Fragment thisFrag = state.xmlFragmentMap.get(group);

		ArrayList<Element> unsaturators = new ArrayList<Element>();
		ArrayList<Element> heteroatoms = new ArrayList<Element>();
		ArrayList<Element> hydrogenElements = new ArrayList<Element>();
		ArrayList<Element> dehydroElements = new ArrayList<Element>();

		Elements children =subOrRoot.getChildElements();
		for (int i = 0; i < children.size(); i++) {
			Element currentEl =children.get(i);
			String elName =currentEl.getLocalName();
			if (elName.equals(UNSATURATOR_EL)){
				unsaturators.add(currentEl);
			}
			else if (elName.equals(HETEROATOM_EL)){
				heteroatoms.add(currentEl);
			}
			else if (elName.equals(HYDRO_EL)){
				if (!currentEl.getValue().equals("dehydro")){
					hydrogenElements.add(currentEl);
				}
				else{
					dehydroElements.add(currentEl);
				}
			}
			else if (elName.equals(HYDROGEN_EL)){
				hydrogenElements.add(currentEl);
			}
			else if (elName.equals(INDICATEDHYDROGEN_EL)){
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
				Atom a =thisFrag.getAtomByLocantOrThrow(locant);
				if (a.hasSpareValency()){
					a.setSpareValency(false);
				}
				else{
					throw new StructureBuildingException("hydrogen addition at locant: " + locant +" was requested, but this atom is not unsaturated");
				}
				hydrogenElements.remove(hydrogen);
				hydrogen.detach();
			}
		}
		
		List<Atom> atomsToFormTripleBondsBetween =new ArrayList<Atom>();//dehydro on a double/aromatic bond forms a triple bond
		for(int i=dehydroElements.size() -1;i >= 0;i--) {
			Element dehydro = dehydroElements.get(i);
			String locant = getLocant(dehydro);
			if(!locant.equals("0")) {
				Atom a =thisFrag.getAtomByLocantOrThrow(locant);
				if (!a.hasSpareValency()){
					a.setSpareValency(true);
				}
				else{
					atomsToFormTripleBondsBetween.add(a);
				}
				hydrogenElements.remove(dehydro);
				dehydro.detach();
			}
			else{
				throw new StructureBuildingException("locants are assumed to be required for the use of dehydro to be unambiguous");
			}
		}
		addDehydroInducedTripleBonds(thisFrag, atomsToFormTripleBondsBetween);

		for(int i=unsaturators.size() -1;i >= 0;i--) {
			Element unsaturator = unsaturators.get(i);
			String locant = getLocant(unsaturator);
			int bondOrder = Integer.parseInt(unsaturator.getAttributeValue(VALUE_ATR));
			if(bondOrder <= 1) {
				unsaturator.detach();
				continue;
			}
			if(!locant.equals("0")){
				unsaturators.remove(unsaturator);
				/*
				 * Is the locant a compound locant e.g. 1(6) 
				 * This would indicated unsaturation between the atoms with locants 1 and 6
				 */
	            Matcher matcher = matchCompoundLocant.matcher(locant);
	            if (matcher.find()) {
	            	String compoundLocant = matcher.group(1);
	            	locant = matcher.replaceAll("");
	            	Integer idOfFirstAtomInMultipleBond=thisFrag.getIDFromLocantOrThrow(locant);
	            	FragmentTools.unsaturate(idOfFirstAtomInMultipleBond, compoundLocant, bondOrder, thisFrag);
	            }
	            else{
	            	Integer idOfFirstAtomInMultipleBond=thisFrag.getIDFromLocantOrThrow(locant);
	            	FragmentTools.unsaturate(idOfFirstAtomInMultipleBond, bondOrder, thisFrag);
	            }
				unsaturator.detach();
			}
		}

		for(int i=heteroatoms.size() -1;i >= 0;i--) {
			Element heteroatom = heteroatoms.get(i);
			String locant = getLocant(heteroatom);
			if(!locant.equals("0")) {
				String atomSMILES = heteroatom.getAttributeValue(VALUE_ATR);
				state.fragManager.makeHeteroatom(thisFrag.getAtomByLocantOrThrow(locant), atomSMILES, true);
				if (heteroatom.getAttribute(LAMBDA_ATR)!=null){
					thisFrag.getAtomByLocantOrThrow(locant).setLambdaConventionValency(Integer.parseInt(heteroatom.getAttributeValue(LAMBDA_ATR)));
				}
				heteroatoms.remove(heteroatom);
				heteroatom.detach();
			}
		}
	}

	/**
	 * Attempts to form triple bond between the atoms in atomsToFormTripleBondsBetween
	 * Throws an exception if the list contains duplicates or atoms with no adjacent atom in the list
	 * @param thisFrag
	 * @param atomsToFormTripleBondsBetween
	 * @throws StructureBuildingException
	 */
	private static void addDehydroInducedTripleBonds(Fragment thisFrag,List<Atom> atomsToFormTripleBondsBetween) throws StructureBuildingException {
		if (atomsToFormTripleBondsBetween.size()>0){
			Set<Atom> atoms = new HashSet<Atom>(atomsToFormTripleBondsBetween);
			if (atomsToFormTripleBondsBetween.size() != atoms.size()){
				throw new StructureBuildingException("locants specified for dehydro specify the same atom too many times");
			}
			atomLoop: for (int i = atomsToFormTripleBondsBetween.size()-1; i >=0; i = i-2) {//two atoms will have a triple bond formed betwen them
				Atom a  = atomsToFormTripleBondsBetween.get(i);
				List<Atom> neighbours = a.getAtomNeighbours();
				for (Atom neighbour : neighbours) {
					if (atomsToFormTripleBondsBetween.contains(neighbour)){
						atomsToFormTripleBondsBetween.remove(i);
						atomsToFormTripleBondsBetween.remove(neighbour);
						Bond b = thisFrag.findBondOrThrow(a, neighbour);
						b.setOrder(3);
						a.setSpareValency(false);
						neighbour.setSpareValency(false);
						continue atomLoop;
					}
				}
				throw new StructureBuildingException("dehydro indicated atom should form a triple bond but no adjacent atoms also had hydrogen removed!");
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
		Elements groups = subOrRoot.getChildElements(GROUP_EL);
		if (groups.size()!=1){
			throw new StructureBuildingException("Each sub or root should only have one group element. This indicates a bug in OPSIN");
		}
		Element group = groups.get(0);
		Fragment thisFrag = state.xmlFragmentMap.get(group);

		ArrayList<Element> unsaturators = new ArrayList<Element>();
		ArrayList<Element> heteroatoms = new ArrayList<Element>();
		ArrayList<Element> hydrogenElements = new ArrayList<Element>();

		Elements children =subOrRoot.getChildElements();
		for (int i = 0; i < children.size(); i++) {
			Element currentEl =children.get(i);
			String elName =currentEl.getLocalName();
			if (elName.equals(UNSATURATOR_EL)){
				unsaturators.add(currentEl);
			}
			else if (elName.equals(HETEROATOM_EL)){
				heteroatoms.add(currentEl);
			}
			else if (elName.equals(HYDRO_EL)){
				if (!currentEl.getValue().equals("dehydro")){
					hydrogenElements.add(currentEl);
				}
			}
			else if (elName.equals(HYDROGEN_EL)){
				hydrogenElements.add(currentEl);
			}
			else if (elName.equals(INDICATEDHYDROGEN_EL)){
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
				if (atom.getType().equals(SUFFIX_TYPE_VAL)){
					break;
				}
				atom.ensureSVIsConsistantWithValency(false);//doesn't take into account suffixes
				if (atom.hasSpareValency()){
					if (atomWillHaveSVImplicitlyRemoved(atom)){
						atomsWhichImplicitlyWillHaveTheirSVRemoved.add(atom);
					}
					else{
						atomsWithSV.add(atom);
					}
				}
			}
			atomsWithSV.addAll(atomsWhichImplicitlyWillHaveTheirSVRemoved);//these end up at the end of the list
			boolean saturateAllAtoms =false;
			for (Element hydrogenElement : hydrogenElements) {
				if (hydrogenElement.getValue().equals("perhydro")){
					saturateAllAtoms =true;
					hydrogenElement.detach();
				}
			}
			if (saturateAllAtoms){
				if (hydrogenElements.size() != 1){
					throw new StructureBuildingException("Unexpected indication of hydrogen when perhydro makes such indication redundnant");
				}
	            for (Atom atomToReduceSpareValencyOn : atomsWithSV) {
	                atomToReduceSpareValencyOn.setSpareValency(false);
	            }
			}
			else{
				if (hydrogenElements.size()> atomsWithSV.size()){
					throw new StructureBuildingException("Cannot find atom to add hydrogen to (" +
							hydrogenElements.size() + " hydrogen adding tags but only " +  atomsWithSV.size() +" positions that can be hydrogenated)" );
				}
	            for (Element hydrogenElement : hydrogenElements) {
	                Atom atomToReduceSpareValencyOn = atomsWithSV.removeFirst();
	                atomToReduceSpareValencyOn.setSpareValency(false);
	                hydrogenElement.detach();
	            }
			}
		}

		int idOfFirstAtomInFragment= thisFrag.getIdOfFirstAtom();
		//defaultId corresponds to an id on thisFrag; whenever it is used it is incremented
		//all non-suffix atoms in the fragment are eventually checked if the defaultId does not correspond to a suitable atom
		//e.g. if the id was 3 in a 5 atom fragment which had ids 1-5, atoms would be checked for suitability in the order 3,4,5,1,2
		int defaultId = idOfFirstAtomInFragment;

        for (Element unsaturator : unsaturators) {
            int bondOrder = Integer.parseInt(unsaturator.getAttributeValue(VALUE_ATR));
            if (bondOrder <= 1) {
            	unsaturator.detach();
                continue;
            }
            //checks if both atoms can accept an extra bond (if double bond) or two extra bonds (if triple bond)

            Atom currentAtom = thisFrag.getAtomByIDOrThrow(defaultId);
            Atom nextAtom = thisFrag.getAtomByIDOrThrow(defaultId + 1);
            while (currentAtom.hasSpareValency() || !ValencyChecker.checkValencyAvailableForBond(currentAtom, bondOrder - 1 + currentAtom.getOutValency()) ||
                    nextAtom.hasSpareValency() || !ValencyChecker.checkValencyAvailableForBond(nextAtom, bondOrder - 1 + nextAtom.getOutValency())) {
                defaultId++;
                currentAtom = thisFrag.getAtomByIDOrThrow(defaultId);
                nextAtom = thisFrag.getAtomByIDOrThrow(defaultId + 1);
                if (currentAtom.getType().equals(SUFFIX_TYPE_VAL) || nextAtom.getType().equals(SUFFIX_TYPE_VAL)) {
                    throw new StructureBuildingException("No suitable atom found");
                }
            }
            Integer idOfFirstAtomInMultipleBond = currentAtom.getID();
            FragmentTools.unsaturate(idOfFirstAtomInMultipleBond, bondOrder, thisFrag);
            defaultId = idOfFirstAtomInMultipleBond + 2;
            unsaturator.detach();
        }
		defaultId = idOfFirstAtomInFragment;

        for (Element heteroatom : heteroatoms) {
            String atomSMILES = heteroatom.getAttributeValue(VALUE_ATR);
            //finds an atom for which changing it to the specified heteroatom will not cause valency to be violated
            Atom atomToReplaceWithHeteroAtom = thisFrag.getAtomByIDOrThrow(defaultId);
            while (!ValencyChecker.checkValencyAvailableForReplacementByHeteroatom(atomToReplaceWithHeteroAtom, atomSMILES)) {
                defaultId++;
                atomToReplaceWithHeteroAtom = thisFrag.getAtomByIDOrThrow(defaultId);
                if (atomToReplaceWithHeteroAtom.getType().equals(SUFFIX_TYPE_VAL)) {
                    throw new StructureBuildingException("No suitable atom found");
                }
            }
            state.fragManager.makeHeteroatom(atomToReplaceWithHeteroAtom, atomSMILES, true);
            if (heteroatom.getAttribute(LABELS_ATR) != null) {
                atomToReplaceWithHeteroAtom.setLambdaConventionValency(Integer.parseInt(heteroatom.getAttributeValue(LABELS_ATR)));
            }
            defaultId++;
            heteroatom.detach();
        }
		defaultId = idOfFirstAtomInFragment;

		if (thisFrag.getOutAtoms().size()>0){//assign any outAtoms that have not been set to a specific atom to a specific atom
			for (OutAtom outAtom : thisFrag.getOutAtoms()) {
				if (!outAtom.isSetExplicitly()){
					defaultId=outAtom.getAtom().getID();
					Atom atomToAssociateOutAtomWith=thisFrag.getAtomByIDOrThrow(defaultId);
					while (!ValencyChecker.checkValencyAvailableForBond(atomToAssociateOutAtomWith, (atomToAssociateOutAtomWith.hasSpareValency() ? 1 : 0) + atomToAssociateOutAtomWith.getOutValency() + outAtom.getValency())){
						defaultId++;
						atomToAssociateOutAtomWith=thisFrag.getAtomByID(defaultId);
						if (atomToAssociateOutAtomWith == null){
							throw new StructureBuildingException("Failed to assign all unlocanted radicals to actual atoms without violating valency");
						}
						if (atomToAssociateOutAtomWith.getType().equals(SUFFIX_TYPE_VAL)){
							throw new StructureBuildingException("No suitable atom found");
						}
					}
					outAtom.setAtom(atomToAssociateOutAtomWith);
					outAtom.setSetExplicitly(true);
					defaultId = idOfFirstAtomInFragment;
				}
			}
		}
	}

	private static boolean atomWillHaveSVImplicitlyRemoved(Atom atom) throws StructureBuildingException {
		boolean canFormDoubleBond =false;
		for(Atom aa : atom.getFrag().getIntraFragmentAtomNeighbours(atom)) {
			if(aa.hasSpareValency()){
				canFormDoubleBond=true;
			}
		}
		if (!canFormDoubleBond){
			return true;
		}
		
		if (atom.hasSpareValency()){
			atom.ensureSVIsConsistantWithValency(true);
			if (!atom.hasSpareValency()){
				atom.setSpareValency(true);
				return true;
			}
		}
		return false;
	}

	private static void performAdditiveOperations(BuildState state, Element subBracketOrRoot) throws StructureBuildingException {
		if (subBracketOrRoot.getAttribute(LOCANT_ATR)!=null){//additive nomenclature does not employ locants
			return;
		}
		Element group;
		if (subBracketOrRoot.getLocalName().equals(BRACKET_EL)){
			group =findRightMostGroupInBracket(subBracketOrRoot);
		}
		else{
			group =subBracketOrRoot.getFirstChildElement(GROUP_EL);
		}
		if (group.getAttribute(RESOLVED_ATR)!=null){
			return;
		}
		Fragment frag = state.xmlFragmentMap.get(group);
		int outAtomCount = frag.getOutAtoms().size();
		if (outAtomCount >=1){
			if (subBracketOrRoot.getAttribute(MULTIPLIER_ATR) ==null){
				Element nextSiblingEl = (Element) XOMTools.getNextSibling(subBracketOrRoot);
				if (nextSiblingEl.getAttribute(MULTIPLIER_ATR)!=null &&
						(outAtomCount >= Integer.parseInt(nextSiblingEl.getAttributeValue(MULTIPLIER_ATR)) || //probably multiplicative nomenclature, should be as many outAtoms as the multiplier
						outAtomCount==1 && frag.getOutAtom(0).getValency()==Integer.parseInt(nextSiblingEl.getAttributeValue(MULTIPLIER_ATR))) &&
						hasRootLikeOrMultiRadicalGroup(state, nextSiblingEl)){
					if (outAtomCount==1){//special case e.g. 4,4'-(benzylidene)dianiline
						FragmentTools.splitOutAtomIntoValency1OutAtoms(frag.getOutAtom(0));
						//special case where something like benzylidene is being used as if it meant benzdiyl for multiplicative nomenclature
						//this is allowed in the IUPAC 79 recommendations but not recommended in the current recommendations
					}
					performMultiplicativeOperations(state, group, nextSiblingEl);
				}
				else if (group.getAttribute(ISAMULTIRADICAL_ATR)!=null){//additive nomenclature e.g. ethyleneoxy
					Fragment nextFrag = getNextInScopeMultiValentFragment(state, subBracketOrRoot);
					if (nextFrag!=null){
						Element nextMultiRadicalGroup = state.xmlFragmentMap.getElement(nextFrag);
						Element parentSubOrRoot = (Element) nextMultiRadicalGroup.getParent();
						if (state.currentWordRule != WordRule.polymer){//imino does not behave like a substituent in polymers only as a linker
							if (nextMultiRadicalGroup.getAttribute(IMINOLIKE_ATR)!=null){//imino/methylene can just act as normal substituents, should an additive bond really be made???
								Fragment adjacentFrag =state.xmlFragmentMap.get(OpsinTools.getNextGroup(subBracketOrRoot));
								
								if (nextFrag !=adjacentFrag){//imino is not the absolute next frag
									if (potentiallyCanSubstitute((Element) nextMultiRadicalGroup.getParent()) || potentiallyCanSubstitute((Element) nextMultiRadicalGroup.getParent().getParent())){
										return;
									}
								}
							}
							if (group.getAttribute(IMINOLIKE_ATR)!=null && levelsToWordEl(group) > levelsToWordEl(nextMultiRadicalGroup)){
								return;//e.g. imino substitutes ((chloroimino)ethylene)dibenzene
							}
						}
						if (parentSubOrRoot.getAttribute(MULTIPLIER_ATR)!=null){
							throw new StructureBuildingException("Attempted to form additive bond to a multiplied component");
						}
						group.addAttribute(new Attribute(RESOLVED_ATR, "yes"));
						joinFragmentsAdditively(state, frag, nextFrag);
					}
				}
				else {//e.g. chlorocarbonyl or hydroxy(sulfanyl)phosphoryl
					List<Fragment> siblingFragments = findAlternativeFragments(state, subBracketOrRoot);
					if (siblingFragments.size()>0){
						Fragment nextFrag = siblingFragments.get(siblingFragments.size()-1);
						Element nextGroup = state.xmlFragmentMap.getElement(nextFrag);
						if (nextGroup.getAttribute(ACCEPTSADDITIVEBONDS_ATR)!=null && nextFrag!=null && nextGroup.getAttribute(ISAMULTIRADICAL_ATR)!=null  && (nextFrag.getOutAtoms().size()>1|| nextGroup.getAttribute(RESOLVED_ATR)!=null && nextFrag.getOutAtoms().size()>=1 )){
							Atom toAtom = nextFrag.getOutAtom(0).getAtom();
							if (calculateSubstitutableHydrogenAtoms(toAtom) ==0){
								group.addAttribute(new Attribute(RESOLVED_ATR, "yes"));
								joinFragmentsAdditively(state, frag, nextFrag);//e.g. aminocarbonyl or aminothio
							}
						}
						if (group.getAttribute(RESOLVED_ATR)==null && siblingFragments.size()>1){
							for (int i = 0; i< siblingFragments.size()-1; i++) {
								Fragment lastFrag = siblingFragments.get(i);
								Element lastGroup = state.xmlFragmentMap.getElement(lastFrag);
								if (lastGroup.getAttribute(ACCEPTSADDITIVEBONDS_ATR)!=null && lastFrag!=null && lastGroup.getAttribute(ISAMULTIRADICAL_ATR)!=null  && (lastFrag.getOutAtoms().size()>1|| lastGroup.getAttribute(RESOLVED_ATR)!=null && lastFrag.getOutAtoms().size()>=1 )){
									Atom toAtom = lastFrag.getOutAtom(0).getAtom();
									if (calculateSubstitutableHydrogenAtoms(toAtom) ==0){
										group.addAttribute(new Attribute(RESOLVED_ATR, "yes"));
										joinFragmentsAdditively(state, frag, lastFrag);//e.g. hydroxy(sulfanyl)phosphoryl
									}
									break;
								}

								//loop may continue if lastFrag was in fact completely unsubtitutable e.g. hydroxy...phosphoryloxy. The oxy is unsubstituable as the phosphoryl will already have bonded to it
								if (lastFrag.getAtomOrNextSuitableAtom(lastFrag.getDefaultInAtom(), frag.getOutAtom(outAtomCount-1).getValency(), true)!=null){
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
					int multiplier = Integer.parseInt(subBracketOrRoot.getAttributeValue(MULTIPLIER_ATR));
					Fragment nextFrag = siblingFragments.get(siblingFragments.size()-1);
					Element nextGroup = state.xmlFragmentMap.getElement(nextFrag);
					if (nextGroup.getAttribute(ACCEPTSADDITIVEBONDS_ATR)!=null && nextFrag!=null && nextGroup.getAttribute(ISAMULTIRADICAL_ATR)!=null  && (nextFrag.getOutAtoms().size()>=multiplier|| nextGroup.getAttribute(RESOLVED_ATR)!=null && nextFrag.getOutAtoms().size()>=multiplier +1 )){
						Atom toAtom = nextFrag.getOutAtom(0).getAtom();
						if (calculateSubstitutableHydrogenAtoms(toAtom) ==0){
							group.addAttribute(new Attribute(RESOLVED_ATR, "yes"));
							multiplyOutAndAdditivelyBond(state, subBracketOrRoot, nextFrag);//e.g.dihydroxyphosphoryl
						}
					}
					if (group.getAttribute(RESOLVED_ATR)==null && siblingFragments.size()>1){
						for (int i = 0; i< siblingFragments.size()-1; i++) {
							Fragment lastFrag = siblingFragments.get(i);
							Element lastGroup = state.xmlFragmentMap.getElement(lastFrag);
							if (lastGroup.getAttribute(ACCEPTSADDITIVEBONDS_ATR)!=null && lastFrag!=null && lastGroup.getAttribute(ISAMULTIRADICAL_ATR)!=null  &&  (lastFrag.getOutAtoms().size()>=multiplier|| lastGroup.getAttribute(RESOLVED_ATR)!=null && lastFrag.getOutAtoms().size()>=multiplier +1 )){
								Atom toAtom = lastFrag.getOutAtom(0).getAtom();
								if (calculateSubstitutableHydrogenAtoms(toAtom) ==0){
									group.addAttribute(new Attribute(RESOLVED_ATR, "yes"));
									multiplyOutAndAdditivelyBond(state, subBracketOrRoot, lastFrag);//e.g. dihydroxyphosphoryloxy
								}
								break;
							}

							//loop may continue if lastFrag was in fact completely unsubtitutable e.g. hydroxy...phosphoryloxy. The oxy is unsubstituable as the phosphoryl will already have bonded to it
							if (lastFrag.getAtomOrNextSuitableAtom(lastFrag.getDefaultInAtom(), frag.getOutAtom(outAtomCount-1).getValency(), true)!=null){
								break;
							}
						}
					}
				}
			}
		}
	}

	/**
	 * Searches the input for something that either is a multiRadical or has no outAtoms i.e. not dimethyl
	 * @param state
	 * @param subBracketOrRoot
	 * @return
	 */
	private static boolean hasRootLikeOrMultiRadicalGroup(BuildState state, Element subBracketOrRoot) {
		List<Element> groups = XOMTools.getDescendantElementsWithTagName(subBracketOrRoot, GROUP_EL);
		if (subBracketOrRoot.getAttribute(INLOCANTS_ATR)!=null){
			return true;// a terminus with specified inLocants
		}
		for (Element group : groups) {
			Fragment frag =state.xmlFragmentMap.get(group);
			int outAtomCount =frag.getOutAtoms().size();
			if (group.getAttribute(ISAMULTIRADICAL_ATR)!=null){
				if (outAtomCount >=1 ){
					return true;//a multi radical
				}
			}
			else if (outAtomCount ==0 && group.getAttribute(RESOLVED_ATR)==null){
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
		int multiplier = Integer.parseInt(subOrBracket.getAttributeValue(MULTIPLIER_ATR));
		subOrBracket.removeAttribute(subOrBracket.getAttribute(MULTIPLIER_ATR));
		List<Element> clonedElements = new ArrayList<Element>();
		List<Element> elementsNotToBeMultiplied = new ArrayList<Element>();//anything before the multiplier in the sub/bracket
		for (int i = multiplier -1; i >=0; i--) {
			Element currentElement;
			if (i!=0){
				currentElement = state.fragManager.cloneElement(state, subOrBracket, StringTools.multiplyString("'", i));
				addPrimesToLocantedStereochemistryElements(currentElement, StringTools.multiplyString("'", i));//Stereochemistry elements with locants will need to have their locants primed (stereochemistry is only processed after structure building)
				clonedElements.add(currentElement);
			}
			else{
				currentElement = subOrBracket;
				Element multiplierEl = subOrBracket.getFirstChildElement(MULTIPLIER_EL);
				if (multiplierEl ==null){
					throw new StructureBuildingException("Multiplier not found where multiplier expected");
				}
				for (int j = subOrBracket.indexOf(multiplierEl) -1 ; j >=0 ; j--) {
					Element el = (Element) subOrBracket.getChild(j);
					el.detach();
					elementsNotToBeMultiplied.add(el);
				}
				multiplierEl.detach();
			}
			Element group;
			if (currentElement.getLocalName().equals(BRACKET_EL)){
				group = findRightMostGroupInBracket(currentElement);
			}
			else{
				group = currentElement.getFirstChildElement(GROUP_EL);
			}
			Fragment frag = state.xmlFragmentMap.get(group);
			if (frag.getOutAtoms().size() !=1 ){
				throw new StructureBuildingException("Additive bond formation failure: Fragment expected to have one OutAtom in this case but had: "+ frag.getOutAtoms().size());
			}
			joinFragmentsAdditively(state, frag, fragToAdditivelyBondTo);
		}
		for (Element clone : clonedElements) {//make sure cloned substituents don't substitute onto each other!
			XOMTools.insertAfter(subOrBracket, clone);
		}
		for (Element el : elementsNotToBeMultiplied) {//re-add anything before multiplier to original subOrBracket
			subOrBracket.insertChild(el, 0);
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
		int multiplier = Integer.parseInt(multipliedParent.getAttributeValue(MULTIPLIER_ATR));
		if (multiplier > multiRadicalBR.getOutAtomCount()){
			throw new StructureBuildingException("Multiplication bond formation failure: number of outAtoms disagree with multiplier(multiplier: " + multiplier + ", outAtom count: " + multiRadicalBR.getOutAtomCount()+ ") , this is an OPSIN bug");
		}
		if (state.debug){System.out.println(multiplier +" multiplicative bonds to be formed");}
		multipliedParent.removeAttribute(multipliedParent.getAttribute(MULTIPLIER_ATR));
		List<String> inLocants = null;
		if (multipliedParent.getAttribute(INLOCANTS_ATR)!=null){//true for the root of a multiplicative name
			String inLocantsString = multipliedParent.getAttributeValue(INLOCANTS_ATR);
			if (inLocantsString.equals(INLOCANTS_DEFAULT)){
				inLocants = new ArrayList<String>(multiplier);
				for (int i = 0; i < multiplier; i++) {
					inLocants.add(INLOCANTS_DEFAULT);
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
			if (currentElement.getLocalName().equals(BRACKET_EL)){
				group =getFirstMultiValentGroup(state, currentElement);
				if (group == null){//root will not have a multivalent group
					group = findRightMostGroupInBracket(currentElement);
				}
			}
			else{
				group = currentElement.getFirstChildElement(GROUP_EL);
			}
			Fragment frag = state.xmlFragmentMap.get(group);
			if (inLocants !=null){
				if (group.getAttribute(ISAMULTIRADICAL_ATR)!=null){//e.g. methylenedisulfonyl dichloride
					if (!multipliedParent.getAttributeValue(INLOCANTS_ATR).equals(INLOCANTS_DEFAULT)){
						throw new StructureBuildingException("inLocants should not be specified for a multiradical parent in multiplicative nomenclature");
					}
					Element rightMostGroup;
					if (currentElement.getLocalName().equals(BRACKET_EL)){
						rightMostGroup = findRightMostGroupInBracket(currentElement);//this is the only one that can accept inAtoms
					}
					else{
						rightMostGroup = currentElement.getFirstChildElement(GROUP_EL);
					}
					rightMostGroup.addAttribute(new Attribute(RESOLVED_ATR, "yes"));
				}
				else{
					group.addAttribute(new Attribute(RESOLVED_ATR, "yes"));
					boolean inAtomAdded =false;
					for (int j = inLocants.size() -1; j >=0; j--) {
						String locant = inLocants.get(j);
						if (locant.equals(INLOCANTS_DEFAULT)){//note that if one entry in inLocantArray is default then they all are "default"
							frag.addInAtom(frag.getAtomOrNextSuitableAtom(frag.getDefaultInAtom(), 1, true), 1);
							inAtomAdded=true;
							inLocants.remove(j);
							break;
						}
						else{
							Atom inAtom = frag.getAtomByLocant(locant);
							if (inAtom!=null){
								frag.addInAtom(inAtom, 1);
								inAtomAdded=true;
								inLocants.remove(j);
								break;
							}
						}
					}
					if (!inAtomAdded){
						throw new StructureBuildingException("Locants for inAtoms on the root were either misassigned to the root or were invalid: " + inLocants.toString() +" could not be assigned!");
					}
				}
			}
			if (frag.getInAtoms().size()!=1 && frag.getOutAtoms().size() ==0 ){
				throw new StructureBuildingException("Multiplication bond formation failure: OPSIN bug, input to joinFragmentsMultiplicatively was unexpected");
			}

			Element multiRadicalGroup =state.xmlFragmentMap.getElement(multiRadicalBR.getOutAtom(i).getAtom().getFrag());
			if (multiRadicalGroup.getAttribute(RESOLVED_ATR)==null){
				resolveUnLocantedFeatures(state, (Element) multiRadicalGroup.getParent());//the addition of unlocanted unsaturators can effect the position of radicals e.g. diazenyl
				multiRadicalGroup.addAttribute(new Attribute(RESOLVED_ATR, "yes"));
			}
			joinFragmentsAdditively(state, multiRadicalBR.getOutAtom(i).getAtom().getFrag(), frag);
			if (currentElement.getLocalName().equals(BRACKET_EL)){
				recursivelyResolveUnLocantedFeatures(state, currentElement);//there may be outAtoms that are involved in unlocanted substitution, these can be safely used now e.g. ...bis((3-hydroxy-4-methoxyphenyl)methylene) where (3-hydroxy-4-methoxyphenyl)methylene is the currentElement
			}

			if (inLocants ==null){
				//currentElement is not a root element. Need to build up a new BuildResults so as to call performMultiplicativeOperations again
				//at this stage an outAtom has been removed from the fragment within currentElement through an additive bond
				newBr.mergeBuildResults(new BuildResults(state, currentElement));
			}
		}

		if (newBr.getFragmentCount()==1){
			throw new StructureBuildingException("Multiplicative nomenclature cannot yield only one temporary terminal fragment");
		}
		if (newBr.getFragmentCount()>=2){
			List<Element> siblings = XOMTools.getNextSiblingsOfTypes(multipliedParent, new String[]{SUBSTITUENT_EL, BRACKET_EL, ROOT_EL});
			if (siblings.size()==0){
				Element parentOfMultipliedEl = (Element) multipliedParent.getParent();
				if (parentOfMultipliedEl.getLocalName().equals(BRACKET_EL)){//brackets are allowed
					siblings = XOMTools.getNextSiblingsOfTypes(parentOfMultipliedEl, new String[]{SUBSTITUENT_EL, BRACKET_EL, ROOT_EL});
					if (siblings.get(0).getAttribute(MULTIPLIER_ATR)==null){
						throw new StructureBuildingException("Multiplier not found where multiplier was expected for succesful multiplicative nomenclature");
					}
					performMultiplicativeOperations(state, newBr, siblings.get(0));
				}
				else{
					throw new StructureBuildingException("Could not find suitable element to continue multiplicative nomenclature");
				}
			}
			else{
				if (siblings.get(0).getAttribute(MULTIPLIER_ATR)==null){
					throw new StructureBuildingException("Multiplier not found where multiplier was expected for successful multiplicative nomenclature");
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
		if (!substituentOrBracket.getLocalName().equals(SUBSTITUENT_EL) && !substituentOrBracket.getLocalName().equals(BRACKET_EL)){
			throw new StructureBuildingException("Input to this function should be a substituent or bracket");
		}
		if (substituentOrBracket.getParent()==null){
			throw new StructureBuildingException("substituent did not have a parent!");
		}
		Element parent =(Element) substituentOrBracket.getParent();

		List<Element> children = XOMTools.getChildElementsWithTagNames(parent, new String[]{SUBSTITUENT_EL, BRACKET_EL, ROOT_EL});//will be returned in index order
		int indexOfSubstituent =parent.indexOf(substituentOrBracket);
		for (Element child : children) {
			if (parent.indexOf(child) <=indexOfSubstituent){//only want things after the input
				continue;
			}
			if (child.getAttribute(MULTIPLIER_ATR)!=null){
				continue;
			}
			List<Element> childDescendants;
			if (child.getLocalName().equals(BRACKET_EL)){
				childDescendants = XOMTools.getDescendantElementsWithTagNames(child, new String[]{SUBSTITUENT_EL, ROOT_EL});//will be returned in depth-first order
			}
			else{
				childDescendants =new ArrayList<Element>();
				childDescendants.add(child);
			}
			for (Element descendantChild : childDescendants) {
				Element group = descendantChild.getFirstChildElement(GROUP_EL);
				if (group == null){
					throw new StructureBuildingException("substituent/root is missing its group");
				}
				Fragment possibleFrag = state.xmlFragmentMap.get(group);
				if (group.getAttribute(ISAMULTIRADICAL_ATR)!=null &&
						(possibleFrag.getOutAtoms().size() >=2 || (possibleFrag.getOutAtoms().size() >=1 && group.getAttribute(RESOLVED_ATR)!=null ))){
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
		if (!bracket.getLocalName().equals(BRACKET_EL)){
			throw new StructureBuildingException("Input to this function should be a bracket");
		}

		List<Element> groups = XOMTools.getDescendantElementsWithTagName(bracket, GROUP_EL);//will be returned in index order
		for (Element group : groups) {
			Fragment possibleFrag = state.xmlFragmentMap.get(group);
			if (group.getAttribute(ISAMULTIRADICAL_ATR)!=null &&
					(possibleFrag.getOutAtoms().size() >=2 || (possibleFrag.getOutAtoms().size() >=1 && group.getAttribute(RESOLVED_ATR)!=null ))){
				return group;
			}
		}
		return null;
	}

	private static void joinFragmentsAdditively(BuildState state, Fragment fragToBeJoined, Fragment parentFrag) throws StructureBuildingException {
		Element elOfFragToBeJoined = state.xmlFragmentMap.getElement(fragToBeJoined);
		if (EPOXYLIKE_SUBTYPE_VAL.equals(elOfFragToBeJoined.getAttributeValue(SUBTYPE_ATR))){
			List<OutAtom> outAtoms = fragToBeJoined.getOutAtoms();
			for (OutAtom outAtom : outAtoms) {
				if (outAtom.getLocant()!=null){
					throw new StructureBuildingException("Inappropriate use of " + elOfFragToBeJoined.getValue());
				}
			}
		}
		int outAtomCountOnFragToBeJoined = fragToBeJoined.getOutAtoms().size();
		if (outAtomCountOnFragToBeJoined ==0){
			throw new StructureBuildingException("Additive bond formation failure: Fragment expected to have at least one OutAtom but had none");
		}

		Atom to;
		int bondOrder;
		if (parentFrag.getInAtoms().size()==1){//special case for the parent of multiplicative nomenclature. This is really substitutive nomenclature
			InAtom in = parentFrag.getInAtom(0);
			to = in.getAtom();
			parentFrag.removeInAtom(in);
			bondOrder = fragToBeJoined.getOutAtom(outAtomCountOnFragToBeJoined-1).getValency();
		}
		else{
			List<OutAtom> outAtomsOnParent = parentFrag.getOutAtoms();
			if (outAtomsOnParent.size() ==0){
				throw new StructureBuildingException("Additive bond formation failure: Fragment expected to have at least one OutAtom but had none");
			}
			OutAtom in = null;
			if (outAtomsOnParent.size() >1){
				int firstOutAtomOrder = outAtomsOnParent.get(0).getValency();
				boolean unresolvedAmbiguity =false;
				for (OutAtom outAtom : outAtomsOnParent) {
					if (outAtom.getValency()!=firstOutAtomOrder){
						unresolvedAmbiguity =true;
					}
				}
				if (unresolvedAmbiguity){//not all outAtoms on parent equivalent
					List<OutAtom> outAtomsOnfragToBeJoined = fragToBeJoined.getOutAtoms();
					firstOutAtomOrder = outAtomsOnfragToBeJoined.get(0).getValency();
					unresolvedAmbiguity =false;
					for (OutAtom outAtom : outAtomsOnfragToBeJoined) {
						if (outAtom.getValency()!=firstOutAtomOrder){
							unresolvedAmbiguity =true;
						}
					}
					if (unresolvedAmbiguity && outAtomsOnfragToBeJoined.size()==2){//not all outAtoms on frag to be joined are equivalent either!
						//Solves the specific case of 2,2'-[ethane-1,2-diylbis(azanylylidenemethanylylidene)]diphenol vs 2,2'-[ethane-1,2-diylidenebis(azanylylidenemethanylylidene)]bis(cyclohexan-1-ol)
						//but does not solve the general case as only a single look behind is performed.
						Element previousGroup = (Element) OpsinTools.getPreviousGroup(elOfFragToBeJoined);
						if (previousGroup!=null){
							List<OutAtom> previousOutAtoms =  state.xmlFragmentMap.get(previousGroup).getOutAtoms();
							if (previousOutAtoms.size()>1){
								int previousGroupFirstOutAtomOrder = previousOutAtoms.get(0).getValency();
								unresolvedAmbiguity =false;
								for (OutAtom outAtom : previousOutAtoms) {
									if (outAtom.getValency()!=previousGroupFirstOutAtomOrder){
										unresolvedAmbiguity =true;
									}
								}
								if (!unresolvedAmbiguity && previousGroupFirstOutAtomOrder==outAtomsOnParent.get(0).getValency()){
									for (OutAtom outAtom : outAtomsOnParent) {
										if (outAtom.getValency()!=previousGroupFirstOutAtomOrder){
											in = outAtom;
											break;
										}
									}
								}
							}
						}
					}
					else{
						for (OutAtom outAtom : outAtomsOnParent) {
							if (outAtom.getValency()==firstOutAtomOrder){
								in = outAtom;
								break;
							}
						}
					}
				}
			}
			if (in==null){
				in = parentFrag.getOutAtom(0);
			}
			to = in.getAtom();
			bondOrder = in.getValency();
			if (!in.isSetExplicitly()){//not set explicitly so may be an inappropriate atom
				to = to.getFrag().getAtomOrNextSuitableAtomOrThrow(to, bondOrder, false);
			}
			parentFrag.removeOutAtom(in);
		}

		OutAtom out =null;

		for (int i =outAtomCountOnFragToBeJoined -1; i>=0; i--) {
			if (fragToBeJoined.getOutAtom(i).getValency() == bondOrder){
				out = fragToBeJoined.getOutAtom(i);
				break;
			}
		}

		if (out ==null){
			if (outAtomCountOnFragToBeJoined >=bondOrder){//handles cases like nitrilo needing to be -N= (remove later outAtoms first as per usual)
				int valency =0;
				Atom lastOutAtom = fragToBeJoined.getOutAtom(outAtomCountOnFragToBeJoined -1).getAtom();
				for (int i =outAtomCountOnFragToBeJoined -1; i >= 0; i--) {
					OutAtom nextOutAtom = fragToBeJoined.getOutAtom(i);
					if (nextOutAtom.getAtom() !=lastOutAtom){
						throw new StructureBuildingException("Additive bond formation failure: bond order disagreement");
					}
					valency += nextOutAtom.getValency();
					if (valency==bondOrder){
						nextOutAtom.setValency(valency);
						out = nextOutAtom;
						break;
					}
					fragToBeJoined.removeOutAtom(nextOutAtom);
				}
				if (out==null){
					throw new StructureBuildingException("Additive bond formation failure: bond order disagreement");
				}
			}
			else{
				throw new StructureBuildingException("Additive bond formation failure: bond order disagreement");
			}
		}

		Atom from = out.getAtom();
		if (!out.isSetExplicitly()){//not set explicitly so may be an inappropriate atom
			from=from.getFrag().getAtomOrNextSuitableAtomOrThrow(from, bondOrder, false);
		}
		fragToBeJoined.removeOutAtom(out);

		state.fragManager.createBond(from, to, bondOrder);
		if (state.debug){System.out.println("Additively bonded " + from.getID() + " (" +state.xmlFragmentMap.getElement(from.getFrag()).getValue()+") " + to.getID() + " (" +state.xmlFragmentMap.getElement(to.getFrag()).getValue()+")" );}
	}

	private static void joinFragmentsSubstitutively(BuildState state, Fragment fragToBeJoined, Atom atomToJoinTo) throws StructureBuildingException {
		Element elOfFragToBeJoined = state.xmlFragmentMap.getElement(fragToBeJoined);
		if (EPOXYLIKE_SUBTYPE_VAL.equals(elOfFragToBeJoined.getAttributeValue(SUBTYPE_ATR))){
			formEpoxide(state, fragToBeJoined, atomToJoinTo);
			return;
		}
		int outAtomCount = fragToBeJoined.getOutAtoms().size();
		if (outAtomCount >1){
			throw new StructureBuildingException("Substitutive bond formation failure: Fragment expected to have one OutAtom but had: "+ outAtomCount);
		}
		if (outAtomCount ==0 ){
			throw new StructureBuildingException("Substitutive bond formation failure: Fragment expected to have one OutAtom but had none");
		}
		if (state.xmlFragmentMap.getElement(fragToBeJoined).getAttribute(IMINOLIKE_ATR)!=null){//special case for methylene/imino
			if (fragToBeJoined.getOutAtoms().size()==1 && fragToBeJoined.getOutAtom(0).getValency()==1 ){
				fragToBeJoined.getOutAtom(0).setValency(2);
			}
		}
		OutAtom out = fragToBeJoined.getOutAtom(0);
		Atom from = out.getAtom();
		int bondOrder = out.getValency();
		if (!out.isSetExplicitly()){//not set explicitly so may be an inappropriate atom
			from=from.getFrag().getAtomOrNextSuitableAtomOrThrow(from, bondOrder, false);
		}
		fragToBeJoined.removeOutAtom(out);

		state.fragManager.createBond(from, atomToJoinTo, bondOrder);
		if (state.debug){System.out.println("Substitutively bonded " + from.getID() + " (" +state.xmlFragmentMap.getElement(from.getFrag()).getValue()+") " + atomToJoinTo.getID() + " (" +state.xmlFragmentMap.getElement(atomToJoinTo.getFrag()).getValue()+")");}
	}

	private static void formEpoxide(BuildState state, Fragment fragToBeJoined, Atom atomToJoinTo) throws StructureBuildingException {
		Fragment fragToJoinTo = atomToJoinTo.getFrag();
		List<Atom> atomList = fragToJoinTo.getAtomList();
		Atom firstAtomToJoinTo;
		if (fragToBeJoined.getOutAtom(0).getLocant()!=null){
			firstAtomToJoinTo = fragToJoinTo.getAtomByLocantOrThrow(fragToBeJoined.getOutAtom(0).getLocant());
		}
		else{
			firstAtomToJoinTo = fragToJoinTo.getAtomOrNextSuitableAtomOrThrow(atomList.get(0), 1, true);
		}
		Atom chalcogenAtom1 = fragToBeJoined.getOutAtom(0).getAtom();
		fragToBeJoined.removeOutAtom(0);
		Atom secondAtomToJoinTo;
		if (fragToBeJoined.getOutAtom(0).getLocant()!=null){
			secondAtomToJoinTo = fragToJoinTo.getAtomByLocantOrThrow(fragToBeJoined.getOutAtom(0).getLocant());
		}
		else{
			int index = atomList.indexOf(firstAtomToJoinTo);
			if (index +1 >= atomList.size()){
				throw new StructureBuildingException("Unable to find second suitable atom to form epoxide");
			}
			secondAtomToJoinTo = fragToJoinTo.getAtomOrNextSuitableAtomOrThrow(atomList.get(index+1), 1, true);
		}
		Atom chalcogenAtom2 = fragToBeJoined.getOutAtom(0).getAtom();
		fragToBeJoined.removeOutAtom(0);
		if (firstAtomToJoinTo == secondAtomToJoinTo){
			throw new StructureBuildingException("Epoxides must be formed between two different atoms");
		}
		//In epoxy chalcogenAtom1 will be chalcogenAtom2. Methylenedioxy is also handled by this method
		state.fragManager.createBond(chalcogenAtom1, firstAtomToJoinTo, 1);
		state.fragManager.createBond(chalcogenAtom2, secondAtomToJoinTo, 1);
	}

	private static Atom findAtomForSubstitution(BuildState state, Element subOrBracket, int bondOrder) throws StructureBuildingException {
		//case where you should actually be substituting onto the previous element e.g. 5-(4-methylphenylcarbonyl)pentane
		Atom to =null;
		List<Fragment> possibleParents =findAlternativeFragments(state, subOrBracket);
		for (Fragment fragment : possibleParents) {
			to = fragment.getAtomOrNextSuitableAtom(fragment.getDefaultInAtom(), bondOrder, true);
			if (to !=null){
				break;
			}
		}
		return to;
	}

	/**
	 * Finds all the groups accessible from the startingElement taking into account brackets
	 * i.e. those that it is feasible that the group of the startingElement could substitute onto
	 * @param state
	 * @param startingElement
	 * @return A list of fragments in the order to try them as possible parent fragments (for substitutive operations)
	 */
	private static List<Fragment> findAlternativeFragments(BuildState state, Element startingElement) {
		Stack<Element> stack = new Stack<Element>();
		stack.add((Element) startingElement.getParent());
		List<Fragment> foundFragments =new ArrayList<Fragment>();
		boolean doneFirstIteration =false;//check on index only done on first iteration to only get elements with an index greater than the starting element
		while (stack.size()>0){
			Element currentElement =stack.pop();
			if (currentElement.getLocalName().equals(GROUP_EL)){
				Fragment groupFrag =state.xmlFragmentMap.get(currentElement);
				foundFragments.add(groupFrag);
				continue;
			}
			List<Element> siblings = XOMTools.getChildElementsWithTagNames(currentElement, new String[]{BRACKET_EL, SUBSTITUENT_EL, ROOT_EL});

			Stack<Element> bracketted = new Stack<Element>();
			for (Element bracketOrSubOrRoot : siblings) {
				if (!doneFirstIteration && currentElement.indexOf(bracketOrSubOrRoot)<=currentElement.indexOf(startingElement)){
					continue;
				}
				if (bracketOrSubOrRoot.getAttribute(MULTIPLIER_ATR)!=null){
					continue;
				}
				if (bracketOrSubOrRoot.getLocalName().equals(BRACKET_EL)){
					if (IMPLICIT_TYPE_VAL.equals(bracketOrSubOrRoot.getAttributeValue(TYPE_ATR))){
						stack.add(bracketOrSubOrRoot);
					}
					else{
						bracketted.add(bracketOrSubOrRoot);
					}
				}
				else{
					Element group = bracketOrSubOrRoot.getFirstChildElement(GROUP_EL);
					stack.add(group);
				}
			}
			stack.addAll(0, bracketted);//locanting into brackets is rarely the desired answer so place at the bottom of the stack
			doneFirstIteration =true;
		}
		return foundFragments;
	}

	/**
	 * Checks through the groups accessible from the currentElement taking into account brackets
	 * i.e. those that it is feasible that the group of the currentElement could substitute onto
	 * @param state
	 * @param startingElement
	 * @param locant: the locant string to check for the presence of
	 * @return The fragment with the locant, or null
	 * @throws StructureBuildingException
	 */
	private static Fragment findFragmentWithLocant(BuildState state, Element startingElement, String locant) throws StructureBuildingException {
		Stack<Element> stack = new Stack<Element>();
		stack.add((Element) startingElement.getParent());
		boolean doneFirstIteration =false;//check on index only done on first iteration to only get elements with an index greater than the starting element
		while (stack.size()>0){
			Element currentElement =stack.pop();
			if (currentElement.getLocalName().equals(GROUP_EL)){
				Fragment groupFrag =state.xmlFragmentMap.get(currentElement);
				if (groupFrag.hasLocant(locant)){
					return groupFrag;
				}
				continue;
			}
			List<Element> siblings = XOMTools.getChildElementsWithTagNames(currentElement, new String[]{BRACKET_EL, SUBSTITUENT_EL, ROOT_EL});

			Stack<Element> bracketted = new Stack<Element>();
			for (Element bracketOrSubOrRoot : siblings) {
				if (!doneFirstIteration && currentElement.indexOf(bracketOrSubOrRoot)<=currentElement.indexOf(startingElement)){
					continue;
				}
				if (bracketOrSubOrRoot.getAttribute(MULTIPLIER_ATR)!=null){
					continue;
				}
				if (bracketOrSubOrRoot.getLocalName().equals(BRACKET_EL)){
					if (IMPLICIT_TYPE_VAL.equals(bracketOrSubOrRoot.getAttributeValue(TYPE_ATR))){
						stack.add(bracketOrSubOrRoot);
					}
					else{
						bracketted.add(bracketOrSubOrRoot);
					}
				}
				else{
					Element group = bracketOrSubOrRoot.getFirstChildElement(GROUP_EL);
					stack.add(group);
				}
			}
			stack.addAll(0, bracketted);//locanting into brackets is rarely the desired answer so place at the bottom of the stack
			doneFirstIteration =true;
		}
		return null;
	}

	static Element findRightMostGroupInBracket(Element bracket) {
		List<Element> subsBracketsAndRoots = XOMTools.getChildElementsWithTagNames(bracket, new String[]{BRACKET_EL, SUBSTITUENT_EL, ROOT_EL});
		while (subsBracketsAndRoots.get(subsBracketsAndRoots.size()-1).getLocalName().equals(BRACKET_EL)){
			subsBracketsAndRoots = XOMTools.getChildElementsWithTagNames(subsBracketsAndRoots.get(subsBracketsAndRoots.size()-1), new String[]{BRACKET_EL, SUBSTITUENT_EL, ROOT_EL});
		}
		return subsBracketsAndRoots.get(subsBracketsAndRoots.size()-1).getFirstChildElement(GROUP_EL);
	}

	private static boolean potentiallyCanSubstitute(Element subBracketOrRoot) {
		Element parent =(Element) subBracketOrRoot.getParent();
		Elements children =parent.getChildElements();
		for (int i = parent.indexOf(subBracketOrRoot) +1 ; i < children.size(); i++) {
			if (!children.get(i).getLocalName().equals(HYPHEN_EL)){
				return true;
			}
		}
		return false;
	}

	/**
	 * In cases such as methylenecyclohexane two outAtoms are combined to form a single outAtom with valency
	 * equal to sum of the valency of the other outAtoms.
	 * This is only allowed on substituents where all the outAtoms are on the same atom
	 * @param frag
	 * @param subType 
	 * @throws StructureBuildingException
	 */
	private static void checkAndApplySpecialCaseWhereOutAtomsCanBeCombinedOrThrow(Fragment frag, String subType) throws StructureBuildingException {
		int outAtomCount = frag.getOutAtoms().size();
		if (outAtomCount<=1){
			return;
		}
		if (EPOXYLIKE_SUBTYPE_VAL.equals(subType)){
			return;
		}
		//special case- all outAtoms on same atom e.g. methylenecyclohexane
		Atom firstOutAtom = frag.getOutAtom(0).getAtom();
		int valencyOfOutAtom =0;
		for (int i = outAtomCount -1; i >=0 ; i--) {//remove all outAtoms and add one with the total valency of all those that have been removed
			OutAtom out = frag.getOutAtom(i);
			if (out.getAtom() !=firstOutAtom){
				throw new StructureBuildingException("Substitutive bond formation failure: Fragment expected to have one OutAtom but had: "+ outAtomCount);
			}
			valencyOfOutAtom +=out.getValency();
			frag.removeOutAtom(i);
		}
		frag.addOutAtom(frag.getFirstAtom(), valencyOfOutAtom, true);
	}

	/**
	 * Calculates the number of substitutable hydrogen by taking into account:
	 * Specified valency if applicable, outAtoms and the lowest valency state that will satisfy these
	 * e.g. thio has 2 outAtoms and no bonds hence -->2 outgoing, lowest stable valency = 2 hence no substitutable hydrogen
	 * e.g. phosphonyl has 2 outAtoms and one double bond -->4 outgoing, lowest stable valency =5 hence 1 substitutable hydrogen
	 * @param atom
	 * @return
	 */
	static int calculateSubstitutableHydrogenAtoms(Atom atom) {
		int valency = atom.determineValency(true);	
		int currentValency =atom.getIncomingValency() + atom.getOutValency();
		return valency-currentValency;
	}

	/**
	 * Stereochemistry terms are assigned right at the end so that checks can be done on whether the indicated atom is in fact chiral.
	 * In the process of multiplication locants are primed. This function adds the appropriate number of primes to any locanted stereochemistry locants
	 * The primesString is the string containing the primes to add to each locant
	 * @param subOrBracket
	 * @param primesString
	 */
	private static void addPrimesToLocantedStereochemistryElements(Element subOrBracket, String primesString) {
		List<Element> stereoChemistryElements =XOMTools.getDescendantElementsWithTagName(subOrBracket, STEREOCHEMISTRY_EL);
		for (Element stereoChemistryElement : stereoChemistryElements) {
			if (!getLocant(stereoChemistryElement).equals("0")){
				stereoChemistryElement.getAttribute(LOCANT_ATR).setValue(getLocant(stereoChemistryElement) + primesString);
			}
		}
	}

	/**
	 * Calculates the number of times getParent() must be called to reach a word element
	 * Returns null if element does not have an enclosing word element.
	 * @param element
	 * @return
	 */
	private static Integer levelsToWordEl(Element element) {
		int count =0;
		while (!element.getLocalName().equals(WORD_EL)){
			element =(Element) element.getParent();
			if (element==null){
				return null;
			}
			count++;
		}	
		return count;
	}

	/**Gets the locant from a group/suffix tag, defaulting to "0"
	 *
	 * @param element
	 * @return The locant on the group/suffix tag.
	 */
	static String getLocant(Element element) {
		String locantStr = element.getAttributeValue(LOCANT_ATR);
		if(locantStr == null) return "0";
		return locantStr;
	}
}
