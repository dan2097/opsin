package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Deque;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.log4j.Logger;

import uk.ac.cam.ch.wwmm.opsin.IsotopeSpecificationParser.IsotopeSpecification;

import static uk.ac.cam.ch.wwmm.opsin.XmlDeclarations.*;
import static uk.ac.cam.ch.wwmm.opsin.OpsinTools.*;

/**
 * Methods for processing the substitutive and additive operations that connect all the fragments together
 * as well as indicated hydrogen/unsaturation/heteroatom replacement
 * @author dl387
 *
 */
class StructureBuildingMethods {
	private static final Logger LOG = Logger.getLogger(StructureBuildingMethods.class);
	private static final Pattern matchCompoundLocant =Pattern.compile("[\\[\\(\\{](\\d+[a-z]?'*)[\\]\\)\\}]");

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
		if (word.getName().equals(WORDRULE_EL)){//already been resolved
			return;
		}
		if (!word.getName().equals(WORD_EL) && !word.getName().equals(BRACKET_EL)){
			throw new StructureBuildingException("A word or bracket is the expected input");
		}
		recursivelyResolveLocantedFeatures(state, word);
		recursivelyResolveUnLocantedFeatures(state, word);
		//TODO check all things that can substitute have outAtoms
		//TOOD think whether you can avoid the need to have a cansubstitute function by only using appropriate group
		List<Element> subsBracketsAndRoots = OpsinTools.getDescendantElementsWithTagNames(word, new String[]{BRACKET_EL, SUBSTITUENT_EL, ROOT_EL});
		for (Element subsBracketsAndRoot : subsBracketsAndRoots) {
			if (subsBracketsAndRoot.getAttribute(MULTIPLIER_ATR) != null) {
				throw new StructureBuildingException("Structure building problem: multiplier on :" + subsBracketsAndRoot.getName() + " was never used");
			}
		}
		List<Element> groups = OpsinTools.getDescendantElementsWithTagName(word, GROUP_EL);
		for (int i = 0; i < groups.size(); i++) {
			Element group = groups.get(i);
			if (group.getAttribute(RESOLVED_ATR)==null && i != groups.size()-1){
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
		if (!word.getName().equals(WORD_EL) && !word.getName().equals(BRACKET_EL)){
			throw new StructureBuildingException("A word or bracket is the expected input");
		}
		List<Element> subsBracketsAndRoots = OpsinTools.getChildElementsWithTagNames(word, new String[]{BRACKET_EL, SUBSTITUENT_EL, ROOT_EL});
		//substitution occurs left to right so by doing this right to left you ensure that any groups that will come into existence
		//due to multipliers being expanded will be in existence
		for (int i =subsBracketsAndRoots.size()-1; i>=0; i--) {
			Element subBracketOrRoot = subsBracketsAndRoots.get(i);
			if (subBracketOrRoot.getName().equals(BRACKET_EL)){
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
		if (!word.getName().equals(WORD_EL) && !word.getName().equals(BRACKET_EL)){
			throw new StructureBuildingException("A word or bracket is the expected input");
		}
		List<Element> subsBracketsAndRoots = OpsinTools.getChildElementsWithTagNames(word, new String[]{BRACKET_EL, SUBSTITUENT_EL, ROOT_EL});
		//substitution occurs left to right so by doing this right to left you ensure that any groups that will come into existence
		//due to multipliers being expanded will be in existence
		for (int i =subsBracketsAndRoots.size()-1; i>=0; i--) {
			Element subBracketOrRoot = subsBracketsAndRoots.get(i);
			if (subBracketOrRoot.getName().equals(BRACKET_EL)){
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
		if (subBracketOrRoot.getName().equals(BRACKET_EL)) {
			group = findRightMostGroupInBracket(subBracketOrRoot);
		}
		else{
			group = subBracketOrRoot.getFirstChildElement(GROUP_EL);
		}
		if (group.getAttribute(RESOLVED_ATR) != null) {
			return;
		}
		Fragment frag = group.getFrag();
		if (frag.getOutAtomCount() >=1 && subBracketOrRoot.getAttribute(LOCANT_ATR) != null){
			String locantString = subBracketOrRoot.getAttributeValue(LOCANT_ATR);
			if (frag.getOutAtomCount() >1){
				checkAndApplySpecialCaseWhereOutAtomsCanBeCombinedOrThrow(frag, group);
			}
			if (subBracketOrRoot.getAttribute(MULTIPLIER_ATR) != null) {//e.g. 1,2-diethyl
				multiplyOutAndSubstitute(state, subBracketOrRoot);
			}
			else{
				Fragment parentFrag = findFragmentWithLocant(subBracketOrRoot, locantString);
				if (parentFrag == null){
					String modifiedLocant = checkForBracketedPrimedLocantSpecialCase(subBracketOrRoot, locantString);
					if (modifiedLocant != null){
						parentFrag = findFragmentWithLocant(subBracketOrRoot, modifiedLocant);
						if (parentFrag != null){
							locantString = modifiedLocant;
						}
					}
				}
				if (parentFrag==null){
					throw new StructureBuildingException("Cannot find in scope fragment with atom with locant " + locantString + ".");
				}
				group.addAttribute(new Attribute(RESOLVED_ATR, "yes"));
				Element groupToAttachTo = parentFrag.getTokenEl();
				if (groupToAttachTo.getAttribute(ACCEPTSADDITIVEBONDS_ATR) != null &&
						parentFrag.getOutAtomCount() > 0 &&
						groupToAttachTo.getAttribute(ISAMULTIRADICAL_ATR) != null &&
						parentFrag.getAtomByLocantOrThrow(locantString).getOutValency() > 0 &&
						frag.getOutAtom(0).getValency() == 1 &&
						parentFrag.getFirstAtom().equals(parentFrag.getAtomByLocantOrThrow(locantString))) {
					//horrible special case to allow C-hydroxycarbonimidoyl and the like
					//If additive nomenclature the first atom should be an out atom
					joinFragmentsAdditively(state, frag, parentFrag);
				}
				else{
					Atom atomToSubstituteAt = parentFrag.getAtomByLocantOrThrow(locantString);
					if (PHOSPHO_SUBTYPE_VAL.equals(group.getAttributeValue(SUBTYPE_ATR)) && frag.getOutAtom(0).getValency() == 1){
						if (atomToSubstituteAt.getElement() != ChemEl.O){
							for (Atom neighbour : atomToSubstituteAt.getAtomNeighbours()) {
								if (neighbour.getElement() == ChemEl.O &&
										neighbour.getBondCount()==1 &&
										neighbour.getFirstBond().getOrder() == 1 &&
										neighbour.getOutValency() == 0 &&
										neighbour.getCharge() == 0){
									atomToSubstituteAt = neighbour;
									break;
								}
							}
						}
					}
					joinFragmentsSubstitutively(state, frag, atomToSubstituteAt);
				}
			}
		}
	}

	private static void performUnLocantedSubstitutiveOperations(BuildState state, Element subBracketOrRoot) throws StructureBuildingException {
		Element group;
		if (subBracketOrRoot.getName().equals(BRACKET_EL)){
			group = findRightMostGroupInBracket(subBracketOrRoot);
		}
		else{
			group = subBracketOrRoot.getFirstChildElement(GROUP_EL);
		}
		if (group.getAttribute(RESOLVED_ATR) != null){
			return;
		}
		Fragment frag = group.getFrag();
		if (frag.getOutAtomCount() >= 1){
			if (subBracketOrRoot.getAttribute(LOCANT_ATR) != null){
				throw new RuntimeException("Substituent has an unused outAtom and has a locant but locanted substitution should already have been performed!");
			}
			if (frag.getOutAtomCount() > 1){
				checkAndApplySpecialCaseWhereOutAtomsCanBeCombinedOrThrow(frag, group);
			}
			if (subBracketOrRoot.getAttribute(MULTIPLIER_ATR) != null) {//e.g. diethyl
				multiplyOutAndSubstitute(state, subBracketOrRoot);
			}
			else{
				if (PERHALOGENO_SUBTYPE_VAL.equals(group.getAttributeValue(SUBTYPE_ATR))) {
					performPerHalogenoSubstitution(state, frag, subBracketOrRoot);
				}
				else{
					List<Atom> atomsToJoinTo = null;
					if (PHOSPHO_SUBTYPE_VAL.equals(group.getAttributeValue(SUBTYPE_ATR)) && frag.getOutAtom(0).getValency() == 1){
						List<Fragment> possibleParents = findAlternativeFragments(subBracketOrRoot);
						for (Fragment fragment : possibleParents) {
							List<Atom> hydroxyAtoms = FragmentTools.findHydroxyGroups(fragment);
							if (hydroxyAtoms.size() >= 1){
								atomsToJoinTo = hydroxyAtoms;
							}
							break;
						}
					}
					if (atomsToJoinTo == null) {
						atomsToJoinTo = findAtomsForSubstitution(subBracketOrRoot, 1, frag.getOutAtom(0).getValency());
					}
					if (atomsToJoinTo == null){
						throw new StructureBuildingException("Unlocanted substitution failed: unable to find suitable atom to bond atom with id:" + frag.getOutAtom(0).getAtom().getID() + " to!");
					}
					if (AmbiguityChecker.isSubstitutionAmbiguous(atomsToJoinTo, 1)) {
						state.addIsAmbiguous("Connection of " + group.getValue() + " to " + atomsToJoinTo.get(0).getFrag().getTokenEl().getValue());
					}
					joinFragmentsSubstitutively(state, frag, atomsToJoinTo.get(0));
				}
				group.addAttribute(new Attribute(RESOLVED_ATR, "yes"));
			}
		}
	}

	/**
	 * Clones the perhalogenFrag sufficiently to replace all in scope hydrogen with halogens.
	 * The cloned fragments are merged into the perhalogenFrag
	 * @param state
	 * @param perhalogenFrag
	 * @param subBracketOrRoot
	 * @throws StructureBuildingException
	 */
	private static void performPerHalogenoSubstitution(BuildState state, Fragment perhalogenFrag, Element subBracketOrRoot) throws StructureBuildingException {
		List<Fragment> fragmentsToAttachTo = findAlternativeFragments(subBracketOrRoot);
		List<Atom> atomsToHalogenate = new ArrayList<Atom>();
		for (Fragment fragment : fragmentsToAttachTo) {
			FragmentTools.convertSpareValenciesToDoubleBonds(fragment);
			for (Atom atom : fragment.getAtomList()) {
				int substitutableHydrogen = calculateSubstitutableHydrogenAtoms(atom);
				if (substitutableHydrogen > 0 && FragmentTools.isCharacteristicAtom(atom)){
					continue;
				}
				for (int i = 0; i < substitutableHydrogen; i++) {
					atomsToHalogenate.add(atom);
				}
			}
		}
		if (atomsToHalogenate.size() == 0){
			throw new RuntimeException("Failed to find any substitutable hydrogen to apply " + perhalogenFrag.getTokenEl().getValue() + " to!");
		}
		List<Fragment> halogens = new ArrayList<Fragment>();
		halogens.add(perhalogenFrag);
		for (int i = 0; i < atomsToHalogenate.size() - 1; i++) {
			halogens.add(state.fragManager.copyFragment(perhalogenFrag));
		}
		for (int i = 0; i < atomsToHalogenate.size(); i++) {
			Fragment halogen = halogens.get(i);
			Atom from = halogen.getOutAtom(0).getAtom();
			halogen.removeOutAtom(0);
			state.fragManager.createBond(from, atomsToHalogenate.get(i), 1);
		}
		for (int i = 1; i < atomsToHalogenate.size(); i++) {
			state.fragManager.incorporateFragment(halogens.get(i), perhalogenFrag);
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
		Attribute multiplierAtr = subOrBracket.getAttribute(MULTIPLIER_ATR);
		int multiplier = Integer.parseInt(multiplierAtr.getValue());
		subOrBracket.removeAttribute(multiplierAtr);
		String[] locants = null;
		String locantsAtrValue = subOrBracket.getAttributeValue(LOCANT_ATR);
		if (locantsAtrValue != null){
			locants = locantsAtrValue.split(",");
		}
		Element parentWordOrBracket = subOrBracket.getParent();
		int indexOfSubOrBracket = parentWordOrBracket.indexOf(subOrBracket);
		subOrBracket.detach();

		List<Element> elementsNotToBeMultiplied = new ArrayList<Element>();//anything before the multiplier in the sub/bracket
		Element multiplierEl = subOrBracket.getFirstChildElement(MULTIPLIER_EL);
		if (multiplierEl == null){
			throw new RuntimeException("Multiplier not found where multiplier expected");
		}
		for (int i = subOrBracket.indexOf(multiplierEl) -1 ; i >=0 ; i--) {
			Element el = subOrBracket.getChild(i);
			el.detach();
			elementsNotToBeMultiplied.add(el);
		}
		multiplierEl.detach();

		List<Element> multipliedElements = new ArrayList<Element>();
		for (int i = multiplier - 1; i >=0; i--) {
			Element currentElement;
			if (i != 0){
				currentElement = state.fragManager.cloneElement(state, subOrBracket, i);
				addPrimesToLocantedStereochemistryElements(currentElement, StringTools.multiplyString("'", i));//Stereochemistry elements with locants will need to have their locants primed (stereochemistry is only processed after structure building)
			}
			else{
				currentElement = subOrBracket;
			}
			multipliedElements.add(currentElement);
			if (locants != null){
				parentWordOrBracket.insertChild(currentElement, indexOfSubOrBracket);
				currentElement.getAttribute(LOCANT_ATR).setValue(locants[i]);
				performLocantedSubstitutiveOperations(state, currentElement);
				currentElement.detach();
			}
		}
		if (locants == null) {
			parentWordOrBracket.insertChild(multipliedElements.get(0), indexOfSubOrBracket);
			performUnlocantedSubstitutiveOperations(state, multipliedElements);
			multipliedElements.get(0).detach();
		}
		for (Element multipliedElement : multipliedElements) {//attach all the multiplied subs/brackets
			parentWordOrBracket.insertChild(multipliedElement, indexOfSubOrBracket);
		}
		for (Element el : elementsNotToBeMultiplied) {//re-add anything before multiplier to original subOrBracket
			subOrBracket.insertChild(el, 0);
		}
	}

	private static void performUnlocantedSubstitutiveOperations(BuildState state, List<Element> multipliedElements) throws StructureBuildingException {
		int numOfSubstituents = multipliedElements.size();
		Element subBracketOrRoot = multipliedElements.get(0);
		Element group;
		if (subBracketOrRoot.getName().equals(BRACKET_EL)){
			group = findRightMostGroupInBracket(subBracketOrRoot);
		}
		else{
			group = subBracketOrRoot.getFirstChildElement(GROUP_EL);
		}
		Fragment frag = group.getFrag();
		if (frag.getOutAtomCount() >= 1){
			if (subBracketOrRoot.getAttribute(LOCANT_ATR) != null){
				throw new RuntimeException("Substituent has an unused outAtom and has a locant but locanted substitution should already been been performed!");
			}
			if (PERHALOGENO_SUBTYPE_VAL.equals(group.getAttributeValue(SUBTYPE_ATR))) {
				throw new StructureBuildingException(group.getValue() + " cannot be multiplied");
			}
			if (frag.getOutAtomCount() > 1){
				checkAndApplySpecialCaseWhereOutAtomsCanBeCombinedOrThrow(frag, group);
			}
			List<Atom> atomsToJoinTo = null;
			if (PHOSPHO_SUBTYPE_VAL.equals(group.getAttributeValue(SUBTYPE_ATR)) && frag.getOutAtom(0).getValency() == 1){
				List<Fragment> possibleParents = findAlternativeFragments(subBracketOrRoot);
				for (Fragment fragment : possibleParents) {
					List<Atom> hydroxyAtoms = FragmentTools.findHydroxyGroups(fragment);
					if (hydroxyAtoms.size() >= numOfSubstituents){
						atomsToJoinTo = hydroxyAtoms;
					}
					break;
				}
			}
			if (atomsToJoinTo == null) {
				atomsToJoinTo = findAtomsForSubstitution(subBracketOrRoot, numOfSubstituents, frag.getOutAtom(0).getValency());
			}
			if (atomsToJoinTo == null) {
				throw new StructureBuildingException("Unlocanted substitution failed: unable to find suitable atom to bond atom with id:" + frag.getOutAtom(0).getAtom().getID() + " to!");
			}
			if (AmbiguityChecker.isSubstitutionAmbiguous(atomsToJoinTo, numOfSubstituents)) {
				state.addIsAmbiguous("Connection of " + group.getValue() + " to " + atomsToJoinTo.get(0).getFrag().getTokenEl().getValue());
				List<Atom> atomsPreferredByEnvironment = AmbiguityChecker.useAtomEnvironmentsToGivePlausibleSubstitution(atomsToJoinTo, numOfSubstituents);
				if (atomsPreferredByEnvironment != null) {
					atomsToJoinTo = atomsPreferredByEnvironment;
				}
			}

			joinFragmentsSubstitutively(state, frag, atomsToJoinTo.get(0));
			group.addAttribute(new Attribute(RESOLVED_ATR, "yes"));
			
			for (int i = 1; i < numOfSubstituents; i++) {
				subBracketOrRoot = multipliedElements.get(i);
				if (subBracketOrRoot.getName().equals(BRACKET_EL)){
					group = findRightMostGroupInBracket(subBracketOrRoot);
				}
				else{
					group = subBracketOrRoot.getFirstChildElement(GROUP_EL);
				}
				frag = group.getFrag();
				if (frag.getOutAtomCount() > 1){//TODO do this prior to multiplication?
					checkAndApplySpecialCaseWhereOutAtomsCanBeCombinedOrThrow(frag, group);
				}
				
				joinFragmentsSubstitutively(state, frag, atomsToJoinTo.get(i));
				group.addAttribute(new Attribute(RESOLVED_ATR, "yes"));
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
		List<Element> groups = subOrRoot.getChildElements(GROUP_EL);
		if (groups.size() != 1){
			throw new StructureBuildingException("Each sub or root should only have one group element. This indicates a bug in OPSIN");
		}
		Element group = groups.get(0);
		Fragment thisFrag = group.getFrag();

		List<Element> unsaturators = new ArrayList<Element>();
		List<Element> heteroatoms = new ArrayList<Element>();
		List<Element> hydrogenElements = new ArrayList<Element>();
		List<Element> subtractivePrefixElements = new ArrayList<Element>();
		List<Element> isotopeSpecifications = new ArrayList<Element>();

		List<Element> children =subOrRoot.getChildElements();
		for (Element currentEl : children) {
			String elName =currentEl.getName();
			if (elName.equals(UNSATURATOR_EL)){
				unsaturators.add(currentEl);
			}
			else if (elName.equals(HETEROATOM_EL)){
				heteroatoms.add(currentEl);
			}
			else if (elName.equals(SUBTRACTIVEPREFIX_EL)){
				subtractivePrefixElements.add(currentEl);
			}
			else if (elName.equals(HYDRO_EL)){
				hydrogenElements.add(currentEl);
			}
			else if (elName.equals(INDICATEDHYDROGEN_EL)){
				hydrogenElements.add(currentEl);
			}
			else if (elName.equals(ADDEDHYDROGEN_EL)){
				hydrogenElements.add(currentEl);
			}
			else if (elName.equals(ISOTOPESPECIFICATION_EL)){
				isotopeSpecifications.add(currentEl);
			}
		}
		/*
		 * Add locanted functionality
		 */

		List<Atom> atomsToDehydro = new ArrayList<Atom>();
		//locanted substitution can be assumed to be irrelevant to subtractive operations hence perform all subtractive operations now
		Map<ChemEl, Integer> unlocantedSubtractivePrefixes = new HashMap<ChemEl, Integer>();

		for(int i = subtractivePrefixElements.size() -1; i >= 0; i--) {
			Element subtractivePrefix = subtractivePrefixElements.get(i);
			String type = subtractivePrefix.getAttributeValue(TYPE_ATR);
			if (type.equals(DEOXY_TYPE_VAL)){
				String locant = subtractivePrefix.getAttributeValue(LOCANT_ATR);
				ChemEl chemEl = ChemEl.valueOf(subtractivePrefix.getAttributeValue(VALUE_ATR));
				if (locant == null) {
					Integer count = unlocantedSubtractivePrefixes.get(chemEl);
					unlocantedSubtractivePrefixes.put(chemEl, count != null ? count + 1 : 1);
				}
				else {
					applySubtractivePrefix(state, thisFrag, chemEl, locant);
				}
			}
			else if (type.equals(ANHYDRO_TYPE_VAL)){
				applyAnhydroPrefix(state, thisFrag, subtractivePrefix);
			}
			else if (type.equals(DEHYDRO_TYPE_VAL)){
				String locant = subtractivePrefix.getAttributeValue(LOCANT_ATR);
				if(locant != null) {
					atomsToDehydro.add(thisFrag.getAtomByLocantOrThrow(locant));
				}
				else{
					throw new StructureBuildingException("locants are assumed to be required for the use of dehydro to be unambiguous");
				}
			}
			else{
				throw new StructureBuildingException("OPSIN bug: Unexpected subtractive prefix type: " + type);
			}
			subtractivePrefix.detach();
		}
		for (Entry<ChemEl, Integer> entry : unlocantedSubtractivePrefixes.entrySet()) {
			applyUnlocantedSubtractivePrefixes(state, thisFrag, entry.getKey(), entry.getValue());
		}
		
		if (atomsToDehydro.size() > 0){
			boolean isCarbohydrateDehydro = false;
			if (group.getAttributeValue(TYPE_ATR).equals(CARBOHYDRATE_TYPE_VAL)){
				Set<Atom> uniquifiedDehydroAtoms = new HashSet<Atom>(atomsToDehydro);
				if (uniquifiedDehydroAtoms.size()==atomsToDehydro.size()){//need to rule out case where dehydro is being used to form triple bonds on carbohydrates
					isCarbohydrateDehydro = true;
				}
			}
			if (isCarbohydrateDehydro){
				for (Atom a : atomsToDehydro) {
					List<Atom> hydroxyAtoms = FragmentTools.findHydroxyLikeTerminalAtoms(a.getAtomNeighbours(), ChemEl.O);
					if (hydroxyAtoms.size() > 0){
						hydroxyAtoms.get(0).getFirstBond().setOrder(2);
					}
					else{
						throw new StructureBuildingException("atom with locant " + a.getFirstLocant() + " did not have a hydroxy group to convert to a ketose");
					}
				}
			}
			else{
				List<Atom> atomsToFormDoubleBonds = new ArrayList<Atom>();
				List<Atom> atomsToFormTripleBondsBetween = new ArrayList<Atom>();//dehydro on a double/aromatic bond forms a triple bond
				
				for (Atom a : atomsToDehydro) {
					if (!a.hasSpareValency()){
						a.setSpareValency(true);
						atomsToFormDoubleBonds.add(a);
					}
					else{
						atomsToFormTripleBondsBetween.add(a);
					}
				}
				
				for (Atom atom : atomsToFormDoubleBonds) {//check that all the dehydro-ed atoms are next to another atom with spare valency
					boolean hasSpareValency =false;
					for (Atom neighbour : atom.getAtomNeighbours()) {
						if (neighbour.hasSpareValency()){
							hasSpareValency = true;
							break;
						}
					}
					if (!hasSpareValency){
						throw new StructureBuildingException("Unexpected use of dehydro; two adjacent atoms were not unsaturated such as to form a double bond");
					}
				}
				addDehydroInducedTripleBonds(atomsToFormTripleBondsBetween);
			}
		}
		
		for(int i=hydrogenElements.size() -1;i >= 0;i--) {
			Element hydrogen = hydrogenElements.get(i);
			String locant = hydrogen.getAttributeValue(LOCANT_ATR);
			if(locant != null) {
				Atom a =thisFrag.getAtomByLocantOrThrow(locant);
				if (a.hasSpareValency()){
					a.setSpareValency(false);
				}
				else{
					if (!acdNameSpiroIndicatedHydrogenBug(group, locant)){
						throw new StructureBuildingException("hydrogen addition at locant: " + locant +" was requested, but this atom is not unsaturated");
					}
				}
				hydrogenElements.remove(i);
				hydrogen.detach();
			}
		}

		for(int i=unsaturators.size() -1;i >= 0;i--) {
			Element unsaturator = unsaturators.get(i);
			String locant = unsaturator.getAttributeValue(LOCANT_ATR);
			int bondOrder = Integer.parseInt(unsaturator.getAttributeValue(VALUE_ATR));
			if(bondOrder <= 1) {
				unsaturator.detach();
				continue;
			}
			if(locant != null) {
				unsaturators.remove(unsaturator);
				/*
				 * Is the locant a compound locant e.g. 1(6) 
				 * This would indicate unsaturation between the atoms with locants 1 and 6
				 */
				Matcher matcher = matchCompoundLocant.matcher(locant);
				if (matcher.find()) {
					String compoundLocant = matcher.group(1);
					locant = matcher.replaceAll("");
					FragmentTools.unsaturate(thisFrag.getAtomByLocantOrThrow(locant), compoundLocant, bondOrder, thisFrag);
				}
				else {
					FragmentTools.unsaturate(thisFrag.getAtomByLocantOrThrow(locant), bondOrder, thisFrag);
				}
				unsaturator.detach();
			}
		}

		for(int i=heteroatoms.size() -1;i >= 0;i--) {
			Element heteroatomEl = heteroatoms.get(i);
			String locant = heteroatomEl.getAttributeValue(LOCANT_ATR);
			if(locant != null) {
				Atom heteroatom = state.fragManager.getHeteroatom(heteroatomEl.getAttributeValue(VALUE_ATR));
				Atom atomToBeReplaced =thisFrag.getAtomByLocantOrThrow(locant);
				if (heteroatom.getElement() == atomToBeReplaced.getElement() && heteroatom.getCharge() == atomToBeReplaced.getCharge()){
					throw new StructureBuildingException("The replacement term " +heteroatomEl.getValue() +" was used on an atom that already is a " + heteroatom.getElement());
				}
				state.fragManager.replaceAtomWithAtom(thisFrag.getAtomByLocantOrThrow(locant), heteroatom, true);
				if (heteroatomEl.getAttribute(LAMBDA_ATR) != null){
					thisFrag.getAtomByLocantOrThrow(locant).setLambdaConventionValency(Integer.parseInt(heteroatomEl.getAttributeValue(LAMBDA_ATR)));
				}
				heteroatoms.remove(heteroatomEl);
				heteroatomEl.detach();
			}
		}
		
		if (isotopeSpecifications.size() > 0) {
			applyIsotopeSpecifications(state, thisFrag, isotopeSpecifications, true);
		}
	}

	/**
	 * ACD/Name has a known bug where it produces names in which a suffixed saturated ring in a polycyclic spiro 
	 * is treated as if it is unsaturated and hence has indicated hydrogens
	 * e.g. 1',3'-dihydro-2H,5H-spiro[imidazolidine-4,2'-indene]-2,5-dione 
	 * @param group
	 * @param indicatedHydrogenLocant
	 * @return
	 */
	private static boolean acdNameSpiroIndicatedHydrogenBug(Element group, String indicatedHydrogenLocant) {
		if (group.getValue().startsWith("spiro")) {
			for (Element suffix : group.getParent().getChildElements(SUFFIX_EL)) {
				String suffixLocant = suffix.getAttributeValue(LOCANT_ATR);
				if (suffixLocant != null && suffixLocant.equals(indicatedHydrogenLocant)) {
					LOG.debug("Indicated hydrogen at " + indicatedHydrogenLocant + " ignored. Known bug in generated IUPAC name");
					return true;
				}
			}
		}
		return false;
	}
	
	/**
	 * Removes a terminal atom of a particular element e.g. oxygen
	 * The locant specifies the atom adjacent to the atom to be removed
	 * Formally the atom is replaced by hydrogen, hence stereochemistry is intentionally preserved
	 * @param state 
	 * @param fragment
	 * @param chemEl
	 * @param locant A locant or null
	 * @throws StructureBuildingException 
	 */
	static void applySubtractivePrefix(BuildState state, Fragment fragment, ChemEl chemEl, String locant) throws StructureBuildingException {
		Atom adjacentAtom = fragment.getAtomByLocantOrThrow(locant);
		List<Atom> applicableTerminalAtoms = FragmentTools.findHydroxyLikeTerminalAtoms(adjacentAtom.getAtomNeighbours(), chemEl);
		if (applicableTerminalAtoms.isEmpty()) {
			throw new StructureBuildingException("Unable to find terminal atom of type: " + chemEl + " at locant "+ locant +" for subtractive nomenclature");
		}
		Atom atomToRemove = applicableTerminalAtoms.get(0);
		if (FragmentTools.isFunctionalAtom(atomToRemove)) {//This can occur with aminoglycosides where the anomeric OH is removed by deoxy
			for (int i = 0, len = fragment.getFunctionalAtomCount(); i < len; i++) {
				if (atomToRemove.equals(fragment.getFunctionalAtom(i).getAtom())) {
					fragment.removeFunctionalAtom(i);
					break;
				}
			}
			fragment.addFunctionalAtom(atomToRemove.getFirstBond().getOtherAtom(atomToRemove));
		}
		FragmentTools.removeTerminalAtom(state, atomToRemove);
	}
	
	/**
	 * Removes terminal atoms of a particular element e.g. oxygen
	 * The number to remove is decided by the count
	 * Formally the atom is replaced by hydrogen, hence stereochemistry is intentionally preserved
	 * @param state
	 * @param fragment
	 * @param chemEl
	 * @param count
	 * @throws StructureBuildingException
	 */
	static void applyUnlocantedSubtractivePrefixes(BuildState state, Fragment fragment, ChemEl chemEl, int count) throws StructureBuildingException {
		List<Atom> applicableTerminalAtoms = FragmentTools.findHydroxyLikeTerminalAtoms(fragment.getAtomList(), chemEl);
		if (applicableTerminalAtoms.isEmpty() || applicableTerminalAtoms.size() < count) {
			throw new StructureBuildingException("Unable to find terminal atom of type: " + chemEl + " for subtractive nomenclature");
		}
		if (AmbiguityChecker.isSubstitutionAmbiguous(applicableTerminalAtoms, count)) {
			state.addIsAmbiguous("Group to remove with subtractive prefix");
		}
		for (int i = 0; i < count; i++) {
			Atom atomToRemove = applicableTerminalAtoms.get(i);
			if (FragmentTools.isFunctionalAtom(atomToRemove)) {//This can occur with aminoglycosides where the anomeric OH is removed by deoxy
				for (int j = 0, len = fragment.getFunctionalAtomCount(); j < len; j++) {
					if (atomToRemove.equals(fragment.getFunctionalAtom(j).getAtom())) {
						fragment.removeFunctionalAtom(j);
						break;
					}
				}
				fragment.addFunctionalAtom(atomToRemove.getFirstBond().getOtherAtom(atomToRemove));
			}
			FragmentTools.removeTerminalAtom(state, atomToRemove);
		}
	}

	private static void applyAnhydroPrefix(BuildState state, Fragment frag, Element subtractivePrefix) throws StructureBuildingException {
		ChemEl chemEl = ChemEl.valueOf(subtractivePrefix.getAttributeValue(VALUE_ATR));
		String locantStr = subtractivePrefix.getAttributeValue(LOCANT_ATR);
		if (locantStr == null) {
			throw new StructureBuildingException("Two locants are required before an anhydro prefix");
		}
		String[] locants = locantStr.split(",");
		Atom backBoneAtom1 = frag.getAtomByLocantOrThrow(locants[0]);
		Atom backBoneAtom2 = frag.getAtomByLocantOrThrow(locants[1]);
		List<Atom> applicableTerminalAtoms = FragmentTools.findHydroxyLikeTerminalAtoms(backBoneAtom1.getAtomNeighbours(), chemEl);
		if (applicableTerminalAtoms.isEmpty()){
			throw new StructureBuildingException("Unable to find terminal atom of type: " + chemEl + " for subtractive nomenclature");
		}
		FragmentTools.removeTerminalAtom(state, applicableTerminalAtoms.get(0));
		
		applicableTerminalAtoms = FragmentTools.findHydroxyLikeTerminalAtoms(backBoneAtom2.getAtomNeighbours(), chemEl);
		if (applicableTerminalAtoms.isEmpty()){
			throw new StructureBuildingException("Unable to find terminal atom of type: " + chemEl + " for subtractive nomenclature");
		}
		state.fragManager.createBond(backBoneAtom1, applicableTerminalAtoms.get(0), 1);
	}

	/**
	 * Attempts to form triple bond between the atoms in atomsToFormTripleBondsBetween
	 * Throws an exception if the list contains duplicates or atoms with no adjacent atom in the list
	 * @param atomsToFormTripleBondsBetween
	 * @throws StructureBuildingException
	 */
	private static void addDehydroInducedTripleBonds(List<Atom> atomsToFormTripleBondsBetween) throws StructureBuildingException {
		if (atomsToFormTripleBondsBetween.size()>0){
			Set<Atom> atoms = new HashSet<Atom>(atomsToFormTripleBondsBetween);
			if (atomsToFormTripleBondsBetween.size() != atoms.size()){
				throw new StructureBuildingException("locants specified for dehydro specify the same atom too many times");
			}
			atomLoop: for (int i = atomsToFormTripleBondsBetween.size()-1; i >=0; i = i-2) {//two atoms will have a triple bond formed betwen them
				Atom a = atomsToFormTripleBondsBetween.get(i);
				List<Atom> neighbours = a.getAtomNeighbours();
				for (Atom neighbour : neighbours) {
					if (atomsToFormTripleBondsBetween.contains(neighbour)){
						atomsToFormTripleBondsBetween.remove(i);
						atomsToFormTripleBondsBetween.remove(neighbour);
						Bond b = a.getBondToAtomOrThrow(neighbour);
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
		List<Element> groups = subOrRoot.getChildElements(GROUP_EL);
		if (groups.size() != 1){
			throw new StructureBuildingException("Each sub or root should only have one group element. This indicates a bug in OPSIN");
		}
		Fragment frag = groups.get(0).getFrag();

		List<Integer> unsaturationBondOrders = new ArrayList<Integer>();
		List<Element> heteroatoms = new ArrayList<Element>();
		List<Element> hydrogenElements = new ArrayList<Element>();
		List<Element> isotopeSpecifications = new ArrayList<Element>();

		List<Element> children = subOrRoot.getChildElements();
		for (Element currentEl : children) {
			String elName = currentEl.getName();
			if (elName.equals(UNSATURATOR_EL)) {
				int bondOrder = Integer.parseInt(currentEl.getAttributeValue(VALUE_ATR));
				if (bondOrder > 1) {
					unsaturationBondOrders.add(bondOrder);
				}
				currentEl.detach();
			}
			else if (elName.equals(HETEROATOM_EL)){
				heteroatoms.add(currentEl);
				currentEl.detach();
			}
			else if (elName.equals(HYDRO_EL) || 
				elName.equals(INDICATEDHYDROGEN_EL) ||
				elName.equals(ADDEDHYDROGEN_EL)){
				hydrogenElements.add(currentEl);
				currentEl.detach();
			}
			else if (elName.equals(ISOTOPESPECIFICATION_EL)){
				isotopeSpecifications.add(currentEl);
			}
		}

		if (hydrogenElements.size() > 0) {
			applyUnlocantedHydro(state, frag, hydrogenElements);
		}

		if (unsaturationBondOrders.size() > 0){
			unsaturateBonds(state, frag, unsaturationBondOrders);
		}

		if (heteroatoms.size() > 0) {
			applyUnlocantedHeteroatoms(state, frag, heteroatoms);
		}
		
		if (isotopeSpecifications.size() > 0) {
			applyIsotopeSpecifications(state, frag, isotopeSpecifications, false);
		}

		if (frag.getOutAtomCount() > 0){//assign any outAtoms that have not been set to a specific atom to a specific atom
			for (int i = 0, l = frag.getOutAtomCount(); i < l; i++) {
				OutAtom outAtom = frag.getOutAtom(i);
				if (!outAtom.isSetExplicitly()){
					outAtom.setAtom(findAtomForUnlocantedRadical(state, frag, outAtom));
					outAtom.setSetExplicitly(true);
				}
			}
		}
	}

	private static void applyUnlocantedHydro(BuildState state, Fragment frag, List<Element> hydrogenElements) throws StructureBuildingException {
		/*
		 * This function is not entirely straightforward as certain atoms definitely should have their spare valency reduced
		 * However names are not consistent as to whether they bother having the hydro tags do this!
		 * The atoms in atomsWithSV are in atom order those that can take a hydro element and then those that shouldn't really take a hydro element as its absence is unambiguous
		 */
		List<Atom> atomsAcceptingHydroPrefix = new ArrayList<Atom>();
		Set<Atom> atomsWhichImplicitlyHadTheirSVRemoved = new HashSet<Atom>();
		List<Atom> atomList = frag.getAtomList();
		for (Atom atom : atomList) {
			if (atom.getType().equals(SUFFIX_TYPE_VAL)){
				continue;
			}
			atom.ensureSVIsConsistantWithValency(false);//doesn't take into account suffixes
			if (atom.hasSpareValency()) {
				atomsAcceptingHydroPrefix.add(atom);
				//if we take into account suffixes is the SV removed
				atom.ensureSVIsConsistantWithValency(true);
				if (!atom.hasSpareValency()) {
					atomsWhichImplicitlyHadTheirSVRemoved.add(atom);
				}
			}
		}
		
		int hydrogenElsCount = hydrogenElements.size();
		for (Element hydrogenElement : hydrogenElements) {
			if (hydrogenElement.getValue().equals("perhydro")) {
				if (hydrogenElsCount != 1){
					throw new StructureBuildingException("Unexpected indication of hydrogen when perhydro makes such indication redundnant");
				}
				for (Atom atom : atomsAcceptingHydroPrefix) {
					atom.setSpareValency(false);
				}
				return;
			}
		}
		
		List<Atom> atomsWithDefiniteSV = new ArrayList<Atom>();
		List<Atom> otherAtomsThatCanHaveHydro = new ArrayList<Atom>();
		for(Atom a : atomsAcceptingHydroPrefix) {
			if (atomsWhichImplicitlyHadTheirSVRemoved.contains(a)) {
				otherAtomsThatCanHaveHydro.add(a);
			}
			else {
				boolean canFormDoubleBond = false;
				for(Atom aa : frag.getIntraFragmentAtomNeighbours(a)) {
					if(aa.hasSpareValency()) {
						canFormDoubleBond = true;
						break;
					}
				}
				if (canFormDoubleBond) {
					atomsWithDefiniteSV.add(a);
				}
				else {
					otherAtomsThatCanHaveHydro.add(a);
				}
			}
		}
		List<Atom> prioritisedAtomsAcceptingHydro = new ArrayList<Atom>(atomsWithDefiniteSV);
		prioritisedAtomsAcceptingHydro.addAll(otherAtomsThatCanHaveHydro);//these end up at the end of the list

		if (hydrogenElsCount > prioritisedAtomsAcceptingHydro.size()) {
			throw new StructureBuildingException("Cannot find atom to add hydrogen to (" +
					hydrogenElsCount + " hydrogens requested but only " + prioritisedAtomsAcceptingHydro.size() +" positions that can be hydrogenated)" );
		}
		
		int svCountAfterRemoval = atomsWithDefiniteSV.size() - hydrogenElsCount;
		if (svCountAfterRemoval > 1) { //ambiguity likely. If it's 1 then an atom will be implicitly hydrogenated
			//NOTE: as hydrogens as added in pairs the unambiguous if one hydrogen is added and allow atoms are identical condition is unlikely to be ever satisfied
			if (!(AmbiguityChecker.allAtomsEquivalent(atomsWithDefiniteSV) &&
					(hydrogenElsCount == 1 || hydrogenElsCount == atomsWithDefiniteSV.size() - 1))) {
				state.addIsAmbiguous("Ambiguous choice of positions to add hydrogen to on " + frag.getTokenEl().getValue());
			}
		}

		for (int i = 0; i < hydrogenElsCount; i++) {
			prioritisedAtomsAcceptingHydro.get(i).setSpareValency(false);
		}
	}

	private static void unsaturateBonds(BuildState state, Fragment frag, List<Integer> unsaturationBondOrders) throws StructureBuildingException {
		int tripleBonds = 0;
		int doublebonds = 0;
		for (Integer bondOrder : unsaturationBondOrders) {
			if (bondOrder == 3) {
				tripleBonds++;
			}
			else if (bondOrder == 2) {
				doublebonds++;
			}
			else {
				throw new RuntimeException("Unexpected unsaturation bon order: " + bondOrder);
			}
		}
		
		if (tripleBonds > 0) {
			unsaturateBonds(state, frag, 3, tripleBonds);
		}
		if (doublebonds > 0) {
			unsaturateBonds(state, frag, 2, doublebonds);
		}
	}

	private static void unsaturateBonds(BuildState state, Fragment frag, int bondOrder, int numToUnsaturate) throws StructureBuildingException {
		List<Bond> bondsThatCouldBeUnsaturated = findBondsToUnSaturate(frag, bondOrder, false);
		List<Bond> alternativeBondsThatCouldBeUnsaturated = Collections.emptyList();
		if (bondsThatCouldBeUnsaturated.size() < numToUnsaturate){
			bondsThatCouldBeUnsaturated = findBondsToUnSaturate(frag, bondOrder, true);
		}
		else {
			alternativeBondsThatCouldBeUnsaturated = findAlternativeBondsToUnSaturate(frag, bondOrder, bondsThatCouldBeUnsaturated);
		}
		if (bondsThatCouldBeUnsaturated.size() < numToUnsaturate){
			throw new StructureBuildingException("Failed to find bond to change to a bond of order: " + bondOrder);
		}
		if (bondsThatCouldBeUnsaturated.size() > numToUnsaturate) {
			//by convention cycloalkanes can have one unsaturation implicitly at the 1 locant
			//terms like oxazoline are formally ambiguous but in practice the lowest locant is the one that will be intended (in this case 2-oxazoline)
			if (!isCycloAlkaneSpecialCase(frag, numToUnsaturate, bondsThatCouldBeUnsaturated) &&
					!HANTZSCHWIDMAN_SUBTYPE_VAL.equals(frag.getSubType())) {
				if (alternativeBondsThatCouldBeUnsaturated.size() >= numToUnsaturate) {
					List<Bond> allBonds = new ArrayList<Bond>(bondsThatCouldBeUnsaturated);
					allBonds.addAll(alternativeBondsThatCouldBeUnsaturated);
					if (!(AmbiguityChecker.allBondsEquivalent(allBonds) &&
							numToUnsaturate == 1 )) {
						state.addIsAmbiguous("Unsaturation of bonds of " + frag.getTokenEl().getValue());
					}
				}
				else {
					if (!(AmbiguityChecker.allBondsEquivalent(bondsThatCouldBeUnsaturated) && 
							(numToUnsaturate == 1 || numToUnsaturate == bondsThatCouldBeUnsaturated.size() - 1))){
						state.addIsAmbiguous("Unsaturation of bonds of " + frag.getTokenEl().getValue());
					}
				}
			}
		}
		for (int i = 0; i < numToUnsaturate; i++) {
			bondsThatCouldBeUnsaturated.get(i).setOrder(bondOrder);
		}
	}

	private static boolean isCycloAlkaneSpecialCase(Fragment frag, int numToUnsaturate, List<Bond> bondsThatCouldBeUnsaturated) {
		if (numToUnsaturate == 1) {
			Bond b = bondsThatCouldBeUnsaturated.get(0);
			Atom a1 = b.getFromAtom();
			Atom a2 = b.getToAtom();
			if ((ALKANESTEM_SUBTYPE_VAL.equals(frag.getSubType()) || HETEROSTEM_SUBTYPE_VAL.equals(frag.getSubType())) && 
					a1.getAtomIsInACycle() && a2.getAtomIsInACycle() &&
					(a1.equals(frag.getFirstAtom()) || a2.equals(frag.getFirstAtom()))) {
				//mono unsaturated cyclo alkanes are unambiguous e.g. cyclohexene
				return true;
			}
		}
		return false;
	}
	
	private static boolean isCycloAlkaneHeteroatomSpecialCase(Fragment frag, int numHeteroatoms, List<Atom> atomsThatCouldBeReplaced) {
		if (numHeteroatoms == 1) {
			if ((ALKANESTEM_SUBTYPE_VAL.equals(frag.getSubType()) || HETEROSTEM_SUBTYPE_VAL.equals(frag.getSubType())) && 
					frag.getFirstAtom().getAtomIsInACycle() && atomsThatCouldBeReplaced.get(0).equals(frag.getFirstAtom())) {
				//single heteroatom implicitly goes to 1 position
				return true;
			}
		}
		return false;
	}
	
	private static class HeteroAtomSmilesAndLambda {
		private final String smiles;
		private final String lambdaConvention;
		
		public HeteroAtomSmilesAndLambda(String smiles, String lambdaConvention) {
			this.smiles = smiles;
			this.lambdaConvention = lambdaConvention;
		}

		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime
					* result
					+ ((lambdaConvention == null) ? 0 : lambdaConvention
							.hashCode());
			result = prime * result
					+ ((smiles == null) ? 0 : smiles.hashCode());
			return result;
		}

		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			HeteroAtomSmilesAndLambda other = (HeteroAtomSmilesAndLambda) obj;
			if (lambdaConvention == null) {
				if (other.lambdaConvention != null)
					return false;
			} else if (!lambdaConvention.equals(other.lambdaConvention))
				return false;
			if (smiles == null) {
				if (other.smiles != null)
					return false;
			} else if (!smiles.equals(other.smiles))
				return false;
			return true;
		}
		
		
	}

	private static void applyUnlocantedHeteroatoms(BuildState state, Fragment frag, List<Element> heteroatoms) throws StructureBuildingException {
		Map<HeteroAtomSmilesAndLambda, Integer> heteroatomDescriptionToCount = new HashMap<HeteroAtomSmilesAndLambda, Integer>();
		for (Element heteroatomEl : heteroatoms) {
			String smiles = heteroatomEl.getAttributeValue(VALUE_ATR);
			String lambdaConvention = heteroatomEl.getAttributeValue(LAMBDA_ATR);
			HeteroAtomSmilesAndLambda desc = new HeteroAtomSmilesAndLambda(smiles, lambdaConvention);
			Integer count = heteroatomDescriptionToCount.get(desc);
			heteroatomDescriptionToCount.put(desc, count != null ? count + 1 : 1);
		}
		List<Atom> atomlist = frag.getAtomList();
		for (Entry<HeteroAtomSmilesAndLambda, Integer> entry : heteroatomDescriptionToCount.entrySet()) {
			HeteroAtomSmilesAndLambda desc = entry.getKey();
			int replacementsRequired = entry.getValue();
			Atom heteroatom = state.fragManager.getHeteroatom(desc.smiles);
			ChemEl heteroatomChemEl = heteroatom.getElement();
			//finds an atom for which changing it to the specified heteroatom will not cause valency to be violated
			List<Atom> atomsThatCouldBeReplaced = new ArrayList<Atom>();
			for (Atom atom : atomlist) {
				if (atom.getType().equals(SUFFIX_TYPE_VAL)) {
					continue;
				}
				if ((heteroatomChemEl.equals(atom.getElement()) && heteroatom.getCharge() == atom.getCharge())){
					continue;//replacement would do nothing
				}
				if(atom.getElement() != ChemEl.C && heteroatomChemEl != ChemEl.C){
					if (atom.getElement() == ChemEl.O && (heteroatomChemEl == ChemEl.S || heteroatomChemEl == ChemEl.Se || heteroatomChemEl == ChemEl.Te)) {
						//by special case allow replacement of oxygen by chalcogen
					}
					else{
						//replacement of heteroatom by another heteroatom
						continue;
					}
				}
				if (ValencyChecker.checkValencyAvailableForReplacementByHeteroatom(atom, heteroatom)) {
					atomsThatCouldBeReplaced.add(atom);
				}
			}
			if (atomsThatCouldBeReplaced.size() < replacementsRequired){
				throw new StructureBuildingException("Cannot find suitable atom for heteroatom replacement");
			}
			
			if (atomsThatCouldBeReplaced.size() > replacementsRequired && !isCycloAlkaneHeteroatomSpecialCase(frag, replacementsRequired, atomsThatCouldBeReplaced)) {
				if (!(AmbiguityChecker.allAtomsEquivalent(atomsThatCouldBeReplaced) &&
						(replacementsRequired == 1 || replacementsRequired == atomsThatCouldBeReplaced.size() - 1))) {
					//by convention cycloalkanes can have one unsaturation implicitly at the 1 locant
					state.addIsAmbiguous("Heteroatom replacement on " + frag.getTokenEl().getValue());
				}
			}
			
			for (int i = 0; i < replacementsRequired; i++) {
				Atom atomToReplaceWithHeteroAtom = atomsThatCouldBeReplaced.get(i);
				state.fragManager.replaceAtomWithAtom(atomToReplaceWithHeteroAtom, heteroatom, true);
				if (desc.lambdaConvention != null) {
					atomToReplaceWithHeteroAtom.setLambdaConventionValency(Integer.parseInt(desc.lambdaConvention));
				}
			}
		}
	}
	
	private static void applyIsotopeSpecifications(BuildState state, Fragment frag, List<Element> isotopeSpecifications, boolean applyLocanted) throws StructureBuildingException {
		for(int i = isotopeSpecifications.size() - 1; i >= 0; i--) {
			Element isotopeSpecification = isotopeSpecifications.get(i);
			IsotopeSpecification isotopeSpec = IsotopeSpecificationParser.parseIsotopeSpecification(isotopeSpecification);
			String[] locants = isotopeSpec.getLocants();
			if(locants != null) {
				if (!applyLocanted) {
					continue;
				}
			}
			else if (applyLocanted) {
				continue;
			}

			ChemEl chemEl = isotopeSpec.getChemEl();
			int isotope = isotopeSpec.getIsotope();
			if(locants != null) {
				if (chemEl == ChemEl.H) {
					for (int j = 0; j < locants.length; j++) {
						Atom atomWithHydrogenIsotope = frag.getAtomByLocantOrThrow(locants[j]);
						Atom hydrogen = state.fragManager.createAtom(isotopeSpec.getChemEl(), frag);
						hydrogen.setIsotope(isotope);
						state.fragManager.createBond(atomWithHydrogenIsotope, hydrogen, 1);
					}
				}
				else {
					for (int j = 0; j < locants.length; j++) {
						Atom atom = frag.getAtomByLocantOrThrow(locants[j]);
						if (chemEl != atom.getElement()) {
							throw new StructureBuildingException("The atom at locant: " + locants[j]  + " was not a " + chemEl.toString() );
						}
						atom.setIsotope(isotope);
					}
				}
			}
			else {
				int multiplier = isotopeSpec.getMultiplier();
				if (chemEl == ChemEl.H) {
					List<Atom> parentAtomsToApplyTo = FragmentTools.findnAtomsForSubstitution(frag, multiplier, 1);
					if (parentAtomsToApplyTo == null){
						throw new StructureBuildingException("Failed to find sufficient hydrogen atoms for unlocanted hydrogen isotope replacement");
					}
					if (AmbiguityChecker.isSubstitutionAmbiguous(parentAtomsToApplyTo, multiplier)) {
						if (!casIsotopeAmbiguitySpecialCase(frag, parentAtomsToApplyTo, multiplier)) {
							state.addIsAmbiguous("Position of hydrogen isotope on " + frag.getTokenEl().getValue());
						}
					}
					for (int j = 0; j < multiplier; j++) {
						Atom atomWithHydrogenIsotope = parentAtomsToApplyTo.get(j);
						Atom hydrogen = state.fragManager.createAtom(isotopeSpec.getChemEl(), frag);
						hydrogen.setIsotope(isotope);
						state.fragManager.createBond(atomWithHydrogenIsotope, hydrogen, 1);
					}
				}
				else {
					List<Atom> parentAtomsToApplyTo = new ArrayList<Atom>();
					for (Atom atom : frag.getAtomList()) {
						if (atom.getElement() == chemEl) {
							parentAtomsToApplyTo.add(atom);
						}
					}
					if (parentAtomsToApplyTo.size() < multiplier) {
						throw new StructureBuildingException("Failed to find sufficient atoms for " + chemEl.toString() + " isotope replacement");
					}
					if (AmbiguityChecker.isSubstitutionAmbiguous(parentAtomsToApplyTo, multiplier)) {
						state.addIsAmbiguous("Position of isotope on " + frag.getTokenEl().getValue());
					}
					for (int j = 0; j < multiplier; j++) {
						parentAtomsToApplyTo.get(j).setIsotope(isotope);
					}
				}
			}
			isotopeSpecification.detach();
		}
	}

	private static boolean casIsotopeAmbiguitySpecialCase(Fragment frag, List<Atom> parentAtomsToApplyTo, int multiplier) throws StructureBuildingException {
		if (multiplier !=1) {
			return false;
		}
		List<Atom> atoms = frag.getAtomList();
		Atom firstAtom = atoms.get(0);
		if (!parentAtomsToApplyTo.get(0).equals(firstAtom)) {
			return false;
		}
		ChemEl firstAtomEl = firstAtom.getElement();
		if (atoms.size() ==2) {
			if (firstAtomEl == atoms.get(1).getElement()) {
				//e.g. ethane
				return true;
			}
		}
		else {
			int intraFragValency = frag.getIntraFragmentIncomingValency(firstAtom);
			boolean spareValency = firstAtom.hasSpareValency();
			if (firstAtom.getAtomIsInACycle()) {
				for (int i = 1; i < atoms.size(); i++) {
					Atom atom = atoms.get(i);
					if (atom.getElement() != firstAtomEl){
						return false;
					}
					if (frag.getIntraFragmentIncomingValency(atom) != intraFragValency){
						return false;
					}
					if (atom.hasSpareValency() != spareValency){
						return false;
					}
				}
				//e.g. benzene
				return true;
			}
		}
		return false;
	}

	static Atom findAtomForUnlocantedRadical(BuildState state, Fragment frag, OutAtom outAtom) throws StructureBuildingException {
		List<Atom> possibleAtoms = FragmentTools.findnAtomsForSubstitution(frag, outAtom.getAtom(), 1, outAtom.getValency(), true);
		if (possibleAtoms == null){
			throw new StructureBuildingException("Failed to assign all unlocanted radicals to actual atoms without violating valency");
		}
		if (!((ALKANESTEM_SUBTYPE_VAL.equals(frag.getSubType()) || HETEROSTEM_SUBTYPE_VAL.equals(frag.getSubType())) && possibleAtoms.get(0).equals(frag.getFirstAtom()))) {
			if (AmbiguityChecker.isSubstitutionAmbiguous(possibleAtoms, 1)) {
				state.addIsAmbiguous("Positioning of radical on: " + frag.getTokenEl().getValue());
			}
		}
		return possibleAtoms.get(0);
	}
	

	private static List<Bond> findAlternativeBondsToUnSaturate(Fragment frag, int bondOrder, Collection<Bond> bondsToIgnore) {
		return findBondsToUnSaturate(frag, bondOrder, false, new HashSet<Bond>(bondsToIgnore));
	}
	
	/**
	 * Finds bond within the fragment that can have their bondOrder increased to the specified bond order
	 * Depending on the value of allowAdjacentUnsaturatedBonds adjacent higher bonds are prevented
	 * @param frag
	 * @param bondOrder
	 * @param allowAdjacentUnsaturatedBonds
	 * @return
	 */
	static List<Bond> findBondsToUnSaturate(Fragment frag, int bondOrder, boolean allowAdjacentUnsaturatedBonds) {
		return findBondsToUnSaturate(frag, bondOrder, allowAdjacentUnsaturatedBonds, Collections.<Bond>emptySet());
	}
	
	private static List<Bond> findBondsToUnSaturate(Fragment frag, int bondOrder, boolean allowAdjacentUnsaturatedBonds, Set<Bond> bondsToIgnore) {
		List<Bond> bondsToUnsaturate = new ArrayList<Bond>();
		mainLoop: for (Atom atom1 : frag.getAtomList()) {
			if (atom1.hasSpareValency() || SUFFIX_TYPE_VAL.equals(atom1.getType()) || atom1.getProperty(Atom.ISALDEHYDE) !=null) {
				continue;
			}
			List<Bond> bonds = atom1.getBonds();
			int incomingValency = 0;
			for (Bond bond : bonds) {
				//don't place implicitly unsaturated bonds next to each other
				if (bond.getOrder() != 1 && !allowAdjacentUnsaturatedBonds) {
					continue mainLoop;
				}
				if (bondsToUnsaturate.contains(bond)) {
					if (!allowAdjacentUnsaturatedBonds) {
						continue mainLoop;
					}
					incomingValency += bondOrder;
				}
				else {
					incomingValency += bond.getOrder();
				}
			}
			
			Integer maxVal = getLambdaValencyOrHwValencyOrMaxValIfCharged(atom1);
			if(maxVal != null && (incomingValency + (bondOrder - 1) + atom1.getOutValency()) > maxVal) {
				continue;
			}
			bondLoop: for (Bond bond : bonds) {
				if (bond.getOrder() == 1 && !bondsToUnsaturate.contains(bond) && !bondsToIgnore.contains(bond)) {
					Atom atom2 = bond.getOtherAtom(atom1);
					if (frag.getAtomByID(atom2.getID()) != null) {//check other atom is actually in the fragment!
						if (atom2.hasSpareValency() || SUFFIX_TYPE_VAL.equals(atom2.getType()) || atom2.getProperty(Atom.ISALDEHYDE) !=null) {
							continue;
						}
						int incomingValency2 = 0;
						for (Bond bond2 : atom2.getBonds()) {
							//don't place implicitly unsaturated bonds next to each other
							if (bond2.getOrder() != 1 && !allowAdjacentUnsaturatedBonds) {
								continue bondLoop;
							}
							if (bondsToUnsaturate.contains(bond2)) {
								if (!allowAdjacentUnsaturatedBonds) {
									continue bondLoop;
								}
								incomingValency2 += bondOrder;
							}
							else {
								incomingValency2 += bond2.getOrder();
							}
						}
						
						Integer maxVal2 = getLambdaValencyOrHwValencyOrMaxValIfCharged(atom2);
						if(maxVal2 != null && (incomingValency2 + (bondOrder - 1) + atom2.getOutValency()) > maxVal2) {
							continue;
						}
						bondsToUnsaturate.add(bond);
						break bondLoop;
					}
				}
			}
		}
		return bondsToUnsaturate;
	}
	
	
	/**
	 * Return the lambda convention derived valency + protons if set
	 * Otherwise if charge is 0 returns {@link ValencyChecker#getHWValency(ChemEl)}
	 * Otherwise return {@link ValencyChecker#getMaximumValency(ChemEl, int)} 
	 * Returns null if the maximum valency is not known
	 * @param a
	 * @return
	 */
	static Integer getLambdaValencyOrHwValencyOrMaxValIfCharged(Atom a) {
		if (a.getLambdaConventionValency() != null) {
			return a.getLambdaConventionValency() + a.getProtonsExplicitlyAddedOrRemoved();
		}
		else if (a.getCharge() == 0){
			return ValencyChecker.getHWValency(a.getElement());
		}
		else {
			return ValencyChecker.getMaximumValency(a.getElement(), a.getCharge());
		}
	}

	private static void performAdditiveOperations(BuildState state, Element subBracketOrRoot) throws StructureBuildingException {
		if (subBracketOrRoot.getAttribute(LOCANT_ATR) != null){//additive nomenclature does not employ locants
			return;
		}
		Element group;
		if (subBracketOrRoot.getName().equals(BRACKET_EL)){
			group =findRightMostGroupInBracket(subBracketOrRoot);
		}
		else{
			group =subBracketOrRoot.getFirstChildElement(GROUP_EL);
		}
		if (group.getAttribute(RESOLVED_ATR) != null){
			return;
		}
		Fragment frag = group.getFrag();
		int outAtomCount = frag.getOutAtomCount();
		if (outAtomCount >=1){
			if (subBracketOrRoot.getAttribute(MULTIPLIER_ATR) ==null){
				Element nextSiblingEl = OpsinTools.getNextSibling(subBracketOrRoot);
				if (nextSiblingEl.getAttribute(MULTIPLIER_ATR) != null &&
						(outAtomCount >= Integer.parseInt(nextSiblingEl.getAttributeValue(MULTIPLIER_ATR)) || //probably multiplicative nomenclature, should be as many outAtoms as the multiplier
						outAtomCount==1 && frag.getOutAtom(0).getValency()==Integer.parseInt(nextSiblingEl.getAttributeValue(MULTIPLIER_ATR))) &&
						hasRootLikeOrMultiRadicalGroup(nextSiblingEl)){
					if (outAtomCount==1){//special case e.g. 4,4'-(benzylidene)dianiline
						FragmentTools.splitOutAtomIntoValency1OutAtoms(frag.getOutAtom(0));
						//special case where something like benzylidene is being used as if it meant benzdiyl for multiplicative nomenclature
						//this is allowed in the IUPAC 79 recommendations but not recommended in the current recommendations
					}
					performMultiplicativeOperations(state, group, nextSiblingEl);
				}
				else if (group.getAttribute(ISAMULTIRADICAL_ATR) != null){//additive nomenclature e.g. ethyleneoxy
					Fragment nextFrag = getNextInScopeMultiValentFragment(subBracketOrRoot);
					if (nextFrag != null){
						Element nextMultiRadicalGroup = nextFrag.getTokenEl();
						Element parentSubOrRoot = nextMultiRadicalGroup.getParent();
						if (state.currentWordRule != WordRule.polymer){//imino does not behave like a substituent in polymers only as a linker
							if (nextMultiRadicalGroup.getAttribute(IMINOLIKE_ATR) != null){//imino/methylene can just act as normal substituents, should an additive bond really be made???
								Fragment adjacentFrag = OpsinTools.getNextGroup(subBracketOrRoot).getFrag();
								
								if (nextFrag != adjacentFrag){//imino is not the absolute next frag
									if (potentiallyCanSubstitute(nextMultiRadicalGroup.getParent()) || potentiallyCanSubstitute(nextMultiRadicalGroup.getParent().getParent())){
										return;
									}
								}
							}
							if (group.getAttribute(IMINOLIKE_ATR) != null && levelsToWordEl(group) > levelsToWordEl(nextMultiRadicalGroup)){
								return;//e.g. imino substitutes ((chloroimino)ethylene)dibenzene
							}
						}
						if (parentSubOrRoot.getAttribute(MULTIPLIER_ATR) != null){
							throw new StructureBuildingException("Attempted to form additive bond to a multiplied component");
						}
						group.addAttribute(new Attribute(RESOLVED_ATR, "yes"));
						joinFragmentsAdditively(state, frag, nextFrag);
					}
				}
				else {//e.g. chlorocarbonyl or hydroxy(sulfanyl)phosphoryl
					List<Fragment> siblingFragments = findAlternativeFragments(subBracketOrRoot);
					if (siblingFragments.size()>0){
						Fragment nextFrag = siblingFragments.get(siblingFragments.size()-1);
						Element nextGroup = nextFrag.getTokenEl();
						if (nextGroup.getAttribute(ACCEPTSADDITIVEBONDS_ATR) != null && nextGroup.getAttribute(ISAMULTIRADICAL_ATR) != null && (nextFrag.getOutAtomCount()>1|| nextGroup.getAttribute(RESOLVED_ATR) != null && nextFrag.getOutAtomCount()>=1 )){
							Atom toAtom = nextFrag.getOutAtom(0).getAtom();
							if (calculateSubstitutableHydrogenAtoms(toAtom) ==0){
								group.addAttribute(new Attribute(RESOLVED_ATR, "yes"));
								joinFragmentsAdditively(state, frag, nextFrag);//e.g. aminocarbonyl or aminothio
							}
						}
						if (group.getAttribute(RESOLVED_ATR)==null && siblingFragments.size()>1){
							for (int i = 0; i< siblingFragments.size()-1; i++) {
								Fragment lastFrag = siblingFragments.get(i);
								Element lastGroup = lastFrag.getTokenEl();
								if (lastGroup.getAttribute(ACCEPTSADDITIVEBONDS_ATR) != null && lastGroup.getAttribute(ISAMULTIRADICAL_ATR) != null && (lastFrag.getOutAtomCount()>1|| lastGroup.getAttribute(RESOLVED_ATR) != null && lastFrag.getOutAtomCount()>=1 )){
									Atom toAtom = lastFrag.getOutAtom(0).getAtom();
									if (calculateSubstitutableHydrogenAtoms(toAtom) ==0){
										group.addAttribute(new Attribute(RESOLVED_ATR, "yes"));
										joinFragmentsAdditively(state, frag, lastFrag);//e.g. hydroxy(sulfanyl)phosphoryl
									}
									break;
								}

								//loop may continue if lastFrag was in fact completely unsubstitutable e.g. hydroxy...phosphoryloxy. The oxy is unsubstituable as the phosphoryl will already have bonded to it
								if (FragmentTools.findSubstituableAtoms(lastFrag, frag.getOutAtom(outAtomCount - 1).getValency()).size() > 0) {
									break;
								}
							}
						}
					}
				}
			}
			else{// e.g. dimethoxyphosphoryl or bis(methylamino)phosphoryl
				List<Fragment> siblingFragments = findAlternativeFragments(subBracketOrRoot);
				if (siblingFragments.size()>0){
					int multiplier = Integer.parseInt(subBracketOrRoot.getAttributeValue(MULTIPLIER_ATR));
					Fragment nextFrag = siblingFragments.get(siblingFragments.size()-1);
					Element nextGroup = nextFrag.getTokenEl();
					if (nextGroup.getAttribute(ACCEPTSADDITIVEBONDS_ATR) != null && nextGroup.getAttribute(ISAMULTIRADICAL_ATR) != null && (nextFrag.getOutAtomCount()>=multiplier|| nextGroup.getAttribute(RESOLVED_ATR) != null && nextFrag.getOutAtomCount()>=multiplier +1 )){
						Atom toAtom = nextFrag.getOutAtom(0).getAtom();
						if (calculateSubstitutableHydrogenAtoms(toAtom) ==0){
							group.addAttribute(new Attribute(RESOLVED_ATR, "yes"));
							multiplyOutAndAdditivelyBond(state, subBracketOrRoot, nextFrag);//e.g.dihydroxyphosphoryl
						}
					}
					if (group.getAttribute(RESOLVED_ATR)==null && siblingFragments.size()>1){
						for (int i = 0; i< siblingFragments.size()-1; i++) {
							Fragment lastFrag = siblingFragments.get(i);
							Element lastGroup = lastFrag.getTokenEl();
							if (lastGroup.getAttribute(ACCEPTSADDITIVEBONDS_ATR) != null && lastGroup.getAttribute(ISAMULTIRADICAL_ATR) != null && (lastFrag.getOutAtomCount()>=multiplier|| lastGroup.getAttribute(RESOLVED_ATR) != null && lastFrag.getOutAtomCount()>=multiplier +1 )){
								Atom toAtom = lastFrag.getOutAtom(0).getAtom();
								if (calculateSubstitutableHydrogenAtoms(toAtom) ==0){
									group.addAttribute(new Attribute(RESOLVED_ATR, "yes"));
									multiplyOutAndAdditivelyBond(state, subBracketOrRoot, lastFrag);//e.g. dihydroxyphosphoryloxy
								}
								break;
							}

							//loop may continue if lastFrag was in fact completely unsubstitutable e.g. hydroxy...phosphoryloxy. The oxy is unsubstituable as the phosphoryl will already have bonded to it
							if (FragmentTools.findSubstituableAtoms(lastFrag, frag.getOutAtom(outAtomCount - 1).getValency()).size() > 0) {
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
	 * @param subBracketOrRoot
	 * @return
	 */
	private static boolean hasRootLikeOrMultiRadicalGroup(Element subBracketOrRoot) {
		List<Element> groups = OpsinTools.getDescendantElementsWithTagName(subBracketOrRoot, GROUP_EL);
		if (subBracketOrRoot.getAttribute(INLOCANTS_ATR) != null){
			return true;// a terminus with specified inLocants
		}
		for (Element group : groups) {
			Fragment frag = group.getFrag();
			int outAtomCount =frag.getOutAtomCount();
			if (group.getAttribute(ISAMULTIRADICAL_ATR) != null){
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
			if (i != 0){
				currentElement = state.fragManager.cloneElement(state, subOrBracket, i);
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
					Element el = subOrBracket.getChild(j);
					el.detach();
					elementsNotToBeMultiplied.add(el);
				}
				multiplierEl.detach();
			}
			Element group;
			if (currentElement.getName().equals(BRACKET_EL)){
				group = findRightMostGroupInBracket(currentElement);
			}
			else{
				group = currentElement.getFirstChildElement(GROUP_EL);
			}
			Fragment frag = group.getFrag();
			if (frag.getOutAtomCount() != 1 ){
				throw new StructureBuildingException("Additive bond formation failure: Fragment expected to have one OutAtom in this case but had: "+ frag.getOutAtomCount());
			}
			joinFragmentsAdditively(state, frag, fragToAdditivelyBondTo);
		}
		for (Element clone : clonedElements) {//make sure cloned substituents don't substitute onto each other!
			OpsinTools.insertAfter(subOrBracket, clone);
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
		BuildResults multiRadicalBR = new BuildResults(group.getParent());
		performMultiplicativeOperations(state, multiRadicalBR, multipliedParent);
	}

	private static void performMultiplicativeOperations(BuildState state, BuildResults multiRadicalBR, Element multipliedParent) throws StructureBuildingException {
		int multiplier = Integer.parseInt(multipliedParent.getAttributeValue(MULTIPLIER_ATR));
		if (multiplier != multiRadicalBR.getOutAtomCount()){
			if (multiRadicalBR.getOutAtomCount() == multiplier*2){
				//TODO substituents like nitrilo can have their outatoms combined
			}
			if (multiplier != multiRadicalBR.getOutAtomCount()){
				throw new StructureBuildingException("Multiplication bond formation failure: number of outAtoms disagree with multiplier(multiplier: " + multiplier + ", outAtom count: " + multiRadicalBR.getOutAtomCount()+ ")");
			}
		}
		if (LOG.isTraceEnabled()){LOG.trace(multiplier +" multiplicative bonds to be formed");}
		multipliedParent.removeAttribute(multipliedParent.getAttribute(MULTIPLIER_ATR));
		List<String> inLocants = null;
		String inLocantsString = multipliedParent.getAttributeValue(INLOCANTS_ATR);
		if (inLocantsString != null){//true for the root of a multiplicative name
			if (inLocantsString.equals(INLOCANTS_DEFAULT)){
				inLocants = new ArrayList<String>(multiplier);
				for (int i = 0; i < multiplier; i++) {
					inLocants.add(INLOCANTS_DEFAULT);
				}
			}
			else{
				inLocants = StringTools.arrayToList(inLocantsString.split(","));
				if (inLocants.size() != multiplier){
					throw new StructureBuildingException("Mismatch between multiplier and number of inLocants in multiplicative nomenclature");
				}
			}
		}
		List<Element> clonedElements = new ArrayList<Element>();
		BuildResults newBr = new BuildResults();
		for (int i = multiplier -1; i >=0; i--) {
			Element multipliedElement;
			if (i != 0){
				multipliedElement = state.fragManager.cloneElement(state, multipliedParent, i);
				addPrimesToLocantedStereochemistryElements(multipliedElement, StringTools.multiplyString("'", i));//Stereochemistry elements with locants will need to have their locants primed (stereochemistry is only processed after structure building)
				clonedElements.add(multipliedElement);
			}
			else{
				multipliedElement = multipliedParent;
			}
			
			//determine group that will be additively bonded to
			Element multipliedGroup;
			if (multipliedElement.getName().equals(BRACKET_EL)) {
				multipliedGroup = getFirstMultiValentGroup(multipliedElement);
				if (multipliedGroup == null){//root will not have a multivalent group
					List<Element> groups = OpsinTools.getDescendantElementsWithTagName(multipliedElement, GROUP_EL);
					if (inLocants == null){
						throw new StructureBuildingException("OPSIN Bug? in locants must be specified for a multiplied root in multiplicative nomenclature");
					}
					if (inLocants.get(0).equals(INLOCANTS_DEFAULT)){
						multipliedGroup = groups.get(groups.size() - 1);
					}
					else{
						groupLoop: for (int j = groups.size()-1; j >=0; j--) {
							Fragment possibleFrag = groups.get(j).getFrag();
							for (String locant : inLocants) {
								if (possibleFrag.hasLocant(locant)){
									multipliedGroup = groups.get(j);
									break groupLoop;
								}
							}
						}
					}
					if (multipliedGroup == null){
						throw new StructureBuildingException("Locants for inAtoms on the root were either misassigned to the root or were invalid: " + inLocants.toString() +" could not be assigned!");
					}
				}
			}
			else{
				multipliedGroup = multipliedElement.getFirstChildElement(GROUP_EL);
			}
			Fragment multipliedFrag = multipliedGroup.getFrag();
			
			OutAtom multiRadicalOutAtom = multiRadicalBR.getOutAtom(i);
			Fragment multiRadicalFrag = multiRadicalOutAtom.getAtom().getFrag();
			Element multiRadicalGroup = multiRadicalFrag.getTokenEl();
			if (multiRadicalGroup.getAttribute(RESOLVED_ATR) == null){
				resolveUnLocantedFeatures(state, multiRadicalGroup.getParent());//the addition of unlocanted unsaturators can effect the position of radicals e.g. diazenyl
				multiRadicalGroup.addAttribute(new Attribute(RESOLVED_ATR, "yes"));
			}

			boolean substitutivelyBondedToRoot = false;
			if (inLocants != null) {
				Element rightMostGroup;
				if (multipliedElement.getName().equals(BRACKET_EL)) {
					rightMostGroup = findRightMostGroupInBracket(multipliedElement);
				}
				else{
					rightMostGroup = multipliedElement.getFirstChildElement(GROUP_EL);
				}
				rightMostGroup.addAttribute(new Attribute(RESOLVED_ATR, "yes"));//this group will not be used further within this word but can in principle be a substituent e.g. methylenedisulfonyl dichloride
				if (multipliedGroup.getAttribute(ISAMULTIRADICAL_ATR) != null) {//e.g. methylenedisulfonyl dichloride
					if (!multipliedParent.getAttributeValue(INLOCANTS_ATR).equals(INLOCANTS_DEFAULT)) {
						throw new StructureBuildingException("inLocants should not be specified for a multiradical parent in multiplicative nomenclature");
					}
				}
				else{
					Atom from = multiRadicalOutAtom.getAtom();
					int bondOrder = multiRadicalOutAtom.getValency();
					//bonding will be substitutive rather additive as this is bonding to a root
					Atom atomToJoinTo = null;
					for (int j = inLocants.size() -1; j >=0; j--) {
						String locant = inLocants.get(j);
						if (locant.equals(INLOCANTS_DEFAULT)){//note that if one entry in inLocantArray is default then they all are "default"
							List<Atom> possibleAtoms = getPossibleAtomsForUnlocantedConnectionToMultipliedRoot(multipliedGroup, bondOrder, i);
							if (possibleAtoms.isEmpty()) {
								throw new StructureBuildingException("No suitable atom found for multiplicative operation");
							}
							if (AmbiguityChecker.isSubstitutionAmbiguous(possibleAtoms, 1)) {
								state.addIsAmbiguous("Connection to multiplied group: " + multipliedGroup.getValue());
							}
							atomToJoinTo = possibleAtoms.get(0);
							inLocants.remove(j);
							break;
						}
						else{
							Atom inAtom = multipliedFrag.getAtomByLocant(locant);
							if (inAtom != null) {
								atomToJoinTo = inAtom;
								inLocants.remove(j);
								break;
							}
						}
					}
					if (atomToJoinTo == null){
						throw new StructureBuildingException("Locants for inAtoms on the root were either misassigned to the root or were invalid: " + inLocants.toString() +" could not be assigned!");
					}

					if (!multiRadicalOutAtom.isSetExplicitly()) {//not set explicitly so may be an inappropriate atom
						from = findAtomForUnlocantedRadical(state, from.getFrag(), multiRadicalOutAtom);
					}
					multiRadicalFrag.removeOutAtom(multiRadicalOutAtom);

					state.fragManager.createBond(from, atomToJoinTo, bondOrder);
					if (LOG.isTraceEnabled()){LOG.trace("Substitutively bonded (multiplicative to root) " + from.getID() + " (" + from.getFrag().getTokenEl().getValue() + ") " + atomToJoinTo.getID() + " (" + atomToJoinTo.getFrag().getTokenEl().getValue() + ")");}
					substitutivelyBondedToRoot = true;
				}
			}
			if (!substitutivelyBondedToRoot) {
				joinFragmentsAdditively(state, multiRadicalFrag, multipliedFrag);
			}
			if (multipliedElement.getName().equals(BRACKET_EL)) {
				recursivelyResolveUnLocantedFeatures(state, multipliedElement);//there may be outAtoms that are involved in unlocanted substitution, these can be safely used now e.g. ...bis((3-hydroxy-4-methoxyphenyl)methylene) where (3-hydroxy-4-methoxyphenyl)methylene is the currentElement
			}

			if (inLocants == null) {
				//currentElement is not a root element. Need to build up a new BuildResults so as to call performMultiplicativeOperations again
				//at this stage an outAtom has been removed from the fragment within currentElement through an additive bond
				newBr.mergeBuildResults(new BuildResults(multipliedElement));
			}
		}

		if (newBr.getFragmentCount() == 1) {
			throw new StructureBuildingException("Multiplicative nomenclature cannot yield only one temporary terminal fragment");
		}
		if (newBr.getFragmentCount() >= 2) {
			List<Element> siblings = OpsinTools.getNextSiblingsOfTypes(multipliedParent, new String[]{SUBSTITUENT_EL, BRACKET_EL, ROOT_EL});
			if (siblings.size() == 0) {
				Element parentOfMultipliedEl = multipliedParent.getParent();
				if (parentOfMultipliedEl.getName().equals(BRACKET_EL)) {//brackets are allowed
					siblings = OpsinTools.getNextSiblingsOfTypes(parentOfMultipliedEl, new String[]{SUBSTITUENT_EL, BRACKET_EL, ROOT_EL});
					if (siblings.get(0).getAttribute(MULTIPLIER_ATR) == null) {
						throw new StructureBuildingException("Multiplier not found where multiplier was expected for succesful multiplicative nomenclature");
					}
					performMultiplicativeOperations(state, newBr, siblings.get(0));
				}
				else{
					throw new StructureBuildingException("Could not find suitable element to continue multiplicative nomenclature");
				}
			}
			else{
				if (siblings.get(0).getAttribute(MULTIPLIER_ATR) == null) {
					throw new StructureBuildingException("Multiplier not found where multiplier was expected for successful multiplicative nomenclature");
				}
				performMultiplicativeOperations(state, newBr, siblings.get(0));
			}
		}

		for (Element clone : clonedElements) {//only insert cloned substituents now so they don't substitute onto each other!
			OpsinTools.insertAfter(multipliedParent, clone);
		}
	}

	/**
	 * Applies special case to prefer the end of chains with the usableAsAJoiner attributes cf. p-phenylenedipropionic acid
	 * Such cases will still be considered to be formally ambiguous
	 * @param multipliedGroup
	 * @param multipliedFrag
	 * @param bondOrder
	 * @param primesAdded 
	 * @return
	 * @throws StructureBuildingException
	 */
	private static List<Atom> getPossibleAtomsForUnlocantedConnectionToMultipliedRoot(Element multipliedGroup, int bondOrder, int primesAdded) throws StructureBuildingException {
		Fragment multipliedFrag = multipliedGroup.getFrag();
		if ("yes".equals(multipliedGroup.getAttributeValue(USABLEASJOINER_ATR)) && multipliedFrag.getDefaultInAtom() == null) {
			Element previous = OpsinTools.getPrevious(multipliedGroup);
			if (previous != null && previous.getName().equals(MULTIPLIER_EL)){
				String locant = getLocantOfEndOfChainIfGreaterThan1(multipliedFrag, primesAdded);
				if (locant != null) {
					Atom preferredAtom = multipliedFrag.getAtomByLocantOrThrow(locant);
					List<Atom> possibleAtoms = FragmentTools.findnAtomsForSubstitution(multipliedFrag.getAtomList(), preferredAtom, 1, bondOrder, true);
					if (possibleAtoms == null) {
						possibleAtoms = Collections.emptyList();
					}
					return possibleAtoms;
				}
			}
		}
		return FragmentTools.findSubstituableAtoms(multipliedFrag, bondOrder);
	}
	
	private static String getLocantOfEndOfChainIfGreaterThan1(Fragment frag, int primes) {
		String primesStr = StringTools.multiplyString("'", primes);
		int length = 0;
		Atom next = frag.getAtomByLocant(Integer.toString(length + 1) + primesStr);
		Atom previous = null;
		while (next != null){
			if (previous != null && previous.getBondToAtom(next) == null){
				break;
			}
			length++;
			previous = next;
			next = frag.getAtomByLocant(Integer.toString(length + 1) + primesStr);
		}
		if (length > 1){
			return Integer.toString(length) + primesStr;
		}
		return null;
	}

	/**
	 * Given a subsituent/bracket finds the next multi valent substituent/root that is in scope and hence its group
	 * e.g. for oxy(dichloromethyl)methylene given oxy substituent the methylene group would be found
	 * for oxy(dichloroethylene) given oxy substituent the ethylene group would be found
	 * for oxy(carbonylimino) given oxy carbonyl would be found
	 * @param substituentOrBracket
	 * @return frag
	 * @throws StructureBuildingException
	 */
	private static Fragment getNextInScopeMultiValentFragment(Element substituentOrBracket) throws StructureBuildingException {
		if (!substituentOrBracket.getName().equals(SUBSTITUENT_EL) && !substituentOrBracket.getName().equals(BRACKET_EL)){
			throw new StructureBuildingException("Input to this function should be a substituent or bracket");
		}
		if (substituentOrBracket.getParent()==null){
			throw new StructureBuildingException("substituent did not have a parent!");
		}
		Element parent = substituentOrBracket.getParent();

		List<Element> children = OpsinTools.getChildElementsWithTagNames(parent, new String[]{SUBSTITUENT_EL, BRACKET_EL, ROOT_EL});//will be returned in index order
		int indexOfSubstituent =parent.indexOf(substituentOrBracket);
		for (Element child : children) {
			if (parent.indexOf(child) <=indexOfSubstituent){//only want things after the input
				continue;
			}
			if (child.getAttribute(MULTIPLIER_ATR) != null){
				continue;
			}
			List<Element> childDescendants;
			if (child.getName().equals(BRACKET_EL)){
				childDescendants = OpsinTools.getDescendantElementsWithTagNames(child, new String[]{SUBSTITUENT_EL, ROOT_EL});//will be returned in depth-first order
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
				Fragment possibleFrag = group.getFrag();
				if (group.getAttribute(ISAMULTIRADICAL_ATR) != null &&
						(possibleFrag.getOutAtomCount() >=2 || (possibleFrag.getOutAtomCount() >=1 && group.getAttribute(RESOLVED_ATR) != null ))){
					return possibleFrag;
				}
			}
		}
		return null;
	}

	/**
	 * Given a bracket searches in a depth first manner for the first multi valent group
	 * @param bracket
	 * @return group
	 * @throws StructureBuildingException
	 */
	private static Element getFirstMultiValentGroup(Element bracket) throws StructureBuildingException {
		if (!bracket.getName().equals(BRACKET_EL)){
			throw new StructureBuildingException("Input to this function should be a bracket");
		}

		List<Element> groups = OpsinTools.getDescendantElementsWithTagName(bracket, GROUP_EL);//will be returned in index order
		for (Element group : groups) {
			Fragment possibleFrag = group.getFrag();
			if (group.getAttribute(ISAMULTIRADICAL_ATR) != null &&
					(possibleFrag.getOutAtomCount() >=2 || (possibleFrag.getOutAtomCount() >=1 && group.getAttribute(RESOLVED_ATR) != null ))){
				return group;
			}
		}
		return null;
	}

	private static void joinFragmentsAdditively(BuildState state, Fragment fragToBeJoined, Fragment parentFrag) throws StructureBuildingException {
		Element elOfFragToBeJoined = fragToBeJoined.getTokenEl();
		if (EPOXYLIKE_SUBTYPE_VAL.equals(elOfFragToBeJoined.getAttributeValue(SUBTYPE_ATR))){
			for (int i = 0, l = fragToBeJoined.getOutAtomCount(); i < l; i++) {
				OutAtom outAtom = fragToBeJoined.getOutAtom(i);
				if (outAtom.getLocant() != null){
					throw new StructureBuildingException("Inappropriate use of " + elOfFragToBeJoined.getValue());
				}
			}
		}
		int outAtomCountOnFragToBeJoined = fragToBeJoined.getOutAtomCount();
		if (outAtomCountOnFragToBeJoined ==0){
			throw new StructureBuildingException("Additive bond formation failure: Fragment expected to have at least one OutAtom but had none");
		}

		if (parentFrag.getOutAtomCount() == 0){
			throw new StructureBuildingException("Additive bond formation failure: Fragment expected to have at least one OutAtom but had none");
		}
		OutAtom in = null;
		if (parentFrag.getOutAtomCount() > 1){
			int firstOutAtomOrder = parentFrag.getOutAtom(0).getValency();
			boolean unresolvedAmbiguity =false;
			for (int i = 1, l = parentFrag.getOutAtomCount(); i < l; i++) {
				OutAtom outAtom = parentFrag.getOutAtom(i);
				if (outAtom.getValency() != firstOutAtomOrder){
					unresolvedAmbiguity =true;
				}
			}
			if (unresolvedAmbiguity){//not all outAtoms on parent equivalent
				firstOutAtomOrder = fragToBeJoined.getOutAtom(0).getValency();
				unresolvedAmbiguity =false;
				for (int i = 1, l = fragToBeJoined.getOutAtomCount(); i < l; i++) {
					OutAtom outAtom = fragToBeJoined.getOutAtom(i);
					if (outAtom.getValency() != firstOutAtomOrder){
						unresolvedAmbiguity =true;
					}
				}
				if (unresolvedAmbiguity && outAtomCountOnFragToBeJoined == 2){//not all outAtoms on frag to be joined are equivalent either!
					//Solves the specific case of 2,2'-[ethane-1,2-diylbis(azanylylidenemethanylylidene)]diphenol vs 2,2'-[ethane-1,2-diylidenebis(azanylylidenemethanylylidene)]bis(cyclohexan-1-ol)
					//but does not solve the general case as only a single look behind is performed.
					Element previousGroup = OpsinTools.getPreviousGroup(elOfFragToBeJoined);
					if (previousGroup != null){
						Fragment previousFrag = previousGroup.getFrag();
						if (previousFrag.getOutAtomCount() > 1){
							int previousGroupFirstOutAtomOrder = previousFrag.getOutAtom(0).getValency();
							unresolvedAmbiguity =false;
							for (int i = 1, l = previousFrag.getOutAtomCount(); i < l; i++) {
								OutAtom outAtom = previousFrag.getOutAtom(i);
								if (outAtom.getValency() != previousGroupFirstOutAtomOrder){
									unresolvedAmbiguity =true;
								}
							}
							if (!unresolvedAmbiguity && previousGroupFirstOutAtomOrder==parentFrag.getOutAtom(0).getValency()){
								for (int i = 1, l = parentFrag.getOutAtomCount(); i < l; i++) {
									OutAtom outAtom = parentFrag.getOutAtom(i);
									if (outAtom.getValency() != previousGroupFirstOutAtomOrder){
										in = outAtom;
										break;
									}
								}
							}
						}
					}
				}
				else{
					for (int i = 0, l = parentFrag.getOutAtomCount(); i < l; i++) {
						OutAtom outAtom = parentFrag.getOutAtom(i);
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
		Atom to = in.getAtom();
		int bondOrder = in.getValency();
		if (!in.isSetExplicitly()){//not set explicitly so may be an inappropriate atom
			to = findAtomForUnlocantedRadical(state, to.getFrag(), in);
		}
		parentFrag.removeOutAtom(in);

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
					if (nextOutAtom.getAtom() != lastOutAtom){
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
			from = findAtomForUnlocantedRadical(state, from.getFrag(), out);
		}
		fragToBeJoined.removeOutAtom(out);

		state.fragManager.createBond(from, to, bondOrder);
		if (LOG.isTraceEnabled()){LOG.trace("Additively bonded " + from.getID() + " (" + from.getFrag().getTokenEl().getValue() + ") " + to.getID() + " (" + to.getFrag().getTokenEl().getValue() + ")" );}
	}

	private static void joinFragmentsSubstitutively(BuildState state, Fragment fragToBeJoined, Atom atomToJoinTo) throws StructureBuildingException {
		Element elOfFragToBeJoined = fragToBeJoined.getTokenEl();
		if (EPOXYLIKE_SUBTYPE_VAL.equals(elOfFragToBeJoined.getAttributeValue(SUBTYPE_ATR))){
			formEpoxide(state, fragToBeJoined, atomToJoinTo);
			return;
		}
		int outAtomCount = fragToBeJoined.getOutAtomCount();
		if (outAtomCount >1){
			throw new StructureBuildingException("Substitutive bond formation failure: Fragment expected to have one OutAtom but had: "+ outAtomCount);
		}
		if (outAtomCount ==0 ){
			throw new StructureBuildingException("Substitutive bond formation failure: Fragment expected to have one OutAtom but had none");
		}
		if (elOfFragToBeJoined.getAttribute(IMINOLIKE_ATR) != null){//special case for methylene/imino
			if (fragToBeJoined.getOutAtomCount()==1 && fragToBeJoined.getOutAtom(0).getValency()==1 ){
				fragToBeJoined.getOutAtom(0).setValency(2);
			}
		}
		OutAtom out = fragToBeJoined.getOutAtom(0);
		Atom from = out.getAtom();
		int bondOrder = out.getValency();
		if (!out.isSetExplicitly()){//not set explicitly so may be an inappropriate atom
			List<Atom> possibleAtoms = FragmentTools.findnAtomsForSubstitution(fragToBeJoined.getAtomList(), from, 1, bondOrder, false);
			if (possibleAtoms == null){
				throw new StructureBuildingException("Failed to assign all unlocanted radicals to actual atoms without violating valency");
			}
			if (!((ALKANESTEM_SUBTYPE_VAL.equals(fragToBeJoined.getSubType()) || HETEROSTEM_SUBTYPE_VAL.equals(fragToBeJoined.getSubType())) && possibleAtoms.get(0).equals(fragToBeJoined.getFirstAtom()))) {
				if (AmbiguityChecker.isSubstitutionAmbiguous(possibleAtoms, 1)) {
					state.addIsAmbiguous("Positioning of radical on: " + fragToBeJoined.getTokenEl().getValue());
				}
			}
			from = possibleAtoms.get(0);
		}
		fragToBeJoined.removeOutAtom(out);

		state.fragManager.createBond(from, atomToJoinTo, bondOrder);
		if (LOG.isTraceEnabled()){LOG.trace("Substitutively bonded " + from.getID() + " (" + from.getFrag().getTokenEl().getValue() + ") " + atomToJoinTo.getID() + " (" + atomToJoinTo.getFrag().getTokenEl().getValue() + ")");}
	}

	/**
	 * Forms a bridge using the given fragment.
	 * The bridgingFragment's outAtoms locants or a combination of the atomToJoinTo and a suitable atom
	 * are used to decide what atoms to form the bridge between
	 * @param state
	 * @param bridgingFragment
	 * @param atomToJoinTo
	 * @return Atoms that the bridgingFragment attached to
	 * @throws StructureBuildingException
	 */
	static Atom[] formEpoxide(BuildState state, Fragment bridgingFragment, Atom atomToJoinTo) throws StructureBuildingException {
		Fragment fragToJoinTo = atomToJoinTo.getFrag();
		List<Atom> atomList = fragToJoinTo.getAtomList();
		if (atomList.size()==1){
			throw new StructureBuildingException("Epoxides must be formed between two different atoms");
		}
		Atom firstAtomToJoinTo;
		if (bridgingFragment.getOutAtom(0).getLocant() != null){
			firstAtomToJoinTo = fragToJoinTo.getAtomByLocantOrThrow(bridgingFragment.getOutAtom(0).getLocant());
		}
		else{
			firstAtomToJoinTo = atomToJoinTo;
		}
		OutAtom outAtom1 = bridgingFragment.getOutAtom(0);
		bridgingFragment.removeOutAtom(0);
		
		//In epoxy chalcogenAtom1 will be chalcogenAtom2. Methylenedioxy is also handled by this method
		state.fragManager.createBond(outAtom1.getAtom(), firstAtomToJoinTo, outAtom1.getValency());
		
		Atom secondAtomToJoinTo;
		if (bridgingFragment.getOutAtom(0).getLocant() != null){
			secondAtomToJoinTo = fragToJoinTo.getAtomByLocantOrThrow(bridgingFragment.getOutAtom(0).getLocant());
		}
		else{
			int index = atomList.indexOf(firstAtomToJoinTo);
			Atom preferredAtom = (index + 1 >= atomList.size()) ? atomList.get(index - 1) : atomList.get(index + 1);
			List<Atom> possibleSecondAtom = FragmentTools.findnAtomsForSubstitution(fragToJoinTo.getAtomList(), preferredAtom, 1, 1, true);
			if (possibleSecondAtom != null) {
				possibleSecondAtom.removeAll(Collections.singleton(firstAtomToJoinTo));
			}
			if (possibleSecondAtom == null || possibleSecondAtom.size() == 0) {
				throw new StructureBuildingException("Unable to find suitable atom to form bridge");
			}
			if (AmbiguityChecker.isSubstitutionAmbiguous(possibleSecondAtom, 1)) {
				state.addIsAmbiguous("Addition of bridge to: "+ fragToJoinTo.getTokenEl().getValue());
			}
			secondAtomToJoinTo = possibleSecondAtom.get(0);
		}
		OutAtom outAtom2 = bridgingFragment.getOutAtom(0);
		bridgingFragment.removeOutAtom(0);
		if (outAtom1.getAtom().equals(outAtom2.getAtom()) && firstAtomToJoinTo == secondAtomToJoinTo){
			throw new StructureBuildingException("Epoxides must be formed between two different atoms");
		}
		int bondValency = outAtom2.getValency();
		if (outAtom2.getAtom().hasSpareValency() && !secondAtomToJoinTo.hasSpareValency()) {
			//bridging groups like azeno are treated as aromatic so that it is not fixed as to which of the two bonds is the double bond
			//if connected to a saturated group though, one of them must be a double bond
			bondValency = 2;
		}
		state.fragManager.createBond(outAtom2.getAtom(), secondAtomToJoinTo, bondValency);
		CycleDetector.assignWhetherAtomsAreInCycles(bridgingFragment);
		return new Atom[]{firstAtomToJoinTo, secondAtomToJoinTo};
	}
	
	/**
	 * Attempts to find an in-scope fragment capable of forming the given numberOfSubstitutions each with the given bondOrder
	 * @param subOrBracket
	 * @param numberOfSubstitutions
	 * @param bondOrder
	 * @return
	 */
	private static List<Atom> findAtomsForSubstitution(Element subOrBracket, int numberOfSubstitutions, int bondOrder) {
		boolean rootHandled = false;
		List<Element> possibleParents = findAlternativeGroups(subOrBracket);
		for (int i = 0, l = possibleParents.size(); i < l; i++) {
			Element possibleParent = possibleParents.get(i);
			Fragment frag = possibleParent.getFrag();
			List<Atom> substitutableAtoms;
			if (possibleParent.getParent().getName().equals(ROOT_EL)){//consider all root groups as if they were one
				if(rootHandled) {
					continue;
				}
				List<Atom> atoms = frag.getAtomList();
				for (int j = i + 1; j < l; j++) {
					Element possibleOtherRoot = possibleParents.get(j);
					if (possibleOtherRoot.getParent().getName().equals(ROOT_EL)) {
						atoms.addAll(possibleOtherRoot.getFrag().getAtomList());
					}
				}
				rootHandled = true;
				substitutableAtoms = FragmentTools.findnAtomsForSubstitution(atoms, frag.getDefaultInAtom(), numberOfSubstitutions, bondOrder, true);
			}
			else{
				substitutableAtoms = FragmentTools.findnAtomsForSubstitution(frag, numberOfSubstitutions, bondOrder);
			}
			if (substitutableAtoms != null){
				return substitutableAtoms;
			}
		}
		return null;
	}

	/**
	 * Finds all the fragments accessible from the startingElement taking into account brackets
	 * i.e. those that it is feasible that the group of the startingElement could substitute onto
	 * @param startingElement
	 * @return A list of fragments in the order to try them as possible parent fragments (for substitutive operations)
	 */
	static List<Fragment> findAlternativeFragments(Element startingElement) {
		List<Fragment> foundFragments = new ArrayList<Fragment>();
		for (Element group : findAlternativeGroups(startingElement)) {
			foundFragments.add(group.getFrag());
		}
		return foundFragments;
	}
	
	/**
	 * Finds all the groups accessible from the startingElement taking into account brackets
	 * i.e. those that it is feasible that the group of the startingElement could substitute onto
	 * @param startingElement
	 * @return A list of groups in the order to try them as possible parent groups (for substitutive operations)
	 */
	static List<Element> findAlternativeGroups(Element startingElement) {
		Deque<Element> stack = new ArrayDeque<Element>();
		stack.add(startingElement.getParent());
		List<Element> foundGroups = new ArrayList<Element>();
		boolean doneFirstIteration = false;//check on index only done on first iteration to only get elements with an index greater than the starting element
		while (stack.size() > 0) {
			Element currentElement =stack.removeLast();
			if (currentElement.getName().equals(GROUP_EL)) {
				foundGroups.add(currentElement);
				continue;
			}
			List<Element> siblings = OpsinTools.getChildElementsWithTagNames(currentElement, new String[]{BRACKET_EL, SUBSTITUENT_EL, ROOT_EL});

			List<Element> bracketted = new ArrayList<Element>();
			for (Element bracketOrSubOrRoot : siblings) {
				if (!doneFirstIteration && currentElement.indexOf(bracketOrSubOrRoot) <= currentElement.indexOf(startingElement)){
					continue;
				}
				if (bracketOrSubOrRoot.getAttribute(MULTIPLIER_ATR) != null){
					continue;
				}
				if (bracketOrSubOrRoot.getName().equals(BRACKET_EL)){
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
			//locanting into brackets is rarely the desired answer so place at the bottom of the stack
			for (int i = bracketted.size() -1; i >=0; i--) {
				stack.addFirst(bracketted.get(i));
			}
			doneFirstIteration = true;
		}
		return foundGroups;
	}

	/**
	 * Checks through the groups accessible from the currentElement taking into account brackets
	 * i.e. those that it is feasible that the group of the currentElement could substitute onto
	 * @param startingElement
	 * @param locant: the locant string to check for the presence of
	 * @return The fragment with the locant, or null
	 * @throws StructureBuildingException
	 */
	private static Fragment findFragmentWithLocant(Element startingElement, String locant) throws StructureBuildingException {
		Deque<Element> stack = new ArrayDeque<Element>();
		stack.add(startingElement.getParent());
		boolean doneFirstIteration = false;//check on index only done on first iteration to only get elements with an index greater than the starting element
		Fragment monoNuclearHydride = null;//e.g. methyl/methane - In this case no locant would be expected as unlocanted substitution is always unambiguous. Hence deprioritise
		while (stack.size() > 0) {
			Element currentElement = stack.removeLast();
			if (currentElement.getName().equals(SUBSTITUENT_EL) || currentElement.getName().equals(ROOT_EL)) {
				Fragment groupFrag = currentElement.getFirstChildElement(GROUP_EL).getFrag();
				if (monoNuclearHydride != null && currentElement.getAttribute(LOCANT_ATR) != null) {//It looks like all groups are locanting onto the monoNuclearHydride e.g. 1-oxo-1-phenyl-sulfanylidene
					return monoNuclearHydride;
				}
				if (groupFrag.hasLocant(locant)) {
					if (locant.equals("1") && groupFrag.getAtomCount() == 1) {
						if (monoNuclearHydride == null) {
							monoNuclearHydride = groupFrag;
						}
					}
					else{
						return groupFrag;
					}
				}
				continue;
			}
			else if (monoNuclearHydride != null) {
				return monoNuclearHydride;
			}
			List<Element> siblings = OpsinTools.getChildElementsWithTagNames(currentElement, new String[]{BRACKET_EL, SUBSTITUENT_EL, ROOT_EL});

			List<Element> bracketted = new ArrayList<Element>();
			if (!doneFirstIteration) {//on the first iteration, ignore elements before the starting element and favour the element directly after the starting element (conditions apply)
				int indexOfStartingEl = currentElement.indexOf(startingElement);
				Element substituentToTryFirst = null;
				for (Element bracketOrSubOrRoot : siblings) {
					int indexOfCurrentEl = currentElement.indexOf(bracketOrSubOrRoot);
					if (indexOfCurrentEl <= indexOfStartingEl) {
						continue;
					}
					if (bracketOrSubOrRoot.getAttribute(MULTIPLIER_ATR) != null) {
						continue;
					}
					
					if (bracketOrSubOrRoot.getName().equals(BRACKET_EL)) {
						if (IMPLICIT_TYPE_VAL.equals(bracketOrSubOrRoot.getAttributeValue(TYPE_ATR)) && bracketOrSubOrRoot.getAttribute(LOCANT_EL) == null) {
							//treat implicit brackets without locants as if they are not there
							for (Element descendent : getChildrenIgnoringLocantlessImplicitBrackets(bracketOrSubOrRoot)) {
								if (descendent.getName().equals(BRACKET_EL)) {
									bracketted.add(descendent);
								}
								else {
									if (substituentToTryFirst == null && descendent.getAttribute(LOCANT_EL) == null && MATCH_NUMERIC_LOCANT.matcher(locant).matches()) {
										substituentToTryFirst = descendent;
									}
									else {
										stack.add(descendent);
									}
								}
							}
						}
						else {
							bracketted.add(bracketOrSubOrRoot);
						}
					}
					else {
						if (substituentToTryFirst == null && bracketOrSubOrRoot.getAttribute(LOCANT_EL) == null && MATCH_NUMERIC_LOCANT.matcher(locant).matches()) {
							substituentToTryFirst = bracketOrSubOrRoot;
						}
						else {
							stack.add(bracketOrSubOrRoot);
						}
					}
				}
				if (substituentToTryFirst != null) {
					stack.add(substituentToTryFirst);
				}
				doneFirstIteration = true;
			}
			else {
				for (Element bracketOrSubOrRoot : siblings) {
					if (bracketOrSubOrRoot.getAttribute(MULTIPLIER_ATR) != null) {
						continue;
					}
					if (bracketOrSubOrRoot.getName().equals(BRACKET_EL)) {
						if (IMPLICIT_TYPE_VAL.equals(bracketOrSubOrRoot.getAttributeValue(TYPE_ATR)) && bracketOrSubOrRoot.getAttribute(LOCANT_EL) == null) {
							//treat implicit brackets without locants as if they are not there
							for (Element descendent : getChildrenIgnoringLocantlessImplicitBrackets(bracketOrSubOrRoot)) {
								if (descendent.getName().equals(BRACKET_EL)) {
									bracketted.add(descendent);
								}
								else {
									stack.add(descendent);
								}
							}
						}
						else {
							bracketted.add(bracketOrSubOrRoot);
						}
					}
					else {
						stack.add(bracketOrSubOrRoot);
					}
				}
			}
			//locanting into brackets is rarely the desired answer so place at the bottom of the stack
			for (int i = bracketted.size() -1; i >=0; i--) {
				stack.addFirst(bracketted.get(i));
			}
		}
		return monoNuclearHydride;
	}

	private static List<Element> getChildrenIgnoringLocantlessImplicitBrackets(Element implicitBracket) {
		List<Element> childrenAndImplicitBracketChildren = new ArrayList<Element>();
		for (Element child : implicitBracket.getChildElements()) {
			if (child.getName().equals(BRACKET_EL) && IMPLICIT_TYPE_VAL.equals(child.getAttributeValue(TYPE_ATR)) && child.getAttribute(LOCANT_EL) == null) {
				childrenAndImplicitBracketChildren.addAll(getChildrenIgnoringLocantlessImplicitBrackets(child));
			}
			else {
				childrenAndImplicitBracketChildren.add(child);
			}
		}
		return childrenAndImplicitBracketChildren;
	}

	static Element findRightMostGroupInBracket(Element bracket) {
		List<Element> subsBracketsAndRoots = OpsinTools.getChildElementsWithTagNames(bracket, new String[]{BRACKET_EL, SUBSTITUENT_EL, ROOT_EL});
		Element lastSubsBracketOrRoot = subsBracketsAndRoots.get(subsBracketsAndRoots.size() - 1);
		while (lastSubsBracketOrRoot.getName().equals(BRACKET_EL)) {
			subsBracketsAndRoots = OpsinTools.getChildElementsWithTagNames(lastSubsBracketOrRoot, new String[]{BRACKET_EL, SUBSTITUENT_EL, ROOT_EL});
			lastSubsBracketOrRoot = subsBracketsAndRoots.get(subsBracketsAndRoots.size() - 1);
		}
		return findRightMostGroupInSubOrRoot(lastSubsBracketOrRoot);
	}
	
	static Element findRightMostGroupInSubBracketOrRoot(Element subBracketOrRoot) {
		if (subBracketOrRoot.getName().equals(BRACKET_EL)) {
			return findRightMostGroupInBracket(subBracketOrRoot);
		}
		else {
			return findRightMostGroupInSubOrRoot(subBracketOrRoot);
		}
	}

	private static Element findRightMostGroupInSubOrRoot(Element subOrRoot) {
		for (int i = subOrRoot.getChildCount() - 1; i >= 0; i--) {
			Element el = subOrRoot.getChild(i);
			if (el.getName().equals(GROUP_EL)) {
				return el;
			}
		}
		return null;
	}

	private static boolean potentiallyCanSubstitute(Element subBracketOrRoot) {
		Element parent = subBracketOrRoot.getParent();
		List<Element> children =parent.getChildElements();
		for (int i = parent.indexOf(subBracketOrRoot) +1 ; i < children.size(); i++) {
			if (!children.get(i).getName().equals(HYPHEN_EL)){
				return true;
			}
		}
		return false;
	}

	static String checkForBracketedPrimedLocantSpecialCase(Element subBracketOrRoot, String locantString) {
		int terminalPrimes = StringTools.countTerminalPrimes(locantString);
		if (terminalPrimes > 0){
			int brackettingDepth = 0;
			Element parent = subBracketOrRoot.getParent();
			while (parent != null && parent.getName().equals(BRACKET_EL)){
				if (!IMPLICIT_TYPE_VAL.equals(parent.getAttributeValue(TYPE_ATR))){
					brackettingDepth++;
				}
				parent = parent.getParent();
			}
			if (terminalPrimes == brackettingDepth){
				return locantString.substring(0, locantString.length() - terminalPrimes);
			}
		}
		return null;
	}

	/**
	 * In cases such as methylenecyclohexane two outAtoms are combined to form a single outAtom with valency
	 * equal to sum of the valency of the other outAtoms.
	 * This is only allowed on substituents where all the outAtoms are on the same atom
	 * @param frag
	 * @param group 
	 * @throws StructureBuildingException
	 */
	private static void checkAndApplySpecialCaseWhereOutAtomsCanBeCombinedOrThrow(Fragment frag, Element group) throws StructureBuildingException {
		int outAtomCount = frag.getOutAtomCount();
		if (outAtomCount <= 1) {
			return;
		}
		if (EPOXYLIKE_SUBTYPE_VAL.equals(group.getAttributeValue(SUBTYPE_ATR))){
			return;
		}
		String groupValue = group.getValue();
		if (groupValue.equals("oxy") || groupValue.equals("thio") || groupValue.equals("seleno") || groupValue.equals("telluro")){//always bivalent
			return;
		}
		//special case: all outAtoms on same atom e.g. methylenecyclohexane
		Atom firstOutAtom = frag.getOutAtom(0).getAtom();
		int valencyOfOutAtom = 0;
		for (int i = outAtomCount - 1; i >=0 ; i--) {//remove all outAtoms and add one with the total valency of all those that have been removed
			OutAtom out = frag.getOutAtom(i);
			if (!out.getAtom().equals(firstOutAtom)){
				throw new StructureBuildingException("Substitutive bond formation failure: Fragment expected to have one OutAtom but had: "+ outAtomCount);
			}
			valencyOfOutAtom += out.getValency();
			frag.removeOutAtom(i);
		}
		frag.addOutAtom(firstOutAtom, valencyOfOutAtom, true);
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
		if (!atom.getImplicitHydrogenAllowed()) {
			return 0;
		}
		int valency = atom.determineValency(true);
		int currentValency = atom.getIncomingValency() + atom.getOutValency();
		int substitutableHydrogen = valency - currentValency;
		return substitutableHydrogen >= 0 ? substitutableHydrogen : 0;
	}

	/**
	 * Stereochemistry terms are assigned right at the end so that checks can be done on whether the indicated atom is in fact chiral.
	 * In the process of multiplication locants are primed. This function adds the appropriate number of primes to any locanted stereochemistry locants
	 * The primesString is the string containing the primes to add to each locant
	 * @param subOrBracket
	 * @param primesString
	 */
	private static void addPrimesToLocantedStereochemistryElements(Element subOrBracket, String primesString) {
		List<Element> stereoChemistryElements =OpsinTools.getDescendantElementsWithTagName(subOrBracket, STEREOCHEMISTRY_EL);
		for (Element stereoChemistryElement : stereoChemistryElements) {
			if (stereoChemistryElement.getAttribute(LOCANT_ATR) != null){
				stereoChemistryElement.getAttribute(LOCANT_ATR).setValue(stereoChemistryElement.getAttributeValue(LOCANT_ATR) + primesString);
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
		while (!element.getName().equals(WORD_EL)){
			element = element.getParent();
			if (element == null){
				return null;
			}
			count++;
		}	
		return count;
	}
}
