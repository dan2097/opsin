package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;

import uk.ac.cam.ch.wwmm.opsin.BondStereo.BondStereoValue;
import uk.ac.cam.ch.wwmm.opsin.ParseWord.WordType;
import uk.ac.cam.ch.wwmm.opsin.StereoAnalyser.StereoBond;
import uk.ac.cam.ch.wwmm.opsin.StereoAnalyser.StereoCentre;

import static uk.ac.cam.ch.wwmm.opsin.XmlDeclarations.*;

import nu.xom.Element;

/**
 * Identifies stereocentres, assigns stereochemistry elements to them and then uses the CIP rules to calculate appropriates atomParity/bondstereo tags
 * @author dl387
 *
 */
class StereochemistryHandler {

	private final static Pattern matchSlash = Pattern.compile("/");
	
	/**
	 * Master method for assigning and processing stereochemistry elements
	 * @param state
	 * @param uniFrag
	 * @param stereoChemistryEls
	 * @throws StructureBuildingException
	 */
	static void processStereochemicalElements(BuildState state, Fragment uniFrag, List<Element> stereoChemistryEls) throws StructureBuildingException {
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
			if (notIn6MemberOrSmallerRing(b)){
				bondStereoBondMap.put(b, stereoBond);
			}
		}
		List<Element> locantedStereoChemistryEls = new ArrayList<Element>();
		List<Element> unlocantedStereoChemistryEls = new ArrayList<Element>();
		for (Element stereoChemistryElement : stereoChemistryEls) {
			if (stereoChemistryElement.getAttributeValue(LOCANT_ATR)!=null){
				locantedStereoChemistryEls.add(stereoChemistryElement);
			}
			else{
				unlocantedStereoChemistryEls.add(stereoChemistryElement);
			}
		}
		//perform locanted before unlocanted to avoid unlocanted elements using the stereocentres a locanted element refers to
		matchStereochemistryToAtomsAndBonds(state, locantedStereoChemistryEls, atomStereoCentreMap, bondStereoBondMap);
		matchStereochemistryToAtomsAndBonds(state, unlocantedStereoChemistryEls, atomStereoCentreMap, bondStereoBondMap);
	}


	/**
	 * Checks whether a bond is part of a 6 member or smaller ring.
	 * This is necessary as such double bonds are assumed to not be capable of having E/Z stereochemistry
	 * @param bond
	 * @return true unless in a 6 member or smaller rings
	 */
	private static boolean notIn6MemberOrSmallerRing(Bond bond) {
		Atom fromAtom =bond.getFromAtom();
		Atom toAtom = bond.getToAtom();
		if (fromAtom.getAtomIsInACycle() && toAtom.getAtomIsInACycle()){//obviously both must be in rings
			//attempt to get from the fromAtom to the toAtom in 6 or fewer steps.
			List<Atom> visitedAtoms = new ArrayList<Atom>();
			LinkedList<Atom> atomsToInvestigate = new LinkedList<Atom>();//A queue is not used as I need to make sure that only up to depth 6 is investigated
			List<Atom> neighbours =fromAtom.getAtomNeighbours();
			neighbours.remove(toAtom);
			for (Atom neighbour : neighbours) {
				atomsToInvestigate.add(neighbour);
			}
			visitedAtoms.add(fromAtom);
			for (int i = 0; i < 5; i++) {//up to 5 bonds from the neighbours of the fromAtom i.e. up to ring size 6
				if (atomsToInvestigate.isEmpty()){
					break;
				}
				LinkedList<Atom> atomsToInvestigateNext = new LinkedList<Atom>();
				while (!atomsToInvestigate.isEmpty()) {
					Atom currentAtom =atomsToInvestigate.removeFirst();
					if (currentAtom == toAtom){
						return false;
					}
					visitedAtoms.add(currentAtom);
					neighbours =currentAtom.getAtomNeighbours();
					for (Atom neighbour : neighbours) {
						if (!visitedAtoms.contains(neighbour) && neighbour.getAtomIsInACycle()){
							atomsToInvestigateNext.add(neighbour);
						}
					}
				}
				atomsToInvestigate = atomsToInvestigateNext;
			}
		}
		return true;
	}


	/**
	 * Finds all stereochemistry elements and attempts to locate a suitable atom/bond for them
	 * Currently only locanted stereochemistry elements are supported
	 * @param state
	 * @param stereoChemistryEls A list of stereochemistry elements
	 * @param atomStereoCentreMap
	 * @param bondStereoBondMap
	 * @throws StructureBuildingException
	 */
	private static void matchStereochemistryToAtomsAndBonds(BuildState state, List<Element> stereoChemistryEls, Map<Atom, StereoCentre> atomStereoCentreMap, Map<Bond, StereoBond> bondStereoBondMap) throws StructureBuildingException {
		for (Element stereoChemistryEl : stereoChemistryEls) {
			String stereoChemistryType =stereoChemistryEl.getAttributeValue(TYPE_ATR);
			if (stereoChemistryType.equals(R_OR_S_TYPE_VAL)){
				assignStereoCentre(state, stereoChemistryEl, atomStereoCentreMap);
			}
			else if (stereoChemistryType.equals(E_OR_Z_TYPE_VAL)){
				assignStereoBond(state, stereoChemistryEl, bondStereoBondMap);
			}
			else if (stereoChemistryType.equals(CISORTRANS_TYPE_VAL)){
				if (!assignCisTransOnRing(state, stereoChemistryEl, atomStereoCentreMap)){
					assignStereoBond(state, stereoChemistryEl, bondStereoBondMap);
				}
			}
			else if (stereoChemistryType.equals(ALPHA_OR_BETA_TYPE_VAL)){
				assignAlphaBetaStereochem(state, stereoChemistryEl, atomStereoCentreMap);
			}
			else{
				throw new StructureBuildingException("Unsupported stereochemistry type: " +stereoChemistryType);
			}
			stereoChemistryEl.detach();
		}

	}

	/**
	 * Handles R/S stereochemistry. r/s is not currently handled
	 * @param state
	 * @param stereoChemistryEl
	 * @param atomStereoCentreMap
	 * @throws StructureBuildingException
	 */
	private static void assignStereoCentre(BuildState state,Element stereoChemistryEl, Map<Atom, StereoCentre> atomStereoCentreMap) throws StructureBuildingException {
		//generally the LAST group in this list will be the appropriate groups e.g. (5S)-5-ethyl-6-methylheptane where the heptane is the appropriate group
		//we use the same algorithm as for unlocanted substitution so as to deprecate assignment into brackets
		Element parentSubBracketOrRoot = (Element) stereoChemistryEl.getParent();
		List<Fragment> possibleFragments = StructureBuildingMethods.findAlternativeFragments(state, parentSubBracketOrRoot);
		List<Element> adjacentGroupEls = XOMTools.getDescendantElementsWithTagName(parentSubBracketOrRoot, GROUP_EL);
		for (int i = adjacentGroupEls.size()-1; i >=0; i--) {
			possibleFragments.add(state.xmlFragmentMap.get(adjacentGroupEls.get(i)));
		}
		String locant = stereoChemistryEl.getAttributeValue(LOCANT_ATR);
		String rOrS = stereoChemistryEl.getAttributeValue(VALUE_ATR);
		for (Fragment fragment : possibleFragments) {
			if (locant ==null){//undefined locant
				List<Atom> atomList = fragment.getAtomList();
				for (Atom potentialStereoAtom : atomList) {
					if (atomStereoCentreMap.containsKey(potentialStereoAtom)){
						applyStereoChemistryToStereoCentre(potentialStereoAtom, atomStereoCentreMap.get(potentialStereoAtom), rOrS);
						atomStereoCentreMap.remove(potentialStereoAtom);
						return;
					}
				}
			}
			else{
				Atom potentialStereoAtom = fragment.getAtomByLocant(locant);
				if (potentialStereoAtom !=null && atomStereoCentreMap.containsKey(potentialStereoAtom)){
					applyStereoChemistryToStereoCentre(potentialStereoAtom, atomStereoCentreMap.get(potentialStereoAtom), rOrS);
					atomStereoCentreMap.remove(potentialStereoAtom);
					return;
				}
			}
		}
		Element possibleWordParent = (Element) parentSubBracketOrRoot.getParent();
		if (possibleWordParent.getLocalName().equals(WORD_EL) && possibleWordParent.getAttributeValue(TYPE_ATR).equals(WordType.substituent.toString())){
			//something like (3R,4R,5R)-ethyl 4-acetamido-5-amino-3-(pentan-3-yloxy)cyclohex-1-enecarboxylate
			//I think this is a violation of the IUPAC rules...but anyway...
			List<Element> words = XOMTools.getChildElementsWithTagNameAndAttribute(((Element)possibleWordParent.getParent()), WORD_EL, TYPE_ATR, WordType.full.toString());
			for (Element word : words) {
				List<Element> possibleGroups = XOMTools.getDescendantElementsWithTagName(word, GROUP_EL);
				for (int i = possibleGroups.size()-1; i >=0; i--) {
					Fragment correspondingFrag = state.xmlFragmentMap.get(possibleGroups.get(i));
					if (locant == null){//undefined locant
						List<Atom> atomList = correspondingFrag.getAtomList();
						for (Atom potentialStereoAtom : atomList) {
							if (atomStereoCentreMap.containsKey(potentialStereoAtom)){
								applyStereoChemistryToStereoCentre(potentialStereoAtom, atomStereoCentreMap.get(potentialStereoAtom), rOrS);
								atomStereoCentreMap.remove(potentialStereoAtom);
								return;
							}
						}
					}
					else{
						Atom potentialStereoAtom = correspondingFrag.getAtomByLocant(locant);
						if (potentialStereoAtom !=null && atomStereoCentreMap.containsKey(potentialStereoAtom)){
							applyStereoChemistryToStereoCentre(potentialStereoAtom, atomStereoCentreMap.get(potentialStereoAtom), rOrS);
							atomStereoCentreMap.remove(potentialStereoAtom);
							return;
						}
					}
				}
			}
		}
		throw new StructureBuildingException("Could not find atom that: " + stereoChemistryEl.toXML() + " appeared to be referring to");
	}


	/**
	 * Assigns atom parity to the given atom in accordance with the CIP rules
	 * @param atom The stereoAtom
	 * @param stereoCentre
	 * @param rOrS The description given in the name
	 * @throws StructureBuildingException
	 */
	private static void applyStereoChemistryToStereoCentre(Atom atom, StereoCentre stereoCentre, String rOrS) throws StructureBuildingException {
		List<Atom> cipOrderedAtoms =stereoCentre.getCipOrderedAtoms();
		if (cipOrderedAtoms.size()!=4){
			throw new StructureBuildingException("Only tetrahedral chirality is currently supported");
		}
		Atom[] atomRefs4 = new Atom[4];
		atomRefs4[0] = cipOrderedAtoms.get(cipOrderedAtoms.size()-1);
		for (int i = 0; i < cipOrderedAtoms.size() -1; i++) {//from highest to lowest (true for S) hence atomParity 1 for S
			atomRefs4[i+1] = cipOrderedAtoms.get(i);
		}
		if (rOrS.equals("R")){
			atom.setAtomParity(atomRefs4, -1);
		}
		else if (rOrS.equals("S")){
			atom.setAtomParity(atomRefs4, 1);
		}
		else{
			throw new StructureBuildingException("Unexpected stereochemistry type: " + rOrS);
		}
	}


	/**
	 * Handles E/Z stereochemistry and cis/trans in cases where cis/trans unambiguously corresponds to E/Z
	 * @param state
	 * @param stereoChemistryEl
	 * @param bondStereoBondMap
	 * @throws StructureBuildingException
	 */
	private static void assignStereoBond(BuildState state, Element stereoChemistryEl, Map<Bond, StereoBond> bondStereoBondMap) throws StructureBuildingException {
		//generally the LAST group in this list will be the appropriate groups e.g. (2Z)-5-ethyl-6-methylhex-2-ene where the hex-2-ene is the appropriate group
		//we use the same algorithm as for unlocanted substitution so as to deprecate assignment into brackets
		Element parentSubBracketOrRoot = (Element) stereoChemistryEl.getParent();
		List<Fragment> possibleFragments = StructureBuildingMethods.findAlternativeFragments(state, parentSubBracketOrRoot);
		List<Element> adjacentGroupEls = XOMTools.getDescendantElementsWithTagName(parentSubBracketOrRoot, GROUP_EL);
		for (int i = adjacentGroupEls.size()-1; i >=0; i--) {
			possibleFragments.add(state.xmlFragmentMap.get(adjacentGroupEls.get(i)));
		}
		String locant = stereoChemistryEl.getAttributeValue(LOCANT_ATR);
		String eOrZ = stereoChemistryEl.getAttributeValue(VALUE_ATR);
		boolean isCisTrans =false;
		if (stereoChemistryEl.getAttributeValue(TYPE_ATR).equals(CISORTRANS_TYPE_VAL)){
			isCisTrans =true;
			String cisOrTrans = stereoChemistryEl.getAttributeValue(VALUE_ATR);
			if (cisOrTrans.equalsIgnoreCase("cis")){
				eOrZ = "Z";
			}
			else if (cisOrTrans.equalsIgnoreCase("trans")){
				eOrZ = "E";
			}
			else{
				throw new StructureBuildingException("Unexpected cis/trans stereochemistry type: " +cisOrTrans);
			}
		}
		for (Fragment fragment : possibleFragments) {
			if (locant == null){//undefined locant
				Set<Bond> bondSet = fragment.getBondSet();
				for (Bond potentialBond : bondSet) {
					if (bondStereoBondMap.containsKey(potentialBond) && (!isCisTrans || cisTranUnambiguousOnBond(potentialBond))){
						applyStereoChemistryToStereoBond(potentialBond, bondStereoBondMap.get(potentialBond), eOrZ);
						bondStereoBondMap.remove(potentialBond);
						return;
					}
				}
				List<Bond> sortedInterFragmentBonds = sortInterFragmentBonds(state.fragManager.getInterFragmentBonds(fragment), fragment);
				for (Bond potentialBond : sortedInterFragmentBonds) {
					if (bondStereoBondMap.containsKey(potentialBond) && (!isCisTrans || cisTranUnambiguousOnBond(potentialBond))){
						applyStereoChemistryToStereoBond(potentialBond, bondStereoBondMap.get(potentialBond), eOrZ);
						bondStereoBondMap.remove(potentialBond);
						return;
					}
				}
			}
			else{
				Atom firstAtomInBond = fragment.getAtomByLocant(locant);
				if (firstAtomInBond !=null){
					Set<Bond> bonds = firstAtomInBond.getBonds();
					for (Bond potentialBond : bonds) {
						if (bondStereoBondMap.containsKey(potentialBond) && (!isCisTrans || cisTranUnambiguousOnBond(potentialBond))){
							applyStereoChemistryToStereoBond(potentialBond, bondStereoBondMap.get(potentialBond), eOrZ);
							bondStereoBondMap.remove(potentialBond);
							return;
						}
					}
				}
			}
		}
		Element possibleWordParent = (Element) parentSubBracketOrRoot.getParent();
		if (possibleWordParent.getLocalName().equals(WORD_EL) && possibleWordParent.getAttributeValue(TYPE_ATR).equals(WordType.substituent.toString())){
			//the element is in front of a substituent and may refer to the full group
			List<Element> words = XOMTools.getChildElementsWithTagNameAndAttribute(((Element)possibleWordParent.getParent()), WORD_EL, TYPE_ATR, WordType.full.toString());
			for (Element word : words) {
				List<Element> possibleGroups = XOMTools.getDescendantElementsWithTagName(word, GROUP_EL);
				for (int i = possibleGroups.size()-1; i >=0; i--) {
					Fragment correspondingFrag = state.xmlFragmentMap.get(possibleGroups.get(i));
					if (locant == null){//undefined locant
						Set<Bond> bondSet = correspondingFrag.getBondSet();
						for (Bond potentialBond : bondSet) {
							if (bondStereoBondMap.containsKey(potentialBond) && (!isCisTrans || cisTranUnambiguousOnBond(potentialBond))){
								applyStereoChemistryToStereoBond(potentialBond, bondStereoBondMap.get(potentialBond), eOrZ);
								bondStereoBondMap.remove(potentialBond);
								return;
							}
						}
						List<Bond> sortedInterFragmentBonds = sortInterFragmentBonds(state.fragManager.getInterFragmentBonds(correspondingFrag), correspondingFrag);
						for (Bond potentialBond : sortedInterFragmentBonds) {
							if (bondStereoBondMap.containsKey(potentialBond) && (!isCisTrans || cisTranUnambiguousOnBond(potentialBond))){
								applyStereoChemistryToStereoBond(potentialBond, bondStereoBondMap.get(potentialBond), eOrZ);
								bondStereoBondMap.remove(potentialBond);
								return;
							}
						}
					}
					else{
						Atom firstAtomInBond = correspondingFrag.getAtomByLocant(locant);
						if (firstAtomInBond !=null){
							Set<Bond> bonds = firstAtomInBond.getBonds();
							for (Bond potentialBond : bonds) {
								if (bondStereoBondMap.containsKey(potentialBond) && (!isCisTrans || cisTranUnambiguousOnBond(potentialBond))){
									applyStereoChemistryToStereoBond(potentialBond, bondStereoBondMap.get(potentialBond), eOrZ);
									bondStereoBondMap.remove(potentialBond);
									return;
								}
							}
						}
					}
				}
			}
		}
		if (isCisTrans){
			throw new StructureBuildingException("Could not find bond that: " + stereoChemistryEl.toXML() + " could refer unambiguously to");
		}
		else{
			throw new StructureBuildingException("Could not find bond that: " + stereoChemistryEl.toXML() + " was referring to");
		}
	}


	/**
	 * Does the stereoBond have a hydrogen connected to both ends of it.
	 * If not it is ambiguous when used in conjunction with cis/trans and E/Z should be used.
	 * @param potentialBond
	 * @return
	 */
	static boolean cisTranUnambiguousOnBond(Bond potentialBond) {
		List<Atom> neighbours1 = potentialBond.getFromAtom().getAtomNeighbours();
		boolean foundHydrogen1 =false;
		for (Atom neighbour : neighbours1) {
			if (neighbour.getElement().equals("H")){
				foundHydrogen1 =true;
			}
		}
		
		List<Atom> neighbours2 = potentialBond.getToAtom().getAtomNeighbours();
		boolean foundHydrogen2 =false;
		for (Atom neighbour : neighbours2) {
			if (neighbour.getElement().equals("H")){
				foundHydrogen2 =true;
			}
		}
		return (foundHydrogen1 && foundHydrogen2);
	}


	/**
	 * Sorts bonds such that those originating from the given fragment are preferred
	 * @param interFragmentBonds A set of interFragment Bonds
	 * @param preferredOriginatingFragment 
	 * @return A sorted list
	 */
	private static List<Bond> sortInterFragmentBonds(Set<Bond> interFragmentBonds, Fragment preferredOriginatingFragment) {
		List<Bond> interFragmentBondList = new ArrayList<Bond>();
		for (Bond bond : interFragmentBonds) {
			if (bond.getFromAtom().getFrag() ==preferredOriginatingFragment){
				interFragmentBondList.add(0, bond);
			}
			else{
				interFragmentBondList.add(bond);
			}
		}
		return interFragmentBondList;
	}


	/**
	 * Assigns bondstereo to the given bond in accordance with the CIP rules
	 * @param bond The stereobond
	 * @param stereoBond
	 * @param eOrZ The stereo description given in the name
	 * @throws StructureBuildingException
	 */
	private static void applyStereoChemistryToStereoBond(Bond bond, StereoBond stereoBond, String eOrZ ) throws StructureBuildingException {
		List<Atom> stereoBondAtoms = stereoBond.getOrderedStereoAtoms();
		//stereoBondAtoms contains the higher priority atom at one end, the two bond atoms and the higher priority atom at the other end
		Atom[] atomRefs4 = new Atom[4];
		atomRefs4[0] = stereoBondAtoms.get(0);
		atomRefs4[1] = stereoBondAtoms.get(1);
		atomRefs4[2] = stereoBondAtoms.get(2);
		atomRefs4[3] = stereoBondAtoms.get(3);
		if (eOrZ.equals("E")){
			bond.setBondStereoElement(atomRefs4, BondStereoValue.TRANS);
		}
		else if (eOrZ.equals("Z")){
			bond.setBondStereoElement(atomRefs4, BondStereoValue.CIS);
		}
		else{
			throw new StructureBuildingException("Unexpected stereochemistry type: " + eOrZ);
		}
	}
	
	/**
	 * Searches for instances of two tetrahedral stereocentres/psuedo-stereocentres
	 * then sets their configuration such that the substituents at these centres are cis or trans to each other
	 * @param state
	 * @param stereoChemistryEl
	 * @param atomStereoCentreMap
	 * @return
	 * @throws StructureBuildingException
	 */
	private static boolean assignCisTransOnRing(BuildState state, Element stereoChemistryEl, Map<Atom, StereoCentre> atomStereoCentreMap) throws StructureBuildingException {
		if (stereoChemistryEl.getAttribute(LOCANT_ATR)!=null){
			return false;
		}
		Element parentSubBracketOrRoot = (Element) stereoChemistryEl.getParent();
		List<Fragment> possibleFragments = StructureBuildingMethods.findAlternativeFragments(state, parentSubBracketOrRoot);
		List<Element> adjacentGroupEls = XOMTools.getDescendantElementsWithTagName(parentSubBracketOrRoot, GROUP_EL);
		for (int i = adjacentGroupEls.size()-1; i >=0; i--) {
			possibleFragments.add(state.xmlFragmentMap.get(adjacentGroupEls.get(i)));
		}
		for (Fragment fragment : possibleFragments) {
			List<Atom> atomList = fragment.getAtomList();
			List<Atom> stereoAtoms = new ArrayList<Atom>();
			for (Atom potentialStereoAtom : atomList) {
				if (potentialStereoAtom.getAtomIsInACycle()){
					List<Atom> neighbours = potentialStereoAtom.getAtomNeighbours();
					if (neighbours.size()==4){
						int hydrogenCount =0;
						int acylicOrNotInFrag =0;
						for (Atom neighbour : neighbours) {
							if (neighbour.getElement().equals("H")){
								hydrogenCount++;
							}
							if (!neighbour.getAtomIsInACycle() || !atomList.contains(neighbour)){
								acylicOrNotInFrag++;
							}
						}
						if (hydrogenCount==1 || (hydrogenCount==0 && acylicOrNotInFrag ==1) ){
							stereoAtoms.add(potentialStereoAtom);
						}
					}
				}
			}
			if (stereoAtoms.size()==2){
				Atom a1 = stereoAtoms.get(0);
				Atom a2 = stereoAtoms.get(1);
				
				List<List<Atom>> paths = CycleDetector.getIntraFragmentPathsBetweenAtoms(a1, a2, fragment);
				paths =findNonOverlappingPaths(paths);
				if (paths.size()!=2 && paths.size()!=3){
					return false;
				}
				if (a1.getAtomParity()!=null && a2.getAtomParity()!=null){//one can have defined stereochemistry but not both
					return false;
				}
				applyStereoChemistryToCisTransOnRing(a1, a2, paths, stereoChemistryEl.getAttributeValue(VALUE_ATR));
				atomStereoCentreMap.remove(stereoAtoms.get(0));
				atomStereoCentreMap.remove(stereoAtoms.get(1));
				return true;
			}
		}
		return false;
	}


	/**
	 * Tries to find the paths that do not overlap.
	 * At least one path for each starting atom must be picked.
	 * If this is not possible no paths are returned
	 * @param paths
	 * @return
	 */
	private static List<List<Atom>> findNonOverlappingPaths(List<List<Atom>> paths) {
		Map<Atom, List<List<Atom>>> sameStartingAtom = new LinkedHashMap<Atom, List<List<Atom>>>();
		List<List<Atom>> zeroLengthPaths = new ArrayList<List<Atom>>();
		for (List<Atom> path : paths) {
			if (path.size()>0){
				Atom firstAtom =path.get(0);
				if (sameStartingAtom.get(firstAtom)==null){
					sameStartingAtom.put(firstAtom, new ArrayList<List<Atom>>());
				}
				sameStartingAtom.get(firstAtom).add(path);
			}
			else{
				zeroLengthPaths.add(path);
			}
		}
		Atom firstKey =sameStartingAtom.keySet().iterator().next();
		List<List<Atom>> possiblePaths = sameStartingAtom.get(firstKey);
		sameStartingAtom.remove(firstKey);
		for (List<Atom> possiblePath : possiblePaths) {
			List<List<Atom>> nonOverlappingPaths = new LinkedList<List<Atom>>();
			nonOverlappingPaths.add(possiblePath);
			for (Atom otherStartingAtom : sameStartingAtom.keySet()) {
				List<List<Atom>> otherPaths = sameStartingAtom.get(otherStartingAtom);
				boolean nonOverlappingPathFound =false;
				for (List<Atom> otherPath : otherPaths) {
					boolean match =false;
					ListLoop: for (List<Atom> nonOverlappingPath : nonOverlappingPaths) {
						for (Atom atom : nonOverlappingPath) {
							if (otherPath.contains(atom)){
								match =true;
								break ListLoop;
							}
						}
					}
					if (!match){
						nonOverlappingPathFound=true;
						nonOverlappingPaths.add(otherPath);
						break;
					}
				}
				if (!nonOverlappingPathFound){
					break;
				}
			}
			if (nonOverlappingPaths.size()==sameStartingAtom.keySet().size()+1){
				nonOverlappingPaths.addAll(0, zeroLengthPaths);
				return nonOverlappingPaths;
			}
		}
		return Collections.emptyList();
	}


	private static void applyStereoChemistryToCisTransOnRing(Atom a1, Atom a2, List<List<Atom>> paths, String cisOrTrans) throws StructureBuildingException {
		List<Atom> a1Neighbours = a1.getAtomNeighbours();
		Atom[] atomRefs4a1 = new Atom[4];
		Atom firstPathAtom = paths.get(0).size()>0 ? paths.get(0).get(0) : a2;
		atomRefs4a1[2] = firstPathAtom;
		Atom secondPathAtom = paths.get(1).size()>0 ? paths.get(1).get(0) : a2;
		atomRefs4a1[3] = secondPathAtom;
		a1Neighbours.remove(firstPathAtom);
		a1Neighbours.remove(secondPathAtom);
		if (paths.size()==3){
			atomRefs4a1[1] = paths.get(2).size()>0 ? paths.get(2).get(0) : a2;
		}
		else{
			for (Atom atom : a1Neighbours) {
				if (atom.getElement().equals("H")){
					atomRefs4a1[1] = atom;
					break;
				}
			}
			if (atomRefs4a1[1] ==null){
				throw new StructureBuildingException("OPSIN Bug: cannot assign cis/trans on ring stereochemistry");
			}
		}
		a1Neighbours.remove(atomRefs4a1[1]);
		atomRefs4a1[0] = a1Neighbours.get(0);
		
		
		List<Atom> a2Neighbours = a2.getAtomNeighbours();
		Atom[] atomRefs4a2 = new Atom[4];
		firstPathAtom = paths.get(0).size()>0 ? paths.get(0).get(paths.get(0).size()-1) : a1;
		atomRefs4a2[2] = firstPathAtom;
		secondPathAtom = paths.get(1).size()>0 ? paths.get(1).get(paths.get(1).size()-1) : a1;
		atomRefs4a2[3] = secondPathAtom;
		a2Neighbours.remove(firstPathAtom);
		a2Neighbours.remove(secondPathAtom);
		if (paths.size()==3){
			atomRefs4a2[1] = paths.get(2).size()>0 ? paths.get(2).get(paths.get(2).size()-1) : a1;
		}
		else{
			for (Atom atom : a2Neighbours) {
				if (atom.getElement().equals("H")){
					atomRefs4a2[1] = atom;
					break;
				}
			}
			if (atomRefs4a1[1] ==null){
				throw new StructureBuildingException("OPSIN Bug: cannot assign cis/trans on ring stereochemistry");
			}
		}
		a2Neighbours.remove(atomRefs4a2[1]);
		atomRefs4a2[0] = a2Neighbours.get(0);
		boolean enantiomer =false;
		if (a1.getAtomParity()!=null){
			if (!checkEquivalencyOfAtomsRefs4AndParity(atomRefs4a1, 1, a1.getAtomParity().getAtomRefs4(), a1.getAtomParity().getParity())){
				enantiomer=true;
			}
		}
		else if (a2.getAtomParity()!=null){
			if (cisOrTrans.equals("cis")){
				if (!checkEquivalencyOfAtomsRefs4AndParity(atomRefs4a2, -1, a2.getAtomParity().getAtomRefs4(), a2.getAtomParity().getParity())){
					enantiomer=true;
				}
			}
			else if (cisOrTrans.equals("trans")){
				if (!checkEquivalencyOfAtomsRefs4AndParity(atomRefs4a2, 1, a2.getAtomParity().getAtomRefs4(), a2.getAtomParity().getParity())){
					enantiomer=true;
				}
			}
		}
		if (enantiomer){
			if (cisOrTrans.equals("cis")){
				a1.setAtomParity(atomRefs4a1, -1);
				a2.setAtomParity(atomRefs4a2, 1);
			}
			else if (cisOrTrans.equals("trans")){
				a1.setAtomParity(atomRefs4a1, -1);
				a2.setAtomParity(atomRefs4a2, -1);
			}
		}
		else{
			if (cisOrTrans.equals("cis")){
				a1.setAtomParity(atomRefs4a1, 1);
				a2.setAtomParity(atomRefs4a2, -1);
			}
			else if (cisOrTrans.equals("trans")){
				a1.setAtomParity(atomRefs4a1, 1);
				a2.setAtomParity(atomRefs4a2, 1);
			}
		}
	}
	
	/**
	 * Handles assignment of alpha and beta stereochemistry to appropriate ring systems
	 * Currently these are only assignable to natural products
	 * @param state
	 * @param stereoChemistryEl
	 * @param atomStereoCentreMap
	 * @throws StructureBuildingException
	 */
	private static void assignAlphaBetaStereochem(BuildState state,Element stereoChemistryEl, Map<Atom, StereoCentre> atomStereoCentreMap) throws StructureBuildingException {
		Element parentSubBracketOrRoot = (Element) stereoChemistryEl.getParent();
		List<Fragment> possibleFragments = StructureBuildingMethods.findAlternativeFragments(state, parentSubBracketOrRoot);
		List<Element> adjacentGroupEls = XOMTools.getDescendantElementsWithTagName(parentSubBracketOrRoot, GROUP_EL);
		for (int i = adjacentGroupEls.size()-1; i >=0; i--) {
			possibleFragments.add(state.xmlFragmentMap.get(adjacentGroupEls.get(i)));
		}
		String locant = stereoChemistryEl.getAttributeValue(LOCANT_ATR);
		String alphaOrBets = stereoChemistryEl.getAttributeValue(VALUE_ATR);
		for (Fragment fragment : possibleFragments) {
			Atom potentialStereoAtom = fragment.getAtomByLocant(locant);
			if (potentialStereoAtom !=null && atomStereoCentreMap.containsKey(potentialStereoAtom)){
				String alphaBetaClockWiseAtomOrdering = state.xmlFragmentMap.getElement(fragment).getAttributeValue(ALPHABETACLOCKWISEATOMORDERING_ATR);
				if (alphaBetaClockWiseAtomOrdering==null){
					throw new StructureBuildingException("Identified fragment is not known to be able to support alpha/beta stereochemistry");
				}
				applyAlphaBetaStereochemistryToStereoCentre(potentialStereoAtom, fragment, alphaBetaClockWiseAtomOrdering, alphaOrBets);
				atomStereoCentreMap.remove(potentialStereoAtom);
				return;
			}
		}
		throw new StructureBuildingException("Could not find atom that: " + stereoChemistryEl.toXML() + " appeared to be referring to");
	}


	/**
	 * Converts the alpha/beta descriptor into an atomRefs4 and parity.
	 * The ordering of atoms in the atomsRefs4 is determined by using the two adjacent atoms along the rings edge as defined by ALPHABETACLOCKWISEATOMORDERING_ATR.
	 * by what atom is also part of the ring or is a hydrogen
	 * and by the substituent atom
	 * @param stereoAtom
	 * @param fragment
	 * @param alphaBetaClockWiseAtomOrdering
	 * @param alphaOrBeta
	 * @throws StructureBuildingException
	 */
	private static void applyAlphaBetaStereochemistryToStereoCentre(Atom stereoAtom, Fragment fragment, String alphaBetaClockWiseAtomOrdering, String alphaOrBeta) throws StructureBuildingException {
		List<String> ringOrder = StringTools.arrayToList(matchSlash.split(alphaBetaClockWiseAtomOrdering));
		int positionInList = ringOrder.indexOf(stereoAtom.getFirstLocant());
		if (stereoAtom.getAtomIsInACycle() && positionInList!=-1){
			Atom[] atomRefs4 = new Atom[4];
			List<Atom> neighbours = stereoAtom.getAtomNeighbours();
			if (neighbours.size()==4){
				int previousIndice = positionInList==0 ? ringOrder.size()-1: positionInList -1;
				int nextindice = positionInList==ringOrder.size()-1? 0: positionInList +1;
				atomRefs4[0] = fragment.getAtomByLocantOrThrow(ringOrder.get(previousIndice));
				atomRefs4[3] = fragment.getAtomByLocantOrThrow(ringOrder.get(nextindice));
				neighbours.remove(atomRefs4[0]);
				neighbours.remove(atomRefs4[3]);
				Atom a1 =neighbours.get(0);
				Atom a2 =neighbours.get(1);
				if ((fragment.getAtomList().contains(a1)&& ringOrder.contains(a1.getFirstLocant()))){
					atomRefs4[1]=a1;
					atomRefs4[2]=a2;
				}
				else if ((fragment.getAtomList().contains(a2)&& ringOrder.contains(a2.getFirstLocant()))){
					atomRefs4[1]=a2;
					atomRefs4[2]=a1;
				}
				else if (a1.getElement().equals("H") && !a2.getElement().equals("H")){
					atomRefs4[1]=a2;
					atomRefs4[2]=a1;
				}
				else if (a2.getElement().equals("H") && !a1.getElement().equals("H")){
					atomRefs4[1]=a1;
					atomRefs4[2]=a2;
				}
				else{
					throw new StructureBuildingException("alpha/beta stereochemistry application is ambiguous at position " +stereoAtom.getFirstLocant());
				}
				if (alphaOrBeta.equals("alpha")){
					stereoAtom.setAtomParity(atomRefs4, 1);
				}
				else if (alphaOrBeta.equals("beta")){
					stereoAtom.setAtomParity(atomRefs4, -1);
				}
				else if (alphaOrBeta.equals("xi")){
					stereoAtom.setAtomParity(null);
				}
				else{
					throw new StructureBuildingException("OPSIN Bug: malformed alpha/beta stereochemistry value");
				}
			}
			else{
				throw new StructureBuildingException("Unsupported stereocentre type for alpha/beta stereochemistry");
			}
		}
		else{
			throw new StructureBuildingException("Unsupported stereocentre type for alpha/beta stereochemistry");
		}
	}


	static int swapsRequiredToSort(Atom[] atomRefs4){
	  Atom[] atomRefs4Copy = atomRefs4.clone();
	  int swapsPerformed = 0;
	  int i,j;

	  for (i=atomRefs4Copy.length; --i >=0;) {
		boolean swapped = false;
		for (j=0; j<i;j++) {
			if (atomRefs4Copy[j].getID() > atomRefs4Copy[j+1].getID()){
				Atom temp = atomRefs4Copy[j+1];
				atomRefs4Copy[j+1] = atomRefs4Copy[j];
				atomRefs4Copy[j] = temp;
				swapsPerformed++;
				swapped=true;
			}
		}
		if (!swapped){
			return swapsPerformed;
		}
	  }
	  return swapsPerformed;
	}
	
	static boolean checkEquivalencyOfAtomsRefs4AndParity(Atom[] atomRefs1, int atomParity1, Atom[] atomRefs2, int atomParity2){
		int swaps1 =swapsRequiredToSort(atomRefs1);
		int swaps2 =swapsRequiredToSort(atomRefs2);
		if (atomParity1<0 && atomParity2>0 || atomParity1>0 && atomParity2<0){
			 swaps1++;
		}
		return swaps1 %2 == swaps2 %2;
	}
}

