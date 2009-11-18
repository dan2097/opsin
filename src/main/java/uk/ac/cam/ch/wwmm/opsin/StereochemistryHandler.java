package uk.ac.cam.ch.wwmm.opsin;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import sea36.chem.core.CMLAtom;
import sea36.chem.stereo.StereoAnalyser;
import sea36.chem.stereo.StereoAnalyser.StereoAnalysis;
import sea36.chem.stereo.StereoAnalyser.StereoBond;
import sea36.chem.stereo.StereoAnalyser.StereoCentre;

import nu.xom.Element;

/**
 * Uses ChemKit to identify stereocentres, assigns stereochemistry elements to them and then uses the CIP rules to calculate appropriates atomParity/bondstereo tags
 * @author dl387
 *
 */
public class StereochemistryHandler {

	/**
	 * Master method for assigning and processing stereochemistry elements
	 * @param state
	 * @param molecule 
	 * @param uniFrag
	 * @param stereoChemistryEls
	 * @throws StructureBuildingException 
	 */
	static void processStereochemicalElements(BuildState state, Element molecule, Fragment uniFrag, List<Element> stereoChemistryEls) throws StructureBuildingException {
		OpsinToChemKitWrapper chemKitMoleculeWrapper = new OpsinToChemKitWrapper(uniFrag);
		StereoAnalysis stereoAnalysis = new StereoAnalyser().findStereoCentres(chemKitMoleculeWrapper.getChemKitMolecule());
	    Map<Atom, StereoCentre> atomStereoCentreMap = new HashMap<Atom, StereoCentre>();//contains all atoms that are stereo centres with a mapping to the corresponding StereoCentre object
		List<StereoCentre> stereoCentres = stereoAnalysis.getStereoCentres();
		for (StereoCentre stereoCentre : stereoCentres) {
			atomStereoCentreMap.put(chemKitMoleculeWrapper.getOpsinAtomFromChemKitAtom(stereoCentre.getStereoAtom()),stereoCentre);
		}
	    Map<Bond, StereoBond> bondStereoBondMap = new HashMap<Bond, StereoBond>();
		List<StereoBond> stereoBonds = stereoAnalysis.getStereoBonds();
		for (StereoBond stereoBond : stereoBonds) {
			Bond b = uniFrag.findBondOrThrow(chemKitMoleculeWrapper.getOpsinAtomFromChemKitAtom(stereoBond.getBond().getAtom0()),chemKitMoleculeWrapper.getOpsinAtomFromChemKitAtom( stereoBond.getBond().getAtom1()));
			bondStereoBondMap.put(b, stereoBond);
		}
		
		matchStereochemistryToAtomsAndBonds(state, molecule, stereoChemistryEls, atomStereoCentreMap, bondStereoBondMap);
	}
	
	
	/**
	 * Finds all stereochemistry elements and attempts to locate a suitable atom/bond for them
	 * Currently only locanted stereochemistry elements are supported
	 * @param state 
	 * @param stereoChemistryEls A list of stereochemistry elements
	 * @param molecule The molecule element
	 * @param bondStereoBondMap 
	 * @param atomStereoCentreMap 
	 * @throws StructureBuildingException 
	 */
	private static void matchStereochemistryToAtomsAndBonds(BuildState state, Element molecule, List<Element> stereoChemistryEls, Map<Atom, StereoCentre> atomStereoCentreMap, Map<Bond, StereoBond> bondStereoBondMap) throws StructureBuildingException {
		for (int i = stereoChemistryEls.size()-1; i >=0 ; i--) {
			Element stereoChemistryEl = stereoChemistryEls.get(i);
			String stereoChemistryType =stereoChemistryEl.getAttributeValue("type");
			if (stereoChemistryType.equals("RorS")){
				assignStereoCentre(state, stereoChemistryEl, atomStereoCentreMap);
			}
			else if (stereoChemistryType.equals("EorZ")){
				assignStereoBond(state, stereoChemistryEl, bondStereoBondMap);
			}
			else if (stereoChemistryType.equals("cisOrTrans")){
				assignStereoBond(state, stereoChemistryEl, bondStereoBondMap);
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
		Element parent = (Element) stereoChemistryEl.getParent().getParent();//want to iterate at the level above the containing substituent or bracket
		//generally the LAST group in this list will be the appropriate groups e.g. (5S)-5-ethyl-6-methylheptane where the heptane is the appropriate group
		List<Element> possibleGroups = XOMTools.getDescendantElementsWithTagName(parent, "group");
		String locant = StructureBuildingMethods.getLocant(stereoChemistryEl);
		String rOrS = stereoChemistryEl.getAttributeValue("value");
		for (int i = possibleGroups.size()-1; i >=0; i--) {//groups further right in scope preferable
			Fragment correspondingFrag = state.xmlFragmentMap.get(possibleGroups.get(i));
			if (locant.equals("0")){//undefined locant
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
		if (parent.getLocalName().equals("word") && parent.getAttributeValue("type").equals("substituent")){
			//something like (3R,4R,5R)-ethyl 4-acetamido-5-amino-3-(pentan-3-yloxy)cyclohex-1-enecarboxylate
			//I think this is a violation of the IUPAC rules...but anyway...
			List<Element> words = XOMTools.getChildElementsWithTagNameAndAttribute(((Element)parent.getParent()), "word", "type", "full");
			for (Element word : words) {
				possibleGroups = XOMTools.getDescendantElementsWithTagName(word, "group");
				for (int i = possibleGroups.size()-1; i >=0; i--) {
					Fragment correspondingFrag = state.xmlFragmentMap.get(possibleGroups.get(i));
					if (locant.equals("0")){//undefined locant
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
		List<CMLAtom> cipOrderedAtoms =stereoCentre.getNeighbours();
		if (cipOrderedAtoms.size()!=4){
			throw new StructureBuildingException("Only tetrahedral chirality is currently supported");
		}
		String atomRefs4= cipOrderedAtoms.get(cipOrderedAtoms.size()-1).getId();//this is "a" + opsin's atom id
		for (int i = 0; i < cipOrderedAtoms.size() -1; i++) {//from highest to lowest (true for S) hence atomParity 1 for S
			atomRefs4 +=" "+cipOrderedAtoms.get(i).getId();
		}
		if (rOrS.equals("R")){
			atom.setAtomParityElement(atomRefs4, -1);
		}
		else if (rOrS.equals("S")){
			atom.setAtomParityElement(atomRefs4, 1);
		}
		else{
			throw new StructureBuildingException("Unexpected stereochemistry type: " + rOrS);
		}
	}


	/**
	 * Handles E/Z stereochemistry
	 * @param state
	 * @param stereoChemistryEl
	 * @param bondStereoBondMap
	 * @throws StructureBuildingException
	 */
	private static void assignStereoBond(BuildState state, Element stereoChemistryEl, Map<Bond, StereoBond> bondStereoBondMap) throws StructureBuildingException {
		Element parent = (Element) stereoChemistryEl.getParent().getParent();//want to iterate at the level above the containing substituent or bracket
		//generally the LAST group in this list will be the appropriate groups e.g. (2Z)-5-ethyl-6-methylhex-2-ene where the hex-2-ene is the appropriate group
		List<Element> possibleGroups = XOMTools.getDescendantElementsWithTagName(parent, "group");
		String locant = StructureBuildingMethods.getLocant(stereoChemistryEl);
		String eOrZ = stereoChemistryEl.getAttributeValue("value");
		for (int i = possibleGroups.size()-1; i >=0; i--) {
			Fragment correspondingFrag = state.xmlFragmentMap.get(possibleGroups.get(i));
			if (locant.equals("0")){//undefined locant
				List<Bond> bondList = correspondingFrag.getBondList();
				for (Bond potentialBond : bondList) {
					if (bondStereoBondMap.containsKey(potentialBond)){
						applyStereoChemistryToStereoBond(potentialBond, bondStereoBondMap.get(potentialBond), eOrZ);
						bondStereoBondMap.remove(potentialBond);
						return;
					}
				}
				Set<Bond> interFragmentBonds = state.fragManager.getInterFragmentBonds(correspondingFrag);
				for (Bond potentialBond : interFragmentBonds) {
					if (bondStereoBondMap.containsKey(potentialBond)){
						applyStereoChemistryToStereoBond(potentialBond, bondStereoBondMap.get(potentialBond), eOrZ);
						bondStereoBondMap.remove(potentialBond);
						return;
					}
				}
			}
			else{
				Atom firstAtomInBond = correspondingFrag.getAtomByLocant(locant);
				if (firstAtomInBond !=null){
					List<Bond> bonds = firstAtomInBond.getBonds();
					for (Bond potentialBond : bonds) {
						if (bondStereoBondMap.containsKey(potentialBond)){
							applyStereoChemistryToStereoBond(potentialBond, bondStereoBondMap.get(potentialBond), eOrZ);
							bondStereoBondMap.remove(potentialBond);
							return;
						}
					}
				}
			}
		}
		if (parent.getLocalName().equals("word") && parent.getAttributeValue("type").equals("substituent")){
			//the element is in front of a substituent and may refer to the full group
			List<Element> words = XOMTools.getChildElementsWithTagNameAndAttribute(((Element)parent.getParent()), "word", "type", "full");
			for (Element word : words) {
				possibleGroups = XOMTools.getDescendantElementsWithTagName(word, "group");
				for (int i = possibleGroups.size()-1; i >=0; i--) {
					Fragment correspondingFrag = state.xmlFragmentMap.get(possibleGroups.get(i));
					if (locant.equals("0")){//undefined locant
						List<Bond> bondList = correspondingFrag.getBondList();
						for (Bond potentialBond : bondList) {
							if (bondStereoBondMap.containsKey(potentialBond)){
								applyStereoChemistryToStereoBond(potentialBond, bondStereoBondMap.get(potentialBond), eOrZ);
								bondStereoBondMap.remove(potentialBond);
								return;
							}
						}
						Set<Bond> interFragmentBonds = state.fragManager.getInterFragmentBonds(correspondingFrag);
						for (Bond potentialBond : interFragmentBonds) {
							if (bondStereoBondMap.containsKey(potentialBond)){
								applyStereoChemistryToStereoBond(potentialBond, bondStereoBondMap.get(potentialBond), eOrZ);
								bondStereoBondMap.remove(potentialBond);
								return;
							}
						}
					}
					else{
						Atom firstAtomInBond = correspondingFrag.getAtomByLocant(locant);
						if (firstAtomInBond !=null){
							List<Bond> bonds = firstAtomInBond.getBonds();
							for (Bond potentialBond : bonds) {
								if (bondStereoBondMap.containsKey(potentialBond)){
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
		throw new StructureBuildingException("Could not find bond that: " + stereoChemistryEl.toXML() + " was referring to");
	}
	

	/**
	 * Assigns bondstereo to the given bond in accordance with the CIP rules
	 * @param bond The stereobond
	 * @param stereoBond
	 * @param eOrZ The stereo description given in the name
	 * @throws StructureBuildingException
	 */
	private static void applyStereoChemistryToStereoBond(Bond bond, StereoBond stereoBond, String eOrZ ) throws StructureBuildingException {
		List<CMLAtom> stereoBondAtoms = stereoBond.getStereoAtoms();
		//stereoBondAtoms contains the higher priority atom at one end, the two bond atoms and the higher priority atom at the other end
		String atomRefs4= stereoBondAtoms.get(0).getId();
		atomRefs4 +=" " + stereoBondAtoms.get(1).getId();
		atomRefs4 +=" " + stereoBondAtoms.get(2).getId();
		atomRefs4 +=" " + stereoBondAtoms.get(3).getId();
		if (eOrZ.equals("E")){
			bond.setBondStereoElement(atomRefs4, "T");
		}
		else if (eOrZ.equals("Z")){
			bond.setBondStereoElement(atomRefs4, "C");
		}
		else{
			throw new StructureBuildingException("Unexpected stereochemistry type: " + bond.getStereochemistry());
		}
	}

}
