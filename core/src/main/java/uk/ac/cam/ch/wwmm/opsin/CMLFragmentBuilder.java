package uk.ac.cam.ch.wwmm.opsin;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import nu.xom.Document;
import nu.xom.Element;
import nu.xom.Elements;

/**A class for building Fragments, based upon a Fragment file in CML.
 *
 * @author ptc24
 *
 */
class CMLFragmentBuilder {

	/**A nu.xom.Document containing the fragments*/
	private Document fragmentDoc;

	/**Initialises the builder.
	 * @param resourceGetter
	 *
	 * @throws Exception If the fragments file can't be found, or the document can't be built.
	 */
	CMLFragmentBuilder(ResourceGetter resourceGetter) throws Exception {
		fragmentDoc = resourceGetter.getXMLDocument("fragments.xml");
	}

	/**
	 * Builds a fragment based on CML in the fragment file.
	 * @param id The name of the fragment to be built.
	 * @param type The type of the fragment, for the purpose of suffix interpretation
	 * @param subType The subType of the fragment, for the purpose of suffix interpretation.
	 * @param fragManager
	 * @return The built fragment.
	 * @throws StructureBuildingException If the fragment cannot be found, or multiple instances are detected,
	 */
	Fragment build(String id, String type, String subType, FragmentManager fragManager) throws StructureBuildingException {
		List<Element> toBuildElements = XOMTools.getChildElementsWithTagNameAndAttribute(fragmentDoc.getRootElement(), "molecule", "id", id);
		if(toBuildElements.size() == 0) throw new StructureBuildingException("Fragment " + id + " could not be found in the fragments file.");
		if(toBuildElements.size() > 1) throw new StructureBuildingException("Fragment " + id + " is not unique in the fragments file.");
		Element cmlMol = toBuildElements.get(0);
		return build(cmlMol, type, subType, fragManager);
	}

	/**Builds a fragment based on a CML molecule element.
	 *
	 * @param cmlMol The CML molecule to build.
	 * @param type The type of the fragment, for the purpose of suffix interpretation.
	 * @param subType The subType of the fragment, for the purpose of suffix interpretation.
	 * @return The built fragment.
	 * @throws StructureBuildingException
	 */
	private Fragment build(Element cmlMol, String type, String subType, FragmentManager fragManager) throws StructureBuildingException {
		Element atomArray = cmlMol.getFirstChildElement("atomArray");
		Element bondArray = cmlMol.getFirstChildElement("bondArray");
		Map<String, Atom> idDict = new HashMap<String, Atom>();// a mapping between CML ids and OPSIN atoms
		Fragment currentFrag = new Fragment(type, subType);
		Elements atoms = atomArray.getChildElements("atom");
		for(int i=0;i<atoms.size();i++) {
			Atom atom = buildAtomFromCML(fragManager, atoms.get(i), currentFrag);
			idDict.put(atoms.get(i).getAttributeValue("id"), atom);
		}
		Elements bonds = bondArray.getChildElements("bond");
		for(int i=0;i<bonds.size();i++) {
			String ar2 = bonds.get(i).getAttributeValue("atomRefs2");
			String[] ar2a = ar2.split(" ");
			Atom fromAtom = idDict.get(ar2a[0]);
			Atom toAtom = idDict.get(ar2a[1]);
			String bondOrder = bonds.get(i).getAttributeValue("order");
			if(bondOrder.equalsIgnoreCase("s")) bondOrder = "1";
			if(bondOrder.equalsIgnoreCase("d")) bondOrder = "2";
			if(bondOrder.equalsIgnoreCase("t")) bondOrder = "3";
			fragManager.createBond(fromAtom, toAtom, Integer.parseInt(bondOrder));
		}
		CycleDetector.assignWhetherAtomsAreInCycles(currentFrag);
		return currentFrag;
	}
	

	/**Builds an Atom from a CML Atom tag. Looks at elementType and formalCharge
	 * attributes, and label tags contained within. id attributes are ignored.
	 *
	 * @param fragManager The current fragment manager
	 * @param cmlAtom The nu.xom.Element for the Atom in CML
	 * @param frag the Fragment to contain the Atom
	 * @throws StructureBuildingException 
	 */
	Atom buildAtomFromCML(FragmentManager fragManager, Element cmlAtom, Fragment frag) throws StructureBuildingException {
		String elementType = cmlAtom.getAttributeValue("elementType");
		if (elementType.equals("H")){
			throw new StructureBuildingException("Explicit hydrogens are not yet supported in OPSIN's CML reading implementation");
		}
		Atom a  =fragManager.createAtom(elementType, frag);
		Elements cmlLocants = cmlAtom.getChildElements("label");
		for(int i=0;i<cmlLocants.size();i++){
			a.addLocant(cmlLocants.get(0).getAttributeValue("value"));
		}
		String chargeStr = cmlAtom.getAttributeValue("formalCharge");
		if(chargeStr != null){
			a.setCharge(Integer.parseInt(chargeStr));
		}
		String hcStr = cmlAtom.getAttributeValue("hydrogenCount");
		if(hcStr != null) {
			throw new StructureBuildingException("Hydrogen count attribute is not yet supported in OPSIN's CML reading implementation");
		}
		return a;
	}
}
