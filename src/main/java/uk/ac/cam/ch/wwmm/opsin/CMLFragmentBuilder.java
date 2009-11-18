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

	/**Builds a fragment based on the fragment file.
	 *
	 * @param id The name of the fragment to be built.
	 * @param type The type of the fragment, for the purpose of suffix interpretation.
	 * @param subType The subType of the fragment, for the purpose of suffix interpretation.
	 * @return The built fragment.
	 * @throws StructureBuildingException If the fragment cannot be found, or multiple instances are detected,
	 */
	Fragment build(String id, String type, String subType, IDManager idManager) throws StructureBuildingException {
		List<Element> toBuildElements = XOMTools.getChildElementsWithTagNameAndAttribute(fragmentDoc.getRootElement(), "molecule", "id", id);
		if(toBuildElements.size() == 0) throw new StructureBuildingException("Fragment " + id + " could not be found in the fragments file.");
		if(toBuildElements.size() > 1) throw new StructureBuildingException("Fragment " + id + " is not unique in the fragments file.");
		Element cmlMol = toBuildElements.get(0);
		return build(cmlMol, type, subType, idManager);
	}

	/**Builds a fragment based on a CML molecule element.
	 *
	 * @param cmlMol The CML molecule to build.
	 * @param type The type of the fragment, for the purpose of suffix interpretation.
	 * @param subType The subType of the fragment, for the purpose of suffix interpretation.
	 * @return The built fragment.
	 * @throws StructureBuildingException
	 */
	Fragment build(Element cmlMol, String type, String subType, IDManager idManager) throws StructureBuildingException {
		Element atomArray = cmlMol.getFirstChildElement("atomArray");
		Element bondArray = cmlMol.getFirstChildElement("bondArray");
		Map<String, Integer> idDict = buildIDDict(atomArray, idManager);
		Fragment currentFrag = new Fragment(type, subType);
		Elements atoms = atomArray.getChildElements("atom");
		for(int i=0;i<atoms.size();i++) {
			Atom atom = new Atom(idDict.get(atoms.get(i).getAttributeValue("id")) , atoms.get(i), currentFrag);
			currentFrag.addAtom(atom);
		}
		Elements bonds = bondArray.getChildElements("bond");
		for(int i=0;i<bonds.size();i++) {
			String ar2 = bonds.get(i).getAttributeValue("atomRefs2");
			String[] ar2a = ar2.split(" ");
			Atom fromAtom = currentFrag.getAtomByIDOrThrow(idDict.get(ar2a[0]));
			Atom toAtom = currentFrag.getAtomByIDOrThrow(idDict.get(ar2a[1]));
			String bondOrder = bonds.get(i).getAttributeValue("order");
			if(bondOrder.equalsIgnoreCase("s")) bondOrder = "1";
			if(bondOrder.equalsIgnoreCase("d")) bondOrder = "2";
			if(bondOrder.equalsIgnoreCase("t")) bondOrder = "3";
			currentFrag.addBond(new Bond(fromAtom, toAtom, Integer.parseInt(bondOrder)));
		}
		List<Atom> AtomList =currentFrag.getAtomList();
		for (Atom atom : AtomList) {
			atom.calculateValencyFromConnectivity();
		}
		return currentFrag;
	}

	/**Builds idDict: A mapping between CML and OPSIN ids for atoms.
	 *
	 * @param atomArray The CML atomArray that contains the CML ids.
	 */
	private Map<String, Integer> buildIDDict(Element atomArray, IDManager idManager) {
		Map<String, Integer> idDict = new HashMap<String, Integer>();
		Elements atoms = atomArray.getChildElements("atom");
		for(int i=0;i<atoms.size();i++) {
			idDict.put(atoms.get(i).getAttributeValue("id"), idManager.getNextID());
		}
		return idDict;
	}
}
