package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;

import nu.xom.Attribute;
import nu.xom.Element;
import nu.xom.Elements;
import nu.xom.Node;


//TODO make this more efficient. You do not need the entire breadthFirstCIPTree to resolve the CIP rules in most cases

/**
 * This function gets the a breadth first tree of atoms from a starting point in accordance with the CIP (Cahn-Ingold-Prelog rules)
 * One can then determine from this the priority of the atoms as given by the CIP system and henc convert from R/S to atom parity.
 * 
 * Differentiation based on isotopic number is not supported
 * @author pm286 (modifed by dl387 to use OPSIN's internal data structure and to be thread safe)
 *
 */
class CIPresolver {
	private final static String ATNUM = "atomicNumber";
	private final static String GHOST = "ghost";
	private final static String ID = "id";
	private final static String NODE = "node";
	private final static String PARENT = "parent";
	private final static String TRUE = "true";
	
	private static HashMap<String, Integer> atomicNumber;

	static {
		atomicNumber = new HashMap<String, Integer>();
		atomicNumber.put("H", 1);
		atomicNumber.put("He", 2);
		atomicNumber.put("Li", 3);
		atomicNumber.put("Be", 4);
		atomicNumber.put("B", 5);
		atomicNumber.put("C", 6);
		atomicNumber.put("N", 7);
		atomicNumber.put("O", 8);
		atomicNumber.put("F", 9);
		atomicNumber.put("Ne", 10);
		atomicNumber.put("Na", 11);
		atomicNumber.put("Mg", 12);
		atomicNumber.put("Al", 13);
		atomicNumber.put("Si", 14);
		atomicNumber.put("P", 15);
		atomicNumber.put("S", 16);
		atomicNumber.put("Cl", 17);
		atomicNumber.put("Ar", 18);
		atomicNumber.put("K", 19);
		atomicNumber.put("Ca", 20);
		atomicNumber.put("Sc", 21);
		atomicNumber.put("Ti", 22);
		atomicNumber.put("V", 23);
		atomicNumber.put("Cr", 24);
		atomicNumber.put("Mn", 25);
		atomicNumber.put("Fe", 26);
		atomicNumber.put("Co", 27);
		atomicNumber.put("Ni", 28);
		atomicNumber.put("Cu", 29);
		atomicNumber.put("Zn", 30);
		atomicNumber.put("Ga", 31);
		atomicNumber.put("Ge", 32);
		atomicNumber.put("As", 33);
		atomicNumber.put("Se", 34);
		atomicNumber.put("Br", 35);
		atomicNumber.put("Kr", 36);
		atomicNumber.put("Rb", 37);
		atomicNumber.put("Sr", 38);
		atomicNumber.put("Y", 39);
		atomicNumber.put("Zr", 40);
		atomicNumber.put("Nb", 41);
		atomicNumber.put("Mo", 42);
		atomicNumber.put("Tc", 43);
		atomicNumber.put("Ru", 44);
		atomicNumber.put("Rh", 45);
		atomicNumber.put("Pd", 46);
		atomicNumber.put("Ag", 47);
		atomicNumber.put("Cd", 48);
		atomicNumber.put("In", 49);
		atomicNumber.put("Sn", 50);
		atomicNumber.put("Sb", 51);
		atomicNumber.put("Te", 52);
		atomicNumber.put("I", 53);
		atomicNumber.put("Xe", 54);
		atomicNumber.put("Cs", 55);
		atomicNumber.put("Ba", 56);
		atomicNumber.put("La", 57);
		atomicNumber.put("Ce", 58);
		atomicNumber.put("Pr", 59);
		atomicNumber.put("Nd", 60);
		atomicNumber.put("Pm", 61);
		atomicNumber.put("Sm", 62);
		atomicNumber.put("Eu", 63);
		atomicNumber.put("Gd", 64);
		atomicNumber.put("Tb", 65);
		atomicNumber.put("Dy", 66);
		atomicNumber.put("Ho", 67);
		atomicNumber.put("Er", 68);
		atomicNumber.put("Tm", 69);
		atomicNumber.put("Yb", 70);
		atomicNumber.put("Lu", 71);
		atomicNumber.put("Hf", 72);
		atomicNumber.put("Ta", 73);
		atomicNumber.put("W", 74);
		atomicNumber.put("Re", 75);
		atomicNumber.put("Os", 76);
		atomicNumber.put("Ir", 77);
		atomicNumber.put("Pt", 78);
		atomicNumber.put("Au", 79);
		atomicNumber.put("Hg", 80);
		atomicNumber.put("Tl", 81);
		atomicNumber.put("Pb", 82);
		atomicNumber.put("Bi", 83);
		atomicNumber.put("Po", 84);
		atomicNumber.put("At", 85);
		atomicNumber.put("Rn", 86);
		atomicNumber.put("Fr", 87);
		atomicNumber.put("Ra", 88);
		atomicNumber.put("Ac", 89);
		atomicNumber.put("Th", 90);
		atomicNumber.put("Pa", 91);
		atomicNumber.put("U", 92);
		atomicNumber.put("Np", 93);
		atomicNumber.put("Pu", 94);
		atomicNumber.put("Am", 95);
		atomicNumber.put("Cm", 96);
		atomicNumber.put("Bk", 97);
		atomicNumber.put("Cf", 98);
		atomicNumber.put("Es", 99);
		atomicNumber.put("Fm", 100);
		atomicNumber.put("Md", 101);
		atomicNumber.put("No", 102);
		atomicNumber.put("Lr", 103);
		atomicNumber.put("Rf", 104);
		atomicNumber.put("Db", 105);
		atomicNumber.put("Sg", 106);
		atomicNumber.put("Bh", 107);
		atomicNumber.put("Hs", 108);
		atomicNumber.put("Mt", 109);
		atomicNumber.put("Ds", 110);
	}

	/** Returns the CIP order of the atoms in startingAtoms. The rootAtom is the stereogenic centre
	 * uses CIP strategy (cf. Wikipedia)
	 * @param molecule
	 * @param rootAtom
	 * @param startingAtoms
	 * @return
	 * @throws StructureBuildingException 
	 */
	public static List<Atom> getCIPOrderedAtoms(Fragment molecule, Atom rootAtom, List<Atom> startingAtoms) throws StructureBuildingException {
		Queue<Element> nodeQueue = new  LinkedList<Element>();
		Element root = makeNode("n/a", Integer.toString(rootAtom.getID()), getAtomicNumber(rootAtom));
		for (Atom atom : startingAtoms) {
			Element newNode =makeNode(Integer.toString(rootAtom.getID()), Integer.toString(atom.getID()), getAtomicNumber(atom));
			nodeQueue.add(newNode);
			root.appendChild(newNode);
		}
		while (!nodeQueue.isEmpty()) {
			Element nextElement = nodeQueue.remove();
			processElement(molecule, nextElement, nodeQueue);
		}
		sortNodesRecursively(root);
		Elements startingEls = root.getChildElements();
		List<Atom> orderedAtoms = new ArrayList<Atom>();
		for (int i = 0; i < startingEls.size(); i++) {
			orderedAtoms.add(molecule.getAtomByID(Integer.parseInt(startingEls.get(i).getAttributeValue(ID))));
		}
		return orderedAtoms;
	}

	/*
	 * explore from atom
	 */
	private static Element processElement(Fragment molecule, Element node, Queue<Element> nodeQueue) throws StructureBuildingException {
		String atomId = node.getAttributeValue(ID);
		String parentId = node.getAttributeValue(PARENT);
		Atom atom = molecule.getAtomByID(Integer.parseInt(atomId));
		List<Atom> ligandList = atom.getAtomNeighbours();
		List<Bond> ligandBondList = atom.getBonds();//ligandList and ligandBondList are in the same order
		for (int i = 0; i < ligandList.size(); i++) {
			Bond ligandBond = ligandBondList.get(i);
			Atom ligandAtom = ligandList.get(i);
			String ligandId = Integer.toString(ligandAtom.getID());
			int ligandAtomicNumber = getAtomicNumber(ligandAtom);
			int order = ligandBond.getOrder();
			// multiple bonds (regardless of whether already in tree)
			for (int j = 2; j <= order; j++) {
				addGhost(node, ligandId, ligandAtomicNumber);
			}
			if (ligandId.equals(parentId)) {
				// no action
			} else if (isAncestorOfLigandId(node, ligandId)) {
				addGhost(node, ligandId, ligandAtomicNumber);
			} else {
				Element newNode = makeNode(atomId, ligandId, ligandAtomicNumber);
				node.appendChild(newNode);
				nodeQueue.add(newNode);
			}
		}
		return node;
	}

	/**
	 * @param ligandAtom
	 * @return
	 */
	private static int getAtomicNumber(Atom ligandAtom) {
		String ligandElementType = ligandAtom.getElement();
		return atomicNumber.get(ligandElementType);
	}
	
	private static boolean isAncestorOfLigandId(Element node, String ligandId) {
		Element currentNode = (Element) node.getParent();
		while (currentNode !=null){
			if (currentNode.getAttributeValue(ID).equals(ligandId)){
				return true;
			}
			currentNode =(Element) currentNode.getParent();
		}
		return false;
	}
	
	private static void addGhost(Element node, String ligandId, int ligandAtomicNumber) {
		Element ghost = makeNode(node.getAttributeValue(ID), ligandId+ "_" +GHOST, ligandAtomicNumber);
		ghost.addAttribute(new Attribute(GHOST, TRUE));
		node.appendChild(ghost);
	}
	
	private static Element makeNode(String parentId, String atomId, int atomAtomicNumber) {
		Element node = new Element(NODE);
		node.addAttribute(new Attribute(PARENT, parentId));
		node.addAttribute(new Attribute(ID, atomId));
		node.addAttribute(new Attribute(ATNUM, ""+atomAtomicNumber));
		return node;
	}

	/** sort nodes on each atom.
	 * ordering is by atomic number, then ghosts
	 * should put nodes with more children first, but I am not sure of this
	 * @param node
	 */
	private static void sortNodesRecursively(Element node) {
		int nNodes = node.getChildCount();
		boolean change = true;
		// sort children first
		for (int i = 0; i < nNodes; i++) {
			sortNodesRecursively((Element) node.getChild(i));
		}
		while (change) {
			change = false;
			for (int i = 0; i < nNodes - 1; i++) {
				Element nodei = (Element) node.getChild(i);
				for (int j = i+1 ; j < nNodes; j++) {
					Element nodej = (Element) node.getChild(j);
					int compare = compareWithGhosts(nodei, nodej);
					if (compare == 0) {
						compare = compareChildrenRecursively(nodei, nodej);
					}
					if (compare < 0) {
						swap(node, i, j);
						change = true;
						break;
					}
				}
			}
		}
	}

	private static int compareChildrenRecursively(Element nodei, Element nodej) {
		int compare = 0;
		// seed comparison
		List<Element> nodesi = new ArrayList<Element>();
		nodesi.add(nodei);
		List<Element> nodesj = new ArrayList<Element>();
		nodesj.add(nodej);
		while (true) {
			// compare at current level of tree
			int nij = Math.min(nodesi.size(), nodesj.size());
			// compare common children
			for (int ij = 0; ij < nij; ij++) {
				compare = compareWithoutGhosts(nodesi.get(ij), nodesj.get(ij));
				if (compare != 0) {
					break;
				}
			}
			// if common children match, are there any children unique to one tree?
			compare = (compare == 0) ? nodesi.size() - nodesj.size() : compare;
			// found inequality
			if (compare != 0) {
				break;
			}
			// recurse to next level of tree
			nodesi = getChildNodes(nodesi);
			nodesj = getChildNodes(nodesj);
			// no more levels in either tree?
			if (nodesi.size() == 0 && nodesj.size() == 0) {
				break;
			}
		}
		return compare;
	}

	private static List<Element> getChildNodes(List<Element> nodes) {
		List<Element> childNodeList = new ArrayList<Element>();
		for (Element node : nodes) {
			Elements childNodes = node.getChildElements("node");
			for (int i = 0; i < childNodes.size(); i++) {
				childNodeList.add((Element) childNodes.get(i));
			}
		}
		return childNodeList;
	}

	/**
	 * compares nodes on atomic numbers
	 * if a tie then elements with ghosts have lower priority
	 * @param nodei
	 * @param nodej
	 * @return atnumi - atnumj
	 * @throws NumberFormatException
	 */
	private static int compareWithGhosts(Element nodei, Element nodej) {
		int compare = compareWithoutGhosts(nodei, nodej);
		if (compare == 0) {
			Attribute ghosti = nodei.getAttribute(GHOST);
			Attribute ghostj = nodej.getAttribute(GHOST);
			// ghosts have lower priority
			if (ghosti == null && ghostj != null) {
				compare = 1;
			}
			if (ghosti != null && ghostj == null) {
				compare = -1;
			}
		}
		return compare;
	}

	/**
	 * @param nodei
	 * @param nodej
	 * @return
	 * @throws NumberFormatException
	 */
	private static int compareWithoutGhosts(Element nodei, Element nodej)
			throws NumberFormatException {
		// sort on atomic numbers
		int compare = 
			Integer.parseInt(nodei.getAttributeValue(ATNUM)) -
			Integer.parseInt(nodej.getAttributeValue(ATNUM));
		return compare;
	}

	// j > i
	private static void swap(Element node, int i, int j) {
		Node nodei = node.getChild(i);
		Node nodej = node.getChild(j);
		nodej.detach();
		node.insertChild(nodej, i);
		nodei.detach();
		node.insertChild(nodei, j);
	}

}
