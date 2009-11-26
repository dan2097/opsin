package uk.ac.cam.ch.wwmm.opsin;

import org.junit.Before;
import org.junit.Test;

import static junit.framework.Assert.*;

public class FragmentTest {

	private Fragment frag;

	@Before
	public void setUp() throws Exception {
		frag = new Fragment();
	}

	@Test
	public void testFragment() {
		assertNotNull("Has atom list", frag.getAtomList());
	}

	//FIXME Argh! I hate namespaces!
	/*public void testtoCMLMolecule() {
		Element elem = frag.toCMLMolecule();
		assertNotNull("Got an Element", elem);
		assertEquals("Element is a cml tag", elem.getLocalName(), "cml");
		frag.addAtom(new Atom(1, 1, "C", frag));
		elem = frag.toCMLMolecule();
		assertEquals("foo", elem.getFirstChildElement("molecule").toXML(), "foo bar");
		assertEquals("Element has 1 atom greatgrandchild", 1, elem.getFirstChildElement("molecule")
				.getFirstChildElement("atomArray").getChildElements("atom").size());
		frag.addAtom(new Atom(2, 2, "C", frag));
		elem = frag.toCMLMolecule();
		assertEquals("Element has 2 atom greatgrandchildren", 2, elem.getFirstChildElement("molecule")
				.getFirstChildElement("atomArray").getChildElements("atom").size());
		frag.addBond(new Bond(1, 2, 1));
		elem = frag.toCMLMolecule();
		assertEquals("Element has 1 bond greatgrandchildren", 1, elem.getFirstChildElement("molecule")
				.getFirstChildElement("bondArray").getChildElements("bond").size());
	}*/

	@Test
	public void testAddAtom() throws StructureBuildingException {
		assertEquals("Has no atoms", 0, frag.getAtomList().size());
		frag.addAtom(new Atom(1, "1", "C", frag));
		assertEquals("Now has one atom", 1, frag.getAtomList().size());
	}

	@Test
	public void testAddBond() throws StructureBuildingException {
		frag.addAtom(new Atom(1, "1", "C", frag));
		frag.addAtom(new Atom(2, "2", "C", frag));
		assertEquals("Has no bonds", 0, frag.getBondList().size());
		frag.addBond(new Bond(frag.getAtomByID(1), frag.getAtomByID(2), 1));
		assertEquals("Now has one bond", 1, frag.getBondList().size());
	}

	@Test
	public void testImportFrag() throws StructureBuildingException {
		frag.addAtom(new Atom(1, "1", "C", frag));
		frag.addAtom(new Atom(2, "2", "C", frag));
		frag.addBond(new Bond(frag.getAtomByID(1), frag.getAtomByID(2), 1));
		Fragment newFrag = new Fragment();
		newFrag.addAtom(new Atom(3, "1", "C", frag));
		newFrag.addAtom(new Atom(4, "2", "C", frag));
		newFrag.addBond(new Bond(newFrag.getAtomByID(3), newFrag.getAtomByID(4), 1));
		assertEquals("Fragment has two atoms", 2, frag.getAtomList().size());
		assertEquals("Fragment has one bond", 1, frag.getBondList().size());
		frag.importFrag(newFrag);
		assertEquals("Fragment now has two atoms", 4, frag.getAtomList().size());
		assertEquals("Fragment now has one bond", 2, frag.getBondList().size());
	}

	@Test
	public void testGetIDFromLocant() throws StructureBuildingException {
		Atom atom = new Atom(10, "1", "C", frag);
		atom.addLocant("a");
		frag.addAtom(atom);
		atom = new Atom(20, "2", "C", frag);
		atom.addLocant("silly");
		frag.addAtom(atom);
		assertEquals("Locant a has ID 10", 10, frag.getIDFromLocant("a"));
		assertEquals("Locant silly has ID 20", 20, frag.getIDFromLocant("silly"));
		assertEquals("Locant 42 is not present", 0, frag.getIDFromLocant("42"));
	}

	@Test
	public void testGetAtomByLocant() throws StructureBuildingException {
		Atom atom1 = new Atom(10, "1", "C", frag);
		atom1.addLocant("a");
		frag.addAtom(atom1);
		Atom atom2 = new Atom(20, "2", "C", frag);
		atom2.addLocant("silly");
		frag.addAtom(atom2);
		assertEquals("Locant a gets atom1", atom1, frag.getAtomByLocant("a"));
		assertEquals("Locant silly gets atom2", atom2, frag.getAtomByLocant("silly"));
		assertNull("Locant 42 is not present", frag.getAtomByLocant("42"));
	}

	@Test
	public void testGetAtomByID() throws StructureBuildingException {
		Atom atom1 = new Atom(10, "1", "C", frag);
		frag.addAtom(atom1);
		Atom atom2 = new Atom(20, "2", "C", frag);
		frag.addAtom(atom2);
		assertEquals("ID 10 gets atom1", atom1, frag.getAtomByID(10));
		assertEquals("ID 20 gets atom2", atom2, frag.getAtomByID(20));
		assertNull("ID 42 is not present", frag.getAtomByID(42));
	}

	@Test
	public void testReduceSpareValency() throws Exception {
		Atom atom = new Atom(10, "1", "C", frag);
		frag.addAtom(atom);
		atom.setSpareValency(1);
		assertEquals("Atom a has spareValency 1", 1, atom.getSpareValency());
		atom.subtractSpareValency(1);
		assertEquals("Atom a has spareValency 0", 0, atom.getSpareValency());
		try {
			atom.subtractSpareValency(1);
			assertFalse("Should throw an exception", true);
		} catch(StructureBuildingException e) {

		}
	}

	@Test
	public void testFindBond() throws StructureBuildingException {
		frag.addAtom(new Atom(1, "1", "C", frag));
		frag.addAtom(new Atom(2, "2", "C", frag));
		frag.addAtom(new Atom(3, "3", "N", frag));
		frag.addAtom(new Atom(4, "4", "O", frag));
		frag.addBond(new Bond(frag.getAtomByID(2), frag.getAtomByID(4), 2));
		frag.addBond(new Bond(frag.getAtomByID(1), frag.getAtomByID(2), 1));
		frag.addBond(new Bond(frag.getAtomByID(1), frag.getAtomByID(3), 3));
		Bond b = frag.findBond(2, 4);
		assertNotNull("Found a bond", b);
		assertEquals("..a double bond", 2, b.getOrder());
		b = frag.findBond(3, 1);
		assertNotNull("Found a different bond", b);
		assertEquals("..a triple bond", 3, b.getOrder());
		b = frag.findBond(2, 3);
		assertNull("Don't find non-existent bonds", b);
	}

	@Test
	public void testGetChainLength() throws StructureBuildingException {
		assertEquals("No chain", 0, frag.getChainLength());
		frag.addAtom(new Atom(1, "1", "C", frag));
		assertEquals("Methane", 1, frag.getChainLength());
		frag.addAtom(new Atom(2, "2", "C", frag));
		frag.addBond(new Bond(frag.getAtomByID(1), frag.getAtomByID(2), 1));
		assertEquals("ethane", 2, frag.getChainLength());
		frag.addAtom(new Atom(3, "3", "C", frag));
		frag.addBond(new Bond(frag.getAtomByID(2), frag.getAtomByID(3), 1));
		assertEquals("propane", 3, frag.getChainLength());

	}


//	@Test
//	public void testRelabelSuffixLocants() throws StructureBuildingException {
//		SMILESFragmentBuilder sBuilder = new SMILESFragmentBuilder();
//		IDManager idManager = new IDManager();
//		Fragment parentFrag = sBuilder.build("CC", idManager);
//		frag = sBuilder.build("N", idManager);
//		assertEquals("Can't find locant N in frag", 0, frag.getIDFromLocant("N"));
//		frag.relabelSuffixLocants(parentFrag);
//		assertEquals("Can find locant N in frag: ID = 3", 3, frag.getIDFromLocant("N"));
//		parentFrag.importFrag(frag);
//		frag = sBuilder.build("N", idManager);
//		assertEquals("Can't find locant N' in frag", 0, frag.getIDFromLocant("N'"));
//		frag.relabelSuffixLocants(parentFrag);
//		assertEquals("Can find locant N' in frag: ID = 4", 4, frag.getIDFromLocant("N'"));
//		frag = sBuilder.build("C", idManager);
//		assertEquals("Can find locant 1 in frag: ID = 5", 5, frag.getIDFromLocant("1"));
//		frag.relabelSuffixLocants(parentFrag);
//		assertEquals("Can't find locant 1 in frag", 0, frag.getIDFromLocant("1"));
//	}

	@Test
	public void testPickUpIndicatedHydrogen() throws StructureBuildingException {
		//TODO: this, properly.
		SMILESFragmentBuilder sBuilder = new SMILESFragmentBuilder();
		Fragment pyrrole = sBuilder.build("NC=CC=C", new IDManager());
		pyrrole.pickUpIndicatedHydrogen();
		//assertEquals("Pyrrole is 1H", "1", pyrrole.)
	}

	@Test
	public void testConvertHighOrderBondsToSpareValencies() throws Exception {
		SMILESFragmentBuilder sBuilder = new SMILESFragmentBuilder();
		Fragment naphthalene = sBuilder.build("C1=CC=CC2=CC=CC=C12", new IDManager());
		CycleDetector.assignWhetherAtomsAreInCycles(naphthalene);
		naphthalene.convertHighOrderBondsToSpareValencies();
		for(Atom a : naphthalene.getAtomList()) {
			assertEquals("All atoms have one sv", 1, a.getSpareValency());
		}
		for(Bond b : naphthalene.getBondList()) {
			assertEquals("All bonds are of order 1", 1, b.getOrder());
		}
	}

	@Test
	public void testConvertSpareValenciesToDoubleBonds() throws Exception {
		SMILESFragmentBuilder sBuilder = new SMILESFragmentBuilder();
		IDManager idManager = new IDManager();
		Fragment dhp = sBuilder.build("C1=CCC=CC1", idManager);
		CycleDetector.assignWhetherAtomsAreInCycles(dhp);
		dhp.convertHighOrderBondsToSpareValencies();
		dhp.convertSpareValenciesToDoubleBonds();
		for(Atom a : dhp.getAtomList()) {
			assertEquals("All atoms have no sv", 0, a.getSpareValency());
		}
		Fragment funnydiene = sBuilder.build("C(=C)C=C", idManager);
		CycleDetector.assignWhetherAtomsAreInCycles(funnydiene);
		funnydiene.convertHighOrderBondsToSpareValencies();
		funnydiene.convertSpareValenciesToDoubleBonds();
		for(Atom a : funnydiene.getAtomList()) {
			assertEquals("All atoms have no sv", 0, a.getSpareValency());
		}
		Fragment naphthalene = sBuilder.build("C1=CC=CC2=CC=CC=C12", idManager);
		CycleDetector.assignWhetherAtomsAreInCycles(naphthalene);
		naphthalene.convertHighOrderBondsToSpareValencies();
		naphthalene.convertSpareValenciesToDoubleBonds();
		for(Atom a : naphthalene.getAtomList()) {
			assertEquals("All atoms have no sv", 0, a.getSpareValency());
		}
		Fragment pentalene = sBuilder.build("C12C(=CC=C1)C=CC=2", idManager);
		CycleDetector.assignWhetherAtomsAreInCycles(pentalene);
		pentalene.convertHighOrderBondsToSpareValencies();
		pentalene.convertSpareValenciesToDoubleBonds();
		for(Atom a : pentalene.getAtomList()) {
			assertEquals("All atoms have no sv", 0, a.getSpareValency());
		}

	}

	@Test
	public void testGetAtomNeighbours() throws Exception {
		SMILESFragmentBuilder sBuilder = new SMILESFragmentBuilder();
		Fragment naphthalene = sBuilder.build("C1=CC=CC2=CC=CC=C12", new IDManager());
		assertEquals("Atom 1 has two neighbours",
				2, naphthalene.getAtomNeighbours(naphthalene.getAtomByID(1)).size());
		assertEquals("Atom 5 has three neighbours",
				3, naphthalene.getAtomNeighbours(naphthalene.getAtomByID(5)).size());
	}

}
