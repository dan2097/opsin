package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;

import org.junit.Before;
import org.junit.Ignore;
import org.junit.Test;

import static junit.framework.Assert.*;
import static org.mockito.Mockito.mock;

public class FragmentTest {

	private Fragment frag;
	private FragmentManager fm = new FragmentManager(new SMILESFragmentBuilder(), mock(CMLFragmentBuilder.class), new IDManager());

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
		frag.addAtom(new Atom(1, "C", frag));
		assertEquals("Now has one atom", 1, frag.getAtomList().size());
	}

	@Test
	public void testAddBond() throws StructureBuildingException {
		frag.addAtom(new Atom(1, "C", frag));
		frag.addAtom(new Atom(2, "C", frag));
		assertEquals("Has no bonds", 0, frag.getBondSet().size());
		fm.createBond(frag.getAtomByID(1), frag.getAtomByID(2), 1);
		assertEquals("Now has one bond", 1, frag.getBondSet().size());
	}

	@Test
	public void testImportFrag() throws StructureBuildingException {
		Fragment frag1 = fm.buildSMILES("CC");
		Fragment frag2 = fm.buildSMILES("CC");
		assertEquals("Fragment has two atoms", 2, frag1.getAtomList().size());
		assertEquals("Fragment has one bond", 1, frag1.getBondSet().size());
		fm.incorporateFragment(frag2, frag1);
		assertEquals("Fragment now has four atoms", 4, frag1.getAtomList().size());
		assertEquals("Fragment now has two bonds", 2, frag1.getBondSet().size());
	}

	@Test
	public void testGetIDFromLocant() throws StructureBuildingException {
		Atom atom = new Atom(10, "C", frag);
		atom.addLocant("a");
		frag.addAtom(atom);
		atom = new Atom(20, "C", frag);
		atom.addLocant("silly");
		frag.addAtom(atom);
		assertEquals("Locant a has ID 10", 10, frag.getIDFromLocant("a"));
		assertEquals("Locant silly has ID 20", 20, frag.getIDFromLocant("silly"));
		assertEquals("Locant 42 is not present", 0, frag.getIDFromLocant("42"));
	}

	@Test
	public void testGetAtomByLocant() throws StructureBuildingException {
		Atom atom1 = new Atom(10, "C", frag);
		atom1.addLocant("a");
		frag.addAtom(atom1);
		Atom atom2 = new Atom(20, "C", frag);
		atom2.addLocant("silly");
		frag.addAtom(atom2);
		assertEquals("Locant a gets atom1", atom1, frag.getAtomByLocant("a"));
		assertEquals("Locant silly gets atom2", atom2, frag.getAtomByLocant("silly"));
		assertNull("Locant 42 is not present", frag.getAtomByLocant("42"));
	}

	@Test
	public void testGetAtomByID() throws StructureBuildingException {
		Atom atom1 = new Atom(10, "C", frag);
		frag.addAtom(atom1);
		Atom atom2 = new Atom(20, "C", frag);
		frag.addAtom(atom2);
		assertEquals("ID 10 gets atom1", atom1, frag.getAtomByID(10));
		assertEquals("ID 20 gets atom2", atom2, frag.getAtomByID(20));
		assertNull("ID 42 is not present", frag.getAtomByID(42));
	}

	@Test
	public void testFindBond() throws StructureBuildingException {
		frag.addAtom(new Atom(1, "C", frag));
		frag.addAtom(new Atom(2, "C", frag));
		frag.addAtom(new Atom(3, "N", frag));
		frag.addAtom(new Atom(4, "O", frag));
		fm.createBond(frag.getAtomByID(2), frag.getAtomByID(4), 2);
		fm.createBond(frag.getAtomByID(1), frag.getAtomByID(2), 1);
		fm.createBond(frag.getAtomByID(1), frag.getAtomByID(3), 3);
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
		Atom a1 =new Atom(1, "C", frag);
		a1.addLocant("1");
		frag.addAtom(a1);
		assertEquals("Methane", 1, frag.getChainLength());
		Atom a2 =new Atom(2, "C", frag);
		a2.addLocant("2");
		frag.addAtom(a2);
		fm.createBond(frag.getAtomByID(1), frag.getAtomByID(2), 1);
		assertEquals("ethane", 2, frag.getChainLength());
		Atom a3 =new Atom(3, "C", frag);
		a1.addLocant("3");
		frag.addAtom(a3);
		fm.createBond(frag.getAtomByID(2), frag.getAtomByID(3), 1);
		assertEquals("propane", 3, frag.getChainLength());

	}


	@Test
	public void testRelabelSuffixLocants() throws StructureBuildingException {
		frag = fm.buildSMILES("C(N)N");
		assertEquals("Can't find locant N in frag", 0, frag.getIDFromLocant("N"));
		assertEquals("Can't find locant N in frag", 0, frag.getIDFromLocant("N'"));
		FragmentTools.assignElementLocants(frag, new ArrayList<Fragment>());
		assertEquals("Can find locant N in frag: ID = 2", 2, frag.getIDFromLocant("N"));
		assertEquals("Can find locant N in frag: ID = 3", 3, frag.getIDFromLocant("N'"));
	}
	
	@Test
	@Ignore
	public void testLabelCarbamimidamido() throws StructureBuildingException {
		frag =  fm.buildSMILES("NC(=N)N-");
		FragmentTools.assignElementLocants(frag, new ArrayList<Fragment>());
		assertEquals("Can find locant N in frag: ID = 4", 4, frag.getIDFromLocant("N"));
		assertEquals("Can find locant N in frag: ID = 1", 1, frag.getIDFromLocant("N'"));
		assertEquals("Can find locant N in frag: ID = 3", 3, frag.getIDFromLocant("N''"));
	}
	
	@Test
	@Ignore
	public void testLabelHydrazonoHydrazide() throws StructureBuildingException {
		frag =  fm.buildSMILES("C(=NN)NN");
		FragmentTools.assignElementLocants(frag, new ArrayList<Fragment>());
		assertEquals("Can find locant N in frag: ID = 4", 4, frag.getIDFromLocant("N"));
		assertEquals("Can find locant N in frag: ID = 5", 5, frag.getIDFromLocant("N'"));
		assertEquals("Can find locant N in frag: ID = 2", 2, frag.getIDFromLocant("N''"));
		assertEquals("Can find locant N in frag: ID = 3", 3, frag.getIDFromLocant("N'''"));
	}

	@Test
	public void testPickUpIndicatedHydrogen() throws StructureBuildingException {
		//TODO: this, properly.
		SMILESFragmentBuilder sBuilder = new SMILESFragmentBuilder();
		Fragment pyrrole = sBuilder.build("NC=CC=C", fm);
		pyrrole.pickUpIndicatedHydrogen();
		//assertEquals("Pyrrole is 1H", "1", pyrrole.)
	}

	@Test
	public void testConvertHighOrderBondsToSpareValencies() throws Exception {
		SMILESFragmentBuilder sBuilder = new SMILESFragmentBuilder();
		Fragment naphthalene = sBuilder.build("C1=CC=CC2=CC=CC=C12", fm);
		CycleDetector.assignWhetherAtomsAreInCycles(naphthalene);
		naphthalene.convertHighOrderBondsToSpareValencies();
		for(Atom a : naphthalene.getAtomList()) {
			assertEquals("All atoms have sv", true, a.hasSpareValency());
		}
		for(Bond b : naphthalene.getBondSet()) {
			assertEquals("All bonds are of order 1", 1, b.getOrder());
		}
	}

	@Test
	public void testConvertSpareValenciesToDoubleBonds() throws Exception {
		SMILESFragmentBuilder sBuilder = new SMILESFragmentBuilder();
		Fragment dhp = sBuilder.build("C1=CCC=CC1", fm);
		CycleDetector.assignWhetherAtomsAreInCycles(dhp);
		dhp.convertHighOrderBondsToSpareValencies();
		dhp.convertSpareValenciesToDoubleBonds();
		for(Atom a : dhp.getAtomList()) {
			assertEquals("All atoms have no sv", false, a.hasSpareValency());
		}
		Fragment funnydiene = sBuilder.build("C(=C)C=C", fm);
		CycleDetector.assignWhetherAtomsAreInCycles(funnydiene);
		funnydiene.convertHighOrderBondsToSpareValencies();
		funnydiene.convertSpareValenciesToDoubleBonds();
		for(Atom a : funnydiene.getAtomList()) {
			assertEquals("All atoms have no sv", false, a.hasSpareValency());
		}
		Fragment naphthalene = sBuilder.build("C1=CC=CC2=CC=CC=C12", fm);
		CycleDetector.assignWhetherAtomsAreInCycles(naphthalene);
		naphthalene.convertHighOrderBondsToSpareValencies();
		naphthalene.convertSpareValenciesToDoubleBonds();
		for(Atom a : naphthalene.getAtomList()) {
			assertEquals("All atoms have no sv", false, a.hasSpareValency());
		}
		Fragment pentalene = sBuilder.build("C12C(=CC=C1)C=CC=2", fm);
		CycleDetector.assignWhetherAtomsAreInCycles(pentalene);
		pentalene.convertHighOrderBondsToSpareValencies();
		pentalene.convertSpareValenciesToDoubleBonds();
		for(Atom a : pentalene.getAtomList()) {
			assertEquals("All atoms have no sv", false, a.hasSpareValency());
		}

	}

	@Test
	public void testGetAtomNeighbours() throws Exception {
		SMILESFragmentBuilder sBuilder = new SMILESFragmentBuilder();
		Fragment naphthalene = sBuilder.build("C1=CC=CC2=CC=CC=C12", fm);
		assertEquals("Atom 1 has two neighbours",
				2, naphthalene.getAtomNeighbours(naphthalene.getAtomByID(1)).size());
		assertEquals("Atom 5 has three neighbours",
				3, naphthalene.getAtomNeighbours(naphthalene.getAtomByID(5)).size());
	}

}
