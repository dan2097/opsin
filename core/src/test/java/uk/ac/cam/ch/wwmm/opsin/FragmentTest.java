package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.List;

import org.junit.Before;
import org.junit.Test;

import static junit.framework.Assert.*;
import static org.mockito.Mockito.mock;
import static uk.ac.cam.ch.wwmm.opsin.XmlDeclarations.*;

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
		a3.addLocant("3");
		frag.addAtom(a3);
		fm.createBond(frag.getAtomByID(2), frag.getAtomByID(3), 1);
		assertEquals("propane", 3, frag.getChainLength());
		Atom a4 =new Atom(4, "C", frag);
		frag.addAtom(a4);
		a4.addLocant("4");
		fm.createBond(frag.getAtomByID(2), frag.getAtomByID(4), 1);
		assertEquals("isobutane", 3, frag.getChainLength());
		fm.removeBond(frag.findBond(a2, a4));
		fm.createBond(a3, a4, 1);
		assertEquals("butane", 4, frag.getChainLength());
	}


	@Test
	public void testRelabelSuffixLocants() throws StructureBuildingException {
		frag = fm.buildSMILES("C(N)N");
		assertEquals("Can't find locant N in frag", 0, frag.getIDFromLocant("N"));
		assertEquals("Can't find locant N' in frag", 0, frag.getIDFromLocant("N'"));
		FragmentTools.assignElementLocants(frag, new ArrayList<Fragment>());
		int idN = frag.getIDFromLocant("N");
		int idNprime = frag.getIDFromLocant("N'");
		if ((idN==2 && idNprime==3) || idN==3 && idNprime==2){
		}
		else{
			fail("Locants misassigned");
		}
	}
	
	@Test
	public void testLabelCarbamimidamido() throws StructureBuildingException {
		frag =  fm.buildSMILES("C(N)(=N)N-", NONCARBOXYLICACID_TYPE_VAL, NONE_LABELS_VAL);
		FragmentTools.assignElementLocants(frag, new ArrayList<Fragment>());
		assertEquals("Can find locant N in frag: ID = 4", 4, frag.getIDFromLocant("N"));
		assertEquals("Can find locant N' in frag: ID = 2", 2, frag.getIDFromLocant("N'"));
		assertEquals("Can find locant N'' in frag: ID = 3", 3, frag.getIDFromLocant("N''"));
	}
	
	@Test
	public void testLabelHydrazonoHydrazide() throws StructureBuildingException {
		frag =  fm.buildSMILES("C(=NN)NN" , NONCARBOXYLICACID_TYPE_VAL, NONE_LABELS_VAL);
		FragmentTools.assignElementLocants(frag, new ArrayList<Fragment>());
		assertEquals("Can find locant N in frag: ID = 4", 4, frag.getIDFromLocant("N"));
		assertEquals("Can find locant N' in frag: ID = 5", 5, frag.getIDFromLocant("N'"));
		assertEquals("Can find locant N'' in frag: ID = 2", 2, frag.getIDFromLocant("N''"));
		assertEquals("Can find locant N''' in frag: ID = 3", 3, frag.getIDFromLocant("N'''"));
	}
	
	
	@Test
	public void testLabelCarbonimidoyl() throws StructureBuildingException {
		frag =  fm.buildSMILES("C(=N)" , ACIDSTEM_TYPE_VAL, NONE_LABELS_VAL);
		frag.addOutAtom(frag.getFirstAtom(), 1, true);
		frag.addOutAtom(frag.getFirstAtom(), 1, true);
		FragmentTools.assignElementLocants(frag, new ArrayList<Fragment>());
		assertEquals("Can find locant N in frag: ID = 2", 2, frag.getIDFromLocant("N"));
		assertEquals("Can find locant N in frag: ID = 1", 1, frag.getIDFromLocant("C"));
	}
	
	@Test
	public void testLabelHydrazonicAmide() throws StructureBuildingException {
		frag =  fm.buildSMILES("C", ACIDSTEM_TYPE_VAL, NONE_LABELS_VAL);
		Fragment suffixfrag =  fm.buildSMILES("[R](N)=NN", SUFFIX_TYPE_VAL, NONE_LABELS_VAL);
		List<Fragment> suffixes = new ArrayList<Fragment>();
		suffixes.add(suffixfrag);
		FragmentTools.assignElementLocants(frag, suffixes);
		fm.incorporateFragment(suffixfrag, frag);
		assertEquals("Can find locant N in frag: ID = 3", 3, frag.getIDFromLocant("N"));
		assertEquals("Can find locant N'' in frag: ID = 5", 5, frag.getIDFromLocant("N''"));
		
		//DEVIATION From systematic behaviour
		assertEquals("Can find locant N' in frag: ID = 5", 5, frag.getIDFromLocant("N'"));
	}
	
	@Test
	public void testLabelHydrazonate() throws StructureBuildingException {
		frag =  fm.buildSMILES("C", ACIDSTEM_TYPE_VAL, NONE_LABELS_VAL);
		Fragment suffixfrag =  fm.buildSMILES("[R]([O-])=NN", SUFFIX_TYPE_VAL, NONE_LABELS_VAL);
		List<Fragment> suffixes = new ArrayList<Fragment>();
		suffixes.add(suffixfrag);
		FragmentTools.assignElementLocants(frag, suffixes);
		fm.incorporateFragment(suffixfrag, frag);
		assertEquals("Can find locant N' in frag: ID = 5", 5, frag.getIDFromLocant("N'"));
		
		//DEVIATION From systematic behaviour
		assertEquals("Can find locant N in frag: ID = 5", 5, frag.getIDFromLocant("N"));
	}
	
	@Test
	public void testLabelHexanDiamide() throws StructureBuildingException {
		frag =  fm.buildSMILES("CCCCCC", CHAIN_TYPE_VAL, "1/2/3/4/5/6");
		Fragment suffixfrag1 =  fm.buildSMILES("[R]N", SUFFIX_TYPE_VAL, NONE_LABELS_VAL);
		Fragment suffixfrag2 =  fm.buildSMILES("[R]N", SUFFIX_TYPE_VAL, NONE_LABELS_VAL);
		List<Fragment> suffixes = new ArrayList<Fragment>();
		suffixes.add(suffixfrag1);
		suffixes.add(suffixfrag2);
		FragmentTools.assignElementLocants(frag, suffixes);
		fm.incorporateFragment(suffixfrag1, frag);
		fm.incorporateFragment(suffixfrag2, frag);
		assertEquals("Can find locant N in frag: ID = 8", 8, frag.getIDFromLocant("N"));
		assertEquals("Can find locant N' in frag: ID = 10", 10, frag.getIDFromLocant("N'"));
	}
	
	@Test
	public void testLabelDiimidooxalicDiamide() throws StructureBuildingException {
		frag =  fm.buildSMILES("CC", ACIDSTEM_TYPE_VAL, "1/2");
		Fragment suffixfrag1 =  fm.buildSMILES("[R](N)=N", SUFFIX_TYPE_VAL, NONE_LABELS_VAL);
		Fragment suffixfrag2 =  fm.buildSMILES("[R](N)=N", SUFFIX_TYPE_VAL, NONE_LABELS_VAL);
		List<Fragment> suffixes = new ArrayList<Fragment>();
		suffixes.add(suffixfrag1);
		suffixes.add(suffixfrag2);
		FragmentTools.assignElementLocants(frag, suffixes);
		fm.incorporateFragment(suffixfrag1, frag);
		fm.incorporateFragment(suffixfrag2, frag);
		assertEquals("Can find locant N in frag: ID = 4", 4, frag.getIDFromLocant("N"));
		assertEquals("Can find locant N' in frag: ID = 7", 7, frag.getIDFromLocant("N'"));
		assertEquals("Can find locant N'' in frag: ID = 5", 5, frag.getIDFromLocant("N''"));
		assertEquals("Can find locant N''' in frag: ID = 8", 8, frag.getIDFromLocant("N'''"));
	}
	
	
	@Test
	public void testLabelHydrazinecarbohydrazide() throws StructureBuildingException {
		frag =  fm.buildSMILES("NN", SIMPLEGROUP_TYPE_VAL, "1/2");
		Fragment suffix =  fm.buildSMILES("[R]C(=O)NN", SUFFIX_TYPE_VAL, "/X///");
		List<Fragment> suffixes = new ArrayList<Fragment>();
		suffixes.add(suffix);
		FragmentTools.assignElementLocants(frag, suffixes);
		fm.incorporateFragment(suffix, frag);
		assertEquals("Can find locant N in frag: ID = 6", 6, frag.getIDFromLocant("N"));
		assertEquals("Can find locant N' in frag: ID = 7", 7, frag.getIDFromLocant("N'"));
		assertEquals("Can't find locant N'' in frag", 0, frag.getIDFromLocant("N''"));
		assertEquals("Can't find locant C in frag", 0, frag.getIDFromLocant("C"));
	}
	
	
	@Test
	public void testLabelCarbonicDihydrazide() throws StructureBuildingException {
		frag =  fm.buildSMILES("C(=O)(NN)NN", NONCARBOXYLICACID_TYPE_VAL, NONE_LABELS_VAL);
		FragmentTools.assignElementLocants(frag, new ArrayList<Fragment>());
		int idN = frag.getIDFromLocant("N");
		int idNprime = frag.getIDFromLocant("N'");
		int idNprime2 = frag.getIDFromLocant("N''");
		int idNprime3 = frag.getIDFromLocant("N'''");
		if ((idN==3 && idNprime==4 && idNprime2==5 && idNprime3==6) || (idN==5 && idNprime==6 && idNprime2==3 && idNprime3==4)){
		}
		else{
			fail("Locants misassigned");
		}
		assertEquals("Can't find locant C in frag", 0, frag.getIDFromLocant("C"));
	}
	
	@Test
	public void testLabelSulfonoThioate() throws StructureBuildingException {
		frag = fm.buildSMILES("C");
		Fragment suffix = fm.buildSMILES("[R]S(=O)(=O)S", SUFFIX_TYPE_VAL, "/X///");
		List<Fragment> suffixes = new ArrayList<Fragment>();
		suffixes.add(suffix);
		FragmentTools.assignElementLocants(frag, suffixes);
		fm.incorporateFragment(suffix, frag);
		assertEquals("Can find locant S in frag: ID = 6", 6, frag.getIDFromLocant("S"));
		assertEquals("Can't find locant S' in frag", 0, frag.getIDFromLocant("S'"));
		int idO = frag.getIDFromLocant("O");
		int idOprime = frag.getIDFromLocant("O'");
		if ((idO==4 && idOprime==5)|| (idO==5 && idOprime==4)){
		}
		else{
			fail("Locants misassigned");
		}
	}
	
	@Test
	public void testLabelAcetoanilide() throws StructureBuildingException {
		frag = fm.buildSMILES("CC");
		Fragment suffix = fm.buildSMILES("[*](=O)Nc1ccccc1", SUFFIX_TYPE_VAL, "///1'/2'/3'/4'/5'/6'");
		List<Fragment> suffixes = new ArrayList<Fragment>();
		suffixes.add(suffix);
		FragmentTools.assignElementLocants(frag, suffixes);
		fm.incorporateFragment(suffix, frag);
		assertEquals("Can find locant N in frag: ID = 5", 5, frag.getIDFromLocant("N"));
	}
	
	
	
	@Test
	public void testLabelPyridine() throws StructureBuildingException {
		frag =  fm.buildSMILES("n1ccccc1", RING_TYPE_VAL, "1/2/3/4/5/6");
		FragmentTools.assignElementLocants(frag, new ArrayList<Fragment>());
		assertEquals("Can find locant N in frag: ID = 1", 1, frag.getIDFromLocant("N"));
		assertEquals("Can't find locant C in frag", 0, frag.getIDFromLocant("C"));
	}
	
	
	@Test
	public void testLabelPiperazine() throws StructureBuildingException {
		frag =  fm.buildSMILES("N1CCNCC1", RING_TYPE_VAL, "1/2/3/4/5/6");
		FragmentTools.assignElementLocants(frag, new ArrayList<Fragment>());
		int idN = frag.getIDFromLocant("N");
		int idNprime = frag.getIDFromLocant("N'");
		if ((idN==1 && idNprime==4) || (idN==4 && idNprime==1)){
		}
		else{
			fail("Locants misassigned");
		}
		assertEquals("Can't find locant C in frag", 0, frag.getIDFromLocant("C"));
	}
	
	@Test
	public void testLabelCarboximidohydrazide() throws StructureBuildingException {
		frag = fm.buildSMILES("c1ccccc1");
		Fragment suffix = fm.buildSMILES("[R]C(=N)NN", SUFFIX_TYPE_VAL, "/X//1'/2'");
		List<Fragment> suffixes = new ArrayList<Fragment>();
		suffixes.add(suffix);
		FragmentTools.assignElementLocants(frag, suffixes);
		fm.incorporateFragment(suffix, frag);
		assertEquals("Can find locant N in frag: ID = 10", 10, frag.getIDFromLocant("N"));
		assertEquals("Can find locant N' in frag: ID = 11", 11, frag.getIDFromLocant("N'"));
		assertEquals("Can find locant N'' in frag: ID = 9", 9, frag.getIDFromLocant("N''"));
	}

	@Test
	public void testPickUpIndicatedHydrogen() throws StructureBuildingException {
		SMILESFragmentBuilder sBuilder = new SMILESFragmentBuilder();
		Fragment pyrrole = sBuilder.build("N1C=CC=C1", fm);
		FragmentTools.convertHighOrderBondsToSpareValencies(pyrrole);
		FragmentTools.pickUpIndicatedHydrogens(pyrrole);
		assertEquals("Pyrrole has 1 indicated hydrogen", 1, pyrrole.getIndicatedHydrogen().size());
		assertEquals("..and the indicated hydrogen is on the nitrogen", pyrrole.getFirstAtom(), pyrrole.getIndicatedHydrogen().get(0));
	}

	@Test
	public void testConvertHighOrderBondsToSpareValencies() throws Exception {
		SMILESFragmentBuilder sBuilder = new SMILESFragmentBuilder();
		Fragment naphthalene = sBuilder.build("C1=CC=CC2=CC=CC=C12", fm);
		FragmentTools.convertHighOrderBondsToSpareValencies(naphthalene);
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
		FragmentTools.convertHighOrderBondsToSpareValencies(dhp);
		FragmentTools.convertSpareValenciesToDoubleBonds(dhp);
		for(Atom a : dhp.getAtomList()) {
			assertEquals("All atoms have no sv", false, a.hasSpareValency());
		}
		Fragment funnydiene = sBuilder.build("C(=C)C=C", fm);
		FragmentTools.convertHighOrderBondsToSpareValencies(funnydiene);
		FragmentTools.convertSpareValenciesToDoubleBonds(funnydiene);
		for(Atom a : funnydiene.getAtomList()) {
			assertEquals("All atoms have no sv", false, a.hasSpareValency());
		}
		Fragment naphthalene = sBuilder.build("C1=CC=CC2=CC=CC=C12", fm);
		FragmentTools.convertHighOrderBondsToSpareValencies(naphthalene);
		FragmentTools.convertSpareValenciesToDoubleBonds(naphthalene);
		for(Atom a : naphthalene.getAtomList()) {
			assertEquals("All atoms have no sv", false, a.hasSpareValency());
		}
		Fragment pentalene = sBuilder.build("C12C(=CC=C1)C=CC=2", fm);
		FragmentTools.convertHighOrderBondsToSpareValencies(pentalene);
		FragmentTools.convertSpareValenciesToDoubleBonds(pentalene);
		for(Atom a : pentalene.getAtomList()) {
			assertEquals("All atoms have no sv", false, a.hasSpareValency());
		}

	}

	@Test
	public void testGetAtomNeighbours() throws Exception {
		SMILESFragmentBuilder sBuilder = new SMILESFragmentBuilder();
		Fragment naphthalene = sBuilder.build("C1=CC=CC2=CC=CC=C12", fm);
		assertEquals("Atom 1 has two neighbours",
				2, naphthalene.getIntraFragmentAtomNeighbours(naphthalene.getAtomByID(1)).size());
		assertEquals("Atom 5 has three neighbours",
				3, naphthalene.getIntraFragmentAtomNeighbours(naphthalene.getAtomByID(5)).size());
	}

}
