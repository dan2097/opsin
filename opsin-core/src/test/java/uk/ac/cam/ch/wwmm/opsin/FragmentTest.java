package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.List;

import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.*;
import static uk.ac.cam.ch.wwmm.opsin.XmlDeclarations.*;

public class FragmentTest {

	private Fragment frag;
	private FragmentManager fm;

	@Before
	public void setUp(){
		IDManager idManager = new IDManager();
		fm = new FragmentManager(new SMILESFragmentBuilder(idManager), idManager);
		try {
			frag = fm.buildSMILES("");
		} catch (StructureBuildingException e) {
			throw new RuntimeException(e);
		}
	}

	@Test
	public void testFragment() {
		assertNotNull("Has atom list", frag.getAtomList());
	}

	@Test
	public void testAddAtom() {
		assertEquals("Has no atoms", 0, frag.getAtomCount());
		frag.addAtom(new Atom(1, ChemEl.C, frag));
		assertEquals("Now has one atom", 1, frag.getAtomCount());
	}

	@Test
	public void testAddBond() {
		frag.addAtom(new Atom(1, ChemEl.C, frag));
		frag.addAtom(new Atom(2, ChemEl.C, frag));
		assertEquals("Has no bonds", 0, frag.getBondSet().size());
		fm.createBond(frag.getAtomByID(1), frag.getAtomByID(2), 1);
		assertEquals("Now has one bond", 1, frag.getBondSet().size());
	}

	@Test
	public void testImportFrag() throws StructureBuildingException {
		Fragment frag1 = fm.buildSMILES("CC");
		Fragment frag2 = fm.buildSMILES("CC");
		assertEquals("Fragment has two atoms", 2, frag1.getAtomCount());
		assertEquals("Fragment has one bond", 1, frag1.getBondSet().size());
		fm.incorporateFragment(frag2, frag1);
		assertEquals("Fragment now has four atoms", 4, frag1.getAtomCount());
		assertEquals("Fragment now has two bonds", 2, frag1.getBondSet().size());
	}
	
	@Test
	public void testImportFragWithIntraFragBonds1() throws StructureBuildingException {
		Fragment frag1 = fm.buildSMILES("C");
		Fragment frag2 = fm.buildSMILES("C");
		fm.createBond(frag1.getFirstAtom(), frag2.getFirstAtom(), 1);
		assertEquals(0, frag1.getBondSet().size());
		assertEquals(0, frag2.getBondSet().size());
		assertEquals(1, fm.getInterFragmentBonds(frag1).size());
		assertEquals(1, fm.getInterFragmentBonds(frag2).size());
		fm.incorporateFragment(frag2, frag1);
		assertEquals(1, frag1.getBondSet().size());
		assertEquals(0, frag2.getBondSet().size());
		assertEquals(0, fm.getInterFragmentBonds(frag1).size());
	}
	
	@Test
	public void testImportFragWithIntraFragBonds2() throws StructureBuildingException {
		Fragment frag1 = fm.buildSMILES("C");
		Fragment frag2 = fm.buildSMILES("C");
		Fragment frag3 = fm.buildSMILES("C");
		fm.createBond(frag2.getFirstAtom(), frag3.getFirstAtom(), 1);
		assertEquals(0, frag1.getBondSet().size());
		assertEquals(0, frag2.getBondSet().size());
		assertEquals(0, frag3.getBondSet().size());
		assertEquals(0, fm.getInterFragmentBonds(frag1).size());
		assertEquals(1, fm.getInterFragmentBonds(frag2).size());
		assertEquals(1, fm.getInterFragmentBonds(frag3).size());
		fm.incorporateFragment(frag2, frag1);
		assertEquals(0, frag1.getBondSet().size());
		assertEquals(0, frag2.getBondSet().size());
		assertEquals(0, frag3.getBondSet().size());
		assertEquals(1, fm.getInterFragmentBonds(frag1).size());
		assertEquals(1, fm.getInterFragmentBonds(frag3).size());
	}

	@Test
	public void testGetIDFromLocant() {
		Atom atom = new Atom(10, ChemEl.C, frag);
		atom.addLocant("a");
		frag.addAtom(atom);
		atom = new Atom(20, ChemEl.C, frag);
		atom.addLocant("silly");
		frag.addAtom(atom);
		assertEquals("Locant a has ID 10", 10, frag.getIDFromLocant("a"));
		assertEquals("Locant silly has ID 20", 20, frag.getIDFromLocant("silly"));
		assertEquals("Locant 42 is not present", 0, frag.getIDFromLocant("42"));
	}

	@Test
	public void testGetAtomByLocant()  {
		Atom atom1 = new Atom(10, ChemEl.C, frag);
		atom1.addLocant("a");
		frag.addAtom(atom1);
		Atom atom2 = new Atom(20, ChemEl.C, frag);
		atom2.addLocant("silly");
		frag.addAtom(atom2);
		assertEquals("Locant a gets atom1", atom1, frag.getAtomByLocant("a"));
		assertEquals("Locant silly gets atom2", atom2, frag.getAtomByLocant("silly"));
		assertNull("Locant 42 is not present", frag.getAtomByLocant("42"));
	}

	@Test
	public void testGetAtomByID() {
		Atom atom1 = new Atom(10, ChemEl.C, frag);
		frag.addAtom(atom1);
		Atom atom2 = new Atom(20, ChemEl.C, frag);
		frag.addAtom(atom2);
		assertEquals("ID 10 gets atom1", atom1, frag.getAtomByID(10));
		assertEquals("ID 20 gets atom2", atom2, frag.getAtomByID(20));
		assertNull("ID 42 is not present", frag.getAtomByID(42));
	}

	@Test
	public void testFindBond() {
		frag.addAtom(new Atom(1, ChemEl.C, frag));
		frag.addAtom(new Atom(2, ChemEl.C, frag));
		frag.addAtom(new Atom(3, ChemEl.N, frag));
		frag.addAtom(new Atom(4, ChemEl.O, frag));
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
	public void testGetChainLength() {
		assertEquals("No chain", 0, frag.getChainLength());
		Atom a1 =new Atom(1, ChemEl.C, frag);
		a1.addLocant("1");
		frag.addAtom(a1);
		assertEquals("Methane", 1, frag.getChainLength());
		Atom a2 =new Atom(2, ChemEl.C, frag);
		a2.addLocant("2");
		frag.addAtom(a2);
		fm.createBond(frag.getAtomByID(1), frag.getAtomByID(2), 1);
		assertEquals("ethane", 2, frag.getChainLength());
		Atom a3 =new Atom(3, ChemEl.C, frag);
		a3.addLocant("3");
		frag.addAtom(a3);
		fm.createBond(frag.getAtomByID(2), frag.getAtomByID(3), 1);
		assertEquals("propane", 3, frag.getChainLength());
		Atom a4 =new Atom(4, ChemEl.C, frag);
		frag.addAtom(a4);
		a4.addLocant("4");
		fm.createBond(frag.getAtomByID(2), frag.getAtomByID(4), 1);
		assertEquals("isobutane", 3, frag.getChainLength());
		fm.removeBond(a2.getBondToAtom(a4));
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
	public void testIndicatedHydrogen() throws StructureBuildingException {
		Fragment pyrrole = fm.buildSMILES("[nH]1cccc1");
		assertEquals("Pyrrole has 1 indicated hydrogen", 1, pyrrole.getIndicatedHydrogen().size());
		assertEquals("..and the indicated hydrogen is on the nitrogen", pyrrole.getFirstAtom(), pyrrole.getIndicatedHydrogen().get(0));
	}

	@Test
	public void testSpareValenciesOnAromaticAtoms() throws StructureBuildingException{
		Fragment naphthalene = fm.buildSMILES("c1cccc2ccccc12");
		for(Atom a : naphthalene.getAtomList()) {
			assertEquals("All atoms have sv", true, a.hasSpareValency());
		}
		for(Bond b : naphthalene.getBondSet()) {
			assertEquals("All bonds are of order 1", 1, b.getOrder());
		}
	}

	@Test
	public void testConvertSpareValenciesToDoubleBonds() throws StructureBuildingException{
		Fragment dhp = fm.buildSMILES("c1cCccC1");
		FragmentTools.convertSpareValenciesToDoubleBonds(dhp);
		for(Atom a : dhp.getAtomList()) {
			assertEquals("All atoms have no sv", false, a.hasSpareValency());
		}
		Fragment funnydiene = fm.buildSMILES("C(=C)C=C");
		FragmentTools.convertSpareValenciesToDoubleBonds(funnydiene);
		for(Atom a : funnydiene.getAtomList()) {
			assertEquals("All atoms have no sv", false, a.hasSpareValency());
		}
		Fragment naphthalene = fm.buildSMILES("c1cccc2ccccc12");
		FragmentTools.convertSpareValenciesToDoubleBonds(naphthalene);
		for(Atom a : naphthalene.getAtomList()) {
			assertEquals("All atoms have no sv", false, a.hasSpareValency());
		}
		Fragment pentalene = fm.buildSMILES("c12c(ccc1)ccc2");
		for(Atom a : pentalene.getAtomList()) {
			assertEquals("All atoms have sv", true, a.hasSpareValency());
		}
		FragmentTools.convertSpareValenciesToDoubleBonds(pentalene);
		for(Atom a : pentalene.getAtomList()) {
			assertEquals("All atoms have no sv", false, a.hasSpareValency());
		}

	}

	@Test
	public void testGetAtomNeighbours() throws StructureBuildingException{
		Fragment naphthalene = fm.buildSMILES("C1=CC=CC2=CC=CC=C12");
		assertEquals("Atom 1 has two neighbours",
				2, naphthalene.getIntraFragmentAtomNeighbours(naphthalene.getAtomByID(1)).size());
		assertEquals("Atom 5 has three neighbours",
				3, naphthalene.getIntraFragmentAtomNeighbours(naphthalene.getAtomByID(5)).size());
	}

}
