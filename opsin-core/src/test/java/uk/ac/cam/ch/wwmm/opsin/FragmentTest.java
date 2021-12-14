package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;
import static uk.ac.cam.ch.wwmm.opsin.XmlDeclarations.*;

public class FragmentTest {

	private Fragment frag;
	private FragmentManager fm;

	@BeforeEach
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
		assertNotNull(frag.getAtomList(), "Has atom list");
	}

	@Test
	public void testAddAtom() {
		assertEquals(0, frag.getAtomCount(), "Has no atoms");
		frag.addAtom(new Atom(1, ChemEl.C, frag));
		assertEquals(1, frag.getAtomCount(), "Now has one atom");
	}

	@Test
	public void testAddBond() {
		frag.addAtom(new Atom(1, ChemEl.C, frag));
		frag.addAtom(new Atom(2, ChemEl.C, frag));
		assertEquals(0, frag.getBondSet().size(), "Has no bonds");
		fm.createBond(frag.getAtomByID(1), frag.getAtomByID(2), 1);
		assertEquals(1, frag.getBondSet().size(), "Now has one bond");
	}

	@Test
	public void testImportFrag() throws StructureBuildingException {
		Fragment frag1 = fm.buildSMILES("CC");
		Fragment frag2 = fm.buildSMILES("CC");
		assertEquals(2, frag1.getAtomCount(), "Fragment has two atoms");
		assertEquals(1, frag1.getBondSet().size(), "Fragment has one bond");
		fm.incorporateFragment(frag2, frag1);
		assertEquals(4, frag1.getAtomCount(), "Fragment now has four atoms");
		assertEquals(2, frag1.getBondSet().size(), "Fragment now has two bonds");
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
		assertEquals(10, frag.getIDFromLocant("a"), "Locant a has ID 10");
		assertEquals(20, frag.getIDFromLocant("silly"), "Locant silly has ID 20");
		assertEquals(0, frag.getIDFromLocant("42"), "Locant 42 is not present");
	}

	@Test
	public void testGetAtomByLocant()  {
		Atom atom1 = new Atom(10, ChemEl.C, frag);
		atom1.addLocant("a");
		frag.addAtom(atom1);
		Atom atom2 = new Atom(20, ChemEl.C, frag);
		atom2.addLocant("silly");
		frag.addAtom(atom2);
		assertEquals(atom1, frag.getAtomByLocant("a"), "Locant a gets atom1");
		assertEquals(atom2, frag.getAtomByLocant("silly"), "Locant silly gets atom2");
		assertNull(frag.getAtomByLocant("42"), "Locant 42 is not present");
	}

	@Test
	public void testGetAtomByID() {
		Atom atom1 = new Atom(10, ChemEl.C, frag);
		frag.addAtom(atom1);
		Atom atom2 = new Atom(20, ChemEl.C, frag);
		frag.addAtom(atom2);
		assertEquals(atom1, frag.getAtomByID(10), "ID 10 gets atom1");
		assertEquals(atom2, frag.getAtomByID(20), "ID 20 gets atom2");
		assertNull(frag.getAtomByID(42), "ID 42 is not present");
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
		assertNotNull(b, "Found a bond");
		assertEquals(2, b.getOrder(), "..a double bond");
		b = frag.findBond(3, 1);
		assertNotNull(b, "Found a different bond");
		assertEquals(3, b.getOrder(), "..a triple bond");
		b = frag.findBond(2, 3);
		assertNull(b, "Don't find non-existent bonds");
	}

	@Test
	public void testGetChainLength() {
		assertEquals(0, frag.getChainLength(), "No chain");
		Atom a1 =new Atom(1, ChemEl.C, frag);
		a1.addLocant("1");
		frag.addAtom(a1);
		assertEquals(1, frag.getChainLength(), "Methane");
		Atom a2 =new Atom(2, ChemEl.C, frag);
		a2.addLocant("2");
		frag.addAtom(a2);
		fm.createBond(frag.getAtomByID(1), frag.getAtomByID(2), 1);
		assertEquals(2, frag.getChainLength(), "ethane");
		Atom a3 =new Atom(3, ChemEl.C, frag);
		a3.addLocant("3");
		frag.addAtom(a3);
		fm.createBond(frag.getAtomByID(2), frag.getAtomByID(3), 1);
		assertEquals(3, frag.getChainLength(), "propane");
		Atom a4 =new Atom(4, ChemEl.C, frag);
		frag.addAtom(a4);
		a4.addLocant("4");
		fm.createBond(frag.getAtomByID(2), frag.getAtomByID(4), 1);
		assertEquals(3, frag.getChainLength(), "isobutane");
		fm.removeBond(a2.getBondToAtom(a4));
		fm.createBond(a3, a4, 1);
		assertEquals(4, frag.getChainLength(), "butane");
	}


	@Test
	public void testRelabelSuffixLocants() throws StructureBuildingException {
		frag = fm.buildSMILES("C(N)N");
		assertEquals(0, frag.getIDFromLocant("N"), "Can't find locant N in frag");
		assertEquals(0, frag.getIDFromLocant("N'"), "Can't find locant N' in frag");
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
		assertEquals(4, frag.getIDFromLocant("N"), "Can find locant N in frag: ID = 4");
		assertEquals(2, frag.getIDFromLocant("N'"), "Can find locant N' in frag: ID = 2");
		assertEquals(3, frag.getIDFromLocant("N''"), "Can find locant N'' in frag: ID = 3");
	}
	
	@Test
	public void testLabelHydrazonoHydrazide() throws StructureBuildingException {
		frag =  fm.buildSMILES("C(=NN)NN" , NONCARBOXYLICACID_TYPE_VAL, NONE_LABELS_VAL);
		FragmentTools.assignElementLocants(frag, new ArrayList<Fragment>());
		assertEquals(4, frag.getIDFromLocant("N"), "Can find locant N in frag: ID = 4");
		assertEquals(5, frag.getIDFromLocant("N'"), "Can find locant N' in frag: ID = 5");
		assertEquals(2, frag.getIDFromLocant("N''"), "Can find locant N'' in frag: ID = 2");
		assertEquals(3, frag.getIDFromLocant("N'''"), "Can find locant N''' in frag: ID = 3");
	}
	
	
	@Test
	public void testLabelCarbonimidoyl() throws StructureBuildingException {
		frag =  fm.buildSMILES("C(=N)" , ACIDSTEM_TYPE_VAL, NONE_LABELS_VAL);
		frag.addOutAtom(frag.getFirstAtom(), 1, true);
		frag.addOutAtom(frag.getFirstAtom(), 1, true);
		FragmentTools.assignElementLocants(frag, new ArrayList<Fragment>());
		assertEquals(2, frag.getIDFromLocant("N"), "Can find locant N in frag: ID = 2");
		assertEquals(1, frag.getIDFromLocant("C"), "Can find locant N in frag: ID = 1");
	}
	
	@Test
	public void testLabelHydrazonicAmide() throws StructureBuildingException {
		frag =  fm.buildSMILES("C", ACIDSTEM_TYPE_VAL, NONE_LABELS_VAL);
		Fragment suffixfrag =  fm.buildSMILES("[R](N)=NN", SUFFIX_TYPE_VAL, NONE_LABELS_VAL);
		List<Fragment> suffixes = new ArrayList<Fragment>();
		suffixes.add(suffixfrag);
		FragmentTools.assignElementLocants(frag, suffixes);
		fm.incorporateFragment(suffixfrag, frag);
		assertEquals(3, frag.getIDFromLocant("N"), "Can find locant N in frag: ID = 3");
		assertEquals(5, frag.getIDFromLocant("N''"), "Can find locant N'' in frag: ID = 5");
		
		//DEVIATION From systematic behaviour
		assertEquals(5, frag.getIDFromLocant("N'"), "Can find locant N' in frag: ID = 5");
	}
	
	@Test
	public void testLabelHydrazonate() throws StructureBuildingException {
		frag =  fm.buildSMILES("C", ACIDSTEM_TYPE_VAL, NONE_LABELS_VAL);
		Fragment suffixfrag =  fm.buildSMILES("[R]([O-])=NN", SUFFIX_TYPE_VAL, NONE_LABELS_VAL);
		List<Fragment> suffixes = new ArrayList<Fragment>();
		suffixes.add(suffixfrag);
		FragmentTools.assignElementLocants(frag, suffixes);
		fm.incorporateFragment(suffixfrag, frag);
		assertEquals(5, frag.getIDFromLocant("N'"), "Can find locant N' in frag: ID = 5");
		
		//DEVIATION From systematic behaviour
		assertEquals(5, frag.getIDFromLocant("N"), "Can find locant N in frag: ID = 5");
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
		assertEquals(8, frag.getIDFromLocant("N"), "Can find locant N in frag: ID = 8");
		assertEquals(10, frag.getIDFromLocant("N'"), "Can find locant N' in frag: ID = 10");
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
		assertEquals(4, frag.getIDFromLocant("N"), "Can find locant N in frag: ID = 4");
		assertEquals(7, frag.getIDFromLocant("N'"), "Can find locant N' in frag: ID = 7");
		assertEquals(5, frag.getIDFromLocant("N''"), "Can find locant N'' in frag: ID = 5");
		assertEquals(8, frag.getIDFromLocant("N'''"), "Can find locant N''' in frag: ID = 8");
	}
	
	
	@Test
	public void testLabelHydrazinecarbohydrazide() throws StructureBuildingException {
		frag =  fm.buildSMILES("NN", SIMPLEGROUP_TYPE_VAL, "1/2");
		Fragment suffix =  fm.buildSMILES("[R]C(=O)NN", SUFFIX_TYPE_VAL, "/X///");
		List<Fragment> suffixes = new ArrayList<Fragment>();
		suffixes.add(suffix);
		FragmentTools.assignElementLocants(frag, suffixes);
		fm.incorporateFragment(suffix, frag);
		assertEquals(6, frag.getIDFromLocant("N"), "Can find locant N in frag: ID = 6");
		assertEquals(7, frag.getIDFromLocant("N'"), "Can find locant N' in frag: ID = 7");
		assertEquals(0, frag.getIDFromLocant("N''"), "Can't find locant N'' in frag");
		assertEquals(0, frag.getIDFromLocant("C"), "Can't find locant C in frag");
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
		assertEquals(0, frag.getIDFromLocant("C"), "Can't find locant C in frag");
	}
	
	@Test
	public void testLabelSulfonoThioate() throws StructureBuildingException {
		frag = fm.buildSMILES("C");
		Fragment suffix = fm.buildSMILES("[R]S(=O)(=O)S", SUFFIX_TYPE_VAL, "/X///");
		List<Fragment> suffixes = new ArrayList<Fragment>();
		suffixes.add(suffix);
		FragmentTools.assignElementLocants(frag, suffixes);
		fm.incorporateFragment(suffix, frag);
		assertEquals(6, frag.getIDFromLocant("S"), "Can find locant S in frag: ID = 6");
		assertEquals(0, frag.getIDFromLocant("S'"), "Can't find locant S' in frag");
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
		assertEquals(5, frag.getIDFromLocant("N"), "Can find locant N in frag: ID = 5");
	}
	
	
	
	@Test
	public void testLabelPyridine() throws StructureBuildingException {
		frag =  fm.buildSMILES("n1ccccc1", RING_TYPE_VAL, "1/2/3/4/5/6");
		FragmentTools.assignElementLocants(frag, new ArrayList<Fragment>());
		assertEquals(1, frag.getIDFromLocant("N"), "Can find locant N in frag: ID = 1");
		assertEquals(0, frag.getIDFromLocant("C"), "Can't find locant C in frag");
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
		assertEquals(0, frag.getIDFromLocant("C"), "Can't find locant C in frag");
	}
	
	@Test
	public void testLabelCarboximidohydrazide() throws StructureBuildingException {
		frag = fm.buildSMILES("c1ccccc1");
		Fragment suffix = fm.buildSMILES("[R]C(=N)NN", SUFFIX_TYPE_VAL, "/X//1'/2'");
		List<Fragment> suffixes = new ArrayList<Fragment>();
		suffixes.add(suffix);
		FragmentTools.assignElementLocants(frag, suffixes);
		fm.incorporateFragment(suffix, frag);
		assertEquals(10, frag.getIDFromLocant("N"), "Can find locant N in frag: ID = 10");
		assertEquals(11, frag.getIDFromLocant("N'"), "Can find locant N' in frag: ID = 11");
		assertEquals(9, frag.getIDFromLocant("N''"), "Can find locant N'' in frag: ID = 9");
	}

	@Test
	public void testIndicatedHydrogen() throws StructureBuildingException {
		Fragment pyrrole = fm.buildSMILES("[nH]1cccc1");
		assertEquals(1, pyrrole.getIndicatedHydrogen().size(), "Pyrrole has 1 indicated hydrogen");
		assertEquals(pyrrole.getFirstAtom(), pyrrole.getIndicatedHydrogen().get(0), "..and the indicated hydrogen is on the nitrogen");
	}

	@Test
	public void testSpareValenciesOnAromaticAtoms() throws StructureBuildingException{
		Fragment naphthalene = fm.buildSMILES("c1cccc2ccccc12");
		for(Atom a : naphthalene) {
			assertEquals(true, a.hasSpareValency(), "All atoms have sv");
		}
		for(Bond b : naphthalene.getBondSet()) {
			assertEquals(1, b.getOrder(), "All bonds are of order 1");
		}
	}

	@Test
	public void testConvertSpareValenciesToDoubleBonds() throws StructureBuildingException{
		Fragment dhp = fm.buildSMILES("c1cCccC1");
		FragmentTools.convertSpareValenciesToDoubleBonds(dhp);
		for(Atom a : dhp) {
			assertEquals(false, a.hasSpareValency(), "All atoms have no sv");
		}
		Fragment funnydiene = fm.buildSMILES("C(=C)C=C");
		FragmentTools.convertSpareValenciesToDoubleBonds(funnydiene);
		for(Atom a : funnydiene) {
			assertEquals(false, a.hasSpareValency(), "All atoms have no sv");
		}
		Fragment naphthalene = fm.buildSMILES("c1cccc2ccccc12");
		FragmentTools.convertSpareValenciesToDoubleBonds(naphthalene);
		for(Atom a : naphthalene) {
			assertEquals(false, a.hasSpareValency(), "All atoms have no sv");
		}
		Fragment pentalene = fm.buildSMILES("c12c(ccc1)ccc2");
		for(Atom a : pentalene) {
			assertEquals(true, a.hasSpareValency(), "All atoms have sv");
		}
		FragmentTools.convertSpareValenciesToDoubleBonds(pentalene);
		for(Atom a : pentalene) {
			assertEquals(false, a.hasSpareValency(), "All atoms have no sv");
		}

	}

	@Test
	public void testGetAtomNeighbours() throws StructureBuildingException {
		Fragment naphthalene = fm.buildSMILES("C1=CC=CC2=CC=CC=C12");
		assertEquals(2, naphthalene.getIntraFragmentAtomNeighbours(naphthalene.getAtomByID(1)).size(),
				"Atom 1 has two neighbours");
		assertEquals(3, naphthalene.getIntraFragmentAtomNeighbours(naphthalene.getAtomByID(5)).size(),
				"Atom 5 has three neighbours");
	}
	
	@Test
	public void testIsCharacteristicAtomSuffix() throws StructureBuildingException{
		Fragment parent = fm.buildSMILES("CC");
		Fragment suffix = fm.buildSMILES("N", SUFFIX_TYPE_VAL, NONE_LABELS_VAL);
		fm.incorporateFragment(suffix, suffix.getFirstAtom(), parent, parent.getFirstAtom(), 1);
		List<Atom> parentAtoms = parent.getAtomList();
		assertFalse(FragmentTools.isCharacteristicAtom(parentAtoms.get(0)));
		assertFalse(FragmentTools.isCharacteristicAtom(parentAtoms.get(1)));
		assertTrue(FragmentTools.isCharacteristicAtom(parentAtoms.get(2)));
	}
	
	@Test
	public void testIsCharacteristicAtomAldehyde() throws StructureBuildingException{
		Fragment parent = fm.buildSMILES("CC");
		Fragment suffix = fm.buildSMILES("O", SUFFIX_TYPE_VAL, NONE_LABELS_VAL);
		fm.incorporateFragment(suffix, suffix.getFirstAtom(), parent, parent.getFirstAtom(), 2);
		List<Atom> parentAtoms = parent.getAtomList();
		parentAtoms.get(1).setProperty(Atom.ISALDEHYDE, true);
		assertFalse(FragmentTools.isCharacteristicAtom(parentAtoms.get(0)));
		assertTrue(FragmentTools.isCharacteristicAtom(parentAtoms.get(1)));
		assertTrue(FragmentTools.isCharacteristicAtom(parentAtoms.get(2)));
	}
	
	@Test
	public void testIsCharacteristicAtomFunctionalAtom() throws StructureBuildingException{
		Fragment parent = fm.buildSMILES("CC(=O)[O-]");
		List<Atom> parentAtoms = parent.getAtomList();
		parent.addFunctionalAtom(parentAtoms.get(3));
		for (int i = 0; i < parentAtoms.size() - 1; i++) {
			assertFalse(FragmentTools.isCharacteristicAtom(parentAtoms.get(i)));
		}
		assertTrue(FragmentTools.isCharacteristicAtom(parentAtoms.get(parentAtoms.size() - 1)));
	}
	
	@Test
	public void testIsCharacteristicAtomHydroxy() throws StructureBuildingException{
		List<Atom> phenolAtoms = fm.buildSMILES("Oc1ccccc1").getAtomList();
		assertTrue(FragmentTools.isCharacteristicAtom(phenolAtoms.get(0)));
		for (int i = 1; i < phenolAtoms.size(); i++) {
			assertFalse(FragmentTools.isCharacteristicAtom(phenolAtoms.get(i)));
		}
	}

}
