package uk.ac.cam.ch.wwmm.opsin;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.fail;

import java.util.Collections;
import java.util.List;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import uk.ac.cam.ch.wwmm.opsin.BondStereo.BondStereoValue;

public class SMILESWriterTest {
	
	private FragmentManager fm;

	@BeforeEach
	public void setup(){
		IDManager idManager = new IDManager();
		fm = new FragmentManager(new SMILESFragmentBuilder(idManager), idManager);
	}

	@Test
	public void testRoundTrip1() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("C");
		fm.makeHydrogensExplicit();
		String smiles = SMILESWriter.generateSmiles(f);
		assertEquals("C", smiles);
	}

	@Test
	public void testRoundTrip2() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("C#N");
		fm.makeHydrogensExplicit();
		String smiles = SMILESWriter.generateSmiles(f);
		assertEquals("C#N", smiles);
	}
	
	@Test
	public void testRoundTrip3() throws StructureBuildingException {
		Fragment f = fm.buildSMILES(StringTools.multiplyString("C",200));
		fm.makeHydrogensExplicit();
		String smiles = SMILESWriter.generateSmiles(f);
		assertEquals(StringTools.multiplyString("C",200), smiles);
	}

	@Test
	public void testRoundTrip4() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("O=C=O");
		fm.makeHydrogensExplicit();
		String smiles = SMILESWriter.generateSmiles(f);
		assertEquals("O=C=O", smiles);
	}

	@Test
	public void testRoundTrip5() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("CCN(CC)CC");
		fm.makeHydrogensExplicit();
		String smiles = SMILESWriter.generateSmiles(f);
		assertEquals("CCN(CC)CC", smiles);
	}

	@Test
	public void testRoundTrip6() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("CC(=O)O");
		fm.makeHydrogensExplicit();
		String smiles = SMILESWriter.generateSmiles(f);
		assertEquals("CC(=O)O", smiles);
	}

	@Test
	public void testRoundTrip7() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("C1CCCCC1");
		fm.makeHydrogensExplicit();
		String smiles = SMILESWriter.generateSmiles(f);
		assertEquals("C1CCCCC1", smiles);
	}
	
	@Test
	public void testRoundTrip8() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("C1=CC=CC=C1");
		fm.makeHydrogensExplicit();
		String smiles = SMILESWriter.generateSmiles(f);
		assertEquals("C1=CC=CC=C1", smiles);
	}

	@Test
	public void testRoundTrip9() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("NC(Cl)(Br)C(=O)O");
		fm.makeHydrogensExplicit();
		String smiles = SMILESWriter.generateSmiles(f);
		assertEquals("NC(Cl)(Br)C(=O)O", smiles);
	}

	@Test
	public void testRoundTrip10() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("[NH4+].[Cl-].F.[He-2]");
		fm.makeHydrogensExplicit();
		String smiles = SMILESWriter.generateSmiles(f);
		assertEquals("[NH4+].[Cl-].F.[He-2]", smiles);
	}
	
	@Test
	public void testRoundTrip11() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("[NH4+].[Cl-].F.[He-2]");
		List<Atom> atomList = f.getAtomList();
		Collections.reverse(atomList);
		f.reorderAtomCollection(atomList);
		fm.makeHydrogensExplicit();
		String smiles = SMILESWriter.generateSmiles(f);
		assertEquals("[He-2].F.[Cl-].[NH4+]", smiles);
	}
	
	@Test
	public void testRoundTrip12() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("CCO.N=O.C#N");
		fm.makeHydrogensExplicit();
		String smiles = SMILESWriter.generateSmiles(f);
		assertEquals("CCO.N=O.C#N", smiles);
	}
	
	@Test
	public void testOrganic1() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("[S]");
		String smiles = SMILESWriter.generateSmiles(f);
		assertEquals("[S]", smiles);
	}
	
	@Test
	public void testOrganic2() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("[S][H]");
		String smiles = SMILESWriter.generateSmiles(f);
		assertEquals("[SH]", smiles);
	}
	
	@Test
	public void testOrganic3() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("[S]([H])[H]");
		String smiles = SMILESWriter.generateSmiles(f);
		assertEquals("S", smiles);
	}
	
	@Test
	public void testOrganic4() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("[S]([H])([H])[H]");
		String smiles = SMILESWriter.generateSmiles(f);
		assertEquals("[SH3]", smiles);
	}
	
	@Test
	public void testOrganic5() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("[S]([H])([H])([H])[H]");
		String smiles = SMILESWriter.generateSmiles(f);
		assertEquals("[SH4]", smiles);
	}
	
	@Test
	public void testOrganic6() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("S(F)(F)(F)F");
		String smiles = SMILESWriter.generateSmiles(f);
		assertEquals("S(F)(F)(F)F", smiles);
	}
	
	@Test
	public void testOrganic7() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("S([H])(F)(F)(F)(F)F");
		String smiles = SMILESWriter.generateSmiles(f);
		assertEquals("S(F)(F)(F)(F)F", smiles);
	}
	
	@Test
	public void testOrganic8() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("S([H])([H])(F)(F)(F)F");
		String smiles = SMILESWriter.generateSmiles(f);
		assertEquals("[SH2](F)(F)(F)F", smiles);
	}
	
	@Test
	public void testOrganic9() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("S(F)(F)(F)(F)(F)(F)F");
		String smiles = SMILESWriter.generateSmiles(f);
		assertEquals("S(F)(F)(F)(F)(F)(F)F", smiles);
	}
	
	@Test
	public void testOrganic10() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("[I]([H])([H])[H]");
		String smiles = SMILESWriter.generateSmiles(f);
		assertEquals("[IH3]", smiles);
	}
	
	@Test
	public void testCharged1() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("[CH3+]");
		fm.makeHydrogensExplicit();
		String smiles = SMILESWriter.generateSmiles(f);
		assertEquals("[CH3+]", smiles);
	}
	
	@Test
	public void testCharged2() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("[Mg+2]");
		fm.makeHydrogensExplicit();
		String smiles = SMILESWriter.generateSmiles(f);
		assertEquals("[Mg+2]", smiles);
	}
	@Test
	public void testCharged3() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("[BH4-]");
		fm.makeHydrogensExplicit();
		String smiles = SMILESWriter.generateSmiles(f);
		assertEquals("[BH4-]", smiles);
	}
	
	@Test
	public void testCharged4() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("[O-2]");
		fm.makeHydrogensExplicit();
		String smiles = SMILESWriter.generateSmiles(f);
		assertEquals("[O-2]", smiles);
	}
	
	@Test
	public void testIsotope() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("[15NH3]");
		fm.makeHydrogensExplicit();
		String smiles = SMILESWriter.generateSmiles(f);
		assertEquals("[15NH3]", smiles);
	}
	
	@Test
	public void testRGroup1() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("[R]CC[R]");
		fm.makeHydrogensExplicit();
		String smiles = SMILESWriter.generateSmiles(f);
		assertEquals("*CC*", smiles);
	}
	
	@Test
	public void testRGroup2() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("[H][R]");
		fm.makeHydrogensExplicit();
		String smiles = SMILESWriter.generateSmiles(f);
		assertEquals("*[H]", smiles);
	}
	
	@Test
	public void testRingOpeningsGreaterThan10() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("C12=C3C4=C5C6=C1C7=C8C9=C1C%10=C%11C(=C29)C3=C2C3=C4C4=C5C5=C9C6=C7C6=C7C8=C1C1=C8C%10=C%10C%11=C2C2=C3C3=C4C4=C5C5=C%11C%12=C(C6=C95)C7=C1C1=C%12C5=C%11C4=C3C3=C5C(=C81)C%10=C23");
		fm.makeHydrogensExplicit();
		String smiles = SMILESWriter.generateSmiles(f);
		assertEquals("C12=C3C4=C5C6=C1C1=C7C8=C9C%10=C%11C(=C28)C3=C3C2=C4C4=C5C5=C8C6=C1C1=C6C7=C9C9=C7C%10=C%10C%11=C3C3=C2C2=C4C4=C5C5=C%11C%12=C(C1=C85)C6=C9C9=C%12C%12=C%11C4=C2C2=C%12C(=C79)C%10=C32", smiles);
	}
	
	@Test
	public void testHydrogenNotBondedToAnyNonHydrogen1() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("[H-].[H+]");
		String smiles = SMILESWriter.generateSmiles(f);
		assertEquals("[H-].[H+]", smiles);
	}
	
	@Test
	public void testHydrogenNotBondedToAnyNonHydrogen2() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("[H][H]");
		String smiles = SMILESWriter.generateSmiles(f);
		assertEquals("[H][H]", smiles);
	}
	
	@Test
	public void testHydrogenNotBondedToAnyNonHydrogen3() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("[2H][H]");
		String smiles = SMILESWriter.generateSmiles(f);
		assertEquals("[2H][H]", smiles);
	}
	
	@Test
	public void testHydrogenNotBondedToAnyNonHydrogen4() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("[H]B1[H]B([H])[H]1");
		String smiles = SMILESWriter.generateSmiles(f);
		assertEquals("B1[H]B[H]1", smiles);
	}
	
	@Test
	public void testTetrahedralChirality1() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("N[C@@H](F)C");
		fm.makeHydrogensExplicit();
		String smiles = SMILESWriter.generateSmiles(f);
		assertEquals("N[C@@H](F)C", smiles);
	}

	@Test
	public void testTetrahedralChirality2() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("N[C@H](F)C");
		fm.makeHydrogensExplicit();
		String smiles = SMILESWriter.generateSmiles(f);
		assertEquals("N[C@H](F)C", smiles);
	}

	@Test
	public void testTetrahedralChirality3() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("C2.N1.F3.[C@@H]231");
		fm.makeHydrogensExplicit();
		String smiles = SMILESWriter.generateSmiles(f);
		assertEquals("C[C@H](F)N", smiles);
	}

	@Test
	public void testTetrahedralChirality4() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("[C@@H]231.C2.N1.F3");
		fm.makeHydrogensExplicit();
		String smiles = SMILESWriter.generateSmiles(f);
		assertEquals("[C@H](C)(N)F", smiles);
	}

	@Test
	public void testTetrahedralChirality5() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("[C@@H](Cl)1[C@H](C)(F).Br1");
		fm.makeHydrogensExplicit();
		String smiles = SMILESWriter.generateSmiles(f);
		assertEquals("[C@H](Cl)([C@H](C)F)Br", smiles);
	}

	@Test
	public void testTetrahedralChirality6() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("I[C@@](Cl)(Br)F");
		fm.makeHydrogensExplicit();
		String smiles = SMILESWriter.generateSmiles(f);
		assertEquals("I[C@@](Cl)(Br)F", smiles);
	}
	

	@Test
	public void testTetrahedralChirality7() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("C[S@](N)=O");
		fm.makeHydrogensExplicit();
		String smiles = SMILESWriter.generateSmiles(f);
		assertEquals("C[S@](N)=O", smiles);
	}
	
	@Test
	public void testDoubleBondSupport1() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("C/C=C/C");
		fm.makeHydrogensExplicit();
		String smiles = SMILESWriter.generateSmiles(f);
		if (!smiles.equals("C/C=C/C") && !smiles.equals("C\\C=C\\C")){
			fail(smiles +" did not correspond to one of the expected SMILES strings");
		}
	}
	
	
	@Test
	public void testDoubleBondSupport2() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("C/C=C\\C");
		fm.makeHydrogensExplicit();
		String smiles = SMILESWriter.generateSmiles(f);
		if (!smiles.equals("C/C=C\\C") && !smiles.equals("C\\C=C/C")){
			fail(smiles +" did not correspond to one of the expected SMILES strings");
		}
	}
	
	
	@Test
	public void testDoubleBondSupport3() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("C/C=C\\C=C/C");
		fm.makeHydrogensExplicit();
		String smiles = SMILESWriter.generateSmiles(f);
		if (!smiles.equals("C/C=C\\C=C/C") && !smiles.equals("C\\C=C/C=C\\C")){
			fail(smiles +" did not correspond to one of the expected SMILES strings");
		}
	}

	@Test
	public void testDoubleBondSupport4() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("ClC(C(=O)[O-])=CC(=CC(=O)[O-])Cl");
		fm.makeHydrogensExplicit();
		f.findBond(2, 6).setBondStereoElement(new Atom[]{f.getAtomByID(1), f.getAtomByID(2), f.getAtomByID(6), f.getAtomByID(7)}, BondStereoValue.TRANS);
		f.findBond(7, 8).setBondStereoElement(new Atom[]{f.getAtomByID(12), f.getAtomByID(7), f.getAtomByID(8), f.getAtomByID(9)}, BondStereoValue.TRANS);
		String smiles = SMILESWriter.generateSmiles(f);
		if (!smiles.equals("Cl\\C(\\C(=O)[O-])=C\\C(=C/C(=O)[O-])\\Cl") && !smiles.equals("Cl/C(/C(=O)[O-])=C/C(=C\\C(=O)[O-])/Cl")){
			fail(smiles +" did not correspond to one of the expected SMILES strings");
		}
	}
	
	@Test
	public void testDoubleBondSupport5() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("C/C=N\\O");
		fm.makeHydrogensExplicit();
		String smiles = SMILESWriter.generateSmiles(f);
		if (!smiles.equals("C/C=N\\O") && !smiles.equals("C\\C=N/O")){
			fail(smiles +" did not correspond to one of the expected SMILES strings");
		}
	}
	
	@Test
	public void testDoubleBondSupport6() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("O=C(/C=C(C(O)=O)\\C=C/C(O)=O)O");
		fm.makeHydrogensExplicit();
		String smiles = SMILESWriter.generateSmiles(f);
		if (!smiles.equals("O=C(/C=C(/C(O)=O)\\C=C/C(O)=O)O") && !smiles.equals("O=C(\\C=C(\\C(O)=O)/C=C\\C(O)=O)O")){
			fail(smiles +" did not correspond to one of the expected SMILES strings");
		}
	}
	
	@Test
	public void testDoubleBondSupport7() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("C(=C(C=CC(=O)O)C(=O)O)C(=O)O");
		fm.makeHydrogensExplicit();
		f.findBond(1, 2).setBondStereoElement(new Atom[]{f.getAtomByID(11), f.getAtomByID(1), f.getAtomByID(2), f.getAtomByID(8)}, BondStereoValue.TRANS);
		f.findBond(3, 4).setBondStereoElement(new Atom[]{f.getAtomByID(2), f.getAtomByID(3), f.getAtomByID(4), f.getAtomByID(5)}, BondStereoValue.TRANS);
		String smiles = SMILESWriter.generateSmiles(f);
		if (!smiles.equals("C(=C(/C=C/C(=O)O)\\C(=O)O)/C(=O)O") && !smiles.equals("C(=C(\\C=C\\C(=O)O)/C(=O)O)\\C(=O)O")){
			fail(smiles +" did not correspond to one of the expected SMILES strings");
		}
	}
	
	@Test
	public void testDoubleBondSupport8() throws StructureBuildingException {
		//hydrogen on the nitrogen must be explicit!
		Fragment f = fm.buildSMILES("[H]/N=C(\\N)/O");
		fm.makeHydrogensExplicit();
		String smiles = SMILESWriter.generateSmiles(f);
		if (!smiles.equals("[H]/N=C(\\N)/O") && !smiles.equals("[H]\\N=C(/N)\\O")){
			fail(smiles +" did not correspond to one of the expected SMILES strings");
		}
	}
	
	@Test
	public void testDoubleBondSupport9() throws StructureBuildingException {
		//Adjacent double bonds need to be processed first to ensure we don't assign inconsistent slashes two the exocyclic double bond 
		Fragment f = fm.buildSMILES("CN/1CCC2=C(\\C(=C/C(/C(=C1)/C(=O)[O-])=C\\OC)\\C1=CC=C(C=C1)C)C=C(C(=C2)OC)OC");
		fm.makeHydrogensExplicit();
		String smiles = SMILESWriter.generateSmiles(f);
		if (!smiles.equals("CN/1CCC2=C(\\C(=C/C(/C(=C1)/C(=O)[O-])=C\\OC)\\C1=CC=C(C=C1)C)C=C(C(=C2)OC)OC")){
			fail(smiles +" did not correspond to one of the expected SMILES strings");
		}
	}
	
	@Test
	public void testLabelling1() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("CCC", "", XmlDeclarations.NONE_LABELS_VAL);
		for (Atom a : f) {
			assertEquals(0, a.getLocants().size());
		}
		
		Fragment f2 = fm.buildSMILES("CCC", "", "");
		for (Atom a : f2) {
			assertEquals(0, a.getLocants().size());
		}
	}
	
	@Test
	public void testLabelling2() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("CCC", "", "1/2,alpha,2'/");
		List<Atom> atoms = f.getAtomList();
		assertEquals(1, atoms.get(0).getLocants().size());
		assertEquals(3, atoms.get(1).getLocants().size());
		assertEquals(0, atoms.get(2).getLocants().size());
		
		assertEquals("1", atoms.get(0).getLocants().get(0));
		assertEquals("2", atoms.get(1).getLocants().get(0));
		assertEquals("alpha", atoms.get(1).getLocants().get(1));
		assertEquals("2'", atoms.get(1).getLocants().get(2));
	}

	@Test
	public void testAtomClassEmitted() throws StructureBuildingException {
		Fragment f = fm.buildSMILES("CCO");
		fm.makeHydrogensExplicit();
		assertEquals("CCO", SMILESWriter.generateSmiles(f));
		for (Atom a : f) {
			a.setProperty(Atom.ATOM_CLASS, a.getID());
		}
		assertEquals("[CH3:1][CH2:2][OH:3]", SMILESWriter.generateSmiles(f));
	}
}
