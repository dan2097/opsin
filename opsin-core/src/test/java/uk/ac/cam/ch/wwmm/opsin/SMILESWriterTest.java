package uk.ac.cam.ch.wwmm.opsin;
import static org.junit.Assert.assertEquals;
import static org.mockito.Mockito.mock;
import static junit.framework.Assert.*;

import java.util.Collections;
import java.util.List;

import org.junit.Before;
import org.junit.Test;

import uk.ac.cam.ch.wwmm.opsin.BondStereo.BondStereoValue;

public class SMILESWriterTest {
	
	BuildState state;
	@Before
	public void setup(){
		state = new BuildState(mock(NameToStructureConfig.class), new SMILESFragmentBuilder(), mock(CMLFragmentBuilder.class));
	}

	@Test
	public void testRoundTrip1() throws StructureBuildingException {
		Fragment f = state.fragManager.buildSMILES("C");
		StructureBuilder.makeHydrogensExplicit(state);
		String smiles = new SMILESWriter(f).generateSmiles();
		assertEquals("C", smiles);
	}

	@Test
	public void testRoundTrip2() throws StructureBuildingException {
		Fragment f = state.fragManager.buildSMILES("C#N");
		StructureBuilder.makeHydrogensExplicit(state);
		String smiles = new SMILESWriter(f).generateSmiles();
		assertEquals("C#N", smiles);
	}
	
	@Test
	public void testRoundTrip3() throws StructureBuildingException {
		Fragment f = state.fragManager.buildSMILES(StringTools.multiplyString("C",200));
		StructureBuilder.makeHydrogensExplicit(state);
		String smiles = new SMILESWriter(f).generateSmiles();
		assertEquals(StringTools.multiplyString("C",200), smiles);
	}

	@Test
	public void testRoundTrip4() throws StructureBuildingException {
		Fragment f = state.fragManager.buildSMILES("O=C=O");
		StructureBuilder.makeHydrogensExplicit(state);
		String smiles = new SMILESWriter(f).generateSmiles();
		assertEquals("O=C=O", smiles);
	}

	@Test
	public void testRoundTrip5() throws StructureBuildingException {
		Fragment f = state.fragManager.buildSMILES("CCN(CC)CC");
		StructureBuilder.makeHydrogensExplicit(state);
		String smiles = new SMILESWriter(f).generateSmiles();
		assertEquals("CCN(CC)CC", smiles);
	}

	@Test
	public void testRoundTrip6() throws StructureBuildingException {
		Fragment f = state.fragManager.buildSMILES("CC(=O)O");
		StructureBuilder.makeHydrogensExplicit(state);
		String smiles = new SMILESWriter(f).generateSmiles();
		assertEquals("CC(=O)O", smiles);
	}

	@Test
	public void testRoundTrip7() throws StructureBuildingException {
		Fragment f = state.fragManager.buildSMILES("C1CCCCC1");
		StructureBuilder.makeHydrogensExplicit(state);
		String smiles = new SMILESWriter(f).generateSmiles();
		assertEquals("C1CCCCC1", smiles);
	}
	
	@Test
	public void testRoundTrip8() throws StructureBuildingException {
		Fragment f = state.fragManager.buildSMILES("C1=CC=CC=C1");
		StructureBuilder.makeHydrogensExplicit(state);
		String smiles = new SMILESWriter(f).generateSmiles();
		assertEquals("C1=CC=CC=C1", smiles);
	}

	@Test
	public void testRoundTrip9() throws StructureBuildingException {
		Fragment f = state.fragManager.buildSMILES("NC(Cl)(Br)C(=O)O");
		StructureBuilder.makeHydrogensExplicit(state);
		String smiles = new SMILESWriter(f).generateSmiles();
		assertEquals("NC(Cl)(Br)C(=O)O", smiles);
	}

	@Test
	public void testRoundTrip10() throws StructureBuildingException {
		Fragment f = state.fragManager.buildSMILES("[NH4+].[Cl-].F.[He-2]");
		StructureBuilder.makeHydrogensExplicit(state);
		String smiles = new SMILESWriter(f).generateSmiles();
		assertEquals("[NH4+].[Cl-].F.[He-2]", smiles);
	}
	
	@Test
	public void testRoundTrip11() throws StructureBuildingException {
		Fragment f = state.fragManager.buildSMILES("[NH4+].[Cl-].F.[He-2]");
		List<Atom> atomList = f.getAtomList();
		Collections.reverse(atomList);
		f.reorderAtomCollection(atomList);
		StructureBuilder.makeHydrogensExplicit(state);
		String smiles = new SMILESWriter(f).generateSmiles();
		assertEquals("[He-2].F.[Cl-].[NH4+]", smiles);
	}
	
	@Test
	public void testRoundTrip12() throws StructureBuildingException {
		Fragment f = state.fragManager.buildSMILES("CCO.N=O.C#N");
		StructureBuilder.makeHydrogensExplicit(state);
		String smiles = new SMILESWriter(f).generateSmiles();
		assertEquals("CCO.N=O.C#N", smiles);
	}
	
	@Test
	public void testOrganic1() throws StructureBuildingException {
		Fragment f = state.fragManager.buildSMILES("[S]");
		String smiles = new SMILESWriter(f).generateSmiles();
		assertEquals("[S]", smiles);
	}
	
	@Test
	public void testOrganic2() throws StructureBuildingException {
		Fragment f = state.fragManager.buildSMILES("[S][H]");
		String smiles = new SMILESWriter(f).generateSmiles();
		assertEquals("[SH]", smiles);
	}
	
	@Test
	public void testOrganic3() throws StructureBuildingException {
		Fragment f = state.fragManager.buildSMILES("[S]([H])[H]");
		String smiles = new SMILESWriter(f).generateSmiles();
		assertEquals("S", smiles);
	}
	
	@Test
	public void testOrganic4() throws StructureBuildingException {
		Fragment f = state.fragManager.buildSMILES("[S]([H])([H])[H]");
		String smiles = new SMILESWriter(f).generateSmiles();
		assertEquals("[SH3]", smiles);
	}
	
	@Test
	public void testOrganic5() throws StructureBuildingException {
		Fragment f = state.fragManager.buildSMILES("[S]([H])([H])([H])[H]");
		String smiles = new SMILESWriter(f).generateSmiles();
		assertEquals("[SH4]", smiles);
	}
	
	@Test
	public void testOrganic6() throws StructureBuildingException {
		Fragment f = state.fragManager.buildSMILES("S(F)(F)(F)F");
		String smiles = new SMILESWriter(f).generateSmiles();
		assertEquals("S(F)(F)(F)F", smiles);
	}
	
	@Test
	public void testOrganic7() throws StructureBuildingException {
		Fragment f = state.fragManager.buildSMILES("S([H])(F)(F)(F)(F)F");
		String smiles = new SMILESWriter(f).generateSmiles();
		assertEquals("S(F)(F)(F)(F)F", smiles);
	}
	
	@Test
	public void testOrganic8() throws StructureBuildingException {
		Fragment f = state.fragManager.buildSMILES("S([H])([H])(F)(F)(F)F");
		String smiles = new SMILESWriter(f).generateSmiles();
		assertEquals("[SH2](F)(F)(F)F", smiles);
	}
	
	@Test
	public void testOrganic9() throws StructureBuildingException {
		Fragment f = state.fragManager.buildSMILES("S(F)(F)(F)(F)(F)(F)F");
		String smiles = new SMILESWriter(f).generateSmiles();
		assertEquals("S(F)(F)(F)(F)(F)(F)F", smiles);
	}
	
	@Test
	public void testOrganic10() throws StructureBuildingException {
		Fragment f = state.fragManager.buildSMILES("[I]([H])([H])[H]");
		String smiles = new SMILESWriter(f).generateSmiles();
		assertEquals("[IH3]", smiles);
	}
	
	@Test
	public void testCharged1() throws StructureBuildingException {
		Fragment f = state.fragManager.buildSMILES("[CH3+]");
		StructureBuilder.makeHydrogensExplicit(state);
		String smiles = new SMILESWriter(f).generateSmiles();
		assertEquals("[CH3+]", smiles);
	}
	
	@Test
	public void testCharged2() throws StructureBuildingException {
		Fragment f = state.fragManager.buildSMILES("[Mg+2]");
		StructureBuilder.makeHydrogensExplicit(state);
		String smiles = new SMILESWriter(f).generateSmiles();
		assertEquals("[Mg+2]", smiles);
	}
	@Test
	public void testCharged3() throws StructureBuildingException {
		Fragment f = state.fragManager.buildSMILES("[BH4-]");
		StructureBuilder.makeHydrogensExplicit(state);
		String smiles = new SMILESWriter(f).generateSmiles();
		assertEquals("[BH4-]", smiles);
	}
	
	@Test
	public void testCharged4() throws StructureBuildingException {
		Fragment f = state.fragManager.buildSMILES("[O-2]");
		StructureBuilder.makeHydrogensExplicit(state);
		String smiles = new SMILESWriter(f).generateSmiles();
		assertEquals("[O-2]", smiles);
	}
	
	@Test
	public void testIsotope() throws StructureBuildingException {
		Fragment f = state.fragManager.buildSMILES("[15NH3]");
		StructureBuilder.makeHydrogensExplicit(state);
		String smiles = new SMILESWriter(f).generateSmiles();
		assertEquals("[15NH3]", smiles);
	}
	
	@Test
	public void testRGroup() throws StructureBuildingException {
		Fragment f = state.fragManager.buildSMILES("[R]CC[R]");
		StructureBuilder.makeHydrogensExplicit(state);
		String smiles = new SMILESWriter(f).generateSmiles();
		assertEquals("*CC*", smiles);
	}
	
	@Test
	public void testRingOpeningsGreaterThan10() throws StructureBuildingException {
		Fragment f = state.fragManager.buildSMILES("C12=C3C4=C5C6=C1C7=C8C9=C1C%10=C%11C(=C29)C3=C2C3=C4C4=C5C5=C9C6=C7C6=C7C8=C1C1=C8C%10=C%10C%11=C2C2=C3C3=C4C4=C5C5=C%11C%12=C(C6=C95)C7=C1C1=C%12C5=C%11C4=C3C3=C5C(=C81)C%10=C23");
		StructureBuilder.makeHydrogensExplicit(state);
		String smiles = new SMILESWriter(f).generateSmiles();
		assertEquals("C12=C3C4=C5C6=C1C1=C7C8=C9C%10=C%11C(=C28)C3=C3C2=C4C4=C5C5=C8C6=C1C1=C6C7=C9C9=C7C%10=C%10C%11=C3C3=C2C2=C4C4=C5C5=C%11C%12=C(C1=C85)C6=C9C9=C%12C%12=C%11C4=C2C2=C%12C(=C79)C%10=C32", smiles);
	}
	
	@Test
	public void testHydrogenNotBondedToAnyNonHydrogen1() throws StructureBuildingException {
		Fragment f = state.fragManager.buildSMILES("[H-].[H+]");
		String smiles = new SMILESWriter(f).generateSmiles();
		assertEquals("[H-].[H+]", smiles);
	}
	
	@Test
	public void testHydrogenNotBondedToAnyNonHydrogen2() throws StructureBuildingException {
		Fragment f = state.fragManager.buildSMILES("[H][H]");
		String smiles = new SMILESWriter(f).generateSmiles();
		assertEquals("[H][H]", smiles);
	}
	
	@Test
	public void testHydrogenNotBondedToAnyNonHydrogen3() throws StructureBuildingException {
		Fragment f = state.fragManager.buildSMILES("[2H][H]");
		String smiles = new SMILESWriter(f).generateSmiles();
		assertEquals("[2H][H]", smiles);
	}
	
	@Test
	public void testHydrogenNotBondedToAnyNonHydrogen4() throws StructureBuildingException {
		Fragment f = state.fragManager.buildSMILES("[H]B1[H]B([H])[H]1");
		String smiles = new SMILESWriter(f).generateSmiles();
		assertEquals("B1[H]B[H]1", smiles);
	}
	
	@Test
	public void testTetrahedralChirality1() throws StructureBuildingException {
		Fragment f = state.fragManager.buildSMILES("N[C@@H](F)C");
		StructureBuilder.makeHydrogensExplicit(state);
		String smiles = new SMILESWriter(f).generateSmiles();
		assertEquals("N[C@@H](F)C", smiles);
	}

	@Test
	public void testTetrahedralChirality2() throws StructureBuildingException {
		Fragment f = state.fragManager.buildSMILES("N[C@H](F)C");
		StructureBuilder.makeHydrogensExplicit(state);
		String smiles = new SMILESWriter(f).generateSmiles();
		assertEquals("N[C@H](F)C", smiles);
	}

	@Test
	public void testTetrahedralChirality3() throws StructureBuildingException {
		Fragment f = state.fragManager.buildSMILES("C2.N1.F3.[C@@H]231");
		StructureBuilder.makeHydrogensExplicit(state);
		String smiles = new SMILESWriter(f).generateSmiles();
		assertEquals("C[C@H](F)N", smiles);
	}

	@Test
	public void testTetrahedralChirality4() throws StructureBuildingException {
		Fragment f = state.fragManager.buildSMILES("[C@@H]231.C2.N1.F3");
		StructureBuilder.makeHydrogensExplicit(state);
		String smiles = new SMILESWriter(f).generateSmiles();
		assertEquals("[C@H](C)(N)F", smiles);
	}

	@Test
	public void testTetrahedralChirality5() throws StructureBuildingException {
		Fragment f = state.fragManager.buildSMILES("[C@@H](Cl)1[C@H](C)(F).Br1");
		StructureBuilder.makeHydrogensExplicit(state);
		String smiles = new SMILESWriter(f).generateSmiles();
		assertEquals("[C@H](Cl)([C@H](C)F)Br", smiles);
	}

	@Test
	public void testTetrahedralChirality6() throws StructureBuildingException {
		Fragment f = state.fragManager.buildSMILES("I[C@@](Cl)(Br)F");
		StructureBuilder.makeHydrogensExplicit(state);
		String smiles = new SMILESWriter(f).generateSmiles();
		assertEquals("I[C@@](Cl)(Br)F", smiles);
	}
	

	@Test
	public void testTetrahedralChirality7() throws StructureBuildingException {
		Fragment f = state.fragManager.buildSMILES("C[S@](N)=O");
		StructureBuilder.makeHydrogensExplicit(state);
		String smiles = new SMILESWriter(f).generateSmiles();
		assertEquals("C[S@](N)=O", smiles);
	}
	
	@Test
	public void testDoubleBondSupport1() throws StructureBuildingException {
		Fragment f = state.fragManager.buildSMILES("C/C=C/C");
		StructureBuilder.makeHydrogensExplicit(state);
		String smiles = new SMILESWriter(f).generateSmiles();
		if (!smiles.equals("C/C=C/C") && !smiles.equals("C\\C=C\\C")){
			fail(smiles +" did not correspond to one of the expected SMILES strings");
		}
	}
	
	
	@Test
	public void testDoubleBondSupport2() throws StructureBuildingException {
		Fragment f = state.fragManager.buildSMILES("C/C=C\\C");
		StructureBuilder.makeHydrogensExplicit(state);
		String smiles = new SMILESWriter(f).generateSmiles();
		if (!smiles.equals("C/C=C\\C") && !smiles.equals("C\\C=C/C")){
			fail(smiles +" did not correspond to one of the expected SMILES strings");
		}
	}
	
	
	@Test
	public void testDoubleBondSupport3() throws StructureBuildingException {
		Fragment f = state.fragManager.buildSMILES("C/C=C\\C=C/C");
		StructureBuilder.makeHydrogensExplicit(state);
		String smiles = new SMILESWriter(f).generateSmiles();
		if (!smiles.equals("C/C=C\\C=C/C") && !smiles.equals("C\\C=C/C=C\\C")){
			fail(smiles +" did not correspond to one of the expected SMILES strings");
		}
	}

	@Test
	public void testDoubleBondSupport4() throws StructureBuildingException {
		Fragment f = state.fragManager.buildSMILES("ClC(C(=O)[O-])=CC(=CC(=O)[O-])Cl");
		StructureBuilder.makeHydrogensExplicit(state);
		f.findBond(2, 6).setBondStereoElement(new Atom[]{f.getAtomByID(1), f.getAtomByID(2), f.getAtomByID(6), f.getAtomByID(7)}, BondStereoValue.TRANS);
		f.findBond(7, 8).setBondStereoElement(new Atom[]{f.getAtomByID(12), f.getAtomByID(7), f.getAtomByID(8), f.getAtomByID(9)}, BondStereoValue.TRANS);
		String smiles = new SMILESWriter(f).generateSmiles();
		if (!smiles.equals("Cl\\C(\\C(=O)[O-])=C\\C(=C/C(=O)[O-])\\Cl") && !smiles.equals("Cl/C(/C(=O)[O-])=C/C(=C\\C(=O)[O-])/Cl")){
			fail(smiles +" did not correspond to one of the expected SMILES strings");
		}
	}
	
	@Test
	public void testDoubleBondSupport5() throws StructureBuildingException {
		Fragment f = state.fragManager.buildSMILES("C/C=N\\O");
		StructureBuilder.makeHydrogensExplicit(state);
		String smiles = new SMILESWriter(f).generateSmiles();
		if (!smiles.equals("C/C=N\\O") && !smiles.equals("C\\C=N/O")){
			fail(smiles +" did not correspond to one of the expected SMILES strings");
		}
	}
	
	@Test
	public void testDoubleBondSupport6() throws StructureBuildingException {
		Fragment f = state.fragManager.buildSMILES("O=C(/C=C(C(O)=O)\\C=C/C(O)=O)O");
		StructureBuilder.makeHydrogensExplicit(state);
		String smiles = new SMILESWriter(f).generateSmiles();
		if (!smiles.equals("O=C(/C=C(/C(O)=O)\\C=C/C(O)=O)O") && !smiles.equals("O=C(\\C=C(\\C(O)=O)/C=C\\C(O)=O)O")){
			fail(smiles +" did not correspond to one of the expected SMILES strings");
		}
	}
	
	@Test
	public void testDoubleBondSupport7() throws StructureBuildingException {
		Fragment f = state.fragManager.buildSMILES("C(=C(C=CC(=O)O)C(=O)O)C(=O)O");
		StructureBuilder.makeHydrogensExplicit(state);
		f.findBond(1, 2).setBondStereoElement(new Atom[]{f.getAtomByID(11), f.getAtomByID(1), f.getAtomByID(2), f.getAtomByID(8)}, BondStereoValue.TRANS);
		f.findBond(3, 4).setBondStereoElement(new Atom[]{f.getAtomByID(2), f.getAtomByID(3), f.getAtomByID(4), f.getAtomByID(5)}, BondStereoValue.TRANS);
		String smiles = new SMILESWriter(f).generateSmiles();
		if (!smiles.equals("C(=C(/C=C/C(=O)O)\\C(=O)O)/C(=O)O") && !smiles.equals("C(=C(\\C=C\\C(=O)O)/C(=O)O)\\C(=O)O")){
			fail(smiles +" did not correspond to one of the expected SMILES strings");
		}
	}
	
	@Test
	public void testDoubleBondSupport8() throws StructureBuildingException {
		//hydrogen on the nitrogen must be explicit!
		Fragment f = state.fragManager.buildSMILES("[H]/N=C(\\N)/O");
		StructureBuilder.makeHydrogensExplicit(state);
		String smiles = new SMILESWriter(f).generateSmiles();
		if (!smiles.equals("[H]/N=C(\\N)/O") && !smiles.equals("[H]\\N=C(/N)\\O")){
			fail(smiles +" did not correspond to one of the expected SMILES strings");
		}
	}
}
