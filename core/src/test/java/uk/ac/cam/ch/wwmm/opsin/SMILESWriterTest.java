package uk.ac.cam.ch.wwmm.opsin;
import static org.mockito.Mockito.mock;
import static junit.framework.Assert.*;

import java.util.Collections;
import java.util.List;

import org.junit.Before;
import org.junit.Test;

public class SMILESWriterTest {
	
	BuildState state;
	@Before
	public void setup(){
		state = new BuildState(mock(NameToStructureConfig.class), new SMILESFragmentBuilder(), mock(CMLFragmentBuilder.class));
	}


	@Test
	public void testRoundTrip1() throws StructureBuildingException {
		Fragment f = state.fragManager.buildSMILES("[NH4+].[Cl-].F.[He-2]");
		StructureBuilder.makeHydrogensExplicit(state);
		String smiles = new SMILESWriter(f).generateSmiles();
		assertEquals("[NH4+].[Cl-].F.[He-2]", smiles);
	}
	
	@Test
	public void testRoundTrip2() throws StructureBuildingException {
		Fragment f = state.fragManager.buildSMILES("[NH4+].[Cl-].F.[He-2]");
		List<Atom> atomList = f.getAtomList();
		Collections.reverse(atomList);
		f.reorderAtomCollection(atomList);
		StructureBuilder.makeHydrogensExplicit(state);
		String smiles = new SMILESWriter(f).generateSmiles();
		assertEquals("[He-2].F.[Cl-].[NH4+]", smiles);
	}
	
	@Test
	public void testRoundTrip3() throws StructureBuildingException {
		Fragment f = state.fragManager.buildSMILES("CCO.N=O.C#N");
		StructureBuilder.makeHydrogensExplicit(state);
		String smiles = new SMILESWriter(f).generateSmiles();
		assertEquals("CCO.N=O.C#N", smiles);
	}
	
	@Test
	public void testRoundTrip4() throws StructureBuildingException {
		Fragment f = state.fragManager.buildSMILES("C1=CC=CC=C1");
		StructureBuilder.makeHydrogensExplicit(state);
		String smiles = new SMILESWriter(f).generateSmiles();
		assertEquals("C1=CC=CC=C1", smiles);
	}
	
	@Test
	public void testRoundTrip5() throws StructureBuildingException {
		Fragment f = state.fragManager.buildSMILES(StringTools.multiplyString("C",200));
		StructureBuilder.makeHydrogensExplicit(state);
		String smiles = new SMILESWriter(f).generateSmiles();
		assertEquals(StringTools.multiplyString("C",200), smiles);
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
	public void testDoubleBondSupport() throws StructureBuildingException {
		Fragment f = state.fragManager.buildSMILES("C/C=C\\C=C/C");
		StructureBuilder.makeHydrogensExplicit(state);
		String smiles = new SMILESWriter(f).generateSmiles();
		if (!smiles.equals("C/C=C\\C=C/C") && !smiles.equals("C\\C=C/C=C\\C")){
			fail(smiles +" did not correspond to one of the expected SMILES strings");
		}
	}
	
	
}
