package uk.ac.cam.ch.wwmm.opsin;
import static junit.framework.Assert.*;
import static org.mockito.Mockito.mock;

import org.junit.Before;
import org.junit.Test;

public class HeteroAtomReplacementTest {

	FragmentManager fragManager;
	Atom a;

	@Before
	public void setUp() throws Exception {
		fragManager = new FragmentManager(new SMILESFragmentBuilder(), mock(CMLFragmentBuilder.class), new IDManager());
		a = new Atom(0, "C", mock(Fragment.class));
	}
	
	@Test
	public void thia() throws StructureBuildingException{
		fragManager.makeHeteroatom(a, "S", false);
		assertEquals(0, a.getCharge());
		assertEquals(0, a.getProtonsExplicitlyAddedOrRemoved());
		assertEquals(2, StructureBuildingMethods.calculateSubstitutableHydrogenAtoms(a));
	}

	@Test
	public void thionia() throws StructureBuildingException{
		fragManager.makeHeteroatom(a, "[SH3+]", false);
		assertEquals(1, a.getCharge());
		assertEquals(1, a.getProtonsExplicitlyAddedOrRemoved());
		assertEquals(3, StructureBuildingMethods.calculateSubstitutableHydrogenAtoms(a));
	}
	
	@Test
	public void sulfanylia() throws StructureBuildingException{
		fragManager.makeHeteroatom(a, "[SH+]", false);
		assertEquals(1, a.getCharge());
		assertEquals(-1, a.getProtonsExplicitlyAddedOrRemoved());
		assertEquals(1, StructureBuildingMethods.calculateSubstitutableHydrogenAtoms(a));
	}
	
	@Test
	public void sulfanida() throws StructureBuildingException{
		fragManager.makeHeteroatom(a, "[SH-]", false);
		assertEquals(-1, a.getCharge());
		assertEquals(-1, a.getProtonsExplicitlyAddedOrRemoved());
		assertEquals(1, StructureBuildingMethods.calculateSubstitutableHydrogenAtoms(a));
	}
	
	@Test
	public void sulfanuida() throws StructureBuildingException{
		fragManager.makeHeteroatom(a, "[SH3-]", false);
		assertEquals(-1, a.getCharge());
		assertEquals(1, a.getProtonsExplicitlyAddedOrRemoved());
		assertEquals(3, StructureBuildingMethods.calculateSubstitutableHydrogenAtoms(a));
	}
	
	@Test
	public void replaceNeutralWithCharged() throws StructureBuildingException{
		Atom a = new Atom(0, "C", mock(Fragment.class));
		fragManager.makeHeteroatom(a, "[NH4+]", false);
		assertEquals(1, a.getCharge());
		assertEquals(1, a.getProtonsExplicitlyAddedOrRemoved());
		assertEquals(4, StructureBuildingMethods.calculateSubstitutableHydrogenAtoms(a));
	}
	
	@Test
	public void replaceChargedWithEquallyCharged() throws StructureBuildingException{
		Atom a = new Atom(0, "C", mock(Fragment.class));
		a.addChargeAndProtons(1, -1);
		fragManager.makeHeteroatom(a, "[NH4+]", false);
		assertEquals(1, a.getCharge());
		assertEquals(1, a.getProtonsExplicitlyAddedOrRemoved());
		assertEquals(4, StructureBuildingMethods.calculateSubstitutableHydrogenAtoms(a));
	}
	
    @Test(expected=StructureBuildingException.class)
	public void replaceChargedWithUnEquallyCharged() throws StructureBuildingException{
		Atom a = new Atom(0, "C", mock(Fragment.class));
		a.addChargeAndProtons(1, -1);
		fragManager.makeHeteroatom(a, "[NH2-]", false);
	}
}
