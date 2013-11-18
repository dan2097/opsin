package uk.ac.cam.ch.wwmm.opsin;
import static junit.framework.Assert.*;
import static org.mockito.Mockito.mock;
import org.junit.Before;
import org.junit.Test;

public class HeteroAtomReplacementTest {

	FragmentManager fragManager;
	Atom a;

	@Before
	public void setUp() {
		IDManager idManager = new IDManager();
		fragManager = new FragmentManager(new SMILESFragmentBuilder(idManager), idManager);
		a = new Atom(0, "C", mock(Fragment.class));
	}
	
	@Test
	public void thia() throws StructureBuildingException{
		fragManager.replaceAtomWithSmiles(a, "S");
		assertEquals(0, a.getCharge());
		assertEquals(0, a.getProtonsExplicitlyAddedOrRemoved());
		assertEquals(2, StructureBuildingMethods.calculateSubstitutableHydrogenAtoms(a));
	}

	@Test
	public void thionia() throws StructureBuildingException{
		fragManager.replaceAtomWithSmiles(a, "[SH3+]");
		assertEquals(1, a.getCharge());
		assertEquals(1, a.getProtonsExplicitlyAddedOrRemoved());
		assertEquals(3, StructureBuildingMethods.calculateSubstitutableHydrogenAtoms(a));
	}
	
	@Test
	public void sulfanylia() throws StructureBuildingException{
		fragManager.replaceAtomWithSmiles(a, "[SH+]");
		assertEquals(1, a.getCharge());
		assertEquals(-1, a.getProtonsExplicitlyAddedOrRemoved());
		assertEquals(1, StructureBuildingMethods.calculateSubstitutableHydrogenAtoms(a));
	}
	
	@Test
	public void sulfanida() throws StructureBuildingException{
		fragManager.replaceAtomWithSmiles(a, "[SH-]");
		assertEquals(-1, a.getCharge());
		assertEquals(-1, a.getProtonsExplicitlyAddedOrRemoved());
		assertEquals(1, StructureBuildingMethods.calculateSubstitutableHydrogenAtoms(a));
	}
	
	@Test
	public void sulfanuida() throws StructureBuildingException{
		fragManager.replaceAtomWithSmiles(a, "[SH3-]");
		assertEquals(-1, a.getCharge());
		assertEquals(1, a.getProtonsExplicitlyAddedOrRemoved());
		assertEquals(3, StructureBuildingMethods.calculateSubstitutableHydrogenAtoms(a));
	}
	
	@Test
	public void replaceNeutralWithCharged() throws StructureBuildingException{
		Atom a = new Atom(0, "C", mock(Fragment.class));
		fragManager.replaceAtomWithSmiles(a, "[NH4+]");
		assertEquals(1, a.getCharge());
		assertEquals(1, a.getProtonsExplicitlyAddedOrRemoved());
		assertEquals(4, StructureBuildingMethods.calculateSubstitutableHydrogenAtoms(a));
	}
	
	@Test
	public void replaceChargedWithEquallyCharged() throws StructureBuildingException{
		Atom a = new Atom(0, "C", mock(Fragment.class));
		a.addChargeAndProtons(1, -1);
		fragManager.replaceAtomWithSmiles(a, "[NH4+]");
		assertEquals(1, a.getCharge());
		assertEquals(1, a.getProtonsExplicitlyAddedOrRemoved());
		assertEquals(4, StructureBuildingMethods.calculateSubstitutableHydrogenAtoms(a));
	}
	
    @Test(expected=StructureBuildingException.class)
	public void replaceChargedWithUnEquallyCharged() throws StructureBuildingException{
		Atom a = new Atom(0, "C", mock(Fragment.class));
		a.addChargeAndProtons(1, -1);
		fragManager.replaceAtomWithSmiles(a, "[NH2-]");
	}
}
