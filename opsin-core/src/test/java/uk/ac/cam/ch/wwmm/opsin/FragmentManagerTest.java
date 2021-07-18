package uk.ac.cam.ch.wwmm.opsin;


import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;

import java.io.IOException;
import java.util.ArrayList;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;


public class FragmentManagerTest {

	FragmentManager fragManager;

	@BeforeEach
	public void setUp() throws IOException{
		IDManager idManager = new IDManager();
		fragManager = new FragmentManager(new SMILESFragmentBuilder(idManager), idManager);
	}

	@Test
	public void testGetUnifiedFrags() throws StructureBuildingException {
		Fragment frag1 = fragManager.buildSMILES("CC");
		Fragment frag2 = fragManager.buildSMILES("CNC");

		fragManager.createBond(frag1.getFirstAtom(), frag2.getFirstAtom(), 1);
		Fragment frag = fragManager.getUnifiedFragment();
		assertEquals(5, frag.getAtomCount(), "Frag has five atoms");
		assertEquals(4, frag.getBondSet().size(), "Frag has four bonds");
	}

	@Test
	public void testRelabelFusedRingSystem() throws StructureBuildingException {
		Fragment naphthalene = fragManager.buildSMILES("C1=CC=CC2=CC=CC=C12");
		FragmentTools.relabelLocantsAsFusedRingSystem(naphthalene.getAtomList());
		assertEquals(1, naphthalene.getIDFromLocant("1"), "Locant 1 = atom 1");
		assertEquals(5, naphthalene.getIDFromLocant("4a"), "Locant 4a = atom 5");
		assertEquals(9, naphthalene.getIDFromLocant("8"), "Locant 8 = atom 9");
		assertEquals(10, naphthalene.getIDFromLocant("8a"), "Locant 8a = atom 10");
		assertEquals(0, naphthalene.getIDFromLocant("9"), "No locant 9");
	}
	
	@Test
	public void testCloneFragment() throws StructureBuildingException {
		Fragment urea = fragManager.buildSMILES("NC(=O)N");
		FragmentTools.assignElementLocants(urea, new ArrayList<Fragment>());
		assertNotNull(urea.getAtomByLocant("N"));
		assertNotNull(urea.getAtomByLocant("N'"));
		assertNull(urea.getAtomByLocant("N''"));
		assertNull(urea.getAtomByLocant("N'''"));
		
		Fragment primedCopy = fragManager.copyAndRelabelFragment(urea, 1);
		assertEquals(4, primedCopy.getAtomCount());
		assertNull(primedCopy.getAtomByLocant("N"));
		assertNull(primedCopy.getAtomByLocant("N'"));
		assertNotNull(primedCopy.getAtomByLocant("N''"));
		assertNotNull(primedCopy.getAtomByLocant("N'''"));
	}
}
