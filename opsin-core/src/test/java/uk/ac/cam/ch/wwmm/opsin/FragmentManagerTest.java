package uk.ac.cam.ch.wwmm.opsin;

import static org.junit.Assert.*;

import java.io.IOException;
import java.util.ArrayList;

import org.junit.Before;
import org.junit.Test;

public class FragmentManagerTest {

	FragmentManager fragManager;

	@Before
	public void setUp() throws IOException{
		IDManager idManager = new IDManager();
		fragManager = new FragmentManager(new SMILESFragmentBuilder(idManager), idManager);
	}

	@Test
	public void testGetUnifiedFrags() throws StructureBuildingException {
		Fragment frag1 = fragManager.buildSMILES("CC");
		Fragment frag2 = fragManager.buildSMILES("CNC");

		fragManager.createBond(frag1.getAtomByLocant("1"), frag2.getAtomByLocant("1"), 1);
		Fragment frag = fragManager.getUnifiedFragment();
		assertEquals("Frag has five atoms", 5, frag.getAtomCount());
		assertEquals("Frag has four bonds", 4, frag.getBondSet().size());
	}

	@Test
	public void testRelabelFusedRingSystem() throws StructureBuildingException {
		Fragment naphthalene = fragManager.buildSMILES("C1=CC=CC2=CC=CC=C12");
		FragmentTools.relabelLocantsAsFusedRingSystem(naphthalene.getAtomList());
		assertEquals("Locant 1 = atom 1", 1, naphthalene.getIDFromLocant("1"));
		assertEquals("Locant 4a = atom 5", 5, naphthalene.getIDFromLocant("4a"));
		assertEquals("Locant 8 = atom 9", 9, naphthalene.getIDFromLocant("8"));
		assertEquals("Locant 8a = atom 10", 10, naphthalene.getIDFromLocant("8a"));
		assertEquals("No locant 9", 0, naphthalene.getIDFromLocant(""));
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
