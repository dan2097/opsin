package uk.ac.cam.ch.wwmm.opsin;

import static junit.framework.Assert.assertEquals;

import org.junit.Before;
import org.junit.Test;


public class FragmentManagerTest {

	FragmentManager fragManager;

	@Before
	public void setUp(){
		fragManager = new FragmentManager(new SMILESFragmentBuilder(), new CMLFragmentBuilder(new ResourceGetter("uk/ac/cam/ch/wwmm/opsin/resources/")), new IDManager());
	}

	@Test
	public void testGetUnifiedFrags() throws StructureBuildingException {
		Fragment frag1 = fragManager.buildSMILES("CC");
		Fragment frag2 = fragManager.buildSMILES("CNC");

		fragManager.createBond(frag1.getAtomByLocant("1"), frag2.getAtomByLocant("1"), 1);
		Fragment frag = fragManager.getUnifiedFragment();
		assertEquals("Frag has five atoms", 5, frag.getAtomList().size());
		assertEquals("Frag has four bonds", 4, frag.getBondSet().size());
	}

	@Test
	public void testRelabelFusedRingSystem() throws StructureBuildingException {
		SMILESFragmentBuilder sBuilder = new SMILESFragmentBuilder();
		Fragment naphthalene = sBuilder.build("C1=CC=CC2=CC=CC=C12", fragManager);
		FragmentTools.relabelFusedRingSystem(naphthalene);
		assertEquals("Locant 1 = atom 1", 1, naphthalene.getIDFromLocant("1"));
		assertEquals("Locant 4a = atom 5", 5, naphthalene.getIDFromLocant("4a"));
		assertEquals("Locant 8 = atom 9", 9, naphthalene.getIDFromLocant("8"));
		assertEquals("Locant 8a = atom 10", 10, naphthalene.getIDFromLocant("8a"));
		assertEquals("No locant 9", 0, naphthalene.getIDFromLocant(""));
	}

}
