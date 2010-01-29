package uk.ac.cam.ch.wwmm.opsin;

import static junit.framework.Assert.assertEquals;

import java.util.List;

import org.junit.BeforeClass;
import org.junit.Test;

import sea36.chem.rings.Ring;
import sea36.chem.rings.SSSRFinder;
import uk.ac.cam.ch.wwmm.opsin.Fragment;
import uk.ac.cam.ch.wwmm.opsin.NameToStructure;
import uk.ac.cam.ch.wwmm.opsin.OpsinToChemKitWrapper;

public class SSSRTest {
	private static NameToStructure n2s;
	@BeforeClass
	public static void setup() throws Exception {
		n2s = NameToStructure.getInstance();
	}
	
	@Test
	public void testFindSSSR() {
		Fragment f = n2s.parseChemicalName("ovalene", false).getStructure();
		OpsinToChemKitWrapper chemKitWrapper  =  new OpsinToChemKitWrapper(f);
		List<Ring> rings = SSSRFinder.findSSSR(chemKitWrapper.getChemKitMolecule());
		assertEquals(10, rings.size());
;	}
}
