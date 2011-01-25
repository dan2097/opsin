package uk.ac.cam.ch.wwmm.opsin;

import static junit.framework.Assert.assertEquals;

import java.util.List;

import org.junit.Test;

import uk.ac.cam.ch.wwmm.opsin.Fragment;
import uk.ac.cam.ch.wwmm.opsin.NameToStructure;

public class SSSRTest {
	
	@Test
	public void testFindSSSR() throws Exception {
		NameToStructure n2s = NameToStructure.getInstance();
		Fragment f = n2s.parseChemicalName("violanthrene", false).getStructure();
		List<Ring> rings = SSSRFinder.getSetOfSmallestRings(f);
		assertEquals(9, rings.size());
		
		f = n2s.parseChemicalName("aceanthrene", false).getStructure();
		rings = SSSRFinder.getSetOfSmallestRings(f);
		assertEquals(4, rings.size());
	}
}
