package uk.ac.cam.ch.wwmm.opsin;


import static org.junit.jupiter.api.Assertions.assertEquals;

import java.util.List;

import org.junit.jupiter.api.Test;

public class SSSRTest {
	
	@Test
	public void testFindSSSR() throws Exception {
		NameToStructure n2s = NameToStructure.getInstance();
		Fragment f = n2s.parseChemicalName("violanthrene").getStructure();
		List<Ring> rings = SSSRFinder.getSetOfSmallestRings(f);
		assertEquals(9, rings.size());
		
		f = n2s.parseChemicalName("aceanthrene").getStructure();
		rings = SSSRFinder.getSetOfSmallestRings(f);
		assertEquals(4, rings.size());
	}
}
