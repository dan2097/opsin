package uk.ac.cam.ch.wwmm.opsin;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;

import org.junit.jupiter.api.Test;

public class PolymerTest {

	@Test
	public void testSimplePolymer() throws ParsingException {
		OpsinResult result = NameToStructure.getInstance().parseChemicalName("poly(oxyethylene)");
		String smiles = result.getSmiles();
		assertNotNull(smiles);
		assertEquals(true, smiles.contains("[*:1]"));
		assertEquals(true, smiles.contains("[*:2]"));
		
		String cml = result.getCml();
		assertEquals(true, cml.contains("alpha"));
		assertEquals(true, cml.contains("omega"));
	}
}
