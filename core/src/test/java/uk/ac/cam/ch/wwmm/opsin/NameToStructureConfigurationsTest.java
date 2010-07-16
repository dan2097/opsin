package uk.ac.cam.ch.wwmm.opsin;
import static org.junit.Assert.*;


import org.junit.Before;
import org.junit.Test;

import uk.ac.cam.ch.wwmm.opsin.OpsinResult.OPSIN_RESULT_STATUS;

public class NameToStructureConfigurationsTest {

		NameToStructure n2s;
		@Before
		public void setUp() throws NameToStructureException {
			n2s = NameToStructure.getInstance();
		}
		
		@Test
		public void testAllowRadicals() throws StructureBuildingException {
			NameToStructureConfig n2sConfig = NameToStructureConfig.getDefaultConfigInstance();
			n2sConfig.setAllowRadicals(false);
			OpsinResult or = n2s.parseChemicalName("methyl", n2sConfig);
			assertEquals(OPSIN_RESULT_STATUS.FAILURE, or.getStatus());
			
			n2sConfig.setAllowRadicals(true);
			or = n2s.parseChemicalName("methyl", n2sConfig);
			assertEquals(OPSIN_RESULT_STATUS.SUCCESS, or.getStatus());
		}
}
