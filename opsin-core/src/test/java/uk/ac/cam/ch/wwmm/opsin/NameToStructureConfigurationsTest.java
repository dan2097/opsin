package uk.ac.cam.ch.wwmm.opsin;
import static org.junit.Assert.*;


import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import uk.ac.cam.ch.wwmm.opsin.OpsinResult.OPSIN_RESULT_STATUS;

public class NameToStructureConfigurationsTest {

		private static NameToStructure n2s;

		@BeforeClass
		public static void setUp() {
			n2s = NameToStructure.getInstance();
		}
		
		@AfterClass
		public static void cleanUp() {
			n2s = null;
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
		
		@Test
		public void testOutputRadicalsAsWildCards() throws StructureBuildingException {
			NameToStructureConfig n2sConfig = NameToStructureConfig.getDefaultConfigInstance();
			n2sConfig.setAllowRadicals(true);
			n2sConfig.setOutputRadicalsAsWildCardAtoms(false);
			OpsinResult or = n2s.parseChemicalName("methyl", n2sConfig);
			assertEquals("[CH3]", or.getSmiles());

			n2sConfig.setOutputRadicalsAsWildCardAtoms(true);
			or = n2s.parseChemicalName("methyl", n2sConfig);
			assertEquals("C*", or.getSmiles());
		}
}
