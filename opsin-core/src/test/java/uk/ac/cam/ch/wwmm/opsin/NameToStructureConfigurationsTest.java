package uk.ac.cam.ch.wwmm.opsin;


import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import uk.ac.cam.ch.wwmm.opsin.OpsinResult.OPSIN_RESULT_STATUS;

public class NameToStructureConfigurationsTest {

		private static NameToStructure n2s;

		@BeforeAll
		public static void setUp() {
			n2s = NameToStructure.getInstance();
		}
		
		@AfterAll
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
			assertEquals("*C", or.getSmiles());
		}

		@Test
		public void testOutputRadicalsAsAttachments() {
			NameToStructureConfig n2sConfig = NameToStructureConfig.getDefaultConfigInstance();
			n2sConfig.setAllowRadicals(true);
			n2sConfig.setOutputRadicalsAsWildCardAtoms(false);
			OpsinResult or = n2s.parseChemicalName("methyl", n2sConfig);
			assertEquals("[CH3]", or.getSmiles(SmilesOptions.CXSMILES_ATOM_LABELS));
			
			n2sConfig.setOutputRadicalsAsWildCardAtoms(true);
			or = n2s.parseChemicalName("methyl", n2sConfig);
			assertEquals("*C |$_AP1$|", or.getSmiles(SmilesOptions.CXSMILES_ATOM_LABELS));
		}
		
		@Test
		public void testInterpretAcidsWithoutTheWordAcid() throws StructureBuildingException {
			NameToStructureConfig n2sConfig = NameToStructureConfig.getDefaultConfigInstance();
			n2sConfig.setInterpretAcidsWithoutTheWordAcid(false);
			OpsinResult or = n2s.parseChemicalName("acetic", n2sConfig);
			assertEquals(OPSIN_RESULT_STATUS.FAILURE, or.getStatus());

			n2sConfig.setInterpretAcidsWithoutTheWordAcid(true);
			or = n2s.parseChemicalName("acetic", n2sConfig);
			assertEquals(OPSIN_RESULT_STATUS.SUCCESS, or.getStatus());
		}

		@Test
		public void testWarnRatherThanFailOnUninterpretableStereochemistry() throws StructureBuildingException {
			NameToStructureConfig n2sConfig = NameToStructureConfig.getDefaultConfigInstance();
			n2sConfig.setWarnRatherThanFailOnUninterpretableStereochemistry(false);
			OpsinResult or = n2s.parseChemicalName("(R)-2,2'-Bis(diphenylphosphino)-1,1'-binaphthyl", n2sConfig);
			assertEquals(OPSIN_RESULT_STATUS.FAILURE, or.getStatus());

			n2sConfig.setWarnRatherThanFailOnUninterpretableStereochemistry(true);
			or = n2s.parseChemicalName("(R)-2,2'-Bis(diphenylphosphino)-1,1'-binaphthyl", n2sConfig);
			assertEquals(OPSIN_RESULT_STATUS.WARNING, or.getStatus());
		}
}
