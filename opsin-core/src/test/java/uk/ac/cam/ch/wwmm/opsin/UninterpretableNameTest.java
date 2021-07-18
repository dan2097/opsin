package uk.ac.cam.ch.wwmm.opsin;

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvFileSource;

import uk.ac.cam.ch.wwmm.opsin.OpsinResult.OPSIN_RESULT_STATUS;

public class UninterpretableNameTest {
	
	private static NameToStructure n2s;

	@BeforeAll
	public static void setUp() {
		n2s = NameToStructure.getInstance();
	}
	
	@AfterAll
	public static void cleanUp(){
		n2s = null;
	}

	@ParameterizedTest
	@CsvFileSource(resources ="uninterpretable.txt",  delimiter='\t')
	public void testNamesThatShoudlBeUninterpretable(String uninterpretablName) {
		OpsinResult result = n2s.parseChemicalName(uninterpretablName);
		assertEquals(OPSIN_RESULT_STATUS.FAILURE, result.getStatus(), uninterpretablName + " should be uninterpretable");
	}
	
}
