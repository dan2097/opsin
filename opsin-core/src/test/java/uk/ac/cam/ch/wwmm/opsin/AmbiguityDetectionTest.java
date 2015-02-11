package uk.ac.cam.ch.wwmm.opsin;

import static org.junit.Assert.assertEquals;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

import org.apache.commons.io.IOUtils;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import uk.ac.cam.ch.wwmm.opsin.OpsinResult.OPSIN_RESULT_STATUS;

public class AmbiguityDetectionTest {
	
	private static NameToStructure n2s;

	@BeforeClass
	public static void setUp() {
		n2s = NameToStructure.getInstance();
	}
	
	@AfterClass
	public static void cleanUp(){
		n2s = null;
	}

	@Test
	public void testNamesThatShouldBeDetectedAsAmbiguous() throws IOException{
		checkNames("ambiguous.txt", OPSIN_RESULT_STATUS.WARNING);
	}
	
	@Test
	public void testUnAmbiguousCounterExamples() throws IOException{
		checkNames("unambiguous.txt", OPSIN_RESULT_STATUS.SUCCESS);
	}
	
	private void checkNames(String file, OPSIN_RESULT_STATUS expectedStatus) throws IOException{
		BufferedReader input = new BufferedReader(new InputStreamReader(this.getClass().getResourceAsStream(file)));
		try {
			String line = null;
			while ((line = input.readLine()) != null) {
				if(line.startsWith("//")){
					continue;
				}
				OpsinResult result = n2s.parseChemicalName(line);
				assertEquals(line + " gave unexpected result", expectedStatus, result.getStatus());
			}
		} finally {
			IOUtils.closeQuietly(input);
		}
	}
	
}
