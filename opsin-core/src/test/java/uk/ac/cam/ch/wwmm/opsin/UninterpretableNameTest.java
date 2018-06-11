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

public class UninterpretableNameTest {
	
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
	public void testNamesThatShoudlBeUninterpretable() throws IOException{
		checkNames("uninterpretable.txt");
	}

	private void checkNames(String file) throws IOException{
		BufferedReader input = new BufferedReader(new InputStreamReader(this.getClass().getResourceAsStream(file)));
		try {
			String line = null;
			while ((line = input.readLine()) != null) {
				if(line.startsWith("//")){
					continue;
				}
				OpsinResult result = n2s.parseChemicalName(line);
				assertEquals(line + " gave unexpected result", OPSIN_RESULT_STATUS.FAILURE, result.getStatus());
			}
		} finally {
			IOUtils.closeQuietly(input);
		}
	}
	
}
