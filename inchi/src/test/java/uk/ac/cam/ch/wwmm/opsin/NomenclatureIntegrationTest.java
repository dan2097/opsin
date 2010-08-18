package uk.ac.cam.ch.wwmm.opsin;
import static org.junit.Assert.*;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;


import org.junit.Before;
import org.junit.Test;

public class NomenclatureIntegrationTest {
	NameToStructure n2s;
	@Before
	public void setUp() throws NameToStructureException {
		n2s = NameToStructure.getInstance();
	}
	
	@Test
	public void testRadicals() throws Exception {
		NameToStructureConfig n2sConfig = NameToStructureConfig.getDefaultConfigInstance();
		n2sConfig.setAllowRadicals(true);
		String file = "radicals.txt";
		checkNamesAgainstInChIs(file, n2sConfig);
	}
	
	@Test
	public void testEpoxyLike() throws Exception {
		NameToStructureConfig n2sConfig = NameToStructureConfig.getDefaultConfigInstance();
		n2sConfig.setAllowRadicals(true);
		String file = "epoxyLike.txt";
		checkNamesAgainstInChIs(file, n2sConfig);
	}

	private void checkNamesAgainstInChIs(String file, NameToStructureConfig n2sConfig) throws IOException{
		BufferedReader input = new BufferedReader(new InputStreamReader(this.getClass().getResourceAsStream(file)));
		try {
			String line = null;
			while ((line = input.readLine()) != null) {
				String[] lineArray = line.split("\t");
				String inchi = NameToInchi.convertResultToInChI(n2s.parseChemicalName(lineArray[0], n2sConfig), false);
				if (inchi!=null){
					String opsinInchi = InchiPruner.mainChargeAndStereochemistryLayers(inchi);
					String referenceInchi = InchiPruner.mainChargeAndStereochemistryLayers(lineArray[1]);
	
					if (!opsinInchi.equals(referenceInchi)){
						fail(lineArray[0] +" was misinterpreted as: " + inchi);
					}
				} else {
					fail(lineArray[0] +" was uninterpretable");
				}
			}
		} finally {
			input.close();
		}
	}
}
