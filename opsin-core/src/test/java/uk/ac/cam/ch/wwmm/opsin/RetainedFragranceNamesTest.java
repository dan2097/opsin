package uk.ac.cam.ch.wwmm.opsin;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;

import uk.ac.cam.ch.wwmm.opsin.OpsinResult.OPSIN_RESULT_STATUS;
import uk.ac.cam.ch.wwmm.opsin.OpsinWarning.OpsinWarningType;

public class RetainedFragranceNamesTest {

	private static NameToStructure n2s;

	@BeforeAll
	public static void setUp() {
		n2s = NameToStructure.getInstance();
	}

	@ParameterizedTest
	@CsvSource(delimiter = ';', value = {
		"ambrox; ambroxide",
		"ambroxan; ambrofix",
		"cetalox; ambrox dl",
		"timberol; norlimbanol",
		"methyl cedryl ketone; acetyl cedrene",
		"vertofix coeur; acetylcedrene",
		"dartanol; levosandol",
		"bacdanol; sandranol",
		"trans-alpha-bergamotol; alpha-trans-bergamotol"
	})
	public void synonymsGiveTheSameStructure(String name1, String name2) {
		assertEquals(smiles(name1), smiles(name2));
	}

	@Test
	public void explicitEZOverridesStoredGeometry() {
		//natural beta-santalol is (Z); (E)-beta-santalol must not silently give (Z)
		String natural = smiles("beta-santalol");
		assertEquals(natural, smiles("(Z)-beta-santalol"));
		assertNotEquals(natural, smiles("(E)-beta-santalol"));
	}

	@Test
	public void opticalRotationIsIgnoredWithWarning() {
		OpsinResult result = n2s.parseChemicalName("(+)-nootkatone");
		assertEquals(OPSIN_RESULT_STATUS.WARNING, result.getStatus());
		assertEquals(OpsinWarningType.STEREOCHEMISTRY_IGNORED, result.getWarnings().get(0).getType());
		assertEquals(smiles("nootkatone"), result.getSmiles());
	}

	@Test
	public void methylCedrylKetoneIsTheCedrene() {
		//the trade product is acetyl cedrene (CAS 32388-55-9), which has a ring C=C conjugated to the ketone
		Fragment frag = n2s.parseChemicalName("methyl cedryl ketone").getStructure();
		int doubleBonds = 0;
		for (Bond b : frag.getBondSet()) {
			if (b.getOrder() == 2) {
				doubleBonds++;
			}
		}
		assertEquals(2, doubleBonds);
	}

	@Test
	public void isocamphylSubstituent() {
		OpsinResult result = n2s.parseChemicalName("3-isocamphylcyclohexanol");
		assertEquals(OPSIN_RESULT_STATUS.SUCCESS, result.getStatus(), result.getMessage());
		//Sandela, C16H28O
		int carbon = 0;
		int hydrogen = 0;
		int oxygen = 0;
		for (Atom a : result.getStructure().getAtomList()) {
			switch (a.getElement()) {
			case C: carbon++; break;
			case H: hydrogen++; break;
			case O: oxygen++; break;
			default: break;
			}
		}
		assertEquals(16, carbon);
		assertEquals(28, hydrogen);
		assertEquals(1, oxygen);
	}

	private static String smiles(String name) {
		OpsinResult result = n2s.parseChemicalName(name);
		assertTrue(result.getStatus() != OPSIN_RESULT_STATUS.FAILURE, name + ": " + result.getMessage());
		return result.getSmiles();
	}
}
