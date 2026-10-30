package uk.ac.cam.ch.wwmm.opsin;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvFileSource;

import uk.ac.cam.ch.wwmm.opsin.OpsinResult.OPSIN_RESULT_STATUS;

public class ChalcogenAcidLocantTest {

	private static NameToStructure n2s;

	@BeforeAll
	public static void setUp() {
		n2s = NameToStructure.getInstance();
	}

	@AfterAll
	public static void cleanUp() {
		n2s = null;
	}

	@ParameterizedTest
	@CsvFileSource(resources = "chalcogenAcidLocants.txt", delimiter = '\t')
	public void testTerminalChalcogens(String name, String expectedTerminalChalcogens) throws StructureBuildingException {
		OpsinResult result = n2s.parseChemicalName(name);
		assertEquals(OPSIN_RESULT_STATUS.SUCCESS, result.getStatus(), name + " should be interpretable: " + result.getMessage());
		assertEquals(expectedTerminalChalcogens, terminalChalcogens(result.getStructure()), name);
	}

	//InChI treats the acidic hydrogens as mobile, so S-acid and O-acid are told apart by the terminal chalcogens e.g. "=O OH SH"
	private static String terminalChalcogens(Fragment frag) throws StructureBuildingException {
		List<String> terminalChalcogens = new ArrayList<>();
		for (Atom atom : frag.getAtomList()) {
			ChemEl chemEl = atom.getElement();
			if (chemEl == ChemEl.O || chemEl == ChemEl.S || chemEl == ChemEl.Se || chemEl == ChemEl.Te) {
				List<Atom> heavyNeighbours = new ArrayList<>();
				for (Atom neighbour : atom.getAtomNeighbours()) {
					if (neighbour.getElement() != ChemEl.H) {
						heavyNeighbours.add(neighbour);
					}
				}
				if (heavyNeighbours.size() == 1) {
					terminalChalcogens.add(atom.getBondToAtomOrThrow(heavyNeighbours.get(0)).getOrder() == 2 ? "=" + chemEl : chemEl + "H");
				}
			}
		}
		Collections.sort(terminalChalcogens);
		return String.join(" ", terminalChalcogens);
	}
}
