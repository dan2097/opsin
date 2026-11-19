package uk.ac.cam.ch.wwmm.opsin;

import uk.ac.cam.ch.wwmm.opsin.OpsinWarning.OpsinWarningType;
import uk.ac.cam.ch.wwmm.opsin.OpsinResult.OPSIN_RESULT_STATUS;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.mockito.Mockito.mock;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

public class HeteroAtomReplacementTest {

	FragmentManager fragManager;
	Atom a;

	@BeforeEach
	public void setUp() {
		IDManager idManager = new IDManager();
		fragManager = new FragmentManager(new SMILESFragmentBuilder(idManager), idManager);
		a = new Atom(0, ChemEl.C, mock(Fragment.class));
	}
	
	@Test
	public void thia() throws StructureBuildingException{
		fragManager.replaceAtomWithSmiles(a, "S");
		assertEquals(0, a.getCharge());
		assertEquals(0, a.getProtonsExplicitlyAddedOrRemoved());
		assertEquals(2, StructureBuildingMethods.calculateSubstitutableHydrogenAtoms(a));
	}

	@Test
	public void thionia() throws StructureBuildingException{
		fragManager.replaceAtomWithSmiles(a, "[SH3+]");
		assertEquals(1, a.getCharge());
		assertEquals(1, a.getProtonsExplicitlyAddedOrRemoved());
		assertEquals(3, StructureBuildingMethods.calculateSubstitutableHydrogenAtoms(a));
	}
	
	@Test
	public void sulfanylia() throws StructureBuildingException{
		fragManager.replaceAtomWithSmiles(a, "[SH+]");
		assertEquals(1, a.getCharge());
		assertEquals(-1, a.getProtonsExplicitlyAddedOrRemoved());
		assertEquals(1, StructureBuildingMethods.calculateSubstitutableHydrogenAtoms(a));
	}
	
	@Test
	public void sulfanida() throws StructureBuildingException{
		fragManager.replaceAtomWithSmiles(a, "[SH-]");
		assertEquals(-1, a.getCharge());
		assertEquals(-1, a.getProtonsExplicitlyAddedOrRemoved());
		assertEquals(1, StructureBuildingMethods.calculateSubstitutableHydrogenAtoms(a));
	}
	
	@Test
	public void sulfanuida() throws StructureBuildingException{
		fragManager.replaceAtomWithSmiles(a, "[SH3-]");
		assertEquals(-1, a.getCharge());
		assertEquals(1, a.getProtonsExplicitlyAddedOrRemoved());
		assertEquals(3, StructureBuildingMethods.calculateSubstitutableHydrogenAtoms(a));
	}
	
	@Test
	public void replaceNeutralWithCharged() throws StructureBuildingException{
		Atom a = new Atom(0, ChemEl.C, mock(Fragment.class));
		fragManager.replaceAtomWithSmiles(a, "[NH4+]");
		assertEquals(1, a.getCharge());
		assertEquals(1, a.getProtonsExplicitlyAddedOrRemoved());
		assertEquals(4, StructureBuildingMethods.calculateSubstitutableHydrogenAtoms(a));
	}
	
	@Test
	public void replaceChargedWithEquallyCharged() throws StructureBuildingException{
		Atom a = new Atom(0, ChemEl.C, mock(Fragment.class));
		a.addChargeAndProtons(1, -1);
		fragManager.replaceAtomWithSmiles(a, "[NH4+]");
		assertEquals(1, a.getCharge());
		assertEquals(1, a.getProtonsExplicitlyAddedOrRemoved());
		assertEquals(4, StructureBuildingMethods.calculateSubstitutableHydrogenAtoms(a));
	}
	
	@Test()
	public void replaceChargedWithUnEquallyCharged() {
		assertThrows(StructureBuildingException.class, () -> {
			Atom a = new Atom(0, ChemEl.C, mock(Fragment.class));
			a.addChargeAndProtons(1, -1);
			fragManager.replaceAtomWithSmiles(a, "[NH2-]");
		});
	}

	@Test
	public void unlocantedReplacementGuessPrefersSaturatedAtom() {
		//sila-linalool is the silanol methyl(4-methylpent-3-enyl)vinylsilanol, not the Si=C compound
		OpsinResult result = NameToStructure.getInstance().parseChemicalName("sila-linalool");
		assertEquals(OPSIN_RESULT_STATUS.WARNING, result.getStatus());
		assertEquals(OpsinWarningType.APPEARS_AMBIGUOUS, result.getWarnings().get(0).getType());
		Atom silicon = null;
		for (Atom atom : result.getStructure().getAtomList()) {
			if (atom.getElement() == ChemEl.Si) {
				silicon = atom;
			}
		}
		assertNotNull(silicon);
		for (Bond bond : silicon.getBonds()) {
			assertEquals(1, bond.getOrder());
		}
		boolean silanol = false;
		for (Atom neighbour : silicon.getAtomNeighbours()) {
			if (neighbour.getElement() == ChemEl.O) {
				silanol = true;
			}
		}
		assertTrue(silanol);
	}

	@Test
	public void unambiguousReplacementUnaffected() {
		OpsinResult result = NameToStructure.getInstance().parseChemicalName("silabenzene");
		assertEquals(OPSIN_RESULT_STATUS.SUCCESS, result.getStatus());
		assertEquals("[SiH]1=CC=CC=C1", result.getSmiles());
	}
}
