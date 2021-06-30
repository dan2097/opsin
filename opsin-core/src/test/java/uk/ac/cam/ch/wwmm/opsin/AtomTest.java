package uk.ac.cam.ch.wwmm.opsin;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.mockito.Mockito.mock;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;


public class AtomTest {

	private Fragment frag;
	private SMILESFragmentBuilder sBuilder = new SMILESFragmentBuilder(new IDManager());
	
	@BeforeEach
	public void setUp() {
		frag = new Fragment(mock(Element.class));
	}
	
	@Test
	public void testAtom() {
		Atom atom = new Atom(10, ChemEl.C, frag);
		assertNotNull(atom, "Got atom");
		assertEquals(10, atom.getID(), "Id = 10");
		assertEquals(ChemEl.C, atom.getElement(), "Element = C");
	}
	
	@Test
	public void testAddLocantHasLocant() {
		Atom atom = new Atom(10, ChemEl.C, frag);
		atom.addLocant("1");
		assertTrue(atom.hasLocant("1"), "Atom has locant '1'");
		assertFalse(atom.hasLocant("C"), "Atom has no locant 'C'");
		atom.addLocant("C");
		assertTrue(atom.hasLocant("C"), "Atom now has locant 'C'");
	}
	
	@Test
	public void testGetIncomingValency() throws StructureBuildingException {
		assertEquals(0, sBuilder.build("C").getFirstAtom().getIncomingValency(), "No bonds");
		assertEquals(1, sBuilder.build("CC").getFirstAtom().getIncomingValency(), "One bond");
		assertEquals(2, sBuilder.build("C(C)C").getFirstAtom().getIncomingValency(), "Two bonds");
		assertEquals(2, sBuilder.build("C=O").getFirstAtom().getIncomingValency(), "Double bond");
		assertEquals(3, sBuilder.build("C#C").getFirstAtom().getIncomingValency(), "Triple bond");
		assertEquals(1, sBuilder.build("CC=CC#N").getFirstAtom().getIncomingValency(), "One bond");
	}
	
}
