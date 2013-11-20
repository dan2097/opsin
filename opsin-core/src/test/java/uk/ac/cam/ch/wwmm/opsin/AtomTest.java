package uk.ac.cam.ch.wwmm.opsin;

import static org.junit.Assert.*;

import org.junit.Before;
import org.junit.Test;

import nu.xom.Element;

public class AtomTest {

	private Fragment frag;
	private SMILESFragmentBuilder sBuilder = new SMILESFragmentBuilder(new IDManager());
	
	@Before
	public void setUp() {
		frag = new Fragment();
	}
	
	@Test
	public void testAtom() {
		Atom atom = new Atom(10, "C", frag);
		assertNotNull("Got atom", atom);
		assertEquals("Id = 10", 10, atom.getID());
		assertEquals("Element = C", "C", atom.getElement());
	}
	
	@Test
	public void testToCMLAtom() {
		Atom atom = new Atom(10, "C", frag);
		atom.addLocant("1");
		Element elem = atom.toCMLAtom();
		assertNotNull("Got XOM Element", elem);
		assertEquals("Correct XML", "<atom xmlns=\"http://www.xml-cml.org/schema\" id=\"a10\" elementType=\"C\" hydrogenCount=\"0\"><label value=\"1\" dictRef=\"cmlDict:locant\" /></atom>", elem.toXML()); 
	}
	
	@Test
	public void testAddLocantHasLocant() {
		Atom atom = new Atom(10, "C", frag);
		atom.addLocant("1");
		assertTrue("Atom has locant '1'", atom.hasLocant("1"));
		assertFalse("Atom has no locant 'C'", atom.hasLocant("C"));
		atom.addLocant("C");
		assertTrue("Atom now has locant 'C'", atom.hasLocant("C"));
	}
	
	@Test
	public void testGetIncomingValency() throws StructureBuildingException {
		assertEquals("No bonds", 0, sBuilder.build("C").getFirstAtom().getIncomingValency());
		assertEquals("One bond", 1, sBuilder.build("CC").getFirstAtom().getIncomingValency());
		assertEquals("Two bonds", 2, sBuilder.build("C(C)C").getFirstAtom().getIncomingValency());
		assertEquals("Double bond", 2, sBuilder.build("C=O").getFirstAtom().getIncomingValency());
		assertEquals("Triple bond", 3, sBuilder.build("C#C").getFirstAtom().getIncomingValency());
		assertEquals("One bond", 1, sBuilder.build("CC=CC#N").getFirstAtom().getIncomingValency());
	}
	
}
