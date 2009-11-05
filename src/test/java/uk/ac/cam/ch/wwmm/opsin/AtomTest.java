package uk.ac.cam.ch.wwmm.opsin;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertFalse;
import static junit.framework.Assert.assertNotNull;
import static junit.framework.Assert.assertTrue;
import nu.xom.Builder;
import nu.xom.Element;

import org.junit.Before;
import org.junit.Test;


public class AtomTest {

	static Builder builder = new Builder();
	private Fragment frag;
	
	@Before
	public void setUp() throws Exception {
		frag = new Fragment();
	}
	
	@Test
	public void testAtom() throws StructureBuildingException {
		Atom atom = new Atom(10, "1", "C", frag);
		assertNotNull("Got atom", atom);
		assertEquals("Id = 10", 10, atom.getID());
		assertEquals("Element = C", "C", atom.getElement());
	}

	@Test
	public void testAtom2() throws Exception {
		Element cmlAtom = builder.build("<atom id=\"a10\" elementType=\"C\" formalCharge=\"-1\">" +
				"<label value=\"1\" /></atom>", "/localhost").getRootElement();
		Atom atom = new Atom(2, cmlAtom, frag);
		assertNotNull("Got atom", atom);
		assertEquals("Id = 2", 2, atom.getID());
		assertEquals("Element = C", "C", atom.getElement());
		assertEquals("Charge = -1", -1, atom.getCharge());
	}
	
	@Test
	public void testToCMLAtom() throws StructureBuildingException {
		Atom atom = new Atom(10, "1", "C", frag);
		Element elem = atom.toCMLAtom();
		assertNotNull("Got XOM Element", elem);
		assertEquals("Correct XML", "<atom id=\"a10\" elementType=\"C\"><label value=\"1\" dictRef=\"cmlDict:locant\" /></atom>", elem.toXML()); 
	}
	
	@Test
	public void testAddLocantHasLocant() throws StructureBuildingException {
		Atom atom = new Atom(10, "1", "C", frag);
		atom.addLocant("1");
		assertTrue("Atom has locant '1'", atom.hasLocant("1"));
		assertFalse("Atom has no locant 'C'", atom.hasLocant("C"));
		atom.addLocant("C");
		assertTrue("Atom now has locant 'C'", atom.hasLocant("C"));
	}
	
	@Test
	public void testGetIncomingValency() throws StructureBuildingException {
		SMILESFragmentBuilder sBuilder = new SMILESFragmentBuilder();
		IDManager idm = new IDManager();
		assertEquals("No bonds", 0, 
				sBuilder.build("C", idm).getAtomList().get(0).getIncomingValency());
		assertEquals("One bond", 1, 
				sBuilder.build("CC", idm).getAtomList().get(0).getIncomingValency());
		assertEquals("Two bonds", 2, 
				sBuilder.build("C(C)C", idm).getAtomList().get(0).getIncomingValency());
		assertEquals("Double bond", 2, 
				sBuilder.build("C=O", idm).getAtomList().get(0).getIncomingValency());
		assertEquals("Triple bond", 3, 
				sBuilder.build("C#C", idm).getAtomList().get(0).getIncomingValency());
		assertEquals("One bond", 1, 
				sBuilder.build("CC=CC#N", idm).getAtomList().get(0).getIncomingValency());
	}
	
}
