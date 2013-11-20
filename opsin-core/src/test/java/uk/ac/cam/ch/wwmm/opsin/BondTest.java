package uk.ac.cam.ch.wwmm.opsin;

import static org.junit.Assert.*;

import org.junit.Test;

import nu.xom.Element;


public class BondTest {
	
	@Test
	public void testBond() {
		Fragment frag = new Fragment();
		Atom a1 = new Atom(1, "C", frag);
		Atom a2 = new Atom(2, "C", frag);
		frag.addAtom(a1);
		frag.addAtom(a2);
		Bond bond = new Bond(a1, a2, 1);
		assertNotNull("Got bond", bond);
		assertEquals("From = 1", 1, bond.getFrom());
		assertEquals("To = 2", 2, bond.getTo());
		assertEquals("Order = 1", 1, bond.getOrder());
	}
	
	@Test
	public void testToCMLBond() {
		Fragment frag = new Fragment();
		Atom a1 = new Atom(1, "C", frag);
		Atom a2 = new Atom(2, "C", frag);
		frag.addAtom(a1);
		frag.addAtom(a2);
		Bond bond = new Bond(a1, a2, 1);
		Element elem = bond.toCMLBond();
		assertNotNull("Got XOM Element", elem);
		assertEquals("Correct XML", "<bond xmlns=\"http://www.xml-cml.org/schema\" id=\"a1_a2\" atomRefs2=\"a1 a2\" order=\"S\" />", elem.toXML());
	}
	
}
