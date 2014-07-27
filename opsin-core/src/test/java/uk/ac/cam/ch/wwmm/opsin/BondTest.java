package uk.ac.cam.ch.wwmm.opsin;

import static org.junit.Assert.*;
import static org.mockito.Mockito.mock;

import org.junit.Test;

import uk.ac.cam.ch.wwmm.opsin.Bond.SMILES_BOND_DIRECTION;
import uk.ac.cam.ch.wwmm.opsin.BondStereo.BondStereoValue;

public class BondTest {
	
	@Test
	public void testBond() {
		Fragment frag = new Fragment(mock(Element.class));
		Atom a1 = new Atom(1, ChemEl.C, frag);
		Atom a2 = new Atom(2, ChemEl.C, frag);
		frag.addAtom(a1);
		frag.addAtom(a2);
		Bond bond = new Bond(a1, a2, 1);
		assertNotNull("Got bond", bond);
		assertEquals("From = 1", 1, bond.getFrom());
		assertEquals("To = 2", 2, bond.getTo());
		assertEquals("Order = 1", 1, bond.getOrder());
		assertEquals(a1, bond.getFromAtom());
		assertEquals(a2, bond.getToAtom());
		assertEquals(a2, bond.getOtherAtom(a1));
		assertEquals(a1, bond.getOtherAtom(a2));
		assertEquals(null, bond.getBondStereo());
		assertEquals(null, bond.getSmilesStereochemistry());
	}

	@Test
	public void testBondMutation() {
		Fragment frag = new Fragment(mock(Element.class));
		Atom a1 = new Atom(1, ChemEl.C, frag);
		Atom a2 = new Atom(2, ChemEl.C, frag);
		Atom a3 = new Atom(3, ChemEl.C, frag);
		Atom a4 = new Atom(4, ChemEl.C, frag);
		frag.addAtom(a1);
		frag.addAtom(a2);
		frag.addAtom(a3);
		frag.addAtom(a4);
		Bond bond = new Bond(a2, a3, 1);
		bond.setOrder(2);
		assertEquals("Order = 2", 2, bond.getOrder());
		BondStereo bondStereo = new BondStereo(new Atom[]{a1,a2,a3,a4}, BondStereoValue.TRANS);
		bond.setBondStereo(bondStereo);
		assertEquals(bondStereo, bond.getBondStereo());
		bond.setSmilesStereochemistry(SMILES_BOND_DIRECTION.LSLASH);
		assertEquals(SMILES_BOND_DIRECTION.LSLASH, bond.getSmilesStereochemistry());
	}
}
