package uk.ac.cam.ch.wwmm.opsin;

import static junit.framework.Assert.assertNotNull;
import static org.junit.Assert.assertEquals;
import static org.mockito.Mockito.mock;

import java.util.List;
import java.util.Set;

import junit.framework.Assert;

import nu.xom.Element;

import org.junit.Before;
import org.junit.Test;

public class SMILESFragmentBuilderTest {

	private Fragment fragment;
	private SMILESFragmentBuilder sBuilder;
	private FragmentManager fm;
	
	@Before
	public void setUp() throws Exception {
		sBuilder = new SMILESFragmentBuilder();
		fm = new FragmentManager(sBuilder, mock(CMLFragmentBuilder.class), new IDManager());
	}

	@Test
	public void testBuild() throws StructureBuildingException {		
		fragment = sBuilder.build("C", fm);
		assertNotNull("Got a fragment", fragment);
	}
	
	@Test
	public void testSimple1() throws StructureBuildingException {		
		fragment = sBuilder.build("CC", fm);
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(2, atomList.size());
		assertEquals("C", atomList.get(0).getElement());
		assertEquals("C", atomList.get(1).getElement());
	}
	
	@Test
	public void testSimple2() throws StructureBuildingException {		
		fragment = sBuilder.build("O=C=O", fm);
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(3, atomList.size());
		assertEquals("O", atomList.get(0).getElement());
		assertEquals("C", atomList.get(1).getElement());
		assertEquals("O", atomList.get(2).getElement());
		Set<Bond> bonds = fragment.getBondSet();
		assertEquals(2, bonds.size());
		for (Bond bond : bonds) {
			assertEquals(2, bond.getOrder());
		}
	}
	
	@Test
	public void testSimple3() throws StructureBuildingException {		
		fragment = sBuilder.build("C#N", fm);
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(2, atomList.size());
		Set<Bond> bonds = fragment.getBondSet();
		assertEquals(1, bonds.size());
		for (Bond bond : bonds) {
			assertEquals(3, bond.getOrder());
		}
	}
	
	@Test
	public void testSimple4() throws StructureBuildingException {		
		fragment = sBuilder.build("CCN(CC)CC", fm);
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(7, atomList.size());
		Atom nitrogen = atomList.get(2);
		assertEquals("N", nitrogen.getElement());
		assertEquals(3, nitrogen.getBonds().size());
		List<Atom> neighbours = nitrogen.getAtomNeighbours();//bonds and hence neighbours come from a linked hash set so the order of the neighbours is deterministic
		assertEquals(3, neighbours.size());
		assertEquals(atomList.get(1), neighbours.get(0));
		assertEquals(atomList.get(3), neighbours.get(1));
		assertEquals(atomList.get(5), neighbours.get(2));
	}
	
	@Test
	public void testSimple5() throws StructureBuildingException {		
		fragment = sBuilder.build("CC(=O)O", fm);
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(4, atomList.size());
		Atom carbon = atomList.get(1);
		List<Atom> neighbours = carbon.getAtomNeighbours();
		assertEquals(3, neighbours.size());
		assertEquals(atomList.get(0), neighbours.get(0));
		assertEquals(atomList.get(2), neighbours.get(1));
		assertEquals(atomList.get(3), neighbours.get(2));
		assertEquals(2, fragment.findBondOrThrow(carbon, atomList.get(2)).getOrder());
	}
	
	@Test
	public void testSimple6() throws StructureBuildingException {		
		fragment = sBuilder.build("C1CCCCC1", fm);
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(6, atomList.size());
		for (Atom atom : atomList) {
			assertEquals(2, atom.getAtomNeighbours().size());
			assertEquals(false, atom.hasSpareValency());
		}
	}
	
	@Test
	public void testSimple7() throws StructureBuildingException {		
		fragment = sBuilder.build("c1ccccc1", fm);
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(6, atomList.size());
		for (Atom atom : atomList) {
			assertEquals(2, atom.getAtomNeighbours().size());
			assertEquals(true, atom.hasSpareValency());
		}
	}
	
	
	@Test
	public void testSimple8() throws StructureBuildingException {		
		fragment = sBuilder.build("[I-].[Na+]", fm);
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(2, atomList.size());
		Atom iodine = atomList.get(0);
		assertEquals(0, iodine.getAtomNeighbours().size());
		assertEquals(-1, iodine.getCharge());
		
		Atom sodium = atomList.get(1);
		assertEquals(0, sodium.getAtomNeighbours().size());
		assertEquals(1, sodium.getCharge());
	}
	
	@Test
	public void testSimple9() throws StructureBuildingException {		
		fragment = sBuilder.build("(C(=O)O)", fm);
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(3, atomList.size());
		Atom carbon = atomList.get(0);
		assertEquals(2, carbon.getAtomNeighbours().size());
	}
	
	@Test
	public void testSimple10() throws StructureBuildingException {		
		fragment = sBuilder.build("C-C-O", fm);
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(3, atomList.size());
	}
	
	@Test
	public void testSimple11() throws StructureBuildingException {		
		fragment = sBuilder.build("NC(Cl)(Br)C(=O)O", fm);
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(7, atomList.size());
		assertEquals("Cl", atomList.get(2).getElement());
	}
	

	@Test
	public void unterminatedRingOpening() {
		try {
			fragment = sBuilder.build("C1CC", fm);
			Assert.fail("Should throw exception for bad smiles");
		} catch (Exception e) {
			;
		}
	}

	@Test
	public void doublePositiveCharge1() throws StructureBuildingException {
		fragment = sBuilder.build("[C++]", fm);
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(1, atomList.size());
		assertEquals(2, atomList.get(0).getCharge());
	}
	
	@Test
	public void doublePositiveCharge2() throws StructureBuildingException {
		fragment = sBuilder.build("[C+2]", fm);
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(1, atomList.size());
		assertEquals(2, atomList.get(0).getCharge());
	}
	
	@Test
	public void doubleNegativeCharge1() throws StructureBuildingException {
		fragment = sBuilder.build("[O--]", fm);
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(1, atomList.size());
		assertEquals(-2, atomList.get(0).getCharge());
	}
	
	@Test
	public void doubleNegativeCharge2() throws StructureBuildingException {
		fragment = sBuilder.build("[O-2]", fm);
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(1, atomList.size());
		assertEquals(-2, atomList.get(0).getCharge());
	}

	@Test
	public void badlyFormedSMILE1() {
		try {
			fragment = sBuilder.build("H5", fm);
			Assert.fail("Should throw exception for bad smiles");
		} catch (Exception e) {
			;
		}
	}

	@Test
	public void badlyFormedSMILE2() {
		try {
			fragment = sBuilder.build("CH4", fm);
			Assert.fail("Should throw exception for bad smiles");
		} catch (Exception e) {
			;
		}
	}

	@Test
	public void badlyFormedSMILE3() {
		try {
			fragment = sBuilder.build("13C", fm);
			Assert.fail("Should throw exception for bad smiles");
		} catch (Exception e) {
			;
		}
	}

	@Test
	public void ringClosureHandling1() throws StructureBuildingException {
		fragment = sBuilder.build("C=1CN1", fm);
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(3, atomList.size());
		assertEquals(2, fragment.findBond(atomList.get(0), atomList.get(2)).getOrder());
	}
	
	@Test
	public void ringClosureHandling2() throws StructureBuildingException {
		fragment = sBuilder.build("C1CN=1", fm);
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(3, atomList.size());
		assertEquals(2, fragment.findBond(atomList.get(0), atomList.get(2)).getOrder());
	}
	
	@Test
	public void ringClosureHandling3() throws StructureBuildingException {
		try {
			fragment = sBuilder.build("C#1CN=1", fm);
			Assert.fail("Should throw exception for bad smiles");
		} catch (Exception e) {
			;
		}
	}
	
	@Test
	public void ringClosureHandling4() throws StructureBuildingException {
		fragment = sBuilder.build("C=1CN=1", fm);
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(3, atomList.size());
		assertEquals(2, fragment.findBond(atomList.get(0), atomList.get(2)).getOrder());
	}
	
	@Test
	public void ringSupportGreaterThan10() throws StructureBuildingException {
		fragment = sBuilder.build("C%10CC%10", fm);
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(3, atomList.size());
		assertEquals(2, atomList.get(0).getAtomNeighbours().size());
	}

	@Test
	public void chiralityTest1() throws StructureBuildingException {
		fragment = sBuilder.build("N[C@@H](F)C", fm);
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(4, atomList.size());
		Atom chiralAtom = atomList.get(1);
		assertEquals(3, chiralAtom.getAtomNeighbours().size());
		Element atomParity = chiralAtom.getAtomParityElement();
		String atomRefs4 = atomParity.getAttributeValue(XmlDeclarations.ATOMREFS4_ATR);
		assertEquals("a1 a2_H a3 a4", atomRefs4);
		assertEquals(1, Integer.parseInt(atomParity.getValue()));
	}
	
	@Test
	public void chiralityTest2() throws StructureBuildingException {
		fragment = sBuilder.build("N[C@H](F)C", fm);
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(4, atomList.size());
		Atom chiralAtom = atomList.get(1);
		assertEquals(3, chiralAtom.getAtomNeighbours().size());
		Element atomParity = chiralAtom.getAtomParityElement();
		String atomRefs4 = atomParity.getAttributeValue(XmlDeclarations.ATOMREFS4_ATR);
		assertEquals("a1 a2_H a3 a4", atomRefs4);
		assertEquals(-1, Integer.parseInt(atomParity.getValue()));
	}
	
	@Test
	public void chiralityTest3() throws StructureBuildingException {
		fragment = sBuilder.build("C2.N1.F3.[C@@H]231", fm);
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(4, atomList.size());
		Atom chiralAtom = atomList.get(3);
		assertEquals(3, chiralAtom.getAtomNeighbours().size());
		Element atomParity = chiralAtom.getAtomParityElement();
		String atomRefs4 = atomParity.getAttributeValue(XmlDeclarations.ATOMREFS4_ATR);
		assertEquals("a4_H a1 a3 a2", atomRefs4);
		assertEquals(1, Integer.parseInt(atomParity.getValue()));
	}
	
	@Test
	public void chiralityTest4() throws StructureBuildingException {
		fragment = sBuilder.build("[C@@H]231.C2.N1.F3", fm);
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(4, atomList.size());
		Atom chiralAtom = atomList.get(0);
		assertEquals(3, chiralAtom.getAtomNeighbours().size());
		Element atomParity = chiralAtom.getAtomParityElement();
		String atomRefs4 = atomParity.getAttributeValue(XmlDeclarations.ATOMREFS4_ATR);
		assertEquals("a1_H a2 a4 a3", atomRefs4);
		assertEquals(1, Integer.parseInt(atomParity.getValue()));
	}
	
	@Test
	public void chiralityTest5() throws StructureBuildingException {
		fragment = sBuilder.build("[C@@H](Cl)1[C@H](C)(F).Br1", fm);
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(6, atomList.size());
		Atom chiralAtom1 = atomList.get(0);
		assertEquals(3, chiralAtom1.getAtomNeighbours().size());
		Element atomParity = chiralAtom1.getAtomParityElement();
		String atomRefs4 = atomParity.getAttributeValue(XmlDeclarations.ATOMREFS4_ATR);
		assertEquals("a1_H a2 a6 a3", atomRefs4);
		assertEquals(1, Integer.parseInt(atomParity.getValue()));
		
		Atom chiralAtom2 = atomList.get(2);
		assertEquals(3, chiralAtom2.getAtomNeighbours().size());
		atomParity = chiralAtom2.getAtomParityElement();
		atomRefs4 = atomParity.getAttributeValue(XmlDeclarations.ATOMREFS4_ATR);
		assertEquals("a1 a3_H a4 a5", atomRefs4);
		assertEquals(-1, Integer.parseInt(atomParity.getValue()));
	}
	
	@Test
	public void chiralityTest6() throws StructureBuildingException {
		fragment = sBuilder.build("I[C@@](Cl)(Br)F", fm);
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(5, atomList.size());
		Atom chiralAtom = atomList.get(1);
		assertEquals(4, chiralAtom.getAtomNeighbours().size());
		Element atomParity = chiralAtom.getAtomParityElement();
		String atomRefs4 = atomParity.getAttributeValue(XmlDeclarations.ATOMREFS4_ATR);
		assertEquals("a1 a3 a4 a5", atomRefs4);
		assertEquals(1, Integer.parseInt(atomParity.getValue()));
	}
	
	@Test
    public void testDoubleBondStereo1() throws StructureBuildingException {
        fragment = sBuilder.build("F/C=C/F", fm);
        Bond b =fragment.findBond(2, 3);
        assertEquals("T", b.getBondStereoElement().getValue());
    }

    @Test
    public void testDoubleBondStereo2() throws StructureBuildingException {
        fragment = sBuilder.build("F\\C=C/F", fm);
        Bond b =fragment.findBond(2, 3);
        assertEquals("C", b.getBondStereoElement().getValue());
    }

    @Test
    public void testDoubleBondStereo3() throws StructureBuildingException {
        fragment = sBuilder.build("C(/F)=C/F", fm);
        Bond b =fragment.findBond(1, 3);
        assertEquals("C", b.getBondStereoElement().getValue());
    }

    @Test
    public void testDoubleBondStereo4() throws StructureBuildingException {
        fragment = sBuilder.build("C(\\F)=C/F", fm);
        Bond b =fragment.findBond(1, 3);
        assertEquals("T", b.getBondStereoElement().getValue());
    }

    @Test
    public void testDoubleBondStereo5a() throws StructureBuildingException {
        fragment = sBuilder.build("CC1=C/F.O\\1", fm);
        Bond b =fragment.findBond(2, 3);
        assertEquals("C", b.getBondStereoElement().getValue());
    }

    @Test
    public void testDoubleBondStereo5b() throws StructureBuildingException {
        fragment = sBuilder.build("CC/1=C/F.O1", fm);
        Bond b =fragment.findBond(2, 3);
        assertEquals("C", b.getBondStereoElement().getValue());
    }

    @Test
    public void testDoubleBondStereo6() throws StructureBuildingException {
        fragment = sBuilder.build("CC1=C/F.O/1", fm);
        Bond b =fragment.findBond(2, 3);
        assertEquals("T", b.getBondStereoElement().getValue());
    }

    @Test
    public void testDoubleBondMulitStereo1() throws StructureBuildingException {
        fragment = sBuilder.build("F/C=C/C=C/C", fm);
        Bond b =fragment.findBond(2, 3);
        assertEquals("T", b.getBondStereoElement().getValue());
        b =fragment.findBond(4, 5);
        assertEquals("T", b.getBondStereoElement().getValue());
    }

    @Test
    public void testDoubleBondMulitStereo2() throws StructureBuildingException {
        fragment = sBuilder.build("F/C=C\\C=C/C", fm);
        Bond b =fragment.findBond(2, 3);
        assertEquals("C", b.getBondStereoElement().getValue());
        b =fragment.findBond(4, 5);
        assertEquals("C", b.getBondStereoElement().getValue());
    }

    @Test
    public void testDoubleBondMulitStereo3() throws StructureBuildingException {
        fragment = sBuilder.build("F/C=C\\C=C\\C", fm);
        Bond b =fragment.findBond(2, 3);
        assertEquals("C", b.getBondStereoElement().getValue());
        b =fragment.findBond(4, 5);
        assertEquals("T", b.getBondStereoElement().getValue());
    }

    @Test
    public void testDoubleBondMulitStereo4() throws StructureBuildingException {
        fragment = sBuilder.build("F/C=C\\C=CC", fm);
        Bond b =fragment.findBond(2, 3);
        assertEquals("C", b.getBondStereoElement().getValue());
        b =fragment.findBond(4, 5);
        assertEquals(null, b.getBondStereoElement());
    } 
}
