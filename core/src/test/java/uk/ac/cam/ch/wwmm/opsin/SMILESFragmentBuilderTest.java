package uk.ac.cam.ch.wwmm.opsin;

import static junit.framework.Assert.assertNotNull;
import static org.junit.Assert.assertEquals;
import static org.mockito.Mockito.mock;

import java.util.List;
import java.util.Set;

import junit.framework.Assert;

import org.junit.Before;
import org.junit.Test;

public class SMILESFragmentBuilderTest {

	private FragmentManager fm;

	@Before
	public void setUp() throws Exception {
		fm = new FragmentManager(new SMILESFragmentBuilder(), mock(CMLFragmentBuilder.class), new IDManager());
	}

	@Test
	public void testBuild() throws StructureBuildingException {
		Fragment fragment = fm.buildSMILES("C");
		assertNotNull("Got a fragment", fragment);
	}

	@Test
	public void testSimple1() throws StructureBuildingException {
		Fragment fragment = fm.buildSMILES("CC");
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(2, atomList.size());
		assertEquals("C", atomList.get(0).getElement());
		assertEquals("C", atomList.get(1).getElement());
	}

	@Test
	public void testSimple2() throws StructureBuildingException {
		Fragment fragment = fm.buildSMILES("O=C=O");
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
		Fragment fragment = fm.buildSMILES("C#N");
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
		Fragment fragment = fm.buildSMILES("CCN(CC)CC");
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
		Fragment fragment = fm.buildSMILES("CC(=O)O");
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
		Fragment fragment = fm.buildSMILES("C1CCCCC1");
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(6, atomList.size());
		for (Atom atom : atomList) {
			assertEquals(2, atom.getAtomNeighbours().size());
			assertEquals(false, atom.hasSpareValency());
		}
	}

	@Test
	public void testSimple7() throws StructureBuildingException {
		Fragment fragment = fm.buildSMILES("c1ccccc1");
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(6, atomList.size());
		for (Atom atom : atomList) {
			assertEquals(2, atom.getAtomNeighbours().size());
			assertEquals(true, atom.hasSpareValency());
		}
	}


	@Test
	public void testSimple8() throws StructureBuildingException {
		Fragment fragment = fm.buildSMILES("[I-].[Na+]");
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
		Fragment fragment = fm.buildSMILES("(C(=O)O)");
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(3, atomList.size());
		Atom carbon = atomList.get(0);
		assertEquals(2, carbon.getAtomNeighbours().size());
	}

	@Test
	public void testSimple10() throws StructureBuildingException {
		Fragment fragment = fm.buildSMILES("C-C-O");
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(3, atomList.size());
	}

	@Test
	public void testSimple11() throws StructureBuildingException {
		Fragment fragment = fm.buildSMILES("NC(Cl)(Br)C(=O)O");
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(7, atomList.size());
		assertEquals("Cl", atomList.get(2).getElement());
	}


	@Test(expected=StructureBuildingException.class)
	public void unterminatedRingOpening() throws StructureBuildingException {
		fm.buildSMILES("C1CC");
		Assert.fail("Should throw exception for bad smiles");
	}

	@Test
	public void doublePositiveCharge1() throws StructureBuildingException {
		Fragment fragment = fm.buildSMILES("[C++]");
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(1, atomList.size());
		assertEquals(2, atomList.get(0).getCharge());
	}

	@Test
	public void doublePositiveCharge2() throws StructureBuildingException {
		Fragment fragment = fm.buildSMILES("[C+2]");
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(1, atomList.size());
		assertEquals(2, atomList.get(0).getCharge());
	}

	@Test
	public void doubleNegativeCharge1() throws StructureBuildingException {
		Fragment fragment = fm.buildSMILES("[O--]");
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(1, atomList.size());
		assertEquals(-2, atomList.get(0).getCharge());
	}

	@Test
	public void doubleNegativeCharge2() throws StructureBuildingException {
		Fragment fragment = fm.buildSMILES("[O-2]");
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(1, atomList.size());
		assertEquals(-2, atomList.get(0).getCharge());
	}

	@Test(expected=StructureBuildingException.class)
	public void badlyFormedSMILE1() throws StructureBuildingException {
		fm.buildSMILES("H5");
		Assert.fail("Should throw exception for bad smiles");
	}

	@Test(expected=StructureBuildingException.class)
	public void badlyFormedSMILE2() throws StructureBuildingException {
		fm.buildSMILES("CH4");
		Assert.fail("Should throw exception for bad smiles");
	}

	@Test(expected=StructureBuildingException.class)
	public void badlyFormedSMILE3() throws StructureBuildingException {
		fm.buildSMILES("13C");
		Assert.fail("Should throw exception for bad smiles");
	}

	@Test
	public void ringClosureHandling1() throws StructureBuildingException {
		Fragment fragment = fm.buildSMILES("C=1CN1");
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(3, atomList.size());
		assertEquals(2, fragment.findBond(atomList.get(0), atomList.get(2)).getOrder());
	}

	@Test
	public void ringClosureHandling2() throws StructureBuildingException {
		Fragment fragment = fm.buildSMILES("C1CN=1");
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(3, atomList.size());
		assertEquals(2, fragment.findBond(atomList.get(0), atomList.get(2)).getOrder());
	}

	@Test(expected=StructureBuildingException.class)
	public void ringClosureHandling3() throws StructureBuildingException {
		fm.buildSMILES("C#1CN=1");
		Assert.fail("Should throw exception for bad smiles");
	}

	@Test
	public void ringClosureHandling4() throws StructureBuildingException {
		Fragment fragment = fm.buildSMILES("C=1CN=1");
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(3, atomList.size());
		assertEquals(2, fragment.findBond(atomList.get(0), atomList.get(2)).getOrder());
	}

	@Test
	public void ringSupportGreaterThan10() throws StructureBuildingException {
		Fragment fragment = fm.buildSMILES("C%10CC%10");
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(3, atomList.size());
		assertEquals(2, atomList.get(0).getAtomNeighbours().size());
	}

	@Test
	public void chiralityTest1() throws StructureBuildingException {
		Fragment fragment = fm.buildSMILES("N[C@@H](F)C");
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(4, atomList.size());
		Atom chiralAtom = atomList.get(1);
		assertEquals(3, chiralAtom.getAtomNeighbours().size());
		AtomParity atomParity  = chiralAtom.getAtomParity();
		Atom[] atomRefs4 = atomParity.getAtomRefs4();
		assertEquals(atomList.get(0), atomRefs4[0]);
		assertEquals(AtomParity.hydrogen, atomRefs4[1]);
		assertEquals(atomList.get(2), atomRefs4[2]);
		assertEquals(atomList.get(3), atomRefs4[3]);
		assertEquals(1, atomParity.getParity());
	}

	@Test
	public void chiralityTest2() throws StructureBuildingException {
		Fragment fragment = fm.buildSMILES("N[C@H](F)C");
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(4, atomList.size());
		Atom chiralAtom = atomList.get(1);
		assertEquals(3, chiralAtom.getAtomNeighbours().size());
		AtomParity atomParity  = chiralAtom.getAtomParity();
		Atom[] atomRefs4 = atomParity.getAtomRefs4();
		assertEquals(atomList.get(0), atomRefs4[0]);
		assertEquals(AtomParity.hydrogen, atomRefs4[1]);
		assertEquals(atomList.get(2), atomRefs4[2]);
		assertEquals(atomList.get(3), atomRefs4[3]);
		assertEquals(-1, atomParity.getParity());
	}

	@Test
	public void chiralityTest3() throws StructureBuildingException {
		Fragment fragment = fm.buildSMILES("C2.N1.F3.[C@@H]231");
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(4, atomList.size());
		Atom chiralAtom = atomList.get(3);
		assertEquals(3, chiralAtom.getAtomNeighbours().size());
		AtomParity atomParity  = chiralAtom.getAtomParity();
		Atom[] atomRefs4 = atomParity.getAtomRefs4();
		assertEquals(AtomParity.hydrogen, atomRefs4[0]);
		assertEquals(atomList.get(0), atomRefs4[1]);
		assertEquals(atomList.get(2), atomRefs4[2]);
		assertEquals(atomList.get(1), atomRefs4[3]);
		assertEquals(1, atomParity.getParity());
	}

	@Test
	public void chiralityTest4() throws StructureBuildingException {
		Fragment fragment = fm.buildSMILES("[C@@H]231.C2.N1.F3");
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(4, atomList.size());
		Atom chiralAtom = atomList.get(0);
		assertEquals(3, chiralAtom.getAtomNeighbours().size());
		AtomParity atomParity  = chiralAtom.getAtomParity();
		Atom[] atomRefs4 = atomParity.getAtomRefs4();
		assertEquals(AtomParity.hydrogen, atomRefs4[0]);
		assertEquals(atomList.get(1), atomRefs4[1]);
		assertEquals(atomList.get(3), atomRefs4[2]);
		assertEquals(atomList.get(2), atomRefs4[3]);
		assertEquals(1, atomParity.getParity());
	}

	@Test
	public void chiralityTest5() throws StructureBuildingException {
		Fragment fragment = fm.buildSMILES("[C@@H](Cl)1[C@H](C)(F).Br1");
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(6, atomList.size());
		Atom chiralAtom1 = atomList.get(0);
		assertEquals(3, chiralAtom1.getAtomNeighbours().size());
		AtomParity atomParity  = chiralAtom1.getAtomParity();
		Atom[] atomRefs4 = atomParity.getAtomRefs4();
		assertEquals(AtomParity.hydrogen, atomRefs4[0]);
		assertEquals(atomList.get(1), atomRefs4[1]);
		assertEquals(atomList.get(5), atomRefs4[2]);
		assertEquals(atomList.get(2), atomRefs4[3]);
		assertEquals(1, atomParity.getParity());

		Atom chiralAtom2 = atomList.get(2);
		assertEquals(3, chiralAtom2.getAtomNeighbours().size());
		atomParity  = chiralAtom2.getAtomParity();
		atomRefs4 = atomParity.getAtomRefs4();
		assertEquals(atomList.get(0), atomRefs4[0]);
		assertEquals(AtomParity.hydrogen, atomRefs4[1]);
		assertEquals(atomList.get(3), atomRefs4[2]);
		assertEquals(atomList.get(4), atomRefs4[3]);
		assertEquals(-1, atomParity.getParity());
	}

	@Test
	public void chiralityTest6() throws StructureBuildingException {
		Fragment fragment = fm.buildSMILES("I[C@@](Cl)(Br)F");
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(5, atomList.size());
		Atom chiralAtom = atomList.get(1);
		assertEquals(4, chiralAtom.getAtomNeighbours().size());
		AtomParity atomParity  = chiralAtom.getAtomParity();
		Atom[] atomRefs4 = atomParity.getAtomRefs4();
		assertEquals(atomList.get(0), atomRefs4[0]);
		assertEquals(atomList.get(2), atomRefs4[1]);
		assertEquals(atomList.get(3), atomRefs4[2]);
		assertEquals(atomList.get(4), atomRefs4[3]);
		assertEquals(1, atomParity.getParity());
	}

	@Test
    public void testDoubleBondStereo1() throws StructureBuildingException {
        Fragment fragment = fm.buildSMILES("F/C=C/F");
        Bond b =fragment.findBond(2, 3);
        assertEquals("T", b.getBondStereoElement().getValue());
    }

    @Test
    public void testDoubleBondStereo2() throws StructureBuildingException {
        Fragment fragment = fm.buildSMILES("F\\C=C/F");
        Bond b =fragment.findBond(2, 3);
        assertEquals("C", b.getBondStereoElement().getValue());
    }

    @Test
    public void testDoubleBondStereo3() throws StructureBuildingException {
        Fragment fragment = fm.buildSMILES("C(/F)=C/F");
        Bond b =fragment.findBond(1, 3);
        assertEquals("C", b.getBondStereoElement().getValue());
    }

    @Test
    public void testDoubleBondStereo4() throws StructureBuildingException {
        Fragment fragment = fm.buildSMILES("C(\\F)=C/F");
        Bond b =fragment.findBond(1, 3);
        assertEquals("T", b.getBondStereoElement().getValue());
    }

    @Test
    public void testDoubleBondStereo5a() throws StructureBuildingException {
        Fragment fragment = fm.buildSMILES("CC1=C/F.O\\1");
        Bond b =fragment.findBond(2, 3);
        assertEquals("C", b.getBondStereoElement().getValue());
    }

    @Test
    public void testDoubleBondStereo5b() throws StructureBuildingException {
        Fragment fragment = fm.buildSMILES("CC/1=C/F.O1");
        Bond b =fragment.findBond(2, 3);
        assertEquals("C", b.getBondStereoElement().getValue());
    }

    @Test
    public void testDoubleBondStereo6() throws StructureBuildingException {
        Fragment fragment = fm.buildSMILES("CC1=C/F.O/1");
        Bond b =fragment.findBond(2, 3);
        assertEquals("T", b.getBondStereoElement().getValue());
    }

    @Test
    public void testDoubleBondMulitStereo1() throws StructureBuildingException {
        Fragment fragment = fm.buildSMILES("F/C=C/C=C/C");
        Bond b =fragment.findBond(2, 3);
        assertEquals("T", b.getBondStereoElement().getValue());
        b =fragment.findBond(4, 5);
        assertEquals("T", b.getBondStereoElement().getValue());
    }

    @Test
    public void testDoubleBondMulitStereo2() throws StructureBuildingException {
        Fragment fragment = fm.buildSMILES("F/C=C\\C=C/C");
        Bond b =fragment.findBond(2, 3);
        assertEquals("C", b.getBondStereoElement().getValue());
        b =fragment.findBond(4, 5);
        assertEquals("C", b.getBondStereoElement().getValue());
    }

    @Test
    public void testDoubleBondMulitStereo3() throws StructureBuildingException {
        Fragment fragment = fm.buildSMILES("F/C=C\\C=C\\C");
        Bond b =fragment.findBond(2, 3);
        assertEquals("C", b.getBondStereoElement().getValue());
        b =fragment.findBond(4, 5);
        assertEquals("T", b.getBondStereoElement().getValue());
    }

    @Test
    public void testDoubleBondMulitStereo4() throws StructureBuildingException {
        Fragment fragment = fm.buildSMILES("F/C=C\\C=CC");
        Bond b =fragment.findBond(2, 3);
        assertEquals("C", b.getBondStereoElement().getValue());
        b =fragment.findBond(4, 5);
        assertEquals(null, b.getBondStereoElement());
    }
}
