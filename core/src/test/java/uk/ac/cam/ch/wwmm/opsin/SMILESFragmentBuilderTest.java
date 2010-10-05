package uk.ac.cam.ch.wwmm.opsin;

import static junit.framework.Assert.assertNotNull;
import static org.junit.Assert.*;
import static org.mockito.Mockito.mock;

import java.util.List;
import java.util.Set;

import junit.framework.Assert;

import org.junit.Before;
import org.junit.Ignore;
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

    @Test(expected=StructureBuildingException.class)
    public void badlyFormedSMILE4() throws StructureBuildingException {
        fm.buildSMILES("C=#C");
        Assert.fail("Should throw exception for bad smiles: is it a double or triple bond?");
    }
    
    @Test(expected=StructureBuildingException.class)
    public void badlyFormedSMILE5() throws StructureBuildingException {
        fm.buildSMILES("C#=C");
        Assert.fail("Should throw exception for bad smiles: is it a double or triple bond?");
    }
    
    @Test(expected=StructureBuildingException.class)
    public void badlyFormedSMILE6() throws StructureBuildingException {
        fm.buildSMILES("F//C=C/F");
        Assert.fail("Should throw exception for bad smiles: bond configuration specified twice");
    }

	
    @Test(expected=StructureBuildingException.class)
    public void badlyFormedSMILE7() throws StructureBuildingException {
        fm.buildSMILES("F/C=C/\\F");
        Assert.fail("Should throw exception for bad smiles: bond configuration specified twice");
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
	public void hydrogenHandling1() throws StructureBuildingException {
		Fragment fragment = fm.buildSMILES("[OH3+]");
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(1, atomList.size());
		assertEquals(1, atomList.get(0).getCharge());
		assertEquals(3, atomList.get(0).determineValency(true));
	}
	
	@Test
	public void hydrogenHandling2() throws StructureBuildingException {
		Fragment fragment = fm.buildSMILES("[CH3][CH2][OH]");
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(3, atomList.size());
		assertEquals(4, atomList.get(0).determineValency(true));
		assertEquals(4, atomList.get(1).determineValency(true));
		assertEquals(2, atomList.get(2).determineValency(true));
	}
	
	@Test
	public void hydrogenHandling3() throws StructureBuildingException {
		Fragment fragment = fm.buildSMILES("[SH2]");
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(1, atomList.size());
		assertEquals(2, atomList.get(0).determineValency(true));
	}
	
	@Test
	public void hydrogenHandling4() throws StructureBuildingException {
		Fragment fragment = fm.buildSMILES("[SH4]");
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(1, atomList.size());
		int minimumVal =atomList.get(0).getMinimumValency();
		assertEquals(4, minimumVal);
		assertEquals(4, atomList.get(0).determineValency(true));
	}
	
	@Test
	public void hydrogenHandling5() throws StructureBuildingException {
		Fragment fragment = fm.buildSMILES("[SH6]");
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(1, atomList.size());
		int minimumVal =atomList.get(0).getMinimumValency();
		assertEquals(6, minimumVal);
		assertEquals(6, atomList.get(0).determineValency(true));
	}
	
	@Test
	public void hydrogenHandling6() throws StructureBuildingException {
		Fragment fragment = fm.buildSMILES("[SH3]");
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(1, atomList.size());
		int minimumVal =atomList.get(0).getMinimumValency();
		assertEquals(3, minimumVal);
		assertEquals(3, atomList.get(0).determineValency(true));
	}
	
	
	@Test
	public void hydrogenHandling7() throws StructureBuildingException {
		Fragment fragment = fm.buildSMILES("[SH3+]");
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(1, atomList.size());
		assertEquals(1, atomList.get(0).getCharge());
		assertEquals(3, atomList.get(0).determineValency(true));
	}
	
	@Test
	public void hydrogenHandling8() throws StructureBuildingException {
		Fragment fragment = fm.buildSMILES("[SH+]");
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(1, atomList.size());
		assertEquals(1, atomList.get(0).getCharge());
		assertEquals(1, atomList.get(0).determineValency(true));
	}
	
	@Test
	public void hydrogenHandling9() throws StructureBuildingException {
		Fragment fragment = fm.buildSMILES("[SH3-]");
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(1, atomList.size());
		assertEquals(-1, atomList.get(0).getCharge());
		assertEquals(3, atomList.get(0).determineValency(true));
	}
	
	@Test
	public void hydrogenHandling10() throws StructureBuildingException {
		Fragment fragment = fm.buildSMILES("[SH-]");
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(1, atomList.size());
		assertEquals(-1, atomList.get(0).getCharge());
		assertEquals(1, atomList.get(0).determineValency(true));
	}
	
	@Test
	public void hydrogenHandling11() throws StructureBuildingException {
		Fragment fragment = fm.buildSMILES("[SH5+]");
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(1, atomList.size());
		int minimumVal =atomList.get(0).getMinimumValency();
		assertEquals(5, minimumVal);
		assertEquals(1, atomList.get(0).getCharge());
		assertEquals(5, atomList.get(0).determineValency(true));
	}
	
	@Test
	public void hydrogenHandling12() throws StructureBuildingException {
		Fragment fragment = fm.buildSMILES("[Li+]");
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(1, atomList.size());
		assertEquals(1, atomList.get(0).getCharge());
		assertEquals(0, atomList.get(0).determineValency(true));
	}
	
	@Test
	public void hydrogenHandling13() throws StructureBuildingException {
		Fragment fragment = fm.buildSMILES("[NaH]");
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(2, atomList.size());
		assertEquals(0, atomList.get(0).getProtonsExplicitlyAddedOrRemoved());
		assertEquals(0, atomList.get(0).getCharge());
		
		assertEquals(0, atomList.get(1).getProtonsExplicitlyAddedOrRemoved());
		assertEquals(0, atomList.get(1).getCharge());
		assertEquals("H", atomList.get(1).getElement());
	}
	
	@Test
	public void hydrogenHandling14() throws StructureBuildingException {
		Fragment fragment = fm.buildSMILES("-[SiH3]");
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(1, atomList.size());
		assertEquals(4, atomList.get(0).determineValency(true));
	}
	
	@Test
	public void hydrogenHandling15() throws StructureBuildingException {
		Fragment fragment = fm.buildSMILES("=[SiH2]");
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(1, atomList.size());
		assertEquals(4, atomList.get(0).determineValency(true));
	}
	
	
	@Test
	public void hydrogenHandling16() throws StructureBuildingException {
		Fragment fragment = fm.buildSMILES("#[SiH]");
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(1, atomList.size());
		assertEquals(4, atomList.get(0).determineValency(true));
	}
	
	@Test
	public void hydrogenHandling17() throws StructureBuildingException {
		Fragment fragment = fm.buildSMILES("[SiH3]-");
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(1, atomList.size());
		assertEquals(4, atomList.get(0).determineValency(true));
	}
	
	@Test
	public void hydrogenHandling18() throws StructureBuildingException {
		Fragment fragment = fm.buildSMILES("[SiH2]=");
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(1, atomList.size());
		assertEquals(4, atomList.get(0).determineValency(true));
	}
	
	@Test
	public void hydrogenHandling19() throws StructureBuildingException {
		Fragment fragment = fm.buildSMILES("[SiH]#");
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(1, atomList.size());
		assertEquals(4, atomList.get(0).determineValency(true));
	}
	
	@Test
	public void hydrogenHandling20() throws StructureBuildingException {
		Fragment fragment = fm.buildSMILES("=[Si]=");
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(1, atomList.size());
		assertEquals(4, atomList.get(0).determineValency(true));
	}
	
	@Test
	public void indicatedHydrogen1() throws StructureBuildingException {
		Fragment fragment = fm.buildSMILES("Nc1[nH]c(=O)c2c(n1)nc[nH]2");
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(11, atomList.size());
		assertEquals(2, fragment.getIndicatedHydrogen().size());
		assertEquals(atomList.get(2), fragment.getIndicatedHydrogen().get(0));
		assertEquals(atomList.get(10),  fragment.getIndicatedHydrogen().get(1));
	}
	
	@Test
	@Ignore
	public void indicatedHydrogen2() throws StructureBuildingException {
		Fragment fragment = fm.buildSMILES("NC(N1)=NC(N=CN2)=C2C1=O");
		List<Atom> atomList = fragment.getAtomList();
		assertEquals(11, atomList.size());
		FragmentTools.convertHighOrderBondsToSpareValencies(fragment);
		assertEquals(2, fragment.getIndicatedHydrogen().size());
		assertEquals(atomList.get(2), fragment.getIndicatedHydrogen().get(0));
		assertEquals(atomList.get(7),  fragment.getIndicatedHydrogen().get(1));
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
    public void testDoubleBondMultiStereo1() throws StructureBuildingException {
        Fragment fragment = fm.buildSMILES("F/C=C/C=C/C");
        Bond b =fragment.findBond(2, 3);
        assertEquals("T", b.getBondStereoElement().getValue());
        b =fragment.findBond(4, 5);
        assertEquals("T", b.getBondStereoElement().getValue());
    }

    @Test
    public void testDoubleBondMultiStereo2() throws StructureBuildingException {
        Fragment fragment = fm.buildSMILES("F/C=C\\C=C/C");
        Bond b =fragment.findBond(2, 3);
        assertEquals("C", b.getBondStereoElement().getValue());
        b =fragment.findBond(4, 5);
        assertEquals("C", b.getBondStereoElement().getValue());
    }

    @Test
    public void testDoubleBondMultiStereo3() throws StructureBuildingException {
        Fragment fragment = fm.buildSMILES("F/C=C\\C=C\\C");
        Bond b =fragment.findBond(2, 3);
        assertEquals("C", b.getBondStereoElement().getValue());
        b =fragment.findBond(4, 5);
        assertEquals("T", b.getBondStereoElement().getValue());
    }

    @Test
    public void testDoubleBondMultiStereo4() throws StructureBuildingException {
        Fragment fragment = fm.buildSMILES("F/C=C\\C=CC");
        Bond b =fragment.findBond(2, 3);
        assertEquals("C", b.getBondStereoElement().getValue());
        b =fragment.findBond(4, 5);
        assertEquals(null, b.getBondStereoElement());
    }
    
    //From http://baoilleach.blogspot.com/2010/09/are-you-on-my-side-or-not-its-ez-part.html
    @Test
    public void testDoubleBondNoela() throws StructureBuildingException {
        Fragment fragment = fm.buildSMILES("C/C=C\\1/NC1");
        Bond b =fragment.findBond(2, 3);
        if ("T".equals( b.getBondStereoElement().getValue())){
            assertEquals("a1 a2 a3 a4", b.getBondStereoElement().getAttributeValue(XmlDeclarations.ATOMREFS4_ATR));
        }
        else{
        	assertEquals("a1 a2 a3 a5", b.getBondStereoElement().getAttributeValue(XmlDeclarations.ATOMREFS4_ATR));
        }
    }
    
    @Test
    public void testDoubleBondNoelb() throws StructureBuildingException {
        Fragment fragment = fm.buildSMILES("C/C=C1/NC1");
        Bond b =fragment.findBond(2, 3);
        assertEquals("T", b.getBondStereoElement().getValue());
        assertEquals("a1 a2 a3 a4", b.getBondStereoElement().getAttributeValue(XmlDeclarations.ATOMREFS4_ATR));
    }
    
    @Test
    public void testDoubleBondNoelc() throws StructureBuildingException {
        Fragment fragment = fm.buildSMILES("C/C=C\\1/NC/1");
        Bond b =fragment.findBond(2, 3);
        if ("T".equals( b.getBondStereoElement().getValue())){
            assertEquals("a1 a2 a3 a4", b.getBondStereoElement().getAttributeValue(XmlDeclarations.ATOMREFS4_ATR));
        }
        else{
        	assertEquals("a1 a2 a3 a5", b.getBondStereoElement().getAttributeValue(XmlDeclarations.ATOMREFS4_ATR));
        }
    }
    
    @Test
    public void testDoubleBondNoeld() throws StructureBuildingException {
        Fragment fragment = fm.buildSMILES("C/C=C1/NC/1");
        Bond b =fragment.findBond(2, 3);
        if ("T".equals( b.getBondStereoElement().getValue())){
            assertEquals("a1 a2 a3 a4", b.getBondStereoElement().getAttributeValue(XmlDeclarations.ATOMREFS4_ATR));
        }
        else{
        	assertEquals("a1 a2 a3 a5", b.getBondStereoElement().getAttributeValue(XmlDeclarations.ATOMREFS4_ATR));
        }
    }
    
    @Test(expected=StructureBuildingException.class)
    public void testDoubleBondNoele() throws StructureBuildingException {
        fm.buildSMILES("C/C=C\\1\\NC1");
        Assert.fail("Should throw exception for bad smiles: contradictory double bond configuration");
    }

    @Test(expected=StructureBuildingException.class)
    public void testDoubleBondNoelf() throws StructureBuildingException {
        fm.buildSMILES("C/C=C\1NC\1");
        Assert.fail("Should throw exception for bad smiles: contradictory double bond configuration");
    }
    
    @Test(expected=StructureBuildingException.class)
    public void testDoubleBondNoelg() throws StructureBuildingException {
        fm.buildSMILES("C/C=C\1/NC\1");
        Assert.fail("Should throw exception for bad smiles: contradictory double bond configuration");
    }
    
    @Test
    public void testDoubleBondNoelLike1() throws StructureBuildingException {
        Fragment fragment = fm.buildSMILES("C\\1NC1=C/C");
        Bond b =fragment.findBond(3, 4);
        assertEquals("C", b.getBondStereoElement().getValue());
        assertEquals("a1 a3 a4 a5", b.getBondStereoElement().getAttributeValue(XmlDeclarations.ATOMREFS4_ATR));
    }
    
    @Test
    public void testDoubleBondNoelLike2() throws StructureBuildingException {
        Fragment fragment = fm.buildSMILES("C1NC/1=C/C");
        Bond b =fragment.findBond(3, 4);
        assertEquals("C", b.getBondStereoElement().getValue());
        assertEquals("a1 a3 a4 a5", b.getBondStereoElement().getAttributeValue(XmlDeclarations.ATOMREFS4_ATR));
    }
}
