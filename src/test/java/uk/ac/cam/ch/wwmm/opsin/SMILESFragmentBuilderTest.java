package uk.ac.cam.ch.wwmm.opsin;

import static junit.framework.Assert.assertNotNull;
import static org.junit.Assert.assertEquals;

import org.junit.Before;
import org.junit.Test;


public class SMILESFragmentBuilderTest {

	Fragment fragment;
	SMILESFragmentBuilder builder;
	
	@Before
	public void setUp() throws Exception {
		builder = new SMILESFragmentBuilder();
	}

	@Test
	public void testBuild() throws StructureBuildingException {		
		fragment = builder.build("CC", new IDManager());
		assertNotNull("Got a fragment", fragment);
	}
	
	@Test
    public void testDoubleBondStereo1() throws StructureBuildingException {
        fragment = builder.build("F/C=C/F", new IDManager());
        Bond b =fragment.findBond(2, 3);
        assertEquals("T", b.getBondStereoElement().getValue());
    }

    @Test
    public void testDoubleBondStereo2() throws StructureBuildingException {
        fragment = builder.build("F\\C=C/F", new IDManager());
        Bond b =fragment.findBond(2, 3);
        assertEquals("C", b.getBondStereoElement().getValue());
    }

    @Test
    public void testDoubleBondStereo3() throws StructureBuildingException {
        fragment = builder.build("C(/F)=C/F", new IDManager());
        Bond b =fragment.findBond(1, 3);
        assertEquals("C", b.getBondStereoElement().getValue());
    }

    @Test
    public void testDoubleBondStereo4() throws StructureBuildingException {
        fragment = builder.build("C(\\F)=C/F", new IDManager());
        Bond b =fragment.findBond(1, 3);
        assertEquals("T", b.getBondStereoElement().getValue());
    }

    @Test
    public void testDoubleBondStereo5a() throws StructureBuildingException {
        fragment = builder.build("CC1=C/F.O\\1", new IDManager());
        Bond b =fragment.findBond(2, 3);
        assertEquals("C", b.getBondStereoElement().getValue());
    }

    @Test
    public void testDoubleBondStereo5b() throws StructureBuildingException {
        fragment = builder.build("CC/1=C/F.O1", new IDManager());
        Bond b =fragment.findBond(2, 3);
        assertEquals("C", b.getBondStereoElement().getValue());
    }

    @Test
    public void testDoubleBondStereo6() throws StructureBuildingException {
        fragment = builder.build("CC1=C/F.O/1", new IDManager());
        Bond b =fragment.findBond(2, 3);
        assertEquals("T", b.getBondStereoElement().getValue());
    }

    @Test
    public void testDoubleBondMulitStereo1() throws StructureBuildingException {
        fragment = builder.build("F/C=C/C=C/C", new IDManager());
        Bond b =fragment.findBond(2, 3);
        assertEquals("T", b.getBondStereoElement().getValue());
        b =fragment.findBond(4, 5);
        assertEquals("T", b.getBondStereoElement().getValue());
    }

    @Test
    public void testDoubleBondMulitStereo2() throws StructureBuildingException {
        fragment = builder.build("F/C=C\\C=C/C", new IDManager());
        Bond b =fragment.findBond(2, 3);
        assertEquals("C", b.getBondStereoElement().getValue());
        b =fragment.findBond(4, 5);
        assertEquals("C", b.getBondStereoElement().getValue());
    }

    @Test
    public void testDoubleBondMulitStereo3() throws StructureBuildingException {
        fragment = builder.build("F/C=C\\C=C\\C", new IDManager());
        Bond b =fragment.findBond(2, 3);
        assertEquals("C", b.getBondStereoElement().getValue());
        b =fragment.findBond(4, 5);
        assertEquals("T", b.getBondStereoElement().getValue());
    }

    @Test
    public void testDoubleBondMulitStereo4() throws StructureBuildingException {
        fragment = builder.build("F/C=C\\C=CC", new IDManager());
        Bond b =fragment.findBond(2, 3);
        assertEquals("C", b.getBondStereoElement().getValue());
        b =fragment.findBond(4, 5);
        assertEquals(null, b.getBondStereoElement());
    } 
}
