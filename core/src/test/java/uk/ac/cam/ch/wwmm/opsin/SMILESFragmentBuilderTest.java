package uk.ac.cam.ch.wwmm.opsin;

import static junit.framework.Assert.assertNotNull;
import static org.junit.Assert.assertEquals;
import static org.mockito.Mockito.mock;

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
		fragment = sBuilder.build("CC", fm);
		assertNotNull("Got a fragment", fragment);
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
