package uk.ac.cam.ch.wwmm.opsin;

import nu.xom.Builder;

import org.junit.Test;


public class StructureBuilderTest {

	static Builder XMLBuilder = new Builder();
	
	//TODO how can this be unit tested in any meaningful way?
	
	@Test
	public void testBuilderFromName() {
		new StructureBuilder();
	}
}
