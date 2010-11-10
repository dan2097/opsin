package uk.ac.cam.ch.wwmm.opsin;
import static junit.framework.Assert.assertEquals;


import org.junit.Test;

public class ParserTest {

	@Test
	public void testConvertStringToComponentRatios1() throws ParsingException {
		String ratio = "(1:2)";
		Integer[] componentRatios = Parser.processStoichometryIndication(ratio);
		assertEquals(2, componentRatios.length);
		for (int i = 0; i < componentRatios.length; i++) {
			if (i==0){
				assertEquals(1,(int) componentRatios[i]);
			}
			if (i==1){
				assertEquals(2,(int) componentRatios[i]);
			}
		}
	}
	
	@Test
	public void testConvertStringToComponentRatios2() throws ParsingException {
		String ratio = "[1/1/2]";
		Integer[] componentRatios = Parser.processStoichometryIndication(ratio);
		assertEquals(3, componentRatios.length);
		for (int i = 0; i < componentRatios.length; i++) {
			if (i==0){
				assertEquals(1,(int) componentRatios[i]);
			}
			if (i==1){
				assertEquals(1,(int) componentRatios[i]);
			}
			if (i==2){
				assertEquals(2,(int) componentRatios[i]);
			}
		}
	}
	
	@Test
	public void testConvertStringToComponentRatios3() throws ParsingException {
		String ratio = "(1:2:?)";
		Integer[] componentRatios = Parser.processStoichometryIndication(ratio);
		assertEquals(3, componentRatios.length);
		for (int i = 0; i < componentRatios.length; i++) {
			if (i==0){
				assertEquals(1,(int) componentRatios[i]);
			}
			if (i==1){
				assertEquals(2,(int) componentRatios[i]);
			}
			if (i==2){
				assertEquals(1,(int) componentRatios[i]);
			}
		}
	}
}
