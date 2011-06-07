package uk.ac.cam.ch.wwmm.opsin;
import org.junit.Test;

import static junit.framework.Assert.*;

public class FusedRingNumbererFunctionsTest {


	@Test
	public void testGetOppositeDirection(){
		assertEquals(4,FusedRingNumberer.getOppositeDirection(0));
		assertEquals(-3,FusedRingNumberer.getOppositeDirection(1));
		assertEquals(-2,FusedRingNumberer.getOppositeDirection(2));
		assertEquals(-1,FusedRingNumberer.getOppositeDirection(3));
		assertEquals(0,FusedRingNumberer.getOppositeDirection(4));
		assertEquals(0,FusedRingNumberer.getOppositeDirection(-4));
		assertEquals(1,FusedRingNumberer.getOppositeDirection(-3));
		assertEquals(2,FusedRingNumberer.getOppositeDirection(-2));
		assertEquals(3,FusedRingNumberer.getOppositeDirection(-1));
	}
//	
//	@Test
//	public void testDetermineAbsoluteDirectionFromPreviousDirection3Membered(){
//		assertEquals(1,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(1, 0, 3));
//		assertEquals(-1,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(-1, 0, 3));
//
//		assertEquals(2,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(1, 1, 3));
//		assertEquals(0,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(-1, 1, 3));
//
//		assertEquals(0,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(1, -1, 3));
//		assertEquals(-2,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(-1, -1, 3));
//	}
//	
//	@Test
//	public void testDetermineAbsoluteDirectionFromPreviousDirection4Membered(){
//		assertEquals(0,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(0, 0, 4));
//		assertEquals(2,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(2, 0, 4));
//		assertEquals(-2,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(-2, 0, 4));
//		
//		assertEquals(2,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(0, 2, 4));
//		assertEquals(4,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(2, 2, 4));
//		assertEquals(0,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(-2, 2, 4));
//		
//		assertEquals(-2,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(0, -2, 4));
//		assertEquals(0,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(2, -2, 4));
//		assertEquals(4,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(-2, -2, 4));
//	}
//	
//	@Test
//	public void testDetermineAbsoluteDirectionFromPreviousDirection5Membered(){
//		assertEquals(0,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(0, 0, 5));
//		assertEquals(1,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(1, 0, 5));
//		assertEquals(2,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(2, 0, 5));
//		assertEquals(3,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(3, 0, 5));
//		assertEquals(-1,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(-1, 0, 5));
//		assertEquals(-2,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(-2, 0, 5));
//		assertEquals(-3,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(-3, 0, 5));
//
//		//assertEquals(,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(0, 1, 5));
//		assertEquals(2,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(1, 1, 5));
//		//assertEquals(3,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(2, 1, 5));
//		assertEquals(4,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(3, 1, 5));
//		assertEquals(0,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(-1, 1, 5));
//		//assertEquals(-1,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(-2, 1, 5));
//		assertEquals(-1,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(-3, 1, 5));
//
//		assertEquals(2,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(0, 2, 5));
//		assertEquals(3,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(1, 2, 5));
//		assertEquals(4,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(2, 2, 5));
//		assertEquals(-3,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(3, 2, 5));
//		assertEquals(1,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(-1, 2, 5));
//		assertEquals(0,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(-2, 2, 5));
//		assertEquals(-1,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(-3, 2, 5));
//
//		//assertEquals(3,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(0, 3, 5));
//		assertEquals(4,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(1, 3, 5));
//		//assertEquals(,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(2, 3, 5));
//		assertEquals(-2,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(3, 3, 5));
//		assertEquals(2,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(-1, 3, 5));
//		assertEquals(1,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(-2, 3, 5));
//		//assertEquals(,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(-3, 3, 5));
////
//		assertEquals(0,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(0, -1, 5));
////		assertEquals(,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(1, -1, 5));
////		assertEquals(,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(2, -1, 5));
////		assertEquals(,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(3, -1, 5));
////		assertEquals(,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(-1, -1, 5));
////		assertEquals(,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(-2, -1, 5));
////		assertEquals(,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(-3, -1, 5));
////
////		assertEquals(-2,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(0, -2, 5));
////		assertEquals(,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(1, -2, 5));
////		assertEquals(,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(2, -2, 5));
////		assertEquals(.FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(3, -2, 5));
////		assertEquals(,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(-1, -2, 5));
////		assertEquals(,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(-2, -2, 5));
////		assertEquals(,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(-3, -2, 5));
////		
////		assertEquals(-3,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(0, -3, 5));
////		assertEquals(,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(1, -3, 5));
////		assertEquals(,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(2, -3, 5));
////		assertEquals(,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(3, -3, 5));
////		assertEquals(,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(-1, -3, 5));
////		assertEquals(,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(-2, -3, 5));
////		assertEquals(,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(-3, -3, 5));
//	}
//	
//	@Test
//	public void testDetermineAbsoluteDirectionFromPreviousDirection8Membered(){
//		assertEquals(0,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(0, 0, 8));
//		assertEquals(1,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(1, 0, 8));
//		assertEquals(3,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(3, 0, 8));
//		assertEquals(-1,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(-1, 0, 8));
//		assertEquals(-3,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(-3, 0, 8));
//		
//		assertEquals(1,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(0, 1, 8));
//		assertEquals(3,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(1, 1, 8));
//		assertEquals(4,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(3, 1, 8));
//		assertEquals(0,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(-1, 1, 8));
//		assertEquals(-1,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(-3, 1, 8));
//		
//		assertEquals(1,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(0, 1, 8));
//		assertEquals(3,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(1, 1, 8));
//		assertEquals(4,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(3, 1, 8));
//		assertEquals(0,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(-1, 1, 8));
//		assertEquals(-1,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(-3, 1, 8));
//		
//		assertEquals(2,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(0, 2, 8));
//		assertEquals(3,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(1, 2, 8));
//		assertEquals(-3,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(3, 2, 8));
//		assertEquals(1,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(-1, 2, 8));
//		assertEquals(-1,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(-3, 2, 8));
//		
//		assertEquals(3,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(0, 3, 8));
//		assertEquals(4,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(1, 3, 8));
//		assertEquals(-3,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(3, 3, 8));
//		assertEquals(1,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(-1, 3, 8));
//		assertEquals(1,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(-1, 3, 8));
//		assertEquals(0,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(-3, 3, 8));
//		
//		assertEquals(-1,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(0, -1, 8));
//		assertEquals(0,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(1, -1, 8));
//		assertEquals(1,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(3, -1, 8));
//		assertEquals(-3,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(-1, -1, 8));
//		assertEquals(4,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(-3, -1, 8));
//
//		assertEquals(-2,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(0, -2, 8));
//		assertEquals(-1,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(1, -2, 8));
//		assertEquals(1,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(3, -2, 8));
//		assertEquals(-3,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(-1, -2, 8));
//		assertEquals(3,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(-3, -2, 8));
//		
//		assertEquals(-3,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(0, -3, 8));
//		assertEquals(-1,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(1, -3, 8));
//		assertEquals(0,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(3, -3, 8));
//		assertEquals(4,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(-1, -3, 8));
//		assertEquals(3,FusedRingNumberer.determineAbsoluteDirectionFromPreviousDirection(-3, -3, 8));
//	}
}
