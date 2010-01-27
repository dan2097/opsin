package uk.ac.cam.ch.wwmm.opsin;

/**
 * Wrapper class for returning three objects
 *
 * @author dl387
 */
final class ThreeReturnValues<F, S, T> {

   public ThreeReturnValues(F f, S s, T t) {
      first = f;
      second = s;
      third = t;
   }

   public F getFirst() {
      return first;
   }

   public S getSecond() {
      return second;
   }
   
   public T getThird() {
	   return third;
   }

   public String toString() {
      return "(" + first.toString() + ", " + second.toString() + ", " + third.toString() + ")";
   }
   private F first;
   private S second;
   private T third;
}
