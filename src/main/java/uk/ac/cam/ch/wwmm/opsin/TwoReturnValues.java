package uk.ac.cam.ch.wwmm.opsin;

/**
 * Wrapper class for returning two objects
 *
 * @author dl387
 */
public final class TwoReturnValues<T, S> {

   public TwoReturnValues(T f, S s) {
      first = f;
      second = s;
   }

   public T getFirst() {
      return first;
   }

   public S getSecond() {
      return second;
   }

   public String toString() {
      return "(" + first.toString() + ", " + second.toString() + ")";
   }
   private T first;
   private S second;
}
