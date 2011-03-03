package uk.ac.cam.ch.wwmm.opsin;

import java.util.List;

/**
 * A wrapper for the results from parsing a chemical name or part of a chemical name
 * through ParseRules
 *
 * @author dl387
 */
public class ParseRulesResults {
   private final List<ParseTokens> parseTokensList;
   private final String uninterpretableName;
   private final String unparseableName;
	
   ParseRulesResults(List<ParseTokens> parseTokensList, String uninterpretableName, String unparseableName) {
	  this.parseTokensList = parseTokensList;
	  this.uninterpretableName = uninterpretableName;
	  this.unparseableName = unparseableName;
   }

   /**
    * One ParseTokens object is returned for each possible interpretation of a chemical name
    * If none of the name can be interpreted this list will be empty
    * @return
    */
   public List<ParseTokens> getParseTokensList() {
      return parseTokensList;
   }

   /**
    * The substring of the name that could not be classified into a substituent/full/functionalTerm
    * e.g. in ethyl-2H-fooarene  "2H-fooarene" will be returned
    * @return
    */
   public String getUninterpretableName() {
      return uninterpretableName;
   }
   
   /**
    * The substring of the name that could not be tokenised at all.
    * This will always be the same or shorter than the uninterpetable substring of name
    * e.g. in ethyl-2H-fooarene  "fooarene" will be returned
    * @return
    */
   public String getUnparseableName() {
	   return unparseableName;
   }

   public String toString() {
      return "(" + parseTokensList.toString() + ", " + uninterpretableName + ", " + unparseableName + ")";
   }

}
