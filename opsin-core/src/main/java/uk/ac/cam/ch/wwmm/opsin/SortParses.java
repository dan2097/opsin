package uk.ac.cam.ch.wwmm.opsin;

import static uk.ac.cam.ch.wwmm.opsin.XmlDeclarations.WORDRULE_ATR;
import static uk.ac.cam.ch.wwmm.opsin.XmlDeclarations.WORDRULE_EL;

import java.util.Comparator;

/**
 * Prefer non-substituent word rules to substituent word rule e.g. ethylene is C=C not -CC-
 * Prefer the parse with the least elements that have 0 children e.g. benzal beats benz al (1 childless element vs 2 childless elements)
 * Prefer less elements e.g. <acryl(acidStem)amide(suffix)> beats <acryl(substituent)><amide(group)>
 */
class SortParses implements Comparator<Element> {
	public int compare(Element el1, Element el2){
		boolean isSubstituent1 = WordRule.substituent.toString().equals(el1.getFirstChildElement(WORDRULE_EL).getAttributeValue(WORDRULE_ATR));
		boolean isSubstituent2 = WordRule.substituent.toString().equals(el2.getFirstChildElement(WORDRULE_EL).getAttributeValue(WORDRULE_ATR));
		if (isSubstituent1 && !isSubstituent2){
			return 1;
		}
		if (!isSubstituent1 && isSubstituent2){
			return -1;
		}
		
		int[] counts1 = OpsinTools.countNumberOfElementsAndNumberOfChildLessElements(el1);
		int[] counts2 = OpsinTools.countNumberOfElementsAndNumberOfChildLessElements(el2);
		int childLessElementsInEl1 = counts1[1];
		int childLessElementsInEl2 = counts2[1];
		if ( childLessElementsInEl1> childLessElementsInEl2){
			return 1;
		}
		else if (childLessElementsInEl1 < childLessElementsInEl2){
			return -1;
		}

		int infixesInEl1 = countInfixes(el1);
		int infixesInEl2 = countInfixes(el2);
		if (infixesInEl1 > infixesInEl2){
			return 1;
		}
		else if (infixesInEl1 < infixesInEl2){
			return -1;
		}

		int elementsInEl1 = counts1[0];
		int elementsInEl2  = counts2[0];
		if ( elementsInEl1> elementsInEl2){
			return 1;
		}
		else if (elementsInEl1 < elementsInEl2){
			return -1;
		}
		else{
			return 0;
		}
	}

	/**
	 * How many infixes the parse uses.
	 *
	 * Tie-broken after the childless element count deliberately. An infix reading is more
	 * compact than the equivalent suffix plus substituent, so on total element count it always
	 * wins, which is how acetamidooxybenzene came to be read as acet + the amid infix applied
	 * to an oxy suffix, giving acetanilide with the oxygen gone rather than acetamido-oxy. Both
	 * readings are structurally legal, so the question is only which to prefer, and the reading
	 * that spends a token on modifying another token is the less likely one when a reading
	 * exists that takes the same tokens at face value.
	 *
	 * It has to come after the childless count rather than before it, so that it only separates
	 * parses that are already using the same number of tokens. chloromethanothiooxybenzene has
	 * an infix-free parse too, but one that needs an extra token and reads the o as a bridging
	 * oxygen, producing a benzocyclopropene; that parse loses on token count before this ever
	 * applies.
	 */
	private static int countInfixes(Element el) {
		int count = 0;
		if (el.getName().equals(XmlDeclarations.INFIX_EL)) {
			count++;
		}
		for (int i = 0, l = el.getChildCount(); i < l; i++) {
			count += countInfixes(el.getChild(i));
		}
		return count;
	}
}