package uk.ac.cam.ch.wwmm.opsin;

import static uk.ac.cam.ch.wwmm.opsin.XmlDeclarations.WORDRULE_ATR;
import static uk.ac.cam.ch.wwmm.opsin.XmlDeclarations.WORDRULE_EL;

import java.util.Comparator;

import nu.xom.Element;
import uk.ac.cam.ch.wwmm.opsin.WordRules.WordRule;

/**
 * Prefer less childless elements e.g. benzal beats benz al
 * Prefer less elements e.g. <acryl(acidStem)amide(suffix)> beats <acryl(substituent)><amide(group)>
 */
class SortParses implements Comparator<Element>{
	public int compare(Element el1, Element el2){
		boolean isSubstituent1 = WordRule.substituent.toString().equals(el1.getFirstChildElement(WORDRULE_EL).getAttributeValue(WORDRULE_ATR));
		boolean isSubstituent2 = WordRule.substituent.toString().equals(el2.getFirstChildElement(WORDRULE_EL).getAttributeValue(WORDRULE_ATR));
		if (isSubstituent1 && !isSubstituent2){
			return 1;
		}
		if (!isSubstituent1 && isSubstituent2){
			return -1;
		}
		
		int[] counts1 = XOMTools.countNumberOfElementsAndNumberOfChildLessElements(el1);
		int[] counts2 = XOMTools.countNumberOfElementsAndNumberOfChildLessElements(el2);
		int childLessElementsInEl1 = counts1[1];
		int childLessElementsInEl2 = counts2[1];
		if ( childLessElementsInEl1> childLessElementsInEl2){
			return 1;
		}
		else if (childLessElementsInEl1 < childLessElementsInEl2){
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
}