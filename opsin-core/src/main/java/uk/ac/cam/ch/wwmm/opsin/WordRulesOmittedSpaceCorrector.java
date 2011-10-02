package uk.ac.cam.ch.wwmm.opsin;

import java.util.List;

import nu.xom.Attribute;
import nu.xom.Element;
import nu.xom.Elements;
import static uk.ac.cam.ch.wwmm.opsin.XmlDeclarations.*;

public class WordRulesOmittedSpaceCorrector {
	private final Element parse;

	public WordRulesOmittedSpaceCorrector(Element parse) {
		this.parse = parse;
	}

	public void correctOmittedSpaces() {
		List<Element> wordRules = XOMTools.getDescendantElementsWithTagName(parse, WORDRULE_EL);
		for (Element wordRule : wordRules) {
			WordRule wordRuleVal = WordRule.valueOf(wordRule.getAttributeValue(WORDRULE_ATR));
			if (wordRuleVal == WordRule.divalentFunctionalGroup){
				checkAndCorrectOmittedSpacesInDivalentFunctionGroupRule(wordRule);
			}
			else if (wordRuleVal == WordRule.ester){
				
			}
		}
	}

	
	/**
	 * Corrects cases like "methylethyl ether" to "methyl ethyl ether"
	 * @param wordRule
	 */
	private void checkAndCorrectOmittedSpacesInDivalentFunctionGroupRule(Element divalentFunctionalGroupWordRule)  {
		List<Element> substituentWords = XOMTools.getChildElementsWithTagNameAndAttribute(divalentFunctionalGroupWordRule, WORD_EL, TYPE_ATR, SUBSTITUENT_TYPE_VAL);
		if (substituentWords.size()==1){//potentially has been "wrongly" interpreted e.g. ethylmethyl ketone is more likely to mean ethyl methyl ketone
			Elements children  =substituentWords.get(0).getChildElements();
			if (children.size()==2){
				Element firstSubstituent =(Element)children.get(0);
				//rule out correct usage e.g. diethyl ether and locanted substituents e.g. 2-methylpropyl ether
				if (firstSubstituent.getAttribute(LOCANT_ATR)==null && firstSubstituent.getAttribute(MULTIPLIER_ATR)==null){
					Element subToMove =children.get(1);
					subToMove.detach();
					Element newWord =new Element(WORD_EL);
					newWord.addAttribute(new Attribute(TYPE_ATR, SUBSTITUENT_TYPE_VAL));
					newWord.appendChild(subToMove);
					XOMTools.insertAfter(substituentWords.get(0), newWord);
				}
			}
		}
	}

}
