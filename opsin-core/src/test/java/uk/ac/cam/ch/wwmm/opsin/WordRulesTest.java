package uk.ac.cam.ch.wwmm.opsin;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static uk.ac.cam.ch.wwmm.opsin.XmlDeclarations.*;

import java.util.List;

import org.junit.jupiter.api.Test;

public class WordRulesTest {
	@Test
	public void testOlateSuffixDoesNotMatchAteGroup() throws Exception {
		for (String spelling : new String[]{"olate", "olat"}) {
			Element suffix = new TokenEl(SUFFIX_EL, spelling);
			suffix.addAttribute(VALUE_ATR, "olate");
			List<Element> wordRules = groupWithAlcohol(suffix);
			assertEquals(2, wordRules.size(), spelling);
			for (Element wordRule : wordRules) {
				assertEquals(WordRule.simple.toString(), wordRule.getAttributeValue(WORDRULE_ATR));
			}
		}
	}

	@Test
	public void testEsterFormingSuffixesMatchAteGroup() throws Exception {
		for (String value : new String[]{"ate", "carboxylate", "ite", "amide", "glycoside"}) {
			Element suffix = new TokenEl(SUFFIX_EL, value);
			suffix.addAttribute(VALUE_ATR, value);
			List<Element> wordRules = groupWithAlcohol(suffix);
			assertEquals(1, wordRules.size(), value);
			assertEquals(WordRule.potentialAlcoholEster.toString(), wordRules.get(0).getAttributeValue(WORDRULE_ATR), value);
		}
	}

	@Test
	public void testPidolateGroupStillMatchesAteGroup() throws Exception {
		for (String spelling : new String[]{"pidolate", "pidolat"}) {
			Element group = new TokenEl(GROUP_EL, spelling);
			group.addAttribute(FUNCTIONALIDS_ATR, "9");
			List<Element> wordRules = groupWithAlcohol(group);
			assertEquals(1, wordRules.size(), spelling);
			assertEquals(WordRule.potentialAlcoholEster.toString(), wordRules.get(0).getAttributeValue(WORDRULE_ATR), spelling);
		}
	}

	private List<Element> groupWithAlcohol(Element endingToken) throws Exception {
		Element molecule = new GroupingEl(MOLECULE_EL);
		molecule.addChild(makeWord("methanol", new TokenEl(GROUP_EL, "methanol")));
		//Use a single token interpretation so alternative tokenisations cannot mask a predicate failure.
		molecule.addChild(makeWord("example", endingToken));
		WordRules rules = new WordRules(new ResourceGetter("uk/ac/cam/ch/wwmm/opsin/resources/"));
		rules.groupWordsIntoWordRules(molecule, NameToStructureConfig.getDefaultConfigInstance(), false, null);
		return molecule.getChildElements(WORDRULE_EL);
	}

	private Element makeWord(String value, Element endingToken) {
		Element word = new GroupingEl(WORD_EL);
		word.addAttribute(TYPE_ATR, WordType.full.toString());
		word.addAttribute(VALUE_ATR, value);
		Element root = new GroupingEl(ROOT_EL);
		root.addChild(endingToken);
		word.addChild(root);
		return word;
	}
}
