package uk.ac.cam.ch.wwmm.opsin;

import static uk.ac.cam.ch.wwmm.opsin.XmlDeclarations.*;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import nu.xom.Document;
import nu.xom.Element;
import nu.xom.Elements;

class SuffixRules {
	
	/**For a given group type what suffixes are applicable. Due to group subTypes altering suffix meaning, the same suffixValue maps to one or more suffixes*/
	private final HashMap<String, HashMap<String, List<Element>>> suffixApplicability;
	/**A mapping between suffix rule names and elements containing the rules for applying the corresponding suffix*/
	private final HashMap<String, Element> suffixRules;
	
	SuffixRules(ResourceGetter resourceGetter) throws IOException{
		suffixApplicability = generateSuffixApplicabilityMap(resourceGetter);
		suffixRules = generateSuffixRulesMap(resourceGetter);
	}

	private HashMap<String, HashMap<String, List<Element>>> generateSuffixApplicabilityMap(ResourceGetter resourceGetter) throws IOException {
		Document suffixApplicabilityDoc = resourceGetter.getXMLDocument("suffixApplicability.xml");
		HashMap<String, HashMap<String, List<Element>>> suffixApplicability = new HashMap<String, HashMap<String,List<Element>>>();
		Elements groupTypes = suffixApplicabilityDoc.getRootElement().getChildElements(SUFFIXAPPLICABILITY_GROUPTYPE_EL);
		for (int i = 0; i < groupTypes.size(); i++) {
			Element groupType =groupTypes.get(i);
			Elements suffixes = groupType.getChildElements(SUFFIXAPPLICABILITY_SUFFIX_EL);
			HashMap<String, List<Element>> suffixToRuleMap= new HashMap<String, List<Element>>();
			for (int j = 0; j < suffixes.size(); j++) {
				Element suffix =suffixes.get(j);
				String suffixValue= suffix.getAttributeValue(SUFFIXAPPLICABILITY_VALUE_ATR);
				if (suffixToRuleMap.get(suffixValue)!=null){//can have multiple entries if subType attribute is set
					suffixToRuleMap.get(suffixValue).add(suffix);
				}
				else{
					List<Element> suffixList =new ArrayList<Element>();
					suffixList.add(suffix);
					suffixToRuleMap.put(suffixValue, suffixList);
				}
			}
			suffixApplicability.put(groupType.getAttributeValue(SUFFIXAPPLICABILITY_TYPE_ATR), suffixToRuleMap);
		}
		return suffixApplicability;
	}

	private HashMap<String, Element> generateSuffixRulesMap(ResourceGetter resourceGetter) throws IOException {
		Document suffixRulesDoc = resourceGetter.getXMLDocument("suffixRules.xml");
		HashMap<String, Element> suffixRules = new HashMap<String, Element>();
		Elements rules = suffixRulesDoc.getRootElement().getChildElements(SUFFIXRULES_RULE_EL);
		for (int i = 0; i < rules.size(); i++) {
			Element rule =rules.get(i);
			String ruleValue=rule.getAttributeValue(SUFFIXRULES_VALUE_ATR);
			if (suffixRules.get(ruleValue)!=null){
				throw new RuntimeException("Suffix: " +ruleValue +" appears multiple times in suffixRules.xml");
			}
			suffixRules.put(ruleValue, rule);
		}
		return suffixRules;
	}
	
	/**
	 * Returns the appropriate suffixRule tags for the given arguments.
	 * The suffix rule tags are the children of the appropriate rule in suffixRules.xml
	 * @param suffixTypeToUse
	 * @param suffixValue
	 * @param subgroupType
	 * @return
	 * @throws ComponentGenerationException
	 */
	Elements getSuffixRuleTags(String suffixTypeToUse, String suffixValue, String subgroupType) throws ComponentGenerationException {
		HashMap<String, List<Element>> groupToSuffixMap = suffixApplicability.get(suffixTypeToUse);
		if (groupToSuffixMap==null){
			throw new ComponentGenerationException("Suffix Type: "+ suffixTypeToUse + " does not have a corresponding groupType entry in suffixApplicability.xml");
		}
		List<Element> potentiallyApplicableSuffixes =groupToSuffixMap.get(suffixValue);
		if(potentiallyApplicableSuffixes==null || potentiallyApplicableSuffixes.size()==0 ) {
			throw new ComponentGenerationException("Suffix: " +suffixValue +" does not apply to the group it was associated with (type: "+  suffixTypeToUse + ")according to suffixApplicability.xml");
		}
		Element chosenSuffix=null;
        for (Element suffix : potentiallyApplicableSuffixes) {
            if (suffix.getAttribute(SUFFIXAPPLICABILITY_SUBTYPE_ATR) != null) {
                if (!suffix.getAttributeValue(SUFFIXAPPLICABILITY_SUBTYPE_ATR).equals(subgroupType)) {
                    continue;
                }
            }
            if (chosenSuffix != null) {
                throw new ComponentGenerationException("Suffix: " + suffixValue + " appears multiple times in suffixApplicability.xml");
            }
            chosenSuffix = suffix;
        }
		if (chosenSuffix==null){
			throw new ComponentGenerationException("Suffix: " +suffixValue +" does not apply to the group it was associated with (type: "+  suffixTypeToUse + ")due to the group's subType: "+ subgroupType +" according to suffixApplicability.xml");
		}
		Element rule =suffixRules.get(chosenSuffix.getValue());
		if(rule ==null) {
			throw new ComponentGenerationException("Suffix: " +chosenSuffix.getValue() +" does not have a rule associated with it in suffixRules.xml");
		}
		return rule.getChildElements();
	}
	
	/**
	 * Does suffixApplicability.xml have an entry for this group type? 
	 * @param groupType
	 * @return
	 */
	boolean isGroupTypeWithSpecificSuffixRules(String groupType){
		return suffixApplicability.containsKey(groupType);
	}
}
