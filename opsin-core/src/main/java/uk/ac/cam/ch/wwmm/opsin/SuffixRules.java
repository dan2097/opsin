package uk.ac.cam.ch.wwmm.opsin;

import static uk.ac.cam.ch.wwmm.opsin.XmlDeclarations.*;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.xml.stream.XMLStreamConstants;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamReader;


class SuffixRules {
	
	/**For a given group type what suffixes are applicable.
	 * Within this group type which are applicable for a given suffixValue
	 * Returns a list as different group subTypes can give different meanings*/
	private final Map<String, Map<String, List<ApplicableSuffix>>> suffixApplicability;

	private static class ApplicableSuffix {

		private final String requiredSubType;
		private final List<SuffixRule> suffixRules;

		public ApplicableSuffix(String requiredSubType, List<SuffixRule> suffixRules) {
			this.requiredSubType = requiredSubType;
			this.suffixRules = suffixRules;
		}
	}
	
	SuffixRules(ResourceGetter resourceGetter) throws IOException {
		Map<String, List<SuffixRule>> suffixRulesMap = generateSuffixRulesMap(resourceGetter);
		suffixApplicability = generateSuffixApplicabilityMap(resourceGetter, suffixRulesMap);
	}
	
	private Map<String, List<SuffixRule>> generateSuffixRulesMap(ResourceGetter resourceGetter) throws IOException {
		Map<String, List<SuffixRule>> suffixRulesMap = new HashMap<>();
		XMLStreamReader reader = resourceGetter.getXMLStreamReader("suffixRules.xml");
		try {
			while (reader.hasNext()) {
				if (reader.next() == XMLStreamConstants.START_ELEMENT && 
						reader.getLocalName().equals(SUFFIXRULES_RULE_EL)) {
					String ruleValue = reader.getAttributeValue(null, SUFFIXRULES_VALUE_ATR);
					if (suffixRulesMap.get(ruleValue) != null) {
						throw new RuntimeException("Suffix: " + ruleValue + " appears multiple times in suffixRules.xml");
					}
					suffixRulesMap.put(ruleValue, processSuffixRules(reader));
				}
			}
		}
		catch (XMLStreamException e) {
			throw new IOException("Parsing exception occurred while reading suffixRules.xml", e);
		}
		finally {
			try {
				reader.close();
			} catch (XMLStreamException e) {
				throw new IOException("Parsing exception occurred while reading suffixRules.xml", e);
			}
		}
		return suffixRulesMap;
	}
	

	private List<SuffixRule> processSuffixRules(XMLStreamReader reader) throws XMLStreamException {
		String startingElName = reader.getLocalName();
		List<SuffixRule> rules = new ArrayList<>();
		while (reader.hasNext()) {
			switch (reader.next()) {
			case XMLStreamConstants.START_ELEMENT:
				String tagName = reader.getLocalName();
				SuffixRuleType type = SuffixRuleType.valueOf(tagName);
				List<Attribute> attributes = new ArrayList<>();
				for (int i = 0, l = reader.getAttributeCount(); i < l; i++) {
					attributes.add(new Attribute(reader.getAttributeLocalName(i), reader.getAttributeValue(i)));
				}
				rules.add(new SuffixRule(type, attributes));
				break;
			case XMLStreamConstants.END_ELEMENT:
				if (reader.getLocalName().equals(startingElName)) {
					return rules;
				}
				break;
			}
		}
		throw new RuntimeException("Malformed suffixRules.xml");
	}

	private Map<String, Map<String, List<ApplicableSuffix>>> generateSuffixApplicabilityMap(ResourceGetter resourceGetter, Map<String, List<SuffixRule>> suffixRulesMap) throws IOException {
		Map<String, Map<String, List<ApplicableSuffix>>> suffixApplicability = new HashMap<>();
		XMLStreamReader reader = resourceGetter.getXMLStreamReader("suffixApplicability.xml");
		try {
			while (reader.hasNext()) {
				if (reader.next() == XMLStreamConstants.START_ELEMENT && 
						reader.getLocalName().equals(SUFFIXAPPLICABILITY_GROUPTYPE_EL)) {
					Map<String, List<ApplicableSuffix>> suffixToRuleMap = new HashMap<>();
					suffixApplicability.put(reader.getAttributeValue(null, SUFFIXAPPLICABILITY_TYPE_ATR), suffixToRuleMap);
					while (reader.hasNext()) {
						int event = reader.next();
						if (event == XMLStreamConstants.START_ELEMENT &&
							reader.getLocalName().equals(SUFFIXAPPLICABILITY_SUFFIX_EL)) {
							String suffixValue = reader.getAttributeValue(null, SUFFIXAPPLICABILITY_VALUE_ATR);
							List<ApplicableSuffix> suffixList = suffixToRuleMap.get(suffixValue);
							//can have multiple entries if subType attribute is set
							if (suffixToRuleMap.get(suffixValue) == null){
								suffixList = new ArrayList<>();
								suffixToRuleMap.put(suffixValue, suffixList);
							}
							String requiredSubType = reader.getAttributeValue(null, SUFFIXAPPLICABILITY_SUBTYPE_ATR);
							String suffixRuleName = reader.getElementText();
							List<SuffixRule> suffixRules = suffixRulesMap.get(suffixRuleName);
							if (suffixRules == null) {
								throw new RuntimeException("Suffix: " + suffixRuleName +" does not have a rule associated with it in suffixRules.xml");
							}
							suffixList.add(new ApplicableSuffix(requiredSubType, suffixRules));
						}
						else if (event == XMLStreamConstants.END_ELEMENT &&
								reader.getLocalName().equals(SUFFIXAPPLICABILITY_GROUPTYPE_EL)) {
							break;
						}
					}
				}
			}
		}
		catch (XMLStreamException e) {
			throw new IOException("Parsing exception occurred while reading suffixApplicability.xml", e);
		}
		finally {
			try {
				reader.close();
			} catch (XMLStreamException e) {
				throw new IOException("Parsing exception occurred while reading suffixApplicability.xml", e);
			}
		}
		return suffixApplicability;
	}
	

	/**
	 * Returns the appropriate suffixRules for the given arguments.
	 * The suffix rules are the children of the appropriate rule in suffixRules.xml
	 * @param suffixTypeToUse
	 * @param suffixValue
	 * @param subgroupType
	 * @return
	 * @throws ComponentGenerationException
	 */
	List<SuffixRule> getSuffixRuleTags(String suffixTypeToUse, String suffixValue, String subgroupType) throws ComponentGenerationException {
		Map<String, List<ApplicableSuffix>> groupToSuffixMap = suffixApplicability.get(suffixTypeToUse);
		if (groupToSuffixMap == null){
			throw new ComponentGenerationException("Suffix Type: " + suffixTypeToUse + " does not have a corresponding groupType entry in suffixApplicability.xml");
		}
		List<ApplicableSuffix> potentiallyApplicableSuffixes = groupToSuffixMap.get(suffixValue);
		if(potentiallyApplicableSuffixes == null || potentiallyApplicableSuffixes.isEmpty() ) {
			throw new ComponentGenerationException("Suffix: " + suffixValue + " does not apply to the group it was associated with (type: " +	suffixTypeToUse + ") according to suffixApplicability.xml");
		}
		List<SuffixRule> suffixRules = null;
		for (ApplicableSuffix suffix : potentiallyApplicableSuffixes) {
			if (suffix.requiredSubType != null) {
				if (!suffix.requiredSubType.equals(subgroupType)) {
					continue;
				}
			}
			if (suffixRules != null) {
				throw new ComponentGenerationException("Suffix: " + suffixValue + " appears multiple times in suffixApplicability.xml");
			}
			suffixRules = suffix.suffixRules;
		}
		if (suffixRules == null){
			throw new ComponentGenerationException("Suffix: " +suffixValue +" does not apply to the group it was associated with (type: "+	suffixTypeToUse + ") due to the group's subType: "+ subgroupType +" according to suffixApplicability.xml");
		}
		return suffixRules;
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
