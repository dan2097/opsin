package uk.ac.cam.ch.wwmm.opsin;

import java.util.List;

 class SuffixRule {

	private final SuffixRuleType type;
	private final List<Attribute> attributes;
	
	SuffixRule(SuffixRuleType type, List<Attribute> attributes) {
		this.type = type;
		this.attributes = attributes;
	}
	
	SuffixRuleType getType() {
		return type;
	}
	
	/**
	 * Returns the value of the attribute with the given name
	 * or null if the attribute doesn't exist
	 * @param name
	 * @return
	 */
	String getAttributeValue(String name) {
		for (Attribute a : attributes) {
			if (a.getName().equals(name)) {
				return a.getValue();
			}
		}
		return null;
	}
}
