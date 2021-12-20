package uk.ac.cam.ch.wwmm.opsin;

class Attribute {

	private final String name;
	private String value;

	Attribute(String name, String value) {
		this.name = name;
		this.value = value;
	}

	/**
	 * Creates a copy
	 * @param attribute
	 */
	Attribute(Attribute attribute) {
		this.name = attribute.getName();
		this.value = attribute.getValue();
	}

	String getValue() {
		return value;
	}

	String getName() {
		return name;
	}
	
	void setValue(String value) {
		this.value = value;
	}

	String toXML() {
		return getName() + "=\"" + OpsinTools.xmlEncode(value) + "\"";
	}

	public String toString() {
		return name +"\t" + value;
	}
}