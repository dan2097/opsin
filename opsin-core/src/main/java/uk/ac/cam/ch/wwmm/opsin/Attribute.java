package uk.ac.cam.ch.wwmm.opsin;

public class Attribute {

	private final String name;
	private String value;

	public Attribute(String name, String value) {
		this.name = name;
		this.value = value;
	}

	/**
	 * Creates a copy
	 * @param attribute
	 */
	public Attribute(Attribute attribute) {
		this.name = attribute.getName();
		this.value = attribute.getValue();
	}

	public String getValue() {
		return value;
	}

	public String getName() {
		return name;
	}
	
	public void setValue(String value) {
		this.value = value;
	}

	public String toXML() {
		return getName() + "=\"" + escapeText(value) + "\"";
	}

	public String toString() {
		return name +"\t" + value;
	}

	private String escapeText(String s) {
		StringBuffer result = new StringBuffer();
		for (int i = 0, l = s.length(); i < l; i++) {
			char c = s.charAt(i);
			switch (c) {
			case '\t':
				result.append("&#x09;");
				break;
			case '\n':
				result.append("&#x0A;");
				break;
			case '\r':
				result.append("&#x0D;");
				break;
			case '"':
				result.append("&quot;");
				break;
			case '&':
				result.append("&amp;");
				break;
			case '<':
				result.append("&lt;");
				break;
			case '>':
				result.append("&gt;");
				break;
			default:
				result.append(c);
			}
		}
		return result.toString();
	}
}