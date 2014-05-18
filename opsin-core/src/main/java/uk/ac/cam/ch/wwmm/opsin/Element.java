package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.List;

abstract class Element {

	protected String name;
	protected Element parent = null;
	protected final List<Attribute> attributes = new ArrayList<Attribute>();

	Element(String name) {
		this.name = name;
	}

	/**
	 * Gets child elements with this name (in iteration order)
	 * @param name
	 * @return
	 */
	abstract List<Element> getChildElements(String name);

	/**
	 * Returns a copy of the child elements
	 * 
	 * @return
	 */
	abstract List<Element> getChildElements();

	/**
	 * Returns the first child element with the specified name
	 * 
	 * @param name
	 * @return
	 */
	abstract Element getFirstChildElement(String name);

	abstract Element getChild(int position);

	abstract int getChildCount();

	void addAttribute(Attribute attribute) {
		attributes.add(attribute);
	}
	
	void addAttribute(String atrName, String atrValue) {
		attributes.add(new Attribute(atrName, atrValue));
	}

	boolean removeAttribute(Attribute attribute) {
		return attributes.remove(attribute);
	}

	/**
	 * Returns the attribute with the given name
	 * or null if the attribute doesn't exist
	 * @param name
	 * @return
	 */
	Attribute getAttribute(String name) {
		for (Attribute a : attributes) {
			if (a.getName().equals(name)) {
				return a;
			}
		}
		return null;
	}

	/**
	 * Returns the value of the attribute with the given name
	 * or null if the attribute doesn't exist
	 * @param name
	 * @return
	 */
	String getAttributeValue(String name) {
		Attribute attribute = getAttribute(name);
		if (attribute != null) {
			return attribute.getValue();
		}
		return null;
	}

	int getAttributeCount() {
		return attributes.size();
	}

	Attribute getAttribute(int index) {
		return attributes.get(index);
	}

	String getName() {
		return name;
	}

	void setName(String name) {
		this.name = name;
	}

	abstract void setValue(String text);
	
	abstract String getValue();
	
	String toXML() {
		return toXML(0).toString();
	}

	private StringBuilder toXML(int indent) {
		StringBuilder result = new StringBuilder();
		for (int i = 0; i < indent; i++) {
			result.append("  ");
		}
		result.append('<');
		result.append(name);
		for (Attribute atr : attributes) {
			result.append(' ');
			result.append(atr.toXML());
		}
		result.append('>');
		if (getChildCount() > 0){
			for (Element child : getChildElements()) {
				result.append('\n');
				result.append(child.toXML(indent + 1));
			}
			result.append('\n');
			for (int i = 0; i < indent; i++) {
				result.append("  ");
			}
		}
		else{
			result.append(getValue());
		}
		result.append("</");
		result.append(name);
		result.append('>');

		return result;
	}

	Element getParent() {
		return this.parent;
	}

	void detach() {
		if (parent != null) {
			parent.removeChild(this);
		}
	}

	abstract void insertChild(Element child, int position);
	abstract void appendChild(Element child);

	abstract int indexOf(Element child);
	
	abstract boolean removeChild(Element child);

	abstract Element removeChild(int i);

	abstract void replaceChild(Element oldChild, Element newChild);
	
	void setParent(Element newParentEl) {
		this.parent = newParentEl;
	}
	
	/**
	 * Creates a deep copy with no parent
	 */
	abstract Element copy();
	
	public String toString() {
		return toXML();
	}

}
