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

	void addAttribute(Attribute attribute) {
		attributes.add(attribute);
	}

	void addAttribute(String atrName, String atrValue) {
		attributes.add(new Attribute(atrName, atrValue));
	}

	/**
	 * Adds a child element
	 * @param child
	 */
	abstract void addChild(Element child);

	/**
	 * Creates a deep copy with no parent
	 */
	abstract Element copy();

	void detach() {
		if (parent != null) {
			parent.removeChild(this);
		}
	}
	
	Attribute getAttribute(int index) {
		return attributes.get(index);
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

	int getAttributeCount() {
		return attributes.size();
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

	/**
	 * Returns the child at the given index in the children list
	 * @param index
	 * @return
	 */
	abstract Element getChild(int index);

	/**
	 * Returns the number of children
	 * @return
	 */
	abstract int getChildCount();

	/**
	 * Returns a copy of the child elements
	 * 
	 * @return
	 */
	abstract List<Element> getChildElements();
	
	/**
	 * Gets child elements with this name (in iteration order)
	 * @param name
	 * @return
	 */
	abstract List<Element> getChildElements(String name);

	/**
	 * Returns the first child element with the specified name
	 * 
	 * @param name
	 * @return
	 */
	abstract Element getFirstChildElement(String name);
	
	/**
	 * Returns the fragment associated with this element (only applicable to group tokens!)
	 * @return
	 */
	Fragment getFrag() {
		throw new UnsupportedOperationException("Only group tokens can have associated fragments");
	}

	String getName() {
		return name;
	}

	Element getParent() {
		return this.parent;
	}

	abstract String getValue();

	/**
	 * Returns the index of the given child in the children list (or -1 if it isn't a child)
	 * @param child
	 * @return
	 */
	abstract int indexOf(Element child);

	/**
	 * Inserts the element at the given index in the children list
	 * @param child
	 * @param index
	 */
	abstract void insertChild(Element child, int index);

	boolean removeAttribute(Attribute attribute) {
		return attributes.remove(attribute);
	}

	/**
	 * Removes the given child element
	 * @param child
	 * @return
	 */
	abstract boolean removeChild(Element child);
	
	/**
	 * Removes the element at the given index in the children list
	 * @param index
	 * @return
	 */
	abstract Element removeChild(int index);
	
	/**
	 * Replaces a child element with another element
	 * @param oldChild
	 * @param newChild
	 */
	abstract void replaceChild(Element oldChild, Element newChild);
	
	/**
	 * Sets the fragment associated with this element (only applicable to group tokens!)
	 * @param frag
	 */
	void setFrag(Fragment frag) {
		throw new UnsupportedOperationException("Only group tokens can have associated fragments");
	}

	void setName(String name) {
		this.name = name;
	}

	void setParent(Element newParentEl) {
		this.parent = newParentEl;
	}

	abstract void setValue(String text);

	public String toString() {
		return toXML();
	}
	
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

}
