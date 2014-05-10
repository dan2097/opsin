package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.List;

public class Element {

	private String localName;
	private String value;
	private Element parent = null;
	private List<Element> children = new ArrayList<Element>();
	private List<Attribute> attributes = new ArrayList<Attribute>();
	
	public Element(String name) {
		this.localName = name;
	}

	/**
	 * Creates a deep copy with no parent
	 * @param element
	 */
	public Element(Element element) {
		this.localName = element.localName;
		this.value = element.value;
		children = new ArrayList<Element>();
		for (Element childEl : element.children) {
			Element newChild = new Element(childEl);
			newChild.setParent(this);
			children.add(newChild);
		}
		
		attributes = new ArrayList<Attribute>();
		for (Attribute atr : element.attributes) {
			attributes.add(new Attribute(atr));
		}
	}

	/**
	 * Gets child elements with this name (in iteration order)
	 * @param name
	 * @return
	 */
	public List<Element> getChildElements(String name) {
		List<Element> elements = new ArrayList<Element>();
		for (Element element : children) {
			if (element.localName.equals(name)) {
				elements.add(element);
			}
		}
		return elements;
	}

	/**
	 * Returns a copy of the child elements
	 * 
	 * @return
	 */
	public List<Element> getChildElements() {
		return new ArrayList<Element>(children);
	}

	/**
	 * Returns the first child element with the specified name
	 * 
	 * @param name
	 * @return
	 */
	public Element getFirstChildElement(String name) {
		for (Element child : children) {
			if (child.getLocalName().equals(name)) {
				return child;
			}
		}
		return null;
	}


	public void addAttribute(Attribute attribute) {
		attributes.add(attribute);
	}
	
	public void addAttribute(String atrName, String atrValue) {
		attributes.add(new Attribute(atrName, atrValue));
	}

	
	public boolean removeAttribute(Attribute attribute) {
		return attributes.remove(attribute);
	}


	/**
	 * Returns the attribute with the given name
	 * or null if the attribute doesn't exist
	 * @param name
	 * @return
	 */
	public Attribute getAttribute(String name) {
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
	public String getAttributeValue(String name) {
		Attribute attribute = getAttribute(name);
		if (attribute != null) {
			return attribute.getValue();
		}
		return null;
	}

	public int getAttributeCount() {
		return attributes.size();
	}


	public Attribute getAttribute(int index) {
		return attributes.get(index);
	}


	public String getLocalName() {
		return localName;
	}

	public void setLocalName(String localName) {
		this.localName = localName;
	}


//TODO this method is inappropriately named
	public void appendChild(String text) {
		this.value = text;
	}

	public void removeChildren() {
		for (Element child : children) {
			child.setParent(null);
		}
		children.clear();
	}


	public String toXML() {
		return localName;

//		StringBuffer result = new StringBuffer(1024);
//		Node current = this;
//		boolean endTag = false;
//		int index = -1;
//		int[] indexes = new int[10];
//		int top = 0;
//		indexes[0] = -1;
//
//		while (true) {
//
//			if (!endTag && current.getChildCount() > 0) {
//				writeStartTag((Element) current, result);
//				current = current.getChild(0);
//				index = 0;
//				top++;
//				indexes = grow(indexes, top);
//				indexes[top] = 0;
//			} else {
//				if (endTag) {
//					writeEndTag((Element) current, result);
//					if (current == this)
//						break;
//				} else if (current.isElement()) {
//					writeStartTag((Element) current, result);
//					if (current == this)
//						break;
//				} else {
//					result.append(current.toXML());
//				}
//				endTag = false;
//				ParentNode parent = current.getParent();
//				if (parent.getChildCount() - 1 == index) {
//					current = parent;
//					top--;
//					if (current != this) {
//						index = indexes[top];
//					}
//					endTag = true;
//				} else {
//					index++;
//					indexes[top] = index;
//					current = parent.getChild(index);
//				}
//
//			}
//
//		}
//
//		return result.toString();

	}


	public String getValue() {
		//TODO should never be null?
		return value != null ? value : "";
	}

	public String toString() {
		return toXML();
	}

	public Element getParent() {
		return this.parent;
	}

	public void detach() {
		if (parent != null) {
			parent.removeChild(this);
		}
	}


	public int getChildCount() {
		return children.size();
	}


	public void insertChild(Element child, int position) {
		child.setParent(this);
		children.add(position, child);
	}


	public void appendChild(Element child) {
		child.setParent(this);
		children.add(child);
	}

	public Element getChild(int position) {
		return children.get(position);
	}


	public int indexOf(Element child) {
		return children.indexOf(child);
	}

	public Element removeChild(int i) {
		Element removed = children.remove(i);
		removed.setParent(null);
		return removed;
	}

	private void setParent(Element newParentEl) {
		this.parent = newParentEl;
	}


	public boolean removeChild(Element child) {
		child.setParent(null);
		return children.remove(child);
	}

	public void replaceChild(Element oldChild, Element newChild) {
		int position = indexOf(oldChild);
		if (position == -1) {
			throw new RuntimeException("oldChild is not a child of this element.");
		}

		removeChild(position);
		insertChild(newChild, position);
	}

}
