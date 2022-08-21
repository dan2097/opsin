package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

class GroupingEl extends Element{
	
	private final List<Element> children = new ArrayList<>();
	
	GroupingEl(String name) {
		super(name);
	}

	@Override
	void addChild(Element child) {
		child.setParent(this);
		children.add(child);
	}
	
	@Override
	Element copy() {
		GroupingEl copy = new GroupingEl(this.name);
		for (Element childEl : this.children) {
			Element newChild = childEl.copy();
			newChild.setParent(copy);
			copy.addChild(newChild);
		}
		for (int i = 0, len = this.attributes.size(); i < len; i++) {
			Attribute atr = this.attributes.get(i);
			copy.addAttribute(new Attribute(atr));
		}
		return copy;
	}
	
	@Override
	Element getChild(int index) {
		return children.get(index);
	}
	
	@Override
	int getChildCount() {
		return children.size();
	}

	@Override
	List<Element> getChildElements() {
		return new ArrayList<>(children);
	}
	
	@Override
	List<Element> getChildElements(String name) {
		List<Element> elements = new ArrayList<>(1);
		for (Element element : children) {
			if (element.name.equals(name)) {
				elements.add(element);
			}
		}
		return elements;
	}

	@Override
	Element getFirstChildElement(String name) {
		for (Element child : children) {
			if (child.getName().equals(name)) {
				return child;
			}
		}
		return null;
	}

	String getValue() {
		int childCount = getChildCount();
		if (childCount == 0) {
			return "";
		}
		StringBuilder result = new StringBuilder();
		for (int i = 0; i < childCount; i++) {
			result.append(children.get(i).getValue());
		}
		return result.toString();
	}
	
	@Override
	int indexOf(Element child) {
		return children.indexOf(child);
	}

	@Override
	void insertChild(Element child, int index) {
		child.setParent(this);
		children.add(index, child);
	}

	@Override
	boolean removeChild(Element child) {
		child.setParent(null);
		return children.remove(child);
	}
	
	@Override
	Element removeChild(int index) {
		Element removed = children.remove(index);
		removed.setParent(null);
		return removed;
	}
	
	@Override
	void replaceChild(Element oldChild, Element newChild) {
		int index = indexOf(oldChild);
		if (index == -1) {
			throw new RuntimeException("oldChild is not a child of this element.");
		}
		removeChild(index);
		insertChild(newChild, index);
	}
	
	void setValue(String text) {
		throw new UnsupportedOperationException("Token groups do not have a value");
	}

	@Override
	public Iterator<Element> iterator() {
		return children.iterator();
	}

}
