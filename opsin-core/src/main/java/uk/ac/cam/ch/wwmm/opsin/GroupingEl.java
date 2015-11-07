package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Deque;
import java.util.List;

class GroupingEl extends Element{
	
	private final List<Element> children = new ArrayList<Element>();
	
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
		for (Attribute atr : this.attributes) {
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
		return new ArrayList<Element>(children);
	}
	
	@Override
	List<Element> getChildElements(String name) {
		List<Element> elements = new ArrayList<Element>(1);
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
		if (children.size() == 0) {
			return "";
		}
		StringBuilder result = new StringBuilder();
		Deque<Element> stack = new ArrayDeque<Element>();
		for (int i = children.size() -1; i >= 0; i--) {
			stack.add(children.get(i));
		}
		while (stack.size() > 0){
			Element currentElement = stack.removeLast();
			result.append(currentElement.getValue());
			for (int i = currentElement.getChildCount() -1; i >= 0; i--) {
				stack.add(currentElement.getChild(i));
			}
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

}
