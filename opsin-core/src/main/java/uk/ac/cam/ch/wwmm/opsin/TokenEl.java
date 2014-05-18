package uk.ac.cam.ch.wwmm.opsin;

import java.util.Collections;
import java.util.List;

class TokenEl extends Element {
	
	private String value;

	TokenEl(String name, String value) {
		super(name);
		this.value = value;
	}
	
	TokenEl(String name) {
		super(name);
		this.value = "";
	}

	@Override
	List<Element> getChildElements(String name) {
		return Collections.emptyList();
	}
	
	@Override
	List<Element> getChildElements() {
		return Collections.emptyList();
	}
	
	@Override
	Element getFirstChildElement(String name) {
		return null;
	}
	
	@Override
	Element getChild(int position) {
		throw new UnsupportedOperationException("Tokens do not have children");
	}

	@Override
	int getChildCount() {
		return 0;
	}

	@Override
	void insertChild(Element child, int position) {
		throw new UnsupportedOperationException("Tokens do not have children");
	}

	@Override
	void appendChild(Element child) {
		throw new UnsupportedOperationException("Tokens do not have children");
	}

	@Override
	int indexOf(Element child) {
		throw new UnsupportedOperationException("Tokens do not have children");
	}
	
	@Override
	boolean removeChild(Element child) {
		throw new UnsupportedOperationException("Tokens do not have children");
	}

	@Override
	Element removeChild(int i) {
		throw new UnsupportedOperationException("Tokens do not have children");
	}

	@Override
	void replaceChild(Element oldChild, Element newChild) {
		throw new UnsupportedOperationException("Tokens do not have children");
	}
	
	void setValue(String text) {
		this.value = text;
	}
	
	String getValue() {
		return value;
	}

	@Override
	Element copy() {
		TokenEl copy = new TokenEl(this.name, this.value);
		for (Attribute atr : this.attributes) {
			copy.addAttribute(new Attribute(atr));
		}
		return copy;
	}

	/**
	 * Creates a copy with no parent
	 * The provided value is used instead of the Element to be copied's value
	 * @param value
	 * @return
	 */
	TokenEl copy(String value) {
		TokenEl copy = new TokenEl(this.name, value);
		for (Attribute atr : this.attributes) {
			copy.addAttribute(new Attribute(atr));
		}
		return copy;
	}

}
