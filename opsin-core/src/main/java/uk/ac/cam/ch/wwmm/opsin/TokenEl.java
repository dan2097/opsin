package uk.ac.cam.ch.wwmm.opsin;

import java.util.Collections;
import java.util.Iterator;
import java.util.List;

class TokenEl extends Element {
	
	private String value;
	private Fragment frag;

	TokenEl(String name) {
		super(name);
		this.value = "";
	}
	
	TokenEl(String name, String value) {
		super(name);
		this.value = value;
	}

	@Override
	void addChild(Element child) {
		throw new UnsupportedOperationException("Tokens do not have children");
	}
	
	@Override
	Element copy() {
		TokenEl copy = new TokenEl(this.name, this.value);
		for (int i = 0, len = this.attributes.size(); i < len; i++) {
			Attribute atr = this.attributes.get(i);
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
		for (int i = 0, len = this.attributes.size(); i < len; i++) {
			Attribute atr = this.attributes.get(i);
			copy.addAttribute(new Attribute(atr));
		}
		return copy;
	}
	
	@Override
	Element getChild(int index) {
		throw new UnsupportedOperationException("Tokens do not have children");
	}

	@Override
	int getChildCount() {
		return 0;
	}

	@Override
	List<Element> getChildElements() {
		return Collections.emptyList();
	}

	@Override
	List<Element> getChildElements(String name) {
		return Collections.emptyList();
	}

	@Override
	Element getFirstChildElement(String name) {
		return null;
	}
	
	@Override
	Element getLastChildElement() {
		return null;
	}
	
	@Override
	Fragment getFrag() {
		return frag;
	}
	
	String getValue() {
		return value;
	}

	@Override
	int indexOf(Element child) {
		return -1;
	}

	@Override
	void insertChild(Element child, int index) {
		throw new UnsupportedOperationException("Tokens do not have children");
	}
	
	@Override
	boolean removeChild(Element child) {
		throw new UnsupportedOperationException("Tokens do not have children");
	}
	
	@Override
	Element removeChild(int index) {
		throw new UnsupportedOperationException("Tokens do not have children");
	}

	@Override
	void replaceChild(Element oldChild, Element newChild) {
		throw new UnsupportedOperationException("Tokens do not have children");
	}
	
	@Override
	void setFrag(Fragment frag) {
		this.frag = frag;
	}

	void setValue(String text) {
		this.value = text;
	}

	@Override
	public Iterator<Element> iterator() {
		return Collections.emptyIterator();
	}

}
