package uk.ac.cam.ch.wwmm.opsin;

import nu.xom.Attribute;
import nu.xom.Element;
import static uk.ac.cam.ch.wwmm.opsin.XmlDeclarations.*;

/**A token in a chemical name. hex, yl, ane, chloro etc.
 * Stores information about the XML element that will be produced for the token.
 *
 * @author ptc24
 * @author dl387
 *
 */
class Token {

	/**A reference copy of the XML element to produce*/
	private final Element elem;

	/**Should this token actually be used. Set to true for meaningless tokens e.g. e, o, endOfSubstituent etc.*/
	private boolean ignoreWhenWritingXML =false;


	/**
	 * Makes a new token using a regexToken Element
	 * @param regexTokenElement
	 */
	Token(Element regexTokenElement) {
		elem = new Element(regexTokenElement.getAttributeValue("tagname"));
		String value = regexTokenElement.getAttributeValue("value");
		String type = regexTokenElement.getAttributeValue("type");
		String subType = regexTokenElement.getAttributeValue("subType");
		if (value != null){
			elem.addAttribute(new Attribute(VALUE_ATR, value));
		}
		if (type != null){
			elem.addAttribute(new Attribute(TYPE_ATR, type));
		}
		if (subType != null){
			elem.addAttribute(new Attribute(SUBTYPE_ATR, subType));
		}
		if ("yes".equals(regexTokenElement.getAttributeValue("ignoreWhenWritingXML"))){
			ignoreWhenWritingXML = true;
		}
	}

	/**Makes a new Token based on reference elements from an XML file.
	 *
	 * @param tokenElement The token element in the XML tokens file.
	 * @param tokenList The tokenList element the token was taken from.
	 */
	Token(Element tokenElement, Element tokenList) {
		elem = OpsinTools.shallowCopy(tokenElement);
		elem.setLocalName(tokenList.getAttributeValue("tagname"));
		String type = tokenList.getAttributeValue("type");
		String subType = tokenList.getAttributeValue("subType");
		if (type != null){
			elem.addAttribute(new Attribute(TYPE_ATR, type));
		}
		if (subType != null){
			elem.addAttribute(new Attribute(SUBTYPE_ATR, subType));
		}
		if ("yes".equals(tokenList.getAttributeValue("ignoreWhenWritingXML"))){
			ignoreWhenWritingXML = true;
		}
	}

	/**Makes an XML element of the token.
	 *
	 * @param text The string to go in the Text node contained within the Element.
	 * @return The element produced.
	 */
	Element makeElement(String text) {
		if (!ignoreWhenWritingXML){
			Element tokenElement = OpsinTools.shallowCopy(elem);
			tokenElement.appendChild(text);
			return tokenElement;
		}
		else{
			return null;
		}
	}
}
