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
		if (regexTokenElement.getAttribute("value")!=null){
			elem.addAttribute(new Attribute(VALUE_ATR, regexTokenElement.getAttributeValue("value")));
		}
		if (regexTokenElement.getAttribute("type")!=null){
			elem.addAttribute(new Attribute(TYPE_ATR, regexTokenElement.getAttributeValue("type")));
		}
		if (regexTokenElement.getAttribute("subType")!=null){
			elem.addAttribute(new Attribute(SUBTYPE_ATR, regexTokenElement.getAttributeValue("subType")));
		}
		if (regexTokenElement.getAttribute("ignoreWhenWritingXML")!=null && regexTokenElement.getAttributeValue("ignoreWhenWritingXML").equals("yes")){
			this.ignoreWhenWritingXML=true;
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
		if(tokenList.getAttribute("type") != null) {
			elem.addAttribute(new Attribute(TYPE_ATR, tokenList.getAttributeValue("type")));
		}
		if(tokenList.getAttribute("subType") != null) {
			elem.addAttribute(new Attribute(SUBTYPE_ATR, tokenList.getAttributeValue("subType")));
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
