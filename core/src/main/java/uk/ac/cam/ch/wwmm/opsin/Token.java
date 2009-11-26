package uk.ac.cam.ch.wwmm.opsin;

import nu.xom.Attribute;
import nu.xom.Element;

/**A token in a chemical name. hex, yl, ane, chloro etc.
 * Stores information about the XML element that will be produced for the token.
 *
 * @author ptc24
 *
 */
class Token {

	/**A reference copy of the XML element to produce*/
	private Element elem;

	/**Should this token actually be used. Set to true for meaningless tokens e.g. e, o, endOfSubstituent etc.*/
	private boolean ignoreWhenWritingXML =false;

	/**Makes a new Token based on tagname, type, and ignoreWhenWritingXML attribute values. To be used for regex tokens.
	 * Type and ignoreWhenWritingXML may be null
	 *
	 * @param tagName
	 * @param type
	 * @param ignoreWhenWritingXML
	 */
	public Token(String tagName, String type, String ignoreWhenWritingXML) {
		elem = new Element(tagName);
		if (type!=null){
			elem.addAttribute(new Attribute("type",type));
		}
		if (ignoreWhenWritingXML!=null){
			if (ignoreWhenWritingXML.equals("yes")){
				this.ignoreWhenWritingXML=true;
			}
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
			elem.addAttribute(new Attribute("type", tokenList.getAttributeValue("type")));
		}
		if(tokenList.getAttribute("subType") != null) {
			elem.addAttribute(new Attribute("subType", tokenList.getAttributeValue("subType")));
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
