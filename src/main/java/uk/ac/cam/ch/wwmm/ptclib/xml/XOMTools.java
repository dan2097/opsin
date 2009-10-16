package uk.ac.cam.ch.wwmm.ptclib.xml;

import java.io.ByteArrayOutputStream;

import nu.xom.Attribute;
import nu.xom.Document;
import nu.xom.Element;
import nu.xom.Elements;
import nu.xom.Node;
import nu.xom.ParentNode;
import nu.xom.Serializer;
import nu.xom.Text;

/** Accessory functions for the manipulation of XOM Elements, Nodes,
 * Documents etc. Originally leaned from Chemname, now added to.
 * 
 * @author Peter Corbett, 2005, based on chemname sources by PMR et al.
 */
public final class XOMTools {
	
    /**Gets the next sibling of a given node.
    *
    * @param node The reference node.
    * @return The next Sibling, or null.
    */	
	public static Node getNextSibling(Node node) {
		ParentNode parent = node.getParent();
		int i = parent.indexOf(node);
		if (i+1 >= parent.getChildCount()) return null;
		return parent.getChild(i+1);
	}
	
	/**Gets the first next sibling of a given node whose tagname matches the given string.
    *
    * @param current The reference node.
    * @param tagName The tagname of a node to look for
    * @return The matched next Sibling, or null.
    */
	public static Node getNextSibling(Node current, String tagName) {
		Element matchedElement =null;
		while (true) {
			Element next = (Element) getNextSibling(current);
			if (next != null) {
				if (next.getLocalName().equals(tagName)){
					matchedElement=next;
					break;
				}
				else{
					current = next;
				}
			} else {
				break;
			}
		}
		return matchedElement;
	}


	/**Gets the previous sibling of a given node.
    *
    * @param node The reference node.
    * @return The previous Sibling, or null.
    */
	public static Node getPreviousSibling(Node node) {
		ParentNode parent = node.getParent();
		int i = parent.indexOf(node);
		if (i==0) return null;
		return parent.getChild(i-1);
	}
	
	
	/**Gets the first previous sibling of a given node whose tagname matches the given string.
    *
    * @param current The reference node.
    * @param tagName The tagname of a node to look for
    * @return The matched previous Sibling, or null.
    */
	public static Node getPreviousSibling(Node current, String tagName) {
		Element matchedElement =null;
		while (true) {
			Element prev = (Element) getPreviousSibling(current);
			if (prev != null) {
				if (prev.getLocalName().equals(tagName)){
					matchedElement=prev;
					break;
				}
				else{
					current = prev;
				}
			} else {
				break;
			}
		}
		return matchedElement;
	}

	/**Inserts a node so that it occurs before a reference node. The new node
     * must not currently have a parent.
     *
     * @param node The reference node.
     * @param newNode The new node to insert.
     */
	public static void insertBefore(Node node, Node newNode) {
		ParentNode parent = node.getParent();
		int i = parent.indexOf(node);
		parent.insertChild(newNode, i);
	}

	/**Inserts a node so that it occurs after a reference node. The new node
     * must not currently have a parent.
     *
     * @param node The reference node.
     * @param newNode The new node to insert.
     */
	public static void insertAfter(Node node, Node newNode) {
		ParentNode parent = node.getParent();
		int i = parent.indexOf(node);
		parent.insertChild(newNode, i+1);
	}
 
	/**Makes a semi-shallow copy of an element, copying the element,
     * the namespace and the attribute, but no other child nodes.
     *
     * @param elem The element to copy.
     * @return The copied element.
     */	
	public static Element shallowCopy(Element elem) {
		Element newElem = new Element(elem.getLocalName());
		newElem.setBaseURI(elem.getBaseURI());
		newElem.setNamespaceURI(elem.getNamespaceURI());
		newElem.setNamespacePrefix(elem.getNamespacePrefix());
		for(int i=0;i<elem.getAttributeCount();i++) {
			newElem.addAttribute((Attribute)elem.getAttribute(i).copy());
		}
		return newElem;
	}
 
	/**Sets the namespace URI of an Element, and all child elements,
	 * recursively.
	 *
	 *@param elem The element to move into the namespace.
	 *@param nameSpace The namespace URI.
	 */
	public static void setNamespaceURIRecursively(Element elem, String nameSpace) {
		elem.setNamespaceURI(nameSpace);
		Elements children =elem.getChildElements();
		for(int i=0, n=children.size();i<n;i++)
			setNamespaceURIRecursively(children.get(i), nameSpace);
	}
	
	/**Removes the namespace URI of an element, and all child elements,
	 * recursively.
	 * 
	 * @param elem The element to remove the namespaces from.
	 */
	public static void clearNamespaces(Element elem) {
		if(elem.getNamespacePrefix() == null || elem.getNamespacePrefix().length() == 0) {
			elem.setNamespaceURI(null);
			Elements children =elem.getChildElements();
			for(int i=0, n=children.size();i<n;i++)
				clearNamespaces(children.get(i));			
		}
	}
	
	/**Removes an Element from a document, putting its child nodes into the gap
	 * left behind. Text nodes are sewn together so as to ensure that no new
	 * instances of a Text node next to another Text node are created.
	 * 
	 * @param elem The element to remove from the document.
	 */
	public static void removeElementPreservingText(Element elem) {
		ParentNode parent = elem.getParent();
		/* Put fingers on relevant nodes */
		Node previous = XOMTools.getPreviousSibling(elem);
		Node next = XOMTools.getNextSibling(elem);
		int index = elem.getParent().indexOf(elem);
		if(elem.getChildCount() > 0) {
			/* Put fingers on relevant inner nodes */
			Node contentsFirst = elem.getChild(0);
			Node contentsLast = elem.getChild(elem.getChildCount()-1);
			/* Transfer inner nodes to parent */
			while(elem.getChildCount() > 0) {
				Node nn = elem.getChild(0);
				nn.detach();
				parent.insertChild(nn, index);
				index++;
			}
			/* Perform Text surgery */
			if((previous instanceof Text) && (contentsFirst instanceof Text)) {
				((Text)previous).setValue(previous.getValue() + contentsFirst.getValue());
				if(contentsFirst == contentsLast) {
					contentsLast = previous;
				}
				contentsFirst.detach();
			} else if (contentsFirst instanceof Text) {
			}
			if((contentsLast instanceof Text) && (next instanceof Text)) {
				((Text)contentsLast).setValue(contentsLast.getValue() + next.getValue());
				next.detach();
			}
			elem.detach();
		} else {
			/* No inner nodes, just discard and maybe perform Text surgery */
			elem.detach();
			if((previous instanceof Text) && (next instanceof Text)) {
				((Text)previous).setValue(previous.getValue() + next.getValue());
				next.detach();
			}
		}
	}
	
	/**Examine all of the child and descendant nodes of this element,
	 * joining consecutive Text nodes together.
	 * 
	 * @param e The element to examine.
	 */
	public static void normalise(Element e) {
		int i = 0;
		while(i < e.getChildCount()) {
			if(e.getChild(i) instanceof Element) {
				normalise((Element)e.getChild(i));
				i++;
			} else if(e.getChild(i) instanceof Text) {
				if(e.getChildCount() > i+1 && e.getChild(i+1) instanceof Text) {
					Text t = (Text)e.getChild(i);
					t.setValue(t.getValue() + e.getChild(i+1).getValue());
					e.getChild(i+1).detach();
				} else {
					i++;
				}
			} else {
				i++;
			}			
		}
	}
	
	/**Copies a Node. Avoids bugs present in older versions of XOM.
	 * 
	 * @param n The node to copy.
	 * @return The copy of the node.
	 */
	public static Node safeCopy(Node n) {
		if(n instanceof Element) {
			Element e = (Element)n;
			Element ee = shallowCopy(e);
			for(int i=0;i<e.getAttributeCount();i++) {
				ee.addAttribute((Attribute)e.getAttribute(i).copy());
			}
			for(int i=0;i<e.getChildCount();i++) {
				ee.appendChild(safeCopy(e.getChild(i)));
			}
			return ee;
		} else {
			return n.copy();
		}
	}
	
	/**Gets the XPoint of a given node, given a root node.
	 * 
	 * @param n The node to make the XPoint for.
	 * @param root The root node.
	 * @return The XPoint
	 */
	public static String getXPointToNode(Node n, Node root) {
		if(n == root) return "/1";
		return getXPointToNode(n.getParent(), root) + "/" + Integer.toString(n.getParent().indexOf(n)+1);		
	}

	/**Finds the node at a given XPoint.
	 * 
	 * @param n The root node that the XPoint is relative to.
	 * @param xPoint The XPoint.
	 * @return The Node.
	 */
	public static Node getNodeAtXPoint(Node n, String xPoint) {
		try {
			if(xPoint.equals("/1")) return n;
			xPoint = xPoint.substring(1);
			xPoint = xPoint.substring(xPoint.indexOf('/'));
			return getNodeAtXPointInternal(n, xPoint);
		} catch (Exception e) {
			return null;
		}
	}
	
	private static Node getNodeAtXPointInternal(Node n, String xPoint) {
		try {
			xPoint = xPoint.substring(1);
			if(xPoint.indexOf('/') == -1) {
				return n.getChild(Integer.parseInt(xPoint)-1);
			}
			String nextXPoint = xPoint.substring(xPoint.indexOf('/'));		
			return getNodeAtXPointInternal(n.getChild(Integer.parseInt(xPoint.substring(0, xPoint.indexOf('/')))-1), nextXPoint);
		} catch (Exception e) {
			return null;
		}
	}

	/**Produces a hash value for an XML document.
	 * 
	 * @param doc The document to get a hash value for.
	 * @return The hash value.
	 */
	public static int documentHash(Document doc) {
		try {
			ByteArrayOutputStream baos = new ByteArrayOutputStream();
			Serializer ser = new Serializer(baos);
			ser.write(doc);
			return baos.toString().hashCode();
		} catch (Exception e) {
			throw new Error(e);
		}
	}
	
}
