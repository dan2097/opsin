package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import uk.ac.cam.ch.wwmm.ptclib.xml.XOMTools;

import nu.xom.Attribute;
import nu.xom.Element;
import nu.xom.Elements;
import nu.xom.Node;
import nu.xom.ParentNode;
import nu.xom.Text;

/**
 * A set of useful methods to assist OPSIN
 * @author dl387
 *
 */
class OpsinTools {

	/**
	 * Returns the next sibling suffix node which is not related to altering charge (ium/ide/id)
	 * @param group
	 */
	public static Element getNextNonChargeSuffix(Element current) {
		Element matchedElement =null;
		while (true) {
			Element next = (Element) XOMTools.getNextSibling(current);
			if (next != null) {
				if (next.getLocalName().equals("suffix")){
					if (next.getAttribute("subType")==null || !next.getAttributeValue("subType").equals("charge")){
						matchedElement=next;
						break;
					}
				}
				current = next;
			} else {
				break;
			}
		}
		return matchedElement;
	}

	/**
	 * Returns an arrayList of elements corresponding to the Elements given
	 * @param elements
	 * @return The new arrayList
	 */
	public static ArrayList<Element> elementsToElementArrayList(Elements elements) {
		ArrayList<Element> elementList =new ArrayList<Element>(elements.size());
		for (int i = 0, n=elements.size(); i < n; i++) {
			elementList.add((Element) elements.get(i));
		}
		return elementList;
	}

	/**
	 * Returns a new list containing the elements of list1 followed by list2
	 * @param list1
	 * @param list2
	 * @return The new list
	 */
	public static ArrayList<Element> combineElementLists(List<Element> list1, List<Element> list2) {
		ArrayList<Element> elementList =new ArrayList<Element>(list1);
		elementList.addAll(list2);
		return elementList;
	}

	/**
	 * Gets the next node. This element need not be a sibling
	 * @param current: starting node
	 * @return
	 */
	public static Node getNext(Node node) {
		Element parent = (Element) node.getParent();
		if (parent == null || parent.getLocalName().equals("molecule")){
			return null;
		}
		int index = parent.indexOf(node);
		if (index +1 >=parent.getChildCount()) return getNext(parent);//reached end of element
		Element next =(Element) parent.getChild(index+1);
		Elements children =next.getChildElements();
		while (children.size()!=0){
			next =children.get(0);
			children =next.getChildElements();
		}
		return next;
	}

	/**
	 * Gets the previous node. This element need not be a sibling
	 * @param current: starting node
	 * @return
	 */
	public static Node getPrevious(Node node) {
		Element parent = (Element) node.getParent();
		if (parent == null || parent.getLocalName().equals("molecule")){
			return null;
		}
		int index = parent.indexOf(node);
		if (index ==0) return getPrevious(parent);//reached beginning of element
		Element previous =(Element) parent.getChild(index-1);
		Elements children =previous.getChildElements();
		while (children.size()!=0){
			previous =children.get(children.size()-1);
			children =previous.getChildElements();
		}
		return previous;
	}

	/**
	 * Gets the next node:
	 * The siblings of the given node will be surveyed.
	 * If no sibling follows null will be returned
	 * Then the first node of the node following the given node will be returned (children are recursively evaluated)
	 * @param current: starting node
	 * @return
	 */
	public static Node getNextAtSameOrLowerLevel(Node node) {
		Element parent = (Element) node.getParent();
		if (parent == null || parent.getLocalName().equals("molecule")){
			return null;
		}
		int index = parent.indexOf(node);
		if (index +1 >=parent.getChildCount()) return null;
		Element next =(Element) parent.getChild(index+1);
		Elements children =next.getChildElements();
		while (children.size()!=0){
			next =children.get(0);
			children =next.getChildElements();
		}
		return next;
	}

	/**
	 * Returns the previous group. This group element need not be a sibling
	 * @param current: starting node
	 * @return
	 */
	public static Node getPreviousGroup(Element current) {
	  if (current.getLocalName().equals("group")){//can start with a group or the sub/root the group is in
		  current=(Element)current.getParent();
	  }
	  Element parent = (Element) current.getParent();
	  if (parent == null || parent.getLocalName().equals("molecule")){
		  return null;
	  }
	  int index = parent.indexOf(current);
	  if (index ==0) return getPreviousGroup(parent);//no group found
	  Element previous =(Element) parent.getChild(index-1);
	  Elements children =previous.getChildElements();
	  while (children.size()!=0){
		  previous =children.get(children.size()-1);
		  children =previous.getChildElements();
	  }
	  Elements groups =((Element)previous.getParent()).getChildElements("group");
	  if (groups.size()==0){
		  return getPreviousGroup(previous);
	  }
	  else{
		  return groups.get(groups.size()-1);//return last group if multiple exist e.g. fused ring
	  }
	}


	/**
	 * If a dash is the last character it is removed
	 * @param locantText
	 * @return
	 */
	public static String removeDashIfPresent(String locantText){
		if(locantText.endsWith("-")) {
			locantText = locantText.substring(0, locantText.length()-1);
		}
		return locantText;
	}

	/**
	 * Sets the first text child of the group to the newName
	 * Throws an exception if the first child is not a Text node
	 * @param group
	 * @param newName
	 * @throws PostProcessingException
	 */
	public static void setTextChild(Element group, String newName) throws PostProcessingException {
		Node textNode =group.getChild(0);
		if (textNode instanceof Text){
			((Text)textNode).setValue(newName);
		}
		else{
			throw new PostProcessingException("No Text Child Found!");
		}
	}

	/**
	 * Finds the word element that encloses the given element.
	 * Returns the word element or throws an exception
	 * @param Element el
	 * @return word Element
	 * @throws PostProcessingException
	 */
	public static Element getParentWord(Element el) throws PostProcessingException {
		Element parent=(Element)el.getParent();
		while(parent !=null && !parent.getLocalName().equals("word")){
			parent =(Element)parent.getParent();
		}
		if (parent==null){
			throw new PostProcessingException("Cannot find enclosing word element");
		}
		else{
			return parent;
		}
	}

	/**
	 * Returns an arrayList containing sibling elements of the given type after the given element.
	 * These elements need not be continuous
	 * @param currentElem: the element to look for following siblings of
	 * @param type: the "localname" of the element type desired
	 * @return
	 * @throws PostProcessingException 
	 */
	public static ArrayList<Element> getLaterSiblingElementsOfType(Element currentElem, String type) throws PostProcessingException {
		ArrayList<Element> laterSiblingElementsOfType= new ArrayList<Element>();
		Element parent =(Element) currentElem.getParent();
		if (parent==null){
			throw new PostProcessingException("Node has no parent!");
		}
		Elements potentialMatches =parent.getChildElements(type);
		int indexOfCurrentElem =parent.indexOf(currentElem);
		for (int i = 0; i < potentialMatches.size(); i++) {
			if (parent.indexOf(potentialMatches.get(i)) > indexOfCurrentElem){
				laterSiblingElementsOfType.add(potentialMatches.get(i));
			}
		}
		return laterSiblingElementsOfType;
	}
	
	/**
	 * Returns an arrayList containing sibling elements of the given type before the given element.
	 * These elements need not be continuous
	 * @param currentElem: the element to look for previous siblings of
	 * @param type: the "localname" of the element type desired
	 * @return
	 * @throws PostProcessingException 
	 */
	public static ArrayList<Element> getEarlierSiblingElementsOfType(Element currentElem, String type) throws PostProcessingException {
		ArrayList<Element> earlierSiblingElementsOfType= new ArrayList<Element>();
		Element parent =(Element) currentElem.getParent();
		if (parent==null){
			throw new PostProcessingException("Node has no parent!");
		}
		Elements potentialMatches =parent.getChildElements(type);
		int indexOfCurrentElem =parent.indexOf(currentElem);
		for (int i = 0; i < potentialMatches.size(); i++) {
			if (parent.indexOf(potentialMatches.get(i)) < indexOfCurrentElem){
				earlierSiblingElementsOfType.add(potentialMatches.get(i));
			}
		}
		return earlierSiblingElementsOfType;
	}

	/**
	 * Gets the next sibling element of the given element. If this element's name is within the elementsToIgnore array this is repeated
	 * If no appropriate element can be found null is returned
	 * @param el
	 * @param elementsToIgnore
	 * @return
	 * @throws PostProcessingException 
	 */
	public static Element getNextIgnoringCertainElements(Element el, String[] elementsToIgnore) throws PostProcessingException {
		ParentNode parent = el.getParent();
		if (parent==null){
			throw new PostProcessingException("Node has no parent!");
		}
		int i = parent.indexOf(el);
		if (i+1 >= parent.getChildCount()) return null;
		Element next =(Element)parent.getChild(i+1);
		String elName =next.getLocalName();
		for (String namesToIgnore : elementsToIgnore) {
			if (elName.equals(namesToIgnore)){
				return getNextIgnoringCertainElements(next, elementsToIgnore);
			}
		}
		return next;
	}

	/**Makes a shallow copy of an element, copying the element
     * and the attribute, but no other child nodes.
     *
     * @param elem The element to copy.
     * @return The copied element.
     */
	public static Element shallowCopy(Element elem) {
		Element newElem = new Element(elem.getLocalName());
		for(int i=0;i<elem.getAttributeCount();i++) {
			newElem.addAttribute((Attribute)elem.getAttribute(i).copy());
		}
		return newElem;
	}

	/**
	 * Finds all descendant elements whose localname matches the given elementName
	 * Equivalent to an xquery of type .//elementName from the startingElement
	 * @param startingElement
	 * @param elementName
	 * @return
	 */
	public static List<Element> findDescendantElementsWithTagName(Element startingElement, String elementName) {
		List<Element> matchingElements = new ArrayList<Element>();
		LinkedList<Element> stack = new LinkedList<Element>();
		Elements children =startingElement.getChildElements();
		for (int i = children.size() -1; i >= 0; i--) {
			stack.add(children.get(i));
		}
		while (stack.size()>0){
			Element currentElement =stack.removeLast();
			if (currentElement.getLocalName().equals(elementName)){
				matchingElements.add(currentElement);
			}
			children =currentElement.getChildElements();
			for (int i = children.size() -1; i >= 0; i--) {
				Element child =children.get(i);
				stack.add(child);
			}
		}
		return matchingElements;
	}
	
	/**
	 * Finds all descendant elements whose localname matches one of the strings in elementNames
	 * Equivalent to an xquery of type .//elementName1|.//elementName2|.//elementName3 from the startingElement
	 * @param startingElement
	 * @param elementNames
	 * @return
	 */
	public static List<Element> findDescendantElementsWithTagNames(Element startingElement, String[] elementNames) {
		List<Element> matchingElements = new ArrayList<Element>();
		LinkedList<Element> stack = new LinkedList<Element>();
		Elements children =startingElement.getChildElements();
		for (int i = children.size() -1; i >= 0; i--) {
			stack.add(children.get(i));
		}
		while (stack.size()>0){
			Element currentElement =stack.removeLast();
			String currentElName=currentElement.getLocalName();
			for (String targetTagName : elementNames) {
				if (currentElName.equals(targetTagName)){
					matchingElements.add(currentElement);
					break;
				}
			}
			children =currentElement.getChildElements();
			for (int i = children.size() -1; i >= 0; i--) {
				Element child =children.get(i);
				stack.add(child);
			}
		}
		return matchingElements;
	}
	
	/**
	 * Finds all child elements whose localname matches one of the strings in elementNames
	 * Equivalent to an xquery of type ./elementName1|./elementName2|./elementName3 from the startingElement
	 * @param startingElement
	 * @param elementNames
	 * @return
	 */
	public static List<Element> findChildElementsWithTagNames(Element startingElement, String[] elementNames) {
		List<Element> matchingElements = new ArrayList<Element>();
		Elements children =startingElement.getChildElements();
		for (int i = 0; i < children.size(); i++) {
			Element child =children.get(i);
			String currentElName=child.getLocalName();
			for (String targetTagName : elementNames) {
				if (currentElName.equals(targetTagName)){
					matchingElements.add(child);
					break;
				}
			}
		}
		return matchingElements;
	}
	
	/**
	 * Finds all descendant elements whose localname matches the given elementName
	 * Additionally the element must have the specified attribute and the value of the attribute must be as specified
	 * Equivalent to an xquery of type .//elementName[@attribute="attributevalue"] from the startingElement
	 * @param startingElement
	 * @param elementName
	 * @return
	 */
	public static List<Element> findDescendantElementsWithTagNameAndAttribute(Element startingElement, String elementName, String attributeName, String attributeValue) {
		List<Element> matchingElements = new ArrayList<Element>();
		LinkedList<Element> stack = new LinkedList<Element>();
		Elements children =startingElement.getChildElements();
		for (int i = children.size() -1; i >= 0; i--) {
			stack.add(children.get(i));
		}
		while (stack.size()>0){
			Element currentElement =stack.removeLast();
			if (currentElement.getLocalName().equals(elementName)){
				if (currentElement.getAttribute(attributeName)!=null && currentElement.getAttributeValue(attributeName).equals(attributeValue)){
					matchingElements.add(currentElement);
				}
			}
			children =currentElement.getChildElements();
			for (int i = children.size() -1; i >= 0; i--) {
				Element child =children.get(i);
				stack.add(child);
			}
		}
		return matchingElements;
	}
	
	/**
	 * Finds all child elements whose localname matches the given elementName
	 * Additionally the element must have the specified attribute and the value of the attribute must be as specified
	 * Equivalent to an xquery of type ./elementName[@attribute="attributevalue"] from the startingElement
	 * @param startingElement
	 * @param elementName
	 * @return
	 */
	public static List<Element> findChildElementsWithTagNameAndAttribute(Element startingElement, String elementName, String attributeName, String attributeValue) {
		List<Element> matchingElements = new ArrayList<Element>();
		Elements children =startingElement.getChildElements();
		for (int i = 0; i < children.size(); i++) {
			Element child =children.get(i);
			if (child.getLocalName().equals(elementName)){
				if (child.getAttribute(attributeName)!=null && child.getAttributeValue(attributeName).equals(attributeValue)){
					matchingElements.add(child);
				}
			}
		}
		return matchingElements;
	}
	
	/**
	 * Finds and returns the count of descendant elements
	 * Equivalent to the size of the nodes array returned from the xquery .//* from the startingElement
	 * @param startingElement
	 * @param elementName
	 * @return
	 */
	public static int countDescendantElements(Element startingElement) {
		int count =0;
		LinkedList<Element> stack = new LinkedList<Element>();
		stack.add(startingElement);
		while (stack.size()>0){
			Element currentElement =stack.removeLast();
			Elements children =currentElement.getChildElements();
			for (int i = 0; i < children.size(); i++) {
				count++;
				stack.add(children.get(i));
			}
		}
		return count;
	}

	/**
	 * Searches in a depth-first manner for a non-suffix atom that has the target locant
	 * Returns either that atom or null if one cannot be found
	 * @param startingAtom
	 * @param targetLocant
	 * @return the matching atom or null
	 * @throws StructureBuildingException 
	 */
	public static Atom depthFirstSearchForNonSuffixAtomWithLocant(Atom startingAtom, String targetLocant) throws StructureBuildingException {
		LinkedList<Atom> stack = new LinkedList<Atom>();
		stack.add(startingAtom);
		Set<Atom> atomsVisited =new HashSet<Atom>();
		while (stack.size() > 0) {
			Atom currentAtom =stack.removeLast();
			atomsVisited.add(currentAtom);
			List<Atom> neighbours = currentAtom.getAtomNeighbours();
			for (Atom neighbour : neighbours) {
				if (atomsVisited.contains(neighbour)){//already visited
					continue;
				}
				List<String> locants = neighbour.getLocants();
				if (!neighbour.getType().equals("suffix")){//A main group atom, would expect to only find one except in something strange like succinimide
					for (String neighbourLocant : locants) {
						if (targetLocant.equals(neighbourLocant)){
							return neighbour;
						}
					}
					continue;
				}
				stack.add(neighbour);
			}
		}
		return null;
	}
}
