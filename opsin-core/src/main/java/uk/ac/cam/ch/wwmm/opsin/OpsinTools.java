package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Deque;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.regex.Pattern;

import static uk.ac.cam.ch.wwmm.opsin.XmlDeclarations.*;

/**
 * A set of useful methods and constants to assist OPSIN
 * @author dl387
 *
 */
class OpsinTools {
	static final Pattern MATCH_COLON = Pattern.compile(":");
	static final Pattern MATCH_COLONORSEMICOLON = Pattern.compile("[:;]");
	static final Pattern MATCH_COMMA = Pattern.compile(",");
	static final Pattern MATCH_DASH = Pattern.compile("-");
	static final Pattern MATCH_SEMICOLON = Pattern.compile(";");
	static final Pattern MATCH_SLASH = Pattern.compile("/");
	static final Pattern MATCH_SPACE =Pattern.compile(" ");
	
	static final Pattern MATCH_AMINOACID_STYLE_LOCANT =Pattern.compile("([A-Z][a-z]?)('*)((\\d+[a-z]?|alpha|beta|gamma|delta|epsilon|zeta|eta|omega)'*)");
	static final Pattern MATCH_ELEMENT_SYMBOL =Pattern.compile("[A-Z][a-z]?");
	static final Pattern MATCH_ELEMENT_SYMBOL_LOCANT =Pattern.compile("[A-Z][a-z]?'*");
	static final Pattern MATCH_NUMERIC_LOCANT =Pattern.compile("(\\d+)[a-z]?'*");
	static final char END_OF_SUBSTITUENT = '\u00e9';
	static final char END_OF_MAINGROUP = '\u00e2';
	static final char END_OF_FUNCTIONALTERM = '\u00FB';
	
	/**
	 * Returns the next sibling suffix node which is not related to altering charge (ium/ide/id)
	 * @param currentEl
	 */
	static Element getNextNonChargeSuffix(Element currentEl) {
		Element next = OpsinTools.getNextSibling(currentEl);
		while (next != null) {
			if (next.getName().equals(SUFFIX_EL) && !CHARGE_TYPE_VAL.equals(next.getAttributeValue(TYPE_ATR))){
				return next;
			}
			next = OpsinTools.getNextSibling(next);
		}
		return null;
	}

	/**
	 * Returns a new list containing the elements of list1 followed by list2
	 * @param list1
	 * @param list2
	 * @return The new list
	 */
	static List<Element> combineElementLists(List<Element> list1, List<Element> list2) {
		List<Element> elementList = new ArrayList<Element>(list1);
		elementList.addAll(list2);
		return elementList;
	}

	/**
	 * Returns the previous group. This group element need not be a sibling
	 * @param current: starting element
	 * @return
	 */
	static Element getPreviousGroup(Element current) {
	  if (current.getName().equals(GROUP_EL)) {//can start with a group or the sub/root the group is in
		  current = current.getParent();
	  }
	  Element parent = current.getParent();
	  if (parent == null || parent.getName().equals(WORDRULE_EL)) {
		  return null;
	  }
	  int index = parent.indexOf(current);
	  if (index ==0) {
		  return getPreviousGroup(parent);//no group found
	  }
	  Element previous = parent.getChild(index - 1);
	  List<Element> children = previous.getChildElements();
	  while (children.size() != 0){
		  previous = children.get(children.size() - 1);
		  children = previous.getChildElements();
	  }
	  List<Element> groups = previous.getParent().getChildElements(GROUP_EL);
	  if (groups.size() == 0){
		  return getPreviousGroup(previous);
	  }
	  else{
		  return groups.get(groups.size() - 1);//return last group if multiple exist e.g. fused ring
	  }
	}
	
	/**
	 * Returns the next group. This group element need not be a sibling
	 * @param current: starting element
	 * @return
	 */
	static Element getNextGroup(Element current) {
	  if (current.getName().equals(GROUP_EL)) {//can start with a group or the sub/root the group is in
		  current = current.getParent();
	  }
	  Element parent = current.getParent();
	  if (parent == null || parent.getName().equals(MOLECULE_EL)) {
		  return null;
	  }
	  int index = parent.indexOf(current);
	  if (index == parent.getChildCount() - 1) {
		  return getNextGroup(parent);//no group found
	  }
	  Element next = parent.getChild(index + 1);
	  List<Element> children = next.getChildElements();
	  while (children.size() != 0){
		  next = children.get(0);
		  children = next.getChildElements();
	  }
	  List<Element> groups = next.getParent().getChildElements(GROUP_EL);
	  if (groups.size() == 0){
		  return getNextGroup(next);
	  }
	  else{
		  return groups.get(0);//return first group if multiple exist e.g. fused ring
	  }
	}

	/**
	 * Finds the wordRule element that encloses the given element.
	 * Returns the wordRule element or throws an exception
	 * @param el
	 * @return wordRule Element
	 */
	static Element getParentWordRule(Element el) {
		Element parent = el.getParent();
		while(parent != null && !parent.getName().equals(WORDRULE_EL)){
			parent = parent.getParent();
		}
		if (parent == null){
			throw new RuntimeException("Cannot find enclosing wordRule element");
		}
		else{
			return parent;
		}
	}

	/**
	 * Searches in a depth-first manner for a non-suffix atom that has the target non element symbol locant
	 * Returns either that atom or null if one cannot be found
	 * @param startingAtom
	 * @param targetLocant
	 * @return the matching atom or null
	 */
	static Atom depthFirstSearchForNonSuffixAtomWithLocant(Atom startingAtom, String targetLocant) {
		Deque<Atom> stack = new ArrayDeque<Atom>();
		stack.add(startingAtom);
		Set<Atom> atomsVisited = new HashSet<Atom>();
		while (stack.size() > 0) {
			Atom currentAtom = stack.removeLast();
			atomsVisited.add(currentAtom);
			List<Atom> neighbours = currentAtom.getAtomNeighbours();
			for (Atom neighbour : neighbours) {
				if (atomsVisited.contains(neighbour)){//already visited
					continue;
				}
				List<String> locants = new ArrayList<String>(neighbour.getLocants());
				locants.removeAll(neighbour.getElementSymbolLocants());

				//A main group atom, would expect to only find one except in something strange like succinimide
				//The locants.size() > 0 condition allows things like terephthalate to work which have an atom between the suffixes and main atoms that has no locant
				if (locants.size() > 0 && !neighbour.getType().equals(SUFFIX_TYPE_VAL)){
					if (locants.contains(targetLocant)){
						return neighbour;
					}
					continue;
				}
				stack.add(neighbour);
			}
		}
		return null;
	}
	
	/**
	 * Searches in a depth-first manner for an atom with a numeric locant
	 * Returns either that atom or null if one cannot be found
	 * @param startingAtom
	 * @return the matching atom or null
	 */
	static Atom depthFirstSearchForAtomWithNumericLocant(Atom startingAtom){
		Deque<Atom> stack = new ArrayDeque<Atom>();
		stack.add(startingAtom);
		Set<Atom> atomsVisited = new HashSet<Atom>();
		while (stack.size() > 0) {
			Atom currentAtom = stack.removeLast();
			atomsVisited.add(currentAtom);
			List<Atom> neighbours = currentAtom.getAtomNeighbours();
			for (Atom neighbour : neighbours) {
				if (atomsVisited.contains(neighbour)){//already visited
					continue;
				}
				List<String> locants = neighbour.getLocants();
				for (String neighbourLocant : locants) {
					if (MATCH_NUMERIC_LOCANT.matcher(neighbourLocant).matches()){
						return neighbour;
					}
				}
				stack.add(neighbour);
			}
		}
		return null;
	}
	
	/**
	 * Given a list of annotations returns the word type as indicated by the final annotation of the list
	 * @param annotations
	 * @return WordType
	 * @throws ParsingException 
	 */
	static WordType determineWordType(List<Character> annotations) throws ParsingException {
		char finalAnnotation = annotations.get(annotations.size() - 1);
		if (finalAnnotation == END_OF_MAINGROUP) {
			return WordType.full;
		}
		else if (finalAnnotation == END_OF_SUBSTITUENT) {
			return WordType.substituent;
		}
		else if (finalAnnotation == END_OF_FUNCTIONALTERM) {
			return WordType.functionalTerm;
		}
		else{
			throw new ParsingException("OPSIN bug: Unable to determine word type!");
		}
		
	}

	/**Gets the next sibling of a given element.
	*
	* @param element The reference element.
	* @return The next Sibling, or null.
	*/
	static Element getNextSibling(Element element) {
		Element parent = element.getParent();
		int i = parent.indexOf(element);
		if (i + 1 >= parent.getChildCount()) {
			return null;
		}
		return parent.getChild(i + 1);
	}

	/**Gets the first next sibling of a given element whose element name matches the given string.
	*
	* @param current The reference element.
	* @param elName The element name to look for
	* @return The matched next Sibling, or null.
	*/
	static Element getNextSibling(Element current, String elName) {
		Element next = getNextSibling(current);
		while (next != null) {
			if (next.getName().equals(elName)){
				return next;
			}
			next = getNextSibling(next);
		}
		return null;
	}

	/**Gets the previous sibling of a given element.
	*
	* @param element The reference element.
	* @return The previous Sibling, or null.
	*/
	static Element getPreviousSibling(Element element) {
		Element parent = element.getParent();
		int i = parent.indexOf(element);
		if (i == 0) {
			return null;
		}
		return parent.getChild(i - 1);
	}

	/**Gets the first previous sibling of a given element whose element name matches the given string.
	*
	* @param current The reference element.
	* @param elName The element name of a element to look for
	* @return The matched previous Sibling, or null.
	*/
	static Element getPreviousSibling(Element current, String elName) {;
		Element prev = getPreviousSibling(current);
		while (prev != null) {
			if (prev.getName().equals(elName)){
				return prev;
			}
			prev = getPreviousSibling(prev);
		}
		return null;
	}

	/**Inserts a element so that it occurs before a reference element. The new element
	 * must not currently have a parent.
	 *
	 * @param element The reference element.
	 * @param newElement The new element to insert.
	 */
	static void insertBefore(Element element, Element newElement) {
		Element parent = element.getParent();
		int i = parent.indexOf(element);
		parent.insertChild(newElement, i);
	}

	/**Inserts an element so that it occurs after a reference element. The new element
	 * must not currently have a parent.
	 *
	 * @param element The reference element.
	 * @param neweElement The new element to insert.
	 */
	static void insertAfter(Element element, Element neweElement) {
		Element parent = element.getParent();
		int i = parent.indexOf(element);
		parent.insertChild(neweElement, i + 1);
	}

	/**
	 * Gets the next element. This element need not be a sibling
	 * @param element: starting element
	 * @return
	 */
	static Element getNext(Element element) {
		Element parent = element.getParent();
		if (parent == null || parent.getName().equals(XmlDeclarations.MOLECULE_EL)){
			return null;
		}
		int index = parent.indexOf(element);
		if (index + 1 >= parent.getChildCount()) {
			return getNext(parent);//reached end of element
		}
		Element next = parent.getChild(index + 1);
		while (next.getChildCount() > 0){
			next = next.getChild(0);
		}
		return next;
	}

	/**
	 * Gets the previous element. This element need not be a sibling
	 * @param element: starting element
	 * @return
	 */
	static Element getPrevious(Element element) {
		Element parent = element.getParent();
		if (parent == null || parent.getName().equals(XmlDeclarations.MOLECULE_EL)){
			return null;
		}
		int index = parent.indexOf(element);
		if (index == 0) {
			return getPrevious(parent);//reached beginning of element
		}
		Element previous = parent.getChild(index - 1);
		while (previous.getChildCount() > 0){
			previous = previous.getChild(previous.getChildCount() - 1);
		}
		return previous;
	}

	/**
	 * Returns a list containing sibling elements with the given element name after the given element.
	 * These elements need not be continuous
	 * @param currentElem: the element to look for following siblings of
	 * @param elName: the name of the elements desired
	 * @return
	 */
	static List<Element> getNextSiblingsOfType(Element currentElem, String elName) {
		List<Element> laterSiblingElementsOfType = new ArrayList<Element>();
		Element parent = currentElem.getParent();
		if (parent == null){
			return laterSiblingElementsOfType;
		}
		int indexOfCurrentElem = parent.indexOf(currentElem);
		for (int i = indexOfCurrentElem + 1; i < parent.getChildCount(); i++) {
			Element child = parent.getChild(i);
			if (child.getName().equals(elName)) {
				laterSiblingElementsOfType.add(child);
			}
		}
		return laterSiblingElementsOfType;
	}

	/**
	 * Returns a list containing sibling elements with the given element name after the given element.
	 * @param currentElem: the element to look for following siblings of
	 * @param elName: the name of the elements desired
	 * @return
	 */
	static List<Element> getNextAdjacentSiblingsOfType(Element currentElem, String elName) {
		List<Element> siblingElementsOfType = new ArrayList<Element>();
		Element parent = currentElem.getParent();
		if (parent == null){
			return siblingElementsOfType;
		}
		Element nextSibling = getNextSibling(currentElem);
		while (nextSibling != null && nextSibling.getName().equals(elName)){
			siblingElementsOfType.add(nextSibling);
			nextSibling = getNextSibling(nextSibling);
		}
		return siblingElementsOfType;
	}

	/**
	 * Returns a list containing sibling elements with the given element names after the given element.
	 * These elements need not be continuous and are returned in the order encountered
	 * @param currentElem: the element to look for following siblings of
	 * @param types: An array of the names of the elements desired
	 * @return
	 */
	static List<Element> getNextSiblingsOfTypes(Element currentElem, String[] elNames){
		List<Element> laterSiblingElementsOfTypes = new ArrayList<Element>();
		currentElem = getNextSibling(currentElem);
		while (currentElem != null){
			String name = currentElem.getName();
			for (String elName : elNames) {
				if (name.equals(elName)){
					laterSiblingElementsOfTypes.add(currentElem);
					break;
				}
			}
			currentElem = getNextSibling(currentElem);
		}
		return laterSiblingElementsOfTypes;
	}

	/**
	 * Returns a list containing sibling elements with the given element name before the given element.
	 * These elements need not be continuous
	 * @param currentElem: the element to look for previous siblings of
	 * @param elName: the name of the elements desired
	 * @return
	 */
	static List<Element> getPreviousSiblingsOfType(Element currentElem, String elName) {
		List<Element> earlierSiblingElementsOfType = new ArrayList<Element>();
		Element parent = currentElem.getParent();
		if (parent == null){
			return earlierSiblingElementsOfType;
		}
		int indexOfCurrentElem = parent.indexOf(currentElem);
		for (int i = 0; i < indexOfCurrentElem - 1; i++) {
			Element child = parent.getChild(i);
			if (child.getName().equals(elName)) {
				earlierSiblingElementsOfType.add(child);
			}
		}
		return earlierSiblingElementsOfType;
	}

	/**
	 * Gets the next sibling element of the given element. If this element's name is within the elementsToIgnore array this is repeated
	 * If no appropriate element can be found null is returned
	 * @param startingEl
	 * @param elNamesToIgnore
	 * @return
	 */
	static Element getNextSiblingIgnoringCertainElements(Element startingEl, String[] elNamesToIgnore){
		Element parent = startingEl.getParent();
		if (parent == null){
			return null;
		}
		int i = parent.indexOf(startingEl);
		if (i + 1 >= parent.getChildCount()) {
			return null;
		}
		Element next = parent.getChild(i + 1);
		String elName = next.getName();
		for (String namesToIgnore : elNamesToIgnore) {
			if (elName.equals(namesToIgnore)){
				return getNextSiblingIgnoringCertainElements(next, elNamesToIgnore);
			}
		}
		return next;
	}

	/**
	 * Gets the previous sibling element of the given element. If this element's name is within the elementsToIgnore array this is repeated
	 * If no appropriate element can be found null is returned
	 * @param startingEl
	 * @param elNamesToIgnore
	 * @return
	 */
	static Element getPreviousSiblingIgnoringCertainElements(Element startingEl, String[] elNamesToIgnore){
		Element parent = startingEl.getParent();
		if (parent == null){
			return null;
		}
		int i = parent.indexOf(startingEl);
		if (i == 0) {
			return null;
		}
		Element previous = parent.getChild(i - 1);
		String elName = previous.getName();
		for (String namesToIgnore : elNamesToIgnore) {
			if (elName.equals(namesToIgnore)){
				return getPreviousSiblingIgnoringCertainElements(previous, elNamesToIgnore);
			}
		}
		return previous;
	}

	/**
	 * Finds all descendant elements whose name  matches the given element name
	 * @param startingElement
	 * @param elementName
	 * @return
	 */
	static List<Element> getDescendantElementsWithTagName(Element startingElement, String elementName) {
		List<Element> matchingElements = new ArrayList<Element>();
		Deque<Element> stack = new ArrayDeque<Element>();
		List<Element> children = startingElement.getChildElements();
		for (int i = children.size() -1; i >= 0; i--) {
			stack.add(children.get(i));
		}
		while (stack.size() > 0){
			Element currentElement = stack.removeLast();
			if (currentElement.getName().equals(elementName)){
				matchingElements.add(currentElement);
			}
			children = currentElement.getChildElements();
			for (int i = children.size() -1; i >= 0; i--) {
				Element child = children.get(i);
				stack.add(child);
			}
		}
		return matchingElements;
	}

	/**
	 * Finds all descendant elements whose element name matches one of the strings in elementNames
	 * @param startingElement
	 * @param elementNames
	 * @return
	 */
	static List<Element> getDescendantElementsWithTagNames(Element startingElement, String[] elementNames) {
		List<Element> matchingElements = new ArrayList<Element>();
		Deque<Element> stack = new ArrayDeque<Element>();
		List<Element> children =startingElement.getChildElements();
		for (int i = children.size() -1; i >= 0; i--) {
			stack.add(children.get(i));
		}
		while (stack.size()>0){
			Element currentElement = stack.removeLast();
			String currentElName = currentElement.getName();
			for (String targetTagName : elementNames) {
				if (currentElName.equals(targetTagName)){
					matchingElements.add(currentElement);
					break;
				}
			}
			children = currentElement.getChildElements();
			for (int i = children.size() -1; i >= 0; i--) {
				Element child =children.get(i);
				stack.add(child);
			}
		}
		return matchingElements;
	}

	/**
	 * Finds all child elements whose element name matches one of the strings in elementNames
	 * @param startingElement
	 * @param elementNames
	 * @return
	 */
	static List<Element> getChildElementsWithTagNames(Element startingElement, String[] elementNames) {
		List<Element> matchingElements = new ArrayList<Element>();
		List<Element> children = startingElement.getChildElements();
		for (Element  child : children) {
			String currentElName = child.getName();
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
	 * Finds all child elements whose localname matches one of the strings in elementNames
	 * Equivalent to an xpath of type ./*[local-name() = 'elementName'] from the startingElement
	 * This is equivalent to XOM's getChildElements(String) other than returning a list
	 * @param startingElement
	 * @param elementName
	 * @return
	 */
	static List<Element> getChildElementsWithTagName(Element startingElement, String elementName) {
		//TODO remove this method as it's redundant
		List<Element> matchingElements = new ArrayList<Element>();
		List<Element> children =startingElement.getChildElements();
		int childCount = children.size();
		for (int i = 0; i < childCount; i++) {
			Element child =children.get(i);
			String currentElName=child.getName();
			if (currentElName.equals(elementName)){
				matchingElements.add(child);
			}
		}
		return matchingElements;
	}

	/**
	 * Finds all descendant elements whose element name matches the given elementName
	 * Additionally the element must have the specified attribute and the value of the attribute must be as specified
	 * @param startingElement
	 * @param elementName
	 * @param attributeName
	 * @param attributeValue
	 * @return
	 */
	static List<Element> getDescendantElementsWithTagNameAndAttribute(Element startingElement, String elementName, String attributeName, String attributeValue) {
		List<Element> matchingElements = new ArrayList<Element>();
		Deque<Element> stack = new ArrayDeque<Element>();
		List<Element> children =startingElement.getChildElements();
		for (int i = children.size() - 1; i >= 0; i--) {
			stack.add(children.get(i));
		}
		while (stack.size() > 0){
			Element currentElement =stack.removeLast();
			if (currentElement.getName().equals(elementName)){
				if (attributeValue.equals(currentElement.getAttributeValue(attributeName))){
					matchingElements.add(currentElement);
				}
			}
			children =currentElement.getChildElements();
			for (int i = children.size() - 1; i >= 0; i--) {
				Element child =children.get(i);
				stack.add(child);
			}
		}
		return matchingElements;
	}

	/**
	 * Finds all child elements whose element name matches the given elementName
	 * Additionally the element must have the specified attribute and the value of the attribute must be as specified
	 * @param startingElement
	 * @param elementName
	 * @return
	 */
	static List<Element> getChildElementsWithTagNameAndAttribute(Element startingElement, String elementName, String attributeName, String attributeValue) {
		List<Element> matchingElements = new ArrayList<Element>();
		List<Element> children = startingElement.getChildElements();
		for (Element child : children) {
			if (child.getName().equals(elementName)){
				if (attributeValue.equals(child.getAttributeValue(attributeName))){
					matchingElements.add(child);
				}
			}
		}
		return matchingElements;
	}

	/**
	 * Finds and returns the number of elements and the number of elements with no children, that are descendants of the startingElement
	 * The 0th position of the returned array is the total number of elements
	 * The 1st position is the number of child less elements
	 * @param startingElement
	 * @return
	 */
	static int[] countNumberOfElementsAndNumberOfChildLessElements(Element startingElement) {
		int[] counts = new int[2];
		Deque<Element> stack = new ArrayDeque<Element>();
		stack.add(startingElement);
		while (stack.size() > 0){
			Element currentElement = stack.removeLast();
			List<Element> children = currentElement.getChildElements();
			int childCount = children.size();
			if (childCount == 0){
				counts[1]++;
			}
			else{
				for (int i = 0; i < childCount; i++) {
					counts[0]++;
					stack.add(children.get(i));
				}
			}
		}
		return counts;
	}

	/**
	 * Find all the later siblings of startingElement up until there is no more siblings or an
	 * element with the given element name  is reached (exclusive of that element)
	 * @param startingEl
	 * @param elName
	 * @return
	 */
	static List<Element> getSiblingsUpToElementWithTagName(Element startingEl, String elName) {
		List<Element> laterSiblings = new ArrayList<Element>();
		Element nextEl = getNextSibling(startingEl);
		while (nextEl != null && !nextEl.getName().equals(elName)){
			laterSiblings.add(nextEl);
			nextEl = getNextSibling(nextEl);
		}
		return laterSiblings;
	}
}
