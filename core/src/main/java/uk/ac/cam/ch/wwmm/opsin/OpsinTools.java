package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.regex.Pattern;

import uk.ac.cam.ch.wwmm.opsin.ParseWord.WordType;


import nu.xom.Attribute;
import nu.xom.Element;
import nu.xom.Elements;
import nu.xom.Node;

import static uk.ac.cam.ch.wwmm.opsin.XmlDeclarations.*;

/**
 * A set of useful methods to assist OPSIN
 * @author dl387
 *
 */
class OpsinTools {
	private final static Pattern matchNumericLocant =Pattern.compile("\\d+[a-z]?'*");
	private final static char endOfSubstituent = '\u00e9';
	private final static char endOfMainGroup = '\u00e2';
	private final static char endOfFunctionalTerm = '\u00FB';
	
	/**
	 * Returns the next sibling suffix node which is not related to altering charge (ium/ide/id)
	 * @param currentEl
	 */
	static Element getNextNonChargeSuffix(Element currentEl) {
		Element next = (Element) XOMTools.getNextSibling(currentEl);
		while (next != null) {
			if (next.getLocalName().equals(SUFFIX_EL) && !CHARGE_TYPE_VAL.equals(next.getAttributeValue(TYPE_ATR))){
				return next;
			}
			next = (Element) XOMTools.getNextSibling(next);
		}
		return null;
	}

	/**
	 * Returns an arrayList of elements corresponding to the Elements given
	 * @param elements
	 * @return The new arrayList
	 */
	static ArrayList<Element> elementsToElementArrayList(Elements elements) {
		ArrayList<Element> elementList =new ArrayList<Element>(elements.size());
		for (int i = 0, n=elements.size(); i < n; i++) {
			elementList.add(elements.get(i));
		}
		return elementList;
	}

	/**
	 * Returns a new list containing the elements of list1 followed by list2
	 * @param list1
	 * @param list2
	 * @return The new list
	 */
	static ArrayList<Element> combineElementLists(List<Element> list1, List<Element> list2) {
		ArrayList<Element> elementList =new ArrayList<Element>(list1);
		elementList.addAll(list2);
		return elementList;
	}

	/**
	 * Returns the previous group. This group element need not be a sibling
	 * @param current: starting node
	 * @return
	 */
	static Node getPreviousGroup(Element current) {
	  if (current.getLocalName().equals(GROUP_EL)){//can start with a group or the sub/root the group is in
		  current=(Element)current.getParent();
	  }
	  Element parent = (Element) current.getParent();
	  if (parent == null || parent.getLocalName().equals(MOLECULE_EL)){
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
	  Elements groups =((Element)previous.getParent()).getChildElements(GROUP_EL);
	  if (groups.size()==0){
		  return getPreviousGroup(previous);
	  }
	  else{
		  return groups.get(groups.size()-1);//return last group if multiple exist e.g. fused ring
	  }
	}
	
	/**
	 * Returns the next group. This group element need not be a sibling
	 * @param current: starting node
	 * @return
	 */
	static Node getNextGroup(Element current) {
	  if (current.getLocalName().equals(GROUP_EL)){//can start with a group or the sub/root the group is in
		  current=(Element)current.getParent();
	  }
	  Element parent = (Element) current.getParent();
	  if (parent == null || parent.getLocalName().equals(MOLECULE_EL)){
		  return null;
	  }
	  int index = parent.indexOf(current);
	  if (index ==parent.getChildElements().size()-1) return getNextGroup(parent);//no group found
	  Element next =(Element) parent.getChild(index +1);
	  Elements children =next.getChildElements();
	  while (children.size()!=0){
		  next =children.get(0);
		  children =next.getChildElements();
	  }
	  Elements groups =((Element)next.getParent()).getChildElements(GROUP_EL);
	  if (groups.size()==0){
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
	 * @throws PostProcessingException
	 */
	static Element getParentWordRule(Element el) throws PostProcessingException {
		Element parent=(Element)el.getParent();
		while(parent !=null && !parent.getLocalName().equals(WORDRULE_EL)){
			parent =(Element)parent.getParent();
		}
		if (parent==null){
			throw new PostProcessingException("Cannot find enclosing wordRule element");
		}
		else{
			return parent;
		}
	}

	/**Makes a shallow copy of an element, copying the element
     * and the attributes, but no other child nodes.
     *
     * @param elem The element to copy.
     * @return The copied element.
     */
	static Element shallowCopy(Element elem) {
		Element newElem = new Element(elem.getLocalName());
		int attributeCount = elem.getAttributeCount();
		for(int i=0; i < attributeCount;i++) {
			newElem.addAttribute(new Attribute(elem.getAttribute(i)));
		}
		return newElem;
	}

	/**
	 * Searches in a depth-first manner for a non-suffix atom that has the target non element symbol locant
	 * Returns either that atom or null if one cannot be found
	 * @param startingAtom
	 * @param targetLocant
	 * @return the matching atom or null
	 */
	static Atom depthFirstSearchForNonSuffixAtomWithLocant(Atom startingAtom, String targetLocant) {
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
				List<String> locants = new ArrayList<String>(neighbour.getLocants());
				locants.removeAll(neighbour.getElementSymbolLocants());

				//A main group atom, would expect to only find one except in something strange like succinimide
				//The locants.size()>0 condition allows things like terephthalate to work which have an atom between the suffixes and main atoms that has no locant
				if (locants.size()>0 && !neighbour.getType().equals(SUFFIX_TYPE_VAL)){
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
	
	/**
	 * Searches in a depth-first manner for an atom with a numeric locant
	 * Returns either that atom or null if one cannot be found
	 * @param startingAtom
	 * @return the matching atom or null
	 */
	static Atom depthFirstSearchForAtomWithNumericLocant(Atom startingAtom){
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
				for (String neighbourLocant : locants) {
					if (matchNumericLocant.matcher(neighbourLocant).matches()){
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
		Character finalAnnotation = annotations.get(annotations.size() -1);
		if (finalAnnotation.equals(endOfMainGroup)){
			return WordType.full;
		}
		else if (finalAnnotation.equals(endOfSubstituent)){
			return WordType.substituent;
		}
		else if (finalAnnotation.equals(endOfFunctionalTerm)){
			return WordType.functionalTerm;
		}
		else{
			throw new ParsingException("OPSIN bug: Unable to determine word type!");
		}
		
	}
}
