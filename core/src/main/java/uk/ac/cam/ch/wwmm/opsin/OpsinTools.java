package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;


import nu.xom.Attribute;
import nu.xom.Element;
import nu.xom.Elements;
import nu.xom.Node;

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
	public static ArrayList<Element> combineElementLists(List<Element> list1, List<Element> list2) {
		ArrayList<Element> elementList =new ArrayList<Element>(list1);
		elementList.addAll(list2);
		return elementList;
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
	 * Finds the wordRule element that encloses the given element.
	 * Returns the wordRule element or throws an exception
	 * @param Element el
	 * @return wordRule Element
	 * @throws PostProcessingException
	 */
	public static Element getParentWordRule(Element el) throws PostProcessingException {
		Element parent=(Element)el.getParent();
		while(parent !=null && !parent.getLocalName().equals("wordRule")){
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

				//A main group atom, would expect to only find one except in something strange like succinimide
				//The locants.size()>0 condition allows things like terephthalate to work which have an atom between the suffixes and main atoms that has no locant
				if (locants.size()>0 && !neighbour.getType().equals("suffix")){
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
