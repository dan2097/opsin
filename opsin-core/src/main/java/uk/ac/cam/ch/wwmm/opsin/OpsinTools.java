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
		Element next = XOMTools.getNextSibling(currentEl);
		while (next != null) {
			if (next.getLocalName().equals(SUFFIX_EL) && !CHARGE_TYPE_VAL.equals(next.getAttributeValue(TYPE_ATR))){
				return next;
			}
			next = XOMTools.getNextSibling(next);
		}
		return null;
	}

	/**
	 * Returns an arrayList of elements corresponding to the Elements given
	 * @param elements
	 * @return The new arrayList
	 */
	static ArrayList<Element> elementsToElementArrayList(List<Element> elements) {
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
	static Element getPreviousGroup(Element current) {
	  if (current.getLocalName().equals(GROUP_EL)){//can start with a group or the sub/root the group is in
		  current = current.getParent();
	  }
	  Element parent = current.getParent();
	  if (parent == null || parent.getLocalName().equals(WORDRULE_EL)){
		  return null;
	  }
	  int index = parent.indexOf(current);
	  if (index ==0) return getPreviousGroup(parent);//no group found
	  Element previous = parent.getChild(index-1);
	  List<Element> children =previous.getChildElements();
	  while (children.size()!=0){
		  previous =children.get(children.size()-1);
		  children =previous.getChildElements();
	  }
	  List<Element> groups = previous.getParent().getChildElements(GROUP_EL);
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
	static Element getNextGroup(Element current) {
	  if (current.getLocalName().equals(GROUP_EL)){//can start with a group or the sub/root the group is in
		  current = current.getParent();
	  }
	  Element parent = current.getParent();
	  if (parent == null || parent.getLocalName().equals(MOLECULE_EL)){
		  return null;
	  }
	  int index = parent.indexOf(current);
	  if (index ==parent.getChildElements().size()-1) return getNextGroup(parent);//no group found
	  Element next = parent.getChild(index +1);
	  List<Element> children =next.getChildElements();
	  while (children.size()!=0){
		  next =children.get(0);
		  children =next.getChildElements();
	  }
	  List<Element> groups = next.getParent().getChildElements(GROUP_EL);
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
	 */
	static Element getParentWordRule(Element el) {
		Element parent = el.getParent();
		while(parent != null && !parent.getLocalName().equals(WORDRULE_EL)){
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
		Deque<Atom> stack = new ArrayDeque<Atom>();
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
		Character finalAnnotation = annotations.get(annotations.size() -1);
		if (finalAnnotation.equals(END_OF_MAINGROUP)){
			return WordType.full;
		}
		else if (finalAnnotation.equals(END_OF_SUBSTITUENT)){
			return WordType.substituent;
		}
		else if (finalAnnotation.equals(END_OF_FUNCTIONALTERM)){
			return WordType.functionalTerm;
		}
		else{
			throw new ParsingException("OPSIN bug: Unable to determine word type!");
		}
		
	}
}
