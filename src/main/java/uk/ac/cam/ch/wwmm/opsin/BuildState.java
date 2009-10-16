package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import nu.xom.Element;

/**
 * Used to pass the current IDManager, FragmentManager and wordRule around as well as a mapping between the XML and fragments
 *
 * Also included is the first element which has multiple outIDs as detected by the preStructureBuilder in a particular word.
 * Usually no groups have multiple outIDs so this hash will have entries
 * If it does have entries this indicates that either the structure is a radical or that it is an example of multiplicative nomenclature
 * Multiplicative nomenclature is typically resolved left to right rather than right to left so from that element onwards in that word
 * the name will be resolved left to right in the structureBuilder. Once it reaches the end of the word it will resolve anything before
 * the detected multi radical in the conventional right to left manner.
 * @author dl387
 *
 */
public class BuildState {
	/**
	 * Wrapper class for returning multiple objects
	 */
	final class BiDirectionalHashMap implements Map<Element, Fragment>{
		HashMap<Element, Fragment> xmlFragmentMap = new HashMap<Element, Fragment>();
		HashMap<Fragment, Element> fragmentXmlMap = new HashMap<Fragment, Element>();
		public void clear() {
			xmlFragmentMap.clear();
			fragmentXmlMap.clear();
		}
		public boolean containsKey(Object key) {
			return xmlFragmentMap.containsKey(key);
		}
		public boolean containsValue(Object value) {
			return xmlFragmentMap.containsValue(value);
		}
		public Set<java.util.Map.Entry<Element, Fragment>> entrySet() {
			return xmlFragmentMap.entrySet();
		}
		public Fragment get(Object key) {
			return xmlFragmentMap.get(key);
		}
		public boolean isEmpty() {
			return xmlFragmentMap.isEmpty();
		}
		public Set<Element> keySet() {
			return xmlFragmentMap.keySet();
		}
		public Fragment put(Element key, Fragment value) {
			fragmentXmlMap.put(value, key);
			return xmlFragmentMap.put(key, value);
		}
		public void putAll(Map<? extends Element, ? extends Fragment> m) {
			for (Element el : m.keySet()) {
				fragmentXmlMap.put(m.get(el), el);
			}
			xmlFragmentMap.putAll(m);
		}
		public Fragment remove(Object key) {
			Fragment f =xmlFragmentMap.remove(key);
			fragmentXmlMap.remove(f);
			return f;
		}
		public int size() {
			return xmlFragmentMap.size();
		}
		public Collection<Fragment> values() {
			return xmlFragmentMap.values();
		}
		public Element getElement(Object key) {
			return fragmentXmlMap.get(key);
		}
	}

	IDManager idManager;
	FragmentManager fragManager;
	String wordRule;
	BiDirectionalHashMap xmlFragmentMap;
	HashMap<Element, ArrayList<Fragment>> xmlSuffixMap;
	HashMap<Element, Element> firstMultiRadical;//hash of word element against substituent/root element with multi radical

	BuildState(SMILESFragmentBuilder sBuilder, CMLFragmentBuilder cmlBuilder) {
		idManager = new IDManager();
		fragManager = new FragmentManager(sBuilder, cmlBuilder, idManager);
		wordRule = null;
		xmlFragmentMap = new BiDirectionalHashMap();
		xmlSuffixMap = new HashMap<Element, ArrayList<Fragment>>();
		firstMultiRadical = new HashMap<Element, Element>();
	}
}
