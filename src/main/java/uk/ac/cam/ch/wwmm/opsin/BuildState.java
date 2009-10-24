package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import uk.ac.cam.ch.wwmm.opsin.PreProcessor.OpsinMode;

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
	final class BiDirectionalHashMap implements Map<Element, List<Fragment>>{
		HashMap<Element, List<Fragment>> xmlFragmentMap = new HashMap<Element, List<Fragment>>();
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
		public Set<java.util.Map.Entry<Element, List<Fragment>>> entrySet() {
			return xmlFragmentMap.entrySet();
		}
		public List<Fragment> get(Object key) {
			return xmlFragmentMap.get(key);
		}
		public boolean isEmpty() {
			return xmlFragmentMap.isEmpty();
		}
		public Set<Element> keySet() {
			return xmlFragmentMap.keySet();
		}
		public List<Fragment> put(Element key, Fragment value) {
			fragmentXmlMap.put(value, key);
			List<Fragment> l = new ArrayList<Fragment>();
			l.add(value);
			return xmlFragmentMap.put(key, l);
		}
		public List<Fragment> put(Element key, List<Fragment> value) {
			for (Fragment frag : value) {
				fragmentXmlMap.put(frag, key);
			}
			return xmlFragmentMap.put(key, value);
		}
		public void putAll(Map<? extends Element, ? extends List<Fragment>> m) {
			for (Element el : m.keySet()) {
				List<Fragment> frags = m.get(el);
				for (Fragment frag : frags) {
					fragmentXmlMap.put(frag, el);
				}
			}
			xmlFragmentMap.putAll(m);
		}
		public List<Fragment> remove(Object key) {
			List<Fragment> f =xmlFragmentMap.remove(key);
			fragmentXmlMap.remove(f);
			return f;
		}
		public int size() {
			return xmlFragmentMap.size();
		}
		public Collection<List<Fragment>> values() {
			return xmlFragmentMap.values();
		}
		public Element getElement(Object key) {
			return fragmentXmlMap.get(key);
		}
		/**
		 * Convenience method to return the first fragment (or null if no fragment was associated)
		 * @param group
		 * @return
		 */
		public Fragment getFirstFragment(Element group) {
			List<Fragment> frags = get(group);
			if (frags ==null){
				return null;
			}
			if (frags.size() >=1){
				return frags.get(0);
			}
			else{
				return null;
			}
		}
	}

	IDManager idManager;
	FragmentManager fragManager;
	String wordRule;
	BiDirectionalHashMap xmlFragmentMap;
	HashMap<Element, ArrayList<Fragment>> xmlSuffixMap;
	OpsinMode mode;

	BuildState(SMILESFragmentBuilder sBuilder, CMLFragmentBuilder cmlBuilder, OpsinMode mode) {
		idManager = new IDManager();
		fragManager = new FragmentManager(sBuilder, cmlBuilder, idManager);
		wordRule = null;
		xmlFragmentMap = new BiDirectionalHashMap();
		xmlSuffixMap = new HashMap<Element, ArrayList<Fragment>>();
		this.mode =mode;
	}
}
