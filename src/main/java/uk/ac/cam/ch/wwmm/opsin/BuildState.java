package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import uk.ac.cam.ch.wwmm.opsin.PreProcessor.OpsinMode;

import nu.xom.Element;

/**
 * Used to pass the current mode, IDManager, FragmentManager and wordRule around as well as a mapping between the XML and fragments
 * Also contains which words appear to feature multiplicative nomenclature to aid in disambiguation
 *
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
	OpsinMode mode;
	HashSet<Element> multiplicativeNomenclaturePresent;
	boolean debug = false;

	BuildState(SMILESFragmentBuilder sBuilder, CMLFragmentBuilder cmlBuilder, OpsinMode mode) {
		idManager = new IDManager();
		fragManager = new FragmentManager(sBuilder, cmlBuilder, idManager);
		wordRule = null;
		xmlFragmentMap = new BiDirectionalHashMap();
		xmlSuffixMap = new HashMap<Element, ArrayList<Fragment>>();
		multiplicativeNomenclaturePresent = new HashSet<Element>();
		this.mode =mode;
	}
}
