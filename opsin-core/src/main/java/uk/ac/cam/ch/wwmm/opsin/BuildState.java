package uk.ac.cam.ch.wwmm.opsin;

import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;


import nu.xom.Element;

/**
 * Used to pass the current mode, IDManager, FragmentManager and wordRule around as well as a mapping between the XML and fragments
 *
 * @author dl387
 *
 */
class BuildState {

	final IDManager idManager;
	final FragmentManager fragManager;
	final BiDirectionalHashMap xmlFragmentMap;
	final HashMap<Element, List<Fragment>> xmlSuffixMap;
	final NameToStructureConfig n2sConfig;
	
	WordRule currentWordRule = null;
	
	private String warningMessage = null;
	
	String getWarningMessage() {
		return warningMessage;
	}

	void addWarningMessage(String warningMessage) {
		if (warningMessage == null){
			this.warningMessage = warningMessage;
		}
		else{
			this.warningMessage += ("\n" + warningMessage);
		}
	}

	/**
	 * Wrapper class for returning multiple objects
	 */
	final static class BiDirectionalHashMap implements Map<Element, Fragment>{
		final HashMap<Element, Fragment> xmlFragmentMap = new HashMap<Element, Fragment>();
		final HashMap<Fragment, Element> fragmentXmlMap = new HashMap<Fragment, Element>();
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
			for (Entry<? extends Element, ? extends Fragment> e : m.entrySet()) {
				fragmentXmlMap.put(e.getValue(), e.getKey());
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
		public Element getElement(Fragment key) {
			return fragmentXmlMap.get(key);
		}
	}

	BuildState(NameToStructureConfig n2sConfig, SMILESFragmentBuilder sBuilder) {
		this.n2sConfig = n2sConfig;
		idManager = new IDManager();
		fragManager = new FragmentManager(sBuilder, idManager);
		xmlFragmentMap = new BiDirectionalHashMap();
		xmlSuffixMap = new HashMap<Element, List<Fragment>>();
	}
}
