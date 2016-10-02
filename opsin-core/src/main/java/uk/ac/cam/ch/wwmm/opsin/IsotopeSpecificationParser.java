package uk.ac.cam.ch.wwmm.opsin;

import static uk.ac.cam.ch.wwmm.opsin.OpsinTools.MATCH_COMMA;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

class IsotopeSpecificationParser {
	
	private static final Pattern matchBoughtonIsotope =Pattern.compile("(?:-([^,]+(?:,[^,]+)*))?-d(\\d+)?");

	static class IsotopeSpecification {
		private final ChemEl chemEl;
		private final int isotope;
		private final int multiplier;
		private final String[] locants;
		
		IsotopeSpecification(ChemEl chemEl, int isotope, int multiplier, String[] locants) {
			this.chemEl = chemEl;
			this.isotope = isotope;
			this.multiplier = multiplier;
			this.locants = locants;
		}

		ChemEl getChemEl() {
			return chemEl;
		}

		int getIsotope() {
			return isotope;
		}

		int getMultiplier() {
			return multiplier;
		}

		String[] getLocants() {
			return locants;
		}
	}
	
	static IsotopeSpecification parseIsotopeSpecification(Element isotopeSpecification) throws StructureBuildingException {
		if (!XmlDeclarations.BOUGHTONSYSTEM_TYPE_VAL.equals(isotopeSpecification.getAttributeValue(XmlDeclarations.TYPE_ATR))) {
			throw new RuntimeException("Unsupported isotope specification syntax");
		}
		String val = isotopeSpecification.getValue();
		Matcher m = matchBoughtonIsotope.matcher(val);
		if (!m.matches()) {
			throw new RuntimeException("Malformed isotope specification: " + val);
		}
		ChemEl chemEl = ChemEl.H;
		int isotope = 2;
		
		int multiplier = 1;
		String multiplierStr = m.group(2);
		if (multiplierStr != null) {
			multiplier = Integer.parseInt(multiplierStr);
		}
		
		String locantsStr = m.group(1);
		String[] locants = null;
		if(locantsStr != null) {
			locants = MATCH_COMMA.split(locantsStr);
			if (locants.length != multiplier) {
				throw new StructureBuildingException("Mismatch between number of locants: " + locants.length + " and number of hydrogen isotopes requested: " + multiplier);
			}
		}
		return new IsotopeSpecification(chemEl, isotope, multiplier, locants);
	}
}
