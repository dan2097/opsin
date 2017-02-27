package uk.ac.cam.ch.wwmm.opsin;

import static uk.ac.cam.ch.wwmm.opsin.OpsinTools.MATCH_COMMA;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

class IsotopeSpecificationParser {
	
	private static final Pattern matchBoughtonIsotope =Pattern.compile("(?:-([^,]+(?:,[^,]+)*))?-d(\\d+)?");
	private static final Pattern matchIupacIsotope =Pattern.compile("(?:([^,]+(?:,[^,]+)*)-)?(\\d+)([A-Z][a-z]?)(\\d+)?");

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
		String type = isotopeSpecification.getAttributeValue(XmlDeclarations.TYPE_ATR);
		if (XmlDeclarations.BOUGHTONSYSTEM_TYPE_VAL.equals(type)) {
			return processBoughtonIsotope(isotopeSpecification);
		}
		else if (XmlDeclarations.IUPACSYSTEM_TYPE_VAL.equals(type)) {
			return processIupacIsotope(isotopeSpecification);
		}
		else {
			throw new RuntimeException("Unsupported isotope specification syntax");
		}
	}

	private static IsotopeSpecification processBoughtonIsotope(Element isotopeSpecification) throws StructureBuildingException {
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
	
	private static IsotopeSpecification processIupacIsotope(Element isotopeSpecification) throws StructureBuildingException {
		String val = isotopeSpecification.getValue();
		Matcher m = matchIupacIsotope.matcher(val);
		if (!m.matches()) {
			throw new RuntimeException("Malformed isotope specification: " + val);
		}
		
		int isotope = Integer.parseInt(m.group(2));
		ChemEl chemEl = ChemEl.valueOf(m.group(3));

		int multiplier = 1;
		String multiplierStr = m.group(4);
		if (multiplierStr != null) {
			multiplier = Integer.parseInt(multiplierStr);
		}
		
		String locantsStr = m.group(1);
		String[] locants = null;
		if(locantsStr != null) {
			locants = MATCH_COMMA.split(locantsStr);
			if (multiplierStr == null) {
				multiplier = locants.length;
			}
			else if (locants.length != multiplier) {
				throw new StructureBuildingException("Mismatch between number of locants: " + locants.length + " and number of " + chemEl.toString() +" isotopes requested: " + multiplier);
			}
		}
		return new IsotopeSpecification(chemEl, isotope, multiplier, locants);
	}
}
