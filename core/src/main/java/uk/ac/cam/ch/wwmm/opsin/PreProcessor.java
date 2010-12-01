package uk.ac.cam.ch.wwmm.opsin;

import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Takes a name:
 * strips leading/trailing white space
 * rejects a few special cases
 * Identifies polymer names
 * @author dl387
 *
 */
class PreProcessor {
	private static final Pattern MATCH_DOLLAR = Pattern.compile("\\$");
	private static final Pattern MATCH_SULPH = Pattern.compile("sulph");
	private static final HashMap<String, String> GREEK_MAP = new HashMap<String, String>();

	private static final String AMIDE = "amide";
	private static final String AMINE = "amine";
	private static final String THIOL = "thiol";

	static {
		GREEK_MAP.put("a", "alpha");
		GREEK_MAP.put("b", "beta");
		GREEK_MAP.put("g", "gamma");
		GREEK_MAP.put("d", "delta");
		GREEK_MAP.put("e", "epsilon");
//		letterGreekMap.put("z", "zeta");
//		letterGreekMap.put("i", "iota");
//		letterGreekMap.put("k", "kappa");
		GREEK_MAP.put("l", "lambda");
//		letterGreekMap.put("m", "mu");
//		letterGreekMap.put("n", "nu");
//		letterGreekMap.put("x", "xi");
//		letterGreekMap.put("p", "pi");
//		letterGreekMap.put("r", "rho");
//		letterGreekMap.put("s", "sigma");
//		letterGreekMap.put("t", "tau");
//		letterGreekMap.put("u", "upsilon");
//		letterGreekMap.put("f", "phi");
//		letterGreekMap.put("o", "omega");
	}

	/**
	 * Master method for PreProcessing
	 * @param chemicalName
	 * @return
	 * @throws PreProcessingException 
	 */
	static String preProcess(String chemicalName) throws PreProcessingException {
		chemicalName=chemicalName.trim();//remove leading and trailing whitespace
		if ("".equals(chemicalName)){
			throw new PreProcessingException("Input chemical name was blank!");
		}

		if(AMIDE.equalsIgnoreCase(chemicalName)) {
			throw new PreProcessingException("Amide is a generic term rather than a specific chemical");//amide
		}
		if(AMINE.equalsIgnoreCase(chemicalName)) {
			throw new PreProcessingException("Amine is a generic term rather than a specific chemical");//trigenericammonia
		}
		if(THIOL.equalsIgnoreCase(chemicalName)) {
			throw new PreProcessingException("Thiol is a generic term rather than a specific chemical");//genericsulfane
		}
		chemicalName = processDollarPrefixedGreeks(chemicalName);
		chemicalName = StringTools.convertNonAsciiAndNormaliseRepresentation(chemicalName);
		chemicalName = MATCH_SULPH.matcher(chemicalName).replaceAll("sulf");//correct British spelling to the IUPAC spelling
		return chemicalName;
	}

	private static String processDollarPrefixedGreeks(String chemicalName) {
		//StringTools.unicodeToLatin(chemicalName);
		Matcher m = MATCH_DOLLAR.matcher(chemicalName);
		while (m.find()){
			if (chemicalName.length()>m.end()){
				String letter = chemicalName.substring(m.end(), m.end()+1).toLowerCase();
				if (GREEK_MAP.containsKey(letter)){
					chemicalName = chemicalName.substring(0, m.end()-1) +GREEK_MAP.get(letter) +  chemicalName.substring(m.end()+1);
					m = MATCH_DOLLAR.matcher(chemicalName);
				}
			}
		}
		return chemicalName;
	}
}
