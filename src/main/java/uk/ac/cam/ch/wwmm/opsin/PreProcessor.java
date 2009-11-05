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
	//private Pattern semiColon =Pattern.compile(";");
	private Pattern matchDollar = Pattern.compile("\\$");

	private HashMap<String, String> letterGreekMap = new HashMap<String, String>();
	
	public PreProcessor() {
		letterGreekMap.put("a", "alpha");
		letterGreekMap.put("b", "beta");
		letterGreekMap.put("g", "gamma");
		letterGreekMap.put("d", "delta");
		letterGreekMap.put("e", "epsilon");
//		letterGreekMap.put("z", "zeta");
//		letterGreekMap.put("i", "iota");
//		letterGreekMap.put("k", "kappa");
		letterGreekMap.put("l", "lambda");
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
	 */
	String preProcess(String chemicalName) {
		chemicalName=chemicalName.trim();//remove leading and trailing whitespace
		if (chemicalName.equals("")){return null;}
		if("amine".equalsIgnoreCase(chemicalName)) return null;//trigenericammonia
		if("thiol".equalsIgnoreCase(chemicalName)) return null;//genericsulfane
		if("carboxylic acid".equalsIgnoreCase(chemicalName)) return null;//genericmethanoic acid
		//Alcohol Aldehyde Alkane Alkene Alkyne Amide Amine Azo compound Benzene derivative Carboxylic acid Cyanate Disulfide Ester Ether Haloalkane Hydrazone Imine Isocyanide Isocyanate Ketone Oxime Nitrile Nitro compound Nitroso compound Peroxide Phosphoric acid Pyridine derivative Sulfone Sulfonic acid Sulfoxide Thioester Thioether Thiol
		
		Matcher m = matchDollar.matcher(chemicalName);
		while (m.find()){
			if (chemicalName.length()>m.end()){
				String letter = chemicalName.substring(m.end(), m.end()+1).toLowerCase();
				if (letterGreekMap.containsKey(letter)){
					chemicalName = chemicalName.substring(0, m.end()-1) +letterGreekMap.get(letter) +  chemicalName.substring(m.end()+1);
					m = matchDollar.matcher(chemicalName);
				}
			}
		}
		
		//chemicalName=semiColon.matcher(chemicalName).replaceAll(" ");

		return chemicalName;
	}
}
