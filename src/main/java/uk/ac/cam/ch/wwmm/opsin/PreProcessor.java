package uk.ac.cam.ch.wwmm.opsin;

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
	/**
	 * Wrapper class for returning multiple objects
	 */
	final class PreProcessorResults {
		final String chemicalName;
		final OpsinMode mode;

		public PreProcessorResults(String chemicalName, OpsinMode mode) {
			this.chemicalName = chemicalName;
			this.mode = mode;
		}

		public String getChemicalName() {
			return chemicalName;
		}

		public OpsinMode getMode() {
			return mode;
		}
	}
	
	/**
	 * OPSIN modes. Chosen by the form of the name as analysed by the Preprocessor
	 */
	enum OpsinMode {
	     normal,
	     poly
	}
	
	//private Pattern semiColon;
	private Pattern matchPoly;
	PreProcessor() {
		//semiColon =Pattern.compile(";");
		matchPoly =Pattern.compile("(?:poly|oligo)[\\[\\(\\{](.+)[\\]\\)\\}]", Pattern.CASE_INSENSITIVE );//poly or oligo followed by a bracket name
	}
	/**
	 * Master method for PreProcessing
	 * @param chemicalName
	 * @return
	 */
	PreProcessorResults preProcess(String chemicalName) {
		OpsinMode mode = OpsinMode.normal;
		chemicalName=chemicalName.trim();//remove leading and trailing whitespace
		if (chemicalName.equals("")){return null;}
		if("amine".equalsIgnoreCase(chemicalName)) return null;//trigenericammonia
		if("thiol".equalsIgnoreCase(chemicalName)) return null;//genericsulfane
		if("carboxylic acid".equalsIgnoreCase(chemicalName)) return null;//genericmethanoic acid
		//Alcohol Aldehyde Alkane Alkene Alkyne Amide Amine Azo compound Benzene derivative Carboxylic acid Cyanate Disulfide Ester Ether Haloalkane Hydrazone Imine Isocyanide Isocyanate Ketone Oxime Nitrile Nitro compound Nitroso compound Peroxide Phosphoric acid Pyridine derivative Sulfone Sulfonic acid Sulfoxide Thioester Thioether Thiol
		
		//chemicalName=semiColon.matcher(chemicalName).replaceAll(" ");

		Matcher polyMatcher = matchPoly.matcher(chemicalName);
		if (polyMatcher.matches()){//Name is a polymer/oligomer name
			chemicalName = polyMatcher.group(1);//strip off starting poly( and trailing closing bracket
			mode = OpsinMode.poly;
		}
		return new PreProcessorResults(chemicalName, mode);
	}
}
