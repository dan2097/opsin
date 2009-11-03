package uk.ac.cam.ch.wwmm.opsin;

/**
 * Takes a name:
 * strips leading/trailing white space
 * rejects a few special cases
 * Identifies polymer names
 * @author dl387
 *
 */
class PreProcessor {
	//private Pattern semiColon;
	
	PreProcessor() {
		//semiColon =Pattern.compile(";");
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
		
		//chemicalName=semiColon.matcher(chemicalName).replaceAll(" ");

		return chemicalName;
	}
}
