package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.List;

import uk.ac.cam.ch.wwmm.opsin.StringTools;

public class InchiPruner {
	/**
	 * Return a modified version of the given InChI where the:
	 * stereochemistry, fixed hydrogen and reconnected layers have been removed
	 * The S indicating standard InChI is also removed
	 * @param inchi
	 * @return
	 */
	public static String mainAndChargeLayers(String inchi){
		String[] inchiLayers = inchi.split("/");
		if (inchiLayers.length < 2){
			return inchi;
		}
		List<String> retainedLayers = new ArrayList<String>();
		if (Character.isLetter(inchiLayers[0].charAt(inchiLayers[0].length() -1))){//remove the S indicating this to be a standard InChI
			inchiLayers[0]=inchiLayers[0].substring(0, inchiLayers[0].length() -1);
		}
		retainedLayers.add(inchiLayers[0]);//version identifier
		retainedLayers.add(inchiLayers[1]);//molecular formula

		for (int i = 2; i < inchiLayers.length; i++) {
			Character c = inchiLayers[i].charAt(0);
			if (c=='c' || c=='h' || c=='q' || c=='p' || c=='i'){
				retainedLayers.add(inchiLayers[i]);
			}
			else if (c!='b' && c!='t' && c!='m' && c!='s'){//ignore stereochemistry but continue as there may be an isotopic layer
				break;
			}
		}
		return StringTools.stringListToString(retainedLayers, "/");
	}
	
	/**
	 * Return a modified version of the given InChI where the:
	 * fixed hydrogen and reconnected layers have been removed
	 * The S indicating standard InChI is also removed
	 * @param inchi
	 * @return
	 */
	public static String mainChargeAndStereochemistryLayers(String inchi){
		String[] inchiLayers = inchi.split("/");
		if (inchiLayers.length < 2){
			return inchi;
		}
		List<String> retainedLayers = new ArrayList<String>();
		if (Character.isLetter(inchiLayers[0].charAt(inchiLayers[0].length() -1))){//remove the S indicating this to be a standard InChI
			inchiLayers[0]=inchiLayers[0].substring(0, inchiLayers[0].length() -1);
		}
		retainedLayers.add(inchiLayers[0]);//version identifier
		retainedLayers.add(inchiLayers[1]);//molecular formula

		for (int i = 2; i < inchiLayers.length; i++) {
			Character c = inchiLayers[i].charAt(0);
			if (c=='c' || c=='h' || c=='q' || c=='p' || c=='b' || c=='t' || c=='m' || c=='s' || c=='i'){
				retainedLayers.add(inchiLayers[i]);
			}
			else{
				break;
			}
		}
		return StringTools.stringListToString(retainedLayers, "/");
	}
}
