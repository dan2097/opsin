package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.HashSet;

import uk.ac.cam.ch.wwmm.opsin.StringTools;

public class InchiPruner {
	/**
	 * removes all bar main layer(connection table and hydrogen) and charge layers. 
	 * @param inchi
	 * @return
	 */
	public static String mainAndChargeLayers(String inchi){
		String[] inchiArray = inchi.split("/");
		ArrayList<String> inchiArrayList =new ArrayList<String>();
		if (inchiArray[0].endsWith("S")){//remove optional S
			inchiArray[0]=inchiArray[0].substring(0, inchiArray[0].length() -1);
		}
		if (inchiArray.length >0) {inchiArrayList.add(inchiArray[0]);}
		if (inchiArray.length >1) {inchiArrayList.add(inchiArray[1]);}
		HashSet<Character> startingChars = new HashSet<Character>();
		for (int i = 2; i < inchiArray.length; i++) {
			Character c = inchiArray[i].charAt(0);
			if (startingChars.contains(c)){
				throw new RuntimeException(inchi + " is invalid or there is a bug in InChI cleaner. Multiple layers starting with " + c + " encountered!");
			}
			if (c=='r' || c=='f' || c=='i'){
				break;
			}
			if (c=='c' || c=='h' || c=='q' || c=='p'){
				inchiArrayList.add(inchiArray[i]);
			}
			startingChars.add(c);
		}
		inchi =StringTools.stringListToString(inchiArrayList, "/");
		
		return inchi;
	}
	
	/**
	 * removes all bar main layer(connection table and hydrogen), charge layers and stereochemistry. 
	 * @param inchi
	 * @return
	 */
	public static String mainChargeAndStereochemistryLayers(String inchi){
		String[] inchiArray = inchi.split("/");
		ArrayList<String> inchiArrayList =new ArrayList<String>();
		if (inchiArray[0].endsWith("S")){//remove optional S
			inchiArray[0]=inchiArray[0].substring(0, inchiArray[0].length() -1);
		}
		if (inchiArray.length >0) {inchiArrayList.add(inchiArray[0]);}
		if (inchiArray.length >1) {inchiArrayList.add(inchiArray[1]);}
		HashSet<Character> startingChars = new HashSet<Character>();
		for (int i = 2; i < inchiArray.length; i++) {
			Character c = inchiArray[i].charAt(0);
			if (startingChars.contains(c)){
				throw new RuntimeException(inchi + " is invalid or there is a bug in InChI cleaner. Multiple layers starting with " + c + " encountered!");
			}
			if (c=='r' || c=='f' || c=='i'){
				break;
			}
			if (c=='c' || c=='h' || c=='q' || c=='p' || c=='m' || c=='t' || c=='s'){
				inchiArrayList.add(inchiArray[i]);
			}
			startingChars.add(c);
		}
		inchi =StringTools.stringListToString(inchiArrayList, "/");
		
		return inchi;
	}
}
