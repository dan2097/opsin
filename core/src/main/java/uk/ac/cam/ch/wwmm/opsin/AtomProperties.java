package uk.ac.cam.ch.wwmm.opsin;

import java.util.HashMap;
import java.util.Map;

/**
 * Holds useful atomic properties
 * @author dl387
 *
 */
class AtomProperties {
	/**
	 * Useful to give an indication of whether a bond is like to be ionic (diff >1.8), polar or covalent (diff < 1.2)
	 */
	static final Map<String, Double> elementToPaulingElectronegativity = new HashMap<String, Double>();
	
	static final Map<String, Integer> elementToAtomicNumber = new HashMap<String, Integer>();
	
	/**
	 * Maps element symbol to the priority of that atom in Hantzch-Widman system. A higher value indicates a higher priority.
	 */
	static final Map<String, Integer> elementToHwPriority = new HashMap<String, Integer>();
	
	static{
		elementToPaulingElectronegativity.put("H", 2.20);
		elementToPaulingElectronegativity.put("Li", 0.98);
		elementToPaulingElectronegativity.put("Be", 1.57);
		elementToPaulingElectronegativity.put("B", 2.04);
		elementToPaulingElectronegativity.put("C", 2.55);
		elementToPaulingElectronegativity.put("N", 3.04);
		elementToPaulingElectronegativity.put("O", 3.44);
		elementToPaulingElectronegativity.put("F", 3.98);
		elementToPaulingElectronegativity.put("Na", 0.93);
		elementToPaulingElectronegativity.put("Mg", 1.31);
		elementToPaulingElectronegativity.put("Al", 1.61);
		elementToPaulingElectronegativity.put("Si", 1.90);
		elementToPaulingElectronegativity.put("P", 2.19);
		elementToPaulingElectronegativity.put("S", 2.58);
		elementToPaulingElectronegativity.put("Cl", 3.16);
		elementToPaulingElectronegativity.put("K", 0.82);
		elementToPaulingElectronegativity.put("Ca", 1.00);
		elementToPaulingElectronegativity.put("Sc", 1.36);
		elementToPaulingElectronegativity.put("Ti", 1.54);
		elementToPaulingElectronegativity.put("V", 1.63);
		elementToPaulingElectronegativity.put("Cr", 1.66);
		elementToPaulingElectronegativity.put("Mn", 1.55);
		elementToPaulingElectronegativity.put("Fe", 1.83);
		elementToPaulingElectronegativity.put("Co", 1.88);
		elementToPaulingElectronegativity.put("Ni", 1.91);
		elementToPaulingElectronegativity.put("Cu", 1.90);
		elementToPaulingElectronegativity.put("Zn", 1.65);
		elementToPaulingElectronegativity.put("Ga", 1.81);
		elementToPaulingElectronegativity.put("Ge", 2.01);
		elementToPaulingElectronegativity.put("As", 2.18);
		elementToPaulingElectronegativity.put("Se", 2.55);
		elementToPaulingElectronegativity.put("Br", 2.96);
		elementToPaulingElectronegativity.put("Kr", 3.00);
		elementToPaulingElectronegativity.put("Rb", 0.82);
		elementToPaulingElectronegativity.put("Sr", 0.95);
		elementToPaulingElectronegativity.put("Y", 1.22);
		elementToPaulingElectronegativity.put("Zr", 1.33);
		elementToPaulingElectronegativity.put("Nb", 1.6);
		elementToPaulingElectronegativity.put("Mo", 2.16);
		elementToPaulingElectronegativity.put("Tc", 1.9);
		elementToPaulingElectronegativity.put("Ru", 2.2);
		elementToPaulingElectronegativity.put("Rh", 2.28);
		elementToPaulingElectronegativity.put("Pd", 2.20);
		elementToPaulingElectronegativity.put("Ag", 1.93);
		elementToPaulingElectronegativity.put("Cd", 1.69);
		elementToPaulingElectronegativity.put("In", 1.78);
		elementToPaulingElectronegativity.put("Sn", 1.96);
		elementToPaulingElectronegativity.put("Sb", 2.05);
		elementToPaulingElectronegativity.put("Te", 2.1);
		elementToPaulingElectronegativity.put("I", 2.66);
		elementToPaulingElectronegativity.put("Xe", 2.60);
		elementToPaulingElectronegativity.put("Cs", 0.79);
		elementToPaulingElectronegativity.put("Ba", 0.89);
		elementToPaulingElectronegativity.put("La", 1.1);
		elementToPaulingElectronegativity.put("Ce", 1.12);
		elementToPaulingElectronegativity.put("Pr", 1.13);
		elementToPaulingElectronegativity.put("Nd", 1.14);
		elementToPaulingElectronegativity.put("Pm", 1.13);
		elementToPaulingElectronegativity.put("Sm", 1.17);
		elementToPaulingElectronegativity.put("Eu", 1.2);
		elementToPaulingElectronegativity.put("Gd", 1.2);
		elementToPaulingElectronegativity.put("Tb", 1.1);
		elementToPaulingElectronegativity.put("Dy", 1.22);
		elementToPaulingElectronegativity.put("Ho", 1.23);
		elementToPaulingElectronegativity.put("Er", 1.24);
		elementToPaulingElectronegativity.put("Tm", 1.25);
		elementToPaulingElectronegativity.put("Yb", 1.1);
		elementToPaulingElectronegativity.put("Lu", 1.27);
		elementToPaulingElectronegativity.put("Hf", 1.3);
		elementToPaulingElectronegativity.put("Ta", 1.5);
		elementToPaulingElectronegativity.put("W", 2.36);
		elementToPaulingElectronegativity.put("Re", 1.9);
		elementToPaulingElectronegativity.put("Os", 2.2);
		elementToPaulingElectronegativity.put("Ir", 2.20);
		elementToPaulingElectronegativity.put("Pt", 2.28);
		elementToPaulingElectronegativity.put("Au", 2.54);
		elementToPaulingElectronegativity.put("Hg", 2.00);
		elementToPaulingElectronegativity.put("Tl", 1.62);
		elementToPaulingElectronegativity.put("Pb", 2.33);
		elementToPaulingElectronegativity.put("Bi", 2.02);
		elementToPaulingElectronegativity.put("Po", 2.0);
		elementToPaulingElectronegativity.put("At", 2.2);
		elementToPaulingElectronegativity.put("Rn", 2.2);
		elementToPaulingElectronegativity.put("Fr", 0.7);
		elementToPaulingElectronegativity.put("Ra", 0.9);
		elementToPaulingElectronegativity.put("Ac", 1.1);
		elementToPaulingElectronegativity.put("Th", 1.3);
		elementToPaulingElectronegativity.put("Pa", 1.5);
		elementToPaulingElectronegativity.put("U", 1.38);
		elementToPaulingElectronegativity.put("Np", 1.36);
		elementToPaulingElectronegativity.put("Pu", 1.28);
		elementToPaulingElectronegativity.put("Am", 1.13);
		elementToPaulingElectronegativity.put("Cm", 1.28);
		elementToPaulingElectronegativity.put("Bk", 1.3);
		elementToPaulingElectronegativity.put("Cf", 1.3);
		elementToPaulingElectronegativity.put("Es", 1.3);
		elementToPaulingElectronegativity.put("Fm", 1.3);
		elementToPaulingElectronegativity.put("Md", 1.3);
		elementToPaulingElectronegativity.put("No", 1.3);
		elementToPaulingElectronegativity.put("Lr", 1.3);
		
		elementToAtomicNumber.put("H", 1);
		elementToAtomicNumber.put("He", 2);
		elementToAtomicNumber.put("Li", 3);
		elementToAtomicNumber.put("Be", 4);
		elementToAtomicNumber.put("B", 5);
		elementToAtomicNumber.put("C", 6);
		elementToAtomicNumber.put("N", 7);
		elementToAtomicNumber.put("O", 8);
		elementToAtomicNumber.put("F", 9);
		elementToAtomicNumber.put("Ne", 10);
		elementToAtomicNumber.put("Na", 11);
		elementToAtomicNumber.put("Mg", 12);
		elementToAtomicNumber.put("Al", 13);
		elementToAtomicNumber.put("Si", 14);
		elementToAtomicNumber.put("P", 15);
		elementToAtomicNumber.put("S", 16);
		elementToAtomicNumber.put("Cl", 17);
		elementToAtomicNumber.put("Ar", 18);
		elementToAtomicNumber.put("K", 19);
		elementToAtomicNumber.put("Ca", 20);
		elementToAtomicNumber.put("Sc", 21);
		elementToAtomicNumber.put("Ti", 22);
		elementToAtomicNumber.put("V", 23);
		elementToAtomicNumber.put("Cr", 24);
		elementToAtomicNumber.put("Mn", 25);
		elementToAtomicNumber.put("Fe", 26);
		elementToAtomicNumber.put("Co", 27);
		elementToAtomicNumber.put("Ni", 28);
		elementToAtomicNumber.put("Cu", 29);
		elementToAtomicNumber.put("Zn", 30);
		elementToAtomicNumber.put("Ga", 31);
		elementToAtomicNumber.put("Ge", 32);
		elementToAtomicNumber.put("As", 33);
		elementToAtomicNumber.put("Se", 34);
		elementToAtomicNumber.put("Br", 35);
		elementToAtomicNumber.put("Kr", 36);
		elementToAtomicNumber.put("Rb", 37);
		elementToAtomicNumber.put("Sr", 38);
		elementToAtomicNumber.put("Y", 39);
		elementToAtomicNumber.put("Zr", 40);
		elementToAtomicNumber.put("Nb", 41);
		elementToAtomicNumber.put("Mo", 42);
		elementToAtomicNumber.put("Tc", 43);
		elementToAtomicNumber.put("Ru", 44);
		elementToAtomicNumber.put("Rh", 45);
		elementToAtomicNumber.put("Pd", 46);
		elementToAtomicNumber.put("Ag", 47);
		elementToAtomicNumber.put("Cd", 48);
		elementToAtomicNumber.put("In", 49);
		elementToAtomicNumber.put("Sn", 50);
		elementToAtomicNumber.put("Sb", 51);
		elementToAtomicNumber.put("Te", 52);
		elementToAtomicNumber.put("I", 53);
		elementToAtomicNumber.put("Xe", 54);
		elementToAtomicNumber.put("Cs", 55);
		elementToAtomicNumber.put("Ba", 56);
		elementToAtomicNumber.put("La", 57);
		elementToAtomicNumber.put("Ce", 58);
		elementToAtomicNumber.put("Pr", 59);
		elementToAtomicNumber.put("Nd", 60);
		elementToAtomicNumber.put("Pm", 61);
		elementToAtomicNumber.put("Sm", 62);
		elementToAtomicNumber.put("Eu", 63);
		elementToAtomicNumber.put("Gd", 64);
		elementToAtomicNumber.put("Tb", 65);
		elementToAtomicNumber.put("Dy", 66);
		elementToAtomicNumber.put("Ho", 67);
		elementToAtomicNumber.put("Er", 68);
		elementToAtomicNumber.put("Tm", 69);
		elementToAtomicNumber.put("Yb", 70);
		elementToAtomicNumber.put("Lu", 71);
		elementToAtomicNumber.put("Hf", 72);
		elementToAtomicNumber.put("Ta", 73);
		elementToAtomicNumber.put("W", 74);
		elementToAtomicNumber.put("Re", 75);
		elementToAtomicNumber.put("Os", 76);
		elementToAtomicNumber.put("Ir", 77);
		elementToAtomicNumber.put("Pt", 78);
		elementToAtomicNumber.put("Au", 79);
		elementToAtomicNumber.put("Hg", 80);
		elementToAtomicNumber.put("Tl", 81);
		elementToAtomicNumber.put("Pb", 82);
		elementToAtomicNumber.put("Bi", 83);
		elementToAtomicNumber.put("Po", 84);
		elementToAtomicNumber.put("At", 85);
		elementToAtomicNumber.put("Rn", 86);
		elementToAtomicNumber.put("Fr", 87);
		elementToAtomicNumber.put("Ra", 88);
		elementToAtomicNumber.put("Ac", 89);
		elementToAtomicNumber.put("Th", 90);
		elementToAtomicNumber.put("Pa", 91);
		elementToAtomicNumber.put("U", 92);
		elementToAtomicNumber.put("Np", 93);
		elementToAtomicNumber.put("Pu", 94);
		elementToAtomicNumber.put("Am", 95);
		elementToAtomicNumber.put("Cm", 96);
		elementToAtomicNumber.put("Bk", 97);
		elementToAtomicNumber.put("Cf", 98);
		elementToAtomicNumber.put("Es", 99);
		elementToAtomicNumber.put("Fm", 100);
		elementToAtomicNumber.put("Md", 101);
		elementToAtomicNumber.put("No", 102);
		elementToAtomicNumber.put("Lr", 103);
		elementToAtomicNumber.put("Rf", 104);
		elementToAtomicNumber.put("Db", 105);
		elementToAtomicNumber.put("Sg", 106);
		elementToAtomicNumber.put("Bh", 107);
		elementToAtomicNumber.put("Hs", 108);
		elementToAtomicNumber.put("Mt", 109);
		elementToAtomicNumber.put("Ds", 110);

		elementToHwPriority.put("F", 23);
		elementToHwPriority.put("Cl", 22);
		elementToHwPriority.put("Br", 21);
		elementToHwPriority.put("I", 20);
		elementToHwPriority.put("O", 19);
		elementToHwPriority.put("S", 18);
		elementToHwPriority.put("Se", 17);
		elementToHwPriority.put("Te", 16);
		elementToHwPriority.put("N", 15);
		elementToHwPriority.put("P", 14);
		elementToHwPriority.put("As", 13);
		elementToHwPriority.put("Sb", 12);
		elementToHwPriority.put("Bi", 11);
		elementToHwPriority.put("Si", 10);
		elementToHwPriority.put("Ge", 9);
		elementToHwPriority.put("Sn", 8);
		elementToHwPriority.put("Pb", 7);
		elementToHwPriority.put("B", 6);
		elementToHwPriority.put("Al", 5);
		elementToHwPriority.put("Ga", 4);
		elementToHwPriority.put("In", 3);
		elementToHwPriority.put("Tl", 2);
		elementToHwPriority.put("Hg", 1);
	}
}
