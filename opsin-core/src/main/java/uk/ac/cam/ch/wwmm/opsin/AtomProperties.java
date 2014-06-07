package uk.ac.cam.ch.wwmm.opsin;

import java.util.EnumMap;
import java.util.Map;

/**
 * Holds useful atomic properties
 * @author dl387
 *
 */
class AtomProperties {

	private static final Map<ChemEl, Double> elementToPaulingElectronegativity = new EnumMap<ChemEl, Double>(ChemEl.class);
	private static final Map<ChemEl, Integer> elementToHwPriority = new EnumMap<ChemEl, Integer>(ChemEl.class);
	
	static{
		elementToPaulingElectronegativity.put(ChemEl.H, 2.20);
		elementToPaulingElectronegativity.put(ChemEl.Li, 0.98);
		elementToPaulingElectronegativity.put(ChemEl.Be, 1.57);
		elementToPaulingElectronegativity.put(ChemEl.B, 2.04);
		elementToPaulingElectronegativity.put(ChemEl.C, 2.55);
		elementToPaulingElectronegativity.put(ChemEl.N, 3.04);
		elementToPaulingElectronegativity.put(ChemEl.O, 3.44);
		elementToPaulingElectronegativity.put(ChemEl.F, 3.98);
		elementToPaulingElectronegativity.put(ChemEl.Na, 0.93);
		elementToPaulingElectronegativity.put(ChemEl.Mg, 1.31);
		elementToPaulingElectronegativity.put(ChemEl.Al, 1.61);
		elementToPaulingElectronegativity.put(ChemEl.Si, 1.90);
		elementToPaulingElectronegativity.put(ChemEl.P, 2.19);
		elementToPaulingElectronegativity.put(ChemEl.S, 2.58);
		elementToPaulingElectronegativity.put(ChemEl.Cl, 3.16);
		elementToPaulingElectronegativity.put(ChemEl.K, 0.82);
		elementToPaulingElectronegativity.put(ChemEl.Ca, 1.00);
		elementToPaulingElectronegativity.put(ChemEl.Sc, 1.36);
		elementToPaulingElectronegativity.put(ChemEl.Ti, 1.54);
		elementToPaulingElectronegativity.put(ChemEl.V, 1.63);
		elementToPaulingElectronegativity.put(ChemEl.Cr, 1.66);
		elementToPaulingElectronegativity.put(ChemEl.Mn, 1.55);
		elementToPaulingElectronegativity.put(ChemEl.Fe, 1.83);
		elementToPaulingElectronegativity.put(ChemEl.Co, 1.88);
		elementToPaulingElectronegativity.put(ChemEl.Ni, 1.91);
		elementToPaulingElectronegativity.put(ChemEl.Cu, 1.90);
		elementToPaulingElectronegativity.put(ChemEl.Zn, 1.65);
		elementToPaulingElectronegativity.put(ChemEl.Ga, 1.81);
		elementToPaulingElectronegativity.put(ChemEl.Ge, 2.01);
		elementToPaulingElectronegativity.put(ChemEl.As, 2.18);
		elementToPaulingElectronegativity.put(ChemEl.Se, 2.55);
		elementToPaulingElectronegativity.put(ChemEl.Br, 2.96);
		elementToPaulingElectronegativity.put(ChemEl.Kr, 3.00);
		elementToPaulingElectronegativity.put(ChemEl.Rb, 0.82);
		elementToPaulingElectronegativity.put(ChemEl.Sr, 0.95);
		elementToPaulingElectronegativity.put(ChemEl.Y, 1.22);
		elementToPaulingElectronegativity.put(ChemEl.Zr, 1.33);
		elementToPaulingElectronegativity.put(ChemEl.Nb, 1.6);
		elementToPaulingElectronegativity.put(ChemEl.Mo, 2.16);
		elementToPaulingElectronegativity.put(ChemEl.Tc, 1.9);
		elementToPaulingElectronegativity.put(ChemEl.Ru, 2.2);
		elementToPaulingElectronegativity.put(ChemEl.Rh, 2.28);
		elementToPaulingElectronegativity.put(ChemEl.Pd, 2.20);
		elementToPaulingElectronegativity.put(ChemEl.Ag, 1.93);
		elementToPaulingElectronegativity.put(ChemEl.Cd, 1.69);
		elementToPaulingElectronegativity.put(ChemEl.In, 1.78);
		elementToPaulingElectronegativity.put(ChemEl.Sn, 1.96);
		elementToPaulingElectronegativity.put(ChemEl.Sb, 2.05);
		elementToPaulingElectronegativity.put(ChemEl.Te, 2.1);
		elementToPaulingElectronegativity.put(ChemEl.I, 2.66);
		elementToPaulingElectronegativity.put(ChemEl.Xe, 2.60);
		elementToPaulingElectronegativity.put(ChemEl.Cs, 0.79);
		elementToPaulingElectronegativity.put(ChemEl.Ba, 0.89);
		elementToPaulingElectronegativity.put(ChemEl.La, 1.1);
		elementToPaulingElectronegativity.put(ChemEl.Ce, 1.12);
		elementToPaulingElectronegativity.put(ChemEl.Pr, 1.13);
		elementToPaulingElectronegativity.put(ChemEl.Nd, 1.14);
		elementToPaulingElectronegativity.put(ChemEl.Pm, 1.13);
		elementToPaulingElectronegativity.put(ChemEl.Sm, 1.17);
		elementToPaulingElectronegativity.put(ChemEl.Eu, 1.2);
		elementToPaulingElectronegativity.put(ChemEl.Gd, 1.2);
		elementToPaulingElectronegativity.put(ChemEl.Tb, 1.1);
		elementToPaulingElectronegativity.put(ChemEl.Dy, 1.22);
		elementToPaulingElectronegativity.put(ChemEl.Ho, 1.23);
		elementToPaulingElectronegativity.put(ChemEl.Er, 1.24);
		elementToPaulingElectronegativity.put(ChemEl.Tm, 1.25);
		elementToPaulingElectronegativity.put(ChemEl.Yb, 1.1);
		elementToPaulingElectronegativity.put(ChemEl.Lu, 1.27);
		elementToPaulingElectronegativity.put(ChemEl.Hf, 1.3);
		elementToPaulingElectronegativity.put(ChemEl.Ta, 1.5);
		elementToPaulingElectronegativity.put(ChemEl.W, 2.36);
		elementToPaulingElectronegativity.put(ChemEl.Re, 1.9);
		elementToPaulingElectronegativity.put(ChemEl.Os, 2.2);
		elementToPaulingElectronegativity.put(ChemEl.Ir, 2.20);
		elementToPaulingElectronegativity.put(ChemEl.Pt, 2.28);
		elementToPaulingElectronegativity.put(ChemEl.Au, 2.54);
		elementToPaulingElectronegativity.put(ChemEl.Hg, 2.00);
		elementToPaulingElectronegativity.put(ChemEl.Tl, 1.62);
		elementToPaulingElectronegativity.put(ChemEl.Pb, 2.33);
		elementToPaulingElectronegativity.put(ChemEl.Bi, 2.02);
		elementToPaulingElectronegativity.put(ChemEl.Po, 2.0);
		elementToPaulingElectronegativity.put(ChemEl.At, 2.2);
		elementToPaulingElectronegativity.put(ChemEl.Rn, 2.2);
		elementToPaulingElectronegativity.put(ChemEl.Fr, 0.7);
		elementToPaulingElectronegativity.put(ChemEl.Ra, 0.9);
		elementToPaulingElectronegativity.put(ChemEl.Ac, 1.1);
		elementToPaulingElectronegativity.put(ChemEl.Th, 1.3);
		elementToPaulingElectronegativity.put(ChemEl.Pa, 1.5);
		elementToPaulingElectronegativity.put(ChemEl.U, 1.38);
		elementToPaulingElectronegativity.put(ChemEl.Np, 1.36);
		elementToPaulingElectronegativity.put(ChemEl.Pu, 1.28);
		elementToPaulingElectronegativity.put(ChemEl.Am, 1.13);
		elementToPaulingElectronegativity.put(ChemEl.Cm, 1.28);
		elementToPaulingElectronegativity.put(ChemEl.Bk, 1.3);
		elementToPaulingElectronegativity.put(ChemEl.Cf, 1.3);
		elementToPaulingElectronegativity.put(ChemEl.Es, 1.3);
		elementToPaulingElectronegativity.put(ChemEl.Fm, 1.3);
		elementToPaulingElectronegativity.put(ChemEl.Md, 1.3);
		elementToPaulingElectronegativity.put(ChemEl.No, 1.3);
		elementToPaulingElectronegativity.put(ChemEl.Lr, 1.3);

		elementToHwPriority.put(ChemEl.F, 23);
		elementToHwPriority.put(ChemEl.Cl, 22);
		elementToHwPriority.put(ChemEl.Br, 21);
		elementToHwPriority.put(ChemEl.I, 20);
		elementToHwPriority.put(ChemEl.O, 19);
		elementToHwPriority.put(ChemEl.S, 18);
		elementToHwPriority.put(ChemEl.Se, 17);
		elementToHwPriority.put(ChemEl.Te, 16);
		elementToHwPriority.put(ChemEl.N, 15);
		elementToHwPriority.put(ChemEl.P, 14);
		elementToHwPriority.put(ChemEl.As, 13);
		elementToHwPriority.put(ChemEl.Sb, 12);
		elementToHwPriority.put(ChemEl.Bi, 11);
		elementToHwPriority.put(ChemEl.Si, 10);
		elementToHwPriority.put(ChemEl.Ge, 9);
		elementToHwPriority.put(ChemEl.Sn, 8);
		elementToHwPriority.put(ChemEl.Pb, 7);
		elementToHwPriority.put(ChemEl.B, 6);
		elementToHwPriority.put(ChemEl.Al, 5);
		elementToHwPriority.put(ChemEl.Ga, 4);
		elementToHwPriority.put(ChemEl.In, 3);
		elementToHwPriority.put(ChemEl.Tl, 2);
		elementToHwPriority.put(ChemEl.Hg, 1);
	}

	/**
	 * Useful to give an indication of whether a bond is like to be ionic (diff >1.8), polar or covalent (diff < 1.2)
	 * @param chemEl
	 * @return
	 */
	static Double getPaulingElectronegativity(ChemEl chemEl) {
		return elementToPaulingElectronegativity.get(chemEl);
	}

	/**
	 * Maps chemEl to the priority of that atom in Hantzch-Widman system. A higher value indicates a higher priority.
	 * @param chemEl
	 * @return
	 */
	static Integer getHwpriority(ChemEl chemEl) {
		return elementToHwPriority.get(chemEl);
	}
}
