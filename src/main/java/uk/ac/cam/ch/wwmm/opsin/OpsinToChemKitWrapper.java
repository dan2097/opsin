package uk.ac.cam.ch.wwmm.opsin;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import nu.xom.Element;
import sea36.chem.base.ChemicalElement;
import sea36.chem.base.Order;
import sea36.chem.core.CMLAtom;
import sea36.chem.core.CMLAtomParity;
import sea36.chem.core.CMLBond;
import sea36.chem.core.CMLBondStereo;
import sea36.chem.core.CMLMolecule;
import sea36.chem.core.CMLAtomParity.Parity;
import sea36.chem.core.CMLBondStereo.Type;

class OpsinToChemKitWrapper {

	private Map<CMLAtom,Atom> chemKitAtomToOpsinAtomMap = new HashMap<CMLAtom, Atom>();
	private Map<Atom, CMLAtom> opsinAtomToChemKitAtomMap = new HashMap<Atom, CMLAtom>();
	private CMLMolecule chemKitMolecule =  new CMLMolecule();

	OpsinToChemKitWrapper(Fragment opsinFragment) {
		List<Atom> atoms = opsinFragment.getAtomList();
		for (Atom atom : atoms) {
			CMLAtom chemKitAtom = opsinAtomToChemKitAtom(atom);
			chemKitAtomToOpsinAtomMap.put(chemKitAtom, atom);
			opsinAtomToChemKitAtomMap.put(atom, chemKitAtom);
			chemKitMolecule.addAtom(chemKitAtom);
		}
		List<Bond> bonds = opsinFragment.getBondList();
		for (Bond bond : bonds) {
			CMLBond chemKitBond = opsinBondToChemKitBond(bond);
			chemKitMolecule.addAndGenerateId(chemKitBond);
		}
	}

	private CMLAtom opsinAtomToChemKitAtom(Atom atom) {
		CMLAtom chemKitAtom = new CMLAtom();
		ChemicalElement ce = ChemicalElement.valueOf(atom.getElement());
		if (ce !=null){
			chemKitAtom.setChemicalElement(ce);
		}
		else{
			chemKitAtom.setChemicalElement(atom.getElement());
		}
		chemKitAtom.setFormalCharge(atom.getCharge());
		chemKitAtom.setId("a" +atom.getID());
		if (atom.getAtomParityElement()!=null){
			Element originalAtomParity = atom.getAtomParityElement();
			CMLAtomParity cmlAtomParity = new CMLAtomParity();
			cmlAtomParity.setAtomRefs4(originalAtomParity.getAttributeValue("atomRefs4"));
			if (Float.parseFloat(originalAtomParity.getValue()) > 0.1){
				cmlAtomParity.setParity(Parity.EVEN);
			}
			else if (Float.parseFloat(originalAtomParity.getValue()) < -0.1){
				cmlAtomParity.setParity(Parity.ODD);
			}
			else{
				cmlAtomParity.setParity(Parity.UNKNOWN);
			}
			chemKitAtom.appendChild(cmlAtomParity);
		}
		return chemKitAtom;
	}

	private CMLBond opsinBondToChemKitBond(Bond bond) {
		CMLAtom atom1 = opsinAtomToChemKitAtomMap.get(bond.getFromAtom());
		CMLAtom atom2 = opsinAtomToChemKitAtomMap.get(bond.getToAtom());
		Order bondOrder;
		if (bond.getOrder()==1){
			bondOrder =Order.SINGLE;
		}
		else if (bond.getOrder()==2){
			bondOrder =Order.DOUBLE;
		}
		else if (bond.getOrder()==3){
			bondOrder =Order.TRIPLE;
		}
		else{
			bondOrder =Order.UNKNOWN;
		}
		CMLBond chemKitBond = new CMLBond(atom1, atom2, bondOrder);
		if (bond.getBondStereoElement()!=null){
			Element originalBondStereo = bond.getBondStereoElement();
			CMLBondStereo cmlBondStereo = new CMLBondStereo();
			cmlBondStereo.setAtomRefs4(originalBondStereo.getAttributeValue("atomRefs4"));
			if (originalBondStereo.getValue().equals("C")){
				cmlBondStereo.setType(Type.CIS);
			}
			else if (originalBondStereo.getValue().equals("T")){
				cmlBondStereo.setType(Type.TRANS);
			}
			else{
				cmlBondStereo.setType(Type.UNKNOWN);
			}
			chemKitBond.appendChild(cmlBondStereo);
		}
		return chemKitBond;
	}

	CMLMolecule getChemKitMolecule() {
		return chemKitMolecule;
	}

	Atom getOpsinAtomFromChemKitAtom(CMLAtom a) {
		return chemKitAtomToOpsinAtomMap.get(a);
	}

	CMLAtom getChemKitAtomFromOpsinAtom(Atom a) {
		return opsinAtomToChemKitAtomMap.get(a);
	}
}
