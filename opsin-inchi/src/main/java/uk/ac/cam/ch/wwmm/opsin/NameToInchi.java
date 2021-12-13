package uk.ac.cam.ch.wwmm.opsin;

import java.util.HashMap;
import java.util.List;
import java.util.Set;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import io.github.dan2097.jnainchi.InchiAtom;
import io.github.dan2097.jnainchi.InchiBond;
import io.github.dan2097.jnainchi.InchiBondType;
import io.github.dan2097.jnainchi.InchiFlag;
import io.github.dan2097.jnainchi.InchiInput;
import io.github.dan2097.jnainchi.InchiKeyOutput;
import io.github.dan2097.jnainchi.InchiOptions;
import io.github.dan2097.jnainchi.InchiStereo;
import io.github.dan2097.jnainchi.InchiOptions.InchiOptionsBuilder;
import io.github.dan2097.jnainchi.InchiOutput;
import io.github.dan2097.jnainchi.InchiStatus;
import io.github.dan2097.jnainchi.InchiStereoParity;
import io.github.dan2097.jnainchi.JnaInchi;

import uk.ac.cam.ch.wwmm.opsin.BondStereo.BondStereoValue;

/**
 * Allows the conversion of OPSIN's output into (Std)InChIs or StdInChIKeys
 * Also can be used, as a convenience method, to directly convert chemical names to (Std)InChIs or StdInChIKeys
 * @author dl387
 *
 */
public class NameToInchi {

	private static final Logger LOG = LogManager.getLogger(NameToInchi.class);
	private NameToStructure n2s;
	public NameToInchi() {
		n2s = NameToStructure.getInstance();
	}

	/**Parses a chemical name, returning an InChI representation of the molecule.
	 *
	 * @param name The chemical name to parse.
	 * @return An InChI string, containing the parsed molecule, or null if the molecule would not parse.
	 */
	public String parseToInchi(String name) {
		OpsinResult result = n2s.parseChemicalName(name);
		return convertResultToInChI(result);
	}
	
	/**Parses a chemical name, returning a StdInChI representation of the molecule.
	 * Note that chemical names typically specify an exact tautomer which is not representable in StdInChI
	 * Use {@link #parseToInchi(String)} if you want to represent the exact tautomer using a fixed hydrogen layer
	 *
	 * @param name The chemical name to parse.
	 * @return A StdInChI string, containing the parsed molecule, or null if the molecule would not parse.
	 */
	public String parseToStdInchi(String name) {
		OpsinResult result = n2s.parseChemicalName(name);
		return convertResultToStdInChI(result);
	}
	
	/**Parses a chemical name, returning a StdInChIKey for the molecule.
	 * Like StdInChI, StdInChIKeys aim to not be tautomer specific
	 *
	 * @param name The chemical name to parse.
	 * @return A StdInChIKey string or null if the molecule would not parse.
	 */
	public String parseToStdInchiKey(String name) {
		OpsinResult result = n2s.parseChemicalName(name);
		return convertResultToStdInChIKey(result);
	}
	
	/**
	 * Converts an OPSIN result to InChI. Null is returned if this conversion fails
	 * @param result
	 * @return String InChI
	 */
	public static String convertResultToInChI(OpsinResult result){
		return convertResultToInChI(result, false);
	}
	
	/**
	 * Converts an OPSIN result to StdInChI. Null is returned if this conversion fails
	 * Note that chemical names typically specify an exact tautomer which is not representable in StdInChI
	 * Use {@link #convertResultToInChI(OpsinResult)} if you want to represent the exact tautomer using a fixed hydrogen layer
	 * @param result
	 * @return String InChI
	 */
	public static String convertResultToStdInChI(OpsinResult result){
		return convertResultToInChI(result, true);
	}
	
	/**
	 * Converts an OPSIN result to a StdInChIKey. Null is returned if this conversion fails
	 * Like StdInChI, StdInChIKeys aim to not be tautomer specific
	 * @param result
	 * @return String InChIKey
	 */
	public static String convertResultToStdInChIKey(OpsinResult result){
		String stdInchi = convertResultToInChI(result, true);
		if (stdInchi != null){
			try {
				InchiKeyOutput key = JnaInchi.inchiToInchiKey(stdInchi);
				return key.getInchiKey();
			} catch (Exception e) {
				if (LOG.isDebugEnabled()){
					LOG.debug(e.getMessage(), e);
				}
				return null;
			}
		}
		return null;
	}
	
	private static String convertResultToInChI(OpsinResult result, boolean produceStdInChI){
		if (result.getStructure() != null){
			String inchi = null;
			try{
				inchi = opsinFragmentToInchi(result.getStructure(), produceStdInChI);
			}
			catch (Exception e) {
				if (LOG.isDebugEnabled()){
					LOG.debug(e.getMessage(), e);
				}
				return null;
			}
			if (inchi ==null){
				//inchi generation failed
				return null;
			}
			if(LOG.isDebugEnabled()){
				LOG.debug(inchi);
			}
			return inchi;
		}
		return null;
	}

	private static String opsinFragmentToInchi(Fragment frag, boolean produceStdInChI) {
		HashMap<Integer, InchiAtom> opsinIdAtomMap = new HashMap<>();
		InchiOptionsBuilder optionsBuilder = new InchiOptions.InchiOptionsBuilder();
		optionsBuilder.withFlag(InchiFlag.AuxNone);

		if (!produceStdInChI){
			optionsBuilder.withFlag(InchiFlag.FixedH);
		}
		InchiInput input = new InchiInput();

		List<Atom> atomList =frag.getAtomList();
		// Generate atoms
		for (Atom atom : atomList) {
			InchiAtom inchiAtom = new InchiAtom(atom.getElement().toString());
			input.addAtom(inchiAtom);
			inchiAtom.setCharge(atom.getCharge());
			Integer isotope = atom.getIsotope();
			if (isotope != null) {
				inchiAtom.setIsotopicMass(isotope);
			}
			opsinIdAtomMap.put(atom.getID(), inchiAtom);
		}
		Set<Bond> bondList = frag.getBondSet();
		for (Bond bond : bondList) {
			input.addBond(new InchiBond(opsinIdAtomMap.get(bond.getFrom()), opsinIdAtomMap.get(bond.getTo()), InchiBondType.of((byte)bond.getOrder())));
		}

		for (Atom atom : atomList) {//add atomParities
			AtomParity atomParity = atom.getAtomParity();
			if (atomParity == null) {
				continue;
			}
			StereoGroup stereoGroupType = atomParity.getStereoGroup();
        	if ((stereoGroupType == StereoGroup.Rac || stereoGroupType == StereoGroup.Rel) &&
					countStereoGroup(atom) == 1) {
        		continue;
        	}
			Atom[] atomRefs4 = atomParity.getAtomRefs4();
			int[] atomRefs4AsInt = new int[4];
			for (int i = 0; i < atomRefs4.length; i++) {
				atomRefs4AsInt[i] = atomRefs4[i].getID();
			}
			InchiStereoParity parity = InchiStereoParity.UNKNOWN;
			if (atomParity.getParity() > 0){
				parity = InchiStereoParity.EVEN;
			}
			else if (atomParity.getParity() < 0){
				parity = InchiStereoParity.ODD;
			}
			input.addStereo(InchiStereo.createTetrahedralStereo(opsinIdAtomMap.get(atom.getID()), opsinIdAtomMap.get(atomRefs4AsInt[0]), opsinIdAtomMap.get(atomRefs4AsInt[1]), opsinIdAtomMap.get(atomRefs4AsInt[2]), opsinIdAtomMap.get(atomRefs4AsInt[3]), parity));
        }

		for (Bond bond : bondList) {//add bondStereos
			BondStereo bondStereo =bond.getBondStereo();
			if (bondStereo != null){
				Atom[] atomRefs4 = bondStereo.getAtomRefs4();
				int[] atomRefs4Ids = new int[4];
				for (int i = 0; i < atomRefs4.length; i++) {
					atomRefs4Ids[i] = atomRefs4[i].getID();
				}
				if (BondStereoValue.CIS.equals(bondStereo.getBondStereoValue())){
					input.addStereo(InchiStereo.createDoubleBondStereo(opsinIdAtomMap.get(atomRefs4Ids[0]), opsinIdAtomMap.get(atomRefs4Ids[1]), opsinIdAtomMap.get(atomRefs4Ids[2]), opsinIdAtomMap.get(atomRefs4Ids[3]), InchiStereoParity.ODD));
				}
				else if (BondStereoValue.TRANS.equals(bondStereo.getBondStereoValue())){
					input.addStereo(InchiStereo.createDoubleBondStereo(opsinIdAtomMap.get(atomRefs4Ids[0]), opsinIdAtomMap.get(atomRefs4Ids[1]), opsinIdAtomMap.get(atomRefs4Ids[2]), opsinIdAtomMap.get(atomRefs4Ids[3]), InchiStereoParity.EVEN));
				}
			}
        }
		InchiOutput output = JnaInchi.toInchi(input, optionsBuilder.build());
    	InchiStatus ret = output.getStatus();
    	if (LOG.isDebugEnabled()){
    		LOG.debug("Inchi generation status: " + ret);
    		if (InchiStatus.SUCCESS != ret){
    			LOG.debug(output.getMessage());
    		}
    	}
    	if (InchiStatus.SUCCESS != ret && InchiStatus.WARNING != ret) {
    		return null;
    	}
    	return output.getInchi();
	}

	private static int countStereoGroup(Atom atom) {
		int count = 0;
		for (Atom a : atom.getFrag().getAtomList()) {
			if (a.getAtomParity() == null)
				continue;
			if (a.getAtomParity().getStereoGroup().equals(atom.getAtomParity().getStereoGroup()) &&
					a.getAtomParity().getStereoGroupNum() == atom.getAtomParity().getStereoGroupNum())
				count++;
		}
		return count;
	}
}
