package uk.ac.cam.ch.wwmm.opsin;

import java.util.*;
import java.util.regex.Pattern;

import uk.ac.cam.ch.wwmm.opsin.Bond.SMILES_BOND_DIRECTION;

import nu.xom.Attribute;
import nu.xom.Element;

/** A builder for fragments specified as SMILES. A custom SMILES dialect is used. With the exception of | to directly set valency and the inclusion of sb/te as being aromatic this is a subset of real SMILES:
 *
 * Allowed:
 * Organic elements B,C,N,O,P,S,F,Cl,Br,I (square brackets not required)
 * Aromatic elements c,n,o,p,s (square brackets not required) as,se,sb,te (square brackets required) Note that the inclusion of sb/te are an unofficial extension
 * =, # for bond orders
 * . for disconnection
 * (, ) for branching
 * [, ] for placing inorganic elements within and specifying charge. Allowed: [Al3+] or [Al+++]
 * 012345679 - ring closures
 * %10 %99 - more ring closures (%100 is ring closure %10 and 0 as in normal SMILES)
 * / and \ to set double bond stereochemistry to cis/trans
 * @ and @@ to set tetrahedral stereochemistry as in SMILES. Note that ONLY in this context is a capital H allowed, in all other cases hydrogen should be implicit
 *
 * |3 |5 etc. can be used to set the valency of an atom e.g. P|5, [Se|2]
 * (| looks slightly like an l as in the lambda convention)
 *
 * Also, an = or # at the start of the string indicates that the group attaches to its parent group via a double or triple bond.
 *
 * A -,=,# on the end indicates that in the absence of locants, other groups attach to
 * *it* via the atom at the end of the string, not at the start of the string with -,=,# meaning single,double or triple bond
 * This behaviour is overridden for certain suffixes to give different meanings to the atom the -,=,# is referring to
 *
 * @author ptc24/dl387
 *
 */
class SMILESFragmentBuilder {

	/**A "struct" to hold information on the parsing stack
	 *
	 * @author ptc24
	 *
	 */
	private class StackFrame {
		/**The Atom currently under consideration.*/
		Atom atom;
		/**The order of the bond about to be formed.*/
		int bondOrder;
		/**Whether the bond is a \ or / bond for use in determining cis/trans.*/
		SMILES_BOND_DIRECTION slash;

		/**Creates a stack frame with given parameters.
		 *
		 * @param a An atom or null
		 * @param bondOrderVal The value for bondOrder.
		 */
		StackFrame(Atom a, int bondOrderVal) {
			atom = a;
			bondOrder = bondOrderVal;
			slash = null;
		}

		/**Creates a copy of an existing StackFrame.
		 *
		 * @param sf The stackframe to copy.
		 */
		StackFrame(StackFrame sf) {
			atom = sf.atom;
			bondOrder = sf.bondOrder;
		}
	}

	/**Upper case letters.*/
	private static String upperLetters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
	/**Lower case letters.*/
	private static String lowerLetters = "abcdefghijklmnopqrstuvwxyz";
	/**Numerical digits.*/
	private static String digits = "0123456789";
	/**Organic Atoms.*/
	private static Set<String> organicAtoms = new HashSet<String>();
	/**Aromatic Atoms.*/
	private static Set<String> aromaticAtoms = new HashSet<String>();

	static {
		organicAtoms.add("B");
		organicAtoms.add("C");
		organicAtoms.add("N");
		organicAtoms.add("O");
		organicAtoms.add("P");
		organicAtoms.add("S");
		organicAtoms.add("F");
		organicAtoms.add("Cl");
		organicAtoms.add("Br");
		organicAtoms.add("I");

		aromaticAtoms.add("c");
		aromaticAtoms.add("n");
		aromaticAtoms.add("o");
		aromaticAtoms.add("p");
		aromaticAtoms.add("s");
		aromaticAtoms.add("as");
		aromaticAtoms.add("se");
		aromaticAtoms.add("sb");
		aromaticAtoms.add("te");
	}

	private static Pattern matchComma =Pattern.compile(",");
	private static Pattern matchSlash =Pattern.compile("/");


	/**Build a Fragment based on a SMILES string, with a null type/subType.
	 *
	 * @param smiles The SMILES string to build from.
	 * @return The built fragment.
	 * @throws StructureBuildingException
	 */
	Fragment build(String smiles, FragmentManager fragManaager) throws StructureBuildingException {
		return build(smiles, "", "", "", fragManaager);
	}

	/**
	 * Build a Fragment based on a SMILES string.
	 * @param smiles The SMILES string to build from.
	 * @param type The type of fragment being built.
	 * @param subType The subtype of fragment being built.
	 * @param labelMapping A string indicating which locants to assign to each atom. Can be a slash delimited list, "" for default numbering or "none"
	 * @param fragManager
	 * @return Fragment The built fragment.
	 * @throws StructureBuildingException
	 */
	Fragment build(String smiles, String type, String subType, String labelMapping, FragmentManager fragManager) throws StructureBuildingException {
		if (smiles==null){
			throw new StructureBuildingException("SMILES specified is null");
		}
		if (type==null){
			throw new StructureBuildingException("type specified is null, use \"\" if a type is not desired ");
		}
		if (subType==null){
			throw new StructureBuildingException("subType specified is null, use \"\" if a subType is not desired ");
		}
		if (labelMapping==null){
			throw new StructureBuildingException("labelMapping is null use \"none\" if you do not want any numbering or \"\" if you would like default numbering");
		}
		List<String> labelMap = null;
		if(!labelMapping.equals("none") && !labelMapping.equals("fusedRing") ) {
			labelMap = new ArrayList<String>();
			String [] mappingTmp = matchSlash.split(labelMapping, -1);
            labelMap.addAll(Arrays.asList(mappingTmp));//place slash delimited labels into arrayList
		}
		int currentNumber = 1;
		Fragment currentFrag = new Fragment(type, subType);
		Stack<StackFrame> stack = new Stack<StackFrame>();
		stack.push(new StackFrame(null, 1));
		HashMap<String, StackFrame> closures = new HashMap<String, StackFrame>();//used for ring closures
		String tmpString = smiles;
		char firstCharacter =tmpString.charAt(0);
		if(firstCharacter == '=' || firstCharacter == '#') {//used by OPSIN to specify the valency with which this fragment connects (1 if not specified)
			tmpString = tmpString.substring(1);
		}
		char lastCharacter =tmpString.charAt(tmpString.length()-1);
		if(lastCharacter == '-' || lastCharacter == '=' || lastCharacter == '#') {//used by OPSIN to specify the valency with which this fragment connects and to indicate it connects via the last atom in the SMILES
			tmpString = tmpString.substring(0, tmpString.length()-1);
		}

		while(tmpString.length() > 0) {
			String nextChar = tmpString.substring(0, 1);
			tmpString = tmpString.substring(1);
			if(nextChar.equals("(")) {
				stack.push(new StackFrame(stack.peek()));
			} else if(nextChar.equals(")")) {
				stack.pop();
			} else if(nextChar.equals("=")){
				stack.peek().bondOrder = 2;
			} else if(nextChar.equals("#")){
				stack.peek().bondOrder = 3;
			} else if(nextChar.equals("/")){
				stack.peek().slash = SMILES_BOND_DIRECTION.RSLASH;
			} else if(nextChar.equals("\\")){
				stack.peek().slash = SMILES_BOND_DIRECTION.LSLASH;
			} else if(nextChar.equals(".")){
				stack.peek().atom = null;
			} else if(upperLetters.contains(nextChar) || lowerLetters.contains(nextChar)) {//organic atoms
		        String elementType = nextChar;
		        boolean spareValency =false;
		        if(upperLetters.contains(elementType)) {//normal atoms
					if(tmpString.length() > 0 && lowerLetters.contains(tmpString.substring(0,1)) && organicAtoms.contains(elementType + tmpString.substring(0,1))) {
						elementType += tmpString.substring(0,1);
						tmpString = tmpString.substring(1);
					}
					else if (!organicAtoms.contains(elementType)){
						throw new StructureBuildingException(elementType +" is not an organic Element. If it is actually an element it should be in square brackets");
					}
		        }
		        else if(lowerLetters.contains(elementType)) {//aromatic atoms
					if (!aromaticAtoms.contains(elementType)){
						throw new StructureBuildingException(elementType +" is not an aromatic Element. If it is actually an element it should not be in lower case");
					}
					elementType = elementType.toUpperCase();
					spareValency =true;
		        }
				Atom atom = fragManager.createAtom(elementType, currentFrag);
				atom.setSpareValency(spareValency);
				if(labelMapping.equals("")) {
					atom.addLocant(Integer.toString(currentNumber));
				} else if (labelMap !=null){
					String labels[] = matchComma.split(labelMap.get(currentNumber-1));
                    for (String label : labels) {
                        if (!label.equals("")) {
                            atom.addLocant(label);
                        }
                    }
				}
				currentFrag.addAtom(atom);
				if(stack.peek().atom !=null) {
					Bond b = fragManager.createBond(stack.peek().atom, atom, stack.peek().bondOrder);
					if (stack.peek().slash!=null){
						b.setSmilesStereochemistry(stack.peek().slash);
						stack.peek().slash = null;
					}
					if (stack.peek().atom.getAtomParityElement()!=null){
						Element atomParityElement = stack.peek().atom.getAtomParityElement();
						Attribute atomRefs4Att = atomParityElement.getAttribute("atomRefs4");
						String atomRefs4 = atomRefs4Att.getValue();
						if (atomRefs4.equals("")){
							atomRefs4 += "a" +atom.getID();
						}
						else{
							atomRefs4 += " a" +atom.getID();
						}
						atomRefs4Att.setValue(atomRefs4);
					}
				}
				stack.peek().atom = atom;
				stack.peek().bondOrder = 1;
				currentNumber += 1;
			} else if(nextChar.equals("[")) {//square brackets- contain non-organic atoms and are used to unambiguously set charge/chirality etc.
				int indexOfRightSquareBracket = tmpString.indexOf("]");
                if (indexOfRightSquareBracket == -1) {
                    throw new StructureBuildingException("[ without matching \"]\"");
                }
                String atomString = tmpString.substring(0, indexOfRightSquareBracket);//the contents of the square bracket
                tmpString = tmpString.substring(indexOfRightSquareBracket +1);
             // isotope
                String isotope = "";
				while(atomString.length() > 0 &&
						digits.contains(atomString.substring(0,1))) {
					isotope += atomString.substring(0,1);
					atomString = atomString.substring(1);
				}
		        if (!isotope.equals("")){
		        	throw new StructureBuildingException("Isotopically labelled atoms in SMILES are not yet supported");
		        }

		        if (atomString.length() > 0){
		        	nextChar = atomString.substring(0,1);
		        	atomString = atomString.substring(1);
		        }
		        else{
		        	throw new StructureBuildingException("No element found in square brackets");
		        }
		// elementType
		        String elementType = nextChar;
		        boolean spareValency = false;
		        if(upperLetters.contains(elementType)) {//normal atoms
					if(atomString.length() > 0 && lowerLetters.contains(atomString.substring(0,1))) {
						elementType += atomString.substring(0,1);
						atomString = atomString.substring(1);
					}
		        }
		        else if(lowerLetters.contains(elementType)) {//aromatic atoms
					if(atomString.length() > 0 && lowerLetters.contains(atomString.substring(0,1))) {
						if (aromaticAtoms.contains(elementType + atomString.substring(0,1))){
							elementType = elementType.toUpperCase() + atomString.substring(0,1);
							atomString = atomString.substring(1);
						}
						else{
							throw new StructureBuildingException(elementType + atomString.substring(0,1) +" is not an aromatic Element. If it is actually an element it should not be in lower case");
						}
					}
					else{
						if (!aromaticAtoms.contains(elementType)){
							throw new StructureBuildingException(elementType +" is not an aromatic Element.");
						}
						elementType = elementType.toUpperCase();
					}
					spareValency =true;
		        }
		        else{
		        	throw new StructureBuildingException(elementType +" is not a valid element type!");
		        }
				Atom atom = fragManager.createAtom(elementType, currentFrag);
				atom.setSpareValency(spareValency);
				if(labelMapping.equals("")) {
					atom.addLocant(Integer.toString(currentNumber));
				} else if (labelMap !=null){
					String labels[] = matchComma.split(labelMap.get(currentNumber-1));
                    for (String label : labels) {
                        if (!label.equals("")) {
                            atom.addLocant(label);
                        }
                    }
				}
				currentFrag.addAtom(atom);
				if(stack.peek().atom != null) {
					Bond b = fragManager.createBond(stack.peek().atom, atom, stack.peek().bondOrder);
					if (stack.peek().slash!=null){
						b.setSmilesStereochemistry(stack.peek().slash);
						stack.peek().slash = null;
					}
					if (stack.peek().atom.getAtomParityElement()!=null){
						Element atomParityElement = stack.peek().atom.getAtomParityElement();
						Attribute atomRefs4Att = atomParityElement.getAttribute("atomRefs4");
						String atomRefs4 = atomRefs4Att.getValue();
						if (atomRefs4.equals("")){
							atomRefs4 += "a" +atom.getID();
						}
						else{
							atomRefs4 += " a" +atom.getID();
						}
						atomRefs4Att.setValue(atomRefs4);
					}
				}
				Atom previousAtom = stack.peek().atom;//needed for setting atomParity elements up
				stack.peek().atom = atom;
				stack.peek().bondOrder = 1;
				currentNumber += 1;

		        int hydrogenCount =0;
		        int charge = 0;
		        Boolean chiralityClockwise = null;
		        while (atomString.length()>0){
		        	nextChar = atomString.substring(0,1);
		        	atomString = atomString.substring(1);
		        	if(nextChar.equals("@")) {// chirality-sets atom parity
		        		if (chiralityClockwise != null){
		        			throw new StructureBuildingException("Atom parity appeared to be specified twice for an atom in a square bracket!");
		        		}
						Element atomParity = new Element("atomParity");
						atom.setAtomParityElement(atomParity);
						chiralityClockwise = false;
						if (atomString.length() > 0 && atomString.substring(0,1).equals("@")){
							chiralityClockwise = true;
							atomString = atomString.substring(1);
						}
						if (chiralityClockwise){
							atomParity.appendChild("1");
						}
						else{
							atomParity.appendChild("-1");
						}
						String atomRefs4 ="";
						if (previousAtom !=null){
							atomRefs4 += "a" +previousAtom.getID();
						}
						if (atomString.length() > 0 && atomString.substring(0,1).equals("H")){
							if (!atomRefs4.equals("")){
								atomRefs4 += " ";
							}
							atomRefs4 += "a" + atom.getID() +"_H";
							atomString = atomString.substring(1);
						}
						atomParity.addAttribute(new Attribute("atomRefs4", atomRefs4));
		 			}
		            else if (nextChar.equals("H")){// hydrogenCount
		            	if (hydrogenCount != 0){
		            		throw new StructureBuildingException("Hydrogen count appeared to be specified twice for an atom in a square bracket!");
		            	}
		            	String hydrogenCountString ="";
    					while(atomString.length() > 0 && digits.contains(atomString.substring(0,1))) {
    						hydrogenCountString += atomString.substring(0,1);
    						atomString = atomString.substring(1);
    					}
    					if (hydrogenCountString.equals("")){
    						hydrogenCount=1;
    					}
    					else{
    						hydrogenCount = Integer.parseInt(hydrogenCountString);
    					}
    					throw new StructureBuildingException("Hydrogen count is currently not supported");
		            }
		            else if(nextChar.equals("+") || nextChar.equals("-")) {// formalCharge
		            	if (charge != 0){
		            		throw new StructureBuildingException("Charge appeared to be specified twice for an atom in a square bracket!");
		            	}
	    				charge = nextChar.equals("+") ? 1 : -1;
    					String changeChargeStr = "";
    					int changeCharge = 1;
    					while(atomString.length() > 0 && digits.contains(atomString.substring(0,1))) {//e.g. [C+2]
    						changeChargeStr+= atomString.substring(0,1);
    						atomString = atomString.substring(1);
    					}
    					if (changeChargeStr.equals("")){
    						while(atomString.length() > 0){//e.g. [C++]
    							nextChar = atomString.substring(0,1);
    							if (nextChar.equals("+")){
    								if (charge != 1){
    									throw new StructureBuildingException("Atom has both positive and negative charges specified!");//e.g. [C+-]
    								}
    							}
    							else if (nextChar.equals("-")){
    								if (charge != -1){
    									throw new StructureBuildingException("Atom has both negative and positive charges specified!");
    								}
    							}
    							else{
    								break;
    							}
    							changeCharge++;
        						atomString = atomString.substring(1);
        					}
    					}
    					changeCharge = changeChargeStr.equals("") ? changeCharge : Integer.parseInt(changeChargeStr);
    					atom.setCharge(atom.getCharge() + (charge * changeCharge) );
		            }
		            else if(nextChar.equals("|")) {
						String lambda = "";
						while(atomString.length() > 0 &&
								digits.contains(atomString.substring(0,1))) {
							lambda += atomString.substring(0,1);
							atomString = atomString.substring(1);
						}
						atom.setLambdaConventionValency(Integer.parseInt(lambda));
					}
		            else{
		            	throw new StructureBuildingException("Unexpected character found in square bracket");
		            }
		        }
			} else if(digits.contains(nextChar) ||
					nextChar.equals("%")) {
				if(nextChar.equals("%")) {
					nextChar = "";
					if (tmpString.length() >=2 && digits.contains(tmpString.substring(0,1)) && digits.contains(tmpString.substring(1,2))) {
						nextChar += tmpString.substring(0,2);
						tmpString = tmpString.substring(2);
					}
					else{
						throw new StructureBuildingException("A ring opening indice after a % must be two digits long");
					}
				}
				if(closures.containsKey(nextChar)) {
					StackFrame sf = closures.remove(nextChar);
					int bondOrder = 1;
					if(sf.bondOrder > 1) {
						if(stack.peek().bondOrder > 1) {
							throw new StructureBuildingException("ring closure should not have bond order specified twice!");
						}
						bondOrder = sf.bondOrder;
					} else if(stack.peek().bondOrder > 1) {
						bondOrder = stack.peek().bondOrder;
					}
					Bond b;
					if (stack.peek().slash ==null){
						b = fragManager.createBond(sf.atom, stack.peek().atom, bondOrder);
					}
					else{
						b = fragManager.createBond(stack.peek().atom, sf.atom, bondOrder);//special case e.g. CC1=C/F.O\1  Bond is done from the O to the the C due to the presence of the \
					}
					if(sf.slash !=null) {
						if(stack.peek().slash !=null) {
							throw new StructureBuildingException("ring closure should not have cis/trans specified twice!");
						}
						b.setSmilesStereochemistry(sf.slash);
					} else if(stack.peek().slash !=null) {
						b.setSmilesStereochemistry(stack.peek().slash);
						stack.peek().slash = null;
					}
					if (stack.peek().atom.getAtomParityElement()!=null){
						Element atomParityElement = stack.peek().atom.getAtomParityElement();
						Attribute atomRefs4Att = atomParityElement.getAttribute("atomRefs4");
						String atomRefs4 = atomRefs4Att.getValue();
						if (atomRefs4.equals("")){
							atomRefs4 += "a" + sf.atom.getID();
						}
						else{
							atomRefs4 += " a" + sf.atom.getID();
						}
						atomRefs4Att.setValue(atomRefs4);
					}
					if (sf.atom.getAtomParityElement()!=null){//replace dummy id with actual id e.g. N[C@@H]1C.F1 where the 1 initially holds a dummy id before being replaced with the id of the F
						Element atomParityElement = sf.atom.getAtomParityElement();
						Attribute atomRefs4Atr = atomParityElement.getAttribute("atomRefs4");
						String atomRefs4 = atomRefs4Atr.getValue();
						atomRefs4 = atomRefs4.replaceFirst("ringclosure" + nextChar, "a" +stack.peek().atom.getID());
						atomRefs4Atr.setValue(atomRefs4);
					}
					stack.peek().bondOrder = 1;
				} else {
					StackFrame sf = new StackFrame(stack.peek());
					if (stack.peek().slash!=null){
						sf.slash = stack.peek().slash;
						stack.peek().slash = null;
					}
					if (sf.atom.getAtomParityElement()!=null){//replace ringclosureX with actual reference to id when it is known
						Element atomParityElement = sf.atom.getAtomParityElement();
						Attribute atomRefs4Att = atomParityElement.getAttribute("atomRefs4");
						String atomRefs4 = atomRefs4Att.getValue();
						if (atomRefs4.equals("")){
							atomRefs4 += "ringclosure" + nextChar;
						}
						else{
							atomRefs4 += " ringclosure" + nextChar;
						}
						atomRefs4Att.setValue(atomRefs4);
					}
					closures.put(nextChar, sf);
					stack.peek().bondOrder = 1;
				}
			}else if(nextChar.equals("|")) {
				String lambda = "";
				while(tmpString.length() > 0 &&
						digits.contains(tmpString.substring(0,1))) {
					lambda += tmpString.substring(0,1);
					tmpString = tmpString.substring(1);
				}
				if(stack.peek().atom !=null) {
					Atom a = stack.peek().atom;
					a.setLambdaConventionValency(Integer.parseInt(lambda));
				}
				else{
					throw new StructureBuildingException("| found in SMILES string at unexpected position");
				}
			}
			else{
				throw new StructureBuildingException(nextChar + " is in an unexpected position. Check this is not a mistake and that this feature of SMILES is supported by OPSIN's SMILES parser");
			}
		}
		if (labelMap != null && labelMap.size() >= currentNumber ){
			throw new StructureBuildingException("Group numbering has been invalidly defined in resource file: labels: " +labelMap.size() + ", atoms: " + (currentNumber -1) );
		}

		if(labelMapping.equals("fusedRing")) {//fragment is a fusedring with atoms in the correct order for fused ring numbering
			//this will do stuff like changing labels from 1,2,3,4,5,6,7,8,9,10->1,2,3,4,4a,5,6,7,8,8a
			FragmentTools.relabelFusedRingSystem(currentFrag);
		}
		Set<Bond> bonds = currentFrag.getBondSet();
		mainLoop: for (Bond centralBond : bonds) {//identify cases of E/Z stereochemistry and add appropriate bondstereo tags
			if (centralBond.getOrder()==2){
				Set<Bond> fromAtomBonds =centralBond.getFromAtom().getBonds();
				for (Bond preceedingBond : fromAtomBonds) {
					if (preceedingBond.getSmilesStereochemistry()!=null){
						Set<Bond> toAtomBonds = centralBond.getToAtom().getBonds();
						for (Bond followingBond : toAtomBonds) {
							if (followingBond.getSmilesStereochemistry()!=null){//now found a double bond surrounded by two bonds with slashs
								Element bondStereoEl = new Element("bondStereo");
								Atom atom2 = centralBond.getFromAtom();
								Atom atom1;
								if (atom2 == preceedingBond.getToAtom()){
									atom1 = preceedingBond.getFromAtom();
								}
								else{
									atom1 = preceedingBond.getToAtom();
								}
								Atom atom3 = centralBond.getToAtom();
								Atom atom4;
								if (atom3 == followingBond.getFromAtom()){
									atom4 = followingBond.getToAtom();
								}
								else{
									atom4 = followingBond.getFromAtom();
								}
								bondStereoEl.addAttribute(new Attribute("atomRefs4", "a" + atom1.getID() +" " + "a" + atom2.getID() + " " + "a" + atom3.getID() +" " + "a" + atom4.getID()));
								Boolean upFirst;
								if (preceedingBond.getSmilesStereochemistry() == SMILES_BOND_DIRECTION.LSLASH){
									upFirst = preceedingBond.getToAtom() == atom2;//in normally constructed SMILES this will be the case but you could write C(/F)=C/F instead of F\C=C/F
								}
								else if (preceedingBond.getSmilesStereochemistry() == SMILES_BOND_DIRECTION.RSLASH){
									upFirst = preceedingBond.getToAtom() != atom2;
								}
								else{
									throw new StructureBuildingException(preceedingBond.getSmilesStereochemistry() + " is not a slash!");
								}

								Boolean upSecond = null;
								if (followingBond.getSmilesStereochemistry() == SMILES_BOND_DIRECTION.LSLASH){
									upSecond = followingBond.getFromAtom() != atom3;
								}
								else if (followingBond.getSmilesStereochemistry() == SMILES_BOND_DIRECTION.RSLASH){
									upSecond = followingBond.getFromAtom() == atom3;
								}
								else{
									throw new StructureBuildingException(followingBond.getSmilesStereochemistry() + " is not a slash!");
								}

								if (upFirst == upSecond){
									bondStereoEl.appendChild("C");
								}
								else{
									bondStereoEl.appendChild("T");
								}
								centralBond.setBondStereoElement(bondStereoEl);
								continue mainLoop;
							}
						}
					}

				}
			}
		}
		for (Bond bond : bonds) {
			bond.setSmilesStereochemistry(null);
		}

		if(lastCharacter == '-' || lastCharacter == '=' || lastCharacter == '#') {
			List<Atom> aList =currentFrag.getAtomList();
			int lastAtomID =aList.get(aList.size()-1).getID();
			if (subType.equals("inSuffix")){
				currentFrag.setDefaultInID(lastAtomID);
			}
			else{
				if (lastCharacter == '#'){
					currentFrag.addOutID(lastAtomID,3, true);
				}
				else if (lastCharacter == '='){
					currentFrag.addOutID(lastAtomID,2, true);
				}
				else{
					currentFrag.addOutID(lastAtomID,1, true);
				}
			}
		}

		if(firstCharacter == '='){
			currentFrag.addOutID(currentFrag.getIdOfFirstAtom(),2, true);
		}
		if (firstCharacter == '#'){
			currentFrag.addOutID(currentFrag.getIdOfFirstAtom(),3, true);
		}

		return currentFrag;
	}

}
