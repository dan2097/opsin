package uk.ac.cam.ch.wwmm.opsin;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import uk.ac.cam.ch.wwmm.ptclib.xml.XOMTools;

import nu.xom.Element;
import nu.xom.Elements;

/**
 * 
 * @author aa593/dl387
 *
 */
public class FusedRingBuilder {


	/**
	 * Sorts by atomSequences by the IUPAC rules for determining the preferred labelling
	 * The most preferred will be sorted to the back (0th position)
	 * @author dl387
	 *
	 */
	class SortAtomSequences implements Comparator<ArrayList<Atom>> {

	    public int compare(ArrayList<Atom> sequenceA, ArrayList<Atom> sequenceB){
	    	if (sequenceA.size() != sequenceB.size()){
	    		//Error in fused ring building. Identified ring sequences not the same lengths!
	    		return 0;
	    	}

	    	int i=0;
	    	int j=0;
	    	//Give low numbers for the heteroatoms as a set.
	    	while(i < sequenceA.size()){
				Atom atomA=sequenceA.get(i);
				boolean isAaHeteroatom =!atomA.getElement().equals("C");


				//bridgehead carbon do not increment numbering
				if (isAaHeteroatom ==false && atomA.getIncomingValency()>=3){
					i++;
					continue;
				}
				
				Atom atomB=sequenceB.get(j);
				boolean isBaHeteroatom =!atomB.getElement().equals("C");
				if (isBaHeteroatom ==false && atomB.getIncomingValency()>=3){
					j++;
					continue;
				}

				if (isAaHeteroatom ==true && isBaHeteroatom ==false){
					return -1;
				}
				if (isBaHeteroatom ==true && isAaHeteroatom ==false){
					return 1;
				}
	    		i++;j++;
	    	}

	    	i=0;
	    	j=0;
	    	//Give low numbers for heteroatoms when considered in the order: O, S, Se, Te, N, P, As, Sb, Bi, Si, Ge, Sn, Pb, B, Hg
	    	while(i < sequenceA.size()){
				Atom atomA=sequenceA.get(i);

				//bridgehead carbon do not increment numbering
				if (atomA.getElement().equals("C")&& atomA.getIncomingValency()>=3){
					i++;
					continue;
				}
				
				Atom atomB=sequenceB.get(j);
				if (atomB.getElement().equals("C") && atomB.getIncomingValency()>=3){
					j++;
					continue;
				}

				int atomAElementValue, atomBElementValue;
				if (heteroAtomValues.containsKey(atomA.getElement())){
					atomAElementValue = heteroAtomValues.get(atomA.getElement());
				}
				else{
					atomAElementValue=0;
				}
				if (heteroAtomValues.containsKey(atomB.getElement())){
					atomBElementValue = heteroAtomValues.get(atomB.getElement());
				}
				else{
					atomBElementValue=0;
				}
				if (atomAElementValue > atomBElementValue){
					return -1;
				}
				if (atomAElementValue < atomBElementValue){
					return 1;
				}
				i++;j++;
	    	}

	    	//Give low numbers to fusion carbon atoms.
	    	for ( i = 0; i < sequenceA.size(); i++) {
				Atom atomA=sequenceA.get(i);
				Atom atomB=sequenceB.get(i);
				if (atomA.getIncomingValency()>=3 && atomA.getElement().equals("C")){
					if (!(atomB.getIncomingValency()>=3 && atomB.getElement().equals("C"))){
						return -1;
					}
				}
				if (atomB.getIncomingValency()>=3 && atomB.getElement().equals("C")){
					if (!(atomA.getIncomingValency()>=3 && atomA.getElement().equals("C"))){
						return 1;
					}
				}
			}
	    	//Note that any sequences still unsorted at this step will have fusion carbon atoms in the same places
	    	//which means you can go through both sequences without constantly looking for fusion carbons i.e. the variable j is no longer needed

	    	//Give low numbers to fusion rather than non-fusion atoms of the same heteroelement.
	    	for (i = 0; i < sequenceA.size(); i++) {
				Atom atomA=sequenceA.get(i);
				Atom atomB=sequenceB.get(i);
				if (atomA.getIncomingValency()>=3){
					if (!(atomB.getIncomingValency()>=3)){
						return -1;
					}
				}
				if (atomB.getIncomingValency()>=3){
					if (!(atomA.getIncomingValency()>=3)){
						return 1;
					}
				}
			}
	    	return 0;
	    }
	}


	/**
	 * Sorts by number, then by letter e.g. 4,3,3b,5,3a,2 -->2,3,3a,3b,4,5
	 * @author dl387
	 *
	 */
	class SortLocants implements Comparator<Atom> {

	    public int compare(Atom atoma, Atom atomb){
	    	String locanta =atoma.getFirstLocant();
	    	String locantb =atomb.getFirstLocant();

	    	Pattern matchdigits = Pattern.compile("(\\d+).*");

	    	Matcher m1  =matchdigits.matcher(locanta);
	    	int locantaNumber=0;
	    	if (m1.matches()){
	    		locantaNumber=Integer.parseInt(m1.group(1));
	    	}
	    	else{
	    		return 0;//invalid locant (could be intentionally invalid)
	    	}

	    	Matcher m2  =matchdigits.matcher(locantb);
	    	int locantbNumber=0;
	    	if (m2.matches()){
	    		locantbNumber=Integer.parseInt(m2.group(1));
	    	}
	    	else{
	    		return 0;//invalid locant (could be intentionally invalid)
	    	}

	        if (locantaNumber >locantbNumber) {
	            return 1;//e.g. 3 vs 2 or 3a vs 2
	        } else if (locantbNumber >locantaNumber) {
	            return -1;//e.g. 2 vs 3 or 2 vs 3a
	        }
	        else{
	        	Pattern matchletters = Pattern.compile(".*([a-z]+)");
	        	m1  =matchletters.matcher(locanta);
	        	String locantaLetter="";
	        	if (m1.matches()){
	        		locantaLetter=m1.group(1);
	        	}
	        	else{
	        		return -1;// e.g. 1 vs 1a
	        	}

	        	m2  =matchletters.matcher(locantb);
	        	String locantbLetter="";
	        	if (m2.matches()){
	        		locantbLetter=m2.group(1);
	        	}
	        	else{
	        		return 1;//e.g. 1a vs 1
	        	}

	            if (locantaLetter.compareTo(locantbLetter)>=1) {
	                return 1;//e.g. 1b vs 1a
	            } else if (locantbLetter.compareTo(locantaLetter)>=1) {
	                return -1;//e.g. 1a vs 1b
	            }
	            return 0;
	        }
	    }
	}

	private HashMap<String, Integer> heteroAtomValues =new HashMap<String, Integer>();
	private Pattern matchSlash = Pattern.compile("/");

	FusedRingBuilder() {
		//unknown heteroatoms or carbon are given a value of 0
		heteroAtomValues.put("Hg",2);
		heteroAtomValues.put("B",3);
		heteroAtomValues.put("Pb",4);
		heteroAtomValues.put("Sn",5);
		heteroAtomValues.put("Ge",6);
		heteroAtomValues.put("Si",7);
		heteroAtomValues.put("Bi",8);
		heteroAtomValues.put("Sb",9);
		heteroAtomValues.put("As",10);
		heteroAtomValues.put("P",12);
		heteroAtomValues.put("N",13);
		heteroAtomValues.put("Te",14);
		heteroAtomValues.put("Se",15);
		heteroAtomValues.put("S",16);
		heteroAtomValues.put("O",17);
	}

	/**
	 * Master method for processing fused rings. If 2 groups are present will attempt to fuse them
	 * Returns the substituent/root with the 2 groups fused together into 1 group
	 * @param state: contains the current id and fragment manager
	 * @param e Element (substituent or root)
	 * @throws PostProcessingException
	 * @throws StructureBuildingException
	 */
	Element processFusedRings(BuildState state, Element e) throws PostProcessingException, StructureBuildingException {
		Elements groups =e.getChildElements("group");
		if (groups.size() < 2){return e;}//nothing to fuse
		if (groups.size() >2){
			throw new PostProcessingException("Unsupported fused ring system: more than 2 groups");
		}
		Elements fusions =e.getChildElements("fusion");
		if (fusions.size() >1){
			throw new PostProcessingException("Unsupported fused ring system: more than 1 fusion");
		}
		ArrayList<Fragment> rings = new ArrayList<Fragment>();
		HashMap<Fragment,HashMap<String, ArrayList<String>>> ringsInformation = new HashMap<Fragment,HashMap<String, ArrayList<String>>>();

		/*
		 * Resolve the Groups
		 */
		for(int i=0;i<groups.size();i++) {
			Element group=groups.get(i);
			Fragment ring =state.xmlFragmentMap.get(group);
			if (group.getAttribute("fusedRingNumbering")!=null && i==groups.size()-1){
				String[] standardNumbering = matchSlash.split(group.getAttributeValue("fusedRingNumbering"),-1);
				List<Atom> atomList =ring.getAtomList();
				for (int j = 0; j < standardNumbering.length; j++) {
					atomList.get(j).replaceLocant(standardNumbering[j]);
				}
			}
			
			rings.add(ring);
		}

		int numbersMissing=0;
		int lettersMissing=0;
		ArrayList<String> numericalLocantsOfChild = new ArrayList<String>();
		ArrayList<String> letterLocantsOfParent = new ArrayList<String>();
		if (fusions.size() > 0){
			Element fusion = fusions.get(0);
			String[] fusionArray=fusion.getValue().split("-");
			if (fusionArray.length ==2){
				String[] locantsOfChildTemp = fusionArray[0].replaceFirst("\\[", "").split(",");
				for (int i = 0; i < locantsOfChildTemp.length; i++) {
					numericalLocantsOfChild.add(locantsOfChildTemp[i]);
				}
				char[] tempLetterLocantsOfParent = fusionArray[1].replaceFirst("\\]", "").toCharArray();
				for (int i = 0; i < tempLetterLocantsOfParent.length; i++) {
					letterLocantsOfParent.add(String.valueOf(tempLetterLocantsOfParent[i]));
				}
			}
			else{
				String tempContents = fusionArray[0].replaceFirst("\\[", "").replaceFirst("\\]", "");
				if (tempContents.contains(",")){//only has digits
					String[] numericalLocantsOfChildTemp =tempContents.split(",");
					for (int i = 0; i < numericalLocantsOfChildTemp.length; i++) {
						numericalLocantsOfChild.add(numericalLocantsOfChildTemp[i]);
					}
					lettersMissing=1;
				}
				else{//only has letters
					char[] tempLetterLocantsOfParentCharArray = tempContents.toCharArray();
					for (int i = 0; i < tempLetterLocantsOfParentCharArray.length; i++) {
						letterLocantsOfParent.add(String.valueOf(tempLetterLocantsOfParentCharArray[i]));
					}
					numbersMissing=1;
				}
			}
		}
		else{
			numbersMissing=1;
			lettersMissing=1;
		}

		if (numbersMissing==1){
			List<Atom> atomlist = rings.get(0).getAtomList();
			int foundCarbon=0;
			String locant="";
			atomlist.add(atomlist.get(0));//TODO use ringIterators
			for (Atom atom : atomlist) {
				if (atom.getElement().matches("C")){
					if (!atom.getFirstLocant().matches("\\d*")){
						continue;
					}
					if (foundCarbon ==1 ){//two in a row ->use this side
						numericalLocantsOfChild.add(locant);
						numericalLocantsOfChild.add(atom.getFirstLocant());
						break;
					}
					foundCarbon =1;
					locant =atom.getFirstLocant();
				}
				else{
					foundCarbon =0;
				}
			}
			atomlist.remove(atomlist.size()-1);
		}

		if (lettersMissing==1){
			List<Atom>  atomlist = rings.get(1).getAtomList();
			int foundCarbon=0;
			String locant ="";
			atomlist.add(0, atomlist.get(atomlist.size() -1));//TODO use ringIterators
			for (int i =atomlist.size() -1; i >=0; i--) {
				Atom atom =atomlist.get(i);
				if (atom.getElement().matches("C")){
					if (!atom.getFirstLocant().matches("\\d*")){
						continue;
					}
					locant =atom.getFirstLocant();
					if (foundCarbon ==1 ){//two in a row ->use this side
						letterLocantsOfParent.add(String.valueOf((char)(Integer.parseInt(locant) +96)));
						break;
					}
					foundCarbon =1;
				}
				else{
					foundCarbon =0;
				}
			}
			atomlist.remove(0);
		}


		HashMap <String, ArrayList<String>> tempHash = new HashMap<String, ArrayList<String>>();
		tempHash.put("numericalLocants",numericalLocantsOfChild);
		ringsInformation.put(rings.get(0), tempHash );
		tempHash = new HashMap<String, ArrayList<String>>();
		tempHash.put("letterLocants",letterLocantsOfParent);
		ringsInformation.put(rings.get(1), tempHash );

		Fragment fusedRing =fuseRings(state, rings, ringsInformation);//fuses the rings using the information contained with ringsInformation
		String fusedRingName=groups.get(0).getValue();
		if (fusedRingName.equals("benz") || fusedRingName.equals("benzo")){
			benzoSpecificAssignHeteroAtomsUsingLocants(state, groups.get(0), fusedRing);
		}
		for(int i=0;i<fusions.size();i++) {
			fusedRingName+=fusions.get(i).getValue();
		}
		fusedRingName+=groups.get(1).getValue();

		Element fusedRingEl =groups.get(groups.size()-1);//reuse this element to save having to remap suffixes...
		fusedRingEl.getAttribute("value").setValue(fusedRingName);
		fusedRingEl.getAttribute("valType").setValue("generatedFragment");
		fusedRingEl.getAttribute("type").setValue("ring");
		fusedRingEl.getAttribute("subType").setValue("fusedRing");
		fusedRingEl.removeChildren();
		fusedRingEl.appendChild(fusedRingName);

		state.xmlFragmentMap.put(fusedRingEl, fusedRing);

		for(int i=0;i<groups.size() -1;i++) {
			groups.get(i).detach();
		}
		for(int i=0;i<fusions.size();i++) {
			groups.get(i).detach();
		}
		return fusedRingEl;
	}

	/**
	 * Performs ring fusions
	 * @param state 
	 * @param ringsInformation 
	 * @param rings 
	 * @return fused Fragment
	 * @throws StructureBuildingException
	 */
	private Fragment fuseRings(BuildState state, ArrayList<Fragment> rings, HashMap<Fragment, HashMap<String, ArrayList<String>>> ringsInformation) throws StructureBuildingException {
		ArrayList<List<Atom>> childAtomsToBeConnectedToParent = new ArrayList<List<Atom>>();
		ArrayList<Integer> parentAtomIds = new ArrayList<Integer>();
		Fragment childRing =rings.get(0);
		Fragment parentRing =rings.get(1);
		ArrayList<String> locantsOfChild = ringsInformation.get(childRing).get("numericalLocants");
		ArrayList<String> letterLocantsOfParent = ringsInformation.get(parentRing).get("letterLocants");
		String firstLocant =locantsOfChild.get(0);
		String secondLocant =locantsOfChild.get(1);

		List<Atom> sortedAtomsInChild =sortFragmentAtomListByLocant(childRing);
		int firstLocantIndex=-1;
		int nonBridgeHeads=0;
		for (int i = 0; i < sortedAtomsInChild.size(); i++) {
			if (firstLocant.equals(sortedAtomsInChild.get(i).getFirstLocant())){
				firstLocantIndex=i;
			}
			if (sortedAtomsInChild.get(i).getFirstLocant().matches("\\d*")){
				nonBridgeHeads++;
			}
		}

		if (Integer.parseInt(secondLocant) - Integer.parseInt(firstLocant) == 1 || (Integer.parseInt(secondLocant) + nonBridgeHeads ) - Integer.parseInt(firstLocant) ==1){ //indexes are ascending
			int endingIndex = firstLocantIndex + sortedAtomsInChild.size();
			int counter=0;
			for (int i = firstLocantIndex; i < endingIndex; i++) {
				int index= i;
				if (index >= sortedAtomsInChild.size()){
					index -=sortedAtomsInChild.size();
				}
				String locant = sortedAtomsInChild.get(index).getFirstLocant();

				if (counter<locantsOfChild.size()){
					if (!locant.equals(locantsOfChild.get(counter))){
						locantsOfChild.add(counter++, locant);
					}
					else{
						counter++;
					}
				}
			}
		}
		else{
			int endingIndex = firstLocantIndex - sortedAtomsInChild.size();
			int counter=0;
			for (int i = firstLocantIndex; i > endingIndex; i--) {
				int index= i;
				if (index < 0){
					index +=sortedAtomsInChild.size();
				}
				String locant = sortedAtomsInChild.get(index).getFirstLocant();

				if (counter<locantsOfChild.size()){
					if (!locant.equals(locantsOfChild.get(counter))){
						locantsOfChild.add(counter++, locant);
					}
					else{
						counter++;
					}
				}
			}
		}

		for (int i = 0; i < locantsOfChild.size(); i++) {
			Atom atom =childRing.getAtomByLocant(locantsOfChild.get(i));
			List<Atom> neighbours = childRing.getAtomNeighbours(atom);

			//remove neighbours that are going to be deleted as they are included in numericalLocantsOfChild
			List<Atom> neighboursToRemove = new ArrayList<Atom>();
			for (Atom neighbour : neighbours) {
				String neighbourLocant=neighbour.getFirstLocant();
				for (int j = 0; j < locantsOfChild.size(); j++) {
					if (neighbourLocant.equals(locantsOfChild.get(j))){
						neighboursToRemove.add(neighbour);
					}
				}
			}

			for (Atom atomToRemove : neighboursToRemove) {//neighbour atoms that are to be subsequently removed are dropped
				neighbours.remove(atomToRemove);
			}

			childAtomsToBeConnectedToParent.add(i, neighbours);
		}
		for (int i = 0; i < locantsOfChild.size(); i++) {
			childRing.removeAtomByLocant(locantsOfChild.get(i), state.fragManager);
		}

		parentAtomIds.add(0, (int)letterLocantsOfParent.get(0).charAt(0) -97);
		for (int i = 0; i < letterLocantsOfParent.size(); i++) {
			int ringEdgeStartAtomNumber =(int)letterLocantsOfParent.get(i).charAt(0) -97;//convert from lower case character through ascii to 0-23
			if (ringEdgeStartAtomNumber +1!=parentRing.getAtomList().size()){
				parentAtomIds.add(parentAtomIds.size(), ringEdgeStartAtomNumber +1);
			}
			else{
				parentAtomIds.add(parentAtomIds.size(), 0);
			}
		}

		List<Atom> sortedAtomsInParent =sortFragmentAtomListByLocant(parentRing);
		parentRing.importFrag(childRing);
		for (int i = 0; i < parentAtomIds.size(); i++) {
			for (Atom atom : childAtomsToBeConnectedToParent.get(i)) {
				//System.out.println("Atom ID " + atom.getID() +" bonded to " +  parentAtomIds.get(i));
				parentRing.addBond(new Bond(atom, sortedAtomsInParent.get(parentAtomIds.get(i)), 1));
			}
		}
		state.fragManager.removeFragment(childRing);
		numberFusedRing(parentRing);//numbers the fused ring;
		Fragment fusedRing =state.fragManager.copyAndRelabel(parentRing);//makes sure the IDs are continuous
		state.fragManager.removeFragment(parentRing);
		return fusedRing;
	}


	/**
	 * Numbers the fused ring
	 * Currently only works for a very limited selection of rings
	 * @param uniFrag
	 * @throws StructureBuildingException
	 */
	public void numberFusedRing(Fragment fusedRing) throws StructureBuildingException {
		
		ArrayList<Ring> rings = (ArrayList<Ring>) getSetOfSmallestRings(fusedRing);
		
		ArrayList<ArrayList<Atom>> atomSequences = new ArrayList<ArrayList<Atom>>();
		
		// Special case when there are only 2 rings. This is expected to be faster than a more thorough analysis		 
		if (rings.size() ==2){
			List<Atom> atomList =fusedRing.getAtomList();
			ArrayList<Atom> bridgeheads =new ArrayList<Atom>();
			for (Atom atom : atomList) {
				if (fusedRing.getAtomNeighbours(atom).size()==3){
					bridgeheads.add(atom);
				}
			}
			for (Atom bridgeheadAtom : bridgeheads) {
				List<Atom>  neighbours =fusedRing.getAtomNeighbours(bridgeheadAtom);
				for (Atom  neighbour :  neighbours) {
					if (!bridgeheads.contains(neighbour)){
						//found starting atom
						ArrayList<Atom> atomsVisited =new ArrayList<Atom>();
						atomsVisited.add(bridgeheadAtom);
	
						Atom nextAtom =neighbour;
						do{
							atomsVisited.add(nextAtom);
							List<Atom> possibleNextInRings =fusedRing.getAtomNeighbours( nextAtom);
							nextAtom=null;
							for (Atom nextInRing:  possibleNextInRings) {
								if (atomsVisited.contains(nextInRing)){
									continue;//already visited
								}
								else{
									nextAtom=nextInRing;
								}
							}
						}
						while (nextAtom != null);
						atomsVisited.remove(bridgeheadAtom);
						atomsVisited.add(bridgeheadAtom);//remove the bridgehead and then re-add it so that it is at the end of the list
						atomSequences.add(atomsVisited);
					}
				}
			}
		}
		else {			
			setFusedRings(rings);
			
			if (checkRingAreInChain(rings, fusedRing)) 
			{
				ArrayList<Ring> tRings = findTerminalRings(rings);
				Ring tRing = tRings.get(0);
				
				ArrayList<Bond> fusedBonds = tRing.getFusedBonds();
				if (fusedBonds == null || fusedBonds.size()<=0) throw new StructureBuildingException("No fused bonds found");
				if (fusedBonds.size()>1) throw new StructureBuildingException("Terminal ring connected to more than 2 rings");
				
				// if there are more, we should go through atom most counterclockwise in the ring segh all the tRings
				
				enumerateRingAtoms(rings, tRing);
				ArrayList<Ring> orderedRings = new ArrayList<Ring>();
				int[] path = getDirectionsPath(rings, fusedBonds.get(0), tRing, orderedRings);
				
				atomSequences = applyRules(path, orderedRings);
			}
			else if(checkRingsAre6Membered(rings))
			{	
				atomSequences = number6MemberRings(rings);
				
				// find missing atoms
				if(atomSequences.size()<=0) throw new StructureBuildingException("No path found");
				ArrayList<Atom> takenAtoms = atomSequences.get(0);
				
				ArrayList<Atom> missingAtoms = new ArrayList<Atom>();
				for(Atom atom : fusedRing.getAtomList()) {
					if(!takenAtoms.contains(atom)) missingAtoms.add(atom);
				}
				// add  missing atoms to each path
				for (ArrayList<Atom> path : atomSequences) {
					for(Atom atom : fusedRing.getAtomList()) {
						if(!path.contains(atom)) path.add(atom);							
					}
				}			
			}
			else {
				int i=1;
				for (Atom atom : fusedRing.getAtomList()) {
					atom.replaceLocant("X" + Integer.toString(i));
					i++;
				}
				return;
			}			
		}
		// find the preferred numbering scheme then relabel with this scheme
		Collections.sort( atomSequences, new SortAtomSequences());
		fusedRing.setDefaultInID(atomSequences.get(0).get(0).getID());
		FragmentManager.relabelFusedRingSystem(fusedRing, atomSequences.get(0));
		fusedRing.reorderAtomCollection(atomSequences.get(0));
				
	}
	
	//*****************************************************************************************************	
	private class ConnectivityTable
	{
		public ArrayList<Ring> col1 = new ArrayList<Ring>();
		public ArrayList<Ring> col2 = new ArrayList<Ring>();
		public ArrayList<Integer> col3 = new ArrayList<Integer>();
		public ArrayList<Ring> usedRings = new ArrayList<Ring>();		
	}
	
	/**
	 * Returns possible enumerations of atoms in a 6-member ring system
	 * @param rings
	 * @return
	 * @throws StructureBuildingException
	 */
	private ArrayList<ArrayList<Atom>> number6MemberRings(ArrayList<Ring> rings) throws StructureBuildingException
	{
		ArrayList<Ring> tRings = findTerminalRings(rings);
		if (tRings == null || tRings.size()<0) throw new StructureBuildingException("Terminal rings not found");
		Ring tRing = tRings.get(0);
		Bond b1 = getNonFusedBond(tRing.getBondSet());
		if(b1 == null) throw new StructureBuildingException("Non-fused bond at termial ring not found");
		// order first bring
		
		ConnectivityTable ct = new ConnectivityTable();
		buildTable(tRing, null, 0, b1, b1.getFromAtom(), ct);
		
		ArrayList<Integer> dirs = findLongestChainDirection(ct);
		
		// add all the paths together and return
		ArrayList<ArrayList<Atom>> paths = new ArrayList<ArrayList<Atom>>();
		for (Integer dir : dirs) {
			ArrayList<ArrayList<Atom>> dirPaths = findPossiblePaths(dir, ct);
			for (ArrayList<Atom> path : dirPaths) {
				paths.add(path);
			}
		}
				
		return paths;				
	}
	
	/**
	 * Finds possible variants of enumerating atoms in a given direction
	 * @param newDir
	 * @param ct
	 * @return
	 * @throws StructureBuildingException
	 */
	private ArrayList<ArrayList<Atom>> findPossiblePaths(int newDir, ConnectivityTable ct) throws StructureBuildingException
	{
		//ArrayList<Integer> col3 = new  ArrayList<Integer>(this.col3);
		if ( ct.col1.size() != ct.col2.size() || ct.col2.size() != ct.col3.size() || ct.col1.size() <= 0) throw new StructureBuildingException("Sizes of arrays are not equal");
		
		int n = ct.col3.size();
		int[] col3 = new int[n];
		int maxx = 0;
		int minx = 0;
		int maxy = 0;
		int miny = 0;
				
		// turn the ring system
		int i=0;
		for(; i<n; i++)
		{			
			col3[i] = changeDirectionWithHistory(ct.col3.get(i), -newDir, ct.col1.get(i).size());
		}
		
		// Find max and min coordinates for ringMap
		// we put the first ring into usedRings to start with it in the connection tbl	
		int nRings = ct.usedRings.size();
		int[][] coordinates = new int[nRings][]; // correspondent to usedRings
		Ring[] takenRings = new Ring[nRings];
		int takenRingsCnt = 0;
		
		takenRings[takenRingsCnt++] = ct.col1.get(0);
		coordinates[0] = new int[]{0,0};
		
		// go through the rings in a system
		// find connected to them and assign coordinates according the directions
		// each time we go to the ring, whose coordinates were already identified.
		for(int tr=0; tr<nRings-1; tr++) 
		{
			Ring c1 = takenRings[tr];
			if (c1 == null) throw new StructureBuildingException();
			
			int ic1 = ct.col1.indexOf(c1);
			int xy[] = coordinates[tr]; // find the correspondent coordinates for the ring
					
			if (ic1 >= 0)
			{	
				for (int j=ic1; j<ct.col1.size(); j++)	{
					if (ct.col1.get(j) == c1)
					{						
						Ring c2 = ct.col2.get(j);
						if (arrayContains(takenRings,c2)) continue;
						
						int[] newxy = new int[2]; 
						newxy[0] = xy[0] + Math.round(2 * countDX(col3[j])); 
						newxy[1] = xy[1] + countDH(col3[j]);
						
						if(takenRingsCnt>takenRings.length) throw new StructureBuildingException("Wrong calculations");
						takenRings[takenRingsCnt] = c2;
						coordinates[takenRingsCnt] = newxy; 
						takenRingsCnt++;
						
						if (newxy[0] > maxx) maxx = newxy[0];
						else if (newxy[0] < minx) { minx = newxy[0]; }
						if (newxy[1] > maxy) maxy = newxy[1];
						else if (newxy[1] < miny) { miny = newxy[1];}
					}
				}
			}
		}
		// the height and the width of the map
		int h = maxy - miny + 1;
		int w = maxx - minx + 1;
		
		Ring[][] ringMap = new Ring[w][h];
		
		// Map rings using coordinates calculated in the previous step, and transform them according to found minx and miny
		
		int ix = -minx;
		int iy = -miny;		
		if (ix >= w || iy >= h) throw new StructureBuildingException("Coordinates are calculated wrongly");		
		ringMap[ix][iy] = ct.col1.get(0);
		
		int curx = 0;
		int cury = 0;
		for (int ti = 0; ti<takenRings.length; ti++)
		{
			int[] xy = coordinates[ti]; 
			curx = xy[0] - minx;
			cury = xy[1] - miny;
			if(curx<0 || curx>w || cury<0 || cury>h) throw new StructureBuildingException("Coordinates are calculated wrongly");
			ringMap[curx][cury] = takenRings[ti];
		}
		
		 
		ArrayList< int[]> chains = findChains(ringMap);
		
		// find candidates for different directions and different quadrants
		
		float[][] chainqs = new float[chains.size()][];
		// here we make array of quadrants for each chain
		for (int c=0; c<chains.size(); c++) {
			int[] chain = chains.get(c);
			int midChain = chain[0] + chain[1] - 1;
			
			float[] qs = countQuadrants(ringMap, chain[0], midChain, chain[2] );
			chainqs[c] = qs;
		}
			
		ArrayList<ArrayList<Atom>> paths = new ArrayList<ArrayList<Atom>> ();
		
		//  order for each right corner candidates for each chain
		ArrayList<String> chainCandidates = new ArrayList<String>(); 
		rulesBCD(chainqs, chainCandidates);
		int c = 0;
		for(String cand : chainCandidates)
		{
			ArrayList<ArrayList<Atom>> chainPaths = new ArrayList<ArrayList<Atom>> ();
			int[] chain = chains.get(c);
			int midChain = chain[0] + chain[1] - 1;
			
			for(int qi=0; qi<cand.length(); qi++) {
				int qr = Integer.parseInt(cand.charAt(qi)+"");
				Ring[][] qRingMap = transformRingWithQuadrant(ringMap, qr);
				boolean inverseAtoms = false;
				if (qr == 1 || qr == 3) inverseAtoms = true;
				ArrayList<Atom> quadrantPath = orderAtoms(qRingMap, midChain, inverseAtoms); 
				chainPaths.add(quadrantPath);
			}
			for (ArrayList<Atom> chainPath : chainPaths) {
				paths.add(chainPath);
			}
			c++;				
		}
		
		return paths;
	}
	
	/**
	 * Enumerates the atoms in a system, first finds the uppermost right ring, takes the next neighbour in the clockwise direction, and so one until the starting atom is reached
	 * @param ringMap
	 * @param midChain
	 * @param inverseAtoms
	 * @return
	 * @throws StructureBuildingException
	 */
	private ArrayList<Atom> orderAtoms(Ring[][] ringMap, int midChain, boolean inverseAtoms) throws StructureBuildingException
	{
		int w = ringMap.length;
		if (w<0 || ringMap[0].length < 0) throw new StructureBuildingException("Mapping results are wrong");
		int h = ringMap[0].length;
		
		ArrayList<Atom> atomPath = new ArrayList<Atom>();
		
		// find upper right ring
		Ring iRing = null;
		for (int i=w-1; i>=0; i--) {
			if (ringMap[i][h-1] != null) { iRing = ringMap[i][h-1]; break; }
		}		
		if (iRing == null) throw new StructureBuildingException("Upper right ring not found");
		
		Ring prevRing = findUpperLeftNeighbor(ringMap, iRing); 		
		Bond prevBond = findFusionBond(iRing, prevRing);
		Bond nextBond = null;
				
		boolean finished = false;
		int stNumber;
		int endNumber;
		int size;
		Ring nextRing = null;
		
		while (!finished) // or nof rings cannot be, cause one ring can be taken 2 times
		{						
			size = iRing.size();
											
			stNumber = iRing.getBondNumber(prevBond) ;
		
			ArrayList<Bond> cbonds = iRing.getCyclicBondSet();
			ArrayList<Bond> fbonds = iRing.getFusedBonds();
			if (!inverseAtoms)
			{
				for(int bi=0; bi<size; bi++)
				{
					int i = (stNumber + bi + 1) % size; // +1 cause we start from the bond next to stBond and end with it
					// if this bond is fused
					Bond bond = cbonds.get(i);
					if(fbonds.contains(bond)) {
						nextBond = bond; break;
					}										
				}
			}
			else 
			{
				for(int bi=0; bi<size; bi++)
				{
					int i = (stNumber - bi -1 + size) % size; // -1 cause we start from the bond next to stBond and end with it
					// if this bond is fused
					Bond bond = cbonds.get(i);
					if(fbonds.contains(bond)) {
						nextBond = bond; break;
					}										
				}
			}
			
			if (nextBond == null) throw new StructureBuildingException();
			// next ring			
			for (Ring ring : nextBond.getFusedRings()) {
				if(ring != iRing) { nextRing = ring; break; }
			}
			
			endNumber = iRing.getBondNumber(nextBond) ;
			
			// Add atoms in order, considering inverse or not inverse
			if (!inverseAtoms)
			{
				Atom atom = null;
				
				// if distance between prev bond and cur bond = 1 (it means that fused bonds are next to each other), but not fused use another scheme				
				// we dont add that atom, cause it was added already
				if ( (endNumber - stNumber + size) % size != 1)
				{
					stNumber = (stNumber + 1) % size;
					endNumber = (endNumber - 1 + size ) % size;
					if (stNumber > endNumber) endNumber += size;
					
					// start from the atom next to fusion								
					for (int j = stNumber; j <= endNumber; j++) // change 4-2
					{
						atom = iRing.getCyclicAtomSet().get(j % size);
						if (atomPath.contains(atom)) { finished = true; break; }
						atomPath.add(atom);
					}
				}
			}			
			else 
			{
				Atom atom = null;
				
				// if distance between prev bond and cur bond = 1 (it means that fused bonds are next to each other), use another scheme				
				if ( ( stNumber - endNumber + size) % size != 1)
				{
					stNumber = (stNumber - 2 + size ) % size;
					endNumber = endNumber % size;				
					if (stNumber < endNumber) stNumber += size;
										
					for ( int j = stNumber; j >= endNumber; j-- ) 
					{				
						atom = iRing.getCyclicAtomSet().get(j % size);
						if (atomPath.contains(atom)) { finished = true; break;}
						atomPath.add(atom);
					}
				}
			}
			prevBond = nextBond;
			prevRing = iRing;
			iRing = nextRing;			
		}
			
		return atomPath;
	}
	
	/**
	 * Finds the neighbour ring, which is the uppermost and on the left side from the given ring. used to find previous bond for the uppermost right ring, from which we start to enumerate 
	 * @param ringMap
	 * @param iRing
	 * @return
	 * @throws StructureBuildingException
	 */	
	private Ring findUpperLeftNeighbor (Ring[][] ringMap, Ring iRing) throws StructureBuildingException
	{
		Ring nRing = null; 
		int minx = Integer.MAX_VALUE;
		int maxy = 0;
				
		for (Ring ring : iRing.getNeighbors())
		{
			// upper left would be previous ring
			int xy[] = findRingPosition(ringMap, ring);
			if (xy==null) throw new StructureBuildingException("Ring is not found on the map");
			
			if (xy[1] > maxy  ||  xy[1] == maxy && xy[0] < minx ) {
				maxy = xy[1];
				minx = xy[0];
				nRing = ring;
			}
		}
		return nRing;
	}

	/**
	 * Finds the position(i,j) of the ring in the map
	 * @param ringMap
	 * @param ring
	 * @return
	 * @throws StructureBuildingException
	 */
	private int[] findRingPosition(Ring[][] ringMap, Ring ring) throws StructureBuildingException
	{
		int w = ringMap.length;
		if (w<0 || ringMap[0].length < 0) throw new StructureBuildingException("Mapping results are wrong");
		int h = ringMap[0].length;
		
		for(int i=0; i<w; i++) {
			for(int j=0; j<h; j++) {
				if (ringMap[i][j] == ring) {
					return new int[]{i,j};
				}					
			}
		}
		
		return null;
	}	
	
	/**
	 * Having upper right corner candidate transform the map to place the candidate to upper right corner 
	 * @param ringMap
	 * @param rq
	 * @return
	 * @throws StructureBuildingException
	 */
	private Ring[][] transformRingWithQuadrant(Ring[][] ringMap, int rq) throws StructureBuildingException 
	{
		int w = ringMap.length;
		if (w<0 || ringMap[0].length < 0) throw new StructureBuildingException("Mapping results are wrong");
		int h = ringMap[0].length;
		
		if (rq == 0) return ringMap.clone();
		
		Ring[][] resMap = new Ring[w][h];
		for (int i=0; i<w; i++) {
			for (int j=0; j<h; j++) {
				if(rq == 1) resMap[w-i-1] [j] = ringMap[i][j];
				else if(rq == 2) resMap[w-i-1] [h-j-1] = ringMap[i][j];
				else if(rq == 3) resMap[i] [h-j-1] = ringMap[i][j];
			}
		}
		
		return resMap;
	}
	
	/**
	 * Finds all the chains and their data:  0-chain length, 1-iChain, 2-jChain, for current direction
	 * @param ringMap
	 * @return
	 * @throws StructureBuildingException
	 */
	private ArrayList< int[]> findChains(Ring[][] ringMap) throws StructureBuildingException
	{
		int w = ringMap.length;
		if (w<0 || ringMap[0].length < 0) throw new StructureBuildingException("Mapping results are wrong");
		int h = ringMap[0].length;
		
		ArrayList< int[]> chains = new ArrayList< int[]>(); // Arraylist containing int arrays with all data for each chain: 0-chain length, 1-iChain, 2-jChain 
		
		int maxChain = 0;
		int chain = 0;
				
		// Find the longest chain
		for (int j=0; j<h; j++)	{
			for (int i=0; i<w; i++)	 {
				if(ringMap[i][j] != null) {
					chain = 1;					
					while( i + 2*chain < w && ringMap[i + 2*chain][j] != null ) chain++; // *2 because along the x axe the step is 2
					if(chain >= maxChain) {
						int[] aChain = new int[]{chain, i, j};
						chains.add(aChain);
						maxChain = chain;
					}
					i += 2*chain;
				}
			}
		}
		
		// remove those chains that were added before we found max
		for (int i=chains.size()-1; i>=0; i--) {
			int[] aChain = chains.get(i);
			if (aChain[0] < maxChain) chains.remove(i);
		}
		
		return chains;
	}
	
	/**
	 * Counts number of rings in each quadrant
	 * @param ringMap
	 * @param chain
	 * @param midChain
	 * @param jChain
	 * @return
	 * @throws StructureBuildingException
	 */
	private float[] countQuadrants(Ring[][] ringMap, int chain, int midChain, int jChain) throws StructureBuildingException
	{		
		float[] qs = new float[4];
		int w = ringMap.length;
		if (w<0 || ringMap[0].length < 0) throw new StructureBuildingException("Mapping results are wrong");
		int h = ringMap[0].length;
				
		//	int midChain = iChain + chain - 1; // actually should be *2/2, because we need the middle of the chain(/2) and each step is equal to 2(*2)
		
		// Count rings in each quadrants 
		for (int i=0; i<w; i++)	 {
			for (int j=0; j<h; j++)	{
				if (ringMap[i][j] == null) continue;
				
				if (i == midChain || j == jChain ) // if the ring is on the axe
				{
					if( i==midChain && j > jChain ) { qs[0]+=0.5; qs[1]+=0.5; }
					else if( i==midChain && j < jChain ) { qs[2]+=0.5; qs[3]+=0.5; }
					else if( i<midChain && j==jChain ) { qs[1]+=0.5; qs[2]+=0.5; }
					else if( i>midChain && j==jChain ) { qs[0]+=0.5; qs[3]+=0.5; }
					// if ( i==midChain && j==jChain ) we dont do anything
				}	
				else if(i>midChain && j>jChain) qs[0]++;
				else if(i<midChain && j>jChain) qs[1]++;
				else if(i<midChain && j<jChain) qs[2]++;
				else if(i>midChain && j<jChain) qs[3]++;
			}
		}
		
		return qs;
	}
	/**
	 * Checks if array contains an object
	 * @param array
	 * @param c2
	 * @return
	 */
	private boolean arrayContains(Object[] array, Object c2)
	{
		for (int i=0; i<array.length; i++) if (c2 == array[i])  return true;
		return false;
	}

	/**
	 * Finds the longest chain of rings in a line, using connectivity table 
	 * @param ct
	 * @return
	 * @throws StructureBuildingException
	 */
	private ArrayList<Integer> findLongestChainDirection(ConnectivityTable ct) throws StructureBuildingException
	{		 
		if (ct.col1.size() != ct.col2.size() || ct.col2.size() != ct.col3.size()) throw new StructureBuildingException("Sizes of arrays are not equal");;
	
		ArrayList<Integer> directions = new  ArrayList<Integer>();
		ArrayList<Integer> lengths = new  ArrayList<Integer>();
		
		// Ring c1;
		Ring c2;		
		int curChain;
		int curDir;
		int maxChain = 0;
		
		for (int i=0; i<ct.col1.size(); i++)
		{				
			c2 = ct.col2.get(i);
			curChain = 1;
			curDir = ct.col3.get(i);
			boolean chainBreak = false;
			
			while (!chainBreak)
			{
				int ic2 = ct.col1.indexOf(c2);
				boolean nextFound = false;
				
				if (ic2 >= 0)
				{	
					for (int j=ic2; j<ct.col1.size(); j++)	{
						if (ct.col1.get(j) == c2 && ct.col3.get(j) == curDir)
						{				
							curChain++;					 
							c2 = ct.col2.get(j);
							nextFound = true;
							break;
						}
					}
				}
				
				if(!nextFound)
				{
					if (curChain >= maxChain ) 
					{							
						maxChain = curChain;
						int oDir = getOppositeDirection(curDir);
						// if we didn't have this direction before, and opposite too, it is the same orientation 
						if(!directions.contains(curDir) && ! directions.contains(oDir)) {								
							directions.add(curDir);
							lengths.add(curChain);
						}
					}
					
					chainBreak = true;
				}
				
			}
			
		}
				
		// take  those with length equal to max
		for (int k = lengths.size()-1; k >= 0; k--) {
			if(lengths.get(k) < maxChain){
				lengths.remove(k);
				directions.remove(k);
			}				
		}
		return directions;						
	}
	
	/**
	 * Recursive function creating the connectivity table of the rings, for each connection includs both directions
	 * @param iRing
	 * @param parent
	 * @param prevDir
	 * @param prevBond
	 * @param atom
	 * @param ct
	 * @throws StructureBuildingException
	 */
	private void buildTable(Ring iRing, Ring parent, int prevDir, Bond prevBond, Atom atom, ConnectivityTable ct) throws StructureBuildingException
	{
		
		
		// order atoms and bonds in the ring		
		iRing.makeCyclicSets(prevBond, atom);
		ct.usedRings.add(iRing);
		
		for (Ring ring : iRing.getNeighbors())
		{
			// go back to ring we come from too, take the connection 2 times
			
			// the rings that are inside are necessary, cause we consider them when counting quadrants.
			// if (ring.size() - ring.getNOFusedBonds() <=0) continue;
			
			// find direction
			Bond curBond = findFusionBond(iRing, ring);
			// calculateRingDirection(iRing, prevBond, curBond, prevDir);
			
			int dir = 0;
			if (ring == parent) {
				dir =getOppositeDirection(prevDir);
			}
			else dir = calculateRingDirection(iRing, prevBond, curBond, prevDir);
				
			
			// place into connectivity table, like graph, rings and there connection
			ct.col1.add(iRing);
			ct.col2.add(ring);
			ct.col3.add(dir);
			
			if (!ct.usedRings.contains(ring))
			{				
				Atom a = getAtomFromBond(iRing, curBond);
				buildTable(ring, iRing, dir, curBond, a, ct);
			}			
		}		
	}
	
	/**
	 * Just returns any non fused bond
	 * @param bondSet
	 * @return
	 */
	private Bond getNonFusedBond(ArrayList<Bond> bondSet)
	{
		for (Bond bond : bondSet) {
			if(bond.getFusedRings() == null || bond.getFusedRings().size() < 1)
				return bond;
		}
		return null;
	}
	
	/**
	 * having the direction of the bond from ring1 to ring2, returns the opposite direction: from ring2 to ring1
	 * @param prevDir
	 * @return
	 */
	private int getOppositeDirection(int prevDir)
	{
		int dir;
		if (prevDir == 0) dir = 4;
		else if (Math.abs(prevDir) == 4) dir =0;
		else if (Math.abs(prevDir) == 1) dir = 3 * (-1) * (int) Math.signum(prevDir);
		else dir = 1 * (-1) * (int) Math.signum(prevDir);
		return dir;
	}
	
	/**
	 * Finds the atom connected to the bond, takes into account the order of the bonds and atoms in the ring
	 * @param ring
	 * @param curBond
	 * @return
	 * @throws StructureBuildingException
	 */
	private Atom getAtomFromBond(Ring ring, Bond curBond) throws StructureBuildingException
	{
		if (ring.getCyclicBondSet() == null) throw new StructureBuildingException("Atoms in the ring are not ordered");
		int i=0;
		for (Bond bond : ring.getCyclicBondSet())	{
			if (bond == curBond) break;
			i++;
		}
		int ai = ( i - 1 + ring.size() ) % ring.size();
		Atom atom = ring.getCyclicAtomSet().get(ai);
		return atom;
	}		
	
	/**
	 * Finds the fusion bond between 2 rings
	 * @param r1
	 * @param r2
	 * @return
	 */
	private Bond findFusionBond (Ring r1, Ring r2)
	{
		ArrayList<Bond> b2 = r2.getBondSet();
		for(Bond bond : r1.getBondSet()) 
			if (b2.contains(bond)) return bond;
		
		return null;
	}
	
	/**
	 * Calculates the direction of the next ring according to the distance between fusing bonds and the previous direction
	 * @param ring
	 * @param prevBond
	 * @param curBond
	 * @param history
	 * @return
	 * @throws StructureBuildingException
	 */
	private int calculateRingDirection(Ring ring, Bond prevBond, Bond curBond, int history) throws StructureBuildingException
	{
		// take the ring fused to one from the previous loop step				
		if ( ring.getCyclicBondSet() == null ) throw new StructureBuildingException();
		int size = ring.size();
									
		int i1 = -1;
		int i2 = -1;
		int cnt = 0;
		for(Bond bond :ring.getCyclicBondSet()) 
		{
			if (bond == prevBond) i1=cnt;
			if (bond == curBond) i2=cnt;
			if (i1>=0 && i2>=0) break;
			cnt++;
		}
		
		int dist = (size + i2 - i1) % size;
		
		if (dist == 0) throw new StructureBuildingException("Distance between bonds is equal to 0");
		
		int dir = getDirectionFromDist(dist, size, history);
				
		return dir;
	}
	/**
	 * Check if all the rings in a system are 6 membered
	 * @param rings
	 * @return
	 */
	private  boolean checkRingsAre6Membered(ArrayList<Ring> rings)
	{
		for (Ring ring : rings) {
			if (ring.size() != 6) return false;
		}
		return true;
	}
	
	//*****************************************************************************************************
	
	
	
	/**
	 * Finds the longest chains and call the function assigning the atoms order.
	 * @param path
	 * @param pRings
	 * @return
	 * @throws StructureBuildingException
	 */	
	private ArrayList<ArrayList<Atom>> applyRules(int[] path, ArrayList<Ring> pRings) throws StructureBuildingException 
	{		
		int i;
		int si=0;
		int maxChain=0;
		
		int curDir = 0;
		ArrayList<Integer> maxDir = new ArrayList<Integer>();		
		
		ArrayList<ArrayList<Atom>> allAtomOrders = new ArrayList<ArrayList<Atom>>();
		
		//int[] path = {0,0,1,1};
		
		// Find the length of the longest chain
		for (i=0; i<path.length; i++)
		{	
			if (path[i] == curDir)  {				
				si++;
			}
			else  // if we switch the direction
			{				
				if (maxChain<si) {
					maxChain = si;					
				}
				
				si = 1; // start to count current symbol
				curDir = path[i];
			}
		}		
		//if sequence ends at the end of array
		if (maxChain<si) {
			maxChain = si;			
		}
		// TODO change here, with delete
		curDir = 0;
		si = 0;
		for (i=0; i<path.length; i++)
		{
			if (path[i] == curDir)	{				
				si++;
			}
			else  // if we switch the direction
			{				
				if (maxChain == si && !maxDir.contains(curDir)) {
					maxDir.add(curDir);
				}
				
				si = 1; // start to count current symbol
				curDir = path[i];
			}
		}
		if (maxChain == si && !maxDir.contains(curDir)) {
			maxDir.add(curDir);
		}
		
		if (maxDir.size()<=0) throw new StructureBuildingException("Chains are not recognized in the molecule"); 
			
		for (int dir : maxDir) {
			ArrayList<ArrayList<Atom>> orders = getOrderInEachDirection(path, pRings, dir, maxChain);
			for (ArrayList<Atom> order : orders) {
				allAtomOrders.add (order);
			}
		}
		
		return allAtomOrders;
	}
	
	
	
	
	/**
	 * get the direction of the main chain, according to it recalculates the path, check the rule, calculating number of rings in each quadrant 
	 * @param path
	 * @param pRings
	 * @param maxDir
	 * @param maxChain
	 * @return
	 * @throws StructureBuildingException
	 */	
	private ArrayList<ArrayList<Atom>> getOrderInEachDirection(int[] path, ArrayList<Ring> pRings, int maxDir, int maxChain) throws StructureBuildingException
	{
		//for further analyses we change the orientation.		
		int i;		
		ArrayList<String> chainVariants = new ArrayList<String>();
		boolean sequence = false;
		path = path.clone(); // not to change real object
		
		//  if 4 than change the path and the order of rings
		if (Math.abs(maxDir) == 4)
		{
			int[] copyPath = path.clone();
			ArrayList<Ring> inversedRings = new ArrayList<Ring>();
			int l = path.length;
			
			if (pRings.size() < path.length) throw new StructureBuildingException("The path does not correspond to array of rings");
			for (i=0; i<l; i++) {
				if ( copyPath[i] == 0 ) path[l-i-1] = 4;
				else if ( Math.abs( copyPath[i] ) == 1 ) path[l-i-1] = 3 * (int) Math.signum(copyPath[i]) * (-1);
				else if ( Math.abs( copyPath[i] ) == 3 ) path[l-i-1] = 1 * (int) Math.signum(copyPath[i]) * (-1);
				else if ( Math.abs( copyPath[i] ) == 4 ) path[l-i-1] = 0;	
				
				inversedRings.add(pRings.get(l-i));				
			} 
			inversedRings.add(pRings.get(0));
			pRings = inversedRings; // change the reference, but not changing initial values of pRings 
		}
		// changed with function
		else if (maxDir != 0 ){			
			for (i=0; i<path.length; i++) {
				path[i] = changeDirectionWithHistory(path[i], -maxDir, pRings.get(i).size());
			}
		}
		
		// Rule A
		ArrayList<Integer> chains = new ArrayList<Integer>();		
		// find the longest chains
		int si=0;
		for (i=0; i<path.length; i++)
		{			
			if (path[i] == 0 ) //|| path[i] == 4)
			{
				if (!sequence) sequence  = true;
				si++;
			}
			else if (sequence) // if we finished the chain
			{
				sequence = false;
				if (maxChain==si) chains.add(i-maxChain); // the begining of the chain
				si=0;
			}			
		}
		// if chain  finishes with the end of path
		if (sequence && maxChain==si) chains.add(i-maxChain);
		
		
		float[][] qs = new float [chains.size()][];
		int c=0;
		
		for (Integer chain : chains) 
		{			
			qs[c] = countQuadrants(chain, maxChain, path);						
			c++;
		}
		
		rulesBCD(qs, chainVariants);
		
		int ichain = 0;
		
		ArrayList<Ring> inversedRings = new ArrayList<Ring>(); 
		for (int k=pRings.size()-1; k>=0; k--) { // make once inversed, dont repeat
			inversedRings.add(pRings.get(k));
		} 
				
		ArrayList<ArrayList<Atom>> atomOrders = new ArrayList<ArrayList<Atom>>();
				
		for (String chain :chainVariants)
		{			
			for (int j=0; j<chain.length(); j++)
			{
				int q = Integer.parseInt(chain.charAt(j)+ "");
				
				boolean inverseAtoms = false; 
				if (q == 1 || q ==3) inverseAtoms = true;
				
				int stChain = chains.get(ichain);
				
				ArrayList<Ring> ringsToPass;
				if (q == 1 || q == 2) {
					stChain = path.length - stChain - maxChain;
					ringsToPass = inversedRings;
				}
				else ringsToPass = new ArrayList<Ring>(pRings); 
				
				ArrayList<Atom> oAtoms = createAtomOrder(ringsToPass, getTransformedPath(path, q), stChain, maxChain, q, inverseAtoms);
				atomOrders.add(oAtoms);
				
			}
			ichain++;
		}
		
		return  atomOrders;
	}
	/**
	 * Applying rules B, C and D for the ring system. The function is used for both types of ring systems. 
	 * @param qs - array with number of ring in each quadrant for each chain.
	 * @param chainVariants 
	 * @throws StructureBuildingException
	 */
	private void rulesBCD(float[][] qs, ArrayList<String> chainVariants) throws StructureBuildingException
	{
		// Analyse quadrants
		
		// Rule B: Maximum number of rings in upper right quadrant. Upper right corner candidates		
		int variantNumber = 0;
		float qmax = 0;
		int c=0;
		int nchains = qs.length;
		
		for (c=0; c<nchains; c++) 
		{				
			for (int j=0; j<4; j++)	{
				if(qs[c][j]>qmax) qmax = qs[c][j];				
			}			
		}
		
		for (c=0; c<nchains; c++)
		{
			String taken = "";
			for (int j=0; j<4; j++){
				if (qs[c][j]==qmax) { taken += j; variantNumber++; }
			}
			chainVariants.add(taken);
		}
		
		
		
		// Rule C: Minimum number of rings in lower left quadrant
		if (variantNumber > 1)
		{		
			c=0;
			variantNumber = 0;
			float qmin = Integer.MAX_VALUE;
			
			for (String chain : chainVariants) {
				for (int j=0; j<chain.length(); j++)
				{
					int q = Integer.parseInt(chain.charAt(j)+ "");
					int qdiagonal = (q + 2) % 4;
					if (qs[c][qdiagonal]<qmin) qmin = qs[c][qdiagonal];
				}
				c++;
			}
			c=0;
			for (String chain : chainVariants) {
				String taken = "";
				for (int j=0; j<chain.length(); j++)
				{
					int q = Integer.parseInt(chain.charAt(j)+ "");
					int qdiagonal = (q + 2) % 4;
					if (qs[c][qdiagonal]==qmin) { taken += q; variantNumber++;}				
				}
				chainVariants.set(c, taken);
				c++;
			}
		}
		else if (variantNumber <= 0)
			throw new StructureBuildingException("Atom enumeration path not found");
		
		
		// Rule D: Maximum number of rings above the horizontal row
		if (variantNumber > 1)
		{
			c=0;		
			float rmax = 0;
			variantNumber = 0;
			for (String chain : chainVariants) {
				for (int j=0; j<chain.length(); j++)
				{
					int q = Integer.parseInt(chain.charAt(j)+ "");
					int qrow;
					if (q % 2 == 0) qrow = q + 1;
					else qrow = q - 1;
		
					if (qs[c][qrow] + qs[c][q] > rmax) rmax = qs[c][qrow] + qs[c][q];
				}
				c++;
			}
			c=0;
			for (String chain : chainVariants) {
				String taken = "";
				for (int j=0; j<chain.length(); j++)
				{
					int q = Integer.parseInt(chain.charAt(j)+ "");
					int qrow;
					if (q % 2 == 0) qrow = q + 1;
					else qrow = q - 1;
						
					if (qs[c][qrow] + qs[c][q] == rmax) { taken += q; variantNumber++; }				
				}
				chainVariants.set(c, taken);
				c++;
			}
		}
		
		if (variantNumber <= 0)
			throw new StructureBuildingException("Atom enumeration path not found");
		
			
		
	}
	/**
	 * Adds atoms to the array in the order according to the rules,. Finds the uppermost right ring, starting from it adds atoms to the result arraylist 
	 * @param rings
	 * @param path
	 * @param chainStart
	 * @param chainLen
	 * @param quadrant
	 * @param inverseAtoms
	 * @return
	 * @throws StructureBuildingException
	 */
	private ArrayList<Atom> createAtomOrder(ArrayList<Ring> rings, int[] path, int chainStart, int chainLen, int quadrant, boolean inverseAtoms) throws StructureBuildingException
	{
		// atom order is changed when we transform the right corner from the 1st and 3d quadrant (start from 0)
		// rings order is inversed when we transform from 1st and 2d quadrant (start from 0)
		ArrayList<Atom> atomPath = new ArrayList<Atom>();
		
		// find the upper right ring				
		int height=0;
		int maxheight = 0;
		float xdist = (float) chainLen / 2;
		float maxDist = xdist;
		boolean foundAfterChain = true;
		int[] ringHeights = new int[path.length+1];
				
		int  i = chainStart + chainLen;
		int upperRightPath = i-1;// if nothing after chain
		
		for ( ; i<path.length; i++)
		{
			height += countDH(path[i]);			
			xdist += countDX(path[i]);
			ringHeights[i+1] = height;
			
			// take the ring if it is the highest and in the right quadrant, and then take the most right
			if ( (height > maxheight && xdist >= 0) || (height == maxheight &&  xdist > maxDist) )
			{
				maxheight = height;				
				maxDist = xdist;
				upperRightPath = i;
			}			
		}
		
		// if we assume that the path can come from the left side to the right, then we should check the beginning of the path 
		height=0;		
		xdist = -(float) chainLen / 2;
				
		i = chainStart - 1;
				
		for ( ; i>=0; i--)
		{
			height -= countDH(path[i]);
			xdist -= countDX(path[i]);
			ringHeights[i+1] = height;
									
			// take the ring if it is the highest and in the right quadrant, and then take the most right
			if ( (height > maxheight && xdist >= 0) || (height == maxheight &&  xdist > maxDist) )
			{
				maxheight = height;				
				maxDist = xdist;
				upperRightPath = i;
				foundAfterChain = false;
			}			
		}
		//  if we found the ring by backtracing we dont need to decrease
		if (foundAfterChain) upperRightPath++; // because we have 1 less elements in the path array		
		if (upperRightPath<0 || upperRightPath>rings.size()) throw new StructureBuildingException();
		
		
		ArrayList<Bond> fusedBonds;
		Bond prevFusedBond = null;
		Bond curFusedBond = null;
		Ring ring = rings.get(upperRightPath);
		boolean finished = false;
		
		// if only one fused bond - the terminal ring is the first ring
		// otherwise we need to find the "previous" bond		
		fusedBonds = ring.getFusedBonds();
		if (fusedBonds.size()==1) prevFusedBond = null;
		else 
		{
			// the ring should be between 2 rings, otherwise it is terminal
			if (upperRightPath+1 >= ringHeights.length || upperRightPath-1 < 0) throw new StructureBuildingException();
			
			Ring prevRing = null;
			
			if (ringHeights[upperRightPath-1] > ringHeights[upperRightPath+1]) prevRing = rings.get(upperRightPath-1);			
			else prevRing = rings.get(upperRightPath+1);
			
			for (Bond bond : fusedBonds) {
				if (bond.getFusedRings().contains(prevRing)) { prevFusedBond = bond; break;}			
			}
			if (prevFusedBond == null) throw new StructureBuildingException();
		}
		
		
		while (!finished) // we go back from the right corner   // we can also ask if there this ring is equal to the first one
		{			
			fusedBonds = ring.getFusedBonds();
						
			// if there is only one fused bond: EITHER the current would be that one and prev=null (1st iteration) OR current would be equal to prev.
			for (Bond bond : fusedBonds) {
				if (bond != prevFusedBond) {
					curFusedBond = bond; 
				}
			}
			
			if (prevFusedBond==null) prevFusedBond = curFusedBond;
			
			int size = ring.size();
			
			int stNumber = ring.getBondNumber(prevFusedBond) ;
			int endNumber = ring.getBondNumber(curFusedBond) ;
			
			if (!inverseAtoms)
			{
				Atom atom = null;
				
				// if distance between prev bond and cur bond = 1 (it means that fused bonds are next to each other), use another scheme				
				// we dont add that atom, cause it was added already
				if ( (endNumber - stNumber + size) % size != 1)
				{
					stNumber = (stNumber + 1) % size;
					endNumber = (endNumber - 1 + size ) % size;
					if (stNumber > endNumber) endNumber += size;
					
					// start from the atom next to fusion								
					for (int j = stNumber; j <= endNumber; j++) // change 4-2
					{
						atom = ring.getCyclicAtomSet().get(j % size);
						if (atomPath.contains(atom)) { finished = true;  break; } 
						atomPath.add(atom); 
					}
				}
			}			
			else 
			{
				Atom atom = null;
				
				// if distance between prev bond and cur bond = 1 (it means that fused bonds are next to each other), use another scheme				
				if ( ( stNumber - endNumber + size) % size != 1)
				{
					stNumber = (stNumber - 2 + size ) % size;
					endNumber = endNumber % size;				
					if (stNumber < endNumber) stNumber += size;
										
					for ( int j = stNumber; j >= endNumber; j-- ) 
					{				
						atom = ring.getCyclicAtomSet().get(j % size);
						if (atomPath.contains(atom)) { finished = true; break; } 
						atomPath.add(atom);
						 
					}
				}
			}
			
			if (finished) break;
			
			ArrayList<Ring> fusedRings = curFusedBond.getFusedRings();
			for (Ring fRing : fusedRings) {
				if (ring != fRing) { ring = fRing; break;}
			}
			
			prevFusedBond = curFusedBond;
		}
	
		return atomPath;
	}	
	/**
	 * Transforms the given path according to the quadrant proposed as the right corner
	 * @param path
	 * @param urCorner
	 * @return
	 */	
	private int[] getTransformedPath(int[] path, int urCorner)
	{
		int l = path.length;
	
		int[] rPath = new int[l];
		
		for(int i=0; i<l; i++)
		{
			if(urCorner == 1){
				rPath[i] = -path[l-i-1]; // changes ring order
			}
			else if(urCorner == 2){
				rPath[i] = path[l-i-1]; // changes ring order
			}
			else if(urCorner == 3){
				rPath[i] = -path[i];
			}
			else{
				rPath[i] = path[i];
			}
		}
		return rPath;
	}
	
	/**
	 * Counts number of rings in each quadrant
	 * @param start
	 * @param len
	 * @param path
	 * @return
	 */
	private float[] countQuadrants(int start, int len, int[] path)
	{
		int i;
		float[] quadrants = new float[4]; // 0-ur, 1-ul, 2-ll, 3-lr, counter-clockwise
		int height = 0;
		float xdist = -(float) len/ 2;
		
		// count left side
		for (i=start-1; i>=0; i--)
		{
			height -= countDH(path[i]);			
			xdist -= countDX(path[i]);
			incrementQuadrant(height, xdist, quadrants);			
		}
		
		height = 0;
		xdist = (float) len/ 2;
		// count right side
		for (i=start+len; i<path.length; i++)
		{			
			height += countDH(path[i]);
			xdist += countDX(path[i]);
			
			incrementQuadrant(height, xdist, quadrants);
		}
		
		return quadrants;
	}
	/**
	 * Used to add the number of the rings to a quadrant, according to the position of the ring
	 * @param height
	 * @param xdist
	 * @param quadrants
	 */
	private void incrementQuadrant(int height, float xdist, float[] quadrants)
	{
		if (height > 0) 
		{
			if (xdist>0) quadrants[0]+=1; // right side
			else if (xdist<0) quadrants[1]+=1; // left side
			else {  // middle
				 quadrants[0]+=0.5;
				 quadrants[1]+=0.5;
			}
		}
		else if (height < 0)  
		{
			if (xdist>0) quadrants[3]+=1; // right side
			else if (xdist<0) quadrants[2]+=1; // left side
			else {  // middle
				 quadrants[3]+=0.5;
				 quadrants[2]+=0.5;
			}				
		}
		else // height=0
		{
			if (xdist>0) 
			{					
				quadrants[0]+=0.5f; 
				quadrants[3]+=0.5f;
			}
			else if (xdist<0) 
			{
				quadrants[1]+=0.5f; 
				quadrants[2]+=0.5f;					
			}
			else { 
				 // should actually add 0.25 to each, but this doesnt make any change.
			}
		}
	}
	/**
	 * Counts delta x distance between previous and next rings
	 * @param val
	 * @return
	 */
	private float countDX (int val)
	{
		float dx = 0;
		
		if (Math.abs(val) == 1) dx += 0.5f;
		else if (Math.abs(val) == 3) dx -= 0.5f;
		else if (Math.abs(val) == 0) dx += 1f;
		else if (Math.abs(val) == 4) dx -= 1f;
		
		return dx;
	}
	/**
	 * Counts delta height between previous and next rings
	 * @param val
	 * @return
	 */
	
	private int countDH(int val)
	{
		int dh = 0;
		if (Math.abs(val) != 4)
		{
			if (val>0) dh = 1;
			if (val<0) dh = -1;
		}
		return dh;
	}
		
	/**
	 * Finds the rings with the min fused bonds
	 * @param rings
	 * @return
	 */
	private ArrayList<Ring> findTerminalRings(ArrayList<Ring> rings)
	{
		ArrayList<Ring> tRings = new ArrayList<Ring>(); 
		
		int minFusedBonds = Integer.MAX_VALUE;
		for  (Ring ring : rings) 
		{
			if (ring.getNOFusedBonds() < minFusedBonds) minFusedBonds = ring.getNOFusedBonds();
		}
		
		for  (Ring ring : rings) 
		{
			if (ring.getNOFusedBonds() == minFusedBonds) tRings.add(ring);
		}
		return tRings;
	}
	
	/**
	 * Fills the value fusedRings for bonds, calcuclates the number of fused bonds in a ring
	 * @param rings
	 */
	private void setFusedRings(ArrayList<Ring> rings)
	{		
		for (Ring curRing : rings) {
			for(Bond bond : curRing.getBondSet()) { 			// go through all the bonds for the current ring
				if (bond.getFusedRings().size()>=2) continue; 	// it means this bond we already analysed and skip it  
				
				for (Ring ring : rings) {  						// check if this bond belongs to any other ring 
					if (curRing != ring) {
						if (ring.getBondSet().contains(bond)) {
							bond.addFusedRing(ring);			// if so, then add the rings into fusedRing array in the bond
							bond.addFusedRing(curRing);			// and decrease number of free bonds for both rings
							
							ring.incrNOFFusedBonds();
							curRing.incrNOFFusedBonds();
							
							ring.addNeighbor(curRing);
							curRing.addNeighbor(ring);
							continue;
						}					
					}
				}
			}
		}
	}
	
	/**
	 * Checks if the given fragment is chain type ring system
	 * @param rings
	 * @param frag
	 * @return
	 */	
	private boolean checkRingAreInChain(ArrayList<Ring> rings, Fragment frag)
	{
		for (Ring ring : rings) {
			if (ring.getNOFusedBonds() > 2) return false;
			if (ring.size()>9) return false;
		}
		
		for (Atom atom : frag.getAtomList()){
			List<Bond> bonds = atom.getBonds();			
			if (bonds.size()>2){
				int nFused = 0;
				for (Bond bond : bonds) {
					if (bond.getFusedRings().size()>1) nFused++;
				}
				if (nFused>2) return false;
			}
		}		
		return true;
	}	

	/**
	 * Orders atoms in each ring, enumeration depends on enumeration of previous ring
	 * @param rings
	 * @param tRing
	 * @throws StructureBuildingException
	 */	
	private void enumerateRingAtoms(ArrayList<Ring> rings, Ring tRing) throws StructureBuildingException
	{
		if (rings == null || rings.size()<=0) throw new StructureBuildingException();
			
		Ring iRing = tRing;		
		Bond stBond = tRing.getBondSet().get(0);
		Atom stAtom = stBond.getToAtom();
		
		for (int i=0; i<rings.size(); i++) 
		{
			iRing.makeCyclicSets(stBond, stAtom);
			
			if (i==rings.size()-1) break;
				
			// find the bond between current ring and next one
			ArrayList<Bond> fusedBonds = iRing.getFusedBonds();
			if (fusedBonds == null || fusedBonds.size()<=0) throw new StructureBuildingException();			
			for (Bond fBond : fusedBonds) {
				if (stBond != fBond) {stBond = fBond; break;} // we take the bond different from the current
			}
						
			int cnt = 0;
			for (Bond bond : iRing.getCyclicBondSet()) {				
				if (bond == stBond)	{
					cnt--;
					if (cnt<0) cnt = iRing.size()-1;
					stAtom = iRing.getCyclicAtomSet().get(cnt); // so that the enumeration go the same direction we give the previous atom
					break;
				}
				cnt++;
			}

			// take next ring fused to the current			
			ArrayList<Ring> fusedRings = stBond.getFusedRings();
			if (fusedRings == null || fusedRings.size()<2) throw new StructureBuildingException(); 
			if (iRing != fusedRings.get(0)) iRing = fusedRings.get(0);
			else iRing = fusedRings.get(1);		
		}			
		
	}
	
	/**
	 * Creates an array describing mutual position of the rings
	 * @param rings
	 * @param startBond
	 * @param tRing
	 * @return
	 * @throws StructureBuildingException
	 */ 
	private int[] getDirectionsPath(ArrayList<Ring> rings, Bond startBond, Ring tRing, ArrayList<Ring> orderedRings) throws StructureBuildingException 
	{
		//String fPath = tRing.size() + "R"; // because from the first ring we go right 
		int[] path = new int[rings.size()-1];		
		path[0]=0;
		orderedRings.add(tRing);
		
		Ring iRing = tRing;
		Bond stBond = startBond;
		Bond nextBond=null;
		int history=0; // we store here the previous  direction
				
		for (int i=0; i<rings.size()-1; i++) 
		{
			// take the ring fused to one from the previous loop step			
			ArrayList<Ring> fusedRings = stBond.getFusedRings();
			if (fusedRings == null || fusedRings.size()<2) throw new StructureBuildingException(); 
			if (iRing != fusedRings.get(0)) iRing = fusedRings.get(0);
			else iRing = fusedRings.get(1);
			
			int size = iRing.size();
			orderedRings.add(iRing);
			
			// find the next fused bond between current ring and the next			
			ArrayList<Bond> fusedBonds = iRing.getFusedBonds();
			if (fusedBonds == null || fusedBonds.size()<=0) throw new StructureBuildingException();
			if (fusedBonds.size() == 1) break;		// we came to the last ring in the chain
			for (Bond fBond : fusedBonds) {
				if (stBond != fBond) {nextBond = fBond; break; } // we take the bond different from the current
			}
						
			int i1 = -1;
			int i2 = -1;
			int cnt = 0;
			for(Bond bond :iRing.getCyclicBondSet())
			{
				if (bond == stBond) i1=cnt;
				if (bond == nextBond) i2=cnt;
				if (i1>=0 && i2>=0) break;
				cnt++;
			}
			
			int dist = (size + i2 - i1) % size;
			
			if (dist == 0) throw new StructureBuildingException("Distance between fused bonds is equal to 0");
			
			int dir = getDirectionFromDist(dist, size, history);
			
			history = dir;
		
			path[i+1] = dir;
							
			stBond = nextBond;
		}
		
		return path;
	}
	
	// take history! or make just 2 directions
	private int getDirectionFromDist(int dist, int size, int history) throws StructureBuildingException
	{
		// positive val of n - Up
		// negative value - Down
		
		int dir=0;
		
		if (size >= 10) throw new StructureBuildingException("rings with more than 10 members are not recognized"); 
		
		if (size == 3) // 3 member ring
		{
			if (dist == 1) dir = 1;
			else if (dist == 2) dir = -1;
			else throw new StructureBuildingException();
		}
		else if (size == 4) // 4 member ring
		{
			if (dist == 2) dir = 0;
			else if (dist < 2) dir = 2;
			else if (dist > 2) dir = -2;
		}
		
		else if (size % 2 == 0) // even
		{
			if (dist == 1) dir = 3;			
			else if (dist == size-1) dir = -3;
			
			else
			{
				dir = size/2 - dist;			
				// 8 and more neighbours
				if (Math.abs(dir) > 2 && size >= 8) dir = 2 * (int) Math.signum(dir);
			}
		}
		else // odd
		{
			if (dist == size/2 || dist == size/2 + 1)  dir = 0;
			else if (size == 5) dir =2;
			
			else if (dist == size-1) dir = -3;
			else if (dist == 1) dir = 3;
			
			else if (size>=9 && dist == size/2-1) dir = 2; // direction number 2 appears only when 
			else if (size>=9 && dist == size/2+2) dir = 2;
				
			else if(dist < size/2) dir = 2;	
			else if(dist > size/2+1) dir = -2;
		}
		
		dir =changeDirectionWithHistory(dir, history, size);
		return dir;

	}	
	
	private int changeDirectionWithHistory(int dir, int history, int size) throws StructureBuildingException
	{
		int relDir = dir;
		
		if (Math.abs(history) == 4) 
		{
			if (dir == 0) dir = 4;
			else
				dir += 4 * (-1) * Math.signum(dir); // if dir<0 we add 4, if dir>0 we add -4
		}
		else
			dir += history;
		
		if (Math.abs(dir)>4) // Added
		{
			dir = Math.round( (8 - Math.abs(dir)) * Math.signum(dir) * (-1) );			
		}
			
		// 6 member ring does not have direction 2
		if (size == 6 && Math.abs(dir) == 2) 
		{
			//if (history == 1 || history == -3) dir++; 
			//else if (history == 3 || history == -1) dir--;
			// changed
			// if (one of them equal to 1 and another is equal to 3, we decrease absolute value and conserve the sign)			
			if (Math.abs(relDir)==1 && Math.abs(history)==3  ||  Math.abs(relDir)==3 && Math.abs(history)==1) {dir = 1 * (int) Math.signum(dir);}
			// if both are equal to 1
			else if(Math.abs(relDir)==1 && Math.abs(history)==1 ) {dir = 3 * (int) Math.signum(dir);}
			// if both are equal to 3
			else if(Math.abs(relDir)==3 && Math.abs(history)==3 ) {dir = 3 * (int) Math.signum(dir);}
			// else it is correctly 2 // else throw new StructureBuildingException();
		}
		
		if (dir == -4) dir = 4;
		
		return dir;
	}
	/** get set of smallest rings.
	 * crude
	 * @param update replace current list
	 * @return list of rings
	 */
	public List<Ring> getSetOfSmallestRings(Fragment frag) throws StructureBuildingException {
		ArrayList<Atom> atomSet = (ArrayList<Atom>) frag.getAtomList();
				
		List<Ring> ringList = getRings(atomSet);
		
		if (ringList.size() > 1) {
			boolean change = true;
			while (change) {
				for (int i = 0; i < ringList.size(); i++) {
					Ring ring = ringList.get(i);
					change = reduceRingSizes(ring, ringList);
				}
			}
		}
		else throw new StructureBuildingException("Ring perception system did not recognize fused rings");
				
		return ringList;
	}
	
	/** get list of rings.
	 * not necessarily SSSR
	 * @return list of rings
	 * @throws StructureBuildingException 
	 */
	public List<Ring> getRings(ArrayList<Atom> atomSet ) throws StructureBuildingException {
		List<Ring> ringList = new ArrayList<Ring>();
		Set<Atom> usedAtoms = new HashSet<Atom>();
		
		Atom root = atomSet.get(0); 
		Atom parentAtom = null;
		Map<Atom, Atom> atomToParentMap = new HashMap<Atom, Atom>();
		Set<Bond> linkBondSet = new LinkedHashSet<Bond>(); 
		
		expand(root, parentAtom, usedAtoms, atomToParentMap, linkBondSet);
		
		for (Bond bond : linkBondSet) {
			Ring ring = getRing(bond, atomToParentMap);
			ringList.add(ring);
		}
		
		return ringList;
	}
	
	private Ring getRing(Bond bond, Map<Atom, Atom> atomToParentMap) throws StructureBuildingException { 
		Atom atomFrom =  bond.getFromAtom() ;
		Atom atomTo = bond.getToAtom(); 
		ArrayList<Atom>  atomSet0 = getAncestors(atomFrom, atomToParentMap);
		ArrayList<Atom>  atomSet1 = getAncestors(atomTo, atomToParentMap);
		ArrayList<Bond>  bondSet0 = getAncestors1(atomFrom, atomToParentMap, atomSet1);
		ArrayList<Bond>  bondSet1 = getAncestors1(atomTo, atomToParentMap, atomSet0);
		ArrayList<Bond>  mergedBondSet = symmetricDifference (bondSet0, bondSet1); 
		
		mergedBondSet.add(bond);
		Ring ring = new Ring(mergedBondSet);
		return ring;
	}
	
	
	private ArrayList<Atom>  getAncestors(Atom atom, Map<Atom, Atom> atomToParentMap) {
		ArrayList<Atom> newAtomSet = new ArrayList<Atom> ();

		atom = (Atom) atomToParentMap.get(atom);
		if (atom != null && !newAtomSet.contains(atom)) {
			newAtomSet.add(atom);
		}

		return newAtomSet;
	}
	
	private ArrayList<Bond>  getAncestors1(Atom atom, Map<Atom, Atom> atomToParentMap, ArrayList<Atom> atomSet) throws StructureBuildingException {
		Fragment molecule =atom.getFrag();
		ArrayList<Bond> newBondSet = new ArrayList<Bond>();
		while (true) {
			Atom atom1 = (Atom) atomToParentMap.get(atom);
			if (atom1 == null) {
				break;
			}
			Bond bond = molecule.findBondOrThrow(atom, atom1);
			if (newBondSet.contains(bond)) {
				break;
			}
			newBondSet.add(bond);
			atom = atom1;
		}
		return newBondSet;
	}
	
	private void expand(Atom atom, Atom parentAtom, 
			Set<Atom> usedAtoms, Map<Atom, Atom> atomToParentMap,
			Set<Bond> linkBondSet) throws StructureBuildingException {
		
		usedAtoms.add(atom);
		atomToParentMap.put(atom, parentAtom);
		List<Atom> ligandAtomList = atom.getAtomNeighbours();
		Fragment fragment = atom.getFrag();
				
		for (int i = 0; i < ligandAtomList.size(); i++) {
			Atom ligandAtom = ligandAtomList.get(i);
			if (ligandAtom.equals(parentAtom)) {
				// skip existing bond
			} else if (usedAtoms.contains(ligandAtom)) {
				Bond linkBond = fragment.findBondOrThrow(atom, ligandAtom);
				linkBondSet.add(linkBond);
				// already treated
			} else {
				expand(ligandAtom, atom, usedAtoms, atomToParentMap, linkBondSet);
			}
		}
	}
	
	
	/**
	 * @param ring
	 * @throws StructureBuildingException 
	 */
	private boolean reduceRingSizes(Ring ring, List<Ring> newList) throws StructureBuildingException {
		boolean change = false;
		for (int i = 0; i < newList.size(); i++) {
			Ring target = newList.get(i);
			if (target == ring) {
				continue;
			}
			
			ArrayList<Bond> newBondSet = symmetricDifference ( target.getBondSet(), ring.getBondSet() ) ;
			if (newBondSet.size() < target.size()) {
				Ring newRing = new Ring(newBondSet);
				newList.set(i, newRing);
				change = true;
			}
		}
		return change;
	}

	
	 public ArrayList<Bond> symmetricDifference(ArrayList<Bond> bondSet1, ArrayList<Bond> bondSet2) {
	 	 ArrayList<Bond> newBondSet = new ArrayList<Bond>();
	        
        for (int i = 0; i < bondSet1.size(); i++) {
            if (!bondSet2.contains(bondSet1.get(i))) {
                newBondSet.add(bondSet1.get(i));
            }
        }
        for (int i = 0; i < bondSet2.size(); i++) {
            Bond bond = bondSet2.get(i);
            if (!bondSet1.contains(bond)) {
                newBondSet.add(bond);
            }
        }

	    return newBondSet;
	 }

	

	/**
	 * Given a fragment returns it's atom list sorted by locant. e.g. 1,2,3,3a,3b,4
	 * @param fragment
	 * @return
	 */
	private List<Atom> sortFragmentAtomListByLocant(Fragment frag) {
		List<Atom> atomsInFragment =frag.getAtomList();
		Collections.sort(atomsInFragment, new SortLocants());
		return atomsInFragment;
	}
	
	/**
	 * Uses locants in front of the benz/benzo group to assign heteroatoms on the now numbered and used fused ring system
	 * @param state
	 * @param benzoEl
	 * @param fusedRing
	 * @throws StructureBuildingException 
	 * @throws PostProcessingException 
	 */
	private void benzoSpecificAssignHeteroAtomsUsingLocants(BuildState state, Element benzoEl, Fragment fusedRing) throws StructureBuildingException, PostProcessingException {
		Element previous = (Element) XOMTools.getPreviousSibling(benzoEl);
		if (previous!=null && previous.getLocalName().equals("multiplier")){//e.g. dibenzothiophene
			throw new StructureBuildingException("multiple benzo groups cannot currently be fused");
		}
		LinkedList<Element> locants =new LinkedList<Element>();
		while(previous != null && previous.getLocalName().equals("locant")) {
			locants.add(previous);
			previous=(Element) XOMTools.getPreviousSibling(previous);
		}
		if (locants.size() >0){
			Elements suffixes=((Element)benzoEl.getParent()).getChildElements("suffix");
			int suffixesWithoutLocants =0;
			for (int i = 0; i < suffixes.size(); i++) {
				if (suffixes.get(i).getAttribute("locant")==null){
					suffixesWithoutLocants++;
				}
			}
			if (locants.size() != suffixesWithoutLocants){//In preference locants will be assigned to suffixes rather than to this nomenclature
				List<Atom> atomList =fusedRing.getAtomList();
				LinkedList<Atom> heteroatoms =new LinkedList<Atom>();
				LinkedList<String> elementOfHeteroAtom =new LinkedList<String>();
				for (Atom atom : atomList) {//this iterates in the same order as the numbering system
					if (!atom.getElement().equals("C")){
						heteroatoms.add(atom);
						elementOfHeteroAtom.add(atom.getElement());
					}
				}
				if (locants.size() >=heteroatoms.size()){//atleast as many locants as there are heteroatoms to assign
					for (Atom atom : heteroatoms) {
						atom.setElement("C");
					}
					fusedRing.pickUpIndicatedHydrogen();
					for (int i=0; i< heteroatoms.size(); i ++){
						String elementSymbol =elementOfHeteroAtom.removeLast();
						Element locant =locants.removeFirst();
						fusedRing.getAtomByLocantOrThrow(locant.getAttributeValue("value")).setElement(elementSymbol);
						locant.detach();
					}
				}
				else if (locants.size()>1){
					throw new PostProcessingException("Unable to assign all locants to benzo-fused ring or multiplier was mising");
				}
			}
		}
	}
}
