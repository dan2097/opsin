package uk.ac.cam.ch.wwmm.opsin;

import java.util.Iterator;
import java.util.List;

/**
 * Convenience class for iterating over a list of atoms that form a ring
 * Doing next on the final atom in the ring will return the first
 * @author dl387
 *
 */
class RingAtomIterator implements Iterator<Atom> {
	final List<Atom> ringAtoms;
	final int startingPosition;
	private boolean started = false;
	private int indice;

	public RingAtomIterator(List<Atom> ringAtoms, int startingPosition) throws StructureBuildingException {
		if (startingPosition>=ringAtoms.size()){
			throw new StructureBuildingException("Starting position is not within ringAtom list");
		}
		this.startingPosition = startingPosition;
		this.indice = startingPosition -1 ;
		this.ringAtoms = ringAtoms;
	}
	public boolean hasNext() {
		int tempIndice = indice +1;
		if (tempIndice>=ringAtoms.size()){
			tempIndice=0;
		}
		return (tempIndice < ringAtoms.size() && (!started || tempIndice != startingPosition));
	}

	public Atom next() {
		int tempIndice = indice + 1;
		if (tempIndice>=ringAtoms.size()){
			tempIndice=0;
		}
		if (tempIndice >= ringAtoms.size() || (started && tempIndice == startingPosition)){//gone through all positions in ring.
			return null;
		}
		else{
			indice =tempIndice;
			started =true;
		}
		return ringAtoms.get(indice);
	}

	public void remove() {
		ringAtoms.remove(indice);
	}
	
	public void reset() {
		indice = startingPosition -1;
		started = false;
	}
}
