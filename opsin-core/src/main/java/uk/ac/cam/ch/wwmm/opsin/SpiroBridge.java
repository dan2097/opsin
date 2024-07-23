package uk.ac.cam.ch.wwmm.opsin;

class SpiroBridge {

	
	private final int chainLength;
	private final Integer locant;
	private final boolean hasSuperscriptedLocant;
	
	SpiroBridge(int chainLength, Integer locant, boolean hasSuperscriptedLocant) {
		this.chainLength = chainLength;
		this.locant = locant;
		this.hasSuperscriptedLocant = hasSuperscriptedLocant;
	}

	int getChainLength() {
		return chainLength;
	}

	Integer getLocant() {
		return locant;
	}

	boolean hasSuperscriptedLocant() {
		return hasSuperscriptedLocant;
	}

}
