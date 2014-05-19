package uk.ac.cam.ch.wwmm.opsin;

class TokenGroupEl extends TokenEl {
	
	private Fragment frag;

	TokenGroupEl(String name) {
		super(name);
	}
	
	TokenGroupEl(String name, String value) {
		super(name, value);
	}

	@Override
	Fragment getFrag() {
		return frag;
	}

	void setFrag(Fragment frag) {
		this.frag = frag;
	}
}
