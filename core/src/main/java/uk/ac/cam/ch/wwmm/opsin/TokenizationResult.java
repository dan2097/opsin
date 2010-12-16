package uk.ac.cam.ch.wwmm.opsin;

class TokenizationResult {

	private Parse parse;
	private String uninterpretableName;
	private String unparsableName;
	private String unparsedName;

	TokenizationResult(String name){
		this.parse = new Parse(name);
		this.uninterpretableName = "";
		this.unparsableName = "";
		this.unparsedName = name;
	}

	boolean isSuccessfullyTokenized() {
		return unparsedName.isEmpty();
	}

	Parse getParse() {
		return parse;
	}

	String getUninterpretableName() {
		return uninterpretableName;
	}

	boolean isFullyInterpretable() {
		return "".equals(uninterpretableName);
	}

	public void setUninterpretableName(String name) {
		this.uninterpretableName = name;
	}

	String getUnparsableName() {
		return unparsableName;
	}

	void setUnparsableName(String name) {
		this.unparsableName = name;
	}

	String getUnparsedName() {
		return unparsedName;
	}

	String getParsedName() {
		return unparsedName.substring(0, unparsedName.length() - uninterpretableName.length());
	}

	void setUnparsedName(String name) {
		this.unparsedName = name;
	}
}
