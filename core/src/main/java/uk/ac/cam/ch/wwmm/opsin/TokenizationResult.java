package uk.ac.cam.ch.wwmm.opsin;

class TokenizationResult {

	private Parse parse;
	private String workingName;
	private String unparsableName;
	private String unparsedName;
	private String uninterpretableName;
	public String foo;

	TokenizationResult(String name) {
		this.parse = new Parse(name);
		this.workingName = "";
		this.unparsableName = "";
		this.unparsedName = name;
		this.uninterpretableName = "";
	}

	boolean isSuccessfullyTokenized() {
		return unparsedName.isEmpty();
	}

	Parse getParse() {
		return parse;
	}

	void setUninterpretableName(String name) {
		this.uninterpretableName = name;
	}

	String getUninterpretableName() {
		return this.uninterpretableName;
	}

	String getWorkingName() {
		return workingName;
	}

	public void setWorkingName(String name) {
		this.workingName = name;
	}

	boolean isFullyInterpretable() {
		return "".equals(workingName);
	}

	String getUnparsableName() {
		return unparsableName;
	}

	void setUnparsableName(String name) {
		if (!"".equals(this.unparsableName)) {
			if (!this.unparsableName.replace(" ", "").equals(name)) {
						this.unparsableName = name;
			}
		} else {
			this.unparsableName = name;
		}
	}

	String getUnparsedName() {
		return unparsedName;
	}

	String getParsedName() {
		return unparsedName.substring(0, unparsedName.length() - workingName.length());
	}

	void setUnparsedName(String name) {
		this.unparsedName = name;
	}
}
