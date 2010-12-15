package uk.ac.cam.ch.wwmm.opsin;

class TokenizationResult {

	private Parse parse;
	private String uninterpretableName;
	private String unparsableName;
	private String unparsedName;

	public TokenizationResult(String name, boolean allowRemovalOfWhiteSpace) throws ParsingException {
		this.parse = new Parse(name);
		this.uninterpretableName = "";
		this.unparsableName = "";
		this.unparsedName = allowRemovalOfWhiteSpace ? removeWhiteSpaceIfBracketsAreUnbalanced(name) : name;
	}

	boolean isSuccessfullyTokenized() {
		return this.unparsedName.isEmpty();
	}

	Parse getParse() {
		return parse;
	}

	String getUninterpretableName() {
		return uninterpretableName;
	}

	boolean hasUninterpretableName() {
		return "".equals(this.uninterpretableName);
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
		return this.unparsedName.substring(0, this.unparsedName.length() - this.uninterpretableName.length());
	}

	void setUnparsedName(String name) {
		this.unparsedName = name;
	}

	/**
	 * Works left to right removing spaces if there are too many opening brackets
	 * @param name
	 * @return
	 * @throws ParsingException If brackets are unbalanced and cannot be balanced by removing whitespace
	 */
	private String removeWhiteSpaceIfBracketsAreUnbalanced(String name) throws ParsingException {
		int bracketLevel = 0;
		int stringLength = name.length();
		for (int i = 0; i < stringLength; i++) {
			char c = name.charAt(i);
			if (c == '(' || c == '[' || c == '{') {
				bracketLevel++;
			} else if (c == ')' || c == ']' || c == '}') {
				bracketLevel--;
			} else if (c == ' ' && bracketLevel > 0) {//brackets unbalanced and a space has been encountered!
				name = name.substring(0, i) + name.substring(i + 1);
				stringLength = name.length();
				i--;
			}
		}
		if (bracketLevel > 0) {
			throw new ParsingException("Unmatched opening bracket found in :" + name);
		} else if (bracketLevel < 0) {
			throw new ParsingException("Unmatched closing bracket found in :" + name);
		}
		return name;
	}
}
