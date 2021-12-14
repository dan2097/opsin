package uk.ac.cam.ch.wwmm.opsin;

import java.util.Objects;

class StereoGroup implements Comparable<StereoGroup> {

	private final StereoGroupType type;
	private final int number;
	
	StereoGroup(StereoGroupType type) {
		this(type, 1);
	}

	StereoGroup(StereoGroupType type, int number) {
		this.type = type;
		this.number = number;
	}

	@Override
	public boolean equals(Object o) {
		if (this == o)
			return true;
		if (o == null || getClass() != o.getClass())
			return false;
		StereoGroup that = (StereoGroup) o;
		return number == that.number && type == that.type;
	}

	@Override
	public int hashCode() {
		return Objects.hash(type, number);
	}

	public int compareTo(StereoGroup that) {
		int cmp = this.type.compareTo(that.type);
		if (cmp != 0) {
			return cmp;
		}
		return Integer.compare(this.number, that.number);
	}

	StereoGroupType getType() {
		return type;
	}

	int getNumber() {
		return number;
	}

	@Override
	public String toString() {
		return "StereoGroup{" + "type=" + type + ", number=" + number + '}';
	}

}
