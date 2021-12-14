package uk.ac.cam.ch.wwmm.opsin;

/**
 * Indicate if a stereo-center belongs to a certain stereo group type. The group
 * can be specified on the {@link AtomParity}. For formats that support enhanced
 * stereo (currently only CXSMILES) you would normally also have a group number
 * (e.g. Rac1/&amp;1, Rac2/&amp;2, Rac3/&amp;3, etc) however there is currently
 * no way to specify this separation in IUPAC. An extension may be to indicate
 * grouping with parenthesis "(1RS)-,(3RS)-" "((1RS),(3RS))-". However the most
 * common cases of mixed stereo groups are likely to be where part of the
 * structure is known (absolute) and part is racemic/relative which can be
 * specified like this "(1RS,3R)-". The most common cases are all absolute, all
 * racemic (AND enantiomer), all relative (OR enantiomer).
 */
public enum StereoGroupType {
	/**
	 * Absolute stereochemistry, the configuration of the stereo center is known.
	 */
  Abs,
	/**
	 * Racemic stereochemistry, the molecule is a mixture of the stereo center.
	 */
  Rac,
	/**
	 * Relative stereochemistry, the configuration of a stereo center is unknown
	 * but may be known relative to another configuration.
	 */
  Rel,
	/**
	 * Fallback sentinel value to ensure non-null.
	 */
  Unk
}
