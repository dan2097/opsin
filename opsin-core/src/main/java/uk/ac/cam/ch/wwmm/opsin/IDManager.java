package uk.ac.cam.ch.wwmm.opsin;

/**A source of unique integers. Starts at 1 by default.
 *
 * @author ptc24
 *
 */
class IDManager {
	/**the last integer generated, or 0 at first*/
	private int currentID;

	int getCurrentID() {
		return currentID;
	}

	/**Initialises currentID at zero - will give 1 when first called */
	IDManager() {
		currentID = 0;
	}

	/**Generates a new, unique integer. This is one
	 * higher than the previous integer, or 1 if previously uncalled.
	 * @return The generated integer.
	 */
	int getNextID() {
		currentID += 1;
		return currentID;
	}

}
