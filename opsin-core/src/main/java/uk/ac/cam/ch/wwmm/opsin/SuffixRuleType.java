package uk.ac.cam.ch.wwmm.opsin;

enum SuffixRuleType {
	addgroup,
	addSuffixPrefixIfNonePresentAndCyclic,
	setOutAtom,
	changecharge,
	addFunctionalAtomsToHydroxyGroups,
	chargeHydroxyGroups,
	removeTerminalOxygen,
	convertHydroxyGroupsToOutAtoms,
	convertHydroxyGroupsToPositiveCharge,
	setAcidicElement
}
