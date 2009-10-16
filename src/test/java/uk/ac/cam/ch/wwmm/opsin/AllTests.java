package uk.ac.cam.ch.wwmm.opsin;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;
import org.junit.runners.Suite.SuiteClasses;

@RunWith(Suite.class)
@SuiteClasses({AtomTest.class,
	BondTest.class,
	CMLFragmentBuilderTest.class,
	FragmentManagerTest.class,
	FragmentTest.class,
	NameToStructureTest_FromPapers.class,
	ParseRulesTest.class,
	SMILESFragmentBuilderTest.class,
	StructureBuilderTest.class,
	FusedRingBuilderTest.class
})
public class AllTests {

}
