package uk.ac.cam.ch.wwmm.opsin;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;
import org.junit.runners.Suite.SuiteClasses;

@RunWith(Suite.class)
@SuiteClasses({AtomTest.class,
	BondTest.class,
	CMLFragmentBuilderTest.class,
	DtdTest.class,
	FragmentManagerTest.class,
	FragmentTest.class,
	FusedRingBuilderTest.class,
	NameToStructureTest_FromPapers.class,
	SMILESFragmentBuilderTest.class,
	SSSRTest.class,
	StereochemistryTest.class,
	StructureBuilderTest.class,
	TokenizerTest.class,
})
public class AllTests {

}
