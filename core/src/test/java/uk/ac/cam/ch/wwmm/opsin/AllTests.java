package uk.ac.cam.ch.wwmm.opsin;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;
import org.junit.runners.Suite.SuiteClasses;

@RunWith(Suite.class)
@SuiteClasses({AtomTest.class,
	BondTest.class,
	CASToolsTest.class,
	CMLFragmentBuilderTest.class,
	ComponentGeneration_AmbiguitiesAndIrregularities.class,
	ComponentGeneration_ProcesslocantsTest.class,
	ComponentGeneration_StereochemistryTest.class,
	CycleDetectorTest.class,
	DtdTest.class,
	FragmentManagerTest.class,
	FragmentTest.class,
	FusedRingNumbererTest.class,
	HeteroAtomReplacementTest.class,
	NameToStructureConfigurationsTest.class,
	NameToStructureTest.class,
	ParserTest.class,
	PreProcessorTest.class,
	SMILESFragmentBuilderTest.class,
	SMILESWriterTest.class,
	SSSRTest.class,
	StereochemistryTest.class,
	StructureBuilderTest.class,
	TokenizerTest.class,
	VerifyFragmentsTest.class
})
public class AllTests {

}
