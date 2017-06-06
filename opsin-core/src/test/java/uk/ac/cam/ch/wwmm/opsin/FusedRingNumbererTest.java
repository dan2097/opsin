package uk.ac.cam.ch.wwmm.opsin;

import static org.junit.Assert.*;
import static org.mockito.Mockito.mock;

import java.util.List;

import org.junit.Ignore;
import org.junit.Test;

/**
 * Tests that fused ring numbering is working as expected. A heteroatom(n) has been placed at the expected locant 1 to make numbering unambiguous where due to symmetry geometric consideration are insufficient to deduce unique numbering
 * Currently interior atoms are not labelled. As this is not seen as a problem, any tests of compounds with interior atoms have not had locants assigned to the interior atoms
 * @author dl387
 *
 */
public class FusedRingNumbererTest {

	private SMILESFragmentBuilder sBuilder = new SMILESFragmentBuilder(new IDManager());

	@Test
	public void aceanthrene() throws StructureBuildingException {
		compareNumbering("C1Cc2cccc3cc4ccccc4c1c23", "1/2/2a/3/4/5/5a/6/6a/7/8/9/10/10a/10b/10c");
	}

	@Test
	public void acenaphthene() throws StructureBuildingException {
		compareNumbering("C1Cc2cccc3cccc1c23", "1/2/2a/3/4/5/5a/6/7/8/8a/8b");
	}

	@Test
	public void acephenanthrene() throws StructureBuildingException {
		compareNumbering("c1ccc2CCc3cc4ccccc4c1c23", "1/2/3/3a/4/5/5a/6/6a/7/8/9/10/10a/10b/10c");
	}

	@Test
	public void arsanthrene() throws StructureBuildingException {
		compareNumbering("C1=CC=CC=2[As]=C3C=CC=CC3=[As]C12", "1/2/3/4/4a/5/5a/6/7/8/9/9a/10/10a");
	}

	@Test
	public void arsanthridine() throws StructureBuildingException {
		compareNumbering("c1cccc2[as]cc3ccccc3c12", "1/2/3/4/4a/5/6/6a/7/8/9/10/10a/10b");
	}

	@Test
	public void arsindole() throws StructureBuildingException {
		compareNumbering("[as]1ccc2ccccc12", "1/2/3/3a/4/5/6/7/7a");
	}

	@Test
	public void arsindoline() throws StructureBuildingException {
		compareNumbering("[as]1cccc2ccccc12", "1/2/3/4/4a/5/6/7/8/8a");
	}

	@Test
	public void betacarboline() throws StructureBuildingException {
		compareNumbering("c1nccc2c3ccccc3nc12", "1/2/3/4/4a/4b/5/6/7/8/8a/9/9a");
	}

	@Test
	public void boranthrene() throws StructureBuildingException {
		compareNumbering("C1=CC=CC=2[B]=C3C=CC=CC3=[B]C12", "1/2/3/4/4a/5/5a/6/7/8/9/9a/10/10a");
	}

	@Test
	public void cholanthrene() throws StructureBuildingException {
		compareNumbering("C1Cc2cccc3cc4c5ccccc5ccc4c1c23", "1/2/2a/3/4/5/5a/6/6a/6b/7/8/9/10/10a/11/12/12a/12b/12c");
	}

	@Test
	public void thiochromane() throws StructureBuildingException {
		compareNumbering("S1CCCc2ccccc12", "1/2/3/4/4a/5/6/7/8/8a");
	}

	@Test
	public void selenochromane() throws StructureBuildingException {
		compareNumbering("[Se]1CCCc2ccccc12", "1/2/3/4/4a/5/6/7/8/8a");
	}

	@Test
	public void tellurochromane() throws StructureBuildingException {
		compareNumbering("[Te]1CCCc2ccccc12", "1/2/3/4/4a/5/6/7/8/8a");
	}

	@Test
	public void thiochromene() throws StructureBuildingException {
		compareNumbering("s1cccc2ccccc12", "1/2/3/4/4a/5/6/7/8/8a");
	}

	@Test
	public void selenochromene() throws StructureBuildingException {
		compareNumbering("[se]1cccc2ccccc12", "1/2/3/4/4a/5/6/7/8/8a");
	}

	@Test
	public void tellurochromene() throws StructureBuildingException {
		compareNumbering("[Te]1CC=Cc2ccccc12", "1/2/3/4/4a/5/6/7/8/8a");
	}

	@Test
	public void coronene() throws StructureBuildingException {
		compareNumbering("c1cc2ccc3ccc4ccc5ccc6ccc71.c28c3c4c5c6c78", "1/2/2a/3/4/4a/5/6/6a/7/8/8a/9/10/10a/11/12/12a/12b/12c/12d/12e/12f/12g");
	}
	
	@Test
	public void indane() throws StructureBuildingException {
		compareNumbering("C1CCc2ccccc12", "1/2/3/3a/4/5/6/7/7a");
	}

	@Test
	public void isoarsindole() throws StructureBuildingException {
		compareNumbering("c1[as]cc2ccccc12", "1/2/3/3a/4/5/6/7/7a");
	}

	@Test
	public void isoarsinoline() throws StructureBuildingException {
		compareNumbering("c1[as]ccc2ccccc12", "1/2/3/4/4a/5/6/7/8/8a");
	}

	@Test
	public void thioisochromane() throws StructureBuildingException {
		compareNumbering("C1SCCc2ccccc12", "1/2/3/4/4a/5/6/7/8/8a");
	}

	@Test
	public void selenoisochromane() throws StructureBuildingException {
		compareNumbering("C1[Se]CCc2ccccc12", "1/2/3/4/4a/5/6/7/8/8a");
	}

	@Test
	public void telluroisochromane() throws StructureBuildingException {
		compareNumbering("C1[Te]CCc2ccccc12", "1/2/3/4/4a/5/6/7/8/8a");
	}

	@Test
	public void isochromene() throws StructureBuildingException {
		compareNumbering("c1occc2ccccc12", "1/2/3/4/4a/5/6/7/8/8a");
	}

	@Test
	public void thioisochromene() throws StructureBuildingException {
		compareNumbering("c1sccc2ccccc12", "1/2/3/4/4a/5/6/7/8/8a");
	}

	@Test
	public void selenoisochromene() throws StructureBuildingException {
		compareNumbering("c1[se]ccc2ccccc12", "1/2/3/4/4a/5/6/7/8/8a");
	}

	@Test
	public void telluroisochromene() throws StructureBuildingException {
		compareNumbering("C1[Te]C=Cc2ccccc12", "1/2/3/4/4a/5/6/7/8/8a");
	}

	@Test
	public void isophosphindole() throws StructureBuildingException {
		compareNumbering("c1pcc2ccccc12", "1/2/3/3a/4/5/6/7/7a");
	}

	@Test
	public void isophosphinoline() throws StructureBuildingException {
		compareNumbering("c1pccc2ccccc12", "1/2/3/4/4a/5/6/7/8/8a");
	}

	@Test
	public void isoviolanthrene() throws StructureBuildingException {
		compareNumbering("c1cccc2c3ccc4c5ccc6Cc7ccccc7c8ccc9c%10ccc%11Cc12.c3%11c4%10.c59c68", "1/2/3/4/4a/4b/5/6/6a/6b/7/8/8a/9/9a/10/11/12/13/13a/13b/14/15/15a/15b/16/17/17a/18/18a/18b/18c/18d/18e");
	}

	@Test
	public void mercuranthrene() throws StructureBuildingException {
		compareNumbering("c1cccc2[Hg]c3ccccc3[Hg]c12", "1/2/3/4/4a/5/5a/6/7/8/9/9a/10/10a");
	}

	@Test
	@Ignore //Transient BUG in path finding code
	public void ovalene() throws StructureBuildingException {
		compareNumbering("c1cc2ccc3ccc4cc5ccc6ccc7ccc8cc91.c19c2c3c4c9c5c6c7c8c19", "1/2/2a/3/4/4a/5/6/6a/7/7a/8/9/9a/10/11/11a/12/13/13a/14/14a/14b/14c/14d/14e/14f/14g/14h/14i/14j/14k");
	}

	@Test
	public void oxanthrene() throws StructureBuildingException {
		compareNumbering("c1cccc2Oc3ccccc3Oc12", "1/2/3/4/4a/5/5a/6/7/8/9/9a/10/10a");
	}

	@Test
	public void perylene() throws StructureBuildingException {
		compareNumbering("c1ccc2cccc3c4cccc5cccc6c71.c237.c456", "1/2/3/3a/4/5/6/6a/6b/7/8/9/9a/10/11/12/12a/12b/12c/12d");
	}

	@Test
	public void phenanthridine() throws StructureBuildingException {
		compareNumbering("c1cccc2ncc3ccccc3c12", "1/2/3/4/4a/5/6/6a/7/8/9/10/10a/10b");
	}

	@Test
	public void phenomercurine() throws StructureBuildingException {
		compareNumbering("c1cccc2[Hg]c3ccccc3[Hg]c12", "1/2/3/4/4a/5/5a/6/7/8/9/9a/10/10a");
	}

	@Test
	public void phenoxazine() throws StructureBuildingException {
		compareNumbering("c1cccc2Oc3ccccc3nc12", "1/2/3/4/4a/5/5a/6/7/8/9/9a/10/10a");
	}

	@Test
	public void phenothiazine() throws StructureBuildingException {
		compareNumbering("c1cccc2Sc3ccccc3nc12", "1/2/3/4/4a/5/5a/6/7/8/9/9a/10/10a");
	}

	@Test
	public void phenoselenazine() throws StructureBuildingException {
		compareNumbering("c1cccc2[Se]c3ccccc3nc12", "1/2/3/4/4a/5/5a/6/7/8/9/9a/10/10a");
	}

	@Test
	public void phenotellurazine() throws StructureBuildingException {
		compareNumbering("c1cccc2[Te]c3ccccc3nc12", "1/2/3/4/4a/5/5a/6/7/8/9/9a/10/10a");
	}

	@Test
	public void phenophosphazinine() throws StructureBuildingException {
		compareNumbering("c1cccc2nc3ccccc3pc12", "1/2/3/4/4a/5/5a/6/7/8/9/9a/10/10a");
	}

	@Test
	public void phenophosphazine() throws StructureBuildingException {
		compareNumbering("c1cccc2nc3ccccc3pc12", "1/2/3/4/4a/5/5a/6/7/8/9/9a/10/10a");
	}

	@Test
	public void phenarsazinine() throws StructureBuildingException {
		compareNumbering("c1cccc2nc3ccccc3[as]c12", "1/2/3/4/4a/5/5a/6/7/8/9/9a/10/10a");
	}

	@Test
	public void phenoarsazine() throws StructureBuildingException {
		compareNumbering("c1cccc2nc3ccccc3[as]c12", "1/2/3/4/4a/5/5a/6/7/8/9/9a/10/10a");
	}

	@Test
	public void phenomercurazine() throws StructureBuildingException {
		compareNumbering("c1cccc2nc3ccccc3[Hg]c12", "1/2/3/4/4a/5/5a/6/7/8/9/9a/10/10a");
	}

	@Test
	public void phenomercazine() throws StructureBuildingException {
		compareNumbering("c1cccc2nc3ccccc3[Hg]c12", "1/2/3/4/4a/5/5a/6/7/8/9/9a/10/10a");
	}

	@Test
	public void phenoxathiine() throws StructureBuildingException {
		compareNumbering("c1cccc2Oc3ccccc3Sc12", "1/2/3/4/4a/5/5a/6/7/8/9/9a/10/10a");
	}

	@Test
	public void phenoxaselenine() throws StructureBuildingException {
		compareNumbering("c1cccc2Oc3ccccc3[Se]c12", "1/2/3/4/4a/5/5a/6/7/8/9/9a/10/10a");
	}

	@Test
	public void phenoxatellurine() throws StructureBuildingException {
		compareNumbering("c1cccc2Oc3ccccc3[Te]c12", "1/2/3/4/4a/5/5a/6/7/8/9/9a/10/10a");
	}

	@Test
	public void phenoxaphosphinine() throws StructureBuildingException {
		compareNumbering("c1cccc2Oc3ccccc3pc12", "1/2/3/4/4a/5/5a/6/7/8/9/9a/10/10a");
	}

	@Test
	public void phenoxaphosphine() throws StructureBuildingException {
		compareNumbering("c1cccc2Oc3ccccc3pc12", "1/2/3/4/4a/5/5a/6/7/8/9/9a/10/10a");
	}

	@Test
	public void phenoxarsinine() throws StructureBuildingException {
		compareNumbering("c1cccc2Oc3ccccc3[as]c12", "1/2/3/4/4a/5/5a/6/7/8/9/9a/10/10a");
	}

	@Test
	public void phenoxarsine() throws StructureBuildingException {
		compareNumbering("c1cccc2Oc3ccccc3[as]c12", "1/2/3/4/4a/5/5a/6/7/8/9/9a/10/10a");
	}

	@Test
	public void phenoxastibinine() throws StructureBuildingException {
		compareNumbering("c1cccc2Oc3ccccc3[sb]c12", "1/2/3/4/4a/5/5a/6/7/8/9/9a/10/10a");
	}

	@Test
	public void phenoxantimonine() throws StructureBuildingException {
		compareNumbering("c1cccc2Oc3ccccc3[sb]c12", "1/2/3/4/4a/5/5a/6/7/8/9/9a/10/10a");
	}

	@Test
	public void phenothiarsinine() throws StructureBuildingException {
		compareNumbering("c1cccc2Sc3ccccc3[as]c12", "1/2/3/4/4a/5/5a/6/7/8/9/9a/10/10a");
	}

	@Test
	public void phenothiarsine() throws StructureBuildingException {
		compareNumbering("c1cccc2Sc3ccccc3[as]c12", "1/2/3/4/4a/5/5a/6/7/8/9/9a/10/10a");
	}

	@Test
	public void phosphanthrene() throws StructureBuildingException {
		compareNumbering("C1=CC=CC=2P=C3C=CC=CC3=PC12", "1/2/3/4/4a/5/5a/6/7/8/9/9a/10/10a");
	}

	@Test
	public void phosphindole() throws StructureBuildingException {
		compareNumbering("p1ccc2ccccc12", "1/2/3/3a/4/5/6/7/7a");
	}

	@Test
	public void phosphinoline() throws StructureBuildingException {
		compareNumbering("p1cccc2ccccc12", "1/2/3/4/4a/5/6/7/8/8a");
	}

	@Test
	public void picene() throws StructureBuildingException {
		compareNumbering("c1cccc2ccc3c4ccc5ccccc5c4ccc3c21", "1/2/3/4/4a/5/6/6a/6b/7/8/8a/9/10/11/12/12a/12b/13/14/14a/14b");
	}

	@Test
	public void pleiadene() throws StructureBuildingException {
		compareNumbering("c1ccc2cccc3cc4ccccc4cc51.c235", "1/2/3/3a/4/5/6/6a/7/7a/8/9/10/11/11a/12/12a/12b");
	}

	@Test
	public void pyranthrene() throws StructureBuildingException {
		compareNumbering("c1cccc2c3cc4ccc5cc6ccccc6c6cc7ccc8cc12.c38c7c4c56", "1/2/3/4/4a/4b/5/5a/6/7/7a/8/8a/9/10/11/12/12a/12b/13/13a/14/15/15a/16/16a////");
	}

	@Test
	public void pyrrolizine() throws StructureBuildingException {
		compareNumbering("c1ccn2cccc12", "1/2/3/4/5/6/7/7a");
	}

	@Test
	public void quinolizine() throws StructureBuildingException {
		compareNumbering("c1cccn2ccccc12", "1/2/3/4/5/6/7/8/9/9a");
	}

	@Test
	public void rubicene() throws StructureBuildingException {
		compareNumbering("c1ccc2c3ccccc3c3c4cccc5c6ccccc6c6c71.c237.c456", "1/2/3/3a/3b/4/5/6/7/7a/7b/7c/8/9/10/10a/10b/11/12/13/14/14a/14b/14c/14d/14e");
	}

	@Test
	public void silanthrene() throws StructureBuildingException {
		compareNumbering("C1=CC=CC=2[Si]=C3C=CC=CC3=[Si]C12", "1/2/3/4/4a/5/5a/6/7/8/9/9a/10/10a");
	}

	@Test
	public void selenanthrene() throws StructureBuildingException {
		compareNumbering("c1cccc2[Se]c3ccccc3[Se]c12", "1/2/3/4/4a/5/5a/6/7/8/9/9a/10/10a");
	}

	@Test
	public void telluranthrene() throws StructureBuildingException {
		compareNumbering("c1cccc2[Te]c3ccccc3[Te]c12", "1/2/3/4/4a/5/5a/6/7/8/9/9a/10/10a");
	}


	@Test
	public void thianthrene() throws StructureBuildingException {
		compareNumbering("c1cccc2Sc3ccccc3Sc12", "1/2/3/4/4a/5/5a/6/7/8/9/9a/10/10a");
	}

	@Test
	public void trindene() throws StructureBuildingException {
		compareNumbering("c1ccc2c3cccc3c4cccc4c12", "1/2/3/3a/3b/4/5/6/6a/6b/7/8/9/9a/9b");
	}

	@Test
	public void violanthrene() throws StructureBuildingException {
		compareNumbering("c1cccc2Cc3ccc4c5ccc6Cc7ccccc7c8ccc9c%10ccc%11c12.c3%11c4%10.c59c68", "1/2/3/4/4a/5/5a/6/7/7a/7b/8/9/9a/10/10a/11/12/13/14/14a/14b/15/16/16a/16b/17/18/18a/18b////");
	}

	@Test
	public void naphthotetraphene() throws StructureBuildingException {
		compareNumbering("c1cc2ccc3ccc4cc5ccccc5cc4c3c2c6ccccc16", "1/2/2a/3/4/4a/5/6/6a/7/7a/8/9/10/11/11a/12/12a/12b/12c/12d/13/14/15/16/16a");
	}
	
	@Test
	public void anthratetraphene() throws StructureBuildingException {
		compareNumbering("c1cc2ccc3ccc4cc5ccccc5cc4c3c2c6cc7ccccc7cc16", "1/2/2a/3/4/4a/5/6/6a/7/7a/8/9/10/11/11a/12/12a/12b/12c/12d/13/13a/14/15/16/17/17a/18/18a");
	}
	
	@Test
	public void octalenotetraphene() throws StructureBuildingException {
		compareNumbering("c1ccc2ccc3ccc4cc5ccccc5cc4c3c2cc6ccccccc16", "1/2/3/3a/4/5/5a/6/7/7a/8/8a/9/10/11/12/12a/13/13a/13b/13c/14/14a/15/16/17/18/19/20/20a");
	}
	
	@Test
	public void difficultChain() throws StructureBuildingException {
		compareNumbering("C1C2C3C4CC5CC6CCCCCCC6CCC5CC4CC3C12", "1/1a/1b/1c/2/2a/3/3a/4/5/6/7/8/9/9a/10/11/11a/12/12a/13/13a/13b");
	}
	
	@Test
	public void difficultChain2() throws StructureBuildingException {
		compareNumbering("C1C2CCCCC2C3C4C5C6CCC7C8CCCCC8C7CCC56C4CCC13", "1/1a/2/3/4/5/5a/5b/5c/5d/5e/6/7/7a/7b/8/9/10/11/11a/11b/12/13/13a/13b/14/15/15a");
	}
	
	@Test
	public void acrindoline() throws StructureBuildingException {
		compareNumbering("c1cccc2c3ccc4cc5ccccc5nc4c3nc12", "1/2/3/4/4a/4b/5/6/6a/7/7a/8/9/10/11/11a/12/12a/12b/13/13a");
	}
	
	@Test
	public void anthrazine() throws StructureBuildingException {
		compareNumbering("c1cccc2cc3c4nc5ccc6cc7ccccc7cc6c5nc4ccc3cc12", "1/2/3/4/4a/5/5a/5b/6/6a/7/8/8a/9/9a/10/11/12/13/13a/14/14a/14b/15/15a/16/17/17a/18/18a");
	}
	
	@Test
	public void anthyridine() throws StructureBuildingException {
		compareNumbering("n1cccc2cc3cccnc3nc12", "1/2/3/4/4a/5/5a/6/7/8/9/9a/10/10a");
	}
	
	@Test
	public void benzo_cd_azulene() throws StructureBuildingException {
		compareNumbering("c1cc2cccc3ccccc1c23", "1/2/2a/3/4/5/5a/6/7/8/9/9a/9b");
	}
	
	@Test
	public void indeno_7_1_cd_azepine() throws StructureBuildingException {
		compareNumbering("c1nccc2ccc3cccc1c23", "1/2/3/4/4a/5/6/6a/7/8/9/9a/9b");
	}
	
	@Test
	@Ignore
	public void tripleSubstituedSevenMembered() throws StructureBuildingException {
		compareNumbering("C1NCCN2c3ncccc3Cc4ccccc4C12", "1/2/3/4/5/5a/6/7/8/9/9a/10/10a/11/12/13/14/14a/14b");
		compareNumbering("c1cccc2C3CNCCN3c4ncccc4Cc12", "1/2/3/4/5/5a/6/7/8/9/9a/10/10a/11/12/13/14/14a/14b");
	}

	/**
	 * Takes smiles and expected labels for a fused ring. Generates the fused ring, numbers it then compares to the given slash delimited labels
	 * @param smiles
	 * @param labels
	 * @throws StructureBuildingException 
	 */
	private void compareNumbering(String smiles, String labels) throws StructureBuildingException {
		Fragment fusedRing = sBuilder.build(smiles, mock(Element.class), XmlDeclarations.NONE_LABELS_VAL);
		String[] labelArray =labels.split("/", -1);
		FusedRingNumberer.numberFusedRing(fusedRing);
		List<Atom> atomList =fusedRing.getAtomList();
		assertEquals(atomList.size(), labelArray.length);//bug in test if not true!
		for (int i = 0; i < atomList.size(); i++) {
			if (!labelArray[i].equals("")){//exterior atom locant
				assertEquals(labelArray[i],atomList.get(i).getFirstLocant());
			}
		}
	}
}
