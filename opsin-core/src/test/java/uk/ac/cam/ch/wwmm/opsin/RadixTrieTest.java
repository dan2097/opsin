package uk.ac.cam.ch.wwmm.opsin;

import static org.junit.Assert.*;

import java.util.List;

import org.junit.Test;


public class RadixTrieTest {

	@Test
	public void testSimpleAddSimpleGet(){
		OpsinRadixTrie trie = new OpsinRadixTrie();
		trie.addToken("benzene");
		List<Integer> matches= trie.findMatches("benzene", 0);
		assertNotNull(matches);
		assertEquals(1, matches.size());
		assertEquals(7, matches.get(0).intValue());
	}
	
	@Test
	public void testSimpleAddFindPrefix(){
		OpsinRadixTrie trie = new OpsinRadixTrie();
		trie.addToken("phenyl");
		List<Integer> matches= trie.findMatches("phenylbenzene", 0);
		assertNotNull(matches);
		assertEquals(1, matches.size());
		assertEquals(6, matches.get(0).intValue());
	}
	
	@Test
	public void testAddWithBranchFindPrefix(){
		OpsinRadixTrie trie = new OpsinRadixTrie();
		trie.addToken("pyridinyl");
		trie.addToken("phenyl");
		List<Integer> matches= trie.findMatches("phenylbenzene", 0);
		assertNotNull(matches);
		assertEquals(1, matches.size());
		assertEquals(6, matches.get(0).intValue());
	}
	
	@Test
	public void testZeroLengthToken(){
		OpsinRadixTrie trie = new OpsinRadixTrie();
		trie.addToken("");//e.g. end of substituent
		List<Integer> matches= trie.findMatches("phenylbenzene", 0);
		assertNotNull(matches);
		assertEquals(1, matches.size());
		assertEquals(0, matches.get(0).intValue());
	}
	
	@Test
	public void testMultipleHits(){
		OpsinRadixTrie trie = new OpsinRadixTrie();
		trie.addToken("methyl");
		trie.addToken("methylidene");
		List<Integer> matches= trie.findMatches("methylidene", 0);
		assertNotNull(matches);
		assertEquals(2, matches.size());
		assertEquals(6, matches.get(0).intValue());
		assertEquals(11, matches.get(1).intValue());
	}
	
	@Test
	public void testMultipleHits2(){
		OpsinRadixTrie trie = new OpsinRadixTrie();
		trie.addToken("abcdef");
		trie.addToken("a");
		trie.addToken("");
		trie.addToken("acd");
		trie.addToken("ab");
		trie.addToken("abcf");
		List<Integer> matches= trie.findMatches("abc", 0);
		assertNotNull(matches);
		assertEquals(3, matches.size());
		assertEquals(0, matches.get(0).intValue());
		assertEquals(1, matches.get(1).intValue());
		assertEquals(2, matches.get(2).intValue());
	}
	
	
	@Test
	public void testReverseMatching(){
		OpsinRadixTrie trie = new OpsinRadixTrie();
		trie.addToken("enedilyhte");
		trie.addToken("lyhte");
		trie.addToken("");
		trie.addToken("ly");
		trie.addToken("lyhtem");
		List<Integer> matches= trie.findMatchesReadingStringRightToLeft("ethyl", 5);
		assertNotNull(matches);
		
		assertEquals(3, matches.size());
		assertEquals(5, matches.get(0).intValue());
		assertEquals(3, matches.get(1).intValue());
		assertEquals(0, matches.get(2).intValue());
	}
}
	
