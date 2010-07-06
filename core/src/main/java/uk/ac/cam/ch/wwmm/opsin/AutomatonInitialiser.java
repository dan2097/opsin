package uk.ac.cam.ch.wwmm.opsin;

import java.io.FileOutputStream;
import java.io.InputStream;

import dk.brics.automaton.Automaton;
import dk.brics.automaton.RegExp;
import dk.brics.automaton.RunAutomaton;
import dk.brics.automaton.SpecialOperations;

/**
 * Handles storing and retrieving automata to/from files
 * This is highly useful to do as building these deterministic automata from scratch can take many seconds.
 * @author dl387
 *
 */
class AutomatonInitialiser {
	private static final ResourceGetter resourceGetter = new ResourceGetter("uk/ac/cam/ch/wwmm/opsin/resources/serialisedAutomata/");
	/**
	 * In preference serialised automata and their hashes will be looked for in the in your resource folder in your workspace
	 * If it cannot be found there then these files will be looked for in the standard resource folder
	 * (this is actually the standard behaviour of the resourceGetter but I'm reiterating it here as if the stored hash doesn't match
	 * the current hash then the creation of an updated serialised automaton and hash will occur in the workspace resource folder as the standard
	 * resource folder will not typically be writable)
	 * @param automatonName : A name for the automaton so that it can it can be saved/loaded from disk
	 * @param regex : the regex from which to build the RunAutomaton
	 * @param reverseAutomaton : should the automaton be reversed
	 * @param tableize: if true, a transition table is created which makes the run method faster in return of a higher memory usage (adds ~256kb)
	 * @return A RunAutomaton, may have been built from scratch or loaded from a file
	 */
	static RunAutomaton getAutomaton(String automatonName, String regex, boolean tableize, boolean reverseAutomaton) {
		if (reverseAutomaton){
			automatonName+="_reversed_";
		}
		String currentRegexHash =Integer.toString(regex.hashCode());
		try {
			/*
			 * This file indicates the hash used to generate the automaton on the disk
			 * This throws an exception if the file cannot be found 
			 */
			String regexHash = resourceGetter.getString(automatonName + "RegexHash.txt");
			if (regexHash.equals(currentRegexHash)){
				InputStream automatonInput= resourceGetter.getStream(automatonName +"SerialisedAutomaton.txt");
				return RunAutomaton.load(automatonInput);
			}
		}
		catch (Exception e) {
			//automaton could not be loaded either because the regex hash could not be loaded, the hashes did not match or the automaton could not be loaded
		}
		Automaton a = new RegExp(regex).toAutomaton();
		if (reverseAutomaton){
			SpecialOperations.reverse(a);
		}
		RunAutomaton ra = new RunAutomaton(a, tableize);
		try {
			FileOutputStream regexHashFOS =(FileOutputStream) resourceGetter.getOutputStream(automatonName + "RegexHash.txt");
			regexHashFOS.write(currentRegexHash.getBytes());
			ra.store(resourceGetter.getOutputStream(automatonName + "SerialisedAutomaton.txt"));
		} catch (Exception e) {
			System.err.println("WARNING: Could not serialize one of OPSIN's automata to disk. This will not prevent OPSIN from working but is unexpected!");
		}

		return ra;
	}

}
