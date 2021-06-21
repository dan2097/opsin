package uk.ac.cam.ch.wwmm.opsin;

import java.io.BufferedInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.nio.charset.StandardCharsets;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import dk.brics.automaton.Automaton;
import dk.brics.automaton.RegExp;
import dk.brics.automaton.RunAutomaton;
import dk.brics.automaton.SpecialOperations;

/**
 * Handles storing and retrieving automata to/from files
 * This is highly useful to do as building these deterministic automata from scratch can take minutes
 * @author dl387
 *
 */
class AutomatonInitialiser {
	
	private static final Logger LOG = LogManager.getLogger(AutomatonInitialiser.class);
	private final ResourceGetter resourceGetter;
	
	AutomatonInitialiser(String resourcePath) {
		resourceGetter = new ResourceGetter(resourcePath);
	}

	/**
	 * In preference serialised automata and their hashes will be looked for in the resource folder in your working directory
	 * If it cannot be found there then these files will be looked for in the standard resource folder
	 * (this is actually the standard behaviour of the resourceGetter but I'm reiterating it here as if the stored hash doesn't match
	 * the current hash then the creation of an updated serialised automaton and hash will occur in the working directory resource folder as the standard
	 * resource folder will not typically be writable)
	 * @param automatonName : A name for the automaton so that it can it can be saved/loaded from disk
	 * @param regex : the regex from which to build the RunAutomaton
	 * @param reverseAutomaton : should the automaton be reversed
	 * @param tableize: if true, a transition table is created which makes the run method faster in return of a higher memory usage (adds ~256kb)
	 * @return A RunAutomaton, may have been built from scratch or loaded from a file
	 */
	RunAutomaton loadAutomaton(String automatonName, String regex, boolean tableize, boolean reverseAutomaton) {
		if (reverseAutomaton){
			automatonName+="_reversed_";
		}
		try{
			if (isAutomatonCached(automatonName, regex)) {
				return loadCachedAutomaton(automatonName);
			}
		}
		catch (IOException e) {
			LOG.warn("Error loading cached automaton: "+automatonName, e);
		}
		RunAutomaton automaton = createAutomaton(regex, tableize, reverseAutomaton);
		cacheAutomaton(automatonName, automaton, regex);
		return automaton;
	}
	
	private boolean isAutomatonCached(String automatonName, String regex) {
		String currentRegexHash = getRegexHash(regex);
		String cachedRegexHash = getCachedRegexHash(automatonName);
		return currentRegexHash.equals(cachedRegexHash);
	}
	
	private String getRegexHash(String regex) {
		return Integer.toString(regex.hashCode());
	}

	private String getCachedRegexHash(String automatonName) {
		/*This file contains the hashcode of the regex which was used to generate the automaton on the disk */
		return resourceGetter.getFileContentsAsString(automatonName + "RegexHash.txt");
	}
	
	private RunAutomaton loadCachedAutomaton(String automatonName) throws IOException{
		try (InputStream automatonInput = resourceGetter.getInputstreamFromFileName(automatonName +"SerialisedAutomaton.aut")){
			return RunAutomaton.load(new BufferedInputStream(automatonInput));
		} catch (Exception e) {
			IOException ioe = new IOException("Error loading automaton");
			ioe.initCause(e);
			throw ioe;
		}
	}
	
	private static RunAutomaton createAutomaton(String regex, boolean tableize, boolean reverseAutomaton) {
		Automaton a = new RegExp(regex).toAutomaton();
		if (reverseAutomaton){
			SpecialOperations.reverse(a);
		}
		return new RunAutomaton(a, tableize);
	}

	private void cacheAutomaton(String automatonName, RunAutomaton automaton, String regex) {
		try (OutputStream regexHashOutputStream = resourceGetter.getOutputStream(automatonName + "RegexHash.txt")) {
			regexHashOutputStream.write(getRegexHash(regex).getBytes(StandardCharsets.UTF_8));
			try (OutputStream automatonOutputStream = resourceGetter.getOutputStream(automatonName + "SerialisedAutomaton.aut")) {
				automaton.store(automatonOutputStream);
			}
		} catch (IOException e) {
			LOG.warn("Error serialising automaton: "+automatonName, e);
		}
	}
	
}
