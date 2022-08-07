package uk.ac.cam.ch.wwmm.opsin;

import java.io.InputStream;
import java.util.Collections;
import java.util.List;
import java.util.Properties;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import uk.ac.cam.ch.wwmm.opsin.OpsinResult.OPSIN_RESULT_STATUS;

/** The "master" class, to turn a name into a structure.
 *
 * @author ptc24
 * @author dl387
 */
public class NameToStructure {
	
	private static final Logger LOG = LogManager.getLogger(NameToStructure.class);
	
	/**Applies OPSIN's grammar to tokenise and assign meaning to tokens*/
	private ParseRules parseRules;

	/**Parses a chemical name into one (or more in the case of ambiguity) parse trees*/
	private Parser parser;
	
	/**Which suffixes apply to what and what their effects are*/
	private SuffixRules suffixRules;

	private static NameToStructure NTS_INSTANCE;

	public static synchronized NameToStructure getInstance() {
		if (NTS_INSTANCE == null) {
			NTS_INSTANCE = new NameToStructure();
		}
		return NTS_INSTANCE;
	}
	
	/**
	 * Returns the version of the OPSIN library
	 * @return Version number String
	 */
	public static String getVersion() {
		try(InputStream is = NameToStructure.class.getResourceAsStream("opsinbuild.props")) {
			Properties props = new Properties();
			props.load(is);
			return props.getProperty("version");
		}
		catch (Exception e) {
			return null;
		}
	}

	/**Initialises the name-to-structure converter.
	 *
	 * @throws NameToStructureException If the converter cannot be initialised, most likely due to bad or missing data files.
	 */
	private NameToStructure() {
		LOG.debug("Initialising OPSIN... ");
		try {
			/*Initialise all of OPSIN's classes. Some classes are injected as dependencies into subsequent classes*/

			//Allows retrieving of OPSIN resources
			ResourceGetter resourceGetter = new ResourceGetter("uk/ac/cam/ch/wwmm/opsin/resources/");
			ResourceManager resourceManager = new ResourceManager(resourceGetter);
			WordRules wordRules = new WordRules(resourceGetter);
			parseRules = new ParseRules(resourceManager);
			Tokeniser tokeniser = new Tokeniser(parseRules);
			parser = new Parser(wordRules, tokeniser, resourceManager);
			suffixRules = new SuffixRules(resourceGetter);
		} catch (Exception e) {
			throw new NameToStructureException(e.getMessage(), e);
		}
		LOG.debug("OPSIN initialised");
	}

	/**
	 * Convenience method for converting a name to CML with OPSIN's default options
	 * @param name The chemical name to parse.
	 * @return A CML element, containing the parsed molecule, or null if the name was uninterpretable.
	 */
	public String parseToCML(String name) {
		OpsinResult result = parseChemicalName(name);
		String cml = result.getCml();
		if(cml != null && LOG.isDebugEnabled()){
			LOG.debug(cml);
		}
		return cml;
	}

	/**
	 * Convenience method for converting a name to SMILES with OPSIN's default options
	 * @param name The chemical name to parse.
	 * @return A SMILES string describing the parsed molecule, or null if the name was uninterpretable.
	 */
	public String parseToSmiles(String name) {
		OpsinResult result = parseChemicalName(name);
		String smiles = result.getSmiles();
		LOG.debug(smiles);
		return smiles;
	}

	/**Parses a chemical name, returning an OpsinResult which represents the molecule.
	 * This object contains in the status whether the name was parsed successfully
	 * A message which may contain additional information if the status was warning/failure
	 * The OpsinResult has methods to generate a SMILES or CML representation
	 * For InChI, the OpsinResult should be given to the NameToInchi class
	 *
	 * @param name The chemical name to parse.
	 * @return OpsinResult
	 */
	public OpsinResult parseChemicalName(String name) {
		NameToStructureConfig n2sConfig = NameToStructureConfig.getDefaultConfigInstance();
		return parseChemicalName(name, n2sConfig);
	}

	/**Parses a chemical name, returning an OpsinResult which represents the molecule.
	 * This object contains in the status whether the name was parsed successfully
	 * A message which may contain additional information if the status was warning/failure
	 * CML and SMILES representations may be retrieved directly from the object
	 * InChI may be generate using NameToInchi
	 *
	 * @param name The chemical name to parse.
	 * @param n2sConfig Options to control how OPSIN interprets the name.
	 * @return OpsinResult
	 */
	public OpsinResult parseChemicalName(String name, NameToStructureConfig n2sConfig) {
		if (name == null){
			throw new IllegalArgumentException("String given for name was null");
		}
		n2sConfig = n2sConfig.clone();//avoid n2sconfig being modified mid name processing

		List<Element> parses;
		try {
			LOG.debug(name);
			String modifiedName = PreProcessor.preProcess(name);
			parses = parser.parse(n2sConfig, modifiedName);
			Collections.sort(parses, new SortParses());//fewer tokens preferred
		} catch (Exception e) {
			if(LOG.isDebugEnabled()) {
				LOG.debug(e.getMessage(), e);
			}
			String message = e.getMessage() != null ? e.getMessage() : "exception with null message";
			return new OpsinResult(null, OPSIN_RESULT_STATUS.FAILURE, message, name);
		}
		String reasonForFailure = "";
		Fragment fragGeneratedWithWarning = null;
		List<OpsinWarning> warnings = Collections.emptyList();
		for(Element parse : parses) {
			try {
				if (LOG.isDebugEnabled()) {
					LOG.debug(parse.toXML());
				}
				//Performs XML manipulation e.g. nesting bracketing, processing some nomenclatures
				BuildState state = new BuildState(n2sConfig);
				new ComponentGenerator(state).processParse(parse);
				if (LOG.isDebugEnabled()) {
					LOG.debug(parse.toXML());
				}
				//Converts the XML to fragments (handles many different nomenclatueres for describing structure). Assigns locants
				new ComponentProcessor(state, new SuffixApplier(state, suffixRules)).processParse(parse);
				if (LOG.isDebugEnabled()) {
					LOG.debug(parse.toXML());
				}
				//Constructs a single fragment from the fragments generated by the ComponentProcessor. Applies stereochemistry
				Fragment frag = new StructureBuilder(state).buildFragment(parse);
				if (LOG.isDebugEnabled()) {
					LOG.debug(parse.toXML());
				}
				if (state.getWarnings().isEmpty()) {
					return new OpsinResult(frag, OPSIN_RESULT_STATUS.SUCCESS, "", name);
				}
				if (fragGeneratedWithWarning == null) {
					//record first frag that had a warning but try other parses as they may work without a warning
					fragGeneratedWithWarning = frag;
					warnings = state.getWarnings();
				}
			} catch (Exception e) {
				if (reasonForFailure.length() == 0) {
					reasonForFailure = e.getMessage() != null ? e.getMessage() : "exception with null message";
				}
				if (LOG.isDebugEnabled()) {
					LOG.debug(e.getMessage(), e);
				}
			}
		}
		if (fragGeneratedWithWarning != null) {
			return new OpsinResult(fragGeneratedWithWarning, OPSIN_RESULT_STATUS.WARNING, warnings, name);
		}
		return new OpsinResult(null, OPSIN_RESULT_STATUS.FAILURE, reasonForFailure, name);
	}
	
	/**
	 * Returns an OPSIN parser
	 * This can be used to determine whether a word can be interpreted as being part of a chemical name.
	 * Just because a word can be split into tokens does not mean the word constitutes a valid chemical name
	 * e.g. ester is interpretable but is not in itself a chemical name
	 * @return Opsin parser for recognition/parsing of a chemical word
	 */
	public static ParseRules getOpsinParser() {
		NameToStructure n2s = NameToStructure.getInstance();
		return n2s.parseRules;
	}

}
