package uk.ac.cam.ch.wwmm.opsin;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.lang.reflect.Method;
import java.util.Collections;
import java.util.List;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import org.apache.commons.cli.UnrecognizedOptionException;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import uk.ac.cam.ch.wwmm.opsin.OpsinResult.OPSIN_RESULT_STATUS;

import nu.xom.Attribute;
import nu.xom.Element;

/** The "master" class, to turn a name into a structure.
 *
 * @author ptc24
 * @author dl387
 */
public class NameToStructure {
	
	private static final Logger LOG = Logger.getLogger(NameToStructure.class);

	/**Does finite-state non-destructive parsing on chemical names.*/
	private Parser parser;
	
	/**Applies OPSIN's grammar to tokenise and assign meanings tokens.*/
	private ParseRules parseRules;
	
	/**Contains rules on how to interpret suffixes*/
	private SuffixRules suffixRules;

	private static NameToStructure NTS_INSTANCE;

	public static synchronized NameToStructure getInstance() {
		if (NTS_INSTANCE ==null){
			NTS_INSTANCE = new NameToStructure();
		}
		return NTS_INSTANCE;
	}

	/**Initialises the name-to-structure converter.
	 *
	 * @throws NameToStructureException If the converter cannot be initialised, most likely due to bad or missing data files.
	 */
	private NameToStructure() {
		LOG.info("Initialising OPSIN... ");
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
		LOG.info("OPSIN initialised");
	}

	/**
	 * Convenience method for converting a name to CML with OPSIN's default options
	 * @param name The chemical name to parse.
	 * @return A CML element, containing the parsed molecule, or null if the name was uninterpretable.
	 */
	public Element parseToCML(String name) {
		OpsinResult result = parseChemicalName(name);
		Element cml = result.getCml();
		if(cml != null && LOG.isDebugEnabled()){
			LOG.debug(new XOMFormatter().elemToString(result.getCml()));
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
		n2sConfig = n2sConfig.clone();//avoid n2sconfig being modified mid name processing
		if (name==null){
			throw new IllegalArgumentException("String given for name was null");
		}
		String message = "";
		try {
			LOG.debug(name);
			String modifiedName = PreProcessor.preProcess(name);
			List<Element> parses = parser.parse(n2sConfig, modifiedName);
			//if(LOG.isDebugEnabled()) for(Element parse : parses) LOG.debug(new XOMFormatter().elemToString(parse));
			Collections.sort(parses, new SortParses());//fewer tokens preferred
			Fragment fragGeneratedWithWarning = null;
			String warningMessage = null;
			for(Element parse : parses) {
				try {
					if (LOG.isDebugEnabled()){
						LOG.debug(new XOMFormatter().elemToString(parse));
					}
					//Performs XML manipulation e.g. nesting bracketing, processing some nomenclatures
					new ComponentGenerator(n2sConfig).processParse(parse);
					if (LOG.isDebugEnabled()){
						LOG.debug(new XOMFormatter().elemToString(parse));
					}
					BuildState state = new BuildState(n2sConfig);
					//Converts the XML to fragments (handles many different nomenclatueres for describing structure). Assigns locants 
					new ComponentProcessor(suffixRules, state).processParse(parse);
					if (LOG.isDebugEnabled()){
						LOG.debug(new XOMFormatter().elemToString(parse));
					}
					//Constructs a single fragment from the fragments generated by the ComponentProcessor. Applies stereochemistry
					Fragment frag = new StructureBuilder(state).buildFragment(parse);
					if (LOG.isDebugEnabled()){
						LOG.debug(new XOMFormatter().elemToString(parse));
					}
					if (state.getWarningMessage() == null){
						return new OpsinResult(frag, OPSIN_RESULT_STATUS.SUCCESS, "", name);
					}
					if (fragGeneratedWithWarning == null){
						//record first frag that had a warning but try other parses as they may work without a warning
						fragGeneratedWithWarning = frag;
						warningMessage = state.getWarningMessage();
					}
				} catch (Exception e) {
					if (message.length() ==0){
						message = e.getMessage();
					}
					if (LOG.isDebugEnabled()){
						LOG.debug(e.getMessage(),e);
					}
				}
			}
			if (fragGeneratedWithWarning != null){
				return new OpsinResult(fragGeneratedWithWarning, OPSIN_RESULT_STATUS.WARNING, warningMessage, name);
			}
		} catch (Exception e) {
			message += e.getMessage();
			if(LOG.isDebugEnabled()) {
				LOG.debug(e.getMessage(),e);
			}
		}
		return new OpsinResult(null, OPSIN_RESULT_STATUS.FAILURE, message, name);
	}
	
	/**
	 * Returns an OPSIN parser
	 * This can be used to determine whether a word can be interpreted as being part of a chemical name.
	 * Just because a word can be split into tokens does not mean the word constitutes a valid chemical name
	 * e.g. ester is interpretable but is not in itself a chemical name
	 * @return

	 */
	public static ParseRules getOpsinParser() {
		NameToStructure n2s = NameToStructure.getInstance();
		return n2s.parseRules;
	}

	/**Run OPSIN as a command-line application.
	 *
	 * @param args
	 * @throws Exception
	 */
	public static void main(String [] args) throws Exception {
		Options options = buildCommandLineOptions();
		CommandLineParser parser = new PosixParser();
		CommandLine cmd = null;
		try{
			cmd = parser.parse(options, args);
		}
		catch (UnrecognizedOptionException e) {
			System.err.println(e.getMessage());
			System.exit(1);
		}
		if (cmd.hasOption("h")){
			displayUsage(options);
		}
		if(cmd.hasOption("v")){
			Logger.getLogger("uk.ac.cam.ch.wwmm.opsin").setLevel(Level.DEBUG);
		}

		NameToStructureConfig n2sconfig = generateOpsinConfigObjectFromCmd(cmd);
		
		InputStream input;
		OutputStream output;
		String[] unparsedArgs = cmd.getArgs();
		if (unparsedArgs.length == 0){
			input = System.in;
			output = System.out;
		}
		else if (unparsedArgs.length == 1){
			input = new FileInputStream(new File(unparsedArgs[0]));
			output = System.out;
		}
		else if (unparsedArgs.length == 2){
			input = new FileInputStream(new File(unparsedArgs[0]));
			output = new FileOutputStream(new File(unparsedArgs[1]));
		}
		else {
			input = null;
			output = null;
			displayUsage(options);
		}

		System.err.println("Run the jar using the -h flag for help. Enter a chemical name to begin:");
		String outputType = cmd.getOptionValue("o", "cml");
		if (outputType.equalsIgnoreCase("cml")){
			interactiveCmlOutput(input, output, n2sconfig);
		}
		else if (outputType.equalsIgnoreCase("smi") || outputType.equalsIgnoreCase("smiles")){
			interactiveSmilesOutput(input, output, n2sconfig);
		}
		else if (outputType.equalsIgnoreCase("inchi")){
			interactiveInchiOutput(input, output, n2sconfig, false);
		}
		else if (outputType.equalsIgnoreCase("stdinchi")){
			interactiveInchiOutput(input, output, n2sconfig, true);
		}
		else{
			System.err.println("Unrecognised output format: " + outputType);
			System.err.println("Expected output types are \"cml\", \"smi\", \"inchi\" and \"stdinchi\"");
			System.exit(1);
		}
		if (unparsedArgs.length == 1){
			input.close();
		}
		else if (unparsedArgs.length == 2){
			input.close();
			output.close();
		}
	}

	private static void displayUsage(Options options) {
		HelpFormatter formatter = new HelpFormatter();
		formatter.printHelp("java -jar opsin-[version]-jar-with-dependencies.jar [options] [inputfile] [outputfile]\n" +
				"OPSIN converts systematic chemical names to CML, SMILES or InChI/StdInChI\n" +
				"Names should be new line delimited and may be read from stdin (default) or a file and output to stdout (default) or a file", options);
		System.exit(0);
	}

	private static Options buildCommandLineOptions() throws ParseException {
		Options options = new Options();
		OptionBuilder.withArgName("o");
		OptionBuilder.withLongOpt("output");
		OptionBuilder.hasArg();
		OptionBuilder.withDescription("Sets OPSIN's output format (default cml)\n" +
				"Allowed values are:\n" +
				"cml for Chemical Markup Language\n" +
				"smi for SMILES\n" +
				"inchi for InChI\n" +
				"stdinchi for StdInChI");
		options.addOption(OptionBuilder.create("o"));
		options.addOption("h", "help", false, "Displays the allowed command line flags");
		options.addOption("v", "verbose", false, "Enables debugging");
		
		options.addOption("a", "allowAcidsWithoutAcid", false, "Allows interpretation of acids without the word acid e.g. \"acetic\"");
		options.addOption("f", "detailedFailureAnalysis", false, "Enables reverse parsing to more accurately determine why parsing failed");
		options.addOption("r", "allowRadicals", false, "Enables interpretation of radicals");
		options.addOption("s", "allowUninterpretableStereo", false, "Allows stereochemistry uninterpretable by OPSIN to be ignored");
		options.addOption("w", "wildcardRadicals", false, "Radicals are output as wildcard atoms");
		return options;
	}
	
	/**
	 * Uses the command line parameters to configure a new NameToStructureConfig
	 * @param cmd
	 * @return
	 */
	private static NameToStructureConfig generateOpsinConfigObjectFromCmd(CommandLine cmd) {
		NameToStructureConfig n2sconfig = new NameToStructureConfig();
		n2sconfig.setInterpretAcidsWithoutTheWordAcid(cmd.hasOption("a"));
		n2sconfig.setDetailedFailureAnalysis(cmd.hasOption("f"));
		n2sconfig.setAllowRadicals(cmd.hasOption("r"));
		n2sconfig.setWarnRatherThanFailOnUninterpretableStereochemistry(cmd.hasOption("s"));
		n2sconfig.setOutputRadicalsAsWildCardAtoms(cmd.hasOption("w"));
		return n2sconfig;
	}

	private static void interactiveCmlOutput(InputStream input, OutputStream out, NameToStructureConfig n2sconfig) throws IOException {
		NameToStructure nts = NameToStructure.getInstance();
		BufferedReader inputReader = new BufferedReader(new InputStreamReader(input, "UTF-8"));
		StreamSerializer serializer = new StreamSerializer(out);
		serializer.setIndent(2);
		serializer.writeXMLDeclaration();
		Element cml = new Element("cml", XmlDeclarations.CML_NAMESPACE);
		cml.addAttribute(new Attribute("convention","conventions:molecular"));
		cml.addNamespaceDeclaration("conventions", "http://www.xml-cml.org/convention/");
		cml.addNamespaceDeclaration("cmlDict", "http://www.xml-cml.org/dictionary/cml/");
		cml.addNamespaceDeclaration("nameDict", "http://www.xml-cml.org/dictionary/cml/name/");
		serializer.writeStartTag(cml);
		int id =1;
		String name;
		while((name =inputReader.readLine()) != null) {
			OpsinResult result = nts.parseChemicalName(name, n2sconfig);
			Element output = result.getCml();
			if(output == null) {
				System.err.println(result.getMessage());
				Element uninterpretableMolecule = new Element("molecule", XmlDeclarations.CML_NAMESPACE);
				uninterpretableMolecule.addAttribute(new Attribute("id", "m" + id++));
				Element nameEl = new Element("name", XmlDeclarations.CML_NAMESPACE);
				nameEl.appendChild(name);
				nameEl.addAttribute(new Attribute("dictRef", "nameDict:unknown"));
				uninterpretableMolecule.appendChild(nameEl);
				serializer.write(uninterpretableMolecule);
				serializer.flush();
			} else {
				Element molecule = XOMTools.getChildElementsWithTagName(output, "molecule").get(0);
				molecule.getAttribute("id").setValue("m" + id++);
				serializer.write(molecule);
				serializer.flush();
			}
		}
		serializer.writeEndTag(cml);
		serializer.flush();
	}
	
	private static void interactiveSmilesOutput(InputStream input, OutputStream out, NameToStructureConfig n2sconfig) throws IOException {
		NameToStructure nts = NameToStructure.getInstance();
		BufferedReader inputReader = new BufferedReader(new InputStreamReader(input, "UTF-8"));
		BufferedWriter outputWriter = new BufferedWriter(new OutputStreamWriter(out, "UTF-8"));
		String name;
		while((name =inputReader.readLine()) != null) {
			OpsinResult result = nts.parseChemicalName(name, n2sconfig);
			String output = result.getSmiles();
			if(output == null) {
				System.err.println(result.getMessage());
			} else {
				outputWriter.write(output);
			}
			outputWriter.newLine();
			outputWriter.flush();
		}
	}

	private static void interactiveInchiOutput(InputStream input, OutputStream out, NameToStructureConfig n2sconfig, boolean produceStdInChI) throws Exception {
		NameToStructure nts = NameToStructure.getInstance();
		BufferedReader inputReader = new BufferedReader(new InputStreamReader(input, "UTF-8"));
		BufferedWriter outputWriter = new BufferedWriter(new OutputStreamWriter(out, "UTF-8"));
		Class<?> c;
		try {
			c = Class.forName("uk.ac.cam.ch.wwmm.opsin.NameToInchi");
		} catch (ClassNotFoundException e) {
			System.err.println("Could not initialise NameToInChI module. Is it on your classpath?");
			throw new RuntimeException(e);
		}
		Method m;
		if (produceStdInChI){
			m = c.getMethod("convertResultToStdInChI", new Class[]{OpsinResult.class});
		}
		else{
			m = c.getMethod("convertResultToInChI", new Class[]{OpsinResult.class});
		}

		String name;
		while((name =inputReader.readLine()) != null) {
			OpsinResult result = nts.parseChemicalName(name, n2sconfig);
			String output = (String) m.invoke(null, new Object[]{result});
			if(output == null) {
				System.err.println(result.getMessage());
			} else {
				outputWriter.write(output);
			}
			outputWriter.newLine();
			outputWriter.flush();
		}
	}
}
