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
import java.nio.charset.StandardCharsets;
import java.util.ArrayDeque;
import java.util.Deque;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Option.Builder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.UnrecognizedOptionException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;

import com.ctc.wstx.api.WstxOutputProperties;
import com.ctc.wstx.stax.WstxOutputFactory;

public class Cli {

	private enum InchiType {
		inchiWithFixedH, stdInchi, stdInchiKey
	}

	/**
	 * Run OPSIN as a command-line application.
	 *
	 * @param args
	 * @throws Exception
	 */
	public static void main(String[] args) throws Exception {
		Options options = buildCommandLineOptions();
		CommandLineParser parser = new DefaultParser();
		CommandLine cmd = null;
		try {
			cmd = parser.parse(options, args);
		} catch (UnrecognizedOptionException e) {
			System.err.println(e.getMessage());
			System.exit(1);
		}
		if (cmd.hasOption("h")) {
			displayUsage(options);
		}
		if (cmd.hasOption("v")) {
			Configurator.setLevel("uk.ac.cam.ch.wwmm.opsin", Level.DEBUG);
		}

		NameToStructureConfig n2sconfig = generateOpsinConfigObjectFromCmd(cmd);

		InputStream input = System.in;
		OutputStream output = System.out;
		String[] unparsedArgs = cmd.getArgs();
		if (unparsedArgs.length == 0) {
			System.err.println("Run the jar using the -h flag for help. Enter a chemical name to begin:");
		} else if (unparsedArgs.length == 1) {
			input = new FileInputStream(new File(unparsedArgs[0]));
		} else if (unparsedArgs.length == 2) {
			input = new FileInputStream(new File(unparsedArgs[0]));
			output = new FileOutputStream(new File(unparsedArgs[1]));
		} else {
			displayUsage(options);
		}
		try {
			String outputType = cmd.getOptionValue("o", "smi");
			boolean outputName = cmd.hasOption("n");
			int threads = parseThreadCount(cmd);
			if (outputType.equalsIgnoreCase("cml")) {
				interactiveCmlOutput(input, output, n2sconfig);
			} else if (outputType.equalsIgnoreCase("smi") || outputType.equalsIgnoreCase("smiles")) {
				interactiveSmilesOutput(input, output, n2sconfig, false, outputName, threads);
			} else if (outputType.equalsIgnoreCase("inchi")) {
				interactiveInchiOutput(input, output, n2sconfig, InchiType.inchiWithFixedH, outputName, threads);
			} else if (outputType.equalsIgnoreCase("stdinchi")) {
				interactiveInchiOutput(input, output, n2sconfig, InchiType.stdInchi, outputName, threads);
			} else if (outputType.equalsIgnoreCase("stdinchikey")) {
				interactiveInchiOutput(input, output, n2sconfig, InchiType.stdInchiKey, outputName, threads);
			} else if (outputType.equalsIgnoreCase("extendedsmi") || outputType.equalsIgnoreCase("extendedsmiles")
					|| outputType.equalsIgnoreCase("cxsmi") || outputType.equalsIgnoreCase("cxsmiles")) {
				interactiveSmilesOutput(input, output, n2sconfig, true, outputName, threads);
			} else {
				System.err.println("Unrecognised output format: " + outputType);
				System.err.println(
						"Expected output types are \"cml\", \"smi\", \"inchi\", \"stdinchi\" and \"stdinchikey\"");
				System.exit(1);
			}
		} finally {
			if (output != System.out) {
				output.close();
			}
			if (input != System.in) {
				input.close();
			}
		}
	}

	private static void displayUsage(Options options) {
		HelpFormatter formatter = new HelpFormatter();
		String version = NameToStructure.getVersion();
		formatter.printHelp("java -jar opsin-" + (version != null ? version : "[version]")
				+ "-jar-with-dependencies.jar [options] [inputfile] [outputfile]" + OpsinTools.NEWLINE
				+ "OPSIN converts systematic chemical names to CML, SMILES or InChI/StdInChI/StdInChIKey"
				+ OpsinTools.NEWLINE
				+ "Names should be new line delimited and may be read from stdin (default) or a file and output to stdout (default) or a file",
				options);
		System.exit(0);
	}

	private static Options buildCommandLineOptions() {
		Options options = new Options();
		Builder outputBuilder = Option.builder("o");
		outputBuilder.longOpt("output");
		outputBuilder.hasArg();
		outputBuilder.argName("format");
		StringBuilder outputOptionsDesc = new StringBuilder();
		outputOptionsDesc.append("Sets OPSIN's output format (default smi)").append(OpsinTools.NEWLINE);
		outputOptionsDesc.append("Allowed values are:").append(OpsinTools.NEWLINE);
		outputOptionsDesc.append("cml for Chemical Markup Language").append(OpsinTools.NEWLINE);
		outputOptionsDesc.append("smi for SMILES").append(OpsinTools.NEWLINE);
		outputOptionsDesc.append("extendedsmi for Extended SMILES").append(OpsinTools.NEWLINE);
		outputOptionsDesc.append("inchi for InChI (with FixedH)").append(OpsinTools.NEWLINE);
		outputOptionsDesc.append("stdinchi for StdInChI").append(OpsinTools.NEWLINE);
		outputOptionsDesc.append("stdinchikey for StdInChIKey");
		outputBuilder.desc(outputOptionsDesc.toString());
		options.addOption(outputBuilder.build());
		options.addOption("h", "help", false, "Displays the allowed command line flags");
		options.addOption("v", "verbose", false, "Enables debugging");

		options.addOption("a", "allowAcidsWithoutAcid", false,
				"Allows interpretation of acids without the word acid e.g. \"acetic\"");
		options.addOption("f", "detailedFailureAnalysis", false,
				"Enables reverse parsing to more accurately determine why parsing failed");
		options.addOption("n", "name", false, "Include name in SMILES/InChI output (tab delimited)");
		options.addOption("r", "allowRadicals", false, "Enables interpretation of radicals");
		options.addOption("s", "allowUninterpretableStereo", false,
				"Allows stereochemistry uninterpretable by OPSIN to be ignored");
		options.addOption("w", "wildcardRadicals", false, "Radicals are output as wildcard atoms");
		Builder threadsBuilder = Option.builder("t");
		threadsBuilder.longOpt("threads");
		threadsBuilder.hasArg();
		threadsBuilder.argName("count");
		threadsBuilder.desc("Number of names to interpret concurrently (default 1)." + OpsinTools.NEWLINE
				+ "Output order always matches input order. Use 0 for one thread per available processor.");
		options.addOption(threadsBuilder.build());
		return options;
	}

	/**
	 * Uses the command line parameters to configure a new NameToStructureConfig
	 * 
	 * @param cmd
	 * @return The configured NameToStructureConfig
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

	private static void interactiveCmlOutput(InputStream input, OutputStream out, NameToStructureConfig n2sconfig) throws IOException, XMLStreamException {
		NameToStructure nts = NameToStructure.getInstance();
		BufferedReader inputReader = new BufferedReader(new InputStreamReader(input, StandardCharsets.UTF_8));
		XMLOutputFactory factory = new WstxOutputFactory();
		factory.setProperty(WstxOutputProperties.P_OUTPUT_ESCAPE_CR, false);
		XMLStreamWriter writer = factory.createXMLStreamWriter(out, "UTF-8");
		writer = new IndentingXMLStreamWriter(writer, 2);
		writer.writeStartDocument();
		CMLWriter cmlWriter = new CMLWriter(writer);
		cmlWriter.writeCmlStart();
		int id = 1;
		String line;
		while ((line = inputReader.readLine()) != null) {
			int splitPoint = line.indexOf('\t');
			String name = splitPoint >= 0 ? line.substring(0, splitPoint) : line;
			OpsinResult result = nts.parseChemicalName(name, n2sconfig);
			String message = result.getMessage();
			if (!message.isEmpty()) {
				System.err.println(message);
			}
			Fragment structure = result.getStructure();
			cmlWriter.writeMolecule(structure, name, id++);
			writer.flush();
		}
		cmlWriter.writeCmlEnd();
		writer.writeEndDocument();
		writer.flush();
		writer.close();
	}


	/** What OPSIN produced for one name: the rendered output, or null if it could not be
	 * interpreted, together with any message OPSIN emitted while interpreting it. */
	private static final class NameResult {
		private final String output;
		private final String message;

		NameResult(String output, String message) {
			this.output = output;
			this.message = message;
		}
	}

	private interface NameInterpreter {
		NameResult interpret(String name);
	}

	private static int parseThreadCount(CommandLine cmd) {
		String value = cmd.getOptionValue("t");
		if (value == null) {
			return 1;
		}
		int threads;
		try {
			threads = Integer.parseInt(value.trim());
		} catch (NumberFormatException e) {
			System.err.println("Number of threads must be an integer: " + value);
			System.exit(1);
			return 1;
		}
		if (threads < 0) {
			System.err.println("Number of threads may not be negative: " + value);
			System.exit(1);
		}
		return threads == 0 ? Runtime.getRuntime().availableProcessors() : threads;
	}

	private static String nameFromLine(String line) {
		int splitPoint = line.indexOf('\t');
		return splitPoint >= 0 ? line.substring(0, splitPoint) : line;
	}

	private static void writeResult(BufferedWriter outputWriter, String line, boolean outputName,
			NameResult result) throws IOException {
		if (!result.message.isEmpty()) {
			System.err.println(result.message);
		}
		if (result.output != null) {
			outputWriter.write(result.output);
		}
		if (outputName) {
			outputWriter.write('\t');
			outputWriter.write(line);
		}
		outputWriter.newLine();
	}

	/**
	 * Interprets each line of input and writes the results, always in input order.
	 *
	 * With one thread this is the historical behaviour: interpret a name, write it, flush, so
	 * that the jar stays usable interactively. With more than one, names are interpreted
	 * concurrently but a bounded window of pending results is drained in order, which keeps
	 * memory flat no matter how large the input is. NameToStructure is shared across the pool;
	 * parseChemicalName holds no mutable state between calls and clones the config it is given.
	 */
	private static void streamNames(BufferedReader inputReader, BufferedWriter outputWriter,
			boolean outputName, int threads, NameInterpreter interpreter) throws IOException {
		if (threads <= 1) {
			String line;
			while ((line = inputReader.readLine()) != null) {
				writeResult(outputWriter, line, outputName, interpreter.interpret(nameFromLine(line)));
				outputWriter.flush();
			}
			return;
		}
		ExecutorService pool = Executors.newFixedThreadPool(threads);
		int window = threads * 4;
		Deque<String> pendingLines = new ArrayDeque<>(window);
		Deque<Future<NameResult>> pendingResults = new ArrayDeque<>(window);
		try {
			String line;
			while ((line = inputReader.readLine()) != null) {
				final String name = nameFromLine(line);
				pendingLines.addLast(line);
				pendingResults.addLast(pool.submit(new Callable<NameResult>() {
					public NameResult call() {
						return interpreter.interpret(name);
					}
				}));
				if (pendingResults.size() >= window) {
					drainOne(outputWriter, pendingLines, pendingResults, outputName);
				}
			}
			while (!pendingResults.isEmpty()) {
				drainOne(outputWriter, pendingLines, pendingResults, outputName);
			}
			outputWriter.flush();
		} finally {
			pool.shutdownNow();
		}
	}

	private static void drainOne(BufferedWriter outputWriter, Deque<String> pendingLines,
			Deque<Future<NameResult>> pendingResults, boolean outputName) throws IOException {
		String line = pendingLines.removeFirst();
		Future<NameResult> future = pendingResults.removeFirst();
		NameResult result;
		try {
			result = future.get();
		} catch (InterruptedException e) {
			Thread.currentThread().interrupt();
			throw new IOException("Interrupted while interpreting: " + line, e);
		} catch (ExecutionException e) {
			throw new IOException("Failed to interpret: " + line, e.getCause());
		}
		writeResult(outputWriter, line, outputName, result);
	}

	private static void interactiveSmilesOutput(InputStream input, OutputStream out, NameToStructureConfig n2sconfig, boolean extendedSmiles, boolean outputName, int threads) throws IOException {
		final NameToStructure nts = NameToStructure.getInstance();
		final NameToStructureConfig config = n2sconfig;
		final boolean extended = extendedSmiles;
		BufferedReader inputReader = new BufferedReader(new InputStreamReader(input, StandardCharsets.UTF_8));
		BufferedWriter outputWriter = new BufferedWriter(new OutputStreamWriter(out, StandardCharsets.UTF_8));
		streamNames(inputReader, outputWriter, outputName, threads, new NameInterpreter() {
			public NameResult interpret(String name) {
				OpsinResult result = nts.parseChemicalName(name, config);
				return new NameResult(extended ? result.getExtendedSmiles() : result.getSmiles(), result.getMessage());
			}
		});
	}

	private static void interactiveInchiOutput(InputStream input, OutputStream out, NameToStructureConfig n2sconfig, InchiType inchiType, boolean outputName, int threads) throws Exception {
		final NameToStructure nts = NameToStructure.getInstance();
		final NameToStructureConfig config = n2sconfig;
		final InchiType type = inchiType;
		BufferedReader inputReader = new BufferedReader(new InputStreamReader(input, StandardCharsets.UTF_8));
		BufferedWriter outputWriter = new BufferedWriter(new OutputStreamWriter(out, StandardCharsets.UTF_8));
		streamNames(inputReader, outputWriter, outputName, threads, new NameInterpreter() {
			public NameResult interpret(String name) {
				OpsinResult result = nts.parseChemicalName(name, config);
				String output;
				switch (type) {
				case inchiWithFixedH:
					output = NameToInchi.convertResultToInChI(result);
					break;
				case stdInchi:
					output = NameToInchi.convertResultToStdInChI(result);
					break;
				case stdInchiKey:
					output = NameToInchi.convertResultToStdInChIKey(result);
					break;
				default:
					throw new IllegalArgumentException("Unexepected enum value: " + type);
				}
				return new NameResult(output, result.getMessage());
			}
		});
	}
}
