package uk.ac.cam.ch.wwmm.opsin;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;

import org.junit.jupiter.api.Test;

/**
 * Concurrent interpretation must not change what the CLI emits: same names in, same bytes out,
 * in the same order.
 */
public class CliThreadingTest {

	private static final String[] NAMES = {
			"ethanol", "benzene", "not a chemical name at all", "toluene",
			"2,4,6-trinitrotoluene", "", "cyclohexanone", "alpha-pinene",
			"(2R,3S)-2,3-dihydroxybutanedioic acid", "caffeine", "zzzzzz",
			"1,3,7-trimethylpurine-2,6-dione", "glucose", "L-alanyl-L-glutamine",
			"bicyclo[2.2.1]heptane", "phenol", "acetic acid", "naphthalene",
			"pyridine", "furan", "thiophene", "indole", "quinoline", "anthracene",
	};

	private String runCli(String[] args, String input) throws Exception {
		ByteArrayOutputStream out = new ByteArrayOutputStream();
		PrintStream originalErr = System.err;
		PrintStream originalOut = System.out;
		try {
			System.setErr(new PrintStream(new ByteArrayOutputStream(), true, "UTF-8"));
			System.setOut(new PrintStream(out, true, "UTF-8"));
			System.setIn(new ByteArrayInputStream(input.getBytes(StandardCharsets.UTF_8)));
			Cli.main(args);
		} finally {
			System.setErr(originalErr);
			System.setOut(originalOut);
		}
		return out.toString("UTF-8");
	}

	private String input() {
		StringBuilder sb = new StringBuilder();
		//repeated so that the bounded window is exercised several times over
		for (int i = 0; i < 20; i++) {
			for (String name : NAMES) {
				sb.append(name).append('\n');
			}
		}
		return sb.toString();
	}

	@Test
	public void threadedSmilesOutputMatchesSequential() throws Exception {
		String in = input();
		String sequential = runCli(new String[] { "-o", "smi", "-n" }, in);
		for (String threads : new String[] { "1", "2", "3", "8" }) {
			assertEquals(sequential, runCli(new String[] { "-o", "smi", "-n", "-t", threads }, in),
					"output with -t " + threads + " should be identical to sequential output");
		}
	}

	@Test
	public void threadedOutputPreservesInputOrderWithoutNameEcho() throws Exception {
		String in = input();
		String sequential = runCli(new String[] { "-o", "smi" }, in);
		assertEquals(sequential, runCli(new String[] { "-o", "smi", "-t", "6" }, in),
				"blank lines for uninterpretable names must stay in position");
	}

	@Test
	public void threadCountOfZeroUsesAvailableProcessors() throws Exception {
		String in = input();
		assertEquals(runCli(new String[] { "-o", "smi", "-n" }, in),
				runCli(new String[] { "-o", "smi", "-n", "-t", "0" }, in));
	}
}
