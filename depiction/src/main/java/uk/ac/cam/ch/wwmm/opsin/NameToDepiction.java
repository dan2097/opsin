package uk.ac.cam.ch.wwmm.opsin;

import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.awt.image.RenderedImage;
import java.io.BufferedReader;
import java.io.InputStreamReader;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.SwingUtilities;

import sea36.chem.core.CMLMolecule;
import sea36.chem.graphics.layout.Layout2D;
import sea36.chem.graphics.renderer.Renderer;
import sea36.chem.tools.HydrogenCalculator;

public class NameToDepiction {

	private NameToStructure n2s;
	public NameToDepiction() throws Exception {
		n2s = NameToStructure.getInstance();
	}

	/**Parses a chemical name, returning a RenderedImage depiction of the molecule.
	 *
	 * @param name The chemical name to parse.
	 * @param verbose Whether to print lots of debugging information to stdin and stderr or not.
	 * @return A RenderedImage containing the parsed molecule, or null if the molecule would not parse.
	 */
	public RenderedImage parseToDepiction(String name, boolean verbose) {
		OpsinResult result = n2s.parseChemicalName(name, verbose);
		if (result.getStructure() == null){
			return null;
		}
		else{
			RenderedImage depiction = null;
			try{
				depiction = opsinFragmentToDepiction(result.getStructure(), verbose);
			}
			catch (Exception e) {
				if (verbose){
					e.printStackTrace();
				}
				return null;
			}
			if (depiction ==null){return null;}//depiction failed
			return depiction;
		}
	}

	private RenderedImage opsinFragmentToDepiction(Fragment frag, boolean verbose) {
		OpsinToChemKitWrapper chemKitWrapper = new OpsinToChemKitWrapper(frag);
		CMLMolecule mol = chemKitWrapper.getChemKitMolecule();
		//Make explicit hydrogens a property of the atom to which they are attached for easier to understand depictions
		HydrogenCalculator.removeExplicitHydrogens(mol);
		Layout2D layout = new Layout2D();
		layout.generateCoordinates(mol);
		Renderer r = new Renderer(mol);
		return r.createImage(new Dimension(640,480));
	}

	/**Run OPSIN as a standalone component for depiction generation
	 *
	 * @param args
	 * @throws Exception
	 */
	public static void main(String [] args) throws Exception {
		NameToDepiction ntd = new NameToDepiction();
		boolean end = false;
		BufferedReader stdinReader = new BufferedReader(new InputStreamReader(System.in));
		System.out.println("OPSIN Prealpha: enter chemical name:");
		while(!end) {
			String name = stdinReader.readLine();
			if(name == null) {
				System.err.println("Disconnected!");
				end = true;
			} else if(name.equals("END")) {
				end = true;
			} else {
				final RenderedImage output = ntd.parseToDepiction(name, true);
				if(output == null) {
					System.out.println("Did not parse.");
					System.out.flush();
				} else {
			        JPanel panel = new JPanel() {
			            /**
						 *
						 */
						private static final long serialVersionUID = 1L;

						protected void paintComponent(Graphics g) {
			                super.paintComponent(g);
			                ((Graphics2D)g).drawImage((BufferedImage) output, null, 0, 0);
			            }
			        };
			        panel.setPreferredSize(new Dimension(640,480));

			        final JFrame frame = new JFrame();
			        frame.setResizable(false);
			        frame.add(panel);
			        frame.pack();
			        SwingUtilities.invokeAndWait(
			            new Runnable() {
			                public void run() {
			                    frame.setVisible(true);
			                }
			            }
			        );

			        while (frame.isVisible()) {
			            Thread.sleep(1000);
			        }
			        frame.dispose();
					System.out.flush();
				}
			}
		}
	}
}
