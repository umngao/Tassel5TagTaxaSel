/*
 * ProjectionLoadPlugin
 */
package net.maizegenetics.analysis.data;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.io.ProjectionGenotypeIO;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginEvent;
import net.maizegenetics.prefs.TasselPrefs;
import net.maizegenetics.util.Utils;

import javax.swing.*;

import java.awt.Container;
import java.awt.Frame;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ComponentAdapter;
import java.awt.event.ComponentEvent;
import java.io.*;

import org.apache.log4j.Logger;

/**
 *
 * @author Alex Lipka This should enable a used to load a projection alignment
 * using the TASSEL GUI
 */
public class ProjectionLoadPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(ProjectionLoadPlugin.class);
    private String myRecombinationBreakpoints = null;
    private String myHighDensityMarkers = null;
    private String myChromosome = null;
    private boolean hasHetSeparator;
    private boolean hasNucleotides = true;
    private boolean isPhysicalMap = true;
    private String hetSeparator;
    private String missingCharacter;

    /**
     * Creates a new instance of ProjectionLoadPlugin
     */
    public ProjectionLoadPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    public DataSet performFunction(DataSet input) {

        if (isInteractive()) {
            ProjectionPluginDialog theDialog = new ProjectionPluginDialog();
            theDialog.setLocationRelativeTo(getParentFrame());
            theDialog.setVisible(true);
            if (theDialog.isCancel()) {
                return null;
            }
            //hasHetSeparator = theDialog.chkHetsSeparated.isSelected();
            //hasNucleotides = theDialog.chkMarkersAreNucleotides.isSelected();
            //hetSeparator = theDialog.txtHetsep.getText().trim();
            //missingCharacter = theDialog.txtMissing.getText().trim();
            //isPhysicalMap = theDialog.radioPhysical.isSelected();
            theDialog.dispose();
        }

        if ((myRecombinationBreakpoints == null) || (myRecombinationBreakpoints.length() == 0)) {
            return null;
        }
        if ((myHighDensityMarkers == null) || (myHighDensityMarkers.length() == 0)) {
            return null;
        }

        try {
            DataSet result = loadFile(myRecombinationBreakpoints, myHighDensityMarkers);
            return result;
        } catch (Exception e) {
            String msg = "Recombination breakpoints " + myRecombinationBreakpoints + " and high density markers"
                    + myHighDensityMarkers + " failed to load. " + "Make sure the import options are properly set.";
            if (isInteractive()) {
                JOptionPane.showMessageDialog(getParentFrame(), msg, "Error uploading Projection files", JOptionPane.ERROR_MESSAGE);
            } else {
                myLogger.error(msg);
            }
            return null;
        }

    }

    public String getRecombinationBreakpoints() {
        return myRecombinationBreakpoints;
    }

    public String getHighDensityMarkers() {
        return myHighDensityMarkers;
    }

    public String getChromosome() {
        return myChromosome;
    }

    public void setRecombinationBreakpoints(String filename) {
        myRecombinationBreakpoints = filename;
    }

    public void setHighDensityMarkers(String filename) {
        myHighDensityMarkers = filename;
    }

    public void setChromosome(String chromosome) {
        myChromosome = chromosome;
    }

    /**
     * Icon for this plugin to be used in buttons, etc.
     *
     * @return ImageIcon
     */
    public ImageIcon getIcon() {
        return null;
    }

    /**
     * Button name for this plugin to be used in buttons, etc.
     *
     * @return String
     */
    public String getButtonName() {
        return "Load Projection Alignment";
    }

    /**
     * Tool Tip Text for this plugin
     *
     * @return String
     */
    public String getToolTipText() {
        return "Load Projection Alignments";
    }

    public DataSet loadFile(String theRecombinationBreakpoints, String theHighDensityMarkers) {

        GenotypeTable theAlignmentForGenotype = null;
        //Terry - fix this
        try {
            theAlignmentForGenotype = ProjectionGenotypeIO.getInstance(theRecombinationBreakpoints, theHighDensityMarkers);

            //Run ProjectionBuilder() to create the breakpoints
            //ProjectionBuilder theProjectionBuilder =  new ProjectionBuilder(theAlignment);
            //theAlignmentForGenotype = theProjectionBuilder.build();
        } catch (Exception e) {
            e.printStackTrace();
            return null;
        }

        //Begin here, Alex
        Datum td = new Datum(Utils.getFilename(theRecombinationBreakpoints, FileLoadPlugin.FILE_EXT_HAPMAP), theAlignmentForGenotype, null);
        DataSet tds = new DataSet(td, this);
        fireDataSetReturned(new PluginEvent(tds, ProjectionLoadPlugin.class));

        return tds;

    }

    public void hasHetSeparator(boolean hasSeparator) {
        hasHetSeparator = hasSeparator;
    }

    public void hasNucleotides(boolean nucleotides) {
        hasNucleotides = nucleotides;
    }

    public void setHetSeparator(String separator) {
        hetSeparator = separator;
    }

    public void setMissingCharacter(String missing) {
        missingCharacter = missing;
    }

    public void hasPhysicalMap(boolean physicalmap) {
        isPhysicalMap = physicalmap;
    }

    class ProjectionPluginDialog extends JDialog {

        private JPanel main = null;
        private final static int TEXT_FIELD_WIDTH = 30;
        private JTextField myHighDensityMarkersField = null;
        private JTextField myRecombinationBreakpointsField = null;
        private JTextField myChromosomeField = null;
        private boolean myIsCancel = false;
        private JButton myHighDensityMarkersBrowseButton = null;
        private JButton myRecombinationBreakpointsBrowseButton = null;
        private final JCheckBox chkHetsSeparated = new JCheckBox("Heterozygotes are separated by a character (eg, A/T rather than AT)", true);
        private final JCheckBox chkMarkersAreNucleotides = new JCheckBox("Markers are nucleotides (ie, A, C, G, or T)", true);
        private final JTextField txtHetsep = new JTextField(5);
        private final JTextField txtMissing = new JTextField(5);
        private final JRadioButton radioPhysical = new JRadioButton("Physical", true);
        private final JRadioButton radioGenetic = new JRadioButton("Genetic", false);

        public ProjectionPluginDialog() {
            super((Frame) null, "File Loader", true);
            try {
                createDialog();
                pack();
            } catch (Exception ex) {
                ex.printStackTrace();
            }
        }

        private void createDialog() {

            myRecombinationBreakpointsField = new JTextField(TEXT_FIELD_WIDTH);
            myRecombinationBreakpointsField.setText("(Select Recombination Breakpoints File)");

            myHighDensityMarkersField = new JTextField(TEXT_FIELD_WIDTH);
            myHighDensityMarkersField.setText("(Select High Density Map)");

            myRecombinationBreakpointsBrowseButton = new JButton("Browse...");
            myHighDensityMarkersBrowseButton = new JButton("Browse...");
            //txtHetsep.setText("/");
            //txtMissing.setText("-");

            setLocationRelativeTo(getParentFrame());

            setTitle("Load Projection Alignment");
            setDefaultCloseOperation(JDialog.HIDE_ON_CLOSE);
            setUndecorated(false);
            getRootPane().setWindowDecorationStyle(JRootPane.NONE);

            Container contentPane = getContentPane();

            contentPane.setLayout(new GridBagLayout());

            addComponentListener(new ComponentAdapter() {

                public void componentShown(ComponentEvent ce) {
                    myRecombinationBreakpointsField.requestFocusInWindow();
                }
            });

            myRecombinationBreakpointsBrowseButton.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent e) {
                    JFileChooser myFileChooser = new JFileChooser(TasselPrefs.getOpenDir());
                    myFileChooser.setDialogTitle("Open a Recombination Breakpoint File");
                    if (myFileChooser.showOpenDialog(main) == JFileChooser.APPROVE_OPTION) {
                        File file = myFileChooser.getSelectedFile();
                        myRecombinationBreakpointsField.setText(file.getPath());
                        TasselPrefs.putOpenDir(myFileChooser.getCurrentDirectory().getPath());
                    }
                }
            });

            myHighDensityMarkersBrowseButton.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent e) {
                    JFileChooser myFileChooser = new JFileChooser(TasselPrefs.getOpenDir());
                    myFileChooser.setDialogTitle("Open a High-density Map File");
                    if (myFileChooser.showOpenDialog(main) == JFileChooser.APPROVE_OPTION) {
                        File file = myFileChooser.getSelectedFile();
                        myHighDensityMarkersField.setText(file.getPath());
                        TasselPrefs.putOpenDir(myFileChooser.getCurrentDirectory().getPath());
                    }
                }
            });

            GridBagConstraints gbc = new GridBagConstraints();

            gbc.gridx = 0;
            gbc.gridy = 0;
            gbc.insets = new Insets(6, 15, 4, 4);
            gbc.anchor = GridBagConstraints.WEST;
            contentPane.add(new JLabel("Recombination Breakpoint File"), gbc);
            gbc.gridx++;
            gbc.gridwidth = 2;
            gbc.insets.left = 4;
            gbc.anchor = GridBagConstraints.CENTER;
            contentPane.add(myRecombinationBreakpointsField, gbc);
            gbc.gridx += 2;
            gbc.gridwidth = 1;
            gbc.insets.right = 10;
            gbc.anchor = GridBagConstraints.WEST;
            contentPane.add(myRecombinationBreakpointsBrowseButton, gbc);

            gbc.gridy++;
            gbc.gridx = 0;
            gbc.anchor = GridBagConstraints.WEST;
            gbc.insets.left = 15;
            contentPane.add(new JLabel("High-Density Marker File"), gbc);
            gbc.gridx++;
            gbc.gridwidth = 2;
            gbc.insets.left = 4;
            gbc.anchor = GridBagConstraints.CENTER;
            contentPane.add(myHighDensityMarkersField, gbc);
            gbc.gridx += 2;
            gbc.gridwidth = 1;
            gbc.insets.right = 10;
            gbc.anchor = GridBagConstraints.WEST;
            contentPane.add(myHighDensityMarkersBrowseButton, gbc);

            JPanel buttonPanel = new JPanel();
            buttonPanel.setLayout(new BoxLayout(buttonPanel, BoxLayout.X_AXIS));
            buttonPanel.add(getOkButton());
            buttonPanel.add(Box.createHorizontalStrut(30));
            buttonPanel.add(getCancelButton());

            gbc.gridy++;
            gbc.gridx = 0;
            gbc.gridwidth = 6;
            gbc.insets.top = 25;
            gbc.anchor = GridBagConstraints.CENTER;
            gbc.fill = GridBagConstraints.NONE;
            contentPane.add(buttonPanel, gbc);

            pack();
        }

        private JButton getCancelButton() {

            JButton cancelButton = new JButton("Cancel");
            cancelButton.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent e) {
                    myIsCancel = true;
                    setVisible(false);
                }
            });

            return cancelButton;

        }

        private JButton getOkButton() {

            JButton okButton = new JButton();
            okButton.setText("Import");
            okButton.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent e) {
                    myRecombinationBreakpoints = myRecombinationBreakpointsField.getText();
                    myHighDensityMarkers = myHighDensityMarkersField.getText();
//                    myChromosome = myChromosomeField.getText();
                    myIsCancel = false;
                    setVisible(false);
                }
            });

            return okButton;

        }

        public boolean isCancel() {
            return myIsCancel;
        }
    }
}
