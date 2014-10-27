/*
 * ExportPlugin.java
 *
 * Created on December 18, 2009
 *
 */
package net.maizegenetics.analysis.data;

import net.maizegenetics.dna.snp.ExportUtils;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.gui.DialogUtils;
import net.maizegenetics.gui.GenotypeTableMask;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.taxa.distance.DistanceMatrix;
import net.maizegenetics.taxa.distance.WriteDistanceMatrix;
import net.maizegenetics.prefs.TasselPrefs;
import net.maizegenetics.tassel.TASSELMainFrame;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.phenotype.PhenotypeUtils;
import net.maizegenetics.util.*;
import org.apache.log4j.Logger;

import javax.swing.*;
import javax.swing.tree.DefaultMutableTreeNode;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.net.URL;

/**
 *
 * @author Terry Casstevens
 */
public class ExportPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(ExportPlugin.class);
    private FileLoadPlugin.TasselFileType myFileType = FileLoadPlugin.TasselFileType.Hapmap;
    private String mySaveFile = null;
    private boolean myIsDiploid = false;
    private boolean myKeepDepth = true;
    private final JFileChooser myFileChooserSave = new JFileChooser(TasselPrefs.getSaveDir());

    /**
     * Creates a new instance of ExportPlugin
     */
    public ExportPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    public DataSet performFunction(DataSet input) {

        try {

            if (input.getSize() != 1) {
                String message = "Please select one and only one item.";
                if (isInteractive()) {
                    JOptionPane.showMessageDialog(getParentFrame(), message);
                } else {
                    myLogger.error("performFunction: " + message);
                }
                return null;
            }

            String filename = mySaveFile;
            try {
                Object data = input.getData(0).getData();

                if (data instanceof GenotypeTable) {
                    filename = performFunctionForAlignment((GenotypeTable) data);
                } else if (data instanceof Phenotype) {
                    filename = performFunctionForPhenotype((Phenotype) data);
                } else if (data instanceof DistanceMatrix) {
                    filename = performFunctionForDistanceMatrix((DistanceMatrix) data);
                } else if (data instanceof TableReport) {
                    filename = performFunctionForTableReport((TableReport) data);
                } else if (data instanceof Report) {
                    filename = performFunctionForReport((Report) data);
                } else {
                    String message = "Don't know how to export data type: " + data.getClass().getName();
                    if (isInteractive()) {
                        JOptionPane.showMessageDialog(getParentFrame(), message);
                    } else {
                        myLogger.error("performFunction: " + message);
                    }
                    return null;
                }
            } catch (Exception e) {
                e.printStackTrace();
                StringBuilder builder = new StringBuilder();
                builder.append(Utils.shortenStrLineLen(ExceptionUtils.getExceptionCauses(e), 50));
                String str = builder.toString();
                if (isInteractive()) {
                    DialogUtils.showError(str, getParentFrame());
                } else {
                    myLogger.error(str);
                }

                return null;
            }

            if (filename != null) {
                myLogger.info("performFunction: wrote dataset: " + input.getData(0).getName() + " to file: " + filename);
                return new DataSet(new Datum("Filename", filename, null), this);
            } else {
                return null;
            }

        } finally {
            fireProgress(100);
        }

    }

    public String performFunctionForDistanceMatrix(DistanceMatrix input) {

        if (isInteractive()) {
            setSaveFile(getFileByChooser());
        }

        if ((mySaveFile == null) || (mySaveFile.length() == 0)) {
            return null;
        }

        try {
            File theFile = new File(Utils.addSuffixIfNeeded(mySaveFile, ".txt"));
            WriteDistanceMatrix.saveDelimitedDistanceMatrix(input, theFile);
            return theFile.getCanonicalPath();
        } catch (Exception e) {
            e.printStackTrace();
            throw new IllegalStateException("ExportPlugin: performFunctionForDistanceMatrix: Problem writing file: " + mySaveFile);
        }

    }

    public String performFunctionForTableReport(TableReport input) {

        if (isInteractive()) {
            setSaveFile(getFileByChooser());
        }

        if ((mySaveFile == null) || (mySaveFile.length() == 0)) {
            return null;
        }

        try {
            File theFile = new File(Utils.addSuffixIfNeeded(mySaveFile, ".txt"));
            TableReportUtils.saveDelimitedTableReport(input, "\t", theFile);
            return theFile.getCanonicalPath();
        } catch (Exception e) {
            e.printStackTrace();
            throw new IllegalStateException("ExportPlugin: performFunctionForTableReport: Problem writing file: " + mySaveFile);
        }

    }

    public String performFunctionForPhenotype(Phenotype input) {

        if (isInteractive()) {
            setSaveFile(getFileByChooser());
        }

        if ((mySaveFile == null) || (mySaveFile.length() == 0)) {
            return null;
        }

        String filename = "";
        try {
            filename = Utils.addSuffixIfNeeded(mySaveFile, ".txt");
            PhenotypeUtils.write(input, filename);
            return new File(filename).getCanonicalPath();
        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("ExportPlugin: performFunctionForPhenotype: Problem writing file: " + filename);
        }

    }

    public String performFunctionForAlignment(GenotypeTable inputAlignment) {

        if (isInteractive()) {
            ExportPluginDialog theDialog = new ExportPluginDialog();
            theDialog.setLocationRelativeTo(getParentFrame());
            theDialog.setVisible(true);
            if (theDialog.isCancel()) {
                return null;
            }
            myFileType = theDialog.getTasselFileType();
            myKeepDepth = theDialog.keepDepth();

            theDialog.dispose();

            setSaveFile(getFileByChooser());
        }

        if ((mySaveFile == null) || (mySaveFile.length() == 0)) {
            return null;
        }

        String resultFile = mySaveFile;

        if ((myFileType == FileLoadPlugin.TasselFileType.Hapmap) || (myFileType == FileLoadPlugin.TasselFileType.HapmapDiploid)) {
            int n = 0;
            DefaultMutableTreeNode node = null;
            if (isInteractive()) {
                DiploidOptionDialog diploidDialog = new DiploidOptionDialog();
                diploidDialog.setLocationRelativeTo(getParentFrame());
                diploidDialog.setVisible(true);
                myIsDiploid = diploidDialog.getDiploid();
                node = (DefaultMutableTreeNode) ((TASSELMainFrame) this.getParentFrame()).getDataTreePanel().getTree().getLastSelectedPathComponent();
                n = node.getChildCount();
            } else {
                if (myFileType == FileLoadPlugin.TasselFileType.Hapmap) {
                    myIsDiploid = false;
                } else if (myFileType == FileLoadPlugin.TasselFileType.HapmapDiploid) {
                    myIsDiploid = true;
                }
            }

            boolean foundImputed = false;
            if ((n == 0) || (!isInteractive())) {
                resultFile = ExportUtils.writeToHapmap(inputAlignment, myIsDiploid, mySaveFile, '\t', this);
            } else {
                int i = 0;
                while (i < n && !foundImputed) {
                    DefaultMutableTreeNode currentNode = (DefaultMutableTreeNode) node.getChildAt(i);
                    Datum currentDatum = (Datum) currentNode.getUserObject();
                    Object currentMask = currentDatum.getData();
                    if (currentMask instanceof GenotypeTableMask) {
                        GenotypeTableMask.MaskType maskType = ((GenotypeTableMask) currentMask).getMaskType();
                        if (maskType == GenotypeTableMask.MaskType.imputed) {
                            ImputeDisplayOptionDialog imputeOptionDialog = new ImputeDisplayOptionDialog();
                            imputeOptionDialog.setLocationRelativeTo(getParentFrame());
                            imputeOptionDialog.setVisible(true);
                            if (imputeOptionDialog.getDisplayImputed()) {
                                //TODO emailed Terry about whether to keep
//                                resultFile = ExportUtils.writeToHapmap(inputAlignment, (GenotypeTableMask) currentMask, myIsDiploid, mySaveFile, '\t', this);
//                                foundImputed = true;
                            } else if (i == (n - 1)) {
                                resultFile = ExportUtils.writeToHapmap(inputAlignment, myIsDiploid, mySaveFile, '\t', this);
                            }
                        }
                    }
                    i++;
                }
            }
        } else if (myFileType == FileLoadPlugin.TasselFileType.Plink) {
            resultFile = ExportUtils.writeToPlink(inputAlignment, mySaveFile, '\t');
        } else if (myFileType == FileLoadPlugin.TasselFileType.Phylip_Seq) {
            PrintWriter out = null;
            try {
                resultFile = Utils.addSuffixIfNeeded(mySaveFile, ".phy");
                out = new PrintWriter(new FileWriter(resultFile));
                ExportUtils.printSequential(inputAlignment, out);
            } catch (Exception e) {
                throw new IllegalStateException("ExportPlugin: performFunction: Problem writing file: " + mySaveFile);
            } finally {
                out.flush();
                out.close();
            }
        } else if (myFileType == FileLoadPlugin.TasselFileType.Phylip_Inter) {
            PrintWriter out = null;
            try {
                resultFile = Utils.addSuffixIfNeeded(mySaveFile, ".phy");
                out = new PrintWriter(new FileWriter(resultFile));
                ExportUtils.printInterleaved(inputAlignment, out);
            } catch (Exception e) {
                throw new IllegalStateException("ExportPlugin: performFunction: Problem writing file: " + mySaveFile);
            } finally {
                out.flush();
                out.close();
            }
        } else if (myFileType == FileLoadPlugin.TasselFileType.Table) {
            resultFile = ExportUtils.saveDelimitedAlignment(inputAlignment, "\t", mySaveFile);
        } else if (myFileType == FileLoadPlugin.TasselFileType.Serial) {
            resultFile = ExportUtils.writeAlignmentToSerialGZ(inputAlignment, mySaveFile);
        } else if (myFileType == FileLoadPlugin.TasselFileType.HDF5) {
            resultFile = ExportUtils.writeGenotypeHDF5(inputAlignment, mySaveFile, myKeepDepth);
        } else if (myFileType == FileLoadPlugin.TasselFileType.VCF) {
            resultFile = ExportUtils.writeToVCF(inputAlignment, mySaveFile, myKeepDepth);
        } else {
            throw new IllegalStateException("ExportPlugin: performFunction: Unknown Alignment File Format: " + myFileType);
        }

        return resultFile;

    }

    public String performFunctionForReport(Report input) {

        if (isInteractive()) {
            ReportOptionDialog theDialog = new ReportOptionDialog();
            theDialog.setLocationRelativeTo(getParentFrame());
            theDialog.setVisible(true);
            if (theDialog.isCancel()) {
                return null;
            }
            myFileType = theDialog.getTasselFileType();

            theDialog.dispose();

            setSaveFile(getFileByChooser());
        }

        if ((mySaveFile == null) || (mySaveFile.length() == 0)) {
            return null;
        }

        String resultFile = Utils.addSuffixIfNeeded(mySaveFile, ".txt");
        if (myFileType == FileLoadPlugin.TasselFileType.Text) {
            BufferedWriter writer = Utils.getBufferedWriter(resultFile);
            try {
                writer.append(input.toString());
            } catch (Exception e) {
                e.printStackTrace();
                throw new IllegalStateException("ExportPlugin: performFunctionForReport: Problem writing file: " + resultFile);
            } finally {
                try {
                    writer.close();
                } catch (Exception e) {
                    // do nothing
                }
            }
        } else {
            PrintWriter writer = null;
            try {
                writer = new PrintWriter(resultFile);
                input.report(writer);
            } catch (Exception e) {
                e.printStackTrace();
                throw new IllegalStateException("ExportPlugin: performFunctionForReport: Problem writing file: " + resultFile);
            } finally {
                try {
                    writer.close();
                } catch (Exception e) {
                    // do nothing
                }
            }
        }
        return resultFile;

    }

    /**
     * Icon for this plugin to be used in buttons, etc.
     *
     * @return ImageIcon
     */
    public ImageIcon getIcon() {
        URL imageURL = ExportPlugin.class.getResource("/net/maizegenetics/analysis/images/Export16.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    /**
     * Button name for this plugin to be used in buttons, etc.
     *
     * @return String
     */
    public String getButtonName() {
        return "Export";
    }

    /**
     * Tool Tip Text for this plugin
     *
     * @return String
     */
    public String getToolTipText() {
        return "Export data to files on your computer.";
    }

    public String getSaveFile() {
        return mySaveFile;
    }

    public void setSaveFile(String saveFile) {
        mySaveFile = saveFile;
    }

    public void setSaveFile(File saveFile) {

        if (saveFile == null) {
            mySaveFile = null;
        } else {
            mySaveFile = saveFile.getPath();
        }

    }

    public void setAlignmentFileType(FileLoadPlugin.TasselFileType type) {
        myFileType = type;
    }

    public void setIsDiploid(boolean isDiploid) {
        myIsDiploid = isDiploid;
    }

    private File getFileByChooser() {
        myFileChooserSave.setMultiSelectionEnabled(false);
        File result = null;
        int returnVal = myFileChooserSave.showSaveDialog(getParentFrame());
        if (returnVal == JFileChooser.SAVE_DIALOG || returnVal == JFileChooser.APPROVE_OPTION) {
            result = myFileChooserSave.getSelectedFile();
            TasselPrefs.putSaveDir(myFileChooserSave.getCurrentDirectory().getPath());
        }
        return result;
    }

    class ExportPluginDialog extends JDialog {

        private boolean myIsCancel = true;
        private ButtonGroup myButtonGroup = new ButtonGroup();
        private JRadioButton myHapMapRadioButton = new JRadioButton("Write Hapmap");
        private JRadioButton myByteHDF5RadioButton = new JRadioButton("Write HDF5");
        private JRadioButton myVCFRadioButton = new JRadioButton("Write VCF");
        private JRadioButton myPlinkRadioButton = new JRadioButton("Write Plink");
        private JRadioButton myPhylipRadioButton = new JRadioButton("Write Phylip (Sequential)");
        private JRadioButton myPhylipInterRadioButton = new JRadioButton("Write Phylip (Interleaved)");
        private JRadioButton myTabTableRadioButton = new JRadioButton("Write Tab Delimited");

        private JCheckBox myKeepDepthCheck = new JCheckBox("Keep Depth (VCF or HDF5)", true);

        public ExportPluginDialog() {
            super((Frame) null, "Export...", true);
            try {
                jbInit();
                pack();
            } catch (Exception ex) {
                ex.printStackTrace();
            }
        }

        private void jbInit() throws Exception {

            setTitle("Export...");
            setDefaultCloseOperation(JDialog.HIDE_ON_CLOSE);
            setUndecorated(false);
            getRootPane().setWindowDecorationStyle(JRootPane.NONE);
            Container contentPane = getContentPane();
            BoxLayout layout = new BoxLayout(contentPane, BoxLayout.Y_AXIS);
            contentPane.setLayout(layout);
            JPanel main = getMain();
            contentPane.add(main);
            pack();
            setResizable(false);
            myButtonGroup.add(myHapMapRadioButton);
            myButtonGroup.add(myByteHDF5RadioButton);
            myButtonGroup.add(myVCFRadioButton);
            myButtonGroup.add(myPlinkRadioButton);
            myButtonGroup.add(myPhylipRadioButton);
            myButtonGroup.add(myPhylipInterRadioButton);
            myButtonGroup.add(myTabTableRadioButton);
            myHapMapRadioButton.setSelected(true);

        }

        private JPanel getMain() {
            JPanel inputs = new JPanel();
            BoxLayout layout = new BoxLayout(inputs, BoxLayout.Y_AXIS);
            inputs.setLayout(layout);
            inputs.setAlignmentX(JPanel.CENTER_ALIGNMENT);
            inputs.add(Box.createRigidArea(new Dimension(1, 10)));
            inputs.add(getLabel());
            inputs.add(Box.createRigidArea(new Dimension(1, 10)));
            inputs.add(getFileTypePanel());
            inputs.add(Box.createRigidArea(new Dimension(1, 10)));
            inputs.add(getOptionPanel());
            inputs.add(Box.createRigidArea(new Dimension(1, 5)));
            inputs.add(getButtons());
            inputs.add(Box.createRigidArea(new Dimension(1, 10)));
            return inputs;
        }

        private JPanel getLabel() {
            JPanel result = new JPanel();
            BoxLayout layout = new BoxLayout(result, BoxLayout.Y_AXIS);
            result.setLayout(layout);
            result.setAlignmentX(JPanel.CENTER_ALIGNMENT);
            JLabel jLabel1 = new JLabel("Choose File Type to Export.");
            jLabel1.setFont(new Font("Dialog", Font.BOLD, 18));
            result.add(jLabel1);
            return result;
        }

        private JPanel getFileTypePanel() {
            JPanel result = new JPanel();
            BoxLayout layout = new BoxLayout(result, BoxLayout.Y_AXIS);
            result.setLayout(layout);
            result.setAlignmentX(JPanel.CENTER_ALIGNMENT);
            result.setBorder(BorderFactory.createEtchedBorder());

            result.add(myHapMapRadioButton);
            result.add(myByteHDF5RadioButton);
            result.add(myVCFRadioButton);
            result.add(myPlinkRadioButton);
            result.add(myPhylipRadioButton);
            result.add(myPhylipInterRadioButton);
            result.add(myTabTableRadioButton);

            result.add(Box.createRigidArea(new Dimension(1, 20)));

            return result;

        }

        private JPanel getOptionPanel() {
            JPanel result = new JPanel();
            BoxLayout layout = new BoxLayout(result, BoxLayout.Y_AXIS);
            result.setLayout(layout);
            result.setAlignmentX(JPanel.CENTER_ALIGNMENT);
            result.setBorder(BorderFactory.createEtchedBorder());
            result.add(myKeepDepthCheck);
            result.add(Box.createRigidArea(new Dimension(1, 10)));

            return result;

        }

        private JPanel getButtons() {

            JButton okButton = new JButton();
            JButton cancelButton = new JButton();

            cancelButton.setText("Cancel");
            cancelButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    cancelButton_actionPerformed(e);
                }
            });

            okButton.setText("OK");
            okButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    okButton_actionPerformed(e);
                }
            });

            JPanel result = new JPanel(new FlowLayout(FlowLayout.CENTER));

            result.add(okButton);

            result.add(cancelButton);

            return result;

        }

        public FileLoadPlugin.TasselFileType getTasselFileType() {
            if (myHapMapRadioButton.isSelected()) {
                return FileLoadPlugin.TasselFileType.Hapmap;
            }
            if (myByteHDF5RadioButton.isSelected()) {
                return FileLoadPlugin.TasselFileType.HDF5;
            }
            if (myVCFRadioButton.isSelected()) {
                return FileLoadPlugin.TasselFileType.VCF;
            }
            if (myPlinkRadioButton.isSelected()) {
                return FileLoadPlugin.TasselFileType.Plink;
            }
            if (myPhylipRadioButton.isSelected()) {
                return FileLoadPlugin.TasselFileType.Phylip_Seq;
            }
            if (myPhylipInterRadioButton.isSelected()) {
                return FileLoadPlugin.TasselFileType.Phylip_Inter;
            }
            if (myTabTableRadioButton.isSelected()) {
                return FileLoadPlugin.TasselFileType.Table;
            }
            return null;
        }

        public boolean keepDepth() {
            return myKeepDepthCheck.isSelected();
        }

        private void okButton_actionPerformed(ActionEvent e) {
            myIsCancel = false;
            setVisible(false);
        }

        private void cancelButton_actionPerformed(ActionEvent e) {
            myIsCancel = true;
            setVisible(false);
        }

        public boolean isCancel() {
            return myIsCancel;
        }
    }
}

class ImputeDisplayOptionDialog extends JDialog {

    boolean displayImputed = true;
    private JPanel mainPanel = new JPanel();
    private JLabel lbl = new JLabel();
    private JButton yesButton = new JButton();
    private JButton noButton = new JButton();
    private GridBagLayout gridBagLayout = new GridBagLayout();

    public ImputeDisplayOptionDialog() {
        super((Frame) null, "File Loader", true);
        try {
            initUI();
            pack();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    void initUI() throws Exception {

        lbl.setFont(new java.awt.Font("Dialog", 1, 12));
        lbl.setText("Would you like Imputed data to be exported in lower case?");

        mainPanel.setMinimumSize(new Dimension(480, 150));
        mainPanel.setPreferredSize(new Dimension(480, 150));
        mainPanel.setLayout(gridBagLayout);

        yesButton.setMaximumSize(new Dimension(63, 27));
        yesButton.setMinimumSize(new Dimension(63, 27));
        yesButton.setText("Yes");
        yesButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(ActionEvent e) {
                yesButton_actionPerformed(e);
            }
        });

        noButton.setMaximumSize(new Dimension(63, 27));
        noButton.setMinimumSize(new Dimension(63, 27));
        noButton.setText("No");
        noButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(ActionEvent e) {
                noButton_actionPerformed(e);
            }
        });

        JPanel buttonPanel = new JPanel();
        buttonPanel.add(yesButton);
        buttonPanel.add(noButton);

        mainPanel.add(lbl, new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(5, 0, 5, 0), 0, 0));
        mainPanel.add(buttonPanel, new GridBagConstraints(0, 1, 1, 1, 1.0, 1.0, GridBagConstraints.PAGE_END, GridBagConstraints.NONE, new Insets(5, 0, 5, 0), 0, 0));

        this.add(mainPanel, BorderLayout.CENTER);
    }

    private void yesButton_actionPerformed(ActionEvent e) {
        displayImputed = true;
        setVisible(false);
    }

    private void noButton_actionPerformed(ActionEvent e) {
        displayImputed = false;
        setVisible(false);
    }

    public boolean getDisplayImputed() {
        return displayImputed;
    }
}

class DiploidOptionDialog extends JDialog {

    boolean displayDiploid = true;
    private JPanel mainPanel = new JPanel();
    private JLabel lbl = new JLabel();
    private JButton yesButton = new JButton();
    private JButton noButton = new JButton();
    private GridBagLayout gridBagLayout = new GridBagLayout();

    public DiploidOptionDialog() {
        super((Frame) null, "File Loader", true);
        try {
            initUI();
            pack();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    void initUI() throws Exception {

        lbl.setFont(new java.awt.Font("Dialog", 1, 12));
        lbl.setText("Would you like SNPs to be exported as Diploids?");

        mainPanel.setMinimumSize(new Dimension(480, 150));
        mainPanel.setPreferredSize(new Dimension(480, 150));
        mainPanel.setLayout(gridBagLayout);

        yesButton.setMaximumSize(new Dimension(63, 27));
        yesButton.setMinimumSize(new Dimension(63, 27));
        yesButton.setText("Yes");
        yesButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(ActionEvent e) {
                yesButton_actionPerformed(e);
            }
        });

        noButton.setMaximumSize(new Dimension(63, 27));
        noButton.setMinimumSize(new Dimension(63, 27));
        noButton.setText("No");
        noButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(ActionEvent e) {
                noButton_actionPerformed(e);
            }
        });

        JPanel buttonPanel = new JPanel();
        buttonPanel.add(yesButton);
        buttonPanel.add(noButton);

        mainPanel.add(lbl, new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(5, 0, 5, 0), 0, 0));
        mainPanel.add(buttonPanel, new GridBagConstraints(0, 1, 1, 1, 1.0, 1.0, GridBagConstraints.PAGE_END, GridBagConstraints.NONE, new Insets(5, 0, 5, 0), 0, 0));

        this.add(mainPanel, BorderLayout.CENTER);
    }

    private void yesButton_actionPerformed(ActionEvent e) {
        displayDiploid = true;
        setVisible(false);
    }

    private void noButton_actionPerformed(ActionEvent e) {
        displayDiploid = false;
        setVisible(false);
    }

    public boolean getDiploid() {
        return displayDiploid;
    }
}

class ReportOptionDialog extends JDialog {

    private boolean myIsCancel = true;
    private ButtonGroup myButtonGroup = new ButtonGroup();
    private JRadioButton myReportRadioButton = new JRadioButton("Write As Report");
    private JRadioButton myTextRadioButton = new JRadioButton("Write As Text");

    public ReportOptionDialog() {
        super((Frame) null, "Export Report...", true);
        try {
            jbInit();
            pack();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    private void jbInit() throws Exception {

        setTitle("Export Report...");
        setDefaultCloseOperation(JDialog.HIDE_ON_CLOSE);
        setUndecorated(false);
        getRootPane().setWindowDecorationStyle(JRootPane.NONE);

        Container contentPane = getContentPane();

        BoxLayout layout = new BoxLayout(contentPane, BoxLayout.Y_AXIS);
        contentPane.setLayout(layout);

        JPanel main = getMain();

        contentPane.add(main);

        pack();

        setResizable(false);

        myButtonGroup.add(myReportRadioButton);
        myButtonGroup.add(myTextRadioButton);
        myReportRadioButton.setSelected(true);

    }

    private JPanel getMain() {

        JPanel inputs = new JPanel();
        BoxLayout layout = new BoxLayout(inputs, BoxLayout.Y_AXIS);
        inputs.setLayout(layout);
        inputs.setAlignmentX(JPanel.CENTER_ALIGNMENT);

        inputs.add(Box.createRigidArea(new Dimension(1, 10)));

        inputs.add(getLabel());

        inputs.add(Box.createRigidArea(new Dimension(1, 10)));

        inputs.add(getOptionPanel());

        inputs.add(Box.createRigidArea(new Dimension(1, 10)));

        inputs.add(getButtons());

        inputs.add(Box.createRigidArea(new Dimension(1, 10)));

        return inputs;

    }

    private JPanel getLabel() {

        JPanel result = new JPanel();
        BoxLayout layout = new BoxLayout(result, BoxLayout.Y_AXIS);
        result.setLayout(layout);
        result.setAlignmentX(JPanel.CENTER_ALIGNMENT);

        JLabel jLabel1 = new JLabel("Choose File Type to Export.");
        jLabel1.setFont(new Font("Dialog", Font.BOLD, 18));
        result.add(jLabel1);

        return result;

    }

    private JPanel getOptionPanel() {

        JPanel result = new JPanel();
        BoxLayout layout = new BoxLayout(result, BoxLayout.Y_AXIS);
        result.setLayout(layout);
        result.setAlignmentX(JPanel.CENTER_ALIGNMENT);
        result.setBorder(BorderFactory.createEtchedBorder());

        result.add(myReportRadioButton);
        result.add(myTextRadioButton);

        result.add(Box.createRigidArea(new Dimension(1, 20)));

        return result;

    }

    private JPanel getButtons() {

        JButton okButton = new JButton();
        JButton cancelButton = new JButton();

        cancelButton.setText("Cancel");
        cancelButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                cancelButton_actionPerformed(e);
            }
        });

        okButton.setText("OK");
        okButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                okButton_actionPerformed(e);
            }
        });

        JPanel result = new JPanel(new FlowLayout(FlowLayout.CENTER));

        result.add(okButton);

        result.add(cancelButton);

        return result;

    }

    public FileLoadPlugin.TasselFileType getTasselFileType() {
        if (myTextRadioButton.isSelected()) {
            return FileLoadPlugin.TasselFileType.Text;
        }
        return null;
    }

    private void okButton_actionPerformed(ActionEvent e) {
        myIsCancel = false;
        setVisible(false);
    }

    private void cancelButton_actionPerformed(ActionEvent e) {
        myIsCancel = true;
        setVisible(false);
    }

    public boolean isCancel() {
        return myIsCancel;
    }
}
