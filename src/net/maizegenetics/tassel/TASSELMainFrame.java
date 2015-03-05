/*
 * TASSEL - Trait Analysis by a aSSociation Evolution & Linkage
 * Copyright (C) 2003 Ed Buckler
 *
 * This software evaluates linkage disequilibrium nucletide diversity and
 * associations. For more information visit http://www.maizegenetics.net
 *
 * This software is distributed under GNU general public license, version 2 and without
 * any warranty or technical support.
 *
 * You can redistribute and/or modify it under the terms of GNU General
 * public license.
 *
 */
//Title:      TASSELMainFrame
//Version:
//Copyright:  Copyright (c) 1997
//Author:     Ed Buckler
//Company:    NCSU
package net.maizegenetics.tassel;

import net.maizegenetics.analysis.popgen.SequenceDiversityPlugin;
import net.maizegenetics.analysis.data.ProjectionLoadPlugin;
import net.maizegenetics.analysis.distance.KinshipPlugin;
import net.maizegenetics.analysis.chart.TableDisplayPlugin;
import net.maizegenetics.analysis.chart.Grid2dDisplayPlugin;
import net.maizegenetics.analysis.chart.ManhattanDisplayPlugin;
import net.maizegenetics.analysis.chart.QQDisplayPlugin;
import net.maizegenetics.analysis.association.EqtlAssociationPlugin;
import net.maizegenetics.analysis.association.FixedEffectLMPlugin;
import net.maizegenetics.analysis.association.MLMPlugin;
import net.maizegenetics.analysis.data.PlinkLoadPlugin;
import net.maizegenetics.analysis.popgen.LinkageDiseqDisplayPlugin;
import net.maizegenetics.analysis.popgen.LinkageDisequilibriumPlugin;
import net.maizegenetics.analysis.data.MergeGenotypeTablesPlugin;
import net.maizegenetics.analysis.data.PrincipalComponentsPlugin;
import net.maizegenetics.analysis.data.UnionAlignmentPlugin;
import net.maizegenetics.analysis.data.FileLoadPlugin;
import net.maizegenetics.analysis.data.IntersectionAlignmentPlugin;
import net.maizegenetics.analysis.data.GenotypeSummaryPlugin;
import net.maizegenetics.analysis.data.ExportPlugin;
import net.maizegenetics.analysis.data.SeparatePlugin;
import net.maizegenetics.analysis.data.SynonymizerPlugin;
import net.maizegenetics.analysis.filter.FilterTaxaAlignmentPlugin;
import net.maizegenetics.analysis.filter.FilterTaxaPropertiesPlugin;
import net.maizegenetics.analysis.filter.FilterSiteNamePlugin;
import net.maizegenetics.analysis.filter.FilterAlignmentPlugin;
import net.maizegenetics.analysis.filter.FilterTraitsPlugin;
import net.maizegenetics.analysis.tree.CreateTreePlugin;
import net.maizegenetics.analysis.tree.ArchaeopteryxPlugin;
import net.maizegenetics.analysis.chart.ChartDisplayPlugin;
import net.maizegenetics.analysis.association.RidgeRegressionEmmaPlugin;
import net.maizegenetics.analysis.modelfitter.StepwiseOLSModelFitterPlugin;
import net.maizegenetics.analysis.numericaltransform.NumericalGenotypePlugin;
import net.maizegenetics.analysis.numericaltransform.TransformDataPlugin;
import net.maizegenetics.gui.PrintHeapAction;
import net.maizegenetics.plugindef.*;
import net.maizegenetics.prefs.TasselPrefs;
import net.maizegenetics.progress.ProgressPanel;
import net.maizegenetics.util.Utils;

import org.apache.log4j.Logger;

import javax.swing.*;

import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.WindowEvent;
import java.io.*;
import java.net.URL;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import net.maizegenetics.analysis.data.GetPositionListPlugin;
import net.maizegenetics.analysis.data.GetTaxaListPlugin;
import net.maizegenetics.analysis.data.HetsToUnknownPlugin;
import net.maizegenetics.analysis.data.ProjectPcsAndRunModelSelectionPlugin;
import net.maizegenetics.analysis.data.SortGenotypeFilePlugin;
import net.maizegenetics.analysis.distance.DistanceMatrixPlugin;
import net.maizegenetics.analysis.gbs.BinaryToTextPlugin;
import net.maizegenetics.analysis.gbs.DiscoverySNPCallerPlugin;
import net.maizegenetics.analysis.gbs.FastqToTagCountPlugin;
import net.maizegenetics.analysis.gbs.MergeMultipleTagCountPlugin;
import net.maizegenetics.analysis.gbs.ModifyTBTHDF5Plugin;
import net.maizegenetics.analysis.gbs.ProductionSNPCallerPlugin;
import net.maizegenetics.analysis.gbs.SAMConverterPlugin;
import net.maizegenetics.analysis.gbs.SeqToTBTHDF5Plugin;
import net.maizegenetics.analysis.gbs.TagCountToFastqPlugin;
import net.maizegenetics.analysis.gbs.UTagCountToTagPairPlugin;
import net.maizegenetics.analysis.gbs.UTagPairToTOPMPlugin;
import net.maizegenetics.analysis.imputation.FILLINFindHaplotypesPlugin;
import net.maizegenetics.analysis.imputation.FILLINImputationPlugin;
import net.maizegenetics.analysis.imputation.FSFHapImputationPlugin;
import net.maizegenetics.analysis.imputation.RemoveIndelsForBeaglePlugin;
import net.maizegenetics.analysis.numericaltransform.ImputationPlugin;

/**
 * TASSELMainFrame
 *
 */
public class TASSELMainFrame extends JFrame implements ActionListener {

    private static final Logger myLogger = Logger.getLogger(TASSELMainFrame.class);
    public static final String version = "5.2.5";
    public static final String versionDate = "March 5, 2015";
    private DataTreePanel myDataTreePanel;
    private String tasselDataFile = "TasselDataFile";
    //a variable to control when the progress bar was last updated
    private JFileChooser filerSave = new JFileChooser();
    private JFileChooser filerOpen = new JFileChooser();
    private JScrollPane reportPanelScrollPane = new JScrollPane();
    private JTextArea reportPanelTextArea = new JTextArea();
    JScrollPane mainPanelScrollPane = new JScrollPane();
    JPanel mainDisplayPanel = new JPanel();
    private JTextArea mainPanelTextArea = new JTextArea();
    private JTextField myStatusTextField = new JTextField();
    private JMenuItem openCompleteDataTreeMenuItem = new JMenuItem();
    private JMenuItem openDataMenuItem = new JMenuItem();
    private JMenuItem saveCompleteDataTreeMenuItem = new JMenuItem();
    private JMenuItem saveDataTreeAsMenuItem = new JMenuItem();
    private PreferencesDialog thePreferencesDialog;
    private final ProgressPanel myProgressPanel = ProgressPanel.getInstance();
    private HashMap<JMenuItem, Plugin> myMenuItemHash = new HashMap<JMenuItem, Plugin>();

    public TASSELMainFrame() {
        try {
            loadSettings();
            myDataTreePanel = new DataTreePanel(this);
            myDataTreePanel.setToolTipText("Data Tree Panel");
            addMenuBar();
            initializeMyFrame();

            this.setTitle("TASSEL (Trait Analysis by aSSociation, Evolution, and Linkage) " + version);

            TasselLogging.basicLoggingInfo();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    //Component initialization
    private void initializeMyFrame() throws Exception {
        getContentPane().setLayout(new BorderLayout());

        int xDim = TasselPrefs.getXDim();
        int yDim = TasselPrefs.getYDim();
        if ((xDim < 50) || (yDim < 50)) {
            Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
            setSize(new Dimension(screenSize.width * 19 / 20, screenSize.height * 19 / 20));
        } else {
            setSize(xDim, yDim);
        }
        setTitle("TASSEL (Trait Analysis by aSSociation, Evolution, and Linkage)");
        addWindowListener(new java.awt.event.WindowAdapter() {
            public void windowClosing(WindowEvent e) {
                TasselLogging.closeInstance();
                TasselPrefs.putXDim(getWidth());
                TasselPrefs.putYDim(getHeight());
                System.exit(0);
            }
        });

        filerSave.setDialogType(JFileChooser.SAVE_DIALOG);

        JSplitPane dataTreeReportPanelsSplitPanel = new JSplitPane();
        dataTreeReportPanelsSplitPanel.setOrientation(JSplitPane.VERTICAL_SPLIT);

        reportPanelTextArea.setEditable(false);
        reportPanelTextArea.setToolTipText("Report Panel");
        reportPanelTextArea.setLineWrap(true);
        reportPanelTextArea.setWrapStyleWord(true);
        mainPanelTextArea.setDoubleBuffered(true);
        mainPanelTextArea.setEditable(false);
        mainPanelTextArea.setFont(new java.awt.Font("Monospaced", 0, 12));
        mainPanelTextArea.setToolTipText("Main Panel");

        myStatusTextField.setBackground(Color.lightGray);
        myStatusTextField.setBorder(null);
        saveCompleteDataTreeMenuItem.setText("Save Data Tree");
        saveCompleteDataTreeMenuItem.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                saveCompleteDataTreeMenuItem_actionPerformed(e);
            }
        });
        saveDataTreeAsMenuItem.setText("Save Data Tree As ...");
        saveDataTreeAsMenuItem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(ActionEvent e) {
                saveDataTreeMenuItem_actionPerformed(e);
            }
        });
        openCompleteDataTreeMenuItem.setText("Open Data Tree");
        openCompleteDataTreeMenuItem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(ActionEvent e) {

                openCompleteDataTreeMenuItem_actionPerformed(e);
            }
        });

        openDataMenuItem.setText("Open Data Tree...");
        openDataMenuItem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(ActionEvent e) {

                openDataMenuItem_actionPerformed(e);
            }
        });

        JSplitPane dataTreeReportMainPanelsSplitPanel = new JSplitPane();

        getContentPane().add(dataTreeReportMainPanelsSplitPanel, BorderLayout.CENTER);

        dataTreeReportMainPanelsSplitPanel.add(dataTreeReportPanelsSplitPanel, JSplitPane.LEFT);

        dataTreeReportPanelsSplitPanel.add(myDataTreePanel, JSplitPane.TOP);

        JSplitPane reportProgressSplitPane = new JSplitPane(JSplitPane.VERTICAL_SPLIT);

        JPanel reportPanel = new JPanel(new BorderLayout());

        reportProgressSplitPane.add(reportPanel, JSplitPane.TOP);

        reportPanel.add(reportPanelScrollPane, BorderLayout.CENTER);

        reportPanelScrollPane.getViewport().add(reportPanelTextArea, null);

        reportProgressSplitPane.add(new JScrollPane(myProgressPanel), JSplitPane.BOTTOM);

        dataTreeReportPanelsSplitPanel.add(reportProgressSplitPane, JSplitPane.BOTTOM);

        dataTreeReportMainPanelsSplitPanel.add(mainDisplayPanel, JSplitPane.RIGHT);

        mainDisplayPanel.setLayout(new BorderLayout());
        mainPanelScrollPane.getViewport().add(mainPanelTextArea, null);
        mainDisplayPanel.add(mainPanelScrollPane, BorderLayout.CENTER);

        mainPanelScrollPane.getViewport().add(mainPanelTextArea, null);

        getContentPane().add(myStatusTextField, BorderLayout.SOUTH);

        dataTreeReportMainPanelsSplitPanel.setDividerLocation(getSize().width / 4);

        dataTreeReportPanelsSplitPanel.setDividerLocation((int) (getSize().height / 3.5));

        reportProgressSplitPane.setDividerLocation((int) (getSize().height / 3.5));
    }

    private void addMenuBar() {

        JMenuBar jMenuBar = new JMenuBar();

        jMenuBar.add(getFileMenu());
        jMenuBar.add(getDataMenu());
        jMenuBar.add(getImputeMenu());
        jMenuBar.add(getFiltersMenu());
        jMenuBar.add(getAnalysisMenu());
        jMenuBar.add(getResultsMenu());
        jMenuBar.add(getGBSMenu());
        jMenuBar.add(Box.createHorizontalGlue());
        jMenuBar.add(getHelpMenu());

        this.setJMenuBar(jMenuBar);

    }

    //Help | About action performed
    private void helpAbout_actionPerformed(ActionEvent e) {
        AboutBox dlg = new AboutBox(this);
        Dimension dlgSize = dlg.getPreferredSize();
        Dimension frmSize = getSize();
        Point loc = getLocation();
        dlg.setLocation((frmSize.width - dlgSize.width) / 2 + loc.x, (frmSize.height - dlgSize.height) / 2 + loc.y);
        dlg.setModal(true);
        dlg.setVisible(true);
    }

    public void sendMessage(String text) {
        myStatusTextField.setForeground(Color.BLACK);
        myStatusTextField.setText(text);
    }

    public void sendErrorMessage(String text) {
        myStatusTextField.setForeground(Color.RED);
        myStatusTextField.setText(text);
    }

    public void setMainText(String text) {
        mainPanelTextArea.setText(text);
        mainPanelTextArea.setCaretPosition(0);
    }

    public void setNoteText(String text) {
        reportPanelTextArea.setText(text);
        reportPanelTextArea.setCaretPosition(0);
    }

    private void loadSettings() {
        filerOpen.setCurrentDirectory(new File(TasselPrefs.getOpenDir()));
        filerSave.setCurrentDirectory(new File(TasselPrefs.getSaveDir()));
    }

    /**
     * Provides a save filer that remembers the last location something was
     * saved to
     */
    private File getSaveFile() {

        File saveFile = null;
        int returnVal = filerSave.showSaveDialog(this);

        if (returnVal == JFileChooser.APPROVE_OPTION) {
            saveFile = filerSave.getSelectedFile();

            TasselPrefs.putSaveDir(filerSave.getCurrentDirectory().getPath());
        }
        return saveFile;
    }

    /**
     * Provides a open filer that remember the last location something was
     * opened from
     */
    private File getOpenFile() {

        File openFile = null;

        int returnVal = filerOpen.showOpenDialog(this);
        System.out.println("returnVal = " + returnVal);

        System.out.println("JFileChooser.OPEN_DIALOG " + JFileChooser.OPEN_DIALOG);

        if (returnVal == JFileChooser.OPEN_DIALOG || returnVal == JFileChooser.APPROVE_OPTION) {
            openFile = filerOpen.getSelectedFile();
            System.out.println("openFile = " + openFile);
            TasselPrefs.putOpenDir(filerOpen.getCurrentDirectory().getPath());
        }

        return openFile;
    }

    public void addDataSet(DataSet theDataSet, String defaultNode) {
        myDataTreePanel.addDataSet(theDataSet, defaultNode);
    }

    private void saveDataTree(String file) {

        Map dataToSerialize = new LinkedHashMap();
        Map dataFromTree = myDataTreePanel.getDataList();
        StringBuilder builder = new StringBuilder();

        Iterator itr = dataFromTree.keySet().iterator();
        while (itr.hasNext()) {

            Datum currentDatum = (Datum) itr.next();
            String currentNode = (String) dataFromTree.get(currentDatum);

            try {
                ByteArrayOutputStream out = new ByteArrayOutputStream();
                ObjectOutputStream oos = new ObjectOutputStream(out);
                oos.writeObject(currentDatum);
                oos.close();
                if (out.toByteArray().length > 0) {
                    dataToSerialize.put(currentDatum, currentNode);
                }
            } catch (Exception e) {
                myLogger.warn("saveDataTree: object: " + currentDatum.getName() + " type: " + currentDatum.getData().getClass().getName() + " does not serialize.");
                myLogger.warn("saveDataTree: message: " + e.getMessage());
                if (builder.length() == 0) {
                    builder.append("Due to error, these data sets could not\n");
                    builder.append("included in the saved file...");
                }
                builder.append("Data set: ");
                builder.append(currentDatum.getName());
                builder.append(" type: ");
                builder.append(currentDatum.getData().getClass().getName());
                builder.append("\n");
            }

        }

        try {

            File theFile = new File(Utils.addSuffixIfNeeded(file, ".zip"));
            FileOutputStream fos = new FileOutputStream(theFile);
            java.util.zip.ZipOutputStream zos = new ZipOutputStream(fos);

            ZipEntry thisEntry = new ZipEntry("DATA");
            zos.putNextEntry(thisEntry);
            ObjectOutputStream oos = new ObjectOutputStream(zos);
            oos.writeObject(dataToSerialize);
            oos.flush();
            zos.closeEntry();
            fos.close();
            sendMessage("Data saved to " + theFile.getAbsolutePath());

        } catch (Exception ee) {
            sendErrorMessage("Data could not be saved: " + ee);
            ee.printStackTrace();
        }

        if (builder.length() != 0) {
            JOptionPane.showMessageDialog(this, builder.toString(), "These data sets not saved...", JOptionPane.INFORMATION_MESSAGE);
        }

    }

    private boolean readDataTree(String file) {

        String dataTreeLoadFailed = "Unable to open the saved data tree.  The file format of this version is "
                + "incompatible with other versions.";

        boolean loadedDataTreePanel = false;
        try {

            FileInputStream fis = null;
            ObjectInputStream ois = null;
            if (file.endsWith("zip")) {

                fis = new FileInputStream(file);
                java.util.zip.ZipInputStream zis = new java.util.zip.ZipInputStream(fis);
                zis.getNextEntry();
                ois = new ObjectInputStream(zis);

            } else {

                fis = new FileInputStream(file);
                ois = new ObjectInputStream(fis);

            }

            try {
                Map data = (Map) ois.readObject();
                Iterator itr = data.keySet().iterator();
                while (itr.hasNext()) {
                    Datum currentDatum = (Datum) itr.next();
                    String currentNode = (String) data.get(currentDatum);
                    myDataTreePanel.addDatum(currentNode, currentDatum);
                }
                loadedDataTreePanel = true;
            } catch (InvalidClassException ice) {
                JOptionPane.showMessageDialog(this, dataTreeLoadFailed, "Incompatible File Format", JOptionPane.INFORMATION_MESSAGE);
            } finally {
                fis.close();
            }

            if (loadedDataTreePanel) {
                sendMessage("Data loaded.");
            }

        } catch (FileNotFoundException fnfe) {
            JOptionPane.showMessageDialog(this, "File not found: " + file, "File not found", JOptionPane.INFORMATION_MESSAGE);
            sendErrorMessage("Data tree could not be loaded.");
            return false;
        } catch (Exception ee) {
            JOptionPane.showMessageDialog(this, dataTreeLoadFailed + ee, "Incompatible File Format", JOptionPane.INFORMATION_MESSAGE);
            sendErrorMessage("Data tree could not be loaded.");
            return false;
        }

        return loadedDataTreePanel;

    }

    private void helpButton_actionPerformed(ActionEvent e) {
        HelpDialog theHelpDialog = new HelpDialog(this);
        theHelpDialog.setLocationRelativeTo(this);
        theHelpDialog.setVisible(true);
    }

    private void openCompleteDataTreeMenuItem_actionPerformed(ActionEvent e) {

        String dataFileName = tasselDataFile + ".zip";
        File dataFile = new File(dataFileName);
        if (dataFile.exists()) {
            readDataTree(dataFileName);
        } else if (new File("QPGADataFile").exists()) {
            // this exists to maintain backward compatibility with previous versions (pre-v0.99)
            readDataTree("QPGADataFile");
        } else {
            JOptionPane.showMessageDialog(this, "File: " + dataFile.getAbsolutePath() + " does not exist.\n"
                    + "Try using File/Open Data Tree...");
        }
    }

    private void openDataMenuItem_actionPerformed(ActionEvent e) {

        File f = getOpenFile();
        if (f != null) {
            readDataTree(f.getAbsolutePath());
        }
    }

    private void saveDataTreeMenuItem_actionPerformed(ActionEvent e) {

        File f = getSaveFile();
        if (f != null) {
            saveDataTree(f.getAbsolutePath());
        }
    }

    private void saveCompleteDataTreeMenuItem_actionPerformed(ActionEvent e) {
        saveDataTree(tasselDataFile + ".zip");
    }

    private void preferencesMenuItem_actionPerformed(ActionEvent e) {
        if (thePreferencesDialog == null) {
            thePreferencesDialog = new PreferencesDialog();
            thePreferencesDialog.pack();
        }

        thePreferencesDialog.setLocationRelativeTo(this);
        thePreferencesDialog.setVisible(true);
    }

    public void updateMainDisplayPanel(JPanel panel) {
        mainDisplayPanel.removeAll();
        mainDisplayPanel.add(panel, BorderLayout.CENTER);
        mainDisplayPanel.repaint();
        mainDisplayPanel.validate();
    }

    public DataTreePanel getDataTreePanel() {
        return myDataTreePanel;
    }

    public ProgressPanel getProgressPanel() {
        return myProgressPanel;
    }

    private JMenuItem createMenuItem(Plugin theTP) {
        return createMenuItem(theTP, -1);
    }

    private JMenuItem createMenuItem(Plugin theTP, boolean adjustForIcon) {
        return createMenuItem(theTP, -1, adjustForIcon);
    }

    private JMenuItem createMenuItem(Plugin theTP, int mnemonic) {
        return createMenuItem(theTP, mnemonic, true);
    }

    private JMenuItem createMenuItem(Plugin theTP, int mnemonic, boolean adjustForIcon) {
        ImageIcon icon = theTP.getIcon();
        JMenuItem menuItem = new JMenuItem(theTP.getButtonName(), icon);
        if (mnemonic != -1) {
            menuItem.setMnemonic(mnemonic);
        }
        if (adjustForIcon) {
            int pixels = 30;
            if (icon != null) {
                pixels -= icon.getIconWidth();
                pixels /= 2;
            }
            menuItem.setIconTextGap(pixels);
        }
        menuItem.setBackground(Color.white);
        menuItem.setMargin(new Insets(2, 2, 2, 2));
        menuItem.setToolTipText(theTP.getToolTipText());
        menuItem.addActionListener(this);
        theTP.addListener(myDataTreePanel);
        myMenuItemHash.put(menuItem, theTP);
        return menuItem;
    }

    private JMenuItem createMenuItem(Action action, int mnemonic) {

        Icon icon = (Icon) action.getValue("SmallIcon");

        JMenuItem menuItem = new JMenuItem(action.getValue("Name").toString(), icon);
        if (mnemonic != -1) {
            menuItem.setMnemonic(mnemonic);
        }
        int pixels = 30;
        if (icon != null) {
            pixels -= icon.getIconWidth();
            pixels /= 2;
        }
        menuItem.setIconTextGap(pixels);
        menuItem.setBackground(Color.white);
        menuItem.setMargin(new Insets(2, 2, 2, 2));
        menuItem.addActionListener(action);
        return menuItem;

    }

    private JMenuItem createMenuItem(String label, boolean adjustForIcon) {

        JMenuItem menuItem = new JMenuItem(label);
        if (adjustForIcon) {
            int pixels = 30;
            menuItem.setIconTextGap(pixels);
        }
        menuItem.setBackground(Color.white);
        menuItem.setMargin(new Insets(2, 2, 2, 2));
        menuItem.setEnabled(false);
        return menuItem;

    }

    private JMenu getFiltersMenu() {
        JMenu result = new JMenu("Filter");
        result.setMnemonic(KeyEvent.VK_F);
        result.add(createMenuItem(new FilterAlignmentPlugin(this, true)));
        result.add(createMenuItem(new FilterSiteNamePlugin(this, true)));
        result.add(createMenuItem(new FilterTaxaAlignmentPlugin(this, true)));
        result.add(createMenuItem(new FilterTaxaPropertiesPlugin(this, true)));
        result.add(createMenuItem(new FilterTraitsPlugin(this, true)));
        return result;
    }

    private JMenu getDataMenu() {

        JMenu result = new JMenu("Data");
        result.setMnemonic(KeyEvent.VK_D);

        PlinkLoadPlugin plinkLoadPlugin = new PlinkLoadPlugin(this, true);
        plinkLoadPlugin.addListener(myDataTreePanel);

        ProjectionLoadPlugin projectionLoadPlugin = new ProjectionLoadPlugin(this, true);
        projectionLoadPlugin.addListener(myDataTreePanel);

        ProjectPcsAndRunModelSelectionPlugin projectPcsAndRunModelSelectionPlugin
                = new ProjectPcsAndRunModelSelectionPlugin(this, true);
        projectPcsAndRunModelSelectionPlugin.addListener(myDataTreePanel);
        result.add(createMenuItem(new FileLoadPlugin(this, true, plinkLoadPlugin, projectionLoadPlugin, projectPcsAndRunModelSelectionPlugin), KeyEvent.VK_L));
        result.add(createMenuItem(new ExportPlugin(this, true)));
        result.add(createMenuItem(new GetTaxaListPlugin(this, true)));
        result.add(createMenuItem(new GetPositionListPlugin(this, true)));
        result.add(createMenuItem(new SortGenotypeFilePlugin(this, true)));
        result.add(createMenuItem(new SynonymizerPlugin(this, true)));
        result.add(createMenuItem(new IntersectionAlignmentPlugin(this, true)));
        result.add(createMenuItem(new UnionAlignmentPlugin(this, true)));
        result.add(createMenuItem(new MergeGenotypeTablesPlugin(this, true)));
        result.add(createMenuItem(new SeparatePlugin(this, true)));
        result.add(createMenuItem(new HetsToUnknownPlugin(this, true)));
        result.add(createMenuItem(new TransformDataPlugin(this, true)));
        result.add(createMenuItem(new NumericalGenotypePlugin(this, true)));
        result.addSeparator();

        JMenuItem delete = new JMenuItem("Delete Dataset");
        delete.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(ActionEvent e) {
                myDataTreePanel.deleteSelectedNodes();
            }
        });
        delete.setToolTipText("Delete Dataset");
        URL imageURL = TASSELMainFrame.class.getResource("images/trash.gif");
        if (imageURL == null) {
            delete.setIconTextGap(30);
        } else {
            delete.setIcon(new ImageIcon(imageURL));
            delete.setIconTextGap(6);
        }
        delete.setIconTextGap(6);
        result.add(delete);

        return result;
    }
    
    private JMenu getImputeMenu() {

        JMenu result = new JMenu("Impute");
        result.setMnemonic(KeyEvent.VK_I);

        result.add(createMenuItem(new FILLINFindHaplotypesPlugin(this, true)));
        result.add(createMenuItem(new FILLINImputationPlugin(this, true)));
        result.add(createMenuItem(new FSFHapImputationPlugin(this, true)));
        result.add(createMenuItem(new ImputationPlugin(this, true)));
        result.add(createMenuItem(new RemoveIndelsForBeaglePlugin(this, true)));
        return result;
    }

    private JMenu getAnalysisMenu() {

        JMenu result = new JMenu("Analysis");
        result.setMnemonic(KeyEvent.VK_A);

        result.add(createMenuItem(new SequenceDiversityPlugin(this, true)));
        result.add(createMenuItem(new LinkageDisequilibriumPlugin(this, true)));
        result.add(createMenuItem(new DistanceMatrixPlugin(this, true)));
        result.add(createMenuItem(new CreateTreePlugin(this, true)));
        result.add(createMenuItem(new KinshipPlugin(this, true)));
        result.add(createMenuItem(new PrincipalComponentsPlugin(this, true)));
        result.add(createMenuItem(new FixedEffectLMPlugin(this, true)));
        result.add(createMenuItem(new MLMPlugin(this, true)));
        result.add(createMenuItem(new RidgeRegressionEmmaPlugin(this, true)));
        result.add(createMenuItem(new GenotypeSummaryPlugin(this, true)));
        result.add(createMenuItem(new StepwiseOLSModelFitterPlugin(this, true)));
        result.add(createMenuItem(new EqtlAssociationPlugin(this, true)));
        return result;
    }

    private JMenu getResultsMenu() {

        JMenu result = new JMenu("Results");
        result.setMnemonic(KeyEvent.VK_R);

        result.add(createMenuItem(new TableDisplayPlugin(this, true)));
        result.add(createMenuItem(new ArchaeopteryxPlugin(this, true)));
        result.add(createMenuItem(new Grid2dDisplayPlugin(this, true)));
        result.add(createMenuItem(new LinkageDiseqDisplayPlugin(this, true)));
        result.add(createMenuItem(new ChartDisplayPlugin(this, true)));
        result.add(createMenuItem(new QQDisplayPlugin(this, true)));
        result.add(createMenuItem(new ManhattanDisplayPlugin(this, true)));
        return result;

    }

    private JMenu getFileMenu() {
        JMenu fileMenu = new JMenu();
        fileMenu.setText("File");
        fileMenu.add(saveCompleteDataTreeMenuItem);
        fileMenu.add(openCompleteDataTreeMenuItem);
        fileMenu.add(saveDataTreeAsMenuItem);
        fileMenu.add(openDataMenuItem);
        JMenuItem preferencesMenuItem = new JMenuItem();
        preferencesMenuItem.setText("Set Preferences");

        preferencesMenuItem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(ActionEvent e) {

                preferencesMenuItem_actionPerformed(e);

            }
        });
        fileMenu.add(preferencesMenuItem);
        fileMenu.addSeparator();

        JMenuItem exitMenuItem = new JMenuItem("Exit");
        exitMenuItem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(ActionEvent e) {
                TasselLogging.closeInstance();
                TasselPrefs.putXDim(getWidth());
                TasselPrefs.putYDim(getHeight());
                System.exit(0);
            }
        });
        fileMenu.add(exitMenuItem);

        return fileMenu;
    }

    private JMenu getGBSMenu() {

        JMenu result = new JMenu("GBS");
        result.setMnemonic(KeyEvent.VK_G);

        result.add(createMenuItem(new BinaryToTextPlugin(this, true), false));
        result.add(createMenuItem(new FastqToTagCountPlugin(this, true), false));
        result.add(createMenuItem(new MergeMultipleTagCountPlugin(this, true), false));
        result.addSeparator();
        result.add(getGBSReferenceMenu());
        result.addSeparator();
        result.add(getUNEAKMenu());
        result.addSeparator();
        result.add(createMenuItem(new SeqToTBTHDF5Plugin(this, true), false));
        result.add(createMenuItem(new ModifyTBTHDF5Plugin(this, true), false));
        result.add(createMenuItem(new DiscoverySNPCallerPlugin(this, true), false));
        result.add(createMenuItem(new ProductionSNPCallerPlugin(this, true), false));

        return result;
    }

    private JMenu getGBSReferenceMenu() {

        JMenu result = new JMenu("Reference Genome");

        result.add(createMenuItem(new TagCountToFastqPlugin(this, true), false));
        result.add(createMenuItem("Align to Reference", false));
        result.add(createMenuItem(new SAMConverterPlugin(this, true), false));

        return result;
    }

    private JMenu getUNEAKMenu() {

        JMenu result = new JMenu("UNEAK (No Reference)");

        result.add(createMenuItem(new UTagCountToTagPairPlugin(this, true)));
        result.add(createMenuItem(new UTagPairToTOPMPlugin(this, true)));

        return result;
    }

    private JMenu getHelpMenu() {
        JMenu helpMenu = new JMenu();
        helpMenu.setMnemonic(KeyEvent.VK_H);
        helpMenu.setText("Help");

        URL infoImageURL = TASSELMainFrame.class.getResource("images/info.gif");
        Icon infoIcon = null;
        if (infoImageURL != null) {
            infoIcon = new ImageIcon(infoImageURL);
        }
        helpMenu.add(createMenuItem(new AbstractAction("Help Manual", infoIcon) {

            @Override
            public void actionPerformed(ActionEvent e) {
                helpButton_actionPerformed(e);
            }
        }, -1));

        URL aboutImageURL = TASSELMainFrame.class.getResource("images/Tassel_Logo16.png");
        Icon aboutIcon = null;
        if (aboutImageURL != null) {
            aboutIcon = new ImageIcon(aboutImageURL);
        }
        helpMenu.add(createMenuItem(new AbstractAction("About", aboutIcon) {

            @Override
            public void actionPerformed(ActionEvent e) {
                helpAbout_actionPerformed(e);
            }
        }, -1));

        helpMenu.add(createMenuItem(PrintHeapAction.getInstance(this), -1));
        helpMenu.add(createMenuItem(TasselLogging.getInstance(this)));

        return helpMenu;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        JMenuItem theMenuItem = (JMenuItem) e.getSource();
        Plugin theTP = this.myMenuItemHash.get(theMenuItem);
        PluginEvent event = new PluginEvent(myDataTreePanel.getSelectedTasselDataSet());
        ProgressPanel progressPanel = getProgressPanel();
        progressPanel.addPlugin(theTP);
        ThreadedPluginListener thread = new ThreadedPluginListener(theTP, event);
        thread.start();
    }
}
