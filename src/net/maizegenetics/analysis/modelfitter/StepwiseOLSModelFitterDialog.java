package net.maizegenetics.analysis.modelfitter;

import java.awt.Color;
import java.awt.Component;
import java.awt.Container;
import java.awt.Frame;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JRootPane;
import javax.swing.JScrollPane;
import javax.swing.JTextField;


public class StepwiseOLSModelFitterDialog extends JDialog implements ActionListener {
	JCheckBox chkNested = null;
	JList<String> listMainEffects = null;
	JButton btnLimits = new JButton("Enter limits by trait");
	JButton btnOk = new JButton("OK");
	JButton btnCancel = new JButton("Cancel");
	JLabel lblWhich = new JLabel("Which factor?");
	JTextField txtEnter, txtExit, txtMax, txtPerm, txtAlpha;
	String enterlim, exitlim;
	String[] mainEffects = null;
	int indexOfSelectedEffect = 0;
	boolean wasCancelled = true;
	int numberOfMainEffects = 0;

	ButtonGroup modelSelectionButtonGroup = new ButtonGroup();
	JRadioButton PvalueRadioButton = new JRadioButton("Use P-value thresholds");
	JRadioButton mBICRadioButton = new JRadioButton("mBIC model selection");
	JRadioButton BICRadioButton = new JRadioButton("BIC model selection");
	JRadioButton AICRadioButton = new JRadioButton("AIC model selection");

	
	public StepwiseOLSModelFitterDialog(String[] mainEffects, Frame parentFrame) {
		super();
		enterlim = "1e-5";
		exitlim = "2e-5";
		this.mainEffects = mainEffects;
		if (mainEffects != null) numberOfMainEffects = mainEffects.length;
		this.setTitle("Specify Model Parameters");
		this.setModalityType(ModalityType.APPLICATION_MODAL);
        setDefaultCloseOperation(JDialog.HIDE_ON_CLOSE);
        setUndecorated(false);
        getRootPane().setWindowDecorationStyle(JRootPane.NONE);
        Container contentPane = getContentPane();
        JPanel myPanel = new JPanel();
        myPanel.setLayout(new BoxLayout(myPanel, BoxLayout.Y_AXIS));
        myPanel.setBorder(BorderFactory.createEmptyBorder(20, 20, 20, 20));
        JScrollPane myScrollPane = new JScrollPane(myPanel);
       
//        contentPane.add(myPanel);
        contentPane.add(myScrollPane);
        
        modelSelectionButtonGroup.add(PvalueRadioButton);
        modelSelectionButtonGroup.add(mBICRadioButton);
        modelSelectionButtonGroup.add(BICRadioButton);
        modelSelectionButtonGroup.add(AICRadioButton);
        PvalueRadioButton.setSelected(true);
        
            if (numberOfMainEffects > 1) {
                    chkNested = new JCheckBox("Nest markers within a factor", null, false);
                    myPanel.add(chkNested);
                    chkNested.setAlignmentX(CENTER_ALIGNMENT);
                    chkNested.addActionListener(this);
                    chkNested.setActionCommand("nest");


                    lblWhich.setAlignmentX(CENTER_ALIGNMENT);
                    lblWhich.setEnabled(false);
                    myPanel.add(Box.createVerticalStrut(30));
                    myPanel.add(lblWhich);

                    listMainEffects = new JList<String>(mainEffects);
                    listMainEffects.setBorder(BorderFactory.createLineBorder(Color.BLACK));
                    myPanel.add(Box.createVerticalStrut(10));
                    myPanel.add(listMainEffects);
                    listMainEffects.setAlignmentX(Component.CENTER_ALIGNMENT);
                    listMainEffects.setEnabled(false);

            } else if (numberOfMainEffects == 1) {
                    chkNested = new JCheckBox("Nest markers within " + mainEffects[0], null, false);
                    myPanel.add(chkNested);
                    chkNested.setAlignmentX(CENTER_ALIGNMENT);
                    chkNested.addActionListener(this);
                    chkNested.setActionCommand("nest");
            }

            if (numberOfMainEffects > 0) {
                    myPanel.add(Box.createVerticalStrut(20));
                    JLabel sep1 = new JLabel("-----------------------");
                    sep1.setAlignmentX(CENTER_ALIGNMENT);
                    myPanel.add(sep1);
            }

            myPanel.add(Box.createVerticalStrut(20));

            Box enterBox = Box.createHorizontalBox();
            enterBox.add(new JLabel("enter limit"));
            txtEnter = new JTextField(enterlim, 15);
            enterBox.add(txtEnter);
            myPanel.add(enterBox);

            myPanel.add(Box.createVerticalStrut(10));

            Box exitBox = Box.createHorizontalBox();
            exitBox.add(new JLabel("exit limit"));
            txtExit = new JTextField(exitlim, 15);
            exitBox.add(txtExit);
            myPanel.add(exitBox);

            myPanel.add(Box.createVerticalStrut(10));
            myPanel.add(btnLimits);
            btnLimits.setActionCommand("limits");
            btnLimits.setAlignmentX(CENTER_ALIGNMENT);


            myPanel.add(Box.createVerticalStrut(20));
            JLabel sep2 = new JLabel("-----------------------");
            myPanel.add(sep2);
            sep2.setAlignmentX(CENTER_ALIGNMENT);
            myPanel.add(Box.createVerticalStrut(20));

            Box maxBox = Box.createHorizontalBox();
            maxBox.add(new JLabel("Maximum number of markers"));
            txtMax = new JTextField("30", 6);
            maxBox.add(txtMax);
            myPanel.add(maxBox);

            myPanel.add(Box.createVerticalStrut(20));
    JLabel sep3 = new JLabel("-----------------------");
        myPanel.add(sep3);
        sep3.setAlignmentX(CENTER_ALIGNMENT);  
		myPanel.add(Box.createVerticalStrut(20));

                Box permBox = Box.createHorizontalBox();
		permBox.add(new JLabel("Number of Permutations (only for P-value model selection)"));
		txtPerm = new JTextField("0", 6);
		permBox.add(txtPerm);
		myPanel.add(permBox);
                myPanel.add(Box.createVerticalStrut(20));
 
		
                myPanel.add(Box.createVerticalStrut(20));
		JLabel sep4 = new JLabel("------");
		myPanel.add(sep4);
		sep4.setAlignmentX(CENTER_ALIGNMENT);
		myPanel.add(Box.createVerticalStrut(20));
                
                Box alphaBox = Box.createHorizontalBox();
		alphaBox.add(new JLabel("Type I Error Rate for Permutations"));
		txtAlpha = new JTextField("0.05", 6);
		alphaBox.add(txtAlpha);
		myPanel.add(alphaBox);
		
		myPanel.add(Box.createVerticalStrut(20));
                myPanel.add(Box.createVerticalStrut(20));
                JLabel sep5 = new JLabel("-----------------------");
                myPanel.add(sep5);
                sep5.setAlignmentX(CENTER_ALIGNMENT);  

        myPanel.add(PvalueRadioButton);
        myPanel.add(mBICRadioButton);
        myPanel.add(BICRadioButton);
        myPanel.add(AICRadioButton);

        PvalueRadioButton.setAlignmentX(CENTER_ALIGNMENT);
        mBICRadioButton.setAlignmentX(CENTER_ALIGNMENT);
        BICRadioButton.setAlignmentX(CENTER_ALIGNMENT);
        AICRadioButton.setAlignmentX(CENTER_ALIGNMENT);

        myPanel.add(Box.createVerticalStrut(20));
                JLabel sep6 = new JLabel("-----------------------");
                myPanel.add(sep6);
		sep6.setAlignmentX(CENTER_ALIGNMENT);
		myPanel.add(Box.createVerticalStrut(20));
		
		Box buttonBox = Box.createHorizontalBox();
		buttonBox.add(btnOk);
		btnOk.setAlignmentX(CENTER_ALIGNMENT);
		btnOk.addActionListener(this);
		btnOk.setActionCommand("ok");
		buttonBox.add(btnCancel);
		btnCancel.setAlignmentX(CENTER_ALIGNMENT);
		btnCancel.addActionListener(this);
		btnCancel.setActionCommand("cancel");
		myPanel.add(buttonBox);
		buttonBox.setAlignmentX(CENTER_ALIGNMENT);

		myPanel.add(Box.createVerticalStrut(20));
		JLabel sep7 = new JLabel("------Developed by Peter Bradbury and Alex Lipka, 2013-----------------");
		myPanel.add(sep7);
		sep7.setAlignmentX(CENTER_ALIGNMENT);
		JLabel note1 = new JLabel("Note: Because the Stepwise method does not allow missing marker data,");
		JLabel note2 = new JLabel("the marker value N will be treated as a marker class, not missing data.");
		note1.setAlignmentX(CENTER_ALIGNMENT);
		note2.setAlignmentX(CENTER_ALIGNMENT);
		myPanel.add(note1);
		myPanel.add(note2);
                
		pack();
		setVisible(true);
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		System.out.println("Action Command: " + e.getActionCommand());
		if (e.getActionCommand().equals("nest") && numberOfMainEffects > 1) {
			if (chkNested.isSelected() ) {
				listMainEffects.setEnabled(true);
				lblWhich.setEnabled(true);
			} else {
				listMainEffects.setEnabled(false);
				lblWhich.setEnabled(false);
			}
		} else if (e.getActionCommand().equals("ok")) {
			wasCancelled = false;
			setVisible(false);
		} else if (e.getActionCommand().equals("cancel")) {
			wasCancelled = true;
			setVisible(false);
		}
	
	}
	
	public boolean isNested() {
		if (chkNested == null) return false;
		return chkNested.isSelected();
	}
	
	public String getNestedEffect() {
		if (mainEffects == null || mainEffects.length == 0) return null;
		if (mainEffects.length == 1) return mainEffects[0];
		return listMainEffects.getSelectedValue();
	}
	
	public boolean wasCancelled() {
		return wasCancelled;
	}
	
	public double[] getEnterLimits() {
		return StepwiseOLSModelFitterPlugin.parseDoubles(txtEnter.getText(), ",");
	}
	
	public double[] getExitLimits() {
		return StepwiseOLSModelFitterPlugin.parseDoubles(txtExit.getText(), ",");
	}
	
        public int getNumberOfPermutations() {
		return Integer.parseInt(txtPerm.getText());
	}

        public double getAlpha() {
		return Double.parseDouble(txtAlpha.getText());
	}
	
	public int getMaxNumberOfMarkers() {
		return Integer.parseInt(txtMax.getText());
	}
        

	public StepwiseOLSModelFitter.MODEL_TYPE getModelType(){
		if (PvalueRadioButton.isSelected()) {
			return StepwiseOLSModelFitter.MODEL_TYPE.pvalue;
		}
		if (mBICRadioButton.isSelected()) {
			return StepwiseOLSModelFitter.MODEL_TYPE.mbic;
		}
		if (BICRadioButton.isSelected()) {
			return StepwiseOLSModelFitter.MODEL_TYPE.bic;
		}
		if (AICRadioButton.isSelected()) {
			return StepwiseOLSModelFitter.MODEL_TYPE.aic;
		}
		return StepwiseOLSModelFitter.MODEL_TYPE.pvalue;

	}

        
}
