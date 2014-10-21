package net.maizegenetics.analysis.data;

import java.awt.Frame;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.util.ArrayList;
import java.util.List;

import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JScrollBar;
import javax.swing.JScrollPane;
import javax.swing.JTextField;

import net.maizegenetics.phenotype.CategoricalAttribute;
import net.maizegenetics.phenotype.NumericAttribute;
import net.maizegenetics.util.OpenBitSet;

public class TransformDataDialog extends JDialog {
	private final List<NumericAttribute> dataAttributes;
	private final List<CategoricalAttribute> factorAttributes;
	private JComboBox<NumericAttribute> cboTraits;
	private JComboBox<CategoricalAttribute> cboFactors;
	private JCheckBox chkPower;
	private JCheckBox chkLog;
	private JRadioButton radNatural;
	private JRadioButton radBase2;
	private JRadioButton radBase10;
	private JTextField txtPower;
	private JCheckBox chkStandard;
	private JButton btnOk;
	private JButton btnCancel;
	
	public TransformDataDialog(Frame parentFrame, List<NumericAttribute> dataAttributes, List<CategoricalAttribute> factorAttributes) {
		super(parentFrame);
		this.dataAttributes = dataAttributes;
		this.factorAttributes = factorAttributes;
		init();
	}
	
	public static void main(String[] args) {
		List<NumericAttribute> dataAttributes = new ArrayList<>();
		dataAttributes.add(new NumericAttribute("trait1", new float[]{1f,2f,3f}, new OpenBitSet(3)));
		dataAttributes.add(new NumericAttribute("trait2", new float[]{1f,2f,3f}, new OpenBitSet(3)));
		dataAttributes.add(new NumericAttribute("trait3", new float[]{1f,2f,3f}, new OpenBitSet(3)));
		dataAttributes.add(new NumericAttribute("trait4", new float[]{1f,2f,3f}, new OpenBitSet(3)));
		dataAttributes.add(new NumericAttribute("trait5", new float[]{1f,2f,3f}, new OpenBitSet(3)));
		dataAttributes.add(new NumericAttribute("trait6", new float[]{1f,2f,3f}, new OpenBitSet(3)));
		dataAttributes.add(new NumericAttribute("trait7", new float[]{1f,2f,3f}, new OpenBitSet(3)));
		
		List<CategoricalAttribute> factorAttributes = new ArrayList<>();
		factorAttributes.add(new CategoricalAttribute("factor1", new String[]{"one","one","two"}));
		
		TransformDataDialog tdd = new TransformDataDialog(null, dataAttributes, factorAttributes);
		tdd.setVisible(true);
		System.exit(0);
	}
	
	private void init() {
		setTitle("Transform Data");
		setModalityType(ModalityType.APPLICATION_MODAL);
		
		cboTraits = new JComboBox<>(dataAttributes.toArray(new NumericAttribute[0]));
		cboFactors = new JComboBox<>(factorAttributes.toArray(new CategoricalAttribute[0]));
		
		chkPower = new JCheckBox("Power", false);
		txtPower = new JTextField();
		
		chkLog = new JCheckBox("Log", false);
		radNatural = new JRadioButton("Natural", true);
		radBase2 = new JRadioButton("Base 2", false);
		radBase10 = new JRadioButton("Base 10", false);
		ButtonGroup logBaseGroup = new ButtonGroup();
		logBaseGroup.add(radNatural);
		logBaseGroup.add(radBase2);
		logBaseGroup.add(radBase10);
		
		chkStandard = new JCheckBox("Standardize");
		btnOk = new JButton("OK");
		btnCancel = new JButton("Cancel");
		
		
		JPanel mainPanel = new JPanel(new GridBagLayout());
		
		GridBagConstraints gbc = new GridBagConstraints();
		gbc.gridx = 0;
		gbc.gridy = 0;
		gbc.gridwidth = 2;
		mainPanel.add(new JLabel("Select Traits:"), gbc);
		
		gbc.gridy++;
		JScrollPane traitScroller = new JScrollPane(cboTraits);
		mainPanel.add(traitScroller, gbc);
		
		gbc.gridy++;
		mainPanel.add(chkPower,gbc);
		
		gbc.gridy++;
		mainPanel.add(txtPower,gbc);
		
		gbc.gridy++;
		mainPanel.add(chkLog,gbc);
		
		gbc.gridy++;
		mainPanel.add(radNatural,gbc);

		gbc.gridy++;
		mainPanel.add(radBase2,gbc);

		gbc.gridy++;
		mainPanel.add(radBase10,gbc);
		
		gbc.gridy++;
		mainPanel.add(chkStandard,gbc);
		
		gbc.gridy++;
		gbc.gridwidth = 1;
		mainPanel.add(btnOk,gbc);
		
		gbc.gridx++;
		mainPanel.add(btnCancel,gbc);
		
		getContentPane().add(mainPanel);
		pack();
}
	
}
