package net.maizegenetics.analysis.data;

import java.awt.Frame;

import javax.swing.ImageIcon;

import org.apache.log4j.Logger;

import net.maizegenetics.analysis.association.FixedEffectLMPlugin;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.util.LoggingUtils;

public class PrincipalComponentsPlugin extends AbstractPlugin {
    private static final Logger myLogger = Logger.getLogger(FixedEffectLMPlugin.class);

	public PrincipalComponentsPlugin(Frame parentFrame, boolean isInteractive) {
		super(parentFrame, isInteractive);
	}
	
	public DataSet processData(DataSet input){
		//TODO implement
		return null;
	}
	
	@Override
	public ImageIcon getIcon() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public String getButtonName() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public String getToolTipText() {
		// TODO Auto-generated method stub
		return null;
	}

}
