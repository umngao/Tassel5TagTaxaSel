package net.maizegenetics.analysis.data;

import java.awt.Frame;

import javax.swing.ImageIcon;
import javax.swing.JDialog;

import org.apache.log4j.Logger;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.PluginParameter;

public class PhenotypeTransformPlugin extends AbstractPlugin {
	private static final Logger myLogger = Logger.getLogger(PhenotypeTransformPlugin.class);
	public enum TRANSFORM_TYPE {power, log, standardize};
	
    public PhenotypeTransformPlugin(Frame parentFrame, boolean isInteractive){
        super(parentFrame, isInteractive);
    }

	@Override
	public ImageIcon getIcon() {
		return null;
	}

	@Override
	public String getButtonName() {
		return "Transform";
	}

	@Override
	public String getToolTipText() {
		return "Transform trait data to correct for lack of normality";
	}

	class TransformationDialog extends JDialog {
		
		TransformationDialog(Frame parentFrame) {
			super(parentFrame);
			
		}
	}
}


