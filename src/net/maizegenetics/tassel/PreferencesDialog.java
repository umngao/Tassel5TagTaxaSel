package net.maizegenetics.tassel;

import java.awt.Frame;
import java.net.URL;
import javax.swing.*;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.prefs.TasselPrefs;

/**
 * @author Terry Casstevens
 */
public class PreferencesDialog extends AbstractPlugin {

    private PluginParameter<Boolean> myRetainRareAlleles = new PluginParameter.Builder<>("retainRareAlleles", TasselPrefs.ALIGNMENT_RETAIN_RARE_ALLELES_DEFAULT, Boolean.class)
            .description("True if rare alleles should be retained.  This has no effect on Nucleotide Data as all alleles will be retained regardless.")
            .build();

    private PluginParameter<Boolean> mySendLogToConsole = new PluginParameter.Builder<>("sendLogToConsole", TasselPrefs.TASSEL_LOG_SEND_TO_CONSOLE_DEFAULT, Boolean.class)
            .description("Flag whether to send logging to the console.")
            .build();

    public PreferencesDialog(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    protected void preProcessParameters(DataSet input) {
        setParameter(myRetainRareAlleles, TasselPrefs.getAlignmentRetainRareAlleles());
        setParameter(mySendLogToConsole, TasselPrefs.getLogSendToConsole());
    }

    @Override
    public DataSet processData(DataSet input) {
        TasselPrefs.putAlignmentRetainRareAlleles(retainRareAlleles());
        TasselPrefs.putLogSendToConsole(sendLogToConsole());
        TasselLogging.updateLoggingLocation();
        return null;
    }

    /**
     * Retain Rare Alleles
     *
     * @return Retain Rare Alleles
     */
    public Boolean retainRareAlleles() {
        return myRetainRareAlleles.value();
    }

    /**
     * Set Retain Rare Alleles. Retain Rare Alleles
     *
     * @param value Retain Rare Alleles
     *
     * @return this plugin
     */
    public PreferencesDialog retainRareAlleles(Boolean value) {
        myRetainRareAlleles = new PluginParameter<>(myRetainRareAlleles, value);
        return this;
    }

    /**
     * Flag whether to send logging to the console.
     *
     * @return Send Log To Console
     */
    public Boolean sendLogToConsole() {
        return mySendLogToConsole.value();
    }

    /**
     * Set Send Log To Console. Flag whether to send logging to the console.
     *
     * @param value Send Log To Console
     *
     * @return this plugin
     */
    public PreferencesDialog sendLogToConsole(Boolean value) {
        mySendLogToConsole = new PluginParameter<>(mySendLogToConsole, value);
        return this;
    }

    @Override
    public ImageIcon getIcon() {
        URL imageURL = TasselLogging.class.getResource("/net/maizegenetics/analysis/images/preferences.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    @Override
    public String getButtonName() {
        return "Preferences";
    }

    @Override
    public String getToolTipText() {
        return "Preferences";
    }
}
