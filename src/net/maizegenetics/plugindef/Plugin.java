/*
 * Plugin.java
 *
 */
package net.maizegenetics.plugindef;

import net.maizegenetics.util.ProgressListener;

import javax.swing.*;
import java.awt.*;

/**
 *
 * @author terryc
 */
public interface Plugin extends PluginListener, ProgressListener, Runnable {

    /**
     * Returns menu that can be added to main menu bar.
     *
     * @return menu
     */
    public JMenu getMenu();

    /**
     * Performs function of this plugin.
     *
     * @param input input
     *
     * @return resulting data set or null.
     */
    public DataSet performFunction(DataSet input);

    /**
     * For the new Generic Plugin Parameter design, performFunction() will
     * automatically call this. Therefore, coders of Plugins should override
     * this instead of performFunction().
     *
     * @param input input
     *
     * @return resulting data set or null.
     */
    public DataSet processData(DataSet input);

    /**
     * Returns parameter value for given parameter key.
     *
     * @param key key
     *
     * @return value
     */
    public Comparable getParameter(Enum key);
    
    /**
     * Returns parameter value for given parameter key.
     *
     * @param key key
     *
     * @return value
     */
    public Comparable getParameter(String key);

    /**
     * Sets parameter value for a given plugin
     *
     * @param param parameter
     * @param value value
     *
     * @return this plugin
     */
    public Plugin setParameter(PluginParameter param, Object value);

    /**
     * Sets parameter value for given parameter key.
     *
     * @param key key
     * @param value value
     *
     * @return this plugin
     */
    public Plugin setParameter(Enum key, Comparable value);

    /**
     * Sets parameter value for given parameter key.
     *
     * @param key key
     * @param value value
     *
     * @return this plugin
     */
    public Plugin setParameter(Enum key, String value);
    
    /**
     * Sets parameter value for given parameter key.
     *
     * @param key key
     * @param value value
     *
     * @return this plugin
     */
    public Plugin setParameter(String key, Comparable value);

    /**
     * Sets parameter value for given parameter key.
     *
     * @param key key
     * @param value value
     *
     * @return this plugin
     */
    public Plugin setParameter(String key, String value);

    /**
     * Sets up this plugin to receive input from another plugin.
     *
     * @param input input
     */
    public void receiveInput(Plugin input);

    /**
     * GUI Panel for this plugin.
     *
     * @return panel
     */
    public JPanel getPanel();

    /**
     * If interactive = true, the plugin will create dialogs and panels to
     * interacts with the user
     *
     * @return boolean
     */
    public boolean isInteractive();

    /**
     * Parent Frame for this plugin. Can be null.
     *
     * @return frame
     */
    public Frame getParentFrame();

    /**
     * Icon for this plugin to be used in buttons, etc.
     *
     * @return ImageIcon
     */
    public ImageIcon getIcon();

    /**
     * Button name for this plugin to be used in buttons, etc.
     *
     * @return String
     */
    public String getButtonName();

    /**
     * Tool Tip Text for this plugin
     *
     * @return String
     */
    public String getToolTipText();

    /**
     * Adds listener to this plugin.
     *
     * @param listener listener to add
     */
    public void addListener(PluginListener listener);

    /**
     * Set whether this plugin is threaded.
     *
     * @param threaded whether to be threaded.
     */
    public void setThreaded(boolean threaded);

    /**
     * Attempt to cancel processing.
     *
     * @return true if plugin will cancel itself.
     */
    public boolean cancel();

    /**
     * Allows self-describing Plugins to use args to set parameters specific to
     * itself.
     *
     * @param args arguments
     */
    public void setParameters(String[] args);

    /**
     * Returns Citation for this plugin.
     *
     * @return Citation
     */
    public String getCitation();
    
    /**
     * Returns description of the plugin.
     * 
     * @return description
     */
    public String pluginDescription();
}
