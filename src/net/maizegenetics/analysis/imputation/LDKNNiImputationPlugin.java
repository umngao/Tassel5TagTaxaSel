package net.maizegenetics.analysis.imputation;

import net.maizegenetics.dna.snp.*;
import net.maizegenetics.dna.snp.io.ProjectionGenotypeIO;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.util.OpenBitSet;

import javax.swing.*;
import java.awt.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

/**
 * Created by edbuckler on 5/4/15.
 */
public class LDKNNiImputationPlugin extends AbstractPlugin{


    public LDKNNiImputationPlugin() {
        super(null, false);
    }

    public LDKNNiImputationPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public DataSet processData(DataSet input) {
        long time=System.currentTimeMillis();
        System.out.println(time);
        return null;
    }

    @Override
    public String getCitation() {
        return "";
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return null;
    }

    @Override
    public String getToolTipText() {
        return null;
    }
}
