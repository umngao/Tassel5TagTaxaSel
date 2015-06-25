package net.maizegenetics.analysis.modelfitter;

import java.util.List;

public interface ForwardRegression {

    /**
     * @return  a list of terms in the model. Each item in the list is an an array whose first element is a Position. The second element an Integer equal to the site index of the term.
     * The third element is a Double equal to the p-value for the term.
     */
    List<Object[]> fittedModel();

    void fitModelForSubsample(int[] subSample);

    void fitModel();

    public static String[] columnLabels() {
        return new String[]{"trait","SnpID","Chr","Pos", "p-value", "-log10p"};
    }
}