package net.maizegenetics.analysis.modelfitter;

public abstract class AbstractAdditiveSite implements AdditiveSite {
    protected final int siteIndex;
    protected double criterionValue;
    protected final CRITERION selectionCriterion;
    protected final int direction;
    
    public AbstractAdditiveSite(int site, CRITERION selectionCriterion) {
        siteIndex = site;
        this.selectionCriterion = selectionCriterion;
        if (selectionCriterion == CRITERION.ss) direction = 1;
        else direction = -1;
    }

    @Override
    public int siteNumber() {
        return siteIndex;
    }

    @Override
    public double criterionValue() {
        return criterionValue;
    }

    @Override
    public void criterionValue(double value) {
        criterionValue = value;
    }

    @Override
    public int compareTo(AdditiveSite other) {
        return direction * Double.compare(criterionValue, other.criterionValue());
    }
    
}
