package net.maizegenetics.analysis.distance;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.dna.snp.score.ReferenceProbability;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrixFactory;
import net.maizegenetics.phenotype.NumericAttribute;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.phenotype.PhenotypeAttribute;
import net.maizegenetics.phenotype.TaxaAttribute;
import net.maizegenetics.phenotype.Phenotype.ATTRIBUTE_TYPE;
import net.maizegenetics.dna.snp.GenotypeTable.GENOTYPE_TABLE_COMPONENT;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.distance.DistanceMatrix;

import java.awt.Frame;
import java.util.List;
import java.util.function.BiConsumer;
import java.util.stream.IntStream;

/**
 * Created by IntelliJ IDEA. User: Zhiwu Date: Apr 29, 2007 Time: 4:21:44 PM To
 * change this template use File | Settings | File Templates.
 */
public class Kinship extends DistanceMatrix {

    Frame parentFrame;
    DistanceMatrix dm = null;
    GenotypeTable mar;
    Phenotype ped;
    int[][] parents;
    private double kMin = 99999;
    private double kMax = -99999;
    private double kAvg = 0;
    private double kSD = 0;
    private double cutOff = 2;
    private int numSeqs;
    private KINSHIP_TYPE kinshipType = KINSHIP_TYPE.IBS;
    private GENOTYPE_TABLE_COMPONENT myDataType;
    private double[][] distance;
    public static double matrixMultiplier = 2; //scale the numeric matrix produced by the transform function or from probabilities which code phenotypes as {1,0.5,0}

    public enum KINSHIP_TYPE {

        Endelman, IBS
    };

    public Kinship(GenotypeTable mar) {
        this(mar, false, true);
    }

    public Kinship(GenotypeTable mar, boolean areHetsRelated, boolean rescaleKinship) {
        this.mar = mar;
        numSeqs = this.mar.numberOfTaxa();
        buildFromMarker();
    }

    public Kinship(GenotypeTable mar, KINSHIP_TYPE kinshipType, GENOTYPE_TABLE_COMPONENT dataType) {
        this.mar = mar;
        this.kinshipType = kinshipType;
        myDataType = dataType;
        numSeqs = this.mar.numberOfTaxa();
        System.out.println("Starting Kinship.buildFromMarker.");
        long start = System.currentTimeMillis();
        buildFromMarker();
        System.out.printf("Built Kinship in %d millisec.\n", System.currentTimeMillis() - start);
    }

    public Kinship(Phenotype ped) {
        this.ped = ped;
        buildFromPed();
    }

    public Kinship(DistanceMatrix dm) {
        this.dm = dm;
    }

    public void buildFromMarker() {
        if (myDataType == GENOTYPE_TABLE_COMPONENT.Genotype) {
            if (kinshipType == KINSHIP_TYPE.Endelman) {
                calculateKinshipFromMarkers();
            } else {
                IBSDistanceMatrix adm = new IBSDistanceMatrix(mar, 0, true, null, true);
                dm = new DistanceMatrix(adm.getDistances(), mar.taxa());
                toSimilarity();
                getKStatistics();
            	//pullBackExtrem();
                //cutOff();
                rescale();
                System.out.println("Kinship was built from markers");
            }
        } else if (myDataType == GENOTYPE_TABLE_COMPONENT.ReferenceProbability) {
            calculateRelationshipKinshipFromReferenceProbability();
        } else {
            throw new IllegalArgumentException("The supplied data type is not currently supported by the Kinship method.");
        }

    }

    public void buildFromPed() {
        // get data from ped (Phenotype) to parents (int[][]);

        System.out.println("Building Kinship From pedigree");
        TaxaAttribute myTaxaAttribute = ped.taxaAttribute();
        int ntaxa = myTaxaAttribute.size();
        List<PhenotypeAttribute> dataAttributeList = ped.attributeListOfType(ATTRIBUTE_TYPE.data);
        int ntraits = dataAttributeList.size();

        parents = new int[ntaxa][ntraits];
        try {
            for (int row = 0; row < ntaxa; row++) {
                int attrNumber = 0;
                for (PhenotypeAttribute attr : dataAttributeList) {
                    NumericAttribute na = (NumericAttribute) attr;
                    parents[row][attrNumber++] = (int) na.floatValue(row);
                }
            }
        } catch (NumberFormatException e) {
            e.printStackTrace();
        }

        dm = new DistanceMatrix(kinshipRelation(parents), new TaxaListBuilder().addAll(myTaxaAttribute.allTaxaAsList()).build());

        System.out.println("Kinship was built from pedigree");
    }

    public static double[][] kinshipRelation(int[][] ped) {
        int n = ped.length;
        int maleParent;
        int femaleParent;
        double[][] aMatrix = new double[n][n];
        //System.out.println("size of ped: "+n);

        //initial: diagonal 1, 0 otherwise;
        for (int i = 0; i < n; i++) {
            aMatrix[i][i] = 1;
            for (int j = i + 1; j < n; j++) {
                aMatrix[i][j] = 0;
                aMatrix[j][i] = 0;
            }
        }

        System.out.println("initial: diagonal 1, 0 otherwise");
        for (int i = 0; i < n; i++) {
            //diagonals
            femaleParent = ped[i][1];
            maleParent = ped[i][2];
            if ((femaleParent > 0) && (maleParent > 0)) {
                aMatrix[i][i] = aMatrix[i][i] + .5 * aMatrix[maleParent - 1][femaleParent - 1];
            }
            //Off diagonals
            for (int j = i + 1; j < n; j++) {
                femaleParent = ped[j][1];
                maleParent = ped[j][2];

                if ((femaleParent > 0) && (maleParent > 0)) {
                    aMatrix[i][j] = .5 * (aMatrix[i][femaleParent - 1] + aMatrix[i][maleParent - 1]);
                } else if (maleParent > 0) {
                    aMatrix[i][j] = .5 * aMatrix[i][maleParent - 1];
                } else if (femaleParent > 0) {
                    aMatrix[i][j] = .5 * aMatrix[i][femaleParent - 1];
                } else {
                    //do nothing
                }
                aMatrix[j][i] = aMatrix[i][j];
            }
        }
        System.out.println("A matrix finished");

        return aMatrix;
    }

    //Convert distance to similarity
    //By Zhiwu Zhang
    public void toSimilarity() {
        double s;
        System.out.println("toSimilarity " + numSeqs);

        for (int i = 0; i < numSeqs; i++) {
            for (int j = i; j < numSeqs; j++) {
                s = cutOff - dm.getDistance(i, j);
                dm.setDistance(i, j, s);
                dm.setDistance(j, i, s);
            }
        }
        System.out.println("toSimilarity finish" + numSeqs);
    }

    public void getKStatistics() {
        //get average
        double total = 0;
        double totalsq = 0;
        double nk = numSeqs * (numSeqs - 1) / 2;
        for (int i = 0; i < numSeqs - 1; i++) {
            for (int j = i + 1; j < numSeqs; j++) {
                total += dm.getDistance(i, j);
                totalsq += (dm.getDistance(i, j) * dm.getDistance(i, j));
                if (dm.getDistance(i, j) < kMin) {
                    kMin = dm.getDistance(i, j);
                }
                if (dm.getDistance(i, j) > kMax) {
                    kMax = dm.getDistance(i, j);
                }
            }
        }
        kAvg = total / nk;
        kSD = Math.sqrt((totalsq - nk * kAvg * kAvg) / (nk - 1));
        System.out.println(kAvg);
    }

    public void pullBackExtrem() {
        //take values beyond 3 sd from mean back
        //By Zhiwu Zhang
        for (int i = 0; i < numSeqs - 1; i++) {
            for (int j = i + 1; j < numSeqs; j++) {
                if (dm.getDistance(i, j) < kAvg - cutOff * kSD) {
                    dm.setDistance(i, j, kAvg - cutOff * kSD);
                    kMin = dm.getDistance(i, j);
                }
            }
        }
        System.out.println("values beyond 3 sd from mean were pulled back");
    }

    public void cutOff() {
        //Set vale to 0 if below than avg
        //By Zhiwu Zhang
        for (int i = 0; i < numSeqs; i++) {
            for (int j = i + 0; j < numSeqs; j++) {
                if (dm.getDistance(i, j) < kAvg) {
                    dm.setDistance(i, j, kAvg);
                }
            }
        }
        kMin = kAvg;
    }

    public void rescale() {
        //rescale from theMin~2 to 0~2
        //By Zhiwu Zhang
        double s;
        for (int i = 0; i < numSeqs; i++) {
            for (int j = i; j < numSeqs; j++) {
                s = (dm.getDistance(i, j) - kMin) * cutOff / (cutOff - kMin);
                dm.setDistance(i, j, s);
                dm.setDistance(j, i, s);
            }
        }
        System.out.println("K rescaled");

    }

    /**
     * Calculates a kinship matrix from genotypes using the method described in
     * Endelman and Jannink (2012) G3 2:1407-1413. It is best to impute missing
     * data before calculating. However, if data is missing it is replaced by
     * the allele average at that site.
     *
     */
    public void calculateKinshipFromMarkers() {
        //mar is the input genotype table
        byte missingAllele = GenotypeTable.UNKNOWN_ALLELE;

		// from Endelman and Jannink. 2012. G3 2:1407ff
        // A = WW'/[2*sumk(pk*qk)]
        // where W = centered genotype matrix (centered on marker mean value, marker coded as 2,1,0)
        // where marker is multi-allelic, leave one allele out to keep markers independent
        int ntaxa = mar.numberOfTaxa();
        int nsites = mar.numberOfSites();
        distance = new double[ntaxa][ntaxa];
        double sumpi = 0;

        //calculate WW' by summing ww' for each allele, where w is a column vector of centered allele counts {2,1,0}
        for (int s = 0; s < nsites; s++) {
            int[][] alleleFreq = mar.allelesSortedByFrequency(s);
            int nalleles = alleleFreq[0].length;
            int totalAlleleCount = mar.totalGametesNonMissingForSite(s);

            for (int a = 0; a < nalleles - 1; a++) {
                double pi = ((double) alleleFreq[1][a]) / ((double) totalAlleleCount);
                double pix2 = 2 * pi;
                sumpi += pi * (1 - pi);
                double[] scores = new double[ntaxa];

                for (int t = 0; t < ntaxa; t++) {
                    byte[] geno = GenotypeTableUtils.getDiploidValues(mar.genotype(t, s));
                    double thisScore = 0;
                    if (geno[0] != missingAllele) {
                        if (geno[0] == alleleFreq[0][a]) {
                            thisScore++;
                        }
                        if (geno[1] == alleleFreq[0][a]) {
                            thisScore++;
                        }
                        thisScore -= pix2;
                    }
                    scores[t] = thisScore;
                }

                for (int r = 0; r < ntaxa; r++) {
                    double rowval = scores[r];
                    distance[r][r] += rowval * rowval;
                    for (int c = r + 1; c < ntaxa; c++) {
                        distance[r][c] += rowval * scores[c];
                    }
                }
            }
        }

        double sumpk = 2 * sumpi;

        for (int r = 0; r < ntaxa; r++) {
            distance[r][r] /= sumpk;
            for (int c = r + 1; c < ntaxa; c++) {
                distance[r][c] = distance[c][r] = distance[r][c] / sumpk;
            }
        }

        dm = new DistanceMatrix(distance, mar.taxa());
    }

    /**
     * Calculates a kinship matrix from genotypes using the method described in
     * Endelman and Jannink (2012) G3 2:1407-1413. It is best to impute missing
     * data before calculating. However, if data is missing it is replaced by
     * the allele average at that site.
     *
     */
    public void calculateKinshipFromMarkersV2() {
        //mar is the input genotype table
        byte missingAllele = GenotypeTable.UNKNOWN_ALLELE;

		// from Endelman and Jannink. 2012. G3 2:1407ff
        // A = WW'/[2*sumk(pk*qk)]
        // where W = centered genotype matrix (centered on marker mean value, marker coded as 2,1,0)
        // where marker is multi-allelic, leave one allele out to keep markers independent
        int ntaxa = mar.numberOfTaxa();
        int nsites = mar.numberOfSites();
        distance = new double[ntaxa][ntaxa];
        double sumpq = IntStream.range(0, nsites).parallel().mapToDouble(s -> {
            int[][] alleleFreq = mar.allelesSortedByFrequency(s);
            int nalleles = alleleFreq[0].length;
            double totalAlleleCount = mar.totalGametesNonMissingForSite(s);
            double pq = 0;
            for (int a = 0; a < nalleles - 1; a++) {
                double p = alleleFreq[1][a] / totalAlleleCount;
                pq += p * (1 - p);
            }
            return pq;
        }).sum();

        //calculate WW' by summing ww' for each allele, where w is a column vector of centered allele counts {2,1,0}
        BiConsumer<double[][], Integer> siteDistance = (dist, site) -> {
            int s = site.intValue();
            int[][] alleleFreq = mar.allelesSortedByFrequency(s);
            int nalleles = alleleFreq[0].length;
            int totalAlleleCount = mar.totalGametesNonMissingForSite(s);

            for (int a = 0; a < nalleles - 1; a++) {
                double pi = ((double) alleleFreq[1][a]) / ((double) totalAlleleCount);
                double pix2 = 2 * pi;
                double[] scores = new double[ntaxa];

                for (int t = 0; t < ntaxa; t++) {
                    byte[] geno = GenotypeTableUtils.getDiploidValues(mar.genotype(t, s));
                    double thisScore = 0;
                    if (geno[0] != missingAllele) {
                        if (geno[0] == alleleFreq[0][a]) {
                            thisScore++;
                        }
                        if (geno[1] == alleleFreq[0][a]) {
                            thisScore++;
                        }
                        thisScore -= pix2;
                    }
                    scores[t] = thisScore;
                }

                for (int r = 0; r < ntaxa; r++) {
                    double rowval = scores[r];
                    dist[r][r] += rowval * rowval;
                    for (int c = r + 1; c < ntaxa; c++) {
                        dist[r][c] += rowval * scores[c];
                    }
                }
            }
        };

        BiConsumer<double[][], double[][]> mergeDistance = (d1, d2) -> {
            for (int r = 0; r < ntaxa; r++) {
                double[] row1 = d1[r];
                double[] row2 = d2[r];
                for (int c = 0; c < ntaxa; c++) {
                    row1[c] += row2[c];
                }
            }
        };

        distance = IntStream.range(0, nsites).boxed().parallel().collect(() -> new double[ntaxa][ntaxa], siteDistance, mergeDistance);

        double sumpk = 2 * sumpq;

        for (int r = 0; r < ntaxa; r++) {
            distance[r][r] /= sumpk;
            for (int c = r + 1; c < ntaxa; c++) {
                distance[r][c] = distance[c][r] = distance[r][c] / sumpk;
            }
        }

        dm = new DistanceMatrix(distance, mar.taxa());
        System.out.print("");
    }

    public void calculateRelationshipKinshipFromReferenceProbability() {
        ReferenceProbability referenceP = mar.referenceProbability();

        //calculate the column averages and sumpq, center W
        int nrow = referenceP.numTaxa();
        int ncol = referenceP.numSites();
        double[][] W = new double[nrow][ncol];

        for (int r = 0; r < nrow; r++) {
            for (int c = 0; c < ncol; c++) {
                W[r][c] = referenceP.value(r, c) * matrixMultiplier;
            }
        }

        double sumpq = 0;
        for (int c = 0; c < ncol; c++) {
            double colTotal = 0;
            int colCount = 0;
            for (int r = 0; r < nrow; r++) {
                if (!Double.isNaN(W[r][c])) {
                    colTotal += W[r][c];
                    colCount++;
                }
            }

            double pi = colTotal / colCount / 2.0;
            double pix2 = pi * 2;
            sumpq += pi * (1 - pi);
            for (int r = 0; r < nrow; r++) {
                if (Double.isNaN(W[r][c])) {
                    W[r][c] = 0;
                } else {
                    W[r][c] -= pix2;
                }
            }
        }

        DoubleMatrix WWt = DoubleMatrixFactory.DEFAULT.make(W).tcrossproduct();

        double[][] scaledIBS = new double[nrow][nrow];
        for (int r = 0; r < nrow; r++) {
            for (int c = 0; c < nrow; c++) {
                scaledIBS[r][c] = WWt.get(r, c) / sumpq / 2;
            }
        }
        dm = new DistanceMatrix(scaledIBS, mar.taxa());
    }

    public void calculateRelationshipKinshipFromAlleleProbabilities() {
        //TODO implement
    }

    public DistanceMatrix getDm() {
        return dm;
    }

}
