package net.maizegenetics.analysis.imputation;

import java.awt.Frame;
import java.util.Arrays;
import java.util.List;

import javax.swing.ImageIcon;

import org.apache.log4j.Logger;

import net.maizegenetics.dna.snp.CombineGenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.ExportUtils;
import net.maizegenetics.dna.snp.FilterGenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginEvent;

public class WritePopulationAlignmentPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(WritePopulationAlignmentPlugin.class);
    boolean mergeAlignments = false;
    boolean writeParentCalls = true;
    boolean writeNucleotides = true;
    boolean outputDiploid = false;
    double minSnpCoverage = 0.1;
    double maxMafForMono = 0.01;
    String baseFileName;

    public WritePopulationAlignmentPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    @Override
    public DataSet performFunction(DataSet input) {
        List<Datum> theData = input.getDataOfType(PopulationData.class);

        if (theData.size() > 0) {
            if (writeParentCalls) {
                writeOutput(theData, false);
            }
            if (writeNucleotides) {
                writeOutput(theData, true);
            }
            fireDataSetReturned(new PluginEvent(input, input.getCreator().getClass()));
            return input;
        } else {
            return null;
        }

    }

    private void writeOutput(List<Datum> theData, boolean asNucleotides) {
        String filename;
        if (mergeAlignments) {
            if (asNucleotides) {
                filename = baseFileName + ".nuc.hmp.txt";
            } else {
                filename = baseFileName + ".parents.hmp.txt";
            }
            GenotypeTable[] allOfTheAlignments = new GenotypeTable[theData.size()];
            int count = 0;
            for (Datum datum : theData) {
                PopulationData family = (PopulationData) datum.getData();
                allOfTheAlignments[count++] = createOutputAlignment(family, asNucleotides);
            }
            GenotypeTable alignment = CombineGenotypeTable.getInstance(allOfTheAlignments, true);
            ExportUtils.writeToHapmap(alignment, outputDiploid, filename, '\t', null);
        } else {
            for (Datum datum : theData) {
                PopulationData family = (PopulationData) datum.getData();
                String familyName = family.name.replace('/', '.');
                String chrName = family.original.chromosomeName(0);
                if (asNucleotides) {
                    filename = baseFileName + ".family." + familyName + ".chr" + chrName + ".nuc.hmp.txt";
                } else {
                    filename = baseFileName + ".family." + familyName + ".chr" + chrName + ".parents.hmp.txt";
                }
                ExportUtils.writeToHapmap(createOutputAlignment(family, asNucleotides), outputDiploid, filename, '\t', null);
            }
        }
    }

    private GenotypeTable createOutputAlignment(PopulationData popdata, boolean asNucleotides) {
    	GenotypeTable out = null;

        if (!asNucleotides) {
            out = popdata.imputed;
        } else {
            //change the parent calls to original nucleotides
        	GenotypeTable outPoly = NucleotideImputationUtils.convertParentCallsToNucleotides(popdata);

            if (!Double.isNaN(minSnpCoverage) && !Double.isNaN(maxMafForMono)) {
                int nsnps = popdata.original.numberOfSites();
                double ngametes = 2 * popdata.original.numberOfTaxa();
                
                int[] monomorphicSnps = new int[nsnps];
                int snpCount = 0;
                for (int s = 0; s < nsnps; s++) {
                    double coverage = popdata.original.totalGametesNonMissingForSite(s) / ngametes;
                    if (!popdata.snpIndex.fastGet(s) && popdata.original.minorAlleleFrequency(s) <= maxMafForMono && coverage >= minSnpCoverage) {
                        monomorphicSnps[snpCount++] = s;
                    }
                }
                monomorphicSnps = Arrays.copyOf(monomorphicSnps, snpCount);
                GenotypeTable fa = FilterGenotypeTable.getInstance(popdata.original, monomorphicSnps);
                if (fa.numberOfSites() == 0) {	//If there are no monomorphic sites (e.g, have been pre-filtered), just return polymorphic ones
                    out = outPoly;
                } else { //Return both monomorphic and polymorphic sites
                	GenotypeTableBuilder builder = GenotypeTableBuilder.getSiteIncremental(fa.taxa());
                	
                    // fill in all values with the major allele
                    nsnps = fa.numberOfSites();
                    int ntaxa = fa.numberOfTaxa();
                    for (int s = 0; s < nsnps; s++) {
                        byte majorAllele = fa.majorAllele(s);
                        byte major = (byte) ((majorAllele << 4) | majorAllele);
                        byte[] snpgeno = new byte[ntaxa];
                        Arrays.fill(snpgeno, major);
                        builder.addSite(fa.positions().get(s), snpgeno);
                    }
                    out = builder.build();
                }
            }
        }
        return out;
    }

    @Override
    public void setParameters(String[] args) {
        if (args == null || args.length == 0) {
            myLogger.info(getUsage());
            return;
        }

        int narg = args.length;
        for (int i = 0; i < narg; i++) {
            if (args[i].equals("-f") || args[i].equalsIgnoreCase("-file")) {
                baseFileName = args[++i];
            } else if (args[i].equals("-m") || args[i].equalsIgnoreCase("-merge")) {
                String val = args[++i];
                if (val.toUpperCase().startsWith("T")) {
                    mergeAlignments = true;
                } else {
                    mergeAlignments = false;
                }
            } else if (args[i].equals("-p") || args[i].equalsIgnoreCase("-parentCalls")) {
                String val = args[++i];
                if (val.toUpperCase().startsWith("T")) {
                    writeParentCalls = true;
                    writeNucleotides = false;
                } else {
                    writeParentCalls = false;
                    writeNucleotides = true;
                }
            } else if (args[i].equals("-o") || args[i].equalsIgnoreCase("-outputType")) {
                String val = args[++i];
                if (val.toUpperCase().startsWith("P")) {
                    writeParentCalls = true;
                    writeNucleotides = false;
                } else if (val.toUpperCase().startsWith("N")) {
                    writeParentCalls = false;
                    writeNucleotides = true;
                } else if (val.toUpperCase().startsWith("B")) {
                    writeParentCalls = true;
                    writeNucleotides = true;
                } else {
                    writeParentCalls = true;
                    writeNucleotides = false;
                }
            } else if (args[i].equals("-d") || args[i].equalsIgnoreCase("-diploid")) {
                String val = args[++i];
                if (val.toUpperCase().startsWith("T")) {
                    outputDiploid = true;
                } else {
                    outputDiploid = false;
                }
            } else if (args[i].equals("-c") || args[i].equalsIgnoreCase("-minCoverage")) {
                minSnpCoverage = Double.parseDouble(args[++i]);
            } else if (args[i].equals("-x") || args[i].equalsIgnoreCase("-maxMono")) {
                maxMafForMono = Double.parseDouble(args[++i]);
            } else if (args[i].equals("?")) {
                myLogger.info(getUsage());
            }
        }
    }

    public void setMergeAlignments(boolean mergeAlignments) {
        this.mergeAlignments = mergeAlignments;
    }

    public void setWriteParentCalls(boolean writeParentCalls) {
        this.writeParentCalls = writeParentCalls;
    }

    public void setOutputDiploid(boolean outputDiploid) {
        this.outputDiploid = outputDiploid;
    }

    public void setBaseFileName(String baseFileName) {
        this.baseFileName = baseFileName;
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "Write Populations";
    }

    @Override
    public String getToolTipText() {
        return null;
    }

    public String getUsage() {
        StringBuilder usage = new StringBuilder("The WritePopulationAlignmentPlugin requires the following parameter:\n");
        usage.append("-f or -file : The base file name for the ouput. .hmp.txt will be appended.\n");
        usage.append("The following parameters are optional:\n");
        usage.append("-m or -merge : if true families are merged into a single file, if false each family is output to a separate file (default = false)\n");
        usage.append("-o or -outputType : parents = output parent calls, nucleotides = output nucleotides, both = output both\n");
        usage.append("-d or -diploid : if true output is AA/CC/AC, if false output is A/C/M\n");
        usage.append("-c or -minCoverage : the minimum coverage for a monomorphic snp to be included in the nucleotide output (default = 0.1)\n");
        usage.append("-x or -maxMono : the maximum minor allele frequency used to call monomorphic snps (default = 0.01)\n");
        usage.append("if -c or -x equals NaN and merge is true, then missing values at monomorphic sites (within a family) will be left missing\n");
        usage.append("? : print the parameter list.\n");

        return usage.toString();
    }

	public void setWriteNucleotides(boolean writeNucleotides) {
		this.writeNucleotides = writeNucleotides;
	}

	public void setMinSnpCoverage(double minSnpCoverage) {
		this.minSnpCoverage = minSnpCoverage;
	}

	public void setMaxMafForMono(double maxMafForMono) {
		this.maxMafForMono = maxMafForMono;
	}
}
