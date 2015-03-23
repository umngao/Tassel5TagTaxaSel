package net.maizegenetics.analysis.imputation;

import java.awt.Frame;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import javax.swing.ImageIcon;
import javax.swing.JOptionPane;

import org.apache.log4j.Logger;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;

public class RandomGenotypeImputationPlugin extends AbstractPlugin {
    private static final Logger myLogger = Logger.getLogger(RandomGenotypeImputationPlugin.class);

	public RandomGenotypeImputationPlugin(Frame parentFrame, boolean isInteractive) {
		super(parentFrame, isInteractive);
	}

	@Override
	public DataSet processData(DataSet input) {
		List<Datum> genotypeList = input.getDataOfType(GenotypeTable.class);
		if (genotypeList.size() == 0) {
			String errmsg = "Error in random imputation: no appropriate data set was selected.";
			if (isInteractive()) JOptionPane.showMessageDialog(getParentFrame(), errmsg, "Error", JOptionPane.ERROR_MESSAGE);
			else myLogger.error(errmsg);
			return null;
		}
		
		Random ran = new Random();
		byte NN = GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
		List<Datum> outData = new ArrayList<>();
		String comment = "Missing genotypes imputed randomly\nImputed genotypes selected from genotype distribution for each site.";
		for (Datum nextDatum : genotypeList) {
			GenotypeTable myGenotype = (GenotypeTable) nextDatum.getData();
			GenotypeTableBuilder myBuilder = GenotypeTableBuilder.getSiteIncremental(myGenotype.taxa());
			int nsites = myGenotype.numberOfSites();
			int ntaxa = myGenotype.numberOfTaxa();
			for (int s = 0; s < nsites; s++) {
				byte[] sitegeno = myGenotype.genotypeAllTaxa(s);
				Object[] genotypeCounts = byteCounts(sitegeno);
				byte[] siteGenotypes = (byte[]) genotypeCounts[0];
				int[] genoCounts = (int[]) genotypeCounts[1];
				int maxCount = genoCounts[genoCounts.length - 1];
				for (int t = 0; t < ntaxa; t++) {
					if (sitegeno[t] == NN) {
						int ranval = ran.nextInt(maxCount);
						int genoIndex = 0;
						while (ranval > genoCounts[genoIndex]) genoIndex++;
						sitegeno[t] = siteGenotypes[genoIndex];
					}
				}
				myBuilder.addSite(myGenotype.positions().get(s), sitegeno);
			}
			
			Datum outDatum = new Datum("Imputed_" + nextDatum.getName(), myBuilder.build(), comment);
			outData.add(outDatum);
			
		}
		
		return new DataSet(outData, this);
	}

	public static Object[] byteCounts(byte[] genotypes) {
		//Object[0] is a byte[] array of genotypes
		//Object[1] is an int[] array of the cumulative counts of the genotypes
		byte NN = GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
		int n = genotypes.length;
		Map<Byte, Long> byteCounts = IntStream.range(0, n).filter(i -> genotypes[i] != NN).boxed()
				.collect(Collectors.groupingBy(i -> new Byte(genotypes[i]), Collectors.counting()));
		int mapSize = byteCounts.size();
		byte[] geno = new byte[mapSize];
		int[] genocount = new int[mapSize];
		int ndx = 0;
		int sum = 0;
		for (Map.Entry<Byte, Long> me : byteCounts.entrySet()) {
			geno[ndx] = me.getKey().byteValue();
			sum += me.getValue().intValue();
			genocount[ndx] = sum;
			ndx++;
		}
		return new Object[] {geno, genocount};
	}
	
	@Override
	public ImageIcon getIcon() {
		return null;
	}

	@Override
	public String getButtonName() {
		return "Random Imputation";
	}

	@Override
	public String getToolTipText() {
		return "Replace missing genotypes with a value drawn from the site genotype distribution.";
	}

	
}
