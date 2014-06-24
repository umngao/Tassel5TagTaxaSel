package net.maizegenetics.analysis.imputation;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.FilterGenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginEvent;

import org.apache.log4j.Appender;
import org.apache.log4j.FileAppender;
import org.apache.log4j.Logger;
import org.apache.log4j.PatternLayout;
import org.apache.log4j.xml.DOMConfigurator;

import javax.swing.*;
import java.awt.*;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

public class CallParentAllelesPlugin extends AbstractPlugin {
	private static final Logger myLogger = Logger.getLogger(CallParentAllelesPlugin.class);
	private String pedfileName = null;
	private int minAlleleCount = 2; //minimum allele count for the minor allele to consider that a site might be polymorphic
	private int windowSize = 50;  //the number of sites to be used a window for determining the original set of snps in LD
	private double cutHeightSnps = 0.2;  //the tree cut height used to find the largest cluster of correlated SNPs
	private double minRforSnps = 0.2;  //the minimum R used to judge whether a snp is in ld with a test group
	private double maxMissing = 0.9;
	private double minMinorAlleleFrequency = -1.0;
	private boolean useBCFilter = true;
	private boolean useMultipleBCFilter = false;
	private boolean useClusterAlgorithm = false;
	private boolean checkSubPops = false;
	private boolean useHets = true;
	private ArrayList<PopulationData> familyList = null;
	
	public CallParentAllelesPlugin(Frame parentFrame) {
        super(parentFrame, false);
	}
	
	@Override
	public DataSet performFunction(DataSet input) {
		if (pedfileName == null && familyList == null) {
			myLogger.error(getUsage());
			return null;
		}
		
		List<Datum> inputAlignments = input.getDataOfType(GenotypeTable.class);
		LinkedList<Datum> datumList = new LinkedList<Datum>();

		for (Datum d : inputAlignments) {
			GenotypeTable align = (GenotypeTable) d.getData();
			if (familyList == null) familyList = PopulationData.readPedigreeFile(pedfileName);
			for (PopulationData family : familyList) {
				myLogger.info("Calling parent alleles for family " + family.name + ", chromosome " + align.chromosomeName(0) + ".");
				
				String[] ids = new String[family.members.size()];
				family.members.toArray(ids);
				
				myLogger.info("creating family alignment for family " + family.name);
                TaxaList tL=new TaxaListBuilder().addAll(ids).build();
				family.original =  GenotypeTableBuilder.getGenotypeCopyInstance((FilterGenotypeTable)FilterGenotypeTable.getInstance(align, tL, false));
				
				if (!useHets) {
					byte NN = NucleotideAlignmentConstants.getNucleotideDiploidByte('N');
//					MutableNucleotideAlignment mna = MutableNucleotideAlignment.getInstance(family.original);
					GenotypeTableBuilder builder = GenotypeTableBuilder.getSiteIncremental(family.original.taxa());
					int nsites = family.original.numberOfSites();
					int ntaxa = family.original.numberOfTaxa();
					for (int s = 0; s < nsites; s++) {
						byte[] siteGeno = family.original.genotypeAllTaxa(s);
						for (int t = 0; t < ntaxa; t++) {
							if (GenotypeTableUtils.isHeterozygous(siteGeno[t])) siteGeno[t] = NN;
						}
						builder.addSite(family.original.positions().get(s), siteGeno);
					}
					family.original = builder.build();
				}
				
				myLogger.info("family alignment created");
				if (useClusterAlgorithm)  {
					BiparentalHaplotypeFinder hapFinder = new BiparentalHaplotypeFinder(family);
					hapFinder.assignHaplotyes();
					hapFinder.convertGenotypesToParentCalls();
				}
//				if (useClusterAlgorithm)  NucleotideImputationUtils.callParentAllelesUsingClusters(family, maxMissing, minMinorAlleleFrequency, windowSize, checkSubPops);
				else if (useBCFilter && (family.contribution1 == 0.75 || family.contribution1 == 0.25)) NucleotideImputationUtils.callParentAllelesByWindowForBackcrosses(family, maxMissing, minMinorAlleleFrequency, windowSize, minRforSnps);
				else if (useMultipleBCFilter) NucleotideImputationUtils.callParentAllelesByWindowForMultipleBC(family, maxMissing, 1, windowSize);
				else NucleotideImputationUtils.callParentAllelesByWindow(family, maxMissing, minMinorAlleleFrequency, windowSize, minRforSnps);
				String comment = "Parent Calls for family " + family.name + " from " + d.getName() + ".";
				datumList.add(new Datum(family.name, family, comment));
			}
		}
		
		DataSet resultDS =  new DataSet(datumList, this);
		fireDataSetReturned(new PluginEvent(resultDS, CallParentAllelesPlugin.class));
		return resultDS;
	}

	@Override
	public void setParameters(String[] args) {
		if (args == null || args.length == 0) {
			myLogger.error(getUsage());
			return;
		}
		
		int narg = args.length;
		for (int i = 0; i < narg; i++) {
			if (args[i].equals("-p") || args[i].equalsIgnoreCase("-pedigrees")) {
				pedfileName = args[++i];
			}
			else if (args[i].equals("-a") || args[i].equalsIgnoreCase("-minAlleleCount")) {
				minAlleleCount = Integer.parseInt(args[++i]);
			}
			else if (args[i].equals("-w") || args[i].equalsIgnoreCase("-windowSize")) {
				windowSize = Integer.parseInt(args[++i]);
			}
			else if (args[i].equals("-h") || args[i].equalsIgnoreCase("-cutHeight")) {
				cutHeightSnps = Double.parseDouble(args[++i]);
			}
			else if (args[i].equals("-r") || args[i].equalsIgnoreCase("-minR")) {
				minRforSnps = Double.parseDouble(args[++i]);
			}
			else if (args[i].equals("-m") || args[i].equalsIgnoreCase("-maxMissing")) {
				maxMissing = Double.parseDouble(args[++i]);
			}
			else if (args[i].equals("-f") || args[i].equalsIgnoreCase("-minMaf")) {
				minMinorAlleleFrequency = Double.parseDouble(args[++i]);
			}
			else if (args[i].equals("-b") || args[i].equalsIgnoreCase("-bc1")) {
				String param = args[++i];
				if (param.toUpperCase().startsWith("F")) useBCFilter = false;
			}
			else if (args[i].equals("-n") || args[i].equalsIgnoreCase("-bcn")) {
				String param = args[++i];
				if (param.toUpperCase().startsWith("T")) useMultipleBCFilter = true;
			}
			else if (args[i].equals("-l") || args[i].equalsIgnoreCase("-logconfig")) {
				addFileLogger(args[++i]);
			}
			else if (args[i].equals("-logfile")) {
				setFileLogger(args[++i]);
			}
			else if (args[i].startsWith("-clust")) {
				useClusterAlgorithm = true;
			}
			else if (args[i].equals("-subpops")) {
				checkSubPops = true;
			}
			else if (args[i].equals("-nohets")) {
				useHets = false;
			}
			
			else if (args[i].equals("?")) myLogger.info(getUsage());
		}
	}
	
	private void addFileLogger(String filename) {
		try {
			Appender fileAppender = Logger.getRootLogger().getAppender("fileAppender");
			if (fileAppender == null) {
				FileAppender filelog = new FileAppender(new PatternLayout("%d %-5p  [%c{1}] %m %n"), filename, true);
				filelog.setName("fileAppender");
				Logger.getRootLogger().addAppender(filelog);
			}
		} catch(Exception e) {
			myLogger.info("log file could not be instantiated");
			e.printStackTrace();
		}
	}
	
	private void setFileLogger(String filename) {
		try {
			Appender fileAppender = Logger.getRootLogger().getAppender("fileAppender");
			if (fileAppender != null) Logger.getRootLogger().removeAppender(fileAppender);
			FileAppender filelog = new FileAppender(new PatternLayout("%d %-5p  [%c{1}] %m %n"), filename, true);
			filelog.setName("fileAppender");
			Logger.getRootLogger().addAppender(filelog);
		} catch(Exception e) {
			myLogger.info("log file could not be instantiated");
			e.printStackTrace();
		}
	}
	
	public void setPedfileName(String pedfileName) {
		this.pedfileName = pedfileName;
	}

	public void setMinAlleleCount(int minAlleleCount) {
		this.minAlleleCount = minAlleleCount;
	}

	public void setWindowSize(int windowSize) {
		this.windowSize = windowSize;
	}

	public void setCutHeightSnps(double cutHeightSnps) {
		this.cutHeightSnps = cutHeightSnps;
	}

	public void setMinRforSnps(double minRforSnps) {
		this.minRforSnps = minRforSnps;
	}

	@Override
	public ImageIcon getIcon() {
		return null;
	}

	@Override
	public String getButtonName() {
		return "Call Parents";
	}

	@Override
	public String getToolTipText() {
		return null;
	}

	public String getUsage() {
		StringBuilder usage = new StringBuilder("The CallParentAllelesPlugin requires the following parameter:\n");
		usage.append("-p or -pedigrees : a file containing pedigrees of the individuals to be imputed\n");
		usage.append("The following parameters are optional:\n");
		usage.append("-w or -windowSize : the number of SNPs to examine for LD clusters (default = 50)\n");
		usage.append("-r or -minR : minimum R used to filter SNPs on LD (default = 0.2, use 0 for no ld filter)\n");
		usage.append("-m or -maxMissing : maximum proportion of missing data allowed for a SNP (default = 0.9)\n");
		usage.append("-f or -minMaf : minimum minor allele frequency used to filter SNPs. If negative, filters on expected segregation ratio from parental contribution (default = -1)\n");
		usage.append("-b or -bc1 : use BC1 specific filter (default = true)\n");
		usage.append("-n or -bcn : use multipe backcross specific filter (default = false)\n");
		usage.append("-logfile : the name of a file to which all logged messages will be printed.\n");
		usage.append("-cluster : use the cluster algorithm. minMaf defaults to 0.05.\n");
		usage.append("-subpops : filter sites for heterozygosity in subpopulations.\n");
		usage.append("-nohets : delete het calls from original data before imputing.\n");
		usage.append("? : print the parameter list.\n");

		return usage.toString();
	}

	public void setFamilyList(ArrayList<PopulationData> familyList) {
		this.familyList = familyList;
	}
}
