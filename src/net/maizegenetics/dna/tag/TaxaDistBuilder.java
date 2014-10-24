package net.maizegenetics.dna.tag;

import cern.colt.list.ShortArrayList;

import com.google.common.collect.ImmutableMultiset;
import com.google.common.collect.Multiset;
import com.google.common.primitives.Ints;
import com.google.common.primitives.Shorts;
import com.google.common.primitives.UnsignedBytes;

import org.xerial.snappy.Snappy;

import gnu.trove.iterator.TIntShortIterator;
import gnu.trove.iterator.TShortShortIterator;
import gnu.trove.map.hash.TIntShortHashMap;
import gnu.trove.map.hash.TShortShortHashMap;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.util.Arrays;
import java.util.stream.Collectors;


/**
 * Builder for TaxaDistribution.  Deals with the
 */
public class TaxaDistBuilder {
	private TaxaDistBuilder() {}

	/**
	 * Create a TaxaDistribution with only a single taxa with a tag.  Very memory efficient, but needs conversion add
	 * an additional taxon with a tag
	 * @param maxTaxa
	 * @param taxaWithTag
	 * @return
	 */
	public static TaxaDistribution create(int maxTaxa, int taxaWithTag) {
		return new TaxaDistSingleTaxon(maxTaxa,taxaWithTag);
	}

	/**
	 * Create an expandable TaxaDistribution with no initial values
	 * @param maxTaxa
	 * @return
	 */
	public static TaxaDistribution create(int maxTaxa) {
		return new TaxaDistExpandable(maxTaxa);

	}

	/**
	 * Create an fixed TaxaDistribution with set values
	 * @param maxTaxa
	 * @return
	 */
	public static TaxaDistribution create(int maxTaxa, int[] taxaWithTags, int[] depthOfTags) {
		return new TaxaDistFixed(maxTaxa,taxaWithTags,depthOfTags);

	}

	/**
	 * Create an fixed TaxaDistribution with set values
	 * @param encodedTaxaDistribution byte array of encoded TaxaDistribution
	 * @return
	 */
	public static TaxaDistribution create(byte[] encodedTaxaDistribution) {
		int[][] decodedTD=getDepthMatrixForEncodedDepths(encodedTaxaDistribution);
		return new TaxaDistFixed(decodedTD[2][0],decodedTD[0], decodedTD[1]);

	}

	/**
	 * Copies a TaxaDistribution to an expandable TaxaDistribution.  Can be used to convert a single TaxaDistribution
	 * to an expandable version
	 * @param srcTaxaDist
	 * @return
	 */
	public static TaxaDistribution create(TaxaDistribution srcTaxaDist) {
		TaxaDistribution dstTD=new TaxaDistExpandable(srcTaxaDist.maxTaxa());
		int[] depths=srcTaxaDist.depths();
		for (int taxaIndex = 0; taxaIndex < depths.length; taxaIndex++) {
			for (int j = 0; j < depths[taxaIndex]; j++) {
				dstTD.increment(taxaIndex);
			}
		}
		return dstTD;
	}

	/**
	 * Combines a TaxaDistribution to an expandable TaxaDistribution.  Can be used to convert a single TaxaDistribution
	 * to an expandable version
	 * @param taxaDist1
	 * @return
	 */
	public static TaxaDistribution combine(TaxaDistribution taxaDist1, TaxaDistribution taxaDist2) {
		if(taxaDist1.maxTaxa()!=taxaDist2.maxTaxa()) throw new IllegalStateException("TaxaDistributions not of same size");
		int[] depths1=taxaDist1.depths();
		int[] depths2=taxaDist2.depths();
		int taxaWithDepth=0;
		for (int i = 0; i < depths1.length; i++) {
			depths1[i]+=depths2[i];
			if(depths1[i]>0) taxaWithDepth++;
		}
		int[] taxaWithTags=new int[taxaWithDepth];
		int[] depthOfTags=new int[taxaWithDepth];
		taxaWithDepth=0;
		for (int i = 0; i < depths1.length; i++) {
			if(depths1[i]>0) {
				taxaWithTags[taxaWithDepth]=i;
				depthOfTags[taxaWithDepth]=depths1[i];
				taxaWithDepth++;
			}
		}
		return create(taxaDist1.maxTaxa(),taxaWithTags,depthOfTags);
	}

	/**
	 * Create expandable TaxaDistribution from encoded taxa distribution.  The int[] encoding use first three bytes
	 * for taxa index, and last byte for depth as unsigned byte.  Depths greater than 256 just increase.
	 * @param maxTaxa
	 * @param encodeTaxaDepths
	 * @return
	 */
	public static TaxaDistribution create(int maxTaxa, int[] encodeTaxaDepths) {
		TaxaDistribution dstTD=new TaxaDistExpandable(maxTaxa);
		for (int taxaDepth : encodeTaxaDepths) {
			int taxa=taxaDepth>>>8;
		int depth=UnsignedBytes.toInt((byte)taxaDepth);
		for (int i = 0; i < depth; i++) {
			dstTD.increment(taxa);
		}
		}
		return dstTD;
	}

	private static int[][] getDepthMatrixForEncodedDepths(byte[] input) {
		try{
			final int maxValueInInt=UnsignedBytes.toInt(UnsignedBytes.MAX_VALUE);
			ByteBuffer bb=ByteBuffer.wrap(Snappy.uncompress(input));
			int maxTaxa=bb.getInt();
			int taxaWithDepth=bb.getInt();
			int[][] result=new int[3][];
			result[0]=new int[taxaWithDepth];
			result[1]=new int[taxaWithDepth];
			result[2]=new int[]{maxTaxa};
			for (int i = 0; i < taxaWithDepth; i++) {
				int space=(i>0)?result[0][i-1]:0;
				int inc=UnsignedBytes.toInt(bb.get());
				while(inc==maxValueInInt) {
					space+=inc;
					inc=UnsignedBytes.toInt(bb.get());
				}
				space+=inc;
				result[0][i]=space;
			}
			for (int i = 0; i < taxaWithDepth; i++) {
				int depth=0;
				int inc=UnsignedBytes.toInt(bb.get());
				while(inc==maxValueInInt) {
					depth+=inc;
					inc=UnsignedBytes.toInt(bb.get());
				}
				depth+=inc;
				result[1][i]=depth;
			}
			return result;
		} catch (IOException e) {
			e.printStackTrace();
			return null;
		}



	}


}

abstract class AbstractTaxaDistribution implements TaxaDistribution {

	@Override
	public byte[] encodeTaxaDepth() {
		int[][] tds=taxaWithDepths();
		ByteBuffer bb=ByteBuffer.allocate(8+(maxTaxa()/64)+(2*tds[0].length)+(totalDepth()/64));
		bb.putInt(maxTaxa());  //maximum number of taxa with depth
		bb.putInt(tds[0].length);  //number of taxa with depth
		bb.put(UnsignedBytes.checkedCast(tds[0][0]));
		for (int i = 1; i < tds[0].length; i++) {
			int space=tds[0][i]-tds[0][i-1];
			while(space>=0) {
				bb.put(UnsignedBytes.saturatedCast(space));
				space-=255;
			}
		}
		for (int i = 0; i < tds[1].length; i++) {
			int depth=tds[1][i];
			while(depth>=0) {
				bb.put(UnsignedBytes.saturatedCast(depth));
				depth-=255;
			}
		}
		try{
			return Snappy.compress(Arrays.copyOf(bb.array(), bb.position()));
		} catch (IOException e) {
			e.printStackTrace();
		}
		return null;
	}

	@Override
	public Multiset<Integer> taxaDepthMap() {
		int[][] tds=taxaWithDepths();
		ImmutableMultiset.Builder<Integer> taxaCnts = new ImmutableMultiset.Builder<>();
		for (int i = 0; i < tds[0].length; i++) {
			taxaCnts.setCount(tds[0][i],tds[1][i]);
		}
		return taxaCnts.build();
	}

	@Override
	public boolean equals(Object o) {
		if (this == o) return true;
		if (o == null) return false;
		if (!(o instanceof TaxaDistribution)) return false;

		TaxaDistribution that = (TaxaDistribution) o;
		if (maxTaxa() != that.maxTaxa()) return false;
		int[][] thisTDS=taxaWithDepths();
		int[][] thatTDS=that.taxaWithDepths();
		if (!Arrays.equals(thisTDS[0], thatTDS[0])) return false;
		if (!Arrays.equals(thisTDS[1], thatTDS[1])) return false;
		return true;
	}

	@Override
	public String toString() {
		return "TaxaDist{" +
				"taxaWithRead=" + numberOfTaxaWithTag() +
				", totalDepth=" + totalDepth() +
				", maxTaxa=" + maxTaxa() +
				", "+taxaDepthMap().toString()+
				'}';
	}
}

/**
 * This is a specialized multiset for recording the distribution of a single across taxa.
 *
 * HashMap or Multiset or large arrays could be reasonable approaches by they do not scale well with hundreds of
 * taxa scored out of the thousands.
 * @author Ed Buckler
 */
class TaxaDistExpandable extends AbstractTaxaDistribution  {

	private ShortArrayList taxaWithTag;
	private TShortShortHashMap taxaTagMapShort = null;
	private TIntShortHashMap taxaTagMapInt = null;
	private int totalDepth;
	private final int maxTaxa;
	private boolean useIntTaxa = false;

	public TaxaDistExpandable(int maxTaxa) {
		this.maxTaxa=maxTaxa;
		if (maxTaxa > Short.MAX_VALUE) {
			useIntTaxa = true;
		}
		taxaWithTag = new ShortArrayList(1); // defaults to 10 items, we only want size of 1 initially
	}

	@Override
	public synchronized TaxaDistribution increment(int taxaNum) {
		if (totalDepth < 1000 ) {
			taxaWithTag.add((short)taxaNum);
		} else { // add previous values to map
			if (useIntTaxa) {
				if (taxaTagMapInt == null) {
					taxaTagMapInt = new TIntShortHashMap(maxTaxa);
				}
				for (int index = 0; index < taxaWithTag.size(); index++) {
					int tagIndex = taxaWithTag.get(index);
					// if not there, add it with value of 1.  Otherwise increment existing value by 1
					taxaTagMapInt.adjustOrPutValue(tagIndex,(short)1, (short)1);

				}
			} else {
				if (taxaTagMapShort == null) {
					taxaTagMapShort = new TShortShortHashMap(maxTaxa);
				}
				for (int index = 0; index < taxaWithTag.size(); index++) {
					short tagIndex = taxaWithTag.get(index);
					taxaTagMapShort.adjustOrPutValue(tagIndex, (short)1, (short)1);
				}
			}
			taxaWithTag.clear(); // no longer need the array
		}
		totalDepth++;
		return this;
	}

	@Override
	public int[] depths() {
		int[] depths=new int[maxTaxa];
		if (taxaWithTag.size() > 0) {
			for (int tn: taxaWithTag.elements()) {
				depths[tn]++;
			}
		} else {
			if (useIntTaxa) {
				for (TIntShortIterator ist = taxaTagMapInt.iterator(); ist.hasNext();) {
					ist.advance();
					depths[ist.key()] = ist.value();
				}
			} else {
				for (TShortShortIterator sst = taxaTagMapShort.iterator(); sst.hasNext();) {
					sst.advance();
					depths[sst.key()] = sst.value();
				}
			}
		}
		return depths;
	}

	@Override
	public int[][] taxaWithDepths() {
		int[] depths=depths();
		int countNoZero=0;
		for (int depth : depths) if(depth>0) countNoZero++;
		int[][] taxaDepth=new int[2][countNoZero];
		countNoZero=0;
		for (int i = 0; i <depths.length ; i++) {
			if(depths[i]>0) {
				taxaDepth[0][countNoZero]=i;
				taxaDepth[1][countNoZero]=depths[i];
				countNoZero++;
			}
		}
		return taxaDepth;
	}

	@Override
	public int totalDepth(){
		return totalDepth;
	}
	@Override
	public int maxTaxa() {
		return maxTaxa;
	}

	@Override
	public int numberOfTaxaWithTag() {
		return taxaWithDepths()[0].length;
	}

	@Override
	public int memorySize() {
		//minimal size 8 (object) + 12 (outer short array) + 12 (sizeArray) + 4+ 4 = 40
				int size=40;
		size+=(taxaWithTag.size()*12);
		if (taxaTagMapInt != null) {
			size+= taxaTagMapInt.size() * 12; // 8 + 4
		}
		if (taxaTagMapShort != null) {
			size+= taxaTagMapShort.size() * 8; // 4 + 4
		}
		return size;
	}

	private int unSignShort(short v) {
		if(v<0) return -v+1;
		return v;
	}
}


/**
 * This is a specialized multiset for recording the distribution of a single across taxa.
 *
 * HashMap or Multiset or large arrays could be reasonable approaches by they do not scale well with hundreds of
 * taxa scored out of the thousands.
 * @author Ed Buckler
 */
class TaxaDistSingleTaxon extends AbstractTaxaDistribution {
	//minimal size 8 + 4 + 4  = 16
	private final int maxTaxa;
	private final int singleTaxa;


	public TaxaDistSingleTaxon(int maxTaxa, int taxaWithTag) {
		this.maxTaxa=maxTaxa;
		singleTaxa=taxaWithTag;
	}

	@Override
	public synchronized TaxaDistribution increment(int taxaNum) {
		throw new UnsupportedOperationException("TaxaDistSingleTaxon cannot be increment.  Expand first.");
	}

	@Override
	public int[] depths() {
		int[] depths=new int[maxTaxa];
		depths[singleTaxa]=1;
		return depths;
	}

	@Override
	public int[][] taxaWithDepths() {
		return new int[][]{{singleTaxa},{1}};
	}

	@Override
	public int totalDepth(){
		return 1;
	}

	@Override
	public int maxTaxa() {
		return maxTaxa;
	}

	@Override
	public int memorySize() {
		return 16;
	}

	@Override
	public int numberOfTaxaWithTag() {
		return 1;
	}
}

/**
 * This is a specialized multiset for recording the distribution of a single across taxa.
 *
 * HashMap or Multiset or large arrays could be reasonable approaches by they do not scale well with hundreds of
 * taxa scored out of the thousands.
 * @author Ed Buckler
 */
class TaxaDistFixed extends AbstractTaxaDistribution  {
	//minimal size 8 + 12 + 12 + 10 + 12 + 4+ 4 = 66
	//This could be changed for the singletons by just making a new class
	private byte[] compTaxaSize;
	private int numTaxaWithTags;
	private int totalDepth;
	private final int maxTaxa;

	public TaxaDistFixed(int maxTaxa, int[] taxaWithTags, int[] depthOfTags) {
		this.maxTaxa=maxTaxa;
		numTaxaWithTags =taxaWithTags.length;
		totalDepth=0;
		for (int depthOfTag : depthOfTags) {totalDepth+=depthOfTag;}
		try{
			compTaxaSize=Snappy.compress(Ints.concat(taxaWithTags,depthOfTags));
			//System.out.printf("tn:%d\tdepth:%d\ts:%d%n", numTaxaWithTags,totalDepth,compTaxaSize.length);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	@Override
	public synchronized TaxaDistribution increment(int taxaNum) {
		throw new UnsupportedOperationException("TaxaDistFixed cannot be increment.  Change to expandable first first.");
	}

	@Override
	public int[] depths() {
		try {
			int[] depths = new int[maxTaxa];
			int[] taxaSize = Snappy.uncompressIntArray(compTaxaSize);
			for (int i = 0; i < numTaxaWithTags; i++) {
				depths[taxaSize[i]]=taxaSize[i+ numTaxaWithTags];
			}
			return depths;
		}  catch (IOException e) {
			e.printStackTrace();
			return null;
		}
	}

	@Override
	public int[][] taxaWithDepths() {

		try {
			int[][] taxaDepth=new int[2][numTaxaWithTags];
			int[] taxaSize = Snappy.uncompressIntArray(compTaxaSize);
			taxaDepth[0]=Arrays.copyOf(taxaSize, numTaxaWithTags);
			taxaDepth[1]=Arrays.copyOfRange(taxaSize, numTaxaWithTags,taxaSize.length);
			return taxaDepth;
		}  catch (IOException e) {
			e.printStackTrace();
			return null;
		}
	}

	@Override
	public int totalDepth(){
		return totalDepth;
	}

	@Override
	public int maxTaxa() {
		return maxTaxa;
	}

	@Override
	public int numberOfTaxaWithTag() {
		return numTaxaWithTags;
	}

	@Override
	public int memorySize() {
		//minimal size 8 (object) + 12 (sizeArray) + 4+ 4 = 40
				return 28+compTaxaSize.length;
	}
}