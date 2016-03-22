package net.maizegenetics.analysis.association;

import java.util.concurrent.BlockingQueue;

import net.maizegenetics.util.BitSet;

public class DiscreteSiteOlsSolver implements Runnable {
	BlockingQueue myQueue;
	double[] myData;
	BitSet missingObs;
	
	public DiscreteSiteOlsSolver(BlockingQueue<Integer> inputQueue, double[] data, BitSet missing) {
		myQueue = inputQueue;
		myData = data;
		missingObs = missing;
	}
	
	@Override
	public void run() {
		// TODO Auto-generated method stub
		
	}
	
	
}
