package probMass;


import java.io.Serializable;
import java.text.DecimalFormat;
import java.util.Arrays;


public class SpoofingPMF implements ProbMassFunction, Serializable {

	/** The number of entries in p. */
	private int N;

	/** The set of bins. */
	private Bin[] bins;

	/** The Cumulative Mass Function built from the Bins. */
	private double[] binCMF;

	/** The hint table built for the Bin_CMF. */
	private int[] hintTable;


	/**
	 * Build a SpoofingPMF.
	 *
	 * @param lifetime - How many random draws this Spoofing distribution should support.
	 * @param weights - A set of unsorted weights.
	 */
	public SpoofingPMF(long lifetime, double[] weights) {
		
		Util.checkPMFInputArray(weights);

		this.N = weights.length;
		this.bins = Bin.buildBins(lifetime, weights);
		buildCMFandHints();
	}


	/** Build the Cumulative Mass Function and Hint table for the set of Bins. */
	private void buildCMFandHints() {
		
		this.binCMF = new double[this.bins.length];
		binCMF[0] = bins[0].getPSum();
		for (int i = 1; i < bins.length; i++) {
			binCMF[i] = binCMF[i - 1] + bins[i].getPSum();
		}

		this.hintTable = new int[this.bins.length];

		double n = (double) this.bins.length;
		for (int i = 0; i < hintTable.length; i++) {
			hintTable[i] = getBinHint(((double) i) / n);
		}

//		System.out.println("Total Number of Bins :: " + this.bins.length);
	}


	/**
	 * @param uniformDraw - A uniformly distributed random number between 0 and 1
	 *
	 * @return A Random number between 0 and (orginalDist.length - 1)
	 */
	@Override
	public int getSample(double uniformDraw) {

		//get the hint
		int currentNum = hintTable[(int) (hintTable.length * uniformDraw)];

		//walk up the table until you are done
		while (binCMF[currentNum] < uniformDraw) {
			currentNum++;
		}

		//transform the uniformDraw into another U(0,1) quantity
		double binMin = (currentNum > 0) ? binCMF[currentNum - 1] : 0;
		double binMax = binCMF[currentNum];
		double newUniform = (uniformDraw - binMin) / (binMax - binMin);

		return bins[currentNum].sample(newUniform);
	}


	/** Generate a Bin hint for a specific uniform draw. */
	private int getBinHint(double point) {
		int val = Arrays.binarySearch(binCMF, point);

		if (val < 0) {
			return -val - 1;
		} else {
			return val;
		}
	}


	/** Use this table to understand what has been built. */
	public void printBinInformationTable() {
		System.out.println("CMF\t\tinBinP\t\tHint\t\ti");
		DecimalFormat df = new DecimalFormat("#.####");
		for (int i = 0; i < bins.length; i++) {
			System.out.println(
					df.format(binCMF[i])
					+ "\t\t" + df.format(bins[i].getPSum())
					+ "\t\t" + df.format(hintTable[i])
					+ "\t\t" + i);
		}
	}
}
