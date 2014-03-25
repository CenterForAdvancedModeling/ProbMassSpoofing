package probMass;


import java.io.Serializable;
import java.util.Arrays;


/**
 * A ChanAsuaPMF accelerates random draws by precomputing the cumulative mass function of a
 * probability mass function and a set of hints that will frequently lead directly to the correct
 * random sample. See CHAN, H. C. AND ASAU, Y. 1974. "On generating random variates from an
 * empirical distribution" . IIE Transactions, Vol. 6 No. 2, 163 - 166.
 *
 * This implementation uses 1 hint for every entry in the input array of weights.
 */
public class ChanAsuaPMF implements ProbMassFunction, Serializable {

	/** The CMF of the input distribution. */
	private double[] cummulativeMassFunction;

	/** A series of hints that will speed up taking samples. */
	private int[] hintTable;


	/**
	 * Build an object that can quickly (and repeatedly) draw samples from the given set of weights
	 *
	 * @param weights - A set of weights, each number will be drawn with probability (weight[i] /
	 * sumOfWeights)
	 */
	public ChanAsuaPMF(double[] weights) {
		
		this.cummulativeMassFunction = Util.buildCMF(weights);
		
		//build hintTable
		this.hintTable = new int[weights.length];
		double n = (double) weights.length;
		for (int i = 0; i < hintTable.length; i++) {
			hintTable[i] = findHint(((double) i) / n);
		}
	}


	/** Use a binary search to find the unaccelerated sample. */
	private int findHint(double point) {
		int val = Arrays.binarySearch(cummulativeMassFunction, point);

		if (val < 0) {
			return -val - 1;
		} else {
			return val;
		}
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
		while (cummulativeMassFunction[currentNum] < uniformDraw) {
			currentNum++;
		}

		return currentNum;
	}
}
