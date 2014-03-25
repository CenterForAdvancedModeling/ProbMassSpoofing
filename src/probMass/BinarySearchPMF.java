package probMass;


import java.io.Serializable;
import java.util.Arrays;


/**
 * A BinarySearchPMF accelerates random draws by precomputing the cumulative mass function of a
 * probability mass function. The precomputed CMF is used to quickly transform a uniform random
 * number to an array index by performing a binary search for the input uniform random number.
 */
public class BinarySearchPMF implements ProbMassFunction, Serializable {

	/** The CMF of the input distribution. */
	private double[] cummMassFun;


	/**
	 * Build an object that can quickly (and repeatedly) draw samples from the given set of weights
	 *
	 * @param weights - A set of weights, each number will be drawn with prob (wieght[i] / sum)
	 */
	public BinarySearchPMF(double[] weights) {
		
		this.cummMassFun = Util.buildCMF(weights);

		assert(cummMassFun[cummMassFun.length - 1] == 1) : "The Final entry in the CMF should be 1";
	}


	/**
	 * @param uniformDraw - A uniformly distributed random number between 0 and 1
	 *
	 * @return A Random number between 0 and (orginalDist.length - 1)
	 */
	@Override
	public int getSample(double uniformDraw) {
		
		int val = Arrays.binarySearch(cummMassFun, uniformDraw);

		if (val < 0) {
			return -val - 1;
		} else {
			return val;
		}
	}
}
