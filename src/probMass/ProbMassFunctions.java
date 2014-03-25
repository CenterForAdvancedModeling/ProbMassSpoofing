package probMass;


public class ProbMassFunctions {

	/**
	 * Create a ProbMassFunction that uses 12 bytes per entry in the input array and should return a
	 * random sample in constant time.
	 *
	 * @param weights - A set of weights, each number will be drawn with probability (weight[i] /
	 * sumOfWeights)
	 */
	public static ProbMassFunction highSpeedHighMemoryPMF(double[] weights) {
		return new ChanAsuaPMF(weights);
	}


	/**
	 * Create a ProbMassFunction that uses 8 bytes per entry in the input array and returns a sample
	 * in logarithmic time.
	 *
	 * @param weights - A set of weights, each number will be drawn with probability (weight[i] /
	 * sumOfWeights)
	 */
	public static ProbMassFunction mediumSpeedMediumMemoryPMF(double[] weights) {
		return new BinarySearchPMF(weights);
	}


	/**
	 * Create a ProbMassFunction based on the lossy "Spoofing" compression technique. This
	 * ProbMassFunction uses roughly 2 bytes per entry in the input array and should return a random
	 * sample almost as fast as the highSpeedHighMemoryPMF.
	 *
	 * @param weights - A set of weights, each number will be drawn with probability (weight[i] /
	 * sumOfWeights)
	 */
	public static ProbMassFunction compressedPMF(long lifetime, double[] weights) {
		return new SpoofingPMF(lifetime, weights);
	}
}
