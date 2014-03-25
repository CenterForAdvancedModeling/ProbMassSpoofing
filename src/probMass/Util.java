package probMass;


import com.google.common.primitives.Doubles;


class Util {

	/**
	 * @param anArray - An array of doubles.
	 *
	 * @return the sum of an array.
	 */
	static double sum(double[] anArray) {
		double sum = 0;
		for (int i = 0; i < anArray.length; i++) {
			sum += anArray[i];
		}
		return sum;
	}


	/**
	 * Check an array of "weights" that should define a probability mass function.
	 *
	 * @param weights - Each entry in this array must be finite and non-negative. Also, the sum of
	 * the entries in this array must be positive.
	 */
	static void checkPMFInputArray(double[] weights) {

		for (double weight : weights) {

			if (!Doubles.isFinite(weight)) {
				throw new IllegalArgumentException("Non-Finite input weight :: " + weight);
			}

			if (weight < 0) {
				throw new IllegalArgumentException("All Weights Must be non-negative :: " + weight);
			}
		}
	}


	/** Build a Cumulative Mass Function for the given weights (which will be normalized). */
	static double[] buildCMF(double[] weights) {

		checkPMFInputArray(weights);
		
		double sum = Util.sum(weights);

		double[] cmf = new double[weights.length];

		cmf[0] = weights[0];

		for (int i = 1; i < weights.length; i++) {
			cmf[i] = cmf[i - 1] + weights[i];
		}
		for (int i = 0; i < cmf.length; i++) {
			cmf[i] /= sum;
		}

		return cmf;
	}
}
