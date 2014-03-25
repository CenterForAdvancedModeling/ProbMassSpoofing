package demo;


import com.google.common.collect.Multiset;
import com.google.common.collect.TreeMultiset;
import java.util.Arrays;
import java.util.Random;
import probMass.BinarySearchPMF;
import probMass.ChanAsuaPMF;
import probMass.ProbMassFunction;
import probMass.SpoofingPMF;


public class Demonstrate {

	public static void main(String[] args) {

		int N = 2000;	//how large should the distribution be.		
		boolean showSamples = false;

		System.out.println("\n\nStarting curved test");
		groupTest(Weights.gentleCurve(N), showSamples);

		System.out.println("\n\nStarting sloped test");
		groupTest(Weights.singleLine(N), showSamples);

		System.out.println("\n\nStarting high-low test");
		groupTest(Weights.alternating(N), showSamples);
	}


	/** Test several different distributions. */
	private static void groupTest(double[] weights, boolean showSamples) {
		int NUM_DRAWS = 100000000;	//100 Mil
		singleTest(NUM_DRAWS, new BinarySearchPMF(weights), showSamples);
		singleTest(NUM_DRAWS, new ChanAsuaPMF(weights), showSamples);
		singleTest(NUM_DRAWS, new SpoofingPMF((long) NUM_DRAWS, weights), showSamples);
	}


	/**
	 * Time, and show output (possible) from a ProbMassFunction
	 *
	 * @param numDraws - Compute this many random draws
	 * @param dist - The distribution to use
	 * @param showSamples - True if you want to see the samples printed (which is probably
	 * false when a distributions size is huge)
	 */
	private static void singleTest(int numDraws, ProbMassFunction dist, boolean showSamples) {

		long startTime = System.nanoTime();

		//find the samples
		Multiset<Integer> multiset = TreeMultiset.create();
		Random rand = new Random(17L);
		for (int i = 0; i < numDraws; i++) {
			int sample = dist.getSample(rand.nextDouble());
			multiset.add(sample);
		}

		long endTime = System.nanoTime();

		System.out.println("Time(ns) :: " + (endTime - startTime));

		if (showSamples) {
			for (Multiset.Entry<Integer> entry : multiset.entrySet()) {
				System.out.println(entry.getElement() + " appeared " + entry.getCount() + " times");
			}
		}
	}

	static class Weights {

		/**
		 * Create many random weights. Where weight_i = 1.0 / rand.nextDouble()
		 *
		 * @param N - The number of weights simulated
		 *
		 * @return - A series of weights
		 */
		private static double[] gentleCurve(int N) {
			Random rand = new Random(17L);
			double[] weights = new double[N];
			for (int i = 0; i < weights.length; i++) {
				weights[i] = 1.0 / rand.nextDouble();
			}
			Arrays.sort(weights);

			return weights;
		}


		/** @return - A series of weights that form a line. */
		private static double[] singleLine(int N) {
			double[] weights = new double[N];
			for (int i = 0; i < weights.length; i++) {
				weights[i] = 20.0 + i;
			}

			return weights;
		}


		/** @return - A series of weights that alternate between high and low values. */
		private static double[] alternating(int N) {
			double[] weights = new double[N];
			for (int i = 0; i < weights.length; i++) {
				weights[i] = 5 * (i % 2) + 1;
			}

			return weights;
		}
	}
}
