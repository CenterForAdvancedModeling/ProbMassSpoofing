package probMass;


import java.io.Serializable;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedList;


class Bin implements Comparable<Bin>, Serializable {

	/** This variable is nulled out when bin refinement is complete. */
	private double[] pValues;

	/** The number of entries in this bin. */
	private int n;

	/** (n+1)/2 is stored rather than frequently recomputed */
	private double nPlusOneOver2;

	/** The total probability in this bin (needed to compute ChanAsua Distribution). */
	private final double pSum;

	/** The y intercept (kinda) of the interpolating line. */
	private double height;

	/** The slope of the interpolating line. */
	private double slope;

	/** The error in this bins interpolation. */
	private double inBinLambda;

	/** Retains the original order (used in place of this.left when input is unsorted). */
	private short[] entries;

	/** The index of the leftMost value (used in place of this.entries when input is presorted). */
	private short leftMost;


	/**
	 * Build a set of bin that supports spoofing for these weights.
	 * 
	 * @param lifetime - How many draws a Spoofing Distribution should support
	 * @param weights - The initial unsorted weights
	 * 
	 * @return - An acceptable set of Bins
	 */
	static Bin[] buildBins(long lifetime, double[] weights) {

		//lambda = (Phi(alpha)^2)/(2*lifetime)
		//assume error tolerance = .05
		//@todo -- allow error tolerance to be assigned
		//note: 1.6449 num standard deviations when there is .05 error in left and right tail
		double ERROR_LIMIT = Math.pow(1.6449, 2) / (double) (2.0 * lifetime);
//		System.out.println("ERROR_LIMIT :: " + ERROR_LIMIT);

		//ensure all weights are > 0
		for (int i = 0; i < weights.length; i++) {
			if (weights[i] < 0) {
				throw new IllegalArgumentException("Negative Weights not allowed");
			}
		}
		double wSum = Util.sum(weights);

		//see if the input weights are presorted (if so, do not store entries)
		boolean inputWasSorted = true;
		for (int i = 1; i < weights.length; i++) {
			if(weights[i] > weights[i-1]) {
				inputWasSorted = false;
				break;
			}
		}

		Bin initalBin = null;
		
		if(inputWasSorted) {

			System.out.println("Input was presorted");

			//normalize input
			double[] pValues = new double[weights.length];
			for (int i = 0; i < pValues.length; i++) {
				pValues[i] = weights[i]/wSum;
			}

			//create first bin -- don't use entries because the input was ordered correctly
			initalBin = new Bin(pValues, (short) 0);

		} else {
			
			//do inital sorting WHILE KEEPING TRACK OF INITAL ORDER
			Sorter[] sortMe = new Sorter[weights.length];
			for (int i = 0; i < sortMe.length; i++) {
				sortMe[i] = new Sorter(weights[i] / wSum, (short) i);
			}
			Arrays.sort(sortMe);

			double[] pValues = new double[sortMe.length];
			short[] entries = new short[sortMe.length];
			for (int i = 0; i < sortMe.length; i++) {
				pValues[i] = sortMe[i].pValue;
				entries[i] = sortMe[i].entry;
			}
			initalBin = new Bin(pValues, entries);
		}

		LinkedList<Bin> binList = new LinkedList<>();
		binList.add(initalBin);

		double totalError = initalBin.inBinLambda;

		while (totalError > ERROR_LIMIT) {

			Collections.sort(binList);
			Bin highErrorBin = binList.removeFirst();

			Bin[] children = highErrorBin.split();

			for (int i = 0; i < children.length; i++) {
				binList.addLast(children[i]);
			}

			totalError = 0;
			double totalSize = 0;
			for (Bin bin : binList) {
				totalError += bin.inBinLambda;
				totalSize += bin.pSum;
			}

//			System.out.println("totalError :: " + totalError + "\n");
//			System.out.println("totalSize :: " + totalSize);
		}

		//output the bins
		Bin[] bins = binList.toArray(new Bin[0]);

		//remove the unnecessary information to save memory
		for (int i = 0; i < bins.length; i++) {
			bins[i].flatten();
		}
		System.out.println("Final Bin Count :: " + bins.length);

		return bins;
	}


	/**
	 * Build a bin -- automatically apply an interpolation
	 * 
	 * @param pValues - The probabilities being interpolated
	 * @param entries - The numbers that must be returned when sampling
	 */
	private Bin(double[] pValues, short[] entries) {

		for (int i = 1; i < pValues.length; i++) {
			if (pValues[i] > pValues[i - 1]) {
				throw new IllegalArgumentException("pValues must be sorted in descending order");
			}
		}

		this.pValues = pValues;
		this.n = pValues.length;
		this.nPlusOneOver2 = (n + 1.0) / 2.0;
		this.pSum = Util.sum(pValues);
		this.height = pSum / ((double) n);
		setSlope();
		this.inBinLambda = computeError();
		this.entries = entries;

//		System.out.println(
//				"Made bin with error :: " + inBinLambda +
//				"\tsize :: " + n +
//				"\t pSum :: " + pSum +
//				"\t slope :: " + slope);
	}


	/**
	 * Build a bin -- automatically apply an interpolation.  This constructor is intended for use
	 * when the input distribution is already sorted.
	 *
	 * @param pValues - The probabilities being interpolated
	 * @param leftMost - The index of the left most (greatest) pValue
	 */
	private Bin(double[] pValues, int leftMost) {
		
		for (int i = 1; i < pValues.length; i++) {
			if (pValues[i] > pValues[i - 1]) {
				throw new IllegalArgumentException("pValues must be sorted in descending order");
			}
		}

		this.pValues = pValues;
		this.n = pValues.length;
		this.nPlusOneOver2 = (n + 1.0) / 2.0;
		this.pSum = Util.sum(pValues);
		this.height = pSum / ((double) n);
		setSlope();
		this.inBinLambda = computeError();
		this.entries = null;
		this.leftMost = (short) leftMost;

//		System.out.println(
//				"Made bin with error :: " + inBinLambda +
//				"\tsize :: " + n +
//				"\t pSum :: " + pSum +
//				"\t slope :: " + slope);
		
	}

	/** @return - The sum of all pValues in this Bin. */
	double getPSum() {
		return pSum;
	}


	/**
	 * Find, and set, the slope the minimizes inBinLambda
	 * 
	 * Parker and Engler derive the forumala used herin to compute the best possible beta.
	 * 
	 * This method can implement the bisection method (which may be slightly more accurate)
	 * or use the formula for Beta shown in Parker and Engler.
	 */
	private void setSlope() {

		//when a bin has size = 2 we can directly solve
		if (n == 2) {
			this.slope = (pValues[0] - this.height) * (-2.0);

//			System.out.println("slope :: " + slope);
//			
//			for (int i = 0; i < pValues.length; i++) {
//				System.out.println("p[" + i + "] :: " + pValues[i]);
//			}
//			for (int i = 0; i < pValues.length; i++) {
//				System.out.println("q[" + i + "] :: " + this.computeQ(i));
//			}
			return;
		}

		//use these values for the inital points in the binary search
		double lowSlope = 0;	//this slope produces a flat bin
		double highSlope = -height / ((n + 1.0) / 2.0);	//this slope makes q_n = 0

		this.slope = lowSlope;
		double lowSlopeError = this.computeError();
		this.slope = highSlope;
		double highSlopeError = this.computeError();

		
		//if inital error is basically numeric error -- end now
		if (lowSlopeError <= 1.0e-25) {
			this.slope = lowSlope;
			return;
		}

		//if inital error is basically numeric error -- end now
		if (highSlopeError <= 1.0e-25) {
			this.slope = highSlope;
			return;
		}


		//while optimal slope is still changing
		while (Math.abs((lowSlope - highSlope) / Math.min(lowSlope, highSlope)) > 0.00000001) {

			if (lowSlopeError < highSlopeError) {
				//low slope is better -- reset highSlope, recompute its error
				highSlope = (lowSlope + highSlope) / 2.0;
				this.slope = highSlope;
				highSlopeError = this.computeError();
			} else {
				//high slope is better -- reset lowSlope, recompute its error
				lowSlope = (lowSlope + highSlope) / 2.0;
				this.slope = lowSlope;
				lowSlopeError = this.computeError();
			}

//			System.out.println("low :: " + lowSlope);
//			System.out.println("high :: " + highSlope);
		}

		this.slope = (lowSlope + highSlope) / 2.0;
	}


	/**
	 * Parker and Engler show how to compute error
	 * 
	 * @return - The sum of all lambda_i values --> aka 
	 */
	private double computeError() {

		double[] errors = new double[n];
		double thisError;
		for (int i = 0; i < errors.length; i++) {
			double q = this.computeQ(i);
			thisError = pValues[i] / q - 1.0;
			errors[i] = q * thisError * thisError;
		}

		return Util.sum(errors);
	}


	/**
	 * Return a properly distributed original entry from this bin
	 * 
	 * @param uniformDraw - A uniform random number between 0 and 1
	 * @return - An entry in this bin
	 */
	int sample(double uniformDraw) {

		//when Bin is exact don't use stacking
		if (n == 2) {
			double q1 = computeQ(0);
			double q2 = computeQ(1);

			if (uniformDraw < q1 / (q1 + q2)) {
				return (this.entries != null) ? (this.entries[0]) : this.leftMost;
				//return this.entries[0];
			} else {
				return (this.entries != null) ? (this.entries[1]) : this.leftMost + 1;
				//return this.entries[1];
			}
		}

		//when Bin is flat return results directly
		if (this.slope == 0) {
			int s = (int) (uniformDraw * n);
			return (this.entries != null) ? (this.entries[s]) : this.leftMost + s;
//			return this.entries[(int) (uniformDraw * entries.length)];
		}

		//no special case found....		
		int columnNum = (int) (uniformDraw * (n / 2.0));	//btw 0 and floor(n/2)

		int candidate1 = columnNum;
		int candidate2 = n - columnNum - 1;

		double q1 = computeQ(candidate1);
		double q2 = computeQ(candidate2);

		double minDraw = (2.0 * columnNum) / n;
		double maxDraw = Math.min(1.0, 2.0 * (((double) columnNum) + 1.0) / n);
		//the above "min" hands odd cases where the maxDraw might be greater than 1
		//this might not be necessary -- not important enough to figure out

//		System.out.println("colNum :: " + columnNum);
//		System.out.println("1 :: " + candidate1);
//		System.out.println("2 :: " + candidate2);
//		System.out.println("minDraw :: " + minDraw);
//		System.out.println("maxDraw :: " + maxDraw);

		double reRandom = (maxDraw - uniformDraw) / (maxDraw - minDraw);

		if (reRandom < q1 / (q1 + q2)) {
			return (this.entries != null) ? (this.entries[candidate1]) : this.leftMost + candidate1;
			//return this.entries[candidate1];
		} else {
			return (this.entries != null) ? (this.entries[candidate2]) : this.leftMost + candidate2;
			//return this.entries[candidate2];
		}
	}


	/** Use the formula for Q to compute an entries's prob. */
	private double computeQ(int index) {
		return height + slope * (index + 1 - nPlusOneOver2);
	}


	/**
	 * Split this Bin in "half".  Ensure the bin with large p_i entries has a size that is 
	 * divisible by two.  This constraint ensures that bins with an odd number of entries (i.e. 
	 * bins that cannot be "split until perfect") are only created with the smallest p_i values.
	 * 
	 * @return - Two bins that represent the original bin after a split.
	 */
	Bin[] split() {

//		DecimalFormat df = new DecimalFormat("#.####");
//		System.out.println(
//				"Split --> " +
//				"error " + df.format(this.inBinLambda) +
//				"   size " + n + "   pSum " + df.format(this.pSum));

		if (this.n <= 3) {
			throw new IllegalStateException(
					"Cannot split a bin with 3 or fewer entries, current size is :: " + n);
		}

		int index = this.n / 2;
		if (index % 2 == 1) {
			index--;
		}

		double[] frontPs = Arrays.copyOfRange(this.pValues, 0, index);
		double[] backPs = Arrays.copyOfRange(this.pValues, index, n);

		if (this.entries != null) {
			
			short[] frontEntries = Arrays.copyOfRange(this.entries, 0, index);
			short[] backEntries = Arrays.copyOfRange(this.entries, index, n);

			return new Bin[]{
						new Bin(frontPs, frontEntries),
						new Bin(backPs, backEntries)
					};
		} else {
			return new Bin[]{
						new Bin(frontPs, this.leftMost),
						new Bin(backPs, this.leftMost + frontPs.length)
					};
		}
	}


	/** Null out the pValues variable, once the approximation is complete pValues in deadWeight. */
	private void flatten() {
		this.pValues = null;
	}


	/** 
	 * Sort by inBinLambda (in descending order).
	 * 
	 * @param other - A second bin
	 */
	@Override
	public int compareTo(Bin other) {
		if (this.inBinLambda < other.inBinLambda) {
			return 1;
		} else if (this.inBinLambda > other.inBinLambda) {
			return -1;
		} else {
			return 0;
		}
	}


	/** This simple utility class sorts the inital weights and retains the inital order. */
	private static class Sorter implements Comparable<Sorter> {

		double pValue;

		short entry;


		Sorter(double pValue, short entry) {
			this.pValue = pValue;
			this.entry = entry;
		}


		/** Sort in descending order. */
		@Override
		public int compareTo(Sorter other) {
			if (this.pValue < other.pValue) {
				return 1;
			} else if (this.pValue > other.pValue) {
				return -1;
			} else {
				return 0;
			}
		}
	}
}
