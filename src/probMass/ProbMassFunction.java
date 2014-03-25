package probMass;


import java.io.Serializable;


/**
 * A ProbMassFunction enables the user to transform a uniform random number drawn from a standard
 * RNG to a sample from a Probability Mass Function.
 */
public interface ProbMassFunction extends Serializable {

	/**
	 * @param uniformRandomDraw - A uniformly distributed random number between 0 and 1
	 *
	 * @return - A sample from the PMF that corresponds to the input uniform random number.
	 */
	public int getSample(double uniformRandomDraw);
}
