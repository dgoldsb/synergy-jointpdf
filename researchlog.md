# Research log

## Saturday 2017/04/15

### NPEET

* Added the NPEET submodule for KNN mutual information

### resilientpdf.py

* Implemented a class that defines a PDF
* Implemented a class that samples from the PDF, for KNN later on
* Implemented a simple cost function that evaluates the inverse of $I(x^{t}:x^{t+1}_\mathrm{nudged})$ and normalizes with the entropy
    * Uses Runge-Kutta 8 for evaluation
    * ODE parameters are stored in a vector
    * PDFs are nudged by shifting the mean
    * The cost function does stuff it should not do: the mutual information becomes negative and small
* Wrote optimizer, use scipy.minimize and scipy.differential_evolution
    * Optimizes the ODE vector
    * Resamples every training iteration
* The cost landscape is changing too quickly, I think I should train on one sample with set nudges, then continue training a few times with new data
    * If this does not work, the cost function is the problem

## Sunday 2017/04/16

### resilientpdf.py

* Altered the optimizing scheme
    * Samples are now drawn before optimizing only, this makes the cost-landscape easier to traverse
    * Several rounds of optimizing can be done, each time with a different initial sample
* Cost function needs updating
    * Removed entropy normalization, this does not work in the continuous case, as entropy gets a different meaning here (not all assumptions apply)
    * I need to include a separate part for the memory and the resilience
    * I need to normalize both, potentially by using $I(X^{t+1}:X^{t+1})$
* Did a bit of testing
    * Minimize converges nicely for trivial example (one variable, good starting guess)
    * Minimize converges nicely for trivial example (one variable, bad starting guess)
    * Cost function is still not functioning nicely in less trivial examples, maybe I should look into the Jaccard distance

## Friday 2017/04/21

### resilientpdf.py

* Altered the nudging technique: instead of shifting the mean, we just add a normal distribution, making sure that we don't train with a bias
    * Now we train to resist a nudge, not to compensate one specific nudge (which would not result in resilience)
    * Essentially we now sample from a different, noisier normal distribution (adding two normal distributions)
* Added option to take out self-loop
* Changed cost function
    * Compare MI nudged and unnudged in $t + \delta t$
* Added basinhopping after all (untested)
* I give KNN a list of lists now, so it evaluates the entire system at once

## Sunday 2017/04/23

### resilientpdf.py

* Added improved logging to a debug file
* Undid KNN list of lists changed, I now take the max of the costs
    * The problem was that with 2 variables the maximum k is 0
* The results are pretty bad after these changes: the cost function is not performing at all
    * Makes sense: memory is thrown out as is, and resilience is done wrong
    * I should do memory with the MI between t and t + dt of the unnudged state
    * I should do resilience with the MSE of the t + dt nudged and unnudged state
    * The problem is: how do I normalize?

### setup.sh

* Created a setup script
    * Pulls latest commit in master from submodules
    * Touches the \_\_init\_\_.py files where necessarry

## Tuesday 2017/04/25

### resilientpdf.py

* Worked on cost function
    * Added Kullback-Leibler divergence to measure nudge impact
    * Went back to the list-of-lists approach, but now properly (transposed the sampleset)
    * I now take the mutual information of the starting state and the unnudged endstate for memory measurement

## Wednesday 2017/04/26

### resilientpdf.py

* Finished histogram plotting functions
* Added visual debugging mode, with levels 0 (none), 1 (regular output) and 2 (every iteration)
* Fixed NPEET implementation of KL-divergence (I was copying by reference)
* Improved logger: I now have a stream with everything to the file, and display only INFO and above
* I think MI is not great for memory, as it ignores that biological systems return to the same state. I think I need to make the following changes:
    1. Use the KL-divergence to determine memory
        1. This assumes that dt is in arbitrary time units
        2. This assumes that a natural system always returns to a value at some point
    2. Discuss with Jaap Kaandorp what induces limit cycles
    3. Brainstorm about the definition of memory: the KL-divergence seems to ignore the possibility off limit cycles
* Replaced Nelder-Mead with BFGS (recommendation Robert Hana)
* Preliminary BFGS results are amazing, it even went from 0.5 with a weird parameter (outlier) to 16, only to go to 0.08 with better parameters

## Thursday 2017/04/27

### resilientpdf.py

* Added ODE plotting library
* Marked annotate (pyplot) for plotting the system
    * Tried to apply this, found something to use: arrow.py
* Plotting the results shows training does not have the intended effect, using a BFGS run with good convergence

## Tuesday 2017/05/02

### resilientpdf.py

* Not doing arrow.py, because a linear ODE is too simple to approximate over time anyway (it is a linear, local in time approximation)
* Back to MI for memory, it is the correct way to do it, and linear ODEs cannot have steady states
* Added scatterplot for MI

<<<<<<< HEAD
## Thursday 2017/05/11

### setup.sh

* Added .mat file extraction from the Yeast model
=======
## Thursday 2017/05/18

(Had a short hiatus due to personal problems)

### resilientpdf.py

* Follow the advice (glucose is a good place to start)
    1. Start from MatLab files from the paper
    2. See if I can build a config to enter these simple systems in my framework
    3. Enter and answer: is the M and R significantly better than a random system (**make this experiment one in my experiments.py**)
    4. See if it can be further optimized
    5. Make a plot of 3, preferably a 2D/3D plot if it has 2/3 dimensions with a dot for our hopefully outlier
>>>>>>> Forgotten commit

## Next time

* Add configuration file library (JSON), add a save/load function
* Follow the advice (glucose is a good place to start)
    1. Start from MatLab files from the paper
    2. See if I can build a config to enter these simple systems in my framework
    3. Enter and answer: is the M and R significantly better than a random system (**make this experiment one in my experiments.py**)
    4. See if it can be further optimized
    5. Make a plot of 3, preferably a 2D/3D plot if it has 2/3 dimensions with a dot for our hopefully outlier
* Plot and pickle training process
* Add training curve, with both components and a cumulative
* Add MI profile
* Add numerical synergy analysis, exact for the case with 2 systems
* Verken CUDA just in time versie van NPEET
* Ook amount syn met random vergelijken
* Markov systeem deterministic, mutual information meten nu en later is dan wat gek