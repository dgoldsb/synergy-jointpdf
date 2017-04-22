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

* X Altered the nudging technique: instead of shifting the mean, we just add a normal distribution, making sure that we don't train with a bias
    * Now we train to resist a nudge, not to compensate one specific nudge (which would not result in resilience)
    * Essentially we now sample from a different, noisier normal distribution (adding two normal distributions)
* Jaccard
* x Option to take out self-loop
* Changed cost function
    * Compare MI nudged and unnudged in $t + \delta t$
* X Toch basinhopping als optie toevoegen
* x Geef list of lists aan KNN
* Add MSE term

## Saturday 2017/04/22

### resilientpdf.py

* Logging to debug file
* Pickling resulting system with a good dataset generated with these params
* Store optimization results over time in a log pickle