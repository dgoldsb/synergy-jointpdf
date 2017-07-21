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

## Thursday 2017/05/11

### setup.sh

* Added .mat file extraction from the Yeast model

## Thursday 2017/05/18

* (Had a short hiatus due to personal problems)
* Send email to my supervisor
    * How the MI is redundant for deterministic systems
    * I can't open with MatLab
    * Proposal: sit with Gorkhan to address the issue, find a more concrete network I can work with?

## Friday 2017/05/19

* Meeting met Manon was goed, er kwam best veel naar boven uit mezelf qua ideeeen!
    * Afgesproken om op vaste dagen te werken aan mijn thesis, die dagen verder vrij van klussen te houden
    * Aan het einde van die dag stuur ik haar een update
    * Ik houd de planning op SCRUM, dus vrij korte termijn, als ik een lange termijn planning heb komt het niet goed zodra er iets wegvalt
* Meeting met Rick: gebruik van MI op deterministische systemen is geen probleem
    * Wel aangeven waarom dit OK is in het paper

## resilientpdf.py

* Add configuration file library (JSON), add a save/load function
* Plot and pickle training process

# Saturday 2017/05/20

* Put on Godspeed You! Black Emperor
* Did my standup, kicked some ass

## hypothesis_testing.py

* Created for later, when I can start testing things
* Thought of a logical way to store results

## resilientpdf.py

* Finished the JSON configuration saving and loading
* Optimized the variables in the System object
* Rearranged the order to make more sense
* Added a MI profile visualization
* Added a plotting function for the training process, with a distinct color for each cycle

# Thursday 2017/05/25

* Read the paper by Peixoto, as suggested by Rick
    * Quite depressing how much smarter the approach of this paper is
    * See my notes, it uses Boolean networks instead of my continuous approach
    * I should look at the yeast network, and then discuss with Rick
* Got a message from Gokhan, see issue #11
* Synergy is still my focus, which is different from the Peixoto paper, I give myself 2 weeks to see if continuous gives me the results I need, otherwise I will discuss a switch to discrete Boolean systems with Rick

# Friday 2017/05/25

* Made the project for the writing part, added the todos as notes
* I am so fucking lost in this MatLab model

# Saturday 2017/05/25

* Email Gokhan
    * I don't really get the structure
    * I don't know what to look for
    * MatLab is infuriating
    * Explain more what I had in mind, is this in there?
    * Should I adapt my model? What is the standard?
    * RBN better?
* Slack to Rick, ask him more about if RBN is better suited
* CC both to Manon, explain how badly it went
* Wrote on my fucking thesis and fuck programming for today

**From now on I just make my research notes on paper, more efficient for me...**