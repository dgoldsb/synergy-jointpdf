# Implemented to-do items

* Dyadic versus polydyadic, I can say that my correlations are dyadic
* Make sure the reader understands why you are telling something in the literature review from what he read in the introduction
* Predictor/Predicted as common terms
* Terugkomen op: In this work we aim to examine the links between the complicatedness of a complex system its resilience, and the memory of the system
* PDF to PMF
* FOCUS MORE ON THE PROS AND CONS OF EACH SYNERGY MEASURE
* REDUNDANCY: from old to new, see OLBRICH 2015
* First Tononi, old
* Synergy is simply that a set of variables COOPERATING together predict a random variable
* Smax: invariant to duplicate predictors, between 0 and the MI(X:Y), overestimates as it assumes unique information is synergy when multiple predictors have unique information on the target
* WMS synergy: not between 0 and MI(X:Y), negative is redundancy, plus is synergy, definition is the total MI minus the individual MIs, if there is tons of redundancy this becomes negative, underestimates synergy becoming loser with higher n because it n-subtracts n-redundancy (THIS IS WHY IT GETS WORSE)
* Synergistic mutual information: idea is whole minus union, instead of minus sum (the sum causes the duplication!), this union information does not exist, global positivity, self-redundancy, invariant to reordirng, monoticity (an information poorer addition does not upset), target monoticity, right range, between upper and lower bound... but relies on a not analytically solveable formula
    * Key idea; synergy is how much the whole exceeds the UNION of its parts
    * Minimization is the problem, best we can do is a numerical algorithm
* Correlation importance: duplicate predictors decrease synergy, measures something different from synergy, does not always fall in the upper-lower range we define, is between 0 and MI(X:Y)
* Early synergy meaures start from I (upsidedown U), the redundancy measure, and go from there, synergy in turn has an I U measure
* Have a notation section at the start of my thesis, very important! (also add predictor and target to notation section, plus terminology GRN, biologically possible (afkorting bedenken?))
* The venn diagram we discuss is called a partial information diagram
    * a simple pi diagram has a surface equal to the MI between the predicting variables and the predicted variable, 
    * it is possible for a slice of information to me synergistic of 12 and redundant with 3
    * all information is unique, synergistic or redundant
    * synergy can be at different levels as {1 2 3} is a higher emergence level than {1 2}
    * use the term emergence level, let this recur later
* Interest in computational biology is great; many phenotypical traits are not coded by one gene, but by a set of COOPERATING genes
* Griffith gives a nice example for redundant information (replica), synergistic information (XOR) and unique information (individual copy)
* Synergy disappears if it is made redundant, but does not disappear if one of the variables part of the synergy is duplicate

# Finished sources

* griffith2011quantifying

Next williams2010nonnegative.

# Killed

* The general consensus is that complexity is strongly dependent on scale 
* Ideally, this quantification allows for the distinction of regular systems, chaotic systems, and systems that show complex behavior.