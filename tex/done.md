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


% Talk about initial profile (short, this has been improved)
A more advanced complexity profile was conceived by Bar-Yam \cite{bar2004multiscale}.
The proposed complexity function
%
\begin{equation}
C(k) = \sum_{k^\prime = k}^n D(k^\prime) 
\end{equation}
%
represents the amount of information shared by at least $k$ variables.
The complexity function utelizes $D(k)$, the information that has a redundancy of $k$ or lower, defined as 
%
\begin{equation}
D(k) =  n(k) \cdot \frac{\log(M)}{N}
\end{equation}
%
where $n(k)$ is the number of possible system subsets of size $k$, and the total number of possible actions a system can take.
This profile has been applied to real-world problems of varying nature in the following years \cite{bar2013computationally}.
It is simple in the sense that it ignores any kind of probability distribution, but only considers the number of system states the system can evolve into, and the number of possible subsets.




It has been suggested by Quax that a mutual information profile such as those proposed by Tononi, where $< I(X^k;Y) >$ is plotted against the subset size $k$ \cite{QuaxPersonal, tononi1999measures}. 
When normalized, this provides us with a non-decreasing, non-negative profile with a range from 0 to 1, that allows us to detect extreme cases of synergy and redundancy.
If there is no redundancy or synergy, we expect to see a straight line.
If the profile is above the straight line, there is more redundancy than synergy in the system, and vice versa if the profile falls below this line, as synergy expresses itself as negative mutual information.
With this information, it is possible to maximum bounds to synergy and redundancy.
For instance, if the profile instantly reaches the maximum possible value, the is no synergy in the system.
If we do not wish to average out the mutual information for each subset size, but instead want to look at the extremes, we can also decide to examine the $\max [ I(X^k;Y) ]$ .
We hope that, in its application, we are able to not only attach bounds to synergy and redundancy, but also at what subset size-level it occurs.
This can provide valuable information about the way synergy and redundancy are incorporated in the structure of the system.

%TODO mismatch
Intuitively, it seems like smaller ecological systems should be the most stable.
After all, these systems consist of few species that have large, robust populations.
However, all around us we see large ecosystems with many components, where the balance seems easy to upset.
This problem was posed around the time the field of information theory was founded, in the 1950s, but at the time no answer that was well-supported by emperical data was proposed \cite{macarthur1955fluctuations}. %TODO ref
At the time, information theory was picked up as a tool to analyze the complicatedness of ecosystems; the 'evenness variable' $H$, the Shannon-entropy, was used as a measuring device for complicatedness.
The entropy was initially applied on stock measurements, to describe the proportional population sizes of different species in the ecosystem.

Shortly after, the still popular point of view was formulated that with an increased complicatedness there are more pathways to reach a consumer in the food web, and thus a higher stability.
After all, if one link between prey and predator in a food web would disappear, for instance due to a low prey population after a harsh winter, the predator population is relatively unaffected as they switch to a less-preferred yet viable prey.
When first posed by Macarthur, information theory was used again to form a definition of complicatedness \cite{macarthur1955fluctuations},
However, this time, it was used to describe the proportional sizes of biomass flows in the foodweb, not proportional stocks sizes.

%TODO waarom stabieler en ref
Nowadays, the question has solidified itself as a form of paradox.
Theoretical studies generally conclude that smaller systems should be more stable, yet in nature we observe many big an complicated predator-prey networks \cite{kondoh2003foraging}.
The preliminary answers to this question of stability and complicatedness in ecosystems since then have been conflicting, especially between emperical data and theoretical studies. \cite{pimm1984complexity}.
In the past decade, computational studies have been added to the arsenal of ecologists in their attempt to answer this paradox.
For instance, in a recent computational study the idea of stability through complicatedness due to an increased flexibility for predators has been reinvestigated, and found as a plausible explanation in the ecosystem model \cite{kondoh2003foraging}.


I got an invite from Sarah to your Bulgarian Saint Trifon's party, but I am a bit confused from where the invitation came?

