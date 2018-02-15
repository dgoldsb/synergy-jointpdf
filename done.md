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




% section: unexpected results
% TODO: (if still exists) The low impact, high memory is a bit fucked up. % TODO: investigate this

% graveyard
% If I have more time...
% \item An actual GRN motif is optimized for memory and resilience
% If I have even more time...
% \item An actual GRN is at the Pareto boundary of the memory/resilience cost function
% \item Synergy is found at a low level in biological networks, the level of common network motifs 
% \item Synergy is found at a low level in trained random GRNs
% \item The (DJ graph) indicates a level of synergistic control that is greater than random


A solution to this problem could be a switch to a stochastic model, where the influence of the omitted part of the network is exerted in each timestep. %TODO how?
This would avoid a permanent lapse into attractors, and improve the synergy, nudge impact and memory measurements, as all these measurements involve a single timestep. %TODO zou het gebruik van een langere-termijn model niet een betere oplossing zijn dan noise toevoegen?


% check my language, network implies the whole thing
While this was surprising, this is not necessarily concerning if we consider this limitation of the model; as we ignore the omitted part of the network we cannot do multiple timesteps, making the search for cycles that span multiple timesteps inconsequential. 
This happens as the system lapses into the available attractors, leaving a predictable and mostly statistic system.
% Ik dacht hier dat het een idee dat een stochastisch model een oplossing kan bieden (elke tijdstap een invloed van de rest van het netwerk simuleren), maar misschien is toch de enige oplossing een ODE model


We found that the model results were not incredibly sensitive to this setting, as long 1-to-1 edges where the dominant type of edge. %TODO do you show this in appenbdix? reference then... try not to make unfounded statements, so always a reference to others or to your own results (in appendix or not).


At the moment, we use a fully deterministic system in which the system memory rises to 100\% over multiple timesteps.


% Draw a broad conclusion
This suggests that synergy does not function as a mechanism to improve resilience.
We ignored any synergy at a higher level, however, as well as any selective pressure in our model.
As such, we suggest that this synergistic control might be emergent in natural networks at a higher level; not at the level of GRN motifs, but at the level of combinations of GRN motifs.
This is an interesting question to tackle in future research, and due computational limitations not one that can be answered using our methodology.



% Answer (3)
We could not draw a conclusion on a difference in the impact of a single-target nudge versus a multi-target nudge with our current methodology.
However, we did find significant results with the same sign in our experiments, which suggest that if there is a difference it is a matter of magnitude, not whether there is an effect at all.



% I am leaving this out for now, the hypothesis was poorly phrased to begin with; the nudge impact cannot increase linearly with this type of nudgeand impact measurement
% we also already include the nudge width in a table, which proves enough of a point
%\begin{figure}[H]
%    \centering
%    \includegraphics[width=\textwidth]{./../../result_pandas/k=3_l=4_e=0.250000/impacts.pdf}
%    \caption{Hidden layer output}
%    \label{fig:ugh}
%\end{figure}

%\subsection{Spread and Contrasts}
% drop this? not so interesting, we already take this from the 3d plot
% just make a note that this is all available for download




# Abstract

* Multiple coordinated interventions, for me this would be multiple uncoordinated interventions

# Literature review

# Code

* Remove my CIs, no normalit

# Overig

Feed forward loop is most common, die heeft ook weinig synergy, toch?

Aangeven hoe het profiel precies werkt het meet synergy indirect, aangezien op een polyadic niveau I(X:Y) zowel redundancy als synergy bevat
Need to chck math for constistency.

Aangeven in methods dat ik deze kies omdat geen optimization




# Working area 

Add explanation on nudging, refer to DJ possibly.
%This method of nudging is reminiscent of the signal-response curves popular in biological sciences \cite{tyson2010functional}.
%These, however, do not carry any mutual information naturally.

Nice to have: BA-netwerk
Ik kan de experimenten nu supermakkelijk schrijven en doen, en dan improve ik mijn constrained searchspace tot er awesome stuff uitkomt (of niet, wat ik dan report) :D

limiet testen van hoeveel genes (timing test)

Improve random generator preferential attachment


Add stochasticity in boolean function NOT GOING TO
Too many hole in initial distribution FIXED

1. is demo logisch bezig: juiste momenten geen synergy
2. Commit werk
3. is eerste resultaat logisch
4. sturen naar Rick
5. notes hieronder verwerken
6. werk aan tweede resultaat
7. is methods logisch

## New direction

Ik kies wel voor regels van type "X moet met 2 verhoogd worden", en ik neem geen X-OR mee (problem solved) in mijn default dictionary

Ik heb eigenlijk 3 vergelijkingen: hele search space van trans tables, de constrained op biologische dingen (deze zijn based op semi-random netwerken), en echte GRN motieven

Constraints zijn

1. barabasi albers (zeker als ik ze groot kan maken)
2. Limiet op soorten relaties (+/-/++/--/AND/NAND)
3. Ratios van relaties
4. Aantal relaties
5. Bereikbaarheid van states

## In paper

Synergy komt uit amalgaam van pos en neg relaties, niet specifiek uit xor

## Multi valued logic

https://arxiv.org/pdf/1310.5697.pdf

## Bolouri_and_Davidson_BioEssays_2002

@article{bolouri2002modeling,
  title={Modeling transcriptional regulatory networks},
  author={Bolouri, Hamid and Davidson, Eric H},
  journal={BioEssays},
  volume={24},
  number={12},
  pages={1118--1129},
  year={2002},
  publisher={Wiley Online Library}
}

addition to identifying transcrip-
tion factor target genes, regulatory linkage analysis can
reveal whether regulatory inputs incident on a gene are
required in combination (logical
and) or individually (logical
or), and whether they activate (logical
imply) or repress
(logical
inversion, or
not) the expression of the regulated
gene (regulatory logic is discussed in Figs. 1–3). Thus,
linkage analysis can provide both the connectivity structure
of a GRN and a logical description of the interactions
between genes.

> Trancsription factors degrade, so case for deactivation?

not so much more useful

---

## GeardWilladsen

@article{geard2009dynamical,
  title={Dynamical approaches to modeling developmental gene regulatory networks},
  author={Geard, Nicholas and Willadsen, Kai},
  journal={Birth Defects Research Part C: Embryo Today: Reviews},
  volume={87},
  number={2},
  pages={131--142},
  year={2009},
  publisher={Wiley Online Library}
}


As such, GRNs
are an important locus of developmental control:  while epigenetic and environmental factors
play an important role, GRNs act throughout development to ensure that the correct types of
cell are produced in the correct place and at the correct time

The complexity of biological systems means that a major challenge in modelling is deciding
on an appropriate level of detail to include.  Too much detail may result in a complicated model
with  reduced  explanatory  power,  in  which  the  essential  nature  of  the  process  of  interest  is
obscured.  On the other hand, too little detail risks omitting critical processes and mechanisms,
resulting  in  a  model  whose  behaviour  is  not  an  accurate  representation  of  the  real  system.
Deciding how much detail to include in a model will be determined by the resources available
(i.e., data, methodological tools) and by the question that motivates the model.

a useful distinction is between parts lists, topology, control logic and dynamics (Schlitt and Brazma, 2007).

A variety of different dynamical systems modelling approaches have been used to simulate
the behaviour of GRNs, and all of these models share some similarities in their approach and
some common abstractions.  All dynamical systems models focus on describing and simulating:
(a) the state of the system, and (b) changes to the system state

In the general parlance of dynamical systems theory, these elements give us
states
,
tran-
sitions
—changes from one state to another—and
state spaces
.  A system’s state space can be
thought of as the collection of all of the possible states of the system, along with all possible
transitions between states.  While state spaces are an abstract construct, they provide a frame-
work for thinking about the dynamics of a system.  The most notable features of a state space
are
attractors
and
basins  of  attraction
.


In  the  RBN  model,  node  (gene)  activation  is  assumed  to  be  Boolean—a  node  is  always
simply  active  or  inactive,  with  no  intermediate  states—and  regulation  is  a  logical  function
of current node activities.

n the most general formulation,  the network structure in RBNs is randomly generated,
and the logical functions of nodes in the network are randomly selected.  This purely random
nature  is  what  makes  the  RBN  model  a  good  basis  for  comparison  with  more  biologically
plausible networks:  in order to understand what behaviours of a specific network are unusual,
a set of baseline behaviours and properties are required to make a comparison.  A significant
contribution of the RBN model was the realisation that ordered behaviour could be obtained
‘for free’; that is, without being specifically engineered into a system (Kauffman, 1993).

The model demonstrates three distinct behavioural regimes:
stable  (or  ordered),  critical  and  chaotic  (or  disordered),  where  the  ‘critical’  regime  is  best
described as the phase transition point between the stable and chaotic regimes.

The  RBN  model  provides  one  possible  baseline  for  comparison  with  real  genetic  regulatory
systems.  It also provides a theoretical framework in which different null hypotheses can be
formulated. Two notable modifications to the RBN that provide alternative null hypotheses for
regulatory network behaviour are scale-free Boolean networks and canalised Boolean networks.

Generalised logical network models are a more descriptive relative of RBNs that aim to provide
a standard method for describing regulatory interactions (Thomas, 1973) using either Boolean
or multi-valued logic (Thomas and Kaufman, 1995).  These networks are distinguished from
simpler Boolean network models primarily by multi-value logic, asynchronous continuous-time
dynamics and time-delay effects. Generalised logic networks provide a framework for modelling
systems with multiple threshold-dependent effects (rather than the single threshold afforded
by Boolean models) or for which timing effects are important.
Generalised logical networks have been used to study regulatory system behaviour in ab-
stract terms (e.g.,  Mestl et al., 1995; Edwards and Glass, 2000), and also to model specific sys-
tems.  Examples of biological systems models constructed using this paradigm include phage-
λ
(Thieffry and Thomas, 1995) and flower morphogenesis in
Arabidopsis  thaliana
(Espinosa-
Soto et al., 2004).

DE models have several advantages over logical models.  In principle, their more detailed
representation of regulatory interactions provides a more accurate representation of the physi-
cal system under investigation.  Additionally, there is a large body of dynamical systems theory
that can be used to analyse such models (Strogatz, 1994).  For example, bifurcation analysis
provides tools for determining the critical values of parameters at which the behaviour of a
system undergoes a qualitative change (see Figure 3 (b)).  As with logical models, analysing
DE models in terms of their dynamical properties can reveal how switches, oscillators and more
complex behaviours are produced from network-level features such as interacting positive and
negative feedback loops (Tyson et al., 2001, 2003; Angeli et al., 2004).
Compared to logical approaches, a disadvantage of DE models is that they contain a large
number  of  kinetic  parameters,  while  the  number  of  systems  for  which  detailed  parameter
values  are  known  is  very  small,  mostly  restricted  to  very  simple  organisms  such  as  phage-
λ
(Shea and Ackers, 1985).  One approach to dealing with unknown parameter values is to use
numerical analysis or computational learning techniques to fit the models.  This approach has
been successfully adopted in models of cell cycle control in
Xenopus
(Novak and Tyson, 1993)
and  the  segment  polarity  network  in
Drosophila
(von  Dassow  et  al.,  2000).   In  both  cases,
the models resulted in the formation of hypotheses about kinetic parameters or interactions
that  were  later  experimentally  verified  (von  Dassow  and  Odell,  2002;  Tyson  et  al.,  2002).
A  further  discovery  resulting  from  this  approach  was  that  the  dynamical  behaviour  of  the
segment polarity network was remarkably robust to variations in the parameter values (von
Dassow  et  al.,  2000).   Similarly  robust  behaviour  was  observed  for  the  signalling  network
containing the Notch-Delta pathway involved in
Drosophila
neurogenesis (Meir et al., 2002).

STOCHASTICITY

An  implicit  assumption  made  by  many  modelling  approaches  is  that  variation  in  product
concentrations is smooth and control decisions are deterministic.  In reality, the biochemical
reactions in a GRN are subject to noise from both intrinsic and extrinsic sources.  Low concen-
trations of regulatory molecules in a cell can cause reaction rates to fluctuate, and the products
of gene transcription appear not continuously but in probabilistic bursts, leading to intrinsic
noise (McAdams and Arkin, 1997; Thattai and van Oudenaarden, 2001).

Several different stochastic modelling approaches have been proposed using both logical and
DE formalisms.  In the domain of logical models, one criticism of the standard RBN model is
based on its use of deterministic synchronous updating (i.e., all nodes are always updated at
each time step) which can be considered unrealistic.

---

## Referred

@article{espinosa2004gene,
  title={A gene regulatory network model for cell-fate determination during Arabidopsis thaliana flower development that is robust and recovers experimental gene expression profiles},
  author={Espinosa-Soto, Carlos and Padilla-Longoria, Pablo and Alvarez-Buylla, Elena R},
  journal={The Plant Cell},
  volume={16},
  number={11},
  pages={2923--2939},
  year={2004},
  publisher={Am Soc Plant Biol}
}

> Logical tables, basically what I can generate based on my model
They also have their single relationships thing

> I should make a plot like fig 4

Each node, except SEP (redundant SEP1, SEP2, and SEP3 genes), stands for the activity of a single gene involved in floral organ fate determination. Most nodes could assume three levels of expression (on, 1 or 2, and off, 0) to enable different activation thresholds when experimental data was available (Thomas, 1991). The system has a finite number of possible initial conditions equal to 139,968, and each one is represented by a vector of dimension 15 in which each column corresponds to the expression state of each network node at initial conditions in the following order: FT EMF1 TFL1 LFY FUL AP1 AP3 PI AG UFO WUS AP2 SEP LUG CLF. The vector of 15 entries that keeps track of the activity level of each node describes the system at each time point. We updated the state of each node synchronously. Starting on each initial condition we iterated the network until it reached an attractor. Thus, we determined the steady gene activation states described in Table 1. All attractors were fixed point attractors (Kauffman, 1993) in which the activity level of all genes remains the same as in the previous iteration. We also kept track of the number of iterations needed for a given initial condition to attain each steady state for future model developments. The set of initial conditions that lead to each of the system's steady states is the basin of attraction of each attractor. The model can be represented by a set of difference equations in which gene interactions are modeled according to logical rules.

> We can capture these conditional relationships with multi: two weak promoters is an AND gate, XOR are an inhibitor by themselves, but prevent each other from working

> In state table, merge rows with the same output

> Maybe it is better for compatibility to do my updates based on a state table...

> Their genes default to on/unsurpressed

> My objective is not though to have the perfect prediction record

---

plotting > state table DONE
middle ground guess: tussen WMS en MI - max(single vs whole MI) DONE
unique_individual_information (misschien interessant) DONE ADDED
Fix scatterplot function DONE

# Notes for paper

* Identity map: volgens Rick measure maximaal synergy, andere manier van nadenken over synergy (dan PID achtige dingen en WMS)
* Mogelijke bias op basis van mislukte measures
* SRV werkt beter op 3-4 dan op bits


# Nice sources 

IN STUK LITERATURE BESCHRIJVEN DAT WAT LOGISCHE DESIGNKEUZES ZIJN

https://link.springer.com/chapter/10.1007/978-3-642-49321-8_15 paywall look who cite him
https://www.ncbi.nlm.nih.gov/pubmed/20868697 paywall
https://arxiv.org/pdf/1703.06746.pdf simple multivalued logic
http://www.cs.columbia.edu/4761/notes07/chapter7.3-regulation.pdf just nice for me, it seems like AND/OR can be seen as inhibiting an edge
https://hal.inria.fr/inria-00072606/document cooperative integrations, negative and als minimum
https://www.researchgate.net/profile/M_Kaufman/publication/225890380_Dynamical_behaviour_of_biological_regulatory_networks-I_Biological_role_of_feedback_loops_and_practical_use_of_the_concept_of_the_loop-characteristic_state/links/554211ec0cf224a89a3335bd/Dynamical-behaviour-of-biological-regulatory-networks-I-Biological-role-of-feedback-loops-and-practical-use-of-the-concept-of-the-loop-characteristic-state.pdf multivalued strength of interaction EN SELF LOOPS ZIE IK


MINIMUM BIOLOGISCH TE ONDERBOUWEN
OMLAAG GAAN ZONDER ZELFVERSTERKING BIOLOGISCH TE ONDERBOUWEN, DECAY
HOEVEEL LAAT EEN + VAN EEN 2 EEN 0 OMHOOG GAAN?


IN STUK LITERATURE WAAROP IK MODELBOUW BASEER
https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-7-S1-S7 zowel ER als BA is OK
http://home.himolde.no/~molka/in765/NetworkBio_Nature-Rev-Genetics-2004.pdf Barabasi claims BA
Small size, dus alles is eigenlijk wel ok
http://www.ibt.unam.mx/~erueda/Curso_IBT_2011/ng1340.pdf scalefree transcription factors produced, dus outgoing connections? PLUS k a bit over 2

https://en.wikipedia.org/wiki/Boolean_network classic boolean network has with self things has no 
https://arxiv.org/ftp/arxiv/papers/1407/1407.6117.pdf self links always work!
	negative self-links make little sense, but I guess it can do something in some cases (2+ and a - from itself makes +)
	
	
We decided on an indegree between 2 and 4, which seems typical for the smaller GRN networks \cite{lahdesmaki2003learning}.

In general, the indegree of individual genes is low as well, meaning hub formation as in Barabasi-Albers networks is atypical.
The Boolean function often consists predominantly of simple relations, such as gene A downregulates gene B \cite{lahdesmaki2003learning, schlitt2007current}.

---

BELANGRIJK

pakket phaleb, final pickles