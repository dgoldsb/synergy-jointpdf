# Working area

Improve random generator preferential attachment
Add stochasticity in boolean function
Too many hole in initial distribution FIXED

1. is demo logisch bezig: juiste momenten geen synergy
2. Commit werk
3. is eerste resultaat logisch
4. sturen naar Rick
5. notes hieronder verwerken
6. werk aan tweede resultaat
7. is methods logisch

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
