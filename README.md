# Introduction

This project investigates synergy in biological networks.
Biological networks have the characteristics of both a long memory (states/oscillations are remembered over time) and resilience to pertubations (disturbances are forgotten).
We hypothesize that these two constraints lead to an abundance of synergy in biological systems.
This will be tested in gene regulation networks, using a simulation study.
In this study, we will build a gene regulation-like model, and optimize the updating rules to maximize memory and minimize resilience.
We will investigate if this system is indeed synergetic, and if the motifs found in the resulting model resemble gene regulation network motifs.
We will also test the hypothesis that a real network has more synergy than a random network, and that it has more memory and resilience than a random network.

# Phases of the project

## Phase 1: bivariate gene regulation network study

### Description

The bivariate case of a gene regulation network has the advantage that the redundancy can be calculated exactly, preferably with integer-values only.
As such, we can do a small simulation study on an existing gene regulation model, without doing any kind of optimization where we alter the variables.
Because we can easily find the redundancy, we can calculate the synergy exactly.
Putting a number on the amount of synergy in the system allows us to show how important synergy is in biological systems.
This justifies the rest of the paper, and is useful information to apply for research funding in the future.
I decided with Rick to first train a simple version, with a linear ODE system, and show the outcome to the gene regulatory network researcher from Maastricht.

## Phase 2: study of MI profiles with controlled PDF systems

### Description

I plan to identify synergy in multivariate systems using the derivative of a MI plot.
This plot shows the average MI of all subsets of variable set X of size n with the output variable Y.
This MI profile is a not previously used measure, and thus a proof-of-concept should be included.

I want to work out a nice example system, preferably one used by Griffith and Ho (2015).
The important features of this profile are the 'base slope', positive spikes, and negative spikes.
The base slope is the average MI between all individual variables in X with Y.
This slope defines the line you expect when there would be zero synergy and zero redundancy.
Spikes can be used to identify that redundancy or synergy happens at that level.

## Phase 3: multivariate gene regulation network simulation

### Description

I want to build a gene regulation model that starts with a real network (based on a real gene regulation system), and then optimizes the parameters of the ODE system.
The goal is (as discussed with Gorkhan) to see if the real network is more optimized for memory and resilience than a random system, and if it has more synergy than a random system.
It should maximize the memory of the system over time, while minimizing the impact of disturbances.
An example would be a neural network: no single misfire should mess up the entire network.

Start with a vector of continuous scalar numbers, which all are drawn from a normal distribution.
Define a set of starting parameters, and define an ODE system with said parameters (for example: dA/dt = gamma * A + delta * B).
Later we can work with stochastic ODE systems, where an error term is added.
Apply a pertubation to one (or several perhaps in later experiments) of the variables.
Use a method such as Runge-Kutta 4 to find the state of the system at a later point in time.
Train the system to maximize the memory of the system, while minimizing the effect of pertubations.
There should be an element of stochasticity, as taking the MI before and after of a deterministic system makes little sense.
An initial measure for memory would be the mutual information between the initial state, and the later state.
Later, the halftime of the mutual information could be used as an improved measure.
The Kullback-Leibler divergence can be used for disturbance impact.
Training can be done using, for example, a genetic algorithm.

We can repeat the simulation many times with different starting points, to generate enough data to determine MIs and entropies.
The kNN method should be used for determining entropies and mutual informations.
For relevant drawings, see the screenshots of the whiteboard in the binaries directory.

## Phase 4: writing of the paper

### Description

Nothing out of the ordinary, I want to write a paper suited for a master's thesis project.
Optionally, I want to edit this down later for publication.

# Workflow

* Do a standup Trello meeting with myself in the morning
* Every day I want to work for at least on hour, to stay in the writing flow
* At least three times a week I want to sit down with one (any) paper, read it, and write in my Tex files based on its contents
* I make notes in my notebook, this should be enough: I don't need to summarize meetings again, parts that should be in the paper I will put in the appropriate section
* The wiki only contains brief summaries of what is in a paper, the detailed stuff should go straight as references into my Tex files
* Go into the rabbit hole in good sources, read their sources
* Make sure I end each day on a cliffhanger, in the middle of something, makes it easier to start again

# General notes/questions for the paper

* Expand the model to include complex triggers, such as the combination of smoke and seeing birds flying in one direction: do these triggers result in synergy if they are trained into the network? (see email Rick on his original research proposal)
* Discuss the difficulty of finding a good optimizer, too quickly changing cost space when I resample every time, KNN problems
* Later in scriptie voor motieven >3 andere redundancy measure, testen of er dezelfde motieven uitkomen
* Transcriptional bursting (periodic) or the plasmid case (steady state) might be interesting bivariate cases
* Periodicity is fine, the MI measure still works
* How do I train the network on the MI between the first and later state, if to determine this I need to do many reruns?
* For a 7+, have a focused literature study that shows the knowledge gaps
* Talk to Rick about his pertubation: is this not the entropy of Y? Is a measurement of the change in this a good measure for resilience?
* I can work in two directions
    * I can start with real gene networks and find synergy in these
    * I can generate systems that satisfy a high system memory and low impact of pertubations and see if these develop synergy (these are important characteristics of biological systems)
* How do I define the memory for a system if it is not static, but oscillating?

# Digital notes from meetings

## Call with Gokhan 2017-04-25

Paar uitkomsten:

- hij zoekt een netwerk van gist op (meest basic en compleet) en wat pointers voor motieven;
- hij zoekt ook een paper met genregulatiemodel die hij realistisch vindt;
- als we een aantal motieven weten dan kunnen we beter 'zoeken' naar de informatie constraints die daarop uitkomen, en de bijbehorende parameters en keuzes, zie verdere punten;
- keuze: model eerst laten equilibreren en dan pas meten? Maar leidt waarschijnlijk tot fixed points en cycles. Dus dan noise generator op iedere node zetten om 'rest van netwerk/cel' te representeren? (kijken wat werkt) 
- selfloops zijn beetje vreemd biologisch gezien, hoewel het weldegelijk voorkomt, is het voor nu beter om die weg te laten (scheelt ook zoekruimte);
- netwerken van 2 zijn beetje vreemd biologisch gezien, omdat die zichzelf heel moeilijk stabiel kan houden (expressie levels dus). Dus vanaf drie kunnen we pas iets zinnigs zeggen.