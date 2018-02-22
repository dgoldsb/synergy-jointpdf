# Synergy-jointpdf

## Introduction

This project investigates the presence synergy in biological networks.
Biological networks have the characteristics of both a long memory (states/oscillations are remembered over time) and resilience to pertubations (disturbances are forgotten).
We hypothesize that these two constraints lead to an abundance of synergy in biological systems.
This will be tested in gene regulation networks, using a simulation study.
In this study, we build a discrete gene regulation-like model, and optimize the updating rules to maximize memory and minimize resilience.
We will investigate if this system is indeed synergetic, and if the motifs found in the resulting model resemble gene regulation network motifs.
We will also test the hypothesis that a real network has more synergy than a random network, and that it has more memory and resilience than a random network.

## Getting started

To get started, make sure both this repository is cloned, and that the submodules are updated.
In addition, a Python init-file needs to be touched in the submodule to make it importable.
To do this, run 'setup.sh'.

To get a feel for how to run our code, we recommend running 'demo.ipynb'.
Here, we demonstrate how to use the discrete motifs.
The experiments ran for this thesis are stored under 'experiments.py'.
To generate tables and plots form the generated experimental data, run 'results.py'.
