'''
Author: Dylan Goldsborough
Email:  dgoldsb@live.nl

This script provides a class to define a simple system with ODEs.
It can be optimized to minimize nudge effect, and maximize memory.
'''

from __future__ import print_function
import numpy as np
import NPEET.entropy_estimators as entr

__author__ = 'dgoldsb'

class Component(object):
    """
    A component for an ODE system.
    For now there is only normal distribution support.
    """
    def __init__(self, mean, stdev):
        self.mean = mean
        self.stdev = stdev

    def sample(self, sample_size):
        """
        Generates a sample from a component.
        """
        return np.random.normal(loc=self.mean, scale=self.stdev, size=sample_size)

class System(object):
    """
    An ODE system starts empty, variables are added with a defined distribution.
    An ODE system can be optimized using various optimization algorithms.
    It has a build in KNN-based cost function.
    """

    def __init__(self, sample_size=100, error_mean=0, error_sd=1, num_affected=1, delta_t=1):
        self.size = 0
        self.components = []
        self.ode_params = []
        self.sample_size = sample_size
        self.num_affected = num_affected
        self.error_mean = error_mean
        self.error_sd = error_sd
        self.delta_t = delta_t

    def generate_ode_params(self):
        """
        Generates the ODE system.
        Is called as a part of the training.
        Should be done initially, it creates enough parameters and makes some guesses.
        """
        self.ode_params = []
        return 0

    def add_component(self, mean, stdev):
        """
        Adds a single, for now normal distributed component.
        """
        self.components.append(Component(mean, stdev))
        self.size = len(self.components)
        return 0

    def sample_system(self, sample_size):
        """
        Creates a sample at time = 0 for the system.
        """
        state = []
        for component in self.components:
            state.append(component.sample(sample_size))
        return state

    def nudge_system(self, state, error_mean, error_sd, num_affected):
        """
        Nudges a state of the system, usually at t = 0.
        """
        return state

    def evolve_system(self, state, delta_t):
        """
        Evolve the system, usually t = 0 after a nudge.
        """
        state_now = state
        state_future = state
        return state_now, state_future

    def cost(self):
        """
        Generally is used to compare the unnudged state with the outcome.
        Done by calculating the MI between the unnudged state, and the final state.
        """
        state_initial = self.sample_system(self.sample_size)
        state_nudged = self.nudge_system(state_initial, self.error_mean,
                                         self.error_sd, self.num_affected)
        state_end = self.evolve_system(state_nudged, self.delta_t)
        for component_initial in state_initial:
            entr.ennt

        # Cost is the average difference between the entropy of the outcome
        # and the MI between the unnudged start and the end-state
        cost = 0
        return cost