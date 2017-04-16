'''
Author: Dylan Goldsborough
Email:  dgoldsb@live.nl

This script provides a class to define a simple system with ODEs.
It can be optimized to minimize nudge effect, and maximize memory.
'''

from __future__ import print_function
import numpy as np
from deap import creator, base, tools, algorithms
from scipy.optimize import basinhopping
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
        Basically centered around zero with some noise slapped on.
        """
        # We throw the parameters in a single list
        # Total size**2 for all component interactions
        # Add size parameters for the constants
        # The pattern is as follows:
        # a, a>a, b>a, c>a, b, a>b, b>b, c>b, c etc.
        num_parameters = self.size**2 + self.size
        self.ode_params = []
        for _ in range(0, num_parameters):
            # For now going for a set error distribution
            self.ode_params.append(np.random.normal(0, 1))
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

    def evolve_system(self, state, delta_t, parameters):
        """
        Evolve the system, usually t = 0 after a nudge.
        """
        state_now = state
        state_future = state
        return state_now, state_future

    def cost(self, parameters):
        """
        Generally is used to compare the unnudged state with the outcome.
        Done by calculating the MI between the unnudged state, and the final state.
        """
        state_initial = self.sample_system(self.sample_size)
        state_nudged = self.nudge_system(state_initial, self.error_mean,
                                         self.error_sd, self.num_affected)
        state_end = self.evolve_system(state_nudged, self.delta_t, parameters)

        # Cost is the average difference between the entropy of the outcome
        # and the MI between the unnudged start and the end-state
        costs = []
        for component_initial, component_end in state_initial, state_end:
            entropy = entr.entropy(x=component_end)
            mutual_info = entr.mi(x=component_initial, y=component_end)
            costs.append(entropy-mutual_info)
        cost = np.average(costs)

        return cost

    def train(self, method):
        """
        Train the parameters using the deap library (genetic algorithm) or scipy (annealing)
        """
        # Start with the ODE params
        # Check if there is a previous result, if not do a guess
        if method == 'evolutionary':
            # Add noise to guess to create a population
            # TODO: follow the example of https://pypi.python.org/pypi/deap
            outcome = []
            self.ode_params = outcome
        if method == 'annealing':
            func = self.cost
            guess = self.ode_params
            self.ode_params = basinhopping(func, guess)
        return 0
