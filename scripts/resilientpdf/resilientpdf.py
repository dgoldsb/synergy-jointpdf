'''
Author: Dylan Goldsborough
Email:  dgoldsb@live.nl

This script provides a class to define a simple system with ODEs.
It can be optimized to minimize nudge effect, and maximize memory.
'''

from __future__ import print_function
import numpy as np
import random
from scipy import optimize
from scipy.integrate import ode as integrator
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

    def __init__(self, sample_size=1000, error_mean=0, error_sd=0.01, num_nudged=1, delta_t=10):
        self.size = 0
        self.components = []
        self.ode_params = []
        self.sample_size = sample_size
        self.num_nudged = num_nudged
        self.error_mean = error_mean
        self.error_sd = error_sd
        self.delta_t = delta_t
        self.iteration = 0

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
            self.ode_params.append(np.random.normal(0, 0.01))
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

    def nudge_system(self, state, error_mean, error_sd):
        """
        Nudges a state of the system, usually at t = 0.
        Each sample is nudged with a new draw from the same random normal distribution.
        """
        state_nudged = state[:]
        nudges_left = self.num_nudged
        nudge_candidates = [i for i in range(0, self.size)]
        while len(nudge_candidates) > 0 and nudges_left > 0:
            # Do a nudge
            target = random.choice(nudge_candidates)
            nudge_candidates.pop(target)
            nudges_left = nudges_left - 1
            for i in range(0, len(state_nudged[target])):
                state_nudged[target][i] = state_nudged[target][i] \
                                          + np.random.normal(error_mean, error_sd)
        return state

    def f_linear(self, t, ys, parameters):
        """
        Gives derivatives for a generic linear ODE system.
        The parameters are given in one long string, for compatibility with other libs.
        """
        derivatives = []
        for i in range(0, self.size):
            start = (self.size + 1) * i
            end = (self.size + 1) + (self.size + 1) * i
            a = parameters[start:end]
            dy_dt = a[0]
            for j in range(0, self.size):
                dy_dt = dy_dt + a[j+1] * ys[j]
            derivatives.append(dy_dt)
        return derivatives

    def evolve_system(self, state, delta_t, parameters):
        """
        Evolve the system, usually t = 0 after a nudge.
        """
        print(parameters)
        # Transpose the thing to get a list of components per sample, instead of
        # samples per component
        state_now = state[:]
        state_now_t = np.asarray(state_now).T.tolist()
        # Create an empty list to fill with future components per sample
        evolver = integrator(f=self.f_linear).set_integrator('dop853')
        evolver.set_f_params(parameters)
        state_future_t = []

        # Use RK5 to integrate the ODE system
        for sample_now in state_now_t:
            evolver.set_initial_value(sample_now)
            t_zero = 0
            sample_future = evolver.integrate(t_zero+delta_t)
            state_future_t.append(sample_future)

        state_future = np.asarray(state_future_t).T.tolist()
        return state_now, state_future

    def cost(self, parameters):
        """
        Generally is used to compare the unnudged state with the outcome.
        Done by calculating the MI between the unnudged state, and the final state.
        """
        state_initial = self.sample_system(self.sample_size)
        state_nudged = self.nudge_system(state_initial, self.error_mean, self.error_sd)
        _, state_end = self.evolve_system(state_nudged, self.delta_t, parameters)

        # Cost is the average difference between the entropy of the outcome
        # and the MI between the unnudged start and the end-state
        costs = []
        for i in range(0, self.size):
            component_initial = [[elem] for elem in state_initial[i]]
            component_end = [[elem] for elem in state_end[i]]
            entropy = entr.entropy(x=component_end, k=3)
            mutual_info = entr.mi(x=component_initial, y=component_end, k=3)
            # TODO: something like the Jaccard distance?
            # costs.append(entropy-mutual_info)
            costs.append(1/mutual_info)
            print(mutual_info)
        cost = np.average(costs)
        self.iteration = self.iteration + 1
        print('Cost at iteration '+str(self.iteration)+' is '+str(cost))

        return cost

    def train(self, method='evolutionary'):
        """
        Train the parameters using the deap library (genetic algorithm) or scipy (annealing)
        """
        # Start with the ODE params
        # Check if there is a previous result, if not do a guess
        if len(self.ode_params) != self.size**2 + self.size:
            self.generate_ode_params()

        # TODO: the cost landscape is changing too quickly, I think I should train on one sample with
        # set nudges, then continue training a few times with new data
        # if this does not work, the cost function is the problem
        # TODO: the cost function does stuff it should not do: the mutual information becomes negative and small
        if method == 'evolutionary':
            # Add noise to guess to create a population
            func = self.cost
            bounds = [[-2, 2] for _ in range(0, len(self.ode_params))]
            self.ode_params = optimize.differential_evolution(func=func, bounds=bounds)
            outcome = []
            self.ode_params = outcome
        if method == 'minimize':
            func = self.cost
            guess = self.ode_params
            bounds = [[-2, 2] for _ in range(0, len(guess))]
            self.ode_params = optimize.minimize(fun=func, x0=guess, bounds=bounds)
        return 0

def main():
    """
    For basic testing purposes.
    """
    # TODO: add debug info in every part, add verbose qualifier and logfile
    # TODO: store the cost over time with training
    system = System(num_nudged=0)
    system.add_component(4, 1)
    system.add_component(2, 1)
    system.train()
    print(system.ode_params)

if __name__ == '__main__':
    main()
