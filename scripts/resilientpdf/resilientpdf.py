'''
Author: Dylan Goldsborough
Email:  dgoldsb@live.nl

This script provides a class to define a simple system with ODEs.
It can be optimized to minimize nudge effect, and maximize memory.
'''

# TODO: put all below in a research logbook
# TODO: something like the Jaccard distance?
# TODO: store the cost over time with training
# TODO: the cost landscape is changing too quickly, 
# I think I should train on one sample with
# set nudges, then continue training a few times with new data
# if this does not work, the cost function is the problem
# TODO: the cost function does stuff it should not do: 
# the mutual information becomes negative and small
# TODO: de nudge moet wel altijd hetzelfde blijven?
# TODO: because of convergence, we cannot resample every time
# TODO: minimize converges nicely for trivial example (one variable, good starting guess)
# TODO: minimize converges nicely for trivial example (one variable, bad starting guess)
# TODO: should the nudge be done once? JUSTIFY
# TODO: should I nudge all samples the same? ANSWER HERE IS YES, ALL NUDGED WITH THE SAME
# TODO: implement logging with debug levels
# TODO: problem of doing the nudge only once is that the nudge will just be compensated for
# So it will learn to counter that specific nudge, not nd arbitrary nudge

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

    def __init__(self, sample_size=1000, error_mean=0, error_sd=0.1, num_nudged=1, delta_t=1, verbose=True):
        self.size = 0
        self.components = []
        self.ode_params = []
        self.sample_size = sample_size
        self.num_nudged = num_nudged
        self.error_mean = error_mean
        self.error_sd = error_sd
        self.delta_t = delta_t
        self.iteration = 0
        self.verbose = verbose
        self.current_sample = None
        self.state_nudged = None

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
        # For now going for a set error distribution
        num_parameters = self.size**2 + self.size
        self.ode_params = []
        for _ in range(0, num_parameters):
            self.ode_params.append(np.random.normal(0, 1))
        if self.verbose:
            print("Making an initial guess for the ODE parameter vector:")
            #print(self.ode_params)
        return 0

    def add_component(self, mean, stdev):
        """
        Adds a single, for now normal distributed component.
        """
        self.components.append(Component(mean, stdev))
        self.size = len(self.components)
        if self.verbose:
            print("Added component: mean "+str(mean)+" and SD "+str(stdev))
        return 0

    def sample_system(self, sample_size):
        """
        Creates a sample at time = 0 for the system.
        """
        state = []
        for component in self.components:
            state.append(component.sample(sample_size))
        if self.verbose:
            print("Drawn a sample of size "+str(sample_size)+" from the system")
        return state

    def nudge_system(self, state, error_mean, error_sd):
        """
        Nudges a state of the system, usually at t = 0.
        Each sample is nudged with the same draw from the same random normal distribution.
        This shifts the mean of the distribution the samples are drawn from
        """
        state_nudged = state[:]
        nudges_left = self.num_nudged
        nudge_candidates = [i for i in range(0, self.size)]
        while len(nudge_candidates) > 0 and nudges_left > 0:
            # Do a nudge
            target = random.choice(nudge_candidates)
            nudge_candidates.pop(target)
            nudges_left = nudges_left - 1
            state_nudged[target] = state_nudged[target] + np.random.normal(error_mean, error_sd)
        if self.verbose:
            print("Nudged the initial sample (before/after):")
            print([np.average(a) for a in state])
            print([np.average(a) for a in state_nudged])
        return state_nudged

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
        if self.verbose:
            print("Evolved the system from t = 0 to t = "+str(delta_t)+" (before/end)")
            print([np.average(a) for a in state_now])
            print([np.average(a) for a in state_future])
            print("The parameters used are:")
            print(parameters)
        return state_now, state_future

    def cost(self, parameters):
        """
        Generally is used to compare the unnudged state with the outcome.
        Done by calculating the MI between the unnudged state, and the final state.
        """
        _, state_end = self.evolve_system(self.state_nudged[:], self.delta_t, parameters)

        # Cost is the average difference between the entropy of the outcome
        # and the MI between the unnudged start and the end-state
        costs = []
        for i in range(0, self.size):
            component_initial = [[elem] for elem in self.current_sample[i]]
            component_end = [[elem] for elem in state_end[i]]
            #entropy = entr.entropy(x=component_end, k=100)
            mutual_info = entr.mi(x=component_initial, y=component_end, k=100)
            if mutual_info < 0:
                mutual_info = 1e-15
            costs.append(1/mutual_info)
            print(mutual_info)
        cost = np.average(costs)
        self.iteration = self.iteration + 1
        if self.verbose:
            print('Cost at iteration '+str(self.iteration)+' is '+str(cost))
        return cost

    def train(self, method='minimize', cycles=10):
        """
        Train the parameters using the deap library (genetic algorithm) or scipy (annealing)
        """
        # Start with the ODE params
        # Check if there is a previous result, if not do a guess
        if len(self.ode_params) != self.size**2 + self.size:
            self.generate_ode_params()

        # Preperatory things for any optimizer
        func = self.cost
        bounds = [[-2, 2] for _ in range(0, len(self.ode_params))]

        for i in range(0, cycles):
            print("Starting cycle "+str(i+1)+" of "+str(cycles))
            self.iteration = 1
            self.current_sample = self.sample_system(self.sample_size)
            self.state_nudged = self.nudge_system(self.current_sample[:], self.error_mean, self.error_sd)
            if method == 'evolutionary':
                if self.verbose:
                    print("Training the ODE parameters using scipy.optimize.differential_evolution")
                self.ode_params = optimize.differential_evolution(func=func, bounds=bounds)
            if method == 'minimize':
                if self.verbose:
                    print("Training the ODE parameters using scipy.optimize.minimize")
                guess = self.ode_params
                self.ode_params = optimize.minimize(fun=func, x0=guess, bounds=bounds
                                                    , method="Nelder-Mead")
        return 0

def main():
    """
    For basic testing purposes.
    """
    system = System(num_nudged=1)
    system.add_component(10, 1)
    #system.add_component(5, 1)
    system.train(method='minimize', cycles=1)
    print(system.ode_params)

if __name__ == '__main__':
    main()
