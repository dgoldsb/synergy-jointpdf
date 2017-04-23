'''
Author: Dylan Goldsborough
Email:  dgoldsb@live.nl

This script provides a class to define a simple system with ODEs.
It can be optimized to minimize nudge effect, and maximize memory.
'''

from __future__ import print_function
import sys
import random
import os
import logging
import numpy as np
from scipy import optimize
from scipy.integrate import ode as integrator
import NPEET.entropy_estimators as entr
__author__ = 'dgoldsb'

ROOT = os.path.join(os.path.dirname(os.path.realpath(__file__)), '../..')
logging.basicConfig(filename=os.path.join(ROOT, 'log/resilientpdf.log'),
                    level=logging.INFO, format='%(asctime)s %(message)s')

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

    def __init__(self, sample_size=1000, error_mean=0, error_sd=0.1,
                 num_nudged=1, delta_t=1, self_loop=False):
        self.size = 0
        self.components = []
        self.ode_params = []
        self.sample_size = sample_size
        self.num_nudged = num_nudged
        self.error_mean = error_mean
        self.error_sd = error_sd
        self.delta_t = delta_t
        self.iteration = 0
        self.current_sample = None
        self.state_nudged = None
        self.self_loop = self_loop

    def generate_ode_params(self):
        """
        Generates the ODE system.
        Is called as a part of the training.
        Should be done initially, it creates enough parameters and makes some guesses.
        Basically centered around zero with some noise slapped on.
        """
        # We throw the parameters in a single list
        # Total size**2 for all component interactions with self-loop
        # Without
        # Add size parameters for the constants
        # The pattern is as follows:
        # a, a>a, b>a, c>a, b, a>b, b>b, c>b, c etc.
        # For now going for a set error distribution
        num_parameters = self.size**2
        if self.self_loop:
            num_parameters = num_parameters + self.size
        self.ode_params = []
        for _ in range(0, num_parameters):
            self.ode_params.append(np.random.normal(0, 1))
        logging.info("Making an initial guess for the ODE parameter vector: "\
                     +str(self.ode_params))
        return 0

    def add_component(self, mean, stdev):
        """
        Adds a single, for now normal distributed component.
        """
        self.components.append(Component(mean, stdev))
        self.size = len(self.components)
        logging.info("Added component: mean "+str(mean)+" and SD "+str(stdev))
        return 0

    def sample_system(self, sample_size):
        """
        Creates a sample at time = 0 for the system.
        """
        state = []
        for component in self.components:
            state.append(component.sample(sample_size))
        logging.info("Drawn a sample of size "+str(sample_size)+" from the system")
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
            # Do a nudge, addition of a noise normal distribution, effectively increases sigma
            target = random.choice(nudge_candidates)
            nudge_candidates.pop(target)
            nudges_left = nudges_left - 1
            for i in range(0, len(state_nudged[target])):
                state_nudged[target][i] = state_nudged[target][i] +\
                                          np.random.normal(error_mean, error_sd)
        logging.info("Nudged the initial sample (avg before/after):")
        logging.debug([np.average(a) for a in state])
        logging.debug([np.average(a) for a in state_nudged])
        return state_nudged

    def f_linear(self, t, ys, parameters):
        """
        Gives derivatives for a generic linear ODE system.
        The parameters are given in one long string, for compatibility with other libs.
        """
        derivatives = []
        for i in range(0, self.size):
            if self.self_loop:
                start = (self.size + 1) * i
                end = (self.size + 1) + (self.size + 1) * i
                a = parameters[start:end]
                dy_dt = a[0]
                for j in range(0, self.size):
                    dy_dt = dy_dt + a[j+1] * ys[j]
                derivatives.append(dy_dt)
            else:
                start = (self.size) * i
                end = (self.size) + (self.size) * i
                a = parameters[start:end]
                dy_dt = a[0]
                ys_popped = ys[:]
                ys_popped.pop(index=i)
                for j in range(0, self.size-1):
                    dy_dt = dy_dt + a[j+1] * ys_popped[j]
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
        logging.info("Evolved the system from t = 0 to t = "+str(delta_t)+" (avg before/end)")
        logging.debug([np.average(a) for a in state_now])
        logging.debug([np.average(a) for a in state_future])
        logging.info("The parameters used are:")
        logging.debug(parameters)
        return state_now, state_future

    def cost(self, parameters):
        """
        Generally is used to compare the unnudged state with the outcome.
        Done by calculating the MI between the unnudged state, and the final state.
        """
        _, state_end_nudged = self.evolve_system(self.state_nudged[:], self.delta_t, parameters)
        _, state_end_unnudged = self.evolve_system(self.current_sample[:], self.delta_t, parameters)

        # Cost is the average difference between the entropy of the outcome
        # and the MI between the unnudged start and the end-state
        mutual_info = entr.mi(x=state_end_nudged, y=state_end_unnudged, k=100)
        if mutual_info < 0:
            mutual_info = 1e-15
        cost = 1/mutual_info
        logging.debug(mutual_info)
        self.iteration = self.iteration + 1
        logging.info('Cost at iteration '+str(self.iteration)+' is '+str(cost))
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
            logging.info("Starting cycle "+str(i+1)+" of "+str(cycles))
            self.iteration = 1
            self.current_sample = self.sample_system(self.sample_size)
            self.state_nudged = self.nudge_system(self.current_sample[:],
                                                  self.error_mean, self.error_sd)
            if method == 'evolutionary':
                logging.info("Training the ODE pars using scipy.optimize.differential_evolution")
                self.ode_params = optimize.differential_evolution(func=func, bounds=bounds)
            elif method == 'minimize':
                logging.info("Training the ODE pars using scipy.optimize.minimize")
                guess = self.ode_params
                self.ode_params = optimize.minimize(fun=func, x0=guess, bounds=bounds
                                                    , method="Nelder-Mead")
            elif method == 'basinhopping':
                logging.info("Training the ODE pars using scipy.optimize.basinhopping")
                guess = self.ode_params
                self.ode_params = optimize.basinhopping(func=func, x0=guess)
            else:
                logging.error("No optimizer selected, exiting...")
                sys.exit(1)
        return 0

def main():
    """
    For basic testing purposes.
    """
    system = System(num_nudged=1)
    system.add_component(10, 1)
    #system.add_component(5, 1)
    system.train(method='basinhopping', cycles=1)
    print(system.ode_params)

if __name__ == '__main__':
    main()
