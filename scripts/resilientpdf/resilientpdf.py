'''
Author: Dylan Goldsborough
Email:  dgoldsb@live.nl

This script provides a class to define a simple system with ODEs.
It can be optimized to minimize nudge effect, and maximize memory.
'''

from __future__ import print_function
import sys
import random
import math
import os
import copy
import logging
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import optimize
from scipy.integrate import ode as integrator
from mpl_toolkits.mplot3d import Axes3D
import NPEET.entropy_estimators as npeet

__author__ = 'dgoldsb'
sns.set(color_codes=True)
ROOT = os.path.join(os.path.dirname(os.path.realpath(__file__)), '../..')

mylogger = logging.getLogger('mylogger')
handler1 = logging.FileHandler(filename=os.path.join(ROOT, 'log/resilientpdf.log'),\
                               mode='w')
handler1.setLevel(logging.DEBUG)
handler1.setFormatter(logging.Formatter('%(asctime)s '+\
                    '[%(levelname)s] %(message)s'))
mylogger.addHandler(handler1)
handler2 = logging.StreamHandler(sys.stdout)
handler2.setLevel(logging.INFO)
handler2.setFormatter(logging.Formatter('%(asctime)s '+\
                    '[%(levelname)s] %(message)s'))
mylogger.addHandler(handler2)
mylogger.setLevel(logging.DEBUG)

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

    def __init__(self, sample_size=1000, error_mean=0, error_sd=1,
                 num_nudged=1, delta_t=1, self_loop=False, visualmode=1):
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
        self.visualmode = visualmode # 0 nothing, 1 before/after, 2 every training iter

    def plot_ode(self, parameters=None):
        """
        Plots the ODE system over time from the current initial values.
        Uses current parameter setup.
        Supports 3D plot for 3-dimensional system.
        Should go to 3 times the normal time
        """
        # Set up the integrator
        if parameters is None:
            parameters = self.ode_params
        if len(self.ode_params) < self.size**2:
            self.generate_ode_params()
            parameters = self.ode_params
        evolver = integrator(f=self.f_linear).set_integrator('dop853')
        evolver.set_f_params(parameters)

        # Set up the initial sample
        t_zero = 0
        y_current = [i.mean for i in self.components]
        evolver.set_initial_value(y_current)

        # Set up the final sets
        timesteps = 1000
        set_x = [i * ((float(self.delta_t) * 3)/timesteps) for i in range(0, timesteps+1)]
        delta_t = (float(self.delta_t) * 3)/timesteps
        set_y_t = [y_current]

        # Evolve so many times
        for i in range(0, timesteps):
            sample_future = evolver.integrate(t_zero+i*delta_t)
            set_y_t.append(sample_future)

        # Transpose
        set_y = np.asarray(set_y_t).T.tolist()

        if self.size == 3:
            mylogger.info('Do with 3D phase spaceplot')
            axis = sns.plt.axes(projection='3d')
            cm_div = sns.diverging_palette(250, 15, s=75, l=40,
                                           center="dark", n=len(set_x))
            axis.scatter(set_y[0], set_y[1], set_y[2], s=20,
                         color=cm_div, depthshade=False)
            sns.plt.show()
            return 0
        elif self.size == 2:
            mylogger.info('Do a 2D phase space plot')
            cm_div = sns.diverging_palette(250, 15, s=75, l=40,
                                           center="dark", n=len(set_x))
            axis = sns.plt.scatter(set_y[0], set_y[1], s=20,
                                   color=cm_div)
            sns.plt.show()

        # Always do a nice plot per parameter
        mylogger.info('Produce a series of plots over time')
        fig, axs = plt.subplots(self.size, 1, figsize=(10, 20), facecolor='w', edgecolor='k')
        fig.suptitle('Behavior of the system over time', fontsize=16)
        fig.subplots_adjust(hspace=.5, wspace=.2)
        axs = axs.ravel()
        for i in range(0, self.size):
            axs[i].plot(set_x, set_y[i])
            axs[i].set_title(str(i+1))
            axs[i].set_xlim(0, 3*self.delta_t)
        plt.show()
        return 0

    def plot_kld(self, set_x, set_xp):
        """
        Plots the two distributions over each other,
        Also adds the Kullback-Leibler divergence.
        """
        no_per_row = 4
        if self.size < 4:
            no_per_row = self.size
        set_x_t = np.asarray(set_x).T.tolist()
        set_xp_t = np.asarray(set_xp).T.tolist()
        kld = npeet.kldiv(x=set_x_t, xp=set_xp_t, k=10)
        rows = int(math.ceil(float(self.size) / no_per_row))
        fig, axs = plt.subplots(rows, no_per_row, figsize=(10, 10), facecolor='w', edgecolor='k')
        fig.suptitle('Comparisons with KL-divergence = '+str(kld), fontsize=16)
        fig.subplots_adjust(hspace=.5, wspace=.2)
        axs = axs.ravel()
        for i in range(0, self.size):
            axs[i].hist([set_x[i], set_xp[i]], color=['r', 'b'], alpha=0.5, normed=True)
            axs[i].set_ylim(0, 1)
            axs[i].set_title(str(i+1))
        plt.show()
        return 0

    def plot_mi(self, set_x, set_y):
        """
        Plots the two distributions over each other,
        Also adds the mutual information.
        """
        no_per_row = 4
        if self.size < 4:
            no_per_row = self.size
        set_x_t = np.asarray(set_x).T.tolist()
        set_y_t = np.asarray(set_y).T.tolist()
        mutualinfo = npeet.mi(x=set_x_t, y=set_y_t, k=10)
        rows = int(math.ceil(float(self.size) / no_per_row))
        print(rows)
        fig, axs = plt.subplots(rows, no_per_row, figsize=(10, 10), facecolor='w', edgecolor='k')
        fig.suptitle('Comparisons with mutual information = '+str(mutualinfo), fontsize=16)
        fig.subplots_adjust(hspace=.5, wspace=.2)
        axs = axs.ravel()
        for i in range(0, self.size):
            axs[i].hist([set_x[i], set_y[i]], color=['r', 'b'], alpha=0.5, normed=True)
            axs[i].set_ylim(0, 1)
            axs[i].set_title(str(i+1))
        plt.show()
        return 0

    def plot_pdfs(self, sets, title=None):
        """
        Plots a set of samples, drawn from PDFs.
        """
        no_per_row = 4
        if self.size < 4:
            no_per_row = self.size
        rows = int(math.ceil(float(len(sets)) / no_per_row))
        fig, axs = plt.subplots(rows, no_per_row, figsize=(10, 10), facecolor='w', edgecolor='k')
        if title is not None:
            fig.suptitle(title, fontsize=16)
        fig.subplots_adjust(hspace=.5, wspace=.2)
        axs = axs.ravel()
        i = 0
        for set_plot in sets:
            axs[i].hist(set_plot, alpha=0.5)
            axs[i].set_title(str(i+1))
            i = i + 1
        return 0

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
        mylogger.info("Making an initial guess for the ODE parameter vector: "\
                     +str(self.ode_params))
        return 0

    def add_component(self, mean, stdev):
        """
        Adds a single, for now normal distributed component.
        """
        self.components.append(Component(mean, stdev))
        self.size = len(self.components)
        mylogger.info("Added component: mean "+str(mean)+" and SD "+str(stdev))
        return 0

    def sample_system(self, sample_size):
        """
        Creates a sample at time = 0 for the system.
        """
        state = []
        for component in self.components:
            state.append(component.sample(sample_size))
        mylogger.info("Drawn a sample of size "+str(sample_size)+" from the system")
        return state

    def nudge_system(self, state, error_mean, error_sd):
        """
        Nudges a state of the system, usually at t = 0.
        Each sample is nudged with the same draw from the same random normal distribution.
        This shifts the mean of the distribution the samples are drawn from
        """
        state_nudged = copy.deepcopy(state)
        nudges_left = self.num_nudged
        nudge_candidates = [i for i in range(0, self.size)]
        while len(nudge_candidates) > 0 and nudges_left > 0:
            # Do a nudge, addition of a noise normal distribution, effectively increases sigma
            target = random.choice(nudge_candidates)
            nudge_candidates = np.delete(nudge_candidates, list(nudge_candidates).index(target))
            nudges_left = nudges_left - 1
            for i in range(0, len(state_nudged[target])):
                state_nudged[target][i] = state_nudged[target][i] +\
                                          np.random.normal(error_mean, error_sd)

        mylogger.debug("Nudged the initial sample (avg before/after):")
        mylogger.debug([np.average(a) for a in state])
        mylogger.debug([np.average(a) for a in state_nudged])
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
                ys_popped = copy.deepcopy(ys)
                ys_popped = np.delete(ys_popped, i)
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
        state_now = state
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
        mylogger.debug("Evolved the system from t = 0 to t = "+str(delta_t)+" (avg before/end)")
        mylogger.debug([np.average(a) for a in state_now])
        mylogger.debug([np.average(a) for a in state_future])
        mylogger.debug("The parameters used are:")
        mylogger.debug(parameters)
        return state_now, state_future

    def cost(self, parameters):
        """
        Generally is used to compare the unnudged state with the outcome.
        Done by calculating the MI between the unnudged state, and the final state.
        """
        if self.visualmode == 2:
            title = 'Comparison startstates nudged vs. unnudged'
            sets = []
            for i in range(0, self.size):
                sets.append([self.state_nudged[i], self.current_sample[i]])
            self.plot_pdfs(sets, title)
        _, state_end_nudged = self.evolve_system(copy.deepcopy(self.state_nudged), self.delta_t, parameters)
        _, state_end_unnudged = self.evolve_system(copy.deepcopy(self.current_sample), self.delta_t, parameters)
        if self.visualmode == 2:
            title = 'Comparison endstates nudged vs. unnudged'
            sets = []
            for i in range(0, self.size):
                sets.append([state_end_nudged[i], state_end_unnudged[i]])
            self.plot_pdfs(sets, title)

        # Use NPEET to find the KL-divergence
        # This is used to measure memory
        try:
            if self.visualmode == 2:
                self.plot_kld(self.current_sample, state_end_unnudged)
            start_set = np.asarray(self.current_sample).T.tolist()
            end_set = np.asarray(state_end_unnudged).T.tolist()
            kl_div_mem = npeet.kldiv(x=start_set, xp=end_set, k=10)
            if kl_div_mem < 0:
                kl_div_mem = 1e-15
                mylogger.error('KL-div below zero encountered')
        except AssertionError as error:
            mylogger.critical('NPEET failed: '+str(error))
            sys.exit(1)
        mylogger.info('KL divergence memory: '+str(kl_div_mem))

        # Use NPEET to find the Kullbar-Leibler divergence
        try:
            if self.visualmode == 2:
                self.plot_kld(state_end_nudged, state_end_unnudged)
            nudged_set = np.asarray(state_end_nudged).T.tolist()
            unnudged_set = np.asarray(state_end_unnudged).T.tolist()
            kl_div_res = npeet.kldiv(x=nudged_set, xp=unnudged_set, k=10)
            if kl_div_res < 0:
                kl_div_res = 1e-15
                mylogger.error('KL-div below zero encountered')
        except AssertionError as error:
            mylogger.critical('NPEET failed: '+str(error))
            sys.exit(1)
        mylogger.info('KL divergence resilience: '+str(kl_div_res))

        cost = kl_div_mem + kl_div_res
        self.iteration = self.iteration + 1
        mylogger.info('Cost at iteration '+str(self.iteration)+' is '+str(cost))
        return cost

    def train(self, method='minimize', cycles=10):
        """
        Train the parameters using scipy
        """
        # Start with the ODE params
        # Check if there is a previous result, if not do a guess
        if len(self.ode_params) < self.size**2:
            self.generate_ode_params()

        # Preperatory things for any optimizer
        func = self.cost
        bounds = [[-2, 2] for _ in range(0, len(self.ode_params))]

        for i in range(0, cycles):
            mylogger.info("Starting cycle "+str(i+1)+" of "+str(cycles))
            self.iteration = 1
            self.current_sample = self.sample_system(self.sample_size)
            self.state_nudged = self.nudge_system(copy.deepcopy(self.current_sample),
                                                  self.error_mean, self.error_sd)
            if method == 'evolutionary':
                mylogger.info("Training the ODE pars using scipy.optimize.differential_evolution")
                self.ode_params = optimize.differential_evolution(func=func, bounds=bounds)
            elif method == 'minimize':
                mylogger.info("Training the ODE pars using scipy.optimize.minimize")
                guess = self.ode_params
                self.ode_params = optimize.minimize(fun=func, x0=guess, bounds=bounds
                                                    , method="BFGS")
            elif method == 'basinhopping':
                mylogger.info("Training the ODE pars using scipy.optimize.basinhopping")
                guess = self.ode_params
                self.ode_params = optimize.basinhopping(func=func, x0=guess)
            else:
                mylogger.error("No optimizer selected, exiting...")
                sys.exit(1)
        return 0

def main():
    """
    For basic testing purposes.
    """
    system = System(num_nudged=1, error_mean=0, error_sd=1, visualmode=1)
    system.add_component(10, 1)
    system.add_component(8, 1)
    system.plot_ode()
    system.train(method='minimize', cycles=1)
    print(system.ode_params)

if __name__ == '__main__':
    main()
