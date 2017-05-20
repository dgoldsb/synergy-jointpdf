'''
Author: Dylan Goldsborough
Email:  dgoldsb@live.nl

This script provides a class to define a simple system with ODEs.
It can be optimized to minimize nudge effect, and maximize memory.
'''

from __future__ import print_function

import copy
import itertools
import json
import logging
import math
import os
import random
import sys

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import NPEET.entropy_estimators as npeet
import numpy as np
import seaborn as sns
from scipy import io
from scipy import optimize
from scipy.integrate import ode as integrator

__author__ = 'dgoldsb'
sns.set(color_codes=True)
ROOT = os.path.join(os.path.dirname(os.path.realpath(__file__)), '../..')
MYLOGGER = logging.getLogger('mylogger')

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
        """
        Set up a logger and initialize the variables
        """
        # Set up a logger
        handler1 = logging.FileHandler(filename=os.path.join(ROOT, 'log/resilientpdf.log'),\
                                    mode='w')
        handler1.setLevel(logging.DEBUG)
        handler1.setFormatter(logging.Formatter('%(asctime)s '+\
                            '[%(levelname)s] %(message)s'))
        MYLOGGER.addHandler(handler1)
        handler2 = logging.StreamHandler(sys.stdout)
        handler2.setLevel(logging.INFO)
        handler2.setFormatter(logging.Formatter('%(asctime)s '+\
                            '[%(levelname)s] %(message)s'))
        MYLOGGER.addHandler(handler2)
        MYLOGGER.setLevel(logging.DEBUG)

        # Initialize the variables
        self.cycle = 0
        self.iteration = 0
        self.current_sample = None
        self.state_nudged = None

        # Initialize the config dict
        self.config = {}
        self.config["filename"] = None
        self.config["size"] = 0
        self.config["components"] = []
        self.config["ode_params"] = []
        self.config["self_loop"] = self_loop

        # Initialize the setting dict
        self.setting = {}
        self.setting["sample_size"] = sample_size
        self.setting["num_nudged"] = num_nudged
        self.setting["error_mean"] = error_mean
        self.setting["error_sd"] = error_sd
        self.setting["delta_t"] = delta_t
        self.setting["visualmode"] = visualmode # 0 nothing, 1 before/after, 2 every training iter

        # Initialize the result dict
        self.training = {}
        self.training['cycle'] = []
        self.training['iteration'] = []
        self.training['memory'] = []
        self.training['resilience'] = []

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
        num_parameters = self.config["size"]**2
        if self.config["self_loop"]:
            num_parameters = num_parameters + self.config["size"]
        self.config["ode_params"] = []
        for _ in range(0, num_parameters):
            self.config["ode_params"].append(np.random.normal(0, 1))
        MYLOGGER.info("Making an initial guess for the ODE parameter vector: "\
                     +str(self.config["ode_params"]))
        return 0

    def add_component(self, mean, stdev):
        """
        Adds a single, for now normal distributed component.
        """
        self.config["components"].append(Component(mean, stdev))
        self.config["size"] = len(self.config["components"])
        MYLOGGER.info("Added component: mean "+str(mean)+" and SD "+str(stdev))
        return 0

    def sample_system(self, sample_size):
        """
        Creates a sample at time = 0 for the system.
        Returns a list of lists, each list is the sampleset of one input.
        """
        state = []
        for component in self.config["components"]:
            state.append(component.sample(sample_size))
        MYLOGGER.info("Drawn a sample of size "+str(sample_size)+" from the system")
        return state

    def nudge_system(self, state, error_mean, error_sd):
        """
        Nudges a state of the system, usually at t = 0.
        Each sample is nudged with the same draw from the same random normal distribution.
        This shifts the mean of the distribution the samples are drawn from
        """
        state_nudged = copy.deepcopy(state)
        nudges_left = self.setting["num_nudged"]
        nudge_candidates = [i for i in range(0, self.config["size"])]
        while len(nudge_candidates) > 0 and nudges_left > 0:
            # Do a nudge, addition of a noise normal distribution, effectively increases sigma
            target = random.choice(nudge_candidates)
            nudge_candidates = np.delete(nudge_candidates, list(nudge_candidates).index(target))
            nudges_left = nudges_left - 1
            for i in range(0, len(state_nudged[target])):
                state_nudged[target][i] = state_nudged[target][i] +\
                                          np.random.normal(error_mean, error_sd)

        MYLOGGER.debug("Nudged the initial sample (avg before/after):")
        MYLOGGER.debug([np.average(a) for a in state])
        MYLOGGER.debug([np.average(a) for a in state_nudged])
        return state_nudged

    def f_linear(self, t, ys, parameters):
        """
        Gives derivatives for a generic linear ODE system.
        The parameters are given in one long string, for compatibility with other libs.
        """
        derivatives = []
        for i in range(0, self.config["size"]):
            if self.config["self_loop"]:
                start = (self.config["size"] + 1) * i
                end = (self.config["size"] + 1) + (self.config["size"] + 1) * i
                a = parameters[start:end]
                dy_dt = a[0]
                for j in range(0, self.config["size"]):
                    dy_dt = dy_dt + a[j+1] * ys[j]
                derivatives.append(dy_dt)
            else:
                start = (self.config["size"]) * i
                end = (self.config["size"]) + (self.config["size"]) * i
                a = parameters[start:end]
                dy_dt = a[0]
                ys_popped = copy.deepcopy(ys)
                ys_popped = np.delete(ys_popped, i)
                for j in range(0, self.config["size"]-1):
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
        MYLOGGER.debug("Evolved the system from t = 0 to t = "+str(delta_t)+" (avg before/end)")
        MYLOGGER.debug([np.average(a) for a in state_now])
        MYLOGGER.debug([np.average(a) for a in state_future])
        MYLOGGER.debug("The parameters used are:")
        MYLOGGER.debug(parameters)
        return state_now, state_future

    def cost(self, parameters):
        """
        Generally is used to compare the unnudged state with the outcome.
        Done by calculating the MI between the unnudged state, and the final state.
        """
        if self.setting["visualmode"] == 2:
            title = 'Comparison startstates nudged vs. unnudged'
            sets = []
            for i in range(0, self.config["size"]):
                sets.append([self.state_nudged[i], self.current_sample[i]])
            self.plot_pdfs(sets, title)
        _, state_end_nudged = self.evolve_system(copy.deepcopy(self.state_nudged),
                                                 self.setting["delta_t"], parameters)
        _, state_end_unnudged = self.evolve_system(copy.deepcopy(self.current_sample),
                                                   self.setting["delta_t"], parameters)
        if self.setting["visualmode"] == 2:
            title = 'Comparison endstates nudged vs. unnudged'
            sets = []
            for i in range(0, self.config["size"]):
                sets.append([state_end_nudged[i], state_end_unnudged[i]])
            self.plot_pdfs(sets, title)

        # Use NPEET to find the MI
        # This is used to measure memory, should be maximized
        try:
            if self.setting["visualmode"] == 2:
                self.plot_mi(self.current_sample, state_end_unnudged)
            start_set = np.asarray(self.current_sample).T.tolist()
            end_set = np.asarray(state_end_unnudged).T.tolist()
            mi_mem = npeet.mi(x=start_set, y=end_set, k=10)
            if mi_mem < 0:
                mi_mem = 1e-15
                MYLOGGER.error('MI below zero encountered')
        except AssertionError as error:
            MYLOGGER.critical('NPEET failed: '+str(error))
            sys.exit(1)
        MYLOGGER.info('Mutual Information memory: '+str(mi_mem))

        # Use NPEET to find the Kullbar-Leibler divergence
        try:
            if self.setting["visualmode"] == 2:
                self.plot_kld(state_end_nudged, state_end_unnudged)
            nudged_set = np.asarray(state_end_nudged).T.tolist()
            unnudged_set = np.asarray(state_end_unnudged).T.tolist()
            kl_div_res = npeet.kldiv(x=nudged_set, xp=unnudged_set, k=10)
            if kl_div_res < 0:
                kl_div_res = 1e-15
                MYLOGGER.error('KL-div below zero encountered')
        except AssertionError as error:
            MYLOGGER.critical('NPEET failed: '+str(error))
            sys.exit(1)
        MYLOGGER.info('KL divergence resilience: '+str(kl_div_res))

        cost = - mi_mem + kl_div_res
        self.training['cycle'].append(self.cycle)
        self.training['iteration'].append(self.iteration)
        self.training['memory'].append(- mi_mem)
        self.training['resilience'].append(kl_div_res)
        self.iteration = self.iteration + 1
        MYLOGGER.info('Cost at iteration '+str(self.iteration)+' is '+str(cost))
        return cost

    def train(self, method='minimize', cycles=10):
        """
        Train the parameters using scipy
        """
        # Start with the ODE params
        # Check if there is a previous result, if not do a guess
        if len(self.config["ode_params"]) < self.config["size"]**2:
            self.generate_ode_params()

        # Preperatory work for any optimizer
        self.iteration = 1
        func = self.cost
        bounds = [[-2, 2] for _ in range(0, len(self.config["ode_params"]))]

        for i in range(0, cycles):
            MYLOGGER.info("Starting cycle "+str(i+1)+" of "+str(cycles))
            self.cycle = i + 1
            self.current_sample = self.sample_system(self.setting["sample_size"])
            self.state_nudged = self.nudge_system(copy.deepcopy(self.current_sample),
                                                  self.setting["error_mean"],
                                                  self.setting["error_sd"])
            if method == 'evolutionary':
                MYLOGGER.info("Training the ODE pars using scipy.optimize.differential_evolution")
                self.config["ode_params"] = optimize.differential_evolution(func=func,
                                                                            bounds=bounds)
            elif method == 'minimize':
                MYLOGGER.info("Training the ODE pars using scipy.optimize.minimize")
                guess = self.config["ode_params"]
                self.config["ode_params"] = optimize.minimize(fun=func, x0=guess, method="BFGS")
            elif method == 'basinhopping':
                MYLOGGER.info("Training the ODE pars using scipy.optimize.basinhopping")
                guess = self.config["ode_params"]
                self.config["ode_params"] = optimize.basinhopping(func=func, x0=guess)
                # Pickle the parameters and the training so far
            else:
                MYLOGGER.error("No optimizer selected, exiting...")
                sys.exit(1)
        return 0

    def plot_ode(self, parameters=None):
        """
        Plots the ODE system over time from the current initial values.
        Uses current parameter setup.
        Supports 3D plot for 3-dimensional system.
        Should go to 3 times the normal time
        """
        # Set up the integrator
        if parameters is None:
            parameters = self.config["ode_params"]
        if len(self.config["ode_params"]) < self.config["size"]**2:
            self.generate_ode_params()
            parameters = self.config["ode_params"]
        evolver = integrator(f=self.f_linear).set_integrator('dop853')
        evolver.set_f_params(parameters)

        # Set up the initial sample
        t_zero = 0
        y_current = [i.mean for i in self.config["components"]]
        evolver.set_initial_value(y_current)

        # Set up the final sets
        timesteps = 1000
        set_x = [i * ((float(self.setting["delta_t"]))/timesteps) for i in range(0, timesteps+1)]
        delta_t = (float(self.setting["delta_t"]))/timesteps
        set_y_t = [y_current]

        # Evolve so many times
        for i in range(0, timesteps):
            sample_future = evolver.integrate(t_zero+i*delta_t)
            set_y_t.append(sample_future)

        # Transpose
        set_y = np.asarray(set_y_t).T.tolist()

        if self.config["size"] == 3:
            MYLOGGER.info('Do with 3D phase spaceplot')
            axis = sns.plt.axes(projection='3d')
            cm_div = sns.diverging_palette(250, 15, s=75, l=40,
                                           center="dark", n=len(set_x))
            axis.scatter(set_y[0], set_y[1], set_y[2], s=20,
                         color=cm_div, depthshade=False)
            sns.plt.show()
        elif self.config["size"] == 2:
            MYLOGGER.info('Do a 2D phase space plot')
            cm_div = sns.diverging_palette(250, 15, s=75, l=40,
                                           center="dark", n=len(set_x))
            axis = sns.plt.scatter(set_y[0], set_y[1], s=20,
                                   color=cm_div)
            sns.plt.show()

        # Always do a nice plot per parameter
        MYLOGGER.info('Produce a series of plots over time')
        fig, axs = plt.subplots(self.config["size"], 1, figsize=(10, 20), facecolor='w', edgecolor='k')
        fig.suptitle('Behavior of the system over time', fontsize=16)
        fig.subplots_adjust(hspace=.5, wspace=.2)
        axs = axs.ravel()
        for i in range(0, self.config["size"]):
            axs[i].plot(set_x, set_y[i])
            axs[i].set_title(str(i+1))
            axs[i].set_xlim(0, 3*self.setting["delta_t"])
        plt.show()
        return 0

    def plot_kld(self, set_x, set_xp):
        """
        Plots the two distributions over each other.
        Also adds the Kullback-Leibler divergence.
        """
        no_per_row = 4
        if self.config["size"] < 4:
            no_per_row = self.config["size"]
        set_x_t = np.asarray(set_x).T.tolist()
        set_xp_t = np.asarray(set_xp).T.tolist()
        kld = npeet.kldiv(x=set_x_t, xp=set_xp_t, k=10)
        rows = int(math.ceil(float(self.config["size"]) / no_per_row))
        fig, axs = plt.subplots(rows, no_per_row, figsize=(10, 10), facecolor='w', edgecolor='k')
        fig.suptitle('Comparisons with KL-divergence = '+str(kld), fontsize=16)
        fig.subplots_adjust(hspace=.5, wspace=.2)
        axs = axs.ravel()
        for i in range(0, self.config["size"]):
            axs[i].hist([set_x[i], set_xp[i]], color=['r', 'b'], alpha=0.5, normed=True)
            axs[i].set_ylim(0, 1)
            axs[i].set_title(str(i+1))
        plt.show()
        return 0

    def plot_mi(self, set_x, set_y):
        """
        Plots the two distributions over each other.
        Also adds the mutual information.
        """
        no_per_row = 4
        if self.config["size"] < 4:
            no_per_row = self.config["size"]
        set_x_t = np.asarray(set_x).T.tolist()
        set_y_t = np.asarray(set_y).T.tolist()
        mutualinfo = npeet.mi(x=set_x_t, y=set_y_t, k=10)
        rows = int(math.ceil(float(self.config["size"]) / no_per_row))
        fig, axs = plt.subplots(rows, no_per_row, figsize=(10, 10), facecolor='w', edgecolor='k')
        fig.suptitle('Comparisons with mutual information = '+str(mutualinfo), fontsize=16)
        fig.subplots_adjust(hspace=.5, wspace=.2)
        axs = axs.ravel()
        for i in range(0, self.config["size"]):
            axs[i].scatter(set_x[i], set_y[i])
            axs[i].set_title(str(i+1))
            axs[i].set_xlabel("X at t=0")
            axs[i].set_ylabel("X at t=dt")
        plt.show()
        return 0

    def plot_pdfs(self, sets, title=None):
        """
        Plots a set of samples, drawn from PDFs.
        """
        no_per_row = 4
        if self.config["size"] < 4:
            no_per_row = self.config["size"]
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

    def plot_training(self):
        """
        Plots the training cost of the different repetitions.
        """
        # Select colors for the cycles
        cycles = np.max(self.training['cycle']) + 1
        cmap = plt.get_cmap('gist_rainbow')
        colors = [cmap(1.*i/cycles) for i in range(cycles)]

        # Produce a plot
        plt.title('Training visualization')
        for _cycle in range(0, cycles):
            iteration = [val for cycle, val in zip(self.training['cycle'],
                                                   self.training['iteration']) if cycle == _cycle]
            memory = [val for cycle, val in zip(self.training['cycle'],
                                                self.training['memory']) if cycle == _cycle]
            resilience = [val for cycle, val in zip(self.training['cycle'],
                                                    self.training['resilience']) if cycle == _cycle]
            cost = resilience - memory
            plt.plot(iteration, memory, ls="--", color=colors[_cycle])
            plt.plot(iteration, resilience, ls=":", color=colors[_cycle])
            plt.plot(iteration, cost, ls="-", color=colors[_cycle])
        plt.show()
        return 0

    def plot_profile(self, mode='maximum'):
        """
        Plots the MI profile of the system.
        We take increasingly larger subsets of the input sets X.
        The output Y is the system state in the future.
        """
        if mode != 'maximum' and mode != 'average':
            MYLOGGER.error('Select the mode maximum or average!')
            return 1

        # Create empty list to fill
        plot_data = [0]

        # Make the x-axis
        labels = []
        for i in range(0, self.config["size"]+1):
            labels.append(i)

        # Take a sample and do the time evaluation
        state_now = self.sample_system(self.setting["sample_size"])
        _, state_end = self.evolve_system(copy.deepcopy(state_now),
                                          self.setting["delta_t"],
                                          self.config["ode_params"])

        for i in range(1, self.config["size"]+1):
            combinations = itertools.combinations(state_now, r=i)
            entropies = []
            for combination in combinations:
                # Use NPEET to find the MI
                try:
                    start_set = np.asarray(combination).T.tolist()
                    end_set = np.asarray(state_end).T.tolist()
                    mi_mem = npeet.mi(x=start_set, y=end_set, k=10)
                    print(mi_mem)
                    if mi_mem < 0:
                        mi_mem = 1e-15
                        MYLOGGER.error('MI below zero encountered')
                except AssertionError as error:
                    MYLOGGER.critical('NPEET failed: '+str(error))
                    sys.exit(1)
                MYLOGGER.info('Mutual Information memory: '+str(mi_mem))
                entropies.append(mi_mem)
            if mode == 'average':
                plot_data.append(np.average(entropies))
            elif mode == 'maximum':
                plot_data.append(np.max(entropies))

        # Normalize with the final value
        print(float(plot_data[-1]))
        plot_data[:] = [x / float(plot_data[-1]) for x in plot_data]

        # Produce a plot
        y_ax = [x / float(self.config["size"]) for x in labels]
        plt.title('MI profile')
        plt.plot(labels, plot_data)
        plt.plot(labels, y_ax, ls="--")
        plt.show()
        return 0

    def load_config(self, filename):
        """
        Loads the configuration dictionary as a JSON.
        Tricky part is translating the list of components.
        """
        # Load the raw dict
        MYLOGGER.info("Loading the existing save "+filename)
        self.config["filename"] = filename
        raw_dict = load_dict(filename)

        # Process the raw dict
        self.config["size"] = raw_dict["size"]
        self.config["ode_params"] = raw_dict["ode_params"]
        self.config["self_loop"] = raw_dict["self_loop"]
        self.config["components"] = []
        for _component in raw_dict["components"]:
            self.config["components"].append(Component(_component[0], _component[1]))
        return 0

    def save_config(self, filename=None):
        """
        Saves the configuration dictionary as a JSON.
        Tricky part is translating the list of components.
        """
        if filename is None and self.config["filename"] is not None:
            MYLOGGER.info("Overwriting the existing save "+self.config["filename"])
        elif filename is None and self.config["filename"] is None:
            MYLOGGER.error("Cannot overwrite, the config was not loaded from a file!")
            return 1
        else:
            MYLOGGER.info("Saving as "+filename)
            self.config["filename"] = filename

        # Prepare a raw dictionary (no complex structures)
        raw_dict = {}
        raw_dict["size"] = self.config["size"]
        raw_dict["ode_params"] = self.config["ode_params"]
        raw_dict["self_loop"] = self.config["self_loop"]
        raw_dict["components"] = []
        for _component in self.config["components"]:
            raw_component = []
            raw_component.append(_component.mean)
            raw_component.append(_component.stdev)
            raw_dict["components"].append(raw_component)

        # Save the raw dictionary
        success = save_dict(filename, raw_dict)
        return success

def save_dict(filename, dictionary):
    """
    Saves to a JSON.
    """
    dirpath = os.path.join(ROOT, 'config')
    if os.path.isdir(dirpath) is False:
        os.makedirs(dirpath)
    filepath = os.path.join(dirpath, filename)
    if filepath.endswith('.json') is False:
        filepath = filepath + '.json'
    with open('data.json', 'r') as filep:
        json.dump(dictionary, filep)
    return 0


def load_dict(filename):
    """
    Loads a JSON.
    """
    dirpath = os.path.join(ROOT, 'config')
    if os.path.isdir(dirpath) is False:
        os.makedirs(dirpath)
    filepath = os.path.join(dirpath, filename)
    if filepath.endswith('.json') is False:
        filepath = filepath + '.json'
    with open('data.json', 'w') as filep:
        dictionary = json.load(filep)
    return dictionary

def main():
    """
    For basic testing purposes.
    """
    # For NPEET testing
    """
    dist1 = np.random.normal(5,1,10000)
    dist2 = [[i] for i in dist1]
    dist3 = [[i] for i in dist1]
    print(npeet.mi(dist2, dist3))
    sys.exit()
    """

    # For .mat file
    """
    config_file = os.path.join(ROOT,'config/bio_bh.mat')
    obj = io.loadmat(config_file)
    print(obj["bio"].length)
    sys.exit()
    """

    # For actual run
    system = System(num_nudged=3, error_mean=0, error_sd=0.2, visualmode=2)
    system.add_component(10, 0.5)
    system.add_component(5, 0.2)
    system.config["ode_params"] = [-1.38566526, 0.30316272, 0.20758156, -1.12657536]
    system.plot_profile()
    system.plot_ode()
    system.train(method='minimize', cycles=3)
    system.plot_ode()
    print(system.config["ode_params"])

if __name__ == '__main__':
    main()
