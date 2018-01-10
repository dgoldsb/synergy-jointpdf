"""
Test what solution to use for the KL-divergence problem,
as this cannot handle states that are zero.
"""

from __future__ import print_function

from copy import deepcopy

from scipy.stats import entropy

def dirichlet_smoothing(state_0, state_1, constant):
    """
    Dirichlet smoothing with this constant.
    """
    for i in range(0, len(state_0)):
        state_0[i] = state_0[i] + constant
        state_1[i] = state_1[i] + constant

    print(entropy(state_0, state_1))

def non_causal(state_0, state_1, percentage):
    """
    Remove determinism, add uncertainty of %.
    """
    for i in range(0, len(state_0)):
        if (state_1[i] == 0) and (i%2 == 0):
            state_1[i] = percentage * state_1[i + 1]
            state_1[i+ 1] = state_1[i + 1] - percentage * state_1[i + 1]
        elif (state_1[i] == 0) and (i%2 != 0):
            state_1[i] = percentage * state_1[i - 1]
            state_1[i - 1] = state_1[i - 1] - percentage * state_1[i - 1]
    print(entropy(state_0, state_1))

def drop_zeros(state_0, state_1):
    """
    Removes zeros.
    """
    for i in range(0, len(state_0)):
        if state_1[i] == 0:
            state_0[i] = 0
    print(entropy(state_0, state_1))

def averaged_sum(state_0, state_1):
    """
    Removes zeros.
    """
    state_2 = deepcopy(state_1)
    for i in range(0, len(state_0)):
        state_2[i] = (state_0[i] + state_1[i]) / 2
    print(entropy(state_0, state_2))

def main():
    state_0 = [0.2, 0.2, 0.2, 0.4]
    state_1 = [0, 0.4, 0, 0.6]
    percentages = [0.00001, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2]

    for percentage in percentages:
        print("Percentage: " + str(percentage))
        non_causal(deepcopy(state_0), deepcopy(state_1), percentage)
        print("Constant: " + str(percentage))
        dirichlet_smoothing(deepcopy(state_0), deepcopy(state_1), percentage)

    print("Dropping zeros...")
    drop_zeros(deepcopy(state_0), deepcopy(state_1))

    print("Comparing to averaged sum...")
    averaged_sum(deepcopy(state_0), deepcopy(state_1))

    print("Comparing to real thing...")
    print(entropy(state_1, state_0))

    print("In the end I should go for Hellinger distance, MI can handle zeros...")

if __name__ == '__main__':
    main()
