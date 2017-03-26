"""
A brief exploration of the jointpdf framework
"""
from __future__ import print_function
import MIprofile
from jointpdf import jointpdf

def main():
    """
    Exploration happens here
    """

    # Create a single profile and plot it
    system = MIprofile.create_system_discrete()
    plot_data = MIprofile.create_mi_profile(system)
    MIprofile.plot_mi_profile(plot_data)

    # Create several profiles and plot them
    # First do a set of 50 systems that are unbiased
    # TODO: make a set of the classes
    # Now, do it for a set that are a function of each other
    # TODO: do this too

    # For testing: create a system, add a variable, 
    system = jointpdf.JointProbabilityMatrix(numvariables=2, numvalues=2
                                             , joint_probs='unbiased')
    print("The entropy of the system is "+str(system.entropy()))
    system.append_variables(num_added_variables=1, added_joint_probabilities=None)
    print("The entropy of the system with another variable appended is "+str(system.entropy()))
    # TODO: verify that the added variable entropy is more than the difference
    system.append_synergistic_variables(num_synergistic_variables=1
                                        , initial_guess_summed_modulo=False
                                        , verbose=True, subject_variables=None, agnostic_about=None
                                        , num_repeats=1, minimize_method=None)
    print("The entropy of the system with a synergistic variable appended is "
          +str(system.entropy()))
    system.set_labels(['a', 'b', 'c', 'd'])
    print(system.get_labels())
    print(system.mutual_information_labels('a', 'b'))

if __name__ == '__main__':
    main()
