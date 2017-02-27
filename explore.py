import MIprofile
import jointpdf

def main():
    # Maak een profiel
    system = jointpdf.JointProbabilityMatrix(numvariables=2, numvalues=2
                                             , joint_probs='unbiased')

    # Maak een set van profielen

    # Plot die profielen

if __name__ == '__main__':
    main()