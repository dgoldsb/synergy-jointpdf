'''
Author: Dylan Goldsborough
Email:  dgoldsb@live.nl

This script tests the hypotheses.
'''

import resilientpdf

def test_synergy_vs_random(config):
    """
    Tests the synergy of a config versus random systems.
    Saves the name of the config (as saved on disk) with the results in a JSON.
    """
    system = resilientpdf.System()
    system.load_config(config)

    return 0

def main():
    """
    Main.
    """
    return 0

if __name__ == '__main__':
    main()
