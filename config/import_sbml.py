"""
Tool to import SBML files into Python
"""

from libsbml import *

def importSBML(sbml_loc):
    """Returns the read SBML file"""
    reader = SBMLReader()
    document = reader.readSBML(sbml_loc)
    model = document.getModel()
    if model is None:
        raise "No model present..."

    return model, document

def main():
    """Main to do stuff in"""
    model, document = importSBML("yeast_7.xml")
    if (model == None):
        print("No model present." + "\n")
        return 1

if __name__ == "__main__":
    main()
