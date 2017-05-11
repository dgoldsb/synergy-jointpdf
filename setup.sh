# Make sure submodules are up to date and ready to be used
git submodule update --init --recursive
touch ./scripts/jointpdf/__init__.py
touch ./scripts/resilientpdf/NPEET/__init__.py
mkdir log
touch ./log/resilientpdf.log

# Unzip matlab files
mkdir ./bin/MatLab
unzip ./bin/journal.pcbi.1004530.s004.ZIP -d ./bin/MatLab/

# Build config files from MatLab files
