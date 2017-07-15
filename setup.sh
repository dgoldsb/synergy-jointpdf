# Make sure submodules are up to date and ready to be used
git submodule update --init --recursive
touch ./scripts/discrete_nudge/__init__.py
touch ./scripts/continuous_nudge/NPEET/__init__.py
mkdir log
touch ./log/resilientpdf.log

# Unzip matlab files
#mkdir ./bin/MatLab
#unzip ./bin/journal.pcbi.1004530.s004.ZIP -d ./bin/MatLab/
# Extract config files from MatLab files
#mkdir ./config
#mv ./bin/MatLab/S1_File/models/models_as_used/*.mat ./config
# Remove matlab files
#rm -r ./bin/MatLab
