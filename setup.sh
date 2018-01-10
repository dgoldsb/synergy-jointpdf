# Make sure submodules are up to date and ready to be used
git submodule update --recursive
touch ./scripts/discrete_nudge/__init__.py
touch ./scripts/continuous_nudge/NPEET/__init__.py
touch ./scripts/discrete_nudge/jointpdf_DJ/code/__init__.py
ln -s ./scripts/discrete_nudge/jointpdf ./scripts/discrete_nudge/jointpdf_DJ/code/jointpdf
ln -s ./scripts/discrete_nudge/jointpdf_DJ/code ./scripts/discrete_nudge/jointpdf_DJ_code

if [ ! -d "./log" ]; then
    mkdir log
    touch ./log/resilientpdf.log
fi
