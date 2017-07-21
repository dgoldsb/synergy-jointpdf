# Make sure submodules are up to date and ready to be used
git submodule update --recursive
touch ./scripts/discrete_nudge/__init__.py
touch ./scripts/continuous_nudge/NPEET/__init__.py
touch ./script/discrete_nudge/jointpdf_DJ/code/__init__.py
ln -s ./script/discrete_nudge/jointpdf ./script/discrete_nudge/jointpdf_DJ/code/jointpdf

if [ ! -d "./log" ]; then
    mkdir log
    touch ./log/resilientpdf.log
fi
