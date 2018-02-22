# Make sure submodules are up to date and ready to be used
git submodule update --recursive
touch ./scripts/__init__.py
touch ./scripts/jointpdf/__init__.py

if [ ! -d "./log" ]; then
    mkdir log
    touch ./log/resilientpdf.log
fi
