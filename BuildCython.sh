cd ./SimulationCode/ErgodicHarvestingLib/
rm -f *.{c,cpp,so}
python ./CythonSetup.py build_ext --inplace
cd ../../