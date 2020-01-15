source ~/miniconda/etc/profile.d/conda.sh
conda activate helio
echo
echo 'using the helio conda environment'
which python
echo '--------------------------------------------'
echo

~/miniconda/envs/helio/bin/python /home/cmoestl/pycode/heliocats/data_update.py

