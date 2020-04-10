source ~/miniconda/etc/profile.d/conda.sh
conda activate helio
echo
echo 'using the helio conda environment'
which python
echo 'cd to directory:'
cd /nas/helio/realcode/test/heliocats/
pwd
echo '--------------------------------------------'
echo

~/miniconda/envs/helio/bin/python /nas/helio/realcode/test/heliocats/data_update.py
~/miniconda/envs/helio/bin/python /nas/helio/realcode/test/heliocats/web_update.py

