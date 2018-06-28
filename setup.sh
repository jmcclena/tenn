module purge
module load omfit/unstable
which python
export GACODE_ROOT=/fusion/projects/codes/atom/btf/gacode
#!/bin/tcsh source /fusion/projects/codes/atom/btf/gacode/shared/bin/gacode_setup.tcsh
export export BTF_HOST=irisd
python test.py


