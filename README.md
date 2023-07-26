# GUD_closure
closure model

source code and configuration files forimplementation of the **closure models** under the **Darwin-MITgcm**.

Enclosed here are the necessary files to reproduce results in our main paper, that need to be combined with the base [Darwin-MITgcm code](https://gitlab.com/darwinproject/gud). Building and execution are dependent upon architecture used to run the code, and we refer to the [Darwin-MITgcm documentation](https://darwin3.readthedocs.io/en/latest/overview/overview.html) for guidance.  Large input files for the MITGCM offline ECCO simulation and DARWIN model inputs will need to be acquired outside of github.

Code was developed off the verification simulation [monod_eccov3_6+4](https://gitlab.com/darwinproject/gud/-/tree/gud/verification/monod_eccov3_6+4).

Two versions are present

**/ver_Gm_p2z1** contains the diamond model setup.

**/ver_Gm_p2z2** contains the parallel model setup.

Contained within each setup:

**/code** - source code containing necessary modifications to the base Darwin-MITgcm model. Files can be added to the code directory of the above-mentioned Darwin-MITgcm

**/input** - configuration files, can be added to the input directory within a given MITgcm verification.   Switching between linear and quadratic closure is accomplished by replacing the file dat.gud with data.gud.linear or data.gud.quad respectivily.

**/scripts** - scripts for compiling, running, and post-processing. These files are all setup specific.  SLURM site specific submission scripts are included for reference.
