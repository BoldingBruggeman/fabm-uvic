# fabm-uvic
FABM-uvic

This is a FABM port of the oceanic biogeochemical components of the UVic earth system model.

This work was co-funded by the European Union under grant agreement no. 101083922 (OceanICU) and UK Research and Innovation (UKRI) under the UK government’s Horizon Europe funding guarantee [grant number 10054454, 10063673, 10064020, 10059241, 10079684, 10059012, 10048179]. The views, opinions and practices used to produce this software are however those of the author(s) only and do not necessarily reflect those of the European Union or European Research Executive Agency. Neither theEuropean Union nor the granting authority can be held responsible for them.
Model code has been kept as faithful to the original code as possible.

The code has been split up into a number of submodules of which several are optional and can be in- or excluded by adjusting the runtime configuration (fabm.yaml), no code change or recompilation needed. Below is a list of included modules:

detritus.F90 | Remineralization and export of biological production

light.F90 | Light absorption and flux-at-depth

nut_chem.F90 | nitrogen, phospherous and carbon chemistry and gas exchange

phytoplankton.F90 | Phytoplankton

sediment.F90 | Solving sediment dynamics and burial

shared.F90 | Common parameters

solar.F90 | Surface solar radiation model

zooplankton.F90 | Zooplankton production and grazing

How to build
This code must be compiled together with FABM. To do this, provide the following additional arguments to cmake when you build FABM: -DFABM_INSTITUTES=uvic -DFABM_UVIC_BASE=<UVICDIR>

Here, <UVICDIR> is the directory with the FABM-UVIC code (the same directory that contains this readme file). Note that -DFABM_INSTITUTES=uvic will make FABM compile UVic as the only available biogeochemical model. If you additionally want to have access to other biogeochemical models included with FABM, you can set FABM_INSTITUTES to a semi-colon separated list, e.g., -DFABM_INSTITUTES="uvic;ersem" (to prevent the shell from interpreting the semi-colons, you typically have to enclose this list with quotes).

For instance, to use UVicC with the latest stable release of the General Ocean Turbulence Model (GOTM), do the following:

git clone --recurse-submodules -b v6.0 https://github.com/gotm-model/code.git gotm
git clone https://github.com/fabm-model/fabm.git
git clone https://github.com/BoldingBruggeman/fabm-uvic.git
mkdir build
cd build
cmake ../gotm -DFABM_BASE=../fabm -DFABM_INSTITUTES=uvic -DFABM_UVIC_BASE=../fabm-uvic
make install
This will install the GOTM executable with support for UVic at ~/local/gotm/bin/gotm.

How to run a FABM-UVic simulation
A fabm.yaml file with the UVic configuration is provided under <UVICDIR>/testcases. You can drop this file in the working directory of a FABM-compatible model such as GOTM to use it during simulation. Note that for GOTM, you will also need to ensure that fabm/use is set to true in gotm.yaml. Otherwise GOTM would run with physics only.

