# PV6-HPMR #

## What is this repository for? ##

This repository contains Bruker ParaVision 6 hpMR method files and associated MATLAB tools for the design and implementation
of hyperpolarized imaging experiments using snapshot EPI readouts, IDEAL encoding and/or spectral-spatial RF excitations.

Please note, this method has currently been tested only on the CentOS 5.11 operating system.

## How do I get set up? ##

* Copy *hpMR* folder to **${PVHOME}/prog/curdir/${USER}/ParaVision/methods/src**
* Compile *hpMR* method in PV6
* Use tools in *matlab/design* folder for designing gradient and RF waveforms
    * Copy gradient waveform files to **${PVHOME}/prog/curdir/${USER}/ParaVision/exp/lists/gp**
    * Add new readout enums to *parsTypes.h* and corresponding parameters to `KAM_ReadModeRel()` in *parsRelations.c* and `EnforceHpmrDesignParams()` in *backbone.c*
    * Copy RF waveform files to **${PVHOME}/prog/curdir/${USER}/ParaVision/exp/lists/wave**
    * Add new excitation enums to *parsTypes.h* and corresponding parameters to `KAM_ExcModeRel()` in *parsRelations.c* and `EnforceHpmrDesignParams()` in *backbone.c*
    * Recompile *hpMR* method in PV6 to enable new readout and excitation modes
* Use tools in *matlab/recon* folder for analyzing and reconstructing acquired data

## Who do I talk to? ##

Keith Michel\
kamichel at mdanderson.org\
The University of Texas MD Anderson Cancer Center\
Department of Imaging Physics\
[Magnetic Resonance Systems Laboratory](https://www.mdanderson.org/research/departments-labs-institutes/labs/magnetic-resonance-systems.html)
