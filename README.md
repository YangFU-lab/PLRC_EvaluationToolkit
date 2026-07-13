# **Comprehensive Evaluation Toolkit for Photoluminescent Radiative Cooling (PL-RC) Materials**

[![DOI]](https://zenodo.org/badge/DOI/10.5281/zenodo.xxxxxxx.svg)](https://doi.org/10.5281/zenodo.xxxxxxx)

## Overview
This repository contains the **complete suite of** custom MATLAB scripts used for the **data processing and evaluation procedures** described in our *Nature Protocols* manuscript: "Spectral characterization and calorimetric validation of photoluminescent radiative cooling materials". 

**The toolkit provides a standardized workflow for spectral reflectance measurement and thermal validation for PL-RC materials.

## Repository Structure
**- `Procedure1_Main.m`:** Script for reconstructing effective spectral reflectance using spectro-fluorescence-photometry (SFP) method and calculating corresponding ESR.

**- `Procedure1_Input.m`:** Input data for Procedure 1.

**- `Procedure2_Main.m`:** Script for combination of the effective spectral reflectance using multi-source post-monochromator spectrometry and calculating corresponding ESR. 

**- `Procedure2_Input.m`:** Input data for Procedure 2.

**- `Procedure3_Main.m`:** Fitting the effective solar reflectance (ESR) using outdoor thermal measurement.

**- `Procedure3_Input.m`:** Input data for Procedure 3.

**- `Procedure4_Main.m`:** Fitting the effective solar reflectance (ESR) using indoor thermal measurement.

**- `Procedure4_Input.m`:** Input data for Procedure 4.

## Prerequisites
- **MATLAB** (Developed and tested on MATLAB R2023a. The scripts use standard matrix operations and are expected to be fully compatible with MATLAB R2018a and later versions).
- **Optimization Toolbox** (Required for non-linear curve fitting functions **in Procedure 4**).

## Usage Instructions
1. Clone or download this repository to your local machine.
2. Open MATLAB and navigate to the downloaded folder.
3. **Execute the procedures sequentially or individually based on your analysis needs:** 
   Open `ProcedureX_Input.m` and replace the empty arrays with your actual experimental data.
   Open `ProcedureX_Main.m` to run the code and obtain the corresponding results.

## Citation
If you use this code in your research, please cite our protocol and the Zenodo archive:

> [作者 A., 作者 B., 等]. [论文标题]. *Nature Protocols* (2026). DOI: [论文DOI]
> 
> [作者 A., 作者 B., 等]. **PLRC_EvaluationToolkit** (Version v1.0.0). *Zenodo*. https://doi.org/10.5281/zenodo.xxxxxxx

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
