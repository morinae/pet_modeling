# pet_modeling_inprogress
This repository provides the source code for the paper "Exact Parameter Identification in PET Pharmacokinetic Modeling Using the Irreversible Two Tissue Compartment Model" as cited below.
## Requirements
* The code in directory [Codes_Reconstruction](Codes_Reconstruction/) was written and tested with Matlab R2019a under Linux. No dedicated installation is needed for the code, simply download the code and get started.
* The code in directory [pet_data_sim](pet_data_sim/) was written and tested with Python $\color{red}{\text{Todo: Version}}$ under $\color{red}{\text{Todo: Distribution}}$. $\color{red}{\text{Todo: any necessary Python modules.}}$
## Running the code
Get started as follows in the given order:
1. Evaluate forward model by calling
```
tac_sim.m
```
2. Generate 4D brainweb phantoms by calling
```
00sim_data.py
```
3. Simulate and reconstruct PET data by calling
```
01osem.py
```
4. Generate plots for cross check by calling
```
02plot_sim.py
```
## Description of the code

## Authors of the code
* Erion Morina [erion.morina@uni-graz.at](mailto:erion.morina@uni-graz.at)
* Georg Schramm [georg.schramm@kuleuven.be](mailto:georg.schramm@kuleuven.be)
* Martin Holler [martin.holler@uni-graz.at](mailto:martin.holler@uni-graz.at)

EM and MH are currently affiliated with the Department of Mathematics and Scientific Computing, University of Graz, Graz, Austria. GS is currently affiliated with the Department of Imaging and Pathology, KU Leuven, Leuven, Belgium.

## Publication
> [!IMPORTANT]
> If you use this code, please cite the following associated publication.
* M. Holler, E. Morina and G. Schramm. Exact parameter identification in PET pharmacokinetic modeling using the irreversible two tissue compartment model. To appear in Physics in Medicine & Biology, 2024. [arXiv](https://arxiv.org/abs/2305.16989)
  
## License
The code in this project is licensed under the GPLv3 license - see the [LICENSE](LICENSE) file for details.
