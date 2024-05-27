# pet_modeling_inprogress
This repository provides the source code for the paper "Exact Parameter Identification in PET Pharmacokinetic Modeling Using the Irreversible Two Tissue Compartment Model" as cited below.
## Requirements
* The code in directory [Codes_Reconstruction](Codes_Reconstruction/) was written and tested with Matlab R2019a under Linux. No dedicated installation is needed for the code, simply download the code and get started.
* The code in directory [pet_data_sim](pet_data_sim/) was written and tested with Python $\color{red}{\text{Todo: Version}}$ under $\color{red}{\text{Todo: Distribution}}$. $\color{red}{\text{Todo: any necessary Python modules.}}$
## Running the code
Follow the instructions in the given order below to get started:

**Data Simulation:** *see files in [pet_data_sim](pet_data_sim/)*
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
> [!CAUTION]
> Ensure that simulated data is saved in [data](pet_data_sim/data/) in separate subdirectories for different patients and counting settings ($\color{red}{\text{Todo: example}}$ to guarantee compatibility with reconstruction code where data is read in (see [read_data.m](Codes_Reconstruction/read_data.m)).

**Parameter Reconstruction:** *see files in [Codes_Reconstruction](Codes_Reconstruction/)*

5. Reconstruct parameters by calling
```
global_run.m
```
> [!NOTE]
> This will reproduce the results (figures, tables) of our [paper](https://arxiv.org/abs/2305.16989) in the setup discussed there. Note further that step 5 automatically generates a .txt file which includes more detailed information on the reconstructions (compare with paper).
## Description of the code
In the following we briefly describe the purpose of each code file.

**Data Simulation:**
* <code>tac_sim.m</code>
* <code>00sim_data.py</code>
* <code>01osem.py</code>
* <code>02plot_sim.py</code>

**Parameter Reconstruction:**
* <code>read_data.m</code>
* <code>forward_model.m</code>
* <code>check_theoretical_assumptions.m</code>
* <code>global_run.m</code>
* <code>IRGNM_modified.m</code>
* <code>IRGNM_reduced.m</code>
* <code>reconstruction_noise.m</code>
* <code>reconstruction_noiseless.m</code>
* <code>comparison_plot.m</code>

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
