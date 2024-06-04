# pet_modeling_inprogress
This repository provides the source code for the paper "Exact Parameter Identification in PET Pharmacokinetic Modeling Using the Irreversible Two Tissue Compartment Model" as cited below.
## Requirements
* The code in directory [Codes_Reconstruction](Codes_Reconstruction/) was written and tested with Matlab R2019a under Linux.
* The code in directory [pet_data_sim](pet_data_sim/) was written and tested with Python $\color{red}{\text{Todo: Version}}$ under $\color{red}{\text{Todo: Distribution}}$. See [environment.yml](pet_data_sim/environment.yml) for the dependencies which are installed by calling
```
conda env create -f environment.yml
conda activate pet_data_sim
```
  
No further dedicated installation is needed for the code, simply download the code and get started.

## Demo example
To run a quick demo example where a simple reconstruction example is computed (without need of preceeding data generation) call
```
matlab -r "run('demo_recon.m');"
```
Note that in [demo_recon.m](Codes_Reconstruction/demo_recon.m) the noise model is kept simple to enable trying out the reconstruction algorithm immediately. For a realistic physical noise model see the steps on data simulation outlined below.

## Reproducing the results of the paper
Follow the instructions in the given order below to get started:

**Data Simulation:** *see files in [pet_data_sim](pet_data_sim/)*
1. Evaluate forward model by calling
```
matlab -nodisplay -nosplash -nodesktop -r "run('eval_forward.m');exit;"
```
2. Simulate PET data
```
python pet_datagen.py
```
> [!CAUTION]
> Ensure that simulated data is saved in [data](pet_data_sim/data/) in separate subdirectories for different patients and counting settings to guarantee compatibility (both file location and names) with reconstruction code where data is read in (see [read_data.m](Codes_Reconstruction/read_data.m)).

**Parameter Reconstruction:** *see files in [Codes_Reconstruction](Codes_Reconstruction/)*

3. Reconstruct parameters by calling
```
matlab -r "run('global_run.m');"
```
> [!NOTE]
> This will reproduce the results (figures, tables) of our [paper](https://doi.org/10.1088/1361-6560/ad539e) in the setup discussed there. Note further that step 5 automatically generates a .txt file which includes more detailed information on the reconstructions (compare with paper).


> [!WARNING]
> Caution is advised regarding the change of hyperparameters (including regularization parameters, hyperparameters of IRGNM, ...) as the reconstructions are not guaranteed to succeed for different choices. See again our [paper](https://doi.org/10.1088/1361-6560/ad539e) for the concrete choices of the hyperparameters.

## Description of the code
In the following we briefly describe the purpose of each code file.

**Data Simulation:**
* <code>eval_forward.m</code>: evaluates forward model for the ground truth parameters
* <code>pet_datagen.py</code>: main script to simulate PET data
* <code>sim_data.py</code>: generates the 4D brainweb phantom
* <code>osem.py</code>: simulates and reconstructs PET data
* <code>plot_sim.py</code>: generates plots for cross check

**Parameter Reconstruction:**
* <code>demo_recon.m</code>: runs a quick demo example
* <code>global_run.m</code>: main script to recover parameter reconstruction results of paper
* <code>read_data.m</code>: reads in simulated data
* <code>forward_model.m</code>: implements forward model in PET pharmacokinetic modeling
* <code>check_theoretical_assumptions.m</code>: verifies requisites necessary for parameter identifiability posed in Theorem 15 of our [paper](https://doi.org/10.1088/1361-6560/ad539e)
* <code>IRGNM_modified.m</code>: implements Iteratively Regularized Gauss Newton Method for full setup
* <code>IRGNM_reduced.m</code>: implements Iteratively Regularized Gauss Newton Method for reduced setup
* <code>reconstruction_noise.m</code>: computes representative reconstructions and gives respective plots in noisy setup
* <code>reconstruction_noiseless.m</code>: computes representative reconstructions and gives respective plots in noiseless setup
* <code>comparison_plot.m</code>: illustrates reconstructed curves, measured data and ground truth curves for comparability

## Authors of the code
* Erion Morina [erion.morina@uni-graz.at](mailto:erion.morina@uni-graz.at)
* Georg Schramm [georg.schramm@kuleuven.be](mailto:georg.schramm@kuleuven.be)
* Martin Holler [martin.holler@uni-graz.at](mailto:martin.holler@uni-graz.at)

EM and MH are currently affiliated with the Department of Mathematics and Scientific Computing, University of Graz, Graz, Austria. GS is currently affiliated with the Department of Imaging and Pathology, KU Leuven, Leuven, Belgium.

## Publication
> [!IMPORTANT]
> If you use this code, please cite the following associated publication.
* M. Holler, E. Morina and G. Schramm. Exact parameter identification in PET pharmacokinetic modeling using the irreversible two tissue compartment model. Physics in Medicine & Biology, 2024. [https://doi.org/10.1088/1361-6560/ad539e](https://doi.org/10.1088/1361-6560/ad539e)
  
## License
The code in this project is licensed under the GPLv3 license - see the [LICENSE](LICENSE) file for details.
