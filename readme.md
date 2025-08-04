# LAMMPS MD Simulation with Langevin Dynamics

**M.Tech Project: Polymer Chain Dynamics and Structural Properties Analysis**

## Overview

This project investigates the diffusion coefficient scaling (D ∝ N^(-ν)) and structural properties of polymer chains using LAMMPS molecular dynamics simulations with Langevin dynamics. The study focuses on analyzing Mean Squared Displacement (MSD), Radial Distribution Function (RDF), Radius of Gyration (Rg), and Hydrodynamic Radius (Rh) for various polymer chain configurations.

## System Configuration

- **Chain Types**: Two types (A and B) with equal numbers
- **Chain Lengths**: N = 12, 16, 20, 24, 28 monomers
- **System Size**: 32³ lattice
- **Number Density**: 10% w.r.t. lattice volume
- **Potentials**: 
  - FENE (bonded interactions)
  - Lennard-Jones (non-bonded interactions)

## Study Cases

### Case 1: Homopolymers, Uncharged
- All atoms neutral
- LJ interaction parameter: ε_AB = 1.0
- Baseline case for comparison

### Case 2: Homopolymers, Modified Interaction
- All atoms neutral
- Reduced LJ interaction: ε_AB = 0.5
- Studies effect of interaction strength

### Case 3: Homopolymers, Charged
- Charges: +q on type A, -q on type B
- Charge spacing patterns: gap = 2, 3, 4 monomers
- Electrostatic effects on polymer dynamics

### Case 4: Block Copolymers
- Extends Cases 1-3 to block copolymer architectures
- Comparative analysis with homopolymers

## Key Analyses

- **MSD (Mean Squared Displacement)**: Chain dynamics and diffusion behavior
- **RDF (Radial Distribution Function)**: Local structural organization
- **Rg (Radius of Gyration)**: Overall chain size and conformation
- **Rh (Hydrodynamic Radius)**: Effective chain size in solution
- **Diffusion Coefficient Scaling**: D ∝ N^(-ν) relationship

## Repository Structure

```
├── src/
│   ├── data_generator.py      # Python script for LAMMPS data file generation
│   ├── lammps_input.in        # LAMMPS input script template
│   └── analysis.py            # Post-processing and analysis script
├── data/
│   ├── input/                 # Generated LAMMPS data files
│   └── output/                # Simulation output files
├── results/
│   ├── plots/                 # Generated analysis plots
│   └── tables/                # Numerical results and statistics
├── scripts/
│   └── run_simulation.sh      # Batch execution script
├── Makefile                   # Automated workflow management
├── requirements.txt           # Python dependencies
└── README.md                  # This file
```

## Dependencies

### Software Requirements
- **LAMMPS** (with Langevin dynamics support)
- **Python 3.8+**
- **Required Python packages**:
  ```
  numpy
  matplotlib
  scipy
  pandas
  seaborn
  ```

### Installation
```bash
# Install Python dependencies
pip install -r requirements.txt

# Ensure LAMMPS is installed and accessible
# Follow LAMMPS installation guide: https://docs.lammps.org/Install.html
```

## Usage

### Quick Start
```bash
# Run complete workflow for all cases
make all

# Run specific case
make case1  # Homopolymers, uncharged
make case2  # Homopolymers, modified interaction
make case3  # Homopolymers, charged
make case4  # Block copolymers
```

### Manual Execution

1. **Generate Data Files**:
   ```bash
   python src/data_generator.py --chain_lengths 12,16,20,24,28 --case homopolymer_uncharged
   ```

2. **Run LAMMPS Simulation**:
   ```bash
   lmp -in src/lammps_input.in -var seed 12345 -var case_type homopolymer_uncharged
   ```

3. **Analyze Results**:
   ```bash
   python src/analysis.py --input data/output/ --case homopolymer_uncharged
   ```

### Configuration Options

#### Data Generator Parameters
- `--chain_lengths`: Comma-separated list of chain lengths (default: 12,16,20,24,28)
- `--case`: Simulation case type
- `--density`: Number density percentage (default: 0.1)
- `--lattice_size`: Lattice dimensions (default: 32)

#### LAMMPS Simulation Variables
- `seed`: Random seed for reproducibility (run 4 different seeds per case)
- `case_type`: Determines interaction parameters and charge distribution
- `chain_length`: Specific chain length for current simulation

#### Analysis Options
- `--ensemble_average`: Average results over multiple seeds
- `--fit_diffusion`: Perform D ∝ N^(-ν) fitting
- `--plot_format`: Output plot format (png, pdf, svg)

## Results and Analysis

### Expected Outputs

1. **MSD vs Time Plots**: Linear regime identification for diffusion coefficient calculation
2. **RDF Plots**: Structural correlation analysis between different chain types
3. **Scaling Plots**: Rg vs N and Rh vs N relationships
4. **Diffusion Scaling**: D vs N with fitted exponent ν
5. **Comparative Analysis**: Effects of charge distribution and block architecture

### Statistical Analysis
- Ensemble averaging over 4 independent runs with different random seeds
- Error bars representing standard deviation across runs
- Statistical significance testing for scaling relationships

## Reproducibility

All simulations use controlled random seeds for reproducibility:
- Seeds: 12345, 23456, 34567, 45678
- Results averaged and analyzed with error propagation
- Version control for all input parameters and scripts

## Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/new-analysis`)
3. Commit changes (`git commit -am 'Add new analysis method'`)
4. Push to branch (`git push origin feature/new-analysis`)
5. Create Pull Request

## Validation

The simulation methodology has been validated against:
- Known polymer scaling laws
- Literature values for similar systems
- Consistency checks across different chain lengths

## References

1. Polymer Physics Literature (add specific references)
2. LAMMPS Documentation: https://docs.lammps.org/
3. Langevin Dynamics Theory (add theoretical background)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

**Author**: Sameer Wanjari
**Institution**: IIT(BHU) Varanasi
**Email**:sameern.wanjari.cd.phy21@itbhu.ac.in 
**Project Supervisor**: Dr.Awaneesh Singh

## Acknowledgments

- IIT(BHU)Varanasi High Performance Computing facility
- Research group members and collaborators

---

**Note**: This is research code developed for academic purposes. Please cite appropriately if used in publications.
