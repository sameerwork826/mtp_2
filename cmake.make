# Makefile for LAMMPS MD Simulation Project
# M.Tech Project: Polymer Chain Diffusion and Structural Analysis
# Automates data generation, simulation runs, and analysis

# Configuration
PYTHON = python3
LAMMPS = lmp_serial
CHAIN_LENGTHS = 12 16 20 24 28
SEEDS = 1 2 3 4
GAPS = 2 3 4
NPROC = 4

# Directories
DATA_DIR = data
OUTPUT_DIR = output
PLOT_DIR = plots
ANALYSIS_DIR = analysis_data

# Scripts
DATA_GENERATOR = generate_data.py
LAMMPS_SCRIPT = run_simulation.in
ANALYSIS_SCRIPT = analyze_simulation.py

# Default target
.PHONY: all
all: setup homopolymers block_copolymers analysis

# Setup directories
.PHONY: setup
setup:
	@echo "Setting up directories..."
	mkdir -p $(DATA_DIR) $(OUTPUT_DIR) $(PLOT_DIR) $(ANALYSIS_DIR)
	@echo "Setup complete."

# Clean all generated files
.PHONY: clean
clean:
	@echo "Cleaning up generated files..."
	rm -rf $(DATA_DIR)/* $(OUTPUT_DIR)/* $(PLOT_DIR)/* $(ANALYSIS_DIR)/*
	@echo "Cleanup complete."

# Full clean including directories
.PHONY: clean-all
clean-all: clean
	rmdir $(DATA_DIR) $(OUTPUT_DIR) $(PLOT_DIR) $(ANALYSIS_DIR) 2>/dev/null || true

# =============================================================================
# HOMOPOLYMER SIMULATIONS
# =============================================================================

# Case 1: Uncharged homopolymers (ε_AB = 1.0)
.PHONY: case1
case1: setup
	@echo "Running Case 1: Uncharged homopolymers..."
	$(MAKE) run-case CASE=case1_uncharged POLYMER=homo EPSILON=1.0 CHARGED=0

# Case 2: Modified interaction homopolymers (ε_AB = 0.5)
.PHONY: case2
case2: setup
	@echo "Running Case 2: Modified interaction homopolymers..."
	$(MAKE) run-case CASE=case2_modified POLYMER=homo EPSILON=0.5 CHARGED=0

# Case 3: Charged homopolymers (various gaps)
.PHONY: case3
case3: setup case3-gap2 case3-gap3 case3-gap4

.PHONY: case3-gap2
case3-gap2:
	@echo "Running Case 3: Charged homopolymers (gap=2)..."
	$(MAKE) run-case CASE=case3_charged_2 POLYMER=homo EPSILON=1.0 CHARGED=1 GAP=2

.PHONY: case3-gap3
case3-gap3:
	@echo "Running Case 3: Charged homopolymers (gap=3)..."
	$(MAKE) run-case CASE=case3_charged_3 POLYMER=homo EPSILON=1.0 CHARGED=1 GAP=3

.PHONY: case3-gap4
case3-gap4:
	@echo "Running Case 3: Charged homopolymers (gap=4)..."
	$(MAKE) run-case CASE=case3_charged_4 POLYMER=homo EPSILON=1.0 CHARGED=1 GAP=4

# All homopolymer cases
.PHONY: homopolymers
homopolymers: case1 case2 case3

# =============================================================================
# BLOCK COPOLYMER SIMULATIONS
# =============================================================================

# Case 4: Block copolymers (uncharged and charged)
.PHONY: case4
case4: setup case4-uncharged case4-charged

.PHONY: case4-uncharged
case4-uncharged:
	@echo "Running Case 4: Uncharged block copolymers..."
	$(MAKE) run-case CASE=case4_uncharged POLYMER=block EPSILON=1.0 CHARGED=0

.PHONY: case4-charged
case4-charged: case4-charged-gap2 case4-charged-gap3 case4-charged-gap4

.PHONY: case4-charged-gap2
case4-charged-gap2:
	@echo "Running Case 4: Charged block copolymers (gap=2)..."
	$(MAKE) run-case CASE=case4_charged_2 POLYMER=block EPSILON=1.0 CHARGED=1 GAP=2

.PHONY: case4-charged-gap3
case4-charged-gap3:
	@echo "Running Case 4: Charged block copolymers (gap=3)..."
	$(MAKE) run-case CASE=case4_charged_3 POLYMER=block EPSILON=1.0 CHARGED=1 GAP=3

.PHONY: case4-charged-gap4
case4-charged-gap4:
	@echo "Running Case 4: Charged block copolymers (gap=4)..."
	$(MAKE) run-case CASE=case4_charged_4 POLYMER=block EPSILON=1.0 CHARGED=1 GAP=4

# All block copolymer cases
.PHONY: block_copolymers
block_copolymers: case4

# =============================================================================
# GENERIC CASE RUNNER
# =============================================================================

.PHONY: run-case
run-case:
	@echo "Processing case: $(CASE)"
	@for N in $(CHAIN_LENGTHS); do \
		echo "  Chain length N = $$N"; \
		if [ "$(CHARGED)" = "1" ]; then \
			$(PYTHON) $(DATA_GENERATOR) --chain-length $$N --polymer-type $(POLYMER) \
				--charged --gap $(GAP) --output $(DATA_DIR)/$(POLYMER)_charged_N$${N}_gap$(GAP).data; \
		else \
			$(PYTHON) $(DATA_GENERATOR) --chain-length $$N --polymer-type $(POLYMER) \
				--output $(DATA_DIR)/$(POLYMER)_$(CASE)_N$${N}.data; \
		fi; \
		for seed in $(SEEDS); do \
			echo "    Running seed $$seed..."; \
			if [ "$(CHARGED)" = "1" ]; then \
				$(LAMMPS) -var datafile $(DATA_DIR)/$(POLYMER)_charged_N$${N}_gap$(GAP).data \
					-var outputprefix $(OUTPUT_DIR)/$(POLYMER)_charged_N$${N}_gap$(GAP)_seed$${seed} \
					-var epsilon_ab $(EPSILON) -var seed $$seed -var charged 1 \
					-in $(LAMMPS_SCRIPT) > $(OUTPUT_DIR)/$(POLYMER)_charged_N$${N}_gap$(GAP)_seed$${seed}.log 2>&1; \
			else \
				$(LAMMPS) -var datafile $(DATA_DIR)/$(POLYMER)_$(CASE)_N$${N}.data \
					-var outputprefix $(OUTPUT_DIR)/$(POLYMER)_$(CASE)_N$${N}_seed$${seed} \
					-var epsilon_ab $(EPSILON) -var seed $$seed -var charged 0 \
					-in $(LAMMPS_SCRIPT) > $(OUTPUT_DIR)/$(POLYMER)_$(CASE)_N$${N}_seed$${seed}.log 2>&1; \
			fi; \
		done; \
	done

# =============================================================================
# ANALYSIS
# =============================================================================

# Run complete analysis
.PHONY: analysis
analysis:
	@echo "Running analysis on all simulation results..."
	$(PYTHON) $(ANALYSIS_SCRIPT) --base-dir ./ --output-dir $(PLOT_DIR) --analysis-dir $(ANALYSIS_DIR)
	@echo "Analysis complete. Check $(PLOT_DIR)/ and $(ANALYSIS_DIR)/ directories."

# Analyze specific case
.PHONY: analyze-case
analyze-case:
	@if [ -z "$(CASE)" ]; then \
		echo "Error: CASE variable not set. Use: make analyze-case CASE=case1_uncharged"; \
		exit 1; \
	fi
	@echo "Analyzing case: $(CASE)"
	$(PYTHON) $(ANALYSIS_SCRIPT) --base-dir ./ --case $(CASE) --polymer-type $(POLYMER) --output-dir $(PLOT_DIR) --analysis-dir $(ANALYSIS_DIR)

# Generate plots only
.PHONY: plots
plots:
	@echo "Generating plots from existing data..."
	$(PYTHON) $(ANALYSIS_SCRIPT) --base-dir ./ --output-dir $(PLOT_DIR) --analysis-dir $(ANALYSIS_DIR) --plots-only

# =============================================================================
# PARALLEL EXECUTION
# =============================================================================

# Run simulations in parallel (requires GNU parallel)
.PHONY: parallel-sims
parallel-sims: setup
	@echo "Running simulations in parallel..."
	@echo "This requires GNU parallel to be installed."
	@parallel --jobs $(NPROC) $(MAKE) run-case CASE={} POLYMER={1} EPSILON={2} CHARGED={3} GAP={4} ::: \
		case1_uncharged,homo,1.0,0,- \
		case2_modified,homo,0.5,0,- \
		case3_charged_2,homo,1.0,1,2 \
		case3_charged_3,homo,1.0,1,3 \
		case3_charged_4,homo,1.0,1,4 \
		case4_uncharged,block,1.0,0,- \
		case4_charged_2,block,1.0,1,2 \
		case4_charged_3,block,1.0,1,3 \
		case4_charged_4,block,1.0,1,4
	@echo "Parallel simulations complete."

# Run single case in parallel
.PHONY: parallel-case
parallel-case:
	@if [ -z "$(CASE)" ]; then \
		echo "Error: CASE variable not set. Use: make parallel-case CASE=case1_uncharged POLYMER=homo EPSILON=1.0 CHARGED=0"; \
		exit 1; \
	fi
	@echo "Running case $(CASE) in parallel..."
	@for N in $(CHAIN_LENGTHS); do \
		if [ "$(CHARGED)" = "1" ]; then \
			$(PYTHON) $(DATA_GENERATOR) --chain-length $$N --polymer-type $(POLYMER) \
				--charged --gap $(GAP) --output $(DATA_DIR)/$(POLYMER)_charged_N$${N}_gap$(GAP).data; \
		else \
			$(PYTHON) $(DATA_GENERATOR) --chain-length $$N --polymer-type $(POLYMER) \
				--output $(DATA_DIR)/$(POLYMER)_$(CASE)_N$${N}.data; \
		fi; \
		parallel --jobs $(NPROC) $(LAMMPS) -var datafile $(DATA_DIR)/$(POLYMER)_{1}_N$${N}{2}.data \
			-var outputprefix $(OUTPUT_DIR)/$(POLYMER)_{1}_N$${N}{2}_seed{3} \
			-var epsilon_ab $(EPSILON) -var seed {3} -var charged $(CHARGED) \
			-in $(LAMMPS_SCRIPT) ">" $(OUTPUT_DIR)/$(POLYMER)_{1}_N$${N}{2}_seed{3}.log 2>&1 ::: \
			$(if $(CHARGED),charged,gap$(GAP),$(CASE)) ::: $(SEEDS); \
	done
	@echo "Parallel case $(CASE) complete."