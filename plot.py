#!/usr/bin/env python3
"""
LAMMPS MD Simulation Analysis Script for M.Tech Project
Analyzes MSD, RDF, Rg, Rh, and diffusion coefficients for polymer chains
Supports homopolymers and block copolymers with various charge configurations
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import glob
import argparse
from scipy import stats
from scipy.optimize import curve_fit
import seaborn as sns

# Set plotting style
plt.style.use('seaborn-v0_8')
sns.set_palette("husl")

class LAMMPSAnalyzer:
    def __init__(self, base_dir="./"):
        self.base_dir = base_dir
        self.chain_lengths = [12, 16, 20, 24, 28]
        self.seeds = [1, 2, 3, 4]
        self.gaps = [2, 3, 4]
        
        # Results storage
        self.results = {
            'msd': {},
            'rdf': {},
            'rg': {},
            'rh': {},
            'diffusion': {}
        }
        
        # Create output directories
        os.makedirs("plots", exist_ok=True)
        os.makedirs("analysis_data", exist_ok=True)
    
    def read_msd_file(self, filename):
        """Read MSD data from LAMMPS output"""
        try:
            data = np.loadtxt(filename, skiprows=1)
            if data.ndim == 1:
                data = data.reshape(1, -1)
            return data[:, 0], data[:, 1]  # time, msd
        except Exception as e:
            print(f"Error reading MSD file {filename}: {e}")
            return None, None
    
    def read_rdf_file(self, filename):
        """Read RDF data from LAMMPS output"""
        try:
            data = np.loadtxt(filename, skiprows=4)
            if data.ndim == 1:
                data = data.reshape(1, -1)
            return data[:, 1], data[:, 2]  # distance, g(r)
        except Exception as e:
            print(f"Error reading RDF file {filename}: {e}")
            return None, None
    
    def read_gyration_file(self, filename):
        """Read radius of gyration data"""
        try:
            data = np.loadtxt(filename, skiprows=1)
            if data.ndim == 1:
                data = data.reshape(1, -1)
            return data[:, 0], data[:, 1]  # time, Rg
        except Exception as e:
            print(f"Error reading gyration file {filename}: {e}")
            return None, None
    
    def calculate_hydrodynamic_radius(self, rg_values):
        """Calculate hydrodynamic radius from radius of gyration"""
        # Rh ≈ 0.665 * Rg for polymer chains (approximate relation)
        return 0.665 * np.array(rg_values)
    
    def calculate_diffusion_coefficient(self, time, msd):
        """Calculate diffusion coefficient from MSD slope"""
        # D = slope / 6 for 3D diffusion (MSD = 6Dt for long times)
        if len(time) < 10:
            return np.nan
        
        # Use latter half of data for linear regime
        start_idx = len(time) // 2
        slope, intercept, r_value, p_value, std_err = stats.linregress(
            time[start_idx:], msd[start_idx:]
        )
        diffusion_coeff = slope / 6.0
        return diffusion_coeff, r_value**2
    
    def power_law(self, x, a, b):
        """Power law function for fitting D vs N"""
        return a * np.power(x, -b)
    
    def analyze_case(self, case_name, polymer_type="homo"):
        """Analyze a specific case across all chain lengths and seeds"""
        print(f"\nAnalyzing {case_name} ({polymer_type})...")
        
        case_results = {
            'msd': {N: [] for N in self.chain_lengths},
            'rdf': {N: [] for N in self.chain_lengths},
            'rg': {N: [] for N in self.chain_lengths},
            'rh': {N: [] for N in self.chain_lengths},
            'diffusion': {N: [] for N in self.chain_lengths}
        }
        
        for N in self.chain_lengths:
            print(f"  Processing N = {N}...")
            
            msd_data_all = []
            rdf_data_all = []
            rg_values = []
            diffusion_coeffs = []
            
            for seed in self.seeds:
                # Construct filenames based on case and parameters
                if "charged" in case_name.lower():
                    gap = int(case_name.split('_')[-1])  # Extract gap from case name
                    prefix = f"{polymer_type}_charged_N{N}_gap{gap}_seed{seed}"
                else:
                    prefix = f"{polymer_type}_{case_name}_N{N}_seed{seed}"
                
                # Read MSD data
                msd_file = f"{self.base_dir}/output/{prefix}_msd.dat"
                time, msd = self.read_msd_file(msd_file)
                if time is not None:
                    msd_data_all.append((time, msd))
                    # Calculate diffusion coefficient
                    if len(time) > 10:
                        D, r2 = self.calculate_diffusion_coefficient(time, msd)
                        diffusion_coeffs.append(D)
                
                # Read RDF data
                rdf_file = f"{self.base_dir}/output/{prefix}_rdf.dat"
                r, gr = self.read_rdf_file(rdf_file)
                if r is not None:
                    rdf_data_all.append((r, gr))
                
                # Read gyration data
                rg_file = f"{self.base_dir}/output/{prefix}_gyration.dat"
                time_rg, rg = self.read_gyration_file(rg_file)
                if rg is not None:
                    rg_values.append(np.mean(rg[-100:]))  # Average over last 100 steps
            
            # Store ensemble averages
            if msd_data_all:
                case_results['msd'][N] = msd_data_all
            if rdf_data_all:
                case_results['rdf'][N] = rdf_data_all
            if rg_values:
                case_results['rg'][N] = np.mean(rg_values)
                case_results['rh'][N] = self.calculate_hydrodynamic_radius(np.mean(rg_values))
            if diffusion_coeffs:
                case_results['diffusion'][N] = np.mean(diffusion_coeffs)
        
        self.results[case_name] = case_results
        return case_results
    
    def plot_msd_analysis(self, case_name, case_results):
        """Plot MSD vs time for all chain lengths"""
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        axes = axes.flatten()
        
        colors = plt.cm.viridis(np.linspace(0, 1, len(self.chain_lengths)))
        
        for i, N in enumerate(self.chain_lengths):
            ax = axes[i]
            msd_data = case_results['msd'].get(N, [])
            
            for j, (time, msd) in enumerate(msd_data):
                ax.loglog(time, msd, alpha=0.3, color=colors[i])
            
            # Plot ensemble average
            if msd_data:
                # Interpolate all curves to common time grid
                min_len = min(len(time) for time, msd in msd_data)
                if min_len > 1:
                    time_common = msd_data[0][0][:min_len]
                    msd_avg = np.mean([msd[:min_len] for time, msd in msd_data], axis=0)
                    ax.loglog(time_common, msd_avg, 'k-', linewidth=2, 
                             label=f'N={N} (avg)')
            
            ax.set_xlabel('Time')
            ax.set_ylabel('MSD')
            ax.set_title(f'N = {N}')
            ax.legend()
            ax.grid(True, alpha=0.3)
        
        # Remove empty subplot
        fig.delaxes(axes[-1])
        
        plt.suptitle(f'MSD Analysis - {case_name}', fontsize=16)
        plt.tight_layout()
        plt.savefig(f'plots/msd_{case_name}.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    def plot_rdf_analysis(self, case_name, case_results):
        """Plot RDF for all chain lengths"""
        fig, ax = plt.subplots(1, 1, figsize=(10, 6))
        
        colors = plt.cm.viridis(np.linspace(0, 1, len(self.chain_lengths)))
        
        for i, N in enumerate(self.chain_lengths):
            rdf_data = case_results['rdf'].get(N, [])
            
            if rdf_data:
                # Calculate ensemble average RDF
                min_len = min(len(r) for r, gr in rdf_data)
                if min_len > 1:
                    r_common = rdf_data[0][0][:min_len]
                    gr_avg = np.mean([gr[:min_len] for r, gr in rdf_data], axis=0)
                    ax.plot(r_common, gr_avg, color=colors[i], linewidth=2, 
                           label=f'N = {N}')
        
        ax.set_xlabel('Distance r')
        ax.set_ylabel('g(r)')
        ax.set_title(f'Radial Distribution Function - {case_name}')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(f'plots/rdf_{case_name}.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    def plot_scaling_analysis(self):
        """Plot Rg, Rh, and D vs N for all cases"""
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        # Plot Rg vs N
        ax = axes[0, 0]
        for case_name, case_results in self.results.items():
            N_vals = []
            rg_vals = []
            for N in self.chain_lengths:
                if N in case_results['rg'] and not np.isnan(case_results['rg'][N]):
                    N_vals.append(N)
                    rg_vals.append(case_results['rg'][N])
            
            if N_vals:
                ax.loglog(N_vals, rg_vals, 'o-', label=case_name, markersize=6)
        
        ax.set_xlabel('Chain Length N')
        ax.set_ylabel('Radius of Gyration Rg')
        ax.set_title('Rg vs N')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Plot Rh vs N
        ax = axes[0, 1]
        for case_name, case_results in self.results.items():
            N_vals = []
            rh_vals = []
            for N in self.chain_lengths:
                if N in case_results['rh'] and not np.isnan(case_results['rh'][N]):
                    N_vals.append(N)
                    rh_vals.append(case_results['rh'][N])
            
            if N_vals:
                ax.loglog(N_vals, rh_vals, 's-', label=case_name, markersize=6)
        
        ax.set_xlabel('Chain Length N')
        ax.set_ylabel('Hydrodynamic Radius Rh')
        ax.set_title('Rh vs N')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Plot D vs N
        ax = axes[1, 0]
        for case_name, case_results in self.results.items():
            N_vals = []
            D_vals = []
            for N in self.chain_lengths:
                if N in case_results['diffusion'] and not np.isnan(case_results['diffusion'][N]):
                    N_vals.append(N)
                    D_vals.append(case_results['diffusion'][N])
            
            if N_vals and len(N_vals) >= 3:
                ax.loglog(N_vals, D_vals, '^-', label=case_name, markersize=6)
                
                # Fit power law to determine ν
                try:
                    popt, pcov = curve_fit(self.power_law, N_vals, D_vals)
                    nu = popt[1]
                    N_fit = np.linspace(min(N_vals), max(N_vals), 100)
                    D_fit = self.power_law(N_fit, *popt)
                    ax.loglog(N_fit, D_fit, '--', alpha=0.7, 
                             label=f'{case_name} fit (ν={nu:.2f})')
                    
                    # Save fitting results
                    with open(f'analysis_data/diffusion_scaling_{case_name}.txt', 'w') as f:
                        f.write(f"Case: {case_name}\n")
                        f.write(f"Scaling exponent ν: {nu:.4f} ± {np.sqrt(pcov[1,1]):.4f}\n")
                        f.write(f"Prefactor A: {popt[0]:.4e}\n")
                        f.write(f"R² of fit: {stats.pearsonr(D_vals, self.power_law(N_vals, *popt))[0]**2:.4f}\n")
                        f.write("\nData points:\n")
                        for n, d in zip(N_vals, D_vals):
                            f.write(f"N={n}: D={d:.4e}\n")
                
                except Exception as e:
                    print(f"Could not fit power law for {case_name}: {e}")
        
        ax.set_xlabel('Chain Length N')
        ax.set_ylabel('Diffusion Coefficient D')
        ax.set_title('D vs N (D ∝ N^(-ν))')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Summary plot - compare ν values
        ax = axes[1, 1]
        case_names = []
        nu_values = []
        
        for case_name in self.results.keys():
            filename = f'analysis_data/diffusion_scaling_{case_name}.txt'
            if os.path.exists(filename):
                with open(filename, 'r') as f:
                    lines = f.readlines()
                    for line in lines:
                        if 'Scaling exponent ν:' in line:
                            nu = float(line.split(':')[1].split('±')[0].strip())
                            case_names.append(case_name)
                            nu_values.append(nu)
                            break
        
        if case_names:
            bars = ax.bar(range(len(case_names)), nu_values, alpha=0.7)
            ax.set_xlabel('Cases')
            ax.set_ylabel('Scaling Exponent ν')
            ax.set_title('Comparison of Scaling Exponents')
            ax.set_xticks(range(len(case_names)))
            ax.set_xticklabels(case_names, rotation=45, ha='right')
            ax.grid(True, alpha=0.3)
            
            # Add value labels on bars
            for bar, nu in zip(bars, nu_values):
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width()/2., height + 0.01,
                       f'{nu:.2f}', ha='center', va='bottom')
        
        plt.tight_layout()
        plt.savefig('plots/scaling_analysis.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    def save_summary_data(self):
        """Save summary of all results to CSV files"""
        # Rg summary
        rg_data = []
        for case_name, case_results in self.results.items():
            for N in self.chain_lengths:
                if N in case_results['rg']:
                    rg_data.append({
                        'Case': case_name,
                        'Chain_Length': N,
                        'Rg': case_results['rg'][N]
                    })
        
        if rg_data:
            pd.DataFrame(rg_data).to_csv('analysis_data/rg_summary.csv', index=False)
        
        # Diffusion coefficient summary
        diff_data = []
        for case_name, case_results in self.results.items():
            for N in self.chain_lengths:
                if N in case_results['diffusion']:
                    diff_data.append({
                        'Case': case_name,
                        'Chain_Length': N,
                        'Diffusion_Coefficient': case_results['diffusion'][N]
                    })
        
        if diff_data:
            pd.DataFrame(diff_data).to_csv('analysis_data/diffusion_summary.csv', index=False)
    
    def run_analysis(self, cases_to_analyze=None):
        """Run complete analysis for specified cases"""
        if cases_to_analyze is None:
            # Default cases
            cases_to_analyze = [
                ("case1_uncharged", "homo"),
                ("case2_modified", "homo"),
                ("case3_charged_2", "homo"),
                ("case3_charged_3", "homo"),
                ("case3_charged_4", "homo"),
                ("case4_uncharged", "block"),
                ("case4_charged_2", "block"),
                ("case4_charged_3", "block"),
                ("case4_charged_4", "block")
            ]
        
        print("Starting LAMMPS simulation analysis...")
        print(f"Base directory: {self.base_dir}")
        print(f"Chain lengths: {self.chain_lengths}")
        print(f"Seeds: {self.seeds}")
        
        # Analyze each case
        for case_name, polymer_type in cases_to_analyze:
            try:
                case_results = self.analyze_case(case_name, polymer_type)
                
                # Generate plots for this case
                self.plot_msd_analysis(case_name, case_results)
                self.plot_rdf_analysis(case_name, case_results)
                
            except Exception as e:
                print(f"Error analyzing {case_name}: {e}")
                continue
        
        # Generate scaling analysis plots
        if self.results:
            print("\nGenerating scaling analysis...")
            self.plot_scaling_analysis()
            
            print("Saving summary data...")
            self.save_summary_data()
        
        print("\nAnalysis complete!")
        print("Results saved in:")
        print("  - plots/ : All generated plots")
        print("  - analysis_data/ : Summary data and fitting results")

def main():
    parser = argparse.ArgumentParser(description='Analyze LAMMPS MD simulation results')
    parser.add_argument('--base-dir', default='./', 
                       help='Base directory containing output files')
    parser.add_argument('--case', default=None,
                       help='Specific case to analyze (default: all cases)')
    parser.add_argument('--polymer-type', default='homo', choices=['homo', 'block'],
                       help='Polymer type (homo or block)')
    
    args = parser.parse_args()
    
    analyzer = LAMMPSAnalyzer(args.base_dir)
    
    if args.case:
        # Analyze specific case
        cases_to_analyze = [(args.case, args.polymer_type)]
    else:
        # Analyze all cases
        cases_to_analyze = None
    
    analyzer.run_analysis(cases_to_analyze)

if __name__ == "__main__":
    main()