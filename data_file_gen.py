#!/usr/bin/env python3
"""
LAMMPS Data File Generator for Polymer Chain MD Simulations
M.Tech Project: Study of diffusion coefficient and structural properties

Usage:
python generate_data.py --chain_length 12 --case 1 --gap 2 --polymer_type homo --seed 1
"""

import numpy as np
import argparse
import random
from typing import List, Tuple, Dict

class PolymerDataGenerator:
    def __init__(self, lattice_size=32, density=0.1, bond_length=1.0):
        self.lattice_size = lattice_size
        self.density = density
        self.bond_length = bond_length
        self.total_sites = lattice_size**3
        self.occupied_sites = int(self.total_sites * density)
        
    def generate_lattice_positions(self) -> List[Tuple[int, int, int]]:
        """Generate all possible lattice positions"""
        positions = []
        for i in range(self.lattice_size):
            for j in range(self.lattice_size):
                for k in range(self.lattice_size):
                    positions.append((i, j, k))
        return positions
    
    def place_chains_on_lattice(self, chain_length: int, num_chains_per_type: int, 
                               polymer_type: str) -> Tuple[List, List]:
        """Place polymer chains on lattice sites"""
        available_positions = self.generate_lattice_positions()
        random.shuffle(available_positions)
        
        atoms = []
        bonds = []
        atom_id = 1
        bond_id = 1
        mol_id = 1
        
        total_chains = 2 * num_chains_per_type  # A and B chains
        
        for chain_idx in range(total_chains):
            chain_type = 1 if chain_idx < num_chains_per_type else 2  # 1=A, 2=B
            
            if polymer_type == "homo":
                # Homopolymer: all monomers same type
                monomer_types = [chain_type] * chain_length
            else:  # block copolymer
                # Block copolymer: first half A, second half B
                mid = chain_length // 2
                if chain_type == 1:  # A chain: A-block then B-block
                    monomer_types = [1] * mid + [2] * (chain_length - mid)
                else:  # B chain: B-block then A-block
                    monomer_types = [2] * mid + [1] * (chain_length - mid)
            
            # Place chain monomers
            chain_atoms = []
            for monomer_idx in range(chain_length):
                if len(available_positions) == 0:
                    raise ValueError("Not enough lattice sites for all atoms")
                
                pos = available_positions.pop()
                x, y, z = pos[0] + 0.5, pos[1] + 0.5, pos[2] + 0.5  # Center in unit cell
                
                atom_type = monomer_types[monomer_idx]
                charge = 0.0  # Will be set later for charged cases
                
                atoms.append({
                    'id': atom_id,
                    'mol': mol_id,
                    'type': atom_type,
                    'charge': charge,
                    'x': x, 'y': y, 'z': z
                })
                chain_atoms.append(atom_id)
                atom_id += 1
            
            # Create bonds within chain
            for i in range(len(chain_atoms) - 1):
                bonds.append({
                    'id': bond_id,
                    'type': 1,  # FENE bond type
                    'atom1': chain_atoms[i],
                    'atom2': chain_atoms[i + 1]
                })
                bond_id += 1
            
            mol_id += 1
        
        return atoms, bonds
    
    def apply_charges(self, atoms: List, gap: int, polymer_type: str):
        """Apply charges to atoms based on gap spacing"""
        for atom in atoms:
            mol_id = atom['mol']
            atom_type = atom['type']
            
            # Find position in chain (0-indexed)
            chain_atoms = [a for a in atoms if a['mol'] == mol_id]
            chain_atoms.sort(key=lambda x: x['id'])
            pos_in_chain = next(i for i, a in enumerate(chain_atoms) if a['id'] == atom['id'])
            
            # Apply charges with gap spacing
            if pos_in_chain % (gap + 1) == 0:  # Every (gap+1)th monomer
                if polymer_type == "homo":
                    # Homopolymer: A chains get +, B chains get -
                    if atom_type == 1:  # A type
                        atom['charge'] = 1.0
                    else:  # B type
                        atom['charge'] = -1.0
                else:  # block copolymer
                    # Block copolymer: A monomers get +, B monomers get -
                    if atom_type == 1:
                        atom['charge'] = 1.0
                    else:
                        atom['charge'] = -1.0
    
    def write_lammps_data(self, filename: str, atoms: List, bonds: List, 
                         case: int, chain_length: int):
        """Write LAMMPS data file"""
        num_atoms = len(atoms)
        num_bonds = len(bonds)
        num_atom_types = 2
        num_bond_types = 1
        
        with open(filename, 'w') as f:
            f.write(f"LAMMPS data file for Case {case}, N={chain_length}\n\n")
            
            f.write(f"{num_atoms} atoms\n")
            f.write(f"{num_bonds} bonds\n")
            f.write("0 angles\n")
            f.write("0 dihedrals\n")
            f.write("0 impropers\n\n")
            
            f.write(f"{num_atom_types} atom types\n")
            f.write(f"{num_bond_types} bond types\n\n")
            
            # Box dimensions
            f.write(f"0.0 {self.lattice_size:.1f} xlo xhi\n")
            f.write(f"0.0 {self.lattice_size:.1f} ylo yhi\n")
            f.write(f"0.0 {self.lattice_size:.1f} zlo zhi\n\n")
            
            # Masses
            f.write("Masses\n\n")
            f.write("1 1.0\n")
            f.write("2 1.0\n\n")
            
            # Atoms
            f.write("Atoms\n\n")
            for atom in atoms:
                f.write(f"{atom['id']} {atom['mol']} {atom['type']} {atom['charge']:.1f} "
                       f"{atom['x']:.6f} {atom['y']:.6f} {atom['z']:.6f}\n")
            
            # Bonds
            if bonds:
                f.write("\nBonds\n\n")
                for bond in bonds:
                    f.write(f"{bond['id']} {bond['type']} {bond['atom1']} {bond['atom2']}\n")
    
    def generate_data_file(self, chain_length: int, case: int, gap: int = 2, 
                          polymer_type: str = "homo", seed: int = 1):
        """Main function to generate LAMMPS data file"""
        random.seed(seed)
        np.random.seed(seed)
        
        # Calculate number of chains
        total_monomers = self.occupied_sites
        monomers_per_chain = chain_length
        total_chains = total_monomers // monomers_per_chain
        
        if total_chains % 2 != 0:
            total_chains -= 1  # Ensure equal number of A and B chains
        
        num_chains_per_type = total_chains // 2
        
        print(f"Generating system with:")
        print(f"  Chain length: {chain_length}")
        print(f"  Chains per type: {num_chains_per_type}")
        print(f"  Total atoms: {total_chains * chain_length}")
        print(f"  Case: {case}")
        print(f"  Polymer type: {polymer_type}")
        if case == 3 or case == 4:
            print(f"  Charge gap: {gap}")
        
        # Place chains on lattice
        atoms, bonds = self.place_chains_on_lattice(chain_length, num_chains_per_type, polymer_type)
        
        # Apply charges for charged cases
        if case == 3 or (case == 4 and gap > 0):
            self.apply_charges(atoms, gap, polymer_type)
        
        # Generate filename
        if case in [1, 2]:
            filename = f"data_{polymer_type}_case{case}_N{chain_length}_seed{seed}.lammps"
        else:
            filename = f"data_{polymer_type}_case{case}_N{chain_length}_gap{gap}_seed{seed}.lammps"
        
        # Write data file
        self.write_lammps_data(filename, atoms, bonds, case, chain_length)
        print(f"Data file written: {filename}")
        
        return filename

def main():
    parser = argparse.ArgumentParser(description='Generate LAMMPS data file for polymer simulation')
    parser.add_argument('--chain_length', type=int, required=True, choices=[12, 16, 20, 24, 28],
                       help='Chain length (N)')
    parser.add_argument('--case', type=int, required=True, choices=[1, 2, 3, 4],
                       help='Simulation case (1-4)')
    parser.add_argument('--gap', type=int, default=2, choices=[2, 3, 4],
                       help='Charge gap for charged cases (default: 2)')
    parser.add_argument('--polymer_type', type=str, default='homo', choices=['homo', 'block'],
                       help='Polymer type: homo or block (default: homo)')
    parser.add_argument('--seed', type=int, default=1,
                       help='Random seed (default: 1)')
    
    args = parser.parse_args()
    
    generator = PolymerDataGenerator()
    generator.generate_data_file(
        chain_length=args.chain_length,
        case=args.case,
        gap=args.gap,
        polymer_type=args.polymer_type,
        seed=args.seed
    )

if __name__ == "__main__":
    main()