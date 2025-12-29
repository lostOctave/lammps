#!/usr/bin/env python3
"""
Simple visualization script for NextTissUe LAMMPS output
Usage: python visualize_cells.py [dump_file]
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import glob
import os

def parse_lammps_dump(filename, timestep=None):
    """
    Parse LAMMPS custom dump file
    Returns: dict with timesteps as keys, each containing numpy array of atom data
    """
    timesteps = {}
    current_timestep = None
    atoms = []
    reading_atoms = False
    headers = []
    
    print(f"Reading file: {filename}")
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            
            if 'ITEM: TIMESTEP' in line:
                # Save previous timestep data
                if current_timestep is not None and atoms:
                    timesteps[current_timestep] = np.array(atoms)
                    atoms = []
                
                current_timestep = int(next(f).strip())
                reading_atoms = False
                
                # If looking for specific timestep and we passed it, break
                if timestep is not None and current_timestep > timestep:
                    break
                    
            elif 'ITEM: NUMBER OF ATOMS' in line:
                n_atoms = int(next(f).strip())
                
            elif 'ITEM: ATOMS' in line:
                headers = line.split()[2:]  # Get column names
                reading_atoms = True
                
            elif reading_atoms and line:
                try:
                    atom_data = [float(x) if '.' in x or 'e' in x.lower() else int(x) 
                                for x in line.split()]
                    atoms.append(atom_data)
                except:
                    pass
    
    # Save last timestep
    if current_timestep is not None and atoms:
        timesteps[current_timestep] = np.array(atoms)
    
    print(f"Found {len(timesteps)} timesteps")
    return timesteps, headers

def plot_snapshot(data, headers, timestep=0, output_file='cell_snapshot.png'):
    """
    Plot a single snapshot of the cell configuration
    """
    # Find column indices
    try:
        x_idx = headers.index('x')
        y_idx = headers.index('y')
        mol_idx = headers.index('mol')
    except ValueError:
        # Default positions if headers not found
        x_idx = 1
        y_idx = 2
        mol_idx = 3
        print(f"Using default column indices. Available headers: {headers}")
    
    x = data[:, x_idx]
    y = data[:, y_idx]
    mol = data[:, mol_idx]
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, 12))
    
    # Plot cells colored by molecule ID
    scatter = ax.scatter(x, y, c=mol, s=50, cmap='tab20', 
                        alpha=0.7, edgecolors='black', linewidth=0.5)
    
    ax.set_xlabel('x', fontsize=14)
    ax.set_ylabel('y', fontsize=14)
    ax.set_title(f'Cell Configuration - Timestep {timestep}', fontsize=16)
    ax.axis('equal')
    ax.grid(True, alpha=0.3)
    
    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('Molecule ID', fontsize=12)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"Saved figure to {output_file}")
    
    return fig, ax

def create_animation(timesteps_data, headers, output_file='cells_animation.mp4'):
    """
    Create an animation from multiple timesteps (requires ffmpeg)
    """
    try:
        from matplotlib.animation import FuncAnimation, FFMpegWriter
    except ImportError:
        print("Animation requires matplotlib.animation. Install ffmpeg for video export.")
        return
    
    timesteps = sorted(timesteps_data.keys())
    
    # Setup figure
    fig, ax = plt.subplots(figsize=(12, 12))
    
    # Find column indices
    try:
        x_idx = headers.index('x')
        y_idx = headers.index('y')
        mol_idx = headers.index('mol')
    except ValueError:
        x_idx = 1
        y_idx = 2
        mol_idx = 3
    
    # Get data range for consistent axes
    all_data = np.vstack([timesteps_data[t] for t in timesteps])
    x_min, x_max = all_data[:, x_idx].min(), all_data[:, x_idx].max()
    y_min, y_max = all_data[:, y_idx].min(), all_data[:, y_idx].max()
    margin = 0.1 * max(x_max - x_min, y_max - y_min)
    
    def update(frame_num):
        ax.clear()
        timestep = timesteps[frame_num]
        data = timesteps_data[timestep]
        
        x = data[:, x_idx]
        y = data[:, y_idx]
        mol = data[:, mol_idx]
        
        ax.scatter(x, y, c=mol, s=50, cmap='tab20', 
                  alpha=0.7, edgecolors='black', linewidth=0.5)
        ax.set_xlim(x_min - margin, x_max + margin)
        ax.set_ylim(y_min - margin, y_max + margin)
        ax.set_xlabel('x', fontsize=14)
        ax.set_ylabel('y', fontsize=14)
        ax.set_title(f'Timestep {timestep}', fontsize=16)
        ax.axis('equal')
        ax.grid(True, alpha=0.3)
    
    anim = FuncAnimation(fig, update, frames=len(timesteps), interval=200)
    
    try:
        writer = FFMpegWriter(fps=5, bitrate=1800)
        anim.save(output_file, writer=writer)
        print(f"Saved animation to {output_file}")
    except Exception as e:
        print(f"Could not save animation: {e}")
        print("Try installing ffmpeg: sudo apt-get install ffmpeg")
    
    plt.close()

def main():
    # Find dump files
    if len(sys.argv) > 1:
        dump_file = sys.argv[1]
    else:
        # Look for dump files
        dump_files = glob.glob('dumps/*.dat')
        if not dump_files:
            dump_files = glob.glob('*.dat')
        
        if not dump_files:
            print("No dump files found!")
            print("Usage: python visualize_cells.py [dump_file]")
            return
        
        # Use the most recent file
        dump_file = max(dump_files, key=os.path.getmtime)
        print(f"Using most recent dump file: {dump_file}")
    
    if not os.path.exists(dump_file):
        print(f"File not found: {dump_file}")
        return
    
    # Parse dump file
    timesteps_data, headers = parse_lammps_dump(dump_file)
    
    if not timesteps_data:
        print("No data found in dump file!")
        return
    
    # Plot first and last timesteps
    timesteps = sorted(timesteps_data.keys())
    
    print(f"\nCreating snapshots...")
    plot_snapshot(timesteps_data[timesteps[0]], headers, 
                 timestep=timesteps[0], 
                 output_file='cell_snapshot_initial.png')
    
    if len(timesteps) > 1:
        plot_snapshot(timesteps_data[timesteps[-1]], headers, 
                     timestep=timesteps[-1], 
                     output_file='cell_snapshot_final.png')
    
    # Create animation if multiple timesteps
    if len(timesteps) > 1:
        print(f"\nCreating animation with {len(timesteps)} frames...")
        create_animation(timesteps_data, headers, 'cells_animation.mp4')
    
    print("\nVisualization complete!")
    print(f"Generated images: cell_snapshot_initial.png", end="")
    if len(timesteps) > 1:
        print(", cell_snapshot_final.png")
    print()

if __name__ == '__main__':
    main()

