#!/bin/bash
# Run NextTissUe simulation with output directly to Windows-accessible directory

OUTPUT_DIR="/mnt/d/NextTissUe_data"

echo "========================================="
echo "  NextTissUe Simulation Runner"
echo "========================================="
echo "Output directory: $OUTPUT_DIR"
echo ""

# Ensure output directory exists
if [ ! -d "$OUTPUT_DIR" ]; then
    echo "Creating output directory..."
    mkdir -p "$OUTPUT_DIR"
fi

# Move to working directory
cd /home/lost_octave/LAMMPS/NextTissUe

# Set up symbolic links for easy access
echo "Setting up output links..."
ln -sf "$OUTPUT_DIR" output_data

echo "Starting LAMMPS simulation..."
echo "Running with 4 MPI processes"
echo "Output will be written directly to: $OUTPUT_DIR"
echo ""
echo "Press Ctrl+C to stop the simulation"
echo "To run in background: ./run_simulation.sh &"
echo ""

# Run simulation
mpirun -np 4 /home/lost_octave/LAMMPS/build/lmp -in in.demo

echo ""
echo "Simulation completed or stopped"
echo "Results are in: $OUTPUT_DIR"

