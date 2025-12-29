#!/bin/bash
# Check status of LAMMPS simulation

OUTPUT_DIR="/mnt/d/NextTissUe_data"

echo "========================================="
echo "  NextTissUe Status Check"
echo "========================================="
echo ""

# Check LAMMPS simulation
echo "🔬 LAMMPS Simulation:"
LAMMPS_RUNNING=$(ps aux | grep "lmp -in in.demo" | grep -v grep | wc -l)
if [ $LAMMPS_RUNNING -gt 0 ]; then
    echo "   ✅ Running (4 MPI processes)"
    # Get last few lines from log
    if [ -f "log.lammps" ]; then
        LAST_STEP=$(tail -20 log.lammps | grep "^[[:space:]]*[0-9]" | tail -1 | awk '{print $1}')
        if [ ! -z "$LAST_STEP" ]; then
            echo "   📊 Last timestep: $LAST_STEP"
        fi
    fi
else
    echo "   ❌ Not running"
    echo "   To start: ./run_simulation.sh &"
fi
echo ""

# Check output configuration
echo "💾 Output Configuration:"
echo "   📂 Direct output to: $OUTPUT_DIR"
echo "   ✅ No local copying - saves disk space!"
echo ""

# Check data sizes
echo "📊 Data Sizes:"
if [ -d "$OUTPUT_DIR" ]; then
    WINDOWS_TOTAL=$(du -sh $OUTPUT_DIR 2>/dev/null | cut -f1)
    echo "   Windows D:\\NextTissUe_data: $WINDOWS_TOTAL"
    
    # Count dump files
    DUMP_COUNT=$(find $OUTPUT_DIR -name "*.dat" 2>/dev/null | wc -l)
    if [ $DUMP_COUNT -gt 0 ]; then
        echo "   Dump files: $DUMP_COUNT"
    fi
fi

if [ -f "log.lammps" ]; then
    LOG_SIZE=$(du -sh log.lammps | cut -f1)
    echo "   log.lammps: $LOG_SIZE"
fi
echo ""

# File checks
echo "📄 Files in Windows folder:"
if [ -d "$OUTPUT_DIR" ]; then
    if [ -f "$OUTPUT_DIR/cell_snapshot_initial.png" ]; then
        echo "   ✅ cell_snapshot_initial.png"
    fi
    if [ -f "$OUTPUT_DIR/cell_snapshot_final.png" ]; then
        echo "   ✅ cell_snapshot_final.png"
    fi
    if [ -f "$OUTPUT_DIR/cells_animation.mp4" ]; then
        echo "   ✅ cells_animation.mp4"
    fi
    if [ -f "$OUTPUT_DIR/README.txt" ]; then
        echo "   ✅ README.txt"
    fi
    if [ -f "$OUTPUT_DIR/view_results.html" ]; then
        echo "   ✅ view_results.html"
    fi
fi
echo ""

echo "========================================="
echo "Commands:"
echo "  Start simulation: ./run_simulation.sh &"
echo "  View log tail: tail -f log.lammps"
echo "  Stop simulation: pkill -f 'lmp -in in.demo'"
echo "  Windows access: D:\\NextTissUe_data\\"
echo "========================================="
