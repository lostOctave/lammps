#!/bin/bash
# Continuous sync script to copy NextTissUe data to Windows-accessible directory

SOURCE_DIR="/home/lost_octave/LAMMPS/NextTissUe"
DEST_DIR="/mnt/d/NextTissUe_data"

echo "=== NextTissUe Data Sync Script ==="
echo "Syncing from: $SOURCE_DIR"
echo "Syncing to: $DEST_DIR"
echo "Started at: $(date)"
echo ""

# Function to sync data
sync_data() {
    # Sync dump files
    rsync -av --progress "$SOURCE_DIR/dumps/" "$DEST_DIR/dumps/" 2>/dev/null
    
    # Sync visualization files if they exist
    if [ -f "$SOURCE_DIR/cell_snapshot_initial.png" ]; then
        cp -u "$SOURCE_DIR/cell_snapshot_initial.png" "$DEST_DIR/"
    fi
    
    if [ -f "$SOURCE_DIR/cell_snapshot_final.png" ]; then
        cp -u "$SOURCE_DIR/cell_snapshot_final.png" "$DEST_DIR/"
    fi
    
    if [ -f "$SOURCE_DIR/cells_animation.mp4" ]; then
        cp -u "$SOURCE_DIR/cells_animation.mp4" "$DEST_DIR/"
    fi
    
    # Sync log file
    if [ -f "$SOURCE_DIR/log.lammps" ]; then
        cp -u "$SOURCE_DIR/log.lammps" "$DEST_DIR/"
    fi
    
    # Get file sizes
    DUMPS_SIZE=$(du -sh "$DEST_DIR/dumps" 2>/dev/null | cut -f1)
    TOTAL_SIZE=$(du -sh "$DEST_DIR" 2>/dev/null | cut -f1)
    
    echo "[$(date +%H:%M:%S)] Synced - Dumps: $DUMPS_SIZE, Total: $TOTAL_SIZE"
}

# Check if running in continuous mode or one-time sync
if [ "$1" == "--once" ]; then
    echo "Running one-time sync..."
    sync_data
    echo "Sync complete!"
else
    echo "Running continuous sync (Ctrl+C to stop)"
    echo "Syncing every 60 seconds..."
    echo ""
    
    while true; do
        sync_data
        sleep 60
    done
fi

