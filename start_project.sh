#!/bin/bash

# Get the directory of the script
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
# SCRIPT_DIR = /Users/kensu/Desktop/SDP/Multilateration-3D-Positioning

# Start Blender in background mode and run the Python script
/Applications/Blender.app/Contents/MacOS/blender "$SCRIPT_DIR/threeD_positioning.blend" --background --python "$SCRIPT_DIR/tag_response.py"

# Run the Python script
# python "$SCRIPT_DIR/height_calculator.py"
python "$SCRIPT_DIR/multilateration.py"