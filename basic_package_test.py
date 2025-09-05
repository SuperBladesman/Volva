#!/usr/bin/env python3
"""
Basic test to identify the specific import issue
"""

import sys
import os

# Add current directory to path
project_root = os.path.dirname(os.path.abspath(__file__))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

print(f"Testing from: {project_root}")

# Test 1: Check if __init__.py files exist
print("\n=== Checking __init__.py files ===")
init_files = [
    'volva/__init__.py',
    'volva/core/__init__.py',
    'volva/data/__init__.py',
    'volva/peaks/__init__.py',
    'volva/analysis/__init__.py',
    'volva/interfaces/__init__.py',
    'volva/utils/__init__.py'
]

missing_files = []
for init_file in init_files:
    full_path = os.path.join(project_root, init_file)
    if os.path.exists(full_path):
        print(f"Found: {init_file}")
    else:
        print(f"MISSING: {init_file}")
        missing_files.append(full_path)

# Create missing __init__.py files
if missing_files:
    print(f"\n=== Creating {len(missing_files)} missing __init__.py files ===")
    for missing_file in missing_files:
        try:
            os.makedirs(os.path.dirname(missing_file), exist_ok=True)
            with open(missing_file, 'w') as f:
                f.write('# Volva package file\n')
            print(f"Created: {missing_file}")
        except Exception as e:
            print(f"Failed to create {missing_file}: {e}")

# Test 2: Try importing volva step by step
print("\n=== Testing imports step by step ===")

try:
    print("Step 1: Importing volva...")
    import volva
    print("SUCCESS: volva imported")
except Exception as e:
    print(f"FAILED: volva - {e}")
    print("Cannot continue without basic volva import")
    sys.exit(1)

try:
    print("Step 2: Importing volva.core...")
    from volva import core
    print("SUCCESS: volva.core imported")
except Exception as e:
    print(f"FAILED: volva.core - {e}")

try:
    print("Step 3: Importing volva.interfaces...")
    from volva import interfaces
    print("SUCCESS: volva.interfaces imported")
except Exception as e:
    print(f"FAILED: volva.interfaces - {e}")

# Test 3: Try importing VolvaGUI without importing the whole gui module
try:
    print("Step 4: Importing VolvaGUI class...")
    from volva.interfaces.gui import VolvaGUI
    print("SUCCESS: VolvaGUI class imported")
except Exception as e:
    print(f"FAILED: VolvaGUI class - {e}")
    print("This suggests an issue in the gui.py file itself")

try:
    print("Step 5: Importing launch_gui function...")
    from volva.interfaces.gui import launch_gui
    print("SUCCESS: launch_gui function imported")
except Exception as e:
    print(f"FAILED: launch_gui function - {e}")

print("\n=== Test complete ===")