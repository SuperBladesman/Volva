#!/usr/bin/env python3
"""
Test script to debug Volva imports
Run this from your project root directory
"""

import sys
import os

# Add current directory to path
project_root = os.path.dirname(os.path.abspath(__file__))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

print(f"Testing imports from: {project_root}")
print(f"Python path: {sys.path[0]}")

# Test 1: Can we import the main volva package?
try:
    import volva
    print("✅ Main volva package imported")
    print(f"   Available attributes: {[attr for attr in dir(volva) if not attr.startswith('_')]}")
except Exception as e:
    print(f"❌ Main volva package failed: {e}")

# Test 2: Can we import core?
try:
    from volva import core
    print("✅ volva.core imported")
except Exception as e:
    print(f"❌ volva.core failed: {e}")

# Test 3: Can we import data processor?
try:
    from volva.data import XRDDataProcessor
    print("✅ XRDDataProcessor imported")
except Exception as e:
    print(f"❌ XRDDataProcessor failed: {e}")

# Test 4: Can we import interfaces?
try:
    from volva import interfaces
    print("✅ volva.interfaces imported")
    print(f"   Available: {[attr for attr in dir(interfaces) if not attr.startswith('_')]}")
except Exception as e:
    print(f"❌ volva.interfaces failed: {e}")

# Test 5: Can we import VolvaGUI?
try:
    from volva.interfaces.gui import VolvaGUI
    print("✅ VolvaGUI class imported")
except Exception as e:
    print(f"❌ VolvaGUI failed: {e}")

# Test 6: Can we import launch_gui?
try:
    from volva.interfaces.gui import launch_gui
    print("✅ launch_gui function imported")
except Exception as e:
    print(f"❌ launch_gui failed: {e}")

# Test 7: Can we import JupyterVolva?
try:
    from volva.interfaces.jupyter import JupyterVolva
    print("✅ JupyterVolva imported")
except Exception as e:
    print(f"❌ JupyterVolva failed: {e}")

# Test 8: Check if __init__.py files exist
init_files = [
    'volva/__init__.py',
    'volva/core/__init__.py',
    'volva/data/__init__.py',
    'volva/peaks/__init__.py',
    'volva/analysis/__init__.py',
    'volva/interfaces/__init__.py',
    'volva/utils/__init__.py'
]

print("\n📁 Checking __init__.py files:")
for init_file in init_files:
    if os.path.exists(init_file):
        print(f"✅ {init_file}")
    else:
        print(f"❌ {init_file} (MISSING)")

print("\n🔧 If you see missing __init__.py files, create them (they can be empty)")