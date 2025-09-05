#!/usr/bin/env python3
"""
Volva XRD Analysis Tool - Main Runner
Launch the GUI or create demo data
"""

import sys
import os

# Add the project directory to Python path
project_root = os.path.dirname(os.path.abspath(__file__))
if project_root not in sys.path:
    sys.path.insert(0, project_root)


def main():
    """Main entry point with better error handling"""
    try:
        # Import Volva with error handling
        print("🔬 Loading Volva XRD Analysis Tool...")

        try:
            # First try the standard import
            from volva.interfaces.gui import launch_gui
            from volva.interfaces.jupyter import JupyterVolva
            print("✅ Volva loaded successfully!")
        except ImportError as e:
            print(f"❌ Standard import failed: {e}")
            print("🔧 Trying alternative import method...")

            # Try importing from the package directly
            try:
                import volva
                launch_gui = volva.launch_gui
                JupyterVolva = volva.JupyterVolva
                print("✅ Volva loaded with package import!")
            except Exception as e2:
                print(f"❌ Package import failed: {e2}")
                print("🔧 Trying direct module import...")

                # Last resort - direct module import
                volva_path = os.path.join(project_root, 'volva')
                if volva_path not in sys.path:
                    sys.path.insert(0, volva_path)

                from interfaces.gui import launch_gui
                from interfaces.jupyter import JupyterVolva
                print("✅ Volva loaded with direct import!")

    except Exception as e:
        print(f"❌ Failed to load Volva: {e}")
        print(f"\n🔧 Debugging info:")
        print(f"Project root: {project_root}")
        print(f"Python path: {sys.path[:3]}...")
        print(f"Volva directory exists: {os.path.exists(os.path.join(project_root, 'volva'))}")
        print(f"GUI file exists: {os.path.exists(os.path.join(project_root, 'volva', 'interfaces', 'gui.py'))}")

        # Try to identify the specific issue
        try:
            from volva.interfaces import gui
            print("✅ GUI module can be imported")
            if hasattr(gui, 'launch_gui'):
                print("✅ launch_gui function exists")
            else:
                print("❌ launch_gui function is missing from gui module")
                print(f"Available attributes: {[attr for attr in dir(gui) if not attr.startswith('_')]}")
        except Exception as debug_e:
            print(f"❌ Cannot import GUI module: {debug_e}")

        return

    print("=" * 50)
    print("Choose an option:")
    print("1. Launch GUI")
    print("2. Create sample data and launch GUI")
    print("3. Jupyter interface demo")
    print("4. Exit")

    while True:
        try:
            choice = input("\nEnter your choice (1-4): ").strip()

            if choice == '1':
                print("🚀 Launching GUI...")
                launch_gui()
                break
            elif choice == '2':
                print("📊 Creating sample data...")
                create_sample_data()
                print("🚀 Launching GUI with sample data...")
                launch_gui()
                break
            elif choice == '3':
                print("📓 Running Jupyter demo...")
                jupyter_demo(JupyterVolva)
                break
            elif choice == '4':
                print("👋 Goodbye!")
                break
            else:
                print("❌ Invalid choice. Please enter 1-4.")
        except (KeyboardInterrupt, EOFError):
            print("\n👋 Goodbye!")
            break


def create_sample_data():
    """Create sample XRD data for testing"""
    import numpy as np

    # Generate realistic XRD pattern
    two_theta = np.linspace(20, 80, 2000)

    # Background
    background = 200 + 50 * np.exp(-(two_theta - 30) ** 2 / 1000)

    # Multiple peaks with progressive broadening
    peak_data = [
        {'center': 28.4, 'intensity': 1200, 'fwhm': 0.08},
        {'center': 40.5, 'intensity': 900, 'fwhm': 0.12},
        {'center': 50.3, 'intensity': 700, 'fwhm': 0.15},
        {'center': 58.3, 'intensity': 600, 'fwhm': 0.18},
        {'center': 69.1, 'intensity': 400, 'fwhm': 0.22},
    ]

    intensity = background.copy()

    for peak in peak_data:
        sigma = peak['fwhm'] / 2.355
        gaussian = peak['intensity'] * np.exp(-(two_theta - peak['center']) ** 2 / (2 * sigma ** 2))
        intensity += gaussian

    # Add realistic noise
    np.random.seed(42)  # Reproducible results
    noise = np.random.normal(0, np.sqrt(np.maximum(intensity, 1)) * 0.3)
    intensity += noise
    intensity = np.maximum(intensity, background * 0.8)

    # Save sample data
    filename = "sample_xrd_data.txt"
    with open(filename, 'w') as f:
        f.write("# Sample XRD data for Volva v4.0 testing\n")
        f.write("# Generated with 5 peaks showing progressive broadening\n")
        f.write("# 2theta(deg)    Intensity(counts)\n")
        for tt, intens in zip(two_theta, intensity):
            f.write(f"{tt:.4f}    {intens:.1f}\n")

    print(f"✅ Sample data saved to: {filename}")
    print(f"📊 Data contains {len(two_theta)} points from {two_theta.min():.1f}° to {two_theta.max():.1f}°")


def jupyter_demo(JupyterVolva):
    """Quick Jupyter interface demo"""
    try:
        print("🔬 Creating Jupyter analyzer...")
        analyzer = JupyterVolva('Cu Ka')

        # Create and load sample data if it doesn't exist
        if not os.path.exists('sample_xrd_data.txt'):
            create_sample_data()

        print("📁 Loading sample data...")
        success = analyzer.load_data('sample_xrd_data.txt')

        if success:
            print("🔍 Detecting peaks...")
            peaks = analyzer.detect_peaks(prominence_factor=0.15)

            if peaks:
                print(f"✅ Found {len(peaks)} peaks")
                print("🧮 Running Scherrer analysis...")
                result = analyzer.run_scherrer()

                if result:
                    print(f"📊 Results:")
                    print(f"   Average crystallite size: {result['average_size']:.1f} ± {result['std_size']:.1f} nm")
                    print(f"   Median size: {result['median_size']:.1f} nm")
                    print(f"   Peaks used: {result['num_peaks']}")
                    print("💡 For comprehensive analysis, use the GUI!")
                else:
                    print("❌ Analysis failed")
            else:
                print("❌ No peaks detected - try adjusting parameters in GUI")
        else:
            print("❌ Failed to load sample data")

    except Exception as e:
        print(f"❌ Demo failed: {e}")
        print("💡 Try the GUI instead: choose option 1")


if __name__ == "__main__":
    main()