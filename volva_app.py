"""Entry point for the frozen (PyInstaller) build and for double-click launching."""
import multiprocessing
import sys

if __name__ == "__main__":
    multiprocessing.freeze_support()
    from volva.cli import main
    sys.exit(main())
