"""
Command-line interface.

  python -m volva                         launch the GUI
  python -m volva calibrate LaB6.xy -o profile.json [--note "Aeris, 1/4° slit"]
  python -m volva analyse sample.xy --profile profile.json -o results/
  python -m volva analyse *.xy --profile profile.json -o results/     (batch)
"""

from __future__ import annotations

import argparse
import glob
import os
import sys
import warnings

from . import __version__
from . import constants as C


def _wl(args, pattern=None):
    from .pipeline import AnalysisSettings, resolve_wavelength
    s = AnalysisSettings(wavelength=args.wavelength)
    if args.wavelength == "Custom":
        s.custom_ka1 = args.ka1
        s.custom_ka2 = args.ka2 or args.ka1
        s.custom_ratio = args.ratio if args.ka2 else 0.0
    return resolve_wavelength(s, pattern)


def cmd_calibrate(args):
    from .instrument import calibrate_from_standard
    from .io import load_pattern
    pat = load_pattern(args.standard)
    wl = _wl(args, pat)
    prof, log, fits, bg = calibrate_from_standard(
        pat, wl, doublet=not args.no_doublet, standard=args.standard_type, name=args.name or "",
        note=args.note or "", bg_half_width=args.bg_width)
    print("\n".join(log))
    if prof is None:
        print("Calibration failed.", file=sys.stderr)
        return 2
    out = args.output or os.path.splitext(args.standard)[0] + "_profile.json"
    prof.save(out)
    print(f"\nSaved instrument profile: {out}")
    if args.figure:
        from matplotlib.figure import Figure
        from . import plots
        fig = Figure(figsize=(9, 5.5))
        plots.plot_instrument(fig, prof)
        fp = os.path.splitext(out)[0] + ".png"
        fig.savefig(fp, dpi=300, facecolor="white")
        print(f"Saved figure: {fp}")
    return 0


def cmd_analyse(args):
    from .instrument import InstrumentProfile
    from .io import LoadError, load_pattern
    from .pipeline import AnalysisSettings, run_analysis
    from .report import write_bundle
    prof = InstrumentProfile.load(args.profile) if args.profile else None
    if prof is None:
        print("WARNING: no instrument profile (--profile); sizes will be lower bounds.", file=sys.stderr)
    files = []
    for pat in args.samples:
        files.extend(sorted(glob.glob(pat)) or [pat])
    status = 0
    for fpath in files:
        try:
            pat = load_pattern(fpath)
        except LoadError as e:
            print(f"{fpath}: {e}", file=sys.stderr)
            status = 1
            continue
        s = AnalysisSettings(
            wavelength=args.wavelength, doublet=not args.no_doublet,
            size_convention=args.size, k_custom=args.k, bg_half_width=args.bg_width,
            prominence_sigma=args.prominence, min_r2=args.min_r2,
            two_theta_min=args.tt_min, two_theta_max=args.tt_max)
        if args.wavelength == "Custom":
            s.custom_ka1 = args.ka1
            s.custom_ka2 = args.ka2 or args.ka1
            s.custom_ratio = args.ratio if args.ka2 else 0.0
        res = run_analysis(pat, prof, s)
        paths = write_bundle(res, args.output, figures=not args.no_figures)
        print(f"{pat.name}: {res.interpretation.headline}")
        for w in res.warnings:
            print(f"   ⚠ {w}")
        print(f"   → {len(paths)} files in {args.output}")
    return status


def build_parser():
    p = argparse.ArgumentParser(prog="volva", description=f"Völva {__version__}: XRD crystallite size and strain")
    p.add_argument("--version", action="version", version=f"Völva {__version__}")
    sub = p.add_subparsers(dest="cmd")

    def common(sp):
        sp.add_argument("--wavelength", default="Cu Ka",
                        choices=list(C.WAVELENGTHS) + ["From file", "Custom"])
        sp.add_argument("--ka1", type=float, help="custom Kα1 (nm)")
        sp.add_argument("--ka2", type=float, help="custom Kα2 (nm)")
        sp.add_argument("--ratio", type=float, default=0.5, help="custom Kα2/Kα1 ratio")
        sp.add_argument("--no-doublet", action="store_true", help="do not fit Kα2 (monochromatic source)")
        sp.add_argument("--bg-width", type=float, default=None, help="SNIP half-width in degrees (default: auto)")

    c = sub.add_parser("calibrate", help="build an instrument profile from a standard")
    c.add_argument("standard")
    c.add_argument("-o", "--output")
    c.add_argument("--standard-type", default=C.DEFAULT_STANDARD, choices=list(C.STANDARDS))
    c.add_argument("--name")
    c.add_argument("--note", help="optics / configuration description stored in the profile")
    c.add_argument("--figure", action="store_true")
    common(c)

    a = sub.add_parser("analyse", aliases=["analyze"], help="analyse one or more sample patterns")
    a.add_argument("samples", nargs="+")
    a.add_argument("--profile", help="instrument profile JSON from 'calibrate'")
    a.add_argument("-o", "--output", default="volva_results")
    a.add_argument("--size", default=C.DEFAULT_SIZE_CONVENTION, choices=list(C.SIZE_CONVENTIONS))
    a.add_argument("--k", type=float, default=None, help="custom shape constant (overrides --size)")
    a.add_argument("--prominence", type=float, default=C.DEFAULT_PROMINENCE_SIGMA)
    a.add_argument("--min-r2", type=float, default=C.DEFAULT_MIN_R2)
    a.add_argument("--tt-min", type=float)
    a.add_argument("--tt-max", type=float)
    a.add_argument("--no-figures", action="store_true")
    common(a)
    return p


def _safe_streams():
    """Never crash on a character the output cannot encode.

    On Windows, output sent to a pipe or a file is encoded as cp1252 by default,
    which has no →, ≈ or η. Files and pipes get UTF-8; a console keeps its own
    encoding and shows anything it cannot encode as "?".
    """
    for stream in (sys.stdout, sys.stderr):
        try:
            if stream.isatty():
                stream.reconfigure(errors="replace")
            else:
                stream.reconfigure(encoding="utf-8", errors="replace")
        except (AttributeError, ValueError, OSError):
            pass


def main(argv=None):
    _safe_streams()
    warnings.filterwarnings("ignore", category=UserWarning, module="uncertainties")
    argv = sys.argv[1:] if argv is None else argv
    if not argv:
        from .gui.app import main as gui_main
        return gui_main()
    args = build_parser().parse_args(argv)
    if args.cmd == "calibrate":
        return cmd_calibrate(args)
    if args.cmd in ("analyse", "analyze"):
        return cmd_analyse(args)
    build_parser().print_help()
    return 0
