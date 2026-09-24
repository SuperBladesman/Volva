"""
Loading diffraction patterns.

Supported: two- or three-column text (.xy, .xye, .txt, .dat, .csv, .asc) with
any header lines, comma/semicolon/tab/space delimiters and comma decimals;
PANalytical XRDML (.xrdml), including the wavelength recorded by the
instrument.
"""

from __future__ import annotations

import hashlib
import os
import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from typing import Optional

import numpy as np


@dataclass
class Pattern:
    two_theta: np.ndarray
    intensity: np.ndarray
    path: str = ""
    sha256: str = ""
    esd: Optional[np.ndarray] = None
    wavelength: Optional[dict] = None     # {'ka1','ka2','ratio'} in nm, if the file records it
    meta: dict = field(default_factory=dict)

    @property
    def name(self) -> str:
        return os.path.basename(self.path) if self.path else "pattern"

    @property
    def step(self) -> float:
        return float(np.median(np.diff(self.two_theta)))


class LoadError(Exception):
    pass


def file_sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_pattern(path: str) -> Pattern:
    """Load a pattern from disk. Raises LoadError with a readable message."""
    if not os.path.isfile(path):
        raise LoadError(f"File not found: {path}")
    ext = os.path.splitext(path)[1].lower()
    if ext == ".xrdml":
        pat = _load_xrdml(path)
    elif ext in (".raw", ".brml", ".rd", ".sd"):
        raise LoadError(
            f"{ext} is a binary instrument format. Export the scan as .xy or .xrdml first.")
    else:
        pat = _load_text(path)
    pat.path = os.path.abspath(path)
    pat.sha256 = file_sha256(path)
    _validate(pat)
    return pat


# ---------------------------------------------------------------------------
# Text formats
# ---------------------------------------------------------------------------
_NUM = r"[-+]?(?:\d+[.,]?\d*|[.,]\d+)(?:[eE][-+]?\d+)?"


def _split(line: str):
    line = line.strip()
    if not line:
        return None
    # Delimiter detection: semicolon or tab first (these allow comma decimals)
    if ";" in line:
        parts = line.split(";")
        parts = [p.replace(",", ".") for p in parts]
    elif "\t" in line:
        parts = line.split("\t")
        parts = [p.replace(",", ".") if p.count(",") == 1 and "." not in p else p for p in parts]
    elif "," in line and "." not in line and re.search(r"\s", line):
        # "12,5  300" : comma decimals, whitespace delimiter
        parts = line.replace(",", ".").split()
    elif "," in line and not re.search(r"\s", line):
        parts = line.split(",")
    elif "," in line:
        # "12.3, 456" style
        parts = re.split(r"[,\s]+", line)
    else:
        parts = line.split()
    out = []
    for p in parts:
        p = p.strip()
        if p == "":
            continue
        try:
            out.append(float(p))
        except ValueError:
            return None
    return out if len(out) >= 2 else None


def _load_text(path: str) -> Pattern:
    rows = []
    skipped = 0
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for raw in f:
            s = raw.strip()
            if not s or s[0] in "#!;%/*[":
                skipped += 1
                continue
            vals = _split(s)
            if vals is None:
                skipped += 1
                continue
            rows.append(vals[:3])
    if len(rows) < 20:
        raise LoadError(
            f"Could not find at least 20 numeric data rows in {os.path.basename(path)} "
            f"({skipped} non-numeric lines skipped).")
    ncol = min(len(r) for r in rows)
    data = np.array([r[:ncol] for r in rows], float)
    esd = data[:, 2] if ncol >= 3 else None
    return Pattern(two_theta=data[:, 0], intensity=data[:, 1], esd=esd,
                   meta={"skipped_lines": skipped, "format": "text"})


# ---------------------------------------------------------------------------
# XRDML
# ---------------------------------------------------------------------------
def _strip_ns(tag: str) -> str:
    return tag.split("}", 1)[-1]


def _load_xrdml(path: str) -> Pattern:
    try:
        root = ET.parse(path).getroot()
    except ET.ParseError as e:
        raise LoadError(f"XRDML could not be parsed: {e}") from e

    def find_all(node, name):
        return [el for el in node.iter() if _strip_ns(el.tag) == name]

    wl = None
    for uw in find_all(root, "usedWavelength"):
        vals = {}
        for ch in uw:
            t = _strip_ns(ch.tag)
            try:
                vals[t] = float(ch.text)
            except (TypeError, ValueError):
                pass
        if "kAlpha1" in vals:
            # XRDML stores Angstrom
            wl = {"ka1": vals["kAlpha1"] / 10.0,
                  "ka2": vals.get("kAlpha2", vals["kAlpha1"]) / 10.0,
                  "ratio": vals.get("ratioKAlpha2KAlpha1", 0.5)}
        break

    dps = find_all(root, "dataPoints")
    if not dps:
        raise LoadError("XRDML contains no dataPoints block.")
    tt_all, I_all = [], []
    for dp in dps:
        pos = None
        counting_time = None
        for ch in dp:
            t = _strip_ns(ch.tag)
            if t == "positions" and ch.get("axis") == "2Theta":
                lp = [c for c in ch if _strip_ns(c.tag) == "listPositions"]
                sp = [c for c in ch if _strip_ns(c.tag) == "startPosition"]
                ep = [c for c in ch if _strip_ns(c.tag) == "endPosition"]
                if lp:
                    pos = np.array(lp[0].text.split(), float)
                elif sp and ep:
                    pos = (float(sp[0].text), float(ep[0].text))
            elif t == "commonCountingTime":
                try:
                    counting_time = float(ch.text)
                except (TypeError, ValueError):
                    pass
        inten = None
        for ch in dp:
            t = _strip_ns(ch.tag)
            if t in ("counts", "intensities"):
                inten = np.array(ch.text.split(), float)
                unit = (ch.get("unit") or "").lower()
                if t == "intensities" and counting_time and ("cps" in unit or "per second" in unit):
                    inten = inten * counting_time   # cps -> counts (SNIP works in counting statistics)
                break
        if inten is None or pos is None:
            continue
        if isinstance(pos, tuple):
            pos = np.linspace(pos[0], pos[1], len(inten))
        tt_all.append(pos)
        I_all.append(inten)
    if not tt_all:
        raise LoadError("XRDML has no 2Theta scan with counts.")
    tt = np.concatenate(tt_all)
    I = np.concatenate(I_all)
    return Pattern(two_theta=tt, intensity=I, wavelength=wl,
                   meta={"format": "xrdml"})


# ---------------------------------------------------------------------------
def _validate(p: Pattern):
    tt, I = p.two_theta, p.intensity
    ok = np.isfinite(tt) & np.isfinite(I)
    tt, I = tt[ok], I[ok]
    if np.any(tt <= 0) or np.any(tt >= 180):
        raise LoadError("First column does not look like 2θ in degrees (values outside 0–180).")
    order = np.argsort(tt)
    tt, I = tt[order], I[order]
    # Remove duplicate 2theta values
    keep = np.concatenate([[True], np.diff(tt) > 0])
    tt, I = tt[keep], I[keep]
    if p.esd is not None:
        p.esd = p.esd[ok][order][keep]
    if len(tt) < 20:
        raise LoadError("Fewer than 20 valid data points.")
    steps = np.diff(tt)
    if np.max(steps) > 5 * np.median(steps):
        p.meta["warning_step"] = (
            f"Non-uniform step: largest gap {np.max(steps):.3f}° vs median {np.median(steps):.4f}°.")
    if np.median(steps) > 0.05:
        p.meta["warning_coarse"] = (
            f"Coarse step ({np.median(steps):.3f}°). Sharp peaks may be under-sampled; "
            "0.02° or finer is recommended for line-profile analysis.")
    if np.min(I) < 0:
        p.meta["warning_negative"] = "Negative intensities present (already background-subtracted?)."
    p.two_theta, p.intensity = tt, I
