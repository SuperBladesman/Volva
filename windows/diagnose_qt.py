"""Diagnose 'DLL load failed while importing QtCore'. Run inside the volva41 env:
    python windows\\diagnose_qt.py
and paste the output."""
import ctypes
import glob
import os
import sys
from ctypes import wintypes


def ver(path):
    try:
        v = ctypes.WinDLL("version")
        n = v.GetFileVersionInfoSizeW(path, None)
        if not n:
            return "?"
        buf = ctypes.create_string_buffer(n)
        v.GetFileVersionInfoW(path, 0, n, buf)
        p = ctypes.c_void_p()
        ln = wintypes.UINT()
        v.VerQueryValueW(buf, "\\", ctypes.byref(p), ctypes.byref(ln))
        d = ctypes.cast(p, ctypes.POINTER(wintypes.DWORD * 13)).contents
        return f"{d[2] >> 16}.{d[2] & 0xFFFF}.{d[3] >> 16}.{d[3] & 0xFFFF}"
    except Exception as e:
        return f"err {e}"


print("python", sys.version.split()[0], sys.executable)
try:
    import PySide6
    print("PySide6", PySide6.__version__, os.path.dirname(PySide6.__file__))
    from PySide6 import QtCore
    print("QtCore import in env: OK", QtCore.qVersion())
except Exception as e:
    print("QtCore import in env: FAILED ->", e)

pyside = os.path.join(sys.prefix, "Lib", "site-packages", "PySide6")
dist = os.path.join(os.environ.get("USERPROFILE", ""), "Volva_build", "dist", "Volva", "_internal")
places = {"env root": sys.prefix, "env Library\\bin": os.path.join(sys.prefix, "Library", "bin"),
          "PySide6": pyside, "shiboken6": os.path.join(sys.prefix, "Lib", "site-packages", "shiboken6"),
          "System32": os.path.join(os.environ.get("SystemRoot", r"C:\Windows"), "System32"),
          "dist _internal": dist, "dist PySide6": os.path.join(dist, "PySide6")}
for label, d in places.items():
    for pat in ("msvcp140*.dll", "vcruntime140*.dll", "python3*.dll", "Qt6Core.dll"):
        for f in sorted(glob.glob(os.path.join(d, pat))):
            print(f"  {label:15s} {os.path.basename(f):22s} {ver(f)}")
print("Qt6Core.dll / Qt5Core.dll on PATH:")
for d in os.environ.get("PATH", "").split(os.pathsep):
    for n in ("Qt6Core.dll", "Qt5Core.dll", "msvcp140.dll"):
        if d and os.path.isfile(os.path.join(d, n)):
            print("  ", os.path.join(d, n), ver(os.path.join(d, n)))

print("Qt / PySide / shiboken DLLs in the built app and in base Anaconda:")
base_ana = os.path.join(os.environ.get("USERPROFILE", ""), "anaconda3", "Library", "bin")
for root in (os.path.join(os.environ.get("USERPROFILE", ""), "Volva_build", "dist", "Volva"), base_ana):
    for dirpath, _, files in os.walk(root):
        for f in files:
            fl = f.lower()
            if fl.endswith(".dll") and (fl.startswith(("qt6", "pyside6", "shiboken6"))):
                print(f"  {os.path.join(dirpath, f)}  {ver(os.path.join(dirpath, f))}")
        if root == base_ana:
            break
