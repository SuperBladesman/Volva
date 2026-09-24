# PyInstaller spec for Völva 4.1 (one-folder, windowed).
# Build via build_windows.bat, or:  pyinstaller --noconfirm volva.spec
import os
import sys

from PyInstaller.utils.hooks import collect_submodules

hidden = collect_submodules("lmfit") + collect_submodules("asteval") + collect_submodules("uncertainties")
hidden += ["matplotlib.backends.backend_qtagg"]

a = Analysis(
    ["volva_app.py"],
    pathex=["."],
    datas=[("volva/gui/volva.ico", "volva/gui"), ("examples", "examples"),
           ("README.md", "."), ("VALIDATION.md", "."), ("CHANGELOG.md", ".")],
    hiddenimports=hidden,
    excludes=["tkinter", "PyQt5", "PyQt6", "IPython", "notebook", "pytest"],
)


# ---------------------------------------------------------------------------
# Windows / conda: make sure the NEWEST Microsoft C++ runtime is bundled.
# Conda's Python ships an older msvcp140.dll; Qt 6 needs a newer one, and if the
# old copy is bundled the app fails with "DLL load failed while importing
# QtCore: The specified procedure could not be found".
# ---------------------------------------------------------------------------
def _file_version(path):
    try:
        import ctypes
        from ctypes import wintypes
        ver = ctypes.WinDLL("version")
        size = ver.GetFileVersionInfoSizeW(path, None)
        if not size:
            return (0, 0, 0, 0)
        buf = ctypes.create_string_buffer(size)
        ver.GetFileVersionInfoW(path, 0, size, buf)
        p = ctypes.c_void_p()
        n = wintypes.UINT()
        ver.VerQueryValueW(buf, "\\", ctypes.byref(p), ctypes.byref(n))

        class VS_FIXEDFILEINFO(ctypes.Structure):
            _fields_ = [(f, wintypes.DWORD) for f in (
                "sig", "strucver", "fvMS", "fvLS", "pvMS", "pvLS", "mask", "flags", "os", "type", "sub", "dMS", "dLS")]

        info = ctypes.cast(p, ctypes.POINTER(VS_FIXEDFILEINFO)).contents
        return (info.fvMS >> 16, info.fvMS & 0xFFFF, info.fvLS >> 16, info.fvLS & 0xFFFF)
    except Exception:
        return (0, 0, 0, 0)


if sys.platform == "win32":
    import PySide6
    import shiboken6
    runtime = {"msvcp140.dll", "msvcp140_1.dll", "msvcp140_2.dll", "vcruntime140.dll",
               "vcruntime140_1.dll", "concrt140.dll"}
    search = [os.path.dirname(PySide6.__file__), os.path.dirname(shiboken6.__file__), sys.base_prefix,
              os.path.join(sys.base_prefix, "Library", "bin"),
              os.path.join(os.environ.get("SystemRoot", r"C:\Windows"), "System32")]
    new_bins = []
    for entry in a.binaries:
        dest, src, kind = entry
        name = os.path.basename(dest).lower()
        if name in runtime and os.path.dirname(dest) in ("", "."):
            cands = [src] + [os.path.join(d, os.path.basename(dest)) for d in search]
            cands = [c for c in cands if os.path.isfile(c)]
            # newest version wins; on a tie prefer the copy shipped with PySide6 (built for its Qt)
            best = max(cands, key=lambda c: (_file_version(c), "pyside6" in c.lower()))
            print(f"volva.spec: {os.path.basename(dest)} <- {best} {_file_version(best)} "
                  f"(was {src} {_file_version(src)})")
            entry = (dest, best, kind)
        new_bins.append(entry)
    a.binaries = new_bins

# ---------------------------------------------------------------------------
# Never bundle DLLs from OUTSIDE the build environment. If another Python/Qt
# install (e.g. base Anaconda with Qt 6.9) is on PATH, PyInstaller can pick up
# its pyside6/shiboken6/Qt DLLs, which then clash with PySide6 6.11 at run time
# ("The specified procedure could not be found").
# ---------------------------------------------------------------------------
if sys.platform == "win32":
    env = os.path.normcase(os.path.abspath(sys.prefix))
    sysroot = os.path.normcase(os.environ.get("SystemRoot", r"C:\Windows"))
    env_dirs = [os.path.dirname(PySide6.__file__), os.path.dirname(shiboken6.__file__), sys.prefix,
                os.path.join(sys.prefix, "Library", "bin"), os.path.join(sys.prefix, "DLLs")]
    fixed = []
    for dest, src, kind in a.binaries:
        s_n = os.path.normcase(os.path.abspath(src))
        if src.lower().endswith((".dll", ".pyd")) and not (s_n.startswith(env) or s_n.startswith(sysroot)):
            repl = next((os.path.join(d, os.path.basename(src)) for d in env_dirs
                         if os.path.isfile(os.path.join(d, os.path.basename(src)))), None)
            print(f"volva.spec: FOREIGN DLL {src} -> {repl or 'kept (no replacement found)'}")
            if repl:
                src = repl
        fixed.append((dest, src, kind))
    a.binaries = fixed

pyz = PYZ(a.pure)
exe = EXE(pyz, a.scripts, [], exclude_binaries=True, name="Volva", console=False,
          icon="volva/gui/volva.ico")
coll = COLLECT(exe, a.binaries, a.datas, name="Volva")
