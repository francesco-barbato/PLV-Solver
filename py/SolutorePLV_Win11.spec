# -*- mode: python ; coding: utf-8 -*-
# Build Windows 11 (.exe onedir) con PyInstaller:
#   cd py && pyinstaller --noconfirm SolutorePLV_Win11.spec

a = Analysis(
    ['ui_strutture.py'],
    pathex=['.'],
    binaries=[],
    datas=[
        ('solutore_plv.py', '.'),
        ('icona.ico', '.'),
        ('icona.jpg', '.'),
        ('../JSON', 'JSON'),
    ],
    hiddenimports=[
        'matplotlib.backends.backend_qtagg',
        'matplotlib.backends.backend_qt',
        'PyQt6.QtPrintSupport',
        'sympy',
        'numpy',
        'typing_extensions',
    ],
    hookspath=[],
    hooksconfig={
        'matplotlib': {
            'backends': ['Agg', 'QtAgg'],
        },
    },
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='SolutorePLV',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    icon=['icona.ico'],
)
coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=False,
    upx_exclude=[],
    name='SolutorePLV',
)
