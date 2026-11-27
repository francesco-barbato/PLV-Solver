@echo off
setlocal
chcp 65001 >nul
cd /d "%~dp0"

set "PYCMD="

where py >nul 2>nul
if errorlevel 1 goto :check_python
py -3 -c "import sys; raise SystemExit(0 if sys.version_info.major > 3 or sys.version_info.minor >= 9 else 1)" >nul 2>nul
if errorlevel 1 goto :check_python
set "PYCMD=py -3"
goto :build

:check_python
where python >nul 2>nul
if errorlevel 1 goto :no_python
python -c "import sys; raise SystemExit(0 if sys.version_info.major > 3 or sys.version_info.minor >= 9 else 1)" >nul 2>nul
if errorlevel 1 goto :no_python
set "PYCMD=python"
goto :build

:no_python
echo [ERRORE] Python 3.9+ non trovato nel PATH di Windows.
echo Scarica Python da https://www.python.org/downloads/
echo e seleziona "Add python.exe to PATH" durante l'installazione.
pause
exit /b 1

:build
echo Verifica pacchetti richiesti...
%PYCMD% -m pip install --upgrade pyinstaller PyQt6 matplotlib numpy sympy typing-extensions
if errorlevel 1 (
    echo [ERRORE] Installazione dipendenze non riuscita.
    pause
    exit /b 1
)

echo Compilazione SolutorePLV.exe per Windows 11 in corso...
%PYCMD% -m PyInstaller --noconfirm SolutorePLV_Win11.spec
if errorlevel 1 (
    echo [ERRORE] Compilazione PyInstaller fallita.
    pause
    exit /b 1
)

echo.
echo [OK] Eseguibile creato in: %~dp0dist\SolutorePLV\SolutorePLV.exe
pause
