@echo off
setlocal
chcp 65001 >nul
cd /d "%~dp0"

set "PYCMD="
set "PYWCMD="

where py >nul 2>nul
if errorlevel 1 goto :check_python
py -3 -c "import sys; raise SystemExit(0 if sys.version_info.major > 3 or sys.version_info.minor >= 9 else 1)" >nul 2>nul
if errorlevel 1 goto :check_python
set "PYCMD=py -3"
where pyw >nul 2>nul
if not errorlevel 1 set "PYWCMD=pyw -3"
goto :check_deps

:check_python
where python >nul 2>nul
if errorlevel 1 goto :no_python
python -c "import sys; raise SystemExit(0 if sys.version_info.major > 3 or sys.version_info.minor >= 9 else 1)" >nul 2>nul
if errorlevel 1 goto :no_python
set "PYCMD=python"
where pythonw >nul 2>nul
if not errorlevel 1 set "PYWCMD=pythonw"
goto :check_deps

:no_python
echo [ERRORE] Python 3.9 o successivo non trovato su Windows 11.
echo Scaricalo da https://www.python.org/downloads/
echo e spunta "Add python.exe to PATH" prima di installare.
pause
exit /b 1

:check_deps
%PYCMD% -c "import PyQt6, PyQt6.QtPrintSupport, matplotlib, numpy, sympy, typing_extensions" >nul 2>nul
if not errorlevel 1 goto :launch

echo Prima configurazione: installazione librerie necessarie in corso...
%PYCMD% -m pip install PyQt6 matplotlib numpy sympy typing-extensions
if errorlevel 1 (
    echo [ERRORE] Impossibile installare le librerie. Verifica la connessione Internet.
    pause
    exit /b 1
)

:launch
if "%PYWCMD%"=="" (
    %PYCMD% "%~dp0py\ui_strutture.py"
    if errorlevel 1 pause
) else (
    start "" %PYWCMD% "%~dp0py\ui_strutture.py"
)
