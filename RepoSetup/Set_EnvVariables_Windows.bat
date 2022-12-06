:: TODO: This script does not check whether the desired paths have already been added to the env variables
@echo off
setlocal enabledelayedexpansion
set SCRIPT_DIR=%~dp0%
echo %SCRIPT_DIR%
set P3_REPO=%SCRIPT_DIR%..
set MODULE_REL_DIRS=..\BdPyTools ..\RunningSimulations\RFTrack
(for %%m in (%MODULE_REL_DIRS%) do (
    echo %%m
    set dirToAdd=%SCRIPT_DIR%%%m
    echo !dirToAdd!
    set PYTHONPATH=!dirToAdd!;!PYTHONPATH!
    set JUPYTER_PATH=!dirToAdd!;!JUPYTER_PATH!
))
echo Going to set env variable P3_REPO to:
echo %P3_REPO%
setx P3_REPO %P3_REPO%
echo Going to set env variable PYTHONPATH to:
echo %PYTHONPATH%
setx PYTHONPATH %PYTHONPATH%
echo Going to set env variable JUPYTER_PATH to:
echo %JUPYTER_PATH%
setx JUPYTER_PATH %JUPYTER_PATH%
echo ... everything done^^!
pause
