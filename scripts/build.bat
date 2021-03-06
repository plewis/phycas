@ECHO OFF

REM This script expects certain environmental variables (BOOST_ROOT, PYTHON_ROOT,
REM and NCL_INSTALL_DIR) to be defined before it is run. Examples of the necessary
REM definitions are provided below; modify appropriately for your installation.

REM The path to your boost installation (download from http://www.boost.org/)
SET BOOST_ROOT=%HOMEDRIVE%\boost_1_56_0

REM The path to your python installation (download from https://www.python.org/)
SET PYTHON_ROOT=%HOMEDRIVE%\Python27

REM The path to the Visual Studio dynamic link libraries  
SET VCREDIST="C:\Program Files (x86)\Microsoft Visual Studio 12.0\VC\redist\x86\Microsoft.VC120.CRT"

REM Modify PATH so that bjam executable can be found. This assumes bjam is just
REM inside $BOOST_ROOT, which will be the case if bootstrap.bat was run
SET PATH=%PATH%;%BOOST_ROOT%

REM This is needed for bjam to find its way
SET BOOST_BUILD_PATH=%HOMEDRIVE%\boost_1_56_0\tools\build\v2

REM Provide path to preinstalled Nexus Class Library
REM (download from http://sourceforge.net/projects/ncl/)
SET NCL_ALREADY_INSTALLED=
SET PHYCAS_NCL_STATIC=1
SET NCL_INSTALL_DIR=%HOMEDRIVE%\ncl-2.1.18

REM This removes the phycas directory created by previous build
RMDIR /S /Q phycas

REM This command initiates the build
echo Running bjam release...
REM SET BJAMDIR=%HOMEPATH%\Documents\libraries\bin
REM SET VCROOT="C:\Program Files (x86)\Microsoft Visual Studio 12.0\Common7\Tools\..\..\VC\"
%HOMEDRIVE%\boost-build-engine\bin\bjam.exe -q address-model=32 architecture=x86 release

CALL scripts\copy_src_python.bat %VCREDIST%

