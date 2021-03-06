@ECHO OFF

SET VSREDIST=%1%

REM ############################################################################
REM # Copy master __init__.py file, plus LICENCE, CHANGES, README, and INSTALL #
REM ############################################################################

IF NOT EXIST "phycas" (
    ECHO %0 expected to find a directory named "phycas" inside %CWD%
    EXIT /B 1)

CD phycas
ECHO Copying files to %cd%...
COPY ..\src\python\__init__.py .
COPY ..\README .
COPY ..\INSTALL .
COPY ..\CHANGES .
COPY ..\LICENSE .

REM ###############################################
REM # Copy python files for conversions extension #
REM ###############################################

IF NOT EXIST "conversions" (
    ECHO %0 expected to find a directory named "conversions" inside %CWD%
    EXIT /B 1)

CD conversions
ECHO Copying files to %cd%...
DEL /Q *.py libncl.lib
COPY ..\..\src\python\conversions\*.py .
COPY %VSREDIST%\msvcp120.dll .
COPY %VSREDIST%\msvcr120.dll .
COPY %VSREDIST%\vccorlib120.dll .
CD ..

REM ###############################################
REM # Copy Nexus Class Library DLL                #
REM ###############################################

ECHO Copying ncl.dll to %cd%\conversions...
COPY %NCL_INSTALL_DIR%\bin\msvc-12.0\release\threading-multi\ncl.dll conversions

REM ##############################################
REM # Copy python files for datamatrix extension #
REM ##############################################

IF NOT EXIST "datamatrix" (
    ECHO %0 expected to find a directory named "datamatrix" inside %CWD%
    EXIT /B 1)

CD datamatrix
ECHO Copying files to %cd%...
DEL /Q *.py
COPY ..\..\src\python\datamatrix\*.py .
CD ..

REM ##############################################
REM # Copy python files for likelihood extension #
REM ##############################################

IF NOT EXIST "likelihood" (
    ECHO %0 expected to find a directory named "likelihood" inside %CWD%
    EXIT /B 1)

CD likelihood
ECHO Copying files to %cd%...
DEL /Q *.py
COPY ..\..\src\python\likelihood\*.py .
CD ..

REM ##############################################
REM # Copy python files for phylogeny extension #
REM ##############################################

IF NOT EXIST "phylogeny" (
    ECHO %0 expected to find a directory named "phylogeny" inside %CWD%
    EXIT /B 1)

CD phylogeny
ECHO Copying files to %cd%...
DEL /Q *.py
COPY ..\..\src\python\phylogeny\*.py .
CD ..

REM ############################################
REM # Copy python files for probdist extension #
REM ############################################

IF NOT EXIST "probdist" (
    ECHO %0 expected to find a directory named "probdist" inside %CWD%
    EXIT /B 1)

CD probdist
ECHO Copying files to %cd%...
DEL /Q *.py
COPY ..\..\src\python\probdist\*.py .
CD ..

REM #############################################
REM # Copy python files for readnexus extension #
REM #############################################

IF NOT EXIST "readnexus" (
    ECHO %0 expected to find a directory named "readnexus" inside %CWD%
    EXIT /B 1)

CD readnexus
ECHO Copying files to %cd%...
DEL /Q *.py
COPY ..\..\src\python\readnexus\*.py .
CD ..

REM #######################################
REM # Copy python files for pdfgen module #
REM #######################################

IF EXIST "pdfgen" (
    RMDIR /S /Q pdfgen)

MKDIR pdfgen
IF NOT ERRORLEVEL 0 (
	ECHO %0 could not create pdfgen directory inside %CWD%
    EXIT /B 1)

CD pdfgen
ECHO Copying files to %cd%...
COPY ..\..\src\python\pdfgen\*.py .
XCOPY /I /E /S ..\..\src\python\pdfgen\AFM AFM
CD ..

REM ###########################################
REM # Copy python files for treeviewer module #
REM ###########################################

IF EXIST "treeviewer" (
    RMDIR /S /Q treeviewer)

MKDIR treeviewer
IF NOT ERRORLEVEL 0 (
	ECHO %0 could not create treeviewer directory inside %CWD%
    EXIT /B 1)

CD treeviewer
ECHO Copying files to %cd%...
COPY ..\..\src\python\treeviewer\*.py .
CD ..

REM ##########################################
REM # Copy python files for utilities module #
REM ##########################################

IF EXIST "utilities" (
    RMDIR /S /Q utilities)

MKDIR utilities
IF NOT ERRORLEVEL 0 (
	ECHO %0 could not create utilities directory inside %CWD%
    EXIT /B 1)

CD utilities
ECHO Copying files to %cd%...
COPY ..\..\src\python\utilities\*.py .
CD ..

REM #########################################
REM # Copy python files for commands module #
REM #########################################

IF EXIST "commands" (
    RMDIR /S /Q commands)

MKDIR commands
IF NOT ERRORLEVEL 0 (
	ECHO %0 could not create commands directory inside %CWD%
    EXIT /B 1)

CD commands
ECHO Copying files to %cd%...
COPY ..\..\src\python\commands\*.py .
CD ..

REM ##############
REM # Copy tests #
REM ##############

IF EXIST "tests" (
    RMDIR /S /Q tests)

MKDIR tests
IF NOT ERRORLEVEL 0 (
	ECHO %0 could not create tests directory inside %CWD%
    EXIT /B 1)

CD tests
ECHO Copying files to %cd%...
XCOPY /E /S ..\..\tests\* .
CD ..

REM #################
REM # Copy examples #
REM #################

IF EXIST "examples" (
    RMDIR /S /Q examples)

MKDIR examples
IF NOT ERRORLEVEL 0 (
	ECHO %0 could not create examples directory inside %CWD%
    EXIT /B 1)

CD examples
ECHO Copying files to %cd%...
XCOPY /E /S ..\..\examples\* .
CD ..

REM ####################################
REM # Navigate out of phycas directory #
REM ####################################
CD ..
ECHO Navigating to original directory %cd%...
