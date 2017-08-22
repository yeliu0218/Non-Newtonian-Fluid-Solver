@rem ##############################################################
@rem equivalent to "make all" with the unix administration makefile
@rem ##############################################################
@if defined VS90COMNTOOLS goto setEnv
@echo error: environment Visual Studio .NET 2008 not installed
@set retcode=1
@goto quit
:setEnv
@rem ### environment setting for VC++
@call "%VS90COMNTOOLS%\vsvars32.bat"
@rem environnement setting for PELICANS
@set HERE=%CD%
@cd ..\..\..
@set PELICANSHOME=%CD%
@cd %HERE%
:parseCommandLine
@rem Build libpel0 and libpel1
@if "%1" EQU "all" (
	call :build pelicans lib0
	call :build pelicans lib1
	goto testerr )
@rem Build libpel0
@if "%1" EQU "lib0" (
	call :build pelicans lib0
	goto testerr )
@rem Build libpel1
@if "%1" EQU "lib1" (	
	call :build pelicans lib1
	goto testerr )
@rem Build libpel2
@if "%1" EQU "lib2" (
	call :build pelicans lib2
	goto testerr )
@rem Build libpel1 and exe1	
@if "%1" EQU "exe1" (
	call :build pelicans lib1
	call :build examplesofapplication "exe1 (lib1)|Win32"
	goto testerr )
@rem Build libpel2 and exe2	
@if "%1" EQU "exe2" (
	call :build pelicans lib2
	call :build examplesofapplication "exe2 (lib2)|Win32"
	goto testerr )	
@rem Build libpel2 , exe2 then run tests	
@if "%1" EQU "check" (
	call :build pelicans libg
	call :build examplesofapplication "exeg (libg)|Win32"
	if not exist %PELICANSHOME%\tests\windows-vc2008 mkdir %PELICANSHOME%\tests\windows-vc2008
	cd %PELICANSHOME%\tests\windows-vc2008
	set PATH=%PELICANSHOME%\lib\windows-vc2008;%PATH%
	echo "Running unit-tests: it can take several minutes depending on machine"...
	%PELICANSHOME%\tests\lib\windows-vc2008\exeg %PELICANSHOME%\admin\check.pel -Call > resu
	type resu
	cd %HERE%
	goto testerr )
@if "%1" EQU "clean" (
	if exist %PELICANSHOME%\tests\lib\windows-vc2008 del /S /Q %PELICANSHOME%\tests\lib\windows-vc2008
	if exist %PELICANSHOME%\lib\windows-vc2008 del /S /Q %PELICANSHOME%\lib\windows-vc2008
	goto testerr )
@echo Usage: %0 TARGET
@echo With TARGET in:
@echo     all         : build libraries lib0 and lib1
@echo     lib0        : build libraries with internal checking level set to 0
@echo     lib1        : build libraries with internal checking level set to 1
@echo     lib2        : build libraries with internal checking level set to 2
@echo     exe1        : build example of Applications executable
@echo                   with internal checking level set to 1 and linked to lib1
@echo     exe2        : build example of Applications executable
@echo                   with internal checking level set to 2 and linked to lib2
@echo     exeg        : build example of Applications executable with debugging informations
@echo                   and internal checking level set to 2 and linked to libg
@echo     check       : check install of built-in functionalities with exeg
@echo     clean       : delete all generated files

@set retcode=1
@goto quit
:testerr
@if not errorlevel 1 goto normalEnd
@echo "Error in compilation"
@if "%2" EQU "shutdown" (
@shutdown /l /f
)
@set retcode=1
@goto quit

:normalEnd
@echo "That's all, folks!"
@set retcode=0

:quit
@if "%2" EQU "shutdown" (
@shutdown /l /f
)
@exit /b %retcode%

:build
@msbuild %1.vcproj /verbosity:d /t:Build /p:Configuration=%2
@if not errorlevel 1 exit /B 0
@exit /B 1
