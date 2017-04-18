if defined VS71COMNTOOLS goto build_and_smoke_tests
@echo error: environment Visual Studio .NET 2003 not installed
@exit /b 1
:build_and_smoke_tests
call "%VS71COMNTOOLS%\vsvars32.bat"
set VERSION=RELEASE
devenv pelicans.sln /build %VERSION%
call run_unit_tests.bat
