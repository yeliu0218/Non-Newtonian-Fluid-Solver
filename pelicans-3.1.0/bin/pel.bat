
@rem ----- ExeScript Options Begin -----
@rem ScriptType: console
@rem DestDirectory: temp
@rem Icon: default
@rem OutputFile: C:\WINPELICANS\bin\pel.exe
@rem ----- ExeScript Options End -----
@rem = '--*-Perl-*--
@echo off
if "%OS%" == "Windows_NT" goto WinNT
perl -I %PELICANSHOME%\tools\pel -S "%PELICANSHOME%\tools\pel\%0.pl" %1 %2 %3 %4 %5 %6 %7 %8 %9
goto endofperl
:WinNT
perl -I %PELICANSHOME%\tools\pel -S %PELICANSHOME%\tools\pel\%0.pl %*
if NOT "%COMSPEC%" == "%SystemRoot%\system32\cmd.exe" goto endofperl
if %errorlevel% == 9009 echo You do not have Perl in your PATH.
if errorlevel 1 goto script_failed_so_exit_with_non_zero_val 2>nul
goto endofperl
:endofperl
