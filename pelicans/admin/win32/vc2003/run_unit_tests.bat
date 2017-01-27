call init.bat
set OUTPUT=%PELICANSHOME%\tests\windows-VC2003
set EXE=%PELICANSHOME%\tests\lib\windows-VC2003\opt1\pelicans.exe
echo %OUTPUT%
mkdir %OUTPUT%
cd %OUTPUT%
echo PELICANSHOME : %PELICANSHOME% > results.txt
echo current date %date% - %time% >> results.txt
%EXE% -v ..\..\admin\win32\non_regression_tests.pel >> results.txt
cd %here%
