@echo off

rem
rem  Macaulay2 wrapper by Martin Bujnak
rem  (GBsolver subroutine)


c:
cd "!!TODO:c:\Program Files\Macaulay2\bin\"

rem  The slashes in the next line should be /
set M2HOME=!!TODO:c:/Progra~1/Macaulay2/lib/Macaulay2-0.9.2

rem  The slashes in the next line should be \
set M2WINH=!!TODO:c:\Progra~1\Macaulay2\lib\Macaulay2-0.9.2

rem to help find cygwin.dll
set PATH=%PATH%;%M2HOME%\..\..\bin

rem  Don't change this line.
"%M2WINH%\libexec\Macaulay2" -tty -ephase=1 "%M2HOME%/m2/setup.m2" -ephase=0 '-e runStartFunctions()' '-silent' '-e load "!!TODO:FULL PATH TO CODE, e.g. C:/../GBpolysolver/gbsMacaulay/code.m2"'


