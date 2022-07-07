call "C:\Program Files (x86)\Microsoft Visual Studio\2019\Professional\VC\Auxiliary\Build\vcvars64.bat"

cls

set "FLAG=/c /EHsc /MP /GS- /W0 /Gy- /I"." /Zi /O2 /fp:strict /D "ENABLE_RELATEDNESS" /D "_CONSOLE" /D "_UNICODE" /D "UNICODE" /fp:except- /GF /GT /WX- /std:c17 /GR- /arch:SSE2 /Gd /Oy /Oi /MT /std:c++17  /EHsc /nologo  /Ot /openmp -msse3 -msse4.1 -msse4.2 -mlzcnt -mpopcnt -Wwritable-strings"

clang-cl.exe %FLAG% dre.cpp
clang-cl.exe %FLAG% ml.cpp
clang-cl.exe %FLAG% ml2.cpp
clang-cl.exe %FLAG% mlbin.cpp
clang-cl.exe %FLAG% mom.cpp
clang-cl.exe %FLAG% mom2.cpp
clang-cl.exe %FLAG% global.cpp
clang-cl.exe %FLAG% hash.cpp
clang-cl.exe %FLAG% math2.cpp
clang-cl.exe %FLAG% /arch:SSE2 mathSSE.cpp
clang-cl.exe %FLAG% /arch:AVX2 mathAVX.cpp
clang-cl.exe %FLAG% /arch:AVX512 math512.cpp
clang-cl.exe %FLAG% misc.cpp
clang-cl.exe %FLAG% file.cpp
clang-cl.exe %FLAG% string2.cpp
clang-cl.exe %FLAG% statistics.cpp
clang-cl.exe %FLAG% matrix.cpp
clang-cl.exe %FLAG% parameters.cpp
clang-cl.exe %FLAG% function.cpp
clang-cl.exe %FLAG% menu.cpp
clang-cl.exe %FLAG% vcfpop.cpp

clang-cl.exe -MT -DSFML_STATIC vcfpop.obj global.obj misc.obj file.obj hash.obj math2.obj math512.obj mathAVX.obj mathSSE.obj string2.obj statistics.obj matrix.obj parameters.obj function.obj menu.obj ml.obj ml2.obj mlbin.obj mom.obj mom2.obj dre.obj    zlibstat.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib libomp.lib

copy vcfpop.exe ..\bin\
del *.obj
del vcfpop.exe
pause