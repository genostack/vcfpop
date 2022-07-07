call "C:\Program Files (x86)\Microsoft Visual Studio\2019\Professional\VC\Auxiliary\Build\vcvars64.bat"

cls
set "FLAG=/c /EHsc /MP /GS- /TP /Qpar /GL /Gy- /Zc:wchar_t /I"." /Zi /O2 /Ob2 /Zc:inline /fp:strict /D "NDEBUG" /D "_CONSOLE" /D "_UNICODE" /D "UNICODE" /fp:except- /GF /w /W0 /WX- /GR- /Gd /Oy /Oi /MT /openmp /std:c++17 /FC /nologo /Ot"

cl.exe %FLAG% dre.cpp
cl.exe %FLAG% ml.cpp
cl.exe %FLAG% ml2.cpp
cl.exe %FLAG% mlbin.cpp
cl.exe %FLAG% mom.cpp
cl.exe %FLAG% mom2.cpp
cl.exe %FLAG% global.cpp
cl.exe %FLAG% hash.cpp
cl.exe %FLAG% math2.cpp
cl.exe %FLAG% mathSSE.cpp
cl.exe %FLAG% mathAVX.cpp
cl.exe %FLAG% math512.cpp
cl.exe %FLAG% misc.cpp
cl.exe %FLAG% file.cpp
cl.exe %FLAG% string2.cpp
cl.exe %FLAG% statistics.cpp
cl.exe %FLAG% matrix.cpp
cl.exe %FLAG% parameters.cpp
cl.exe %FLAG% function.cpp
cl.exe %FLAG% menu.cpp
cl.exe %FLAG% vcfpop.cpp

cl.exe -MT -DSFML_STATIC /OUT:vcfpop.exe /LTCG:incremental /DYNAMICBASE /NDEBUG /MACHINE:X64 /OPT:REF /INCREMENTAL:NO /OPT:ICF /ERRORREPORT:PROMPT /TLBID:1 vcfpop.obj global.obj misc.obj file.obj hash.obj math2.obj math512.obj mathAVX.obj mathSSE.obj string2.obj statistics.obj matrix.obj parameters.obj function.obj menu.obj ml.obj ml2.obj mlbin.obj mom.obj mom2.obj dre.obj zlibstat.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib libomp.lib

copy vcfpop.exe ..\bin\
del *.obj
del vcfpop.exe
pause