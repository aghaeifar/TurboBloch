clc

build_dir   = fileparts(mfilename('fullpath'));
parent_dir  = fileparts(build_dir);
src1        = [' ' '"' fullfile(build_dir, 'bloch_mex.cpp') '"' ];   
src2        = [' ' '"' fullfile(parent_dir, 'bloch.cpp') '"']; 
out_path    = fullfile(parent_dir, 'lib');
out_path    = build_dir;

if isunix
    flags = ' CXXFLAGS="\$CXXFLAGS -std=c++17 -D__MEASURE_ELAPSED_TIME__ -D__SINGLE_PRECISION__"';
    ink_dir =  ''; % [' -I' '"' '/usr/include/mkl/' '"'];
    % libs = ' /usr/lib/x86_64-linux-gnu/libmkl_intel_lp64.a /usr/lib/x86_64-linux-gnu/libmkl_intel_thread.a /usr/lib/x86_64-linux-gnu/libmkl_core.a -liomp5 -lpthread -lm -ldl -ltbb';
    libs = ' -lm -ldl -ltbb';
    eval(['mex -v' flags ink_dir libs src1 src2 ' -R2018a' ' -outdir ' out_path]); 
elseif ispc
    ink_dir = ''; %[' -I' '"' fullfile(getenv('MKLROOT'), 'include') '"'];
    lnk_dir1 = ''; % [' -L' '"' fullfile(getenv('MKLROOT'), 'lib', 'intel64') '"'];   
    lnk_dir2 = ''; %[' -L' '"' fullfile(getenv('ONEAPI_ROOT'), 'compiler', 'latest', 'windows', 'compiler', 'lib',  'intel64_win') '"']; 
    libs = ''; %' -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -llibiomp5md';
    flags = ' COMPFLAGS="$COMPFLAGS /std:c++17 /D__MEASURE_ELAPSED_TIME__ /D__SINGLE_PRECISION__"';
    eval(['mex ' flags ink_dir lnk_dir1 lnk_dir2 libs src1 src2 ' -R2018a' ' -outdir ' out_path]); 
else
    error('Platform not supported')
end

