%% Mex files creation minGW compiler needed (matlab add-on)!!!

disp('Compilation started (might take some times if it is first time runing minGW)')

mex binarySearch/binarySearch.c; % Compile binary search
run CONVNFFT_Folder/CONVNFFT_Folder/convnfft_install.m; % Compile inplaceprod

if ispc
    run matlab-dtts-1.1/compileDttMex.m; % Compile fast sine transform mex
else % Compiled fast discrete sine transfrom for 2D, not needed for 4D, for mac it has to be compiled manually, see warning.
    warning('2D SFT: For windows, and mac with intel, given files should work (just try to run main2D.m), if recompilation is needed (main2D.m not working), it has to be done manually, please look into matlab-dtts-1.1/compileDttMex.m')
end

disp('Now you can run main2D.m or main4D.m')
