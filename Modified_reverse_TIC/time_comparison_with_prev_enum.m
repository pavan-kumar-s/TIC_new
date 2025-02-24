clear
load('D:\OneDrive_moved_files\Microbial_community_gap_filling\With_AGORA2\AGORA2\Brevibacillus_brevis_FJAT_0809_GLX.mat')

% tic
% [TICs,Direction] = ThermOptEnumMILP_rev(model);
% toc
% % Elapsed time is 680.226309 seconds.
% tic
% [TICs1,Direction1] = ThermOptEnumMILP(model);
% toc
% % Elapsed time is 679.545432 seconds.

tic
[TICs1,Direction1] = ThermOptEnumMILP(model);
toc
tic
[TICs,Direction] = ThermOptEnumMILP_rev(model);
toc
