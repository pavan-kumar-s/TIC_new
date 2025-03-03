clear
load('D:\OneDrive_moved_files\Microbial_community_gap_filling\With_AGORA2\AGORA2\Brevibacillus_brevis_FJAT_0809_GLX.mat')
tic
[TICs1,Direction1] = ThermOptEnumMILP_feas(model);
toc
tic
[TICs,Direction] = ThermOptEnumMILP_rev(model);
toc

[tics1,id1] = sort(cellfun(@(x)strjoin(sort(x)),TICs1,'UniformOutput',false));
[tics,id] = sort(cellfun(@(x)strjoin(sort(x)),TICs,'UniformOutput',false));