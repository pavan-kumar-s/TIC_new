load('./yeast9.mat')
[TICs,Direction,TIC_Rxns,modModel,opt] = ThermOptEnumMILP_feas_yeast(model);
