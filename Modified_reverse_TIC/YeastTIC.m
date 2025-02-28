load('./yeast9.mat')
[TICs,Direction,TIC_Rxns,modModel,opt] = ThermOptEnumMILP_rev_yeast(model);
