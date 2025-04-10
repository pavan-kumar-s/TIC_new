pth = 'D:\OneDrive - smail.iitm.ac.in\SprintCore\TIC\deleteit_2\TIC_new\minTIC\models\';
models = {'STM_v1_0','e_coli_core','iAB_RBC_283','iAF1260','iAF1260b','iAF692','iAF987'};
for m=2:numel(models)
    load([pth,models{m}]);
    tic
    [minFlux, maxFlux] = fluxVariability(model,'allowLoops', 0);
    fva=toc
    a1 = getfeasFluxDir(minFlux,maxFlux);
    tic
    [a,modModel] = ThermOptCC(model,1e-4);
    tocc=toc
end

function a1 = getfeasFluxDir(minF,maxF)
a1={};
for i =1:numel(minF)
    mi = minF(i); ma = maxF(i);
    if mi<0 && ma>0
        a1{i} = 'Reversible';
    elseif mi>=0 && ma>0
        a1{i} = 'Forward';
    elseif mi<0 && ma<=0
        a1{i} = 'Reverse';
    elseif mi==0 && ma==0
        a1{i} = 'Blocked';
    end
end
a1=a1';
end
  
