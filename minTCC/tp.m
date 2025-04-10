pth = 'D:\OneDrive - smail.iitm.ac.in\SprintCore\TIC\deleteit_2\TIC_new\minTIC\models\';
load([pth,'e_coli_core.mat'])
[a,modModel] = ThermOptCC(model,1e-4);
