function [TICs,Direction,TIC_Rxns,modModel,opt] = ThermOptEnumMILP_minimal(model,timeLimit)
% Enumerates all the Thermodynamically infeasible cycles in a given model
%
% USAGE: 
%   [TICs,Direction,TIC_Rxns,modModel,opt] = ThermoOptEnumMILP(model,timeLimit)
%
% INPUTS:
%     model:     COBRA model structure for which TICs has be found
%
% OPTIONAL INPUTS:
%     timeLimit: If the algorithm takes more than this time null set will
%                be returned for all the outputs
%
% OUTPUTS:
%     TICs:       List of all the Thermodynamically infeasible cycles in
%                 the given input model
%     Direction:  Relative flux coefficients for reactions 
%                 in the corresponding TICs
%     TIC_Rxns:   Reaction list that participates in the TICs
%     modModel:   Modified model that has no irreversible reactions that
%                 carry flux in reverse direction
%     opt:        Says whether the provided solution is optimal or not
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras

% checking for reactions irreversible in reverse direction
changeCobraSolverParams('LP','feasTol',1e-9);
[~,n] = size(model.S); % number of metabolites and reactions
IrR = model.ub<=0; % reactions irreversible in reverse direction
temp = model.ub(IrR);
model.S(:,IrR) = -model.S(:,IrR);
model.ub(IrR) = -1*model.lb(IrR);
model.lb(IrR) = -1*temp;
rev = ones(n,1);
rev(model.lb>=0) = 0;
model.rev=rev;
% converting all the positive lower bounds to zero lower bounds
model.lb(model.lb>0)=0;

% normalising the bounds to max of 1000
temp1 = max(abs(model.lb));
temp2 = max(abs(model.ub));
model.lb(abs(model.lb)==temp1) = sign(model.lb(abs(model.lb)==temp1))*1000;
model.ub(abs(model.ub)==temp2) = sign(model.ub(abs(model.ub)==temp2))*1000;

modModel = model;
tol=1;
ind=findExcRxns(model);
% blocking all the exchange reactions
model.lb(ind) = 0; model.ub(ind) = 0;

if exist('sprintcc','file')
    a = sprintcc(model,tol); 
else
    a = fastcc(model,tol);
end

if isempty(a)
    TICs={};Direction={};TIC_Rxns={};opt=1;
    return
end
TIC_Rxns = model.rxns(a);
% model with only TICs
model = removeRxns(model,model.rxns(setdiff([1:numel(model.rxns)],a)));
model_temp =model;
TICs={};Direction={};
order = getOrderOfRxns(model);
RXNS = model.rxns(order);

if exist('timeLimit', 'var') && ~isempty(timeLimit)
     tic
else
    opt=1;
end
TICmat=[];
TICcons={}; % TIC constraints
while ~isempty(model.rxns)
    r = RXNS(1);
    core = ismember(model.rxns,r);
    if exist('timeLimit', 'var') && ~isempty(timeLimit)
        if toc>timeLimit
            opt=0;
            TICs=[];Direction=[];TIC_Rxns=[];modModel=[];
            warning('TimeLimitReached')
            return 
        else
            opt=1;
        end
    end
    
    while 1
        if exist('timeLimit', 'var') && ~isempty(timeLimit)
            if toc>timeLimit
                opt=0;
                TICs=[];Direction=[];TIC_Rxns=[];modModel=[];
                warning('TimeLimitReached')
                return  
            else
                opt=1;
            end
        end
        
        dir=1;
        [~,indices] = ismember(model.rxns,model_temp.rxns);
        [TICmat,TICcons] = getTICmat_TICcons(TICs,Direction,model_temp.rxns,model.rxns);
        [m1,blkdCore,flux] = getTICModel2(model,core,tol,dir,TICcons,TICmat,indices);
        if blkdCore
            break
        end
        TICs{end+1,1} = m1; 
        ids =find(ismember(model.rxns,m1));
        Direction{end+1,1} =flux(ids)/min(abs(flux(ids)));
        
        
    end
    
    while 1
        if exist('timeLimit', 'var') && ~isempty(timeLimit)
            if toc>timeLimit
                opt=0;
                TICs=[];Direction=[];TIC_Rxns=[];modModel=[];
                warning('TimeLimitReached')
                return 
            else
                opt=1;
            end
        end
        dir=-1;
        [~,indices] = ismember(model.rxns,model_temp.rxns);
        [TICmat,TICcons] = getTICmat_TICcons(TICs,Direction,model_temp.rxns,model.rxns);
        [m1,blkdCore,flux] = getTICModel2(model,core,tol,dir,TICcons,TICmat,indices);
        if blkdCore
            break
        end
        TICs{end+1,1} = m1; 
        ids =find(ismember(model.rxns,m1));
        Direction{end+1,1} =flux(ids)/min(abs(flux(ids)));
        
    end
    model.lb(core) = 0;model.ub(core) = 0;
    if exist('sprintcc','file')
        a = sprintcc(model,tol); 
    else
        a = fastcc(model,tol);
    end

    rmRXNS = model.rxns(setdiff([1:numel(model.rxns)],a));
    model = removeRxns(model,rmRXNS);
    RXNS = RXNS(~ismember(RXNS,rmRXNS));
end
end

function [TICmat,TICcons] = getTICmat_TICcons(TICs,Direction,rxns1,rxns2)
TICmat = zeros(numel(rxns1),numel(TICs));
TICcons={};
for t =1:numel(TICs)
    temp=TICs{t};
    temp2 = Direction{t};
    [loca,locb] = ismember(rxns1,temp);
    TICmat(loca,t)= temp2(locb(locb~=0));
    if numel(temp)==sum(ismember(rxns2,temp)) % checking if all the tics are present in the current model
        TICcons{end+1,1} = find(ismember(rxns2,temp));
    end
end
end