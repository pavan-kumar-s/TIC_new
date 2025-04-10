function [a,modModel] = ThermOptCC(model,tol)
% Identifies thermodynamically feasible flux directions for all the
% reactions in the input model
%
% USAGE: 
%   a = ThermOptCC(model,tol)
%
% INPUTS:
%     model:     COBRA model structure for which thermodynamic feasibility
%                of the reactions has to be identified
%     tol:       Tolerance value (User defined non-zero value).
%
% OUTPUTS:
%     a:          A cell describing the thermodynamically feasible direction 
%                 of the reactions in the given input model
%     modModel:   Modified model that has no irreversible reactions that
%                 carry flux in reverse direction. The feasible directions
%                 are reported for this model.
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras


% checking for reactions irreversible in reverse direction
IrR = model.ub<=0; % reactions irreversible in reverse direction
[~,n] = size(model.S);
temp = model.ub(IrR);
model.S(:,IrR) = -model.S(:,IrR);
model.ub(IrR) = -1*model.lb(IrR);
model.lb(IrR) = -1*temp;
rev = ones(n,1);
rev(model.lb>=0) = 0;
model.rev=rev;
modModel = model;
% normalising the bounds to max of 10000
temp1 = max(abs(model.lb));
temp2 = max(abs(model.ub));
model.lb(abs(model.lb)==temp1) = sign(model.lb(abs(model.lb)==temp1))*10000;
model.ub(abs(model.ub)==temp2) = sign(model.ub(abs(model.ub)==temp2))*10000;

% the model should not have any reaction with flux only in reverse
% direction
[modelIrrev, ~, rev2irrev] = convertToIrreversible(model);
[~,n] = size(modelIrrev.S);
[nullity,a] = getNullSpaceAndTICrxns(modelIrrev,tol);

core = 1:n; newCore = [];
while numel(core)~=numel(newCore)
    newCore = core;
    IDS = findConsistentIDS(modelIrrev,core,nullity,a,tol);
    core = setdiff(core,IDS);
end
consIrr = double(~ismember([1:n],core)); % consistent irreversible reactions
a = cellfun(@(x)getConsDir(consIrr,x),rev2irrev,'UniformOutput',0);
end


function consDir = getConsDir(consIrr,x)
temp = consIrr(x);
if sum(temp)==2
    consDir = 'Reversible'; % consistent in both the direction
elseif sum(temp)==0
    consDir = 'Blocked'; % blocked reaction
elseif temp(1)==1
    consDir = 'Forward'; % consistent only in forward direction
elseif temp(2)==1
    consDir = 'Reverse'; % consistent only in reverse direction
end
end

function reacInd = findConsistentIDS(model,core,nullity,a,tol)
[m,n] = size(model.S);
n_TICR = numel(a);
K=1000;
% order of decision variables: v (\in R^n), z (\in R^n), mu (\in
% R^n_TICR), eps (\in B^n_TICR)

% objective
f = [zeros(n,1);double(ismember([1:n]',core));zeros(2*n_TICR,1)];

% equalities
Aeq1 = [model.S, sparse(m,n+2*n_TICR)];
beq1 = zeros(m,1);
csenseeq1 = repmat('E',m,1); % equality

% inequalities
temp1 = -1*speye(n);
temp2 = speye(n);
Aineq1 = [temp1,temp2,sparse(n,2*n_TICR)];
bineq1 = zeros(n,1);
csenseineq1 = repmat('L',n,1); % lesser than

% equalities
k = size(nullity,2);
Aeq2 = [sparse(k,2*n),nullity',sparse(k,n_TICR)];
beq2 = zeros(k,1);
csenseeq2 = repmat('E',k,1); % equality

% inequalities
Aineq2 = [sparse(n_TICR,2*n),speye(n_TICR),(K+1)*speye(n_TICR)];
bineq2 = ones(n_TICR,1);
csenseineq2 = repmat('G',n_TICR,1); % lesser than

% inequalities
Aineq3 = [sparse(n_TICR,2*n),speye(n_TICR),(K+1)*speye(n_TICR)];
bineq3 = K*ones(n_TICR,1);
csenseineq3 = repmat('L',n_TICR,1); % lesser than

% inequalities
temp = -1*speye(n);
temp = temp(a,:);
Aineq4 = [temp,sparse(n_TICR,n+n_TICR),K*speye(n_TICR)];
bineq4 = K*ones(n_TICR,1);
csenseineq4 = repmat('L',n_TICR,1); % lesser than

% inequalities
temp = -1*speye(n);
temp = temp(a,:);
Aineq5 = [temp,sparse(n_TICR,n+n_TICR),K*speye(n_TICR)];
bineq5 = zeros(n_TICR,1);
csenseineq5 = repmat('G',n_TICR,1); % lesser than


% bounds
lb = [model.lb;zeros(n,1);-K*ones(n_TICR,1);zeros(n_TICR,1)];
ub = [model.ub;tol*ones(n,1);K*ones(n_TICR,1);ones(n_TICR,1)];

% Set up MILP problem
MILPproblem.A=[Aeq1;Aeq2;Aineq1;Aineq2;Aineq3;Aineq4;Aineq5];
MILPproblem.b=[beq1;beq2;bineq1;bineq2;bineq3;bineq4;bineq5];
MILPproblem.lb=lb;
MILPproblem.ub=ub;
MILPproblem.c=f;
MILPproblem.osense=-1;%maximise
MILPproblem.vartype = [repmat('C',2*n+n_TICR,1);repmat('B',n_TICR,1)];
MILPproblem.csense = [csenseeq1; csenseeq2; csenseineq1; csenseineq2; csenseineq3; csenseineq4; csenseineq5];
solution = solveCobraMILP(MILPproblem);
stat=solution.stat;
if stat~=1
    x = [];
    reacInd = [];
else
    x=solution.full;
    reacInd = intersect(find(x(1:n)>=0.99*tol),core);
end
end


function [nullity,a] = getNullSpaceAndTICrxns(model,vTol)
model_temp = model;
ind=findExcRxns(model_temp);

% blocking all the exchange reactions
model_temp.lb(ind) = 0; model_temp.ub(ind) = 0;
model_temp.lb(model_temp.lb>0)=0;
if exist('sprintcc','file')
    a = sprintcc(model_temp,1); % the reaction that are involved in TICs (or reversible ractions)
else
    a = fastcc(model_temp,1); % the reaction that are involved in TICs (or reversible ractions)
end

% model with only TICs
model_temp = removeRxns(model,model.rxns(setdiff([1:numel(model.rxns)],a)));
[~,a] = ismember(model_temp.rxns,model.rxns);
% nullity = null(model_temp.S);


nullity = fast_snp(model_temp.S,model_temp.lb,model_temp.ub,vTol);
end