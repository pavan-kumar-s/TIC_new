function [m1,blkdCore,flux,stat] = getTICModel2(model,core,tol,dir,TICcons,TICmat,indices)
[~,n] = size(model.S); % number of metabolites and reactions

direction = zeros(n,1);
if dir ==1
    % for forward direction
    direction(core)=1;
    [reacInd,x,stat] = findConsistentReacIDTIC(model,direction,tol,TICcons,TICmat,indices);
elseif dir==-1
    % for reverse direction
    direction(core)=-1;
    [reacInd,x,stat] = findConsistentReacIDTIC(model,direction,tol,TICcons,TICmat,indices);
end

if stat~=1
    blkdCore=1;flux={};m1={};
else
    blkdCore=0;
    flux = x(1:numel(model.rxns));
    m1 = model.rxns(reacInd);
    model_temp =removeRxns(model,setdiff(model.rxns,m1));
    idx= ismember(model_temp.rxns,model.rxns(core));
    model_temp.lb(idx)=0; model_temp.ub(idx)=0;
    if exist('sprintcc','file')
        a = sprintcc(model_temp,tol); 
    else
        a = fastcc(model_temp,tol);
    end
    if numel(a)~= 0
        model_temp.c= zeros(numel(model_temp.rxns),1); 
        a_temp =a(model_temp.lb(a)==0);
        if isempty(a_temp)
            a_temp=a;
        end
        if numel(a_temp)==1
            idx2=a_temp;
        else
            idx2 = randsample(a_temp,1);
        end
        model_temp.c(idx2)=1;
        sol = optimizeCbModel(model_temp);
        [ids,locb] = ismember(model.rxns,model_temp.rxns(a));
        m1 = model.rxns(ids);
        temp = zeros(numel(model.rxns),1);
        temp2 = sol.x(a);
        temp(ids) = temp2(locb(locb~=0));
        flux = temp;
    end
end
end



function [reacInd,x,stat] = findConsistentReacIDTIC(model,direction,tol,TICcons,TICmat,indices)

[m,n] = size(model.S);
dir0 = direction==0;
n_ = sum(dir0);
% objective
f = [zeros(n,1);ones(n_,1);0];

% equalities
Aeq = [model.S, sparse(m,n_+1)];
beq = zeros(m,1);
csenseeq = repmat('E',m,1); % equality

% inequalities
temp1 = speye(n);
temp2 = -1*spdiag(model.lb(dir0));
Aineq1 = [temp1(dir0,:),temp2,sparse(n_,1)];
bineq1 = zeros(n_,1);
csenseineq1 = repmat('G',n_,1); % greater than

temp2 = -1*spdiag(model.ub(dir0));
Aineq2 = [temp1(dir0,:),temp2,sparse(n_,1)];
bineq2 = zeros(n_,1);
csenseineq2 = repmat('L',n_,1); % lesser than


Aineq3=[];bineq3=[];csenseineq3=[];
if ~isempty(TICcons)
    for i=1:size(TICcons,1)
        temp1 = sparse(1,n);
        Tids = TICcons{i,1};
        temp2 = sparse(1,n); temp2(Tids) = ones(numel(Tids),1);
        temp2 = temp2(dir0);
        Aineq3 = [Aineq3;temp1,temp2,0];
        bineq3 = [bineq3;sum(temp2)-1];
        csenseineq3 = [csenseineq3;'L']; % lesser than
    end
end


if ~isempty(TICmat)
    TICmat_o = orth(TICmat);
    P = sparse(TICmat_o*TICmat_o');
    temp1 = rand(1,size(P,2))*(speye(size(P,2))-P);
    temp2 = sparse(1,n_);

    Aineq5=[temp1(indices),temp2,-1000];bineq5=-1;
    csenseineq5 = 'L';

    Aineq4=[temp1(indices),temp2,-1000];bineq4=tol-1000;
    csenseineq4 = 'G';
else
    Aineq5 =[];bineq5=[];csenseineq5=[];
    Aineq4 =[];bineq4=[];csenseineq4=[];
end



% bounds
lb = model.lb;
lb(direction==1)=max([tol*ones(sum(direction==1),1),lb(direction==1)],[],2);
lb = [lb;zeros(n_+1,1)];
ub = model.ub;
ub(direction==-1)=-tol*ones(sum(direction==-1),1);
ub = [ub;ones(n_+1,1)];

% Set up LP problem
MILPproblem.A=[Aeq;Aineq1;Aineq2;Aineq3;Aineq4;Aineq5];
MILPproblem.b=[beq;bineq1;bineq2;bineq3;bineq4;bineq5];
MILPproblem.lb=lb;
MILPproblem.ub=ub;
MILPproblem.c=f;
MILPproblem.osense=1;%minimise
MILPproblem.vartype = [repmat('C',n,1);repmat('B',n_+1,1)];
MILPproblem.csense = [csenseeq; csenseineq1; csenseineq2; csenseineq3;csenseineq4;csenseineq5];
solution = solveCobraMILP(MILPproblem);
stat=solution.stat;
if stat~=1
    x = [];
    reacInd = [];
else
    x=solution.full;
    reacInd = abs(x(1:n))>=tol*1e-3;
end
end