function [m1,blkdCore,flux,stat] = getTICModel(model,core,tol,dir,TICcons,k)
[~,n] = size(model.S); % number of metabolites and reactions

direction = zeros(n,1);
if dir ==1
    % for forward direction
    direction(core)=1;
    [reacInd,x,stat] = findConsistentReacIDTIC(model,direction,tol,TICcons,k);
elseif dir==-1
    % for reverse direction
    direction(core)=-1;
    [reacInd,x,stat] = findConsistentReacIDTIC(model,direction,tol,TICcons,k);
end

if stat~=1
    blkdCore=1;flux={};m1={};
else
    blkdCore=0;
    flux = x(1:numel(model.rxns));
    m1 = model.rxns(reacInd);
end


end



function [reacInd,x,stat] = findConsistentReacIDTIC(model,direction,tol,TICcons,k)
eps=1e-3;
[m,n] = size(model.S);
dir0 = direction==0;
n_ = sum(dir0);
n_rev = sum(model.rev(dir0)); % number of reversible reactions that are non-core
% objective
f = [zeros(n,1);ones(n_+n_rev,1)];

% equalities
Aeq = [model.S, sparse(m,n_+n_rev)];
beq = zeros(m,1);
csenseeq = repmat('E',m,1); % equality

% inequalities for irreversible reactions
temp1 = speye(n);
temp2 = -eps*spdiag(ones(sum(dir0),1));
try
    Aineq1 = [temp1(dir0 & ~model.rev,:),temp2(~model.rev(dir0),:),sparse(sum(dir0 & ~model.rev),n_rev)];
catch
    asd=1;
end
bineq1 = zeros(sum(dir0 & ~model.rev),1);
csenseineq1 = repmat('G',sum(dir0 & ~model.rev),1); % greater than

temp2 = -1*spdiag(model.ub(dir0));
Aineq2 = [temp1(dir0 & ~model.rev,:),temp2(~model.rev(dir0),:),sparse(sum(dir0 & ~model.rev),n_rev)];
bineq2 = zeros(sum(dir0 & ~model.rev),1);
csenseineq2 = repmat('L',sum(dir0 & ~model.rev),1); % lesser than

% inequalities for reversible reactions
temp1 = speye(n);
temp2 = -eps*spdiag(ones(sum(dir0),1));
temp3 = 1000*spdiag(ones(sum(dir0&model.rev),1));
Aineq3 = [temp1(dir0 & model.rev,:),temp2(model.rev(dir0),:),temp3];
bineq3 = zeros(sum(dir0 & model.rev),1);
csenseineq3 = repmat('G',sum(dir0 & model.rev),1); % greater than

temp2 = -1000*spdiag(ones(sum(dir0),1));
temp3 = eps*spdiag(ones(sum(dir0&model.rev),1));
Aineq4 = [temp1(dir0 & model.rev,:),temp2(model.rev(dir0),:),temp3];
bineq4 = zeros(sum(dir0 & model.rev),1);
csenseineq4 = repmat('L',sum(dir0 & model.rev),1); % greater than

% constraints on the reversible reactions
temp2 = speye(sum(dir0));
Aineq5 = [sparse(n_rev,n),temp2(model.rev(dir0),:),speye(n_rev)];
bineq5 = ones(n_rev,1);
csenseineq5 = repmat('L',n_rev,1);

% TIC constraints
Aineq6=[];bineq6=[];csenseineq6=[];
if ~isempty(TICcons)
    for i=1:size(TICcons,1)
        temp1 = sparse(1,n);
        Tids = TICcons{i,1};
        temp2 = sparse(1,n); temp2(Tids) = ones(numel(Tids),1);
        temp2 = temp2(dir0);
        temp3 = sparse(1,n); temp3(Tids) = ones(numel(Tids),1);
        temp3 = temp3(dir0&model.rev);
        Aineq6 = [Aineq6;temp1,temp2,temp3];
        bineq6 = [bineq6;sum(temp2)-1];
        csenseineq6 = [csenseineq6;'L']; % lesser than
    end
end

% inequality constraint on the total number of reactions that has to be in a TIC
Aineq7 = f';
bineq7 = k-1 ;
csenseineq7 = 'G';

% bounds
lb = model.lb;
lb(direction==1)=max([tol*ones(sum(direction==1),1),lb(direction==1)],[],2);
lb = [lb;zeros(n_+n_rev,1)];
ub = model.ub;
ub(direction==-1)=-tol*ones(sum(direction==-1),1);
ub = [ub;ones(n_+n_rev,1)];

% Set up MILP problem
MILPproblem.A=[Aeq;Aineq1;Aineq2;Aineq3;Aineq4;Aineq5;Aineq6;Aineq7];
MILPproblem.b=[beq;bineq1;bineq2;bineq3;bineq4;bineq5;bineq6;bineq7];
MILPproblem.lb=lb;
MILPproblem.ub=ub;
MILPproblem.c=f;
MILPproblem.osense=1;%minimise
MILPproblem.vartype = [repmat('C',n,1);repmat('B',n_+n_rev,1)];
MILPproblem.csense = [csenseeq; csenseineq1; csenseineq2; csenseineq3; csenseineq4; csenseineq5; csenseineq6; csenseineq7];
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