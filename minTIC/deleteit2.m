p=dir('./models');
p={p(3:end).name}';
tics = {}; Dir= {}; t1=[];
times=[];
for i=1:numel(p)
    load(['./models/',p{i}])
    mm=getTICmodel(model);
    t1(i) = size(null(mm.S),2);
    tic
    [tics{i},Dir{i}] = ThermOptEnumMILP_minimal(model);
    times(i)=toc;
    temp = tics{i};
    for t =1:numel(temp)
        tempS = model.S(:,ismember(model.rxns,temp{t}));
        if numel(temp{t})-rank(tempS)~=1
            sad=1 %% to demtify the TICs which are not minimal
        end
    end
    TICmat = getTICmat(tics{i},Dir{i},model.rxns);
    if rank(TICmat)~=size(TICmat,2)
        sad=1 %% to see if the TIC matrix is full rank or not
    end
    
    if rank(TICmat)~=size(mm.S,2)-rank(mm.S)
        sad=1 %% to see if the TIC matrix and the null space matrix span the same space
    end
end

function model = getTICmodel(model)
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

tol=1;
ind=findExcRxns(model);
% blocking all the exchange reactions
model.lb(ind) = 0; model.ub(ind) = 0;

if exist('sprintcc','file')
    a = sprintcc(model,tol); 
else
    a = fastcc(model,tol);
end

% model with only TICs
model = removeRxns(model,model.rxns(setdiff([1:numel(model.rxns)],a)));
end