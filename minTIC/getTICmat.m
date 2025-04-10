function TICmat = getTICmat(TICs,Direction,rxns)
TICmat = zeros(numel(rxns),numel(TICs));
for t =1:numel(TICs)
    temp=TICs{t};
    temp2 = Direction{t};
    [loca,locb] = ismember(rxns,temp);
    TICmat(loca,t)= temp2(locb(locb~=0));
end
end