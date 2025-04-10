temp = null(mm.S);
TICmat = getTICmat(tics{12},Dir{12},mm.rxns);
TICmat = round(TICmat,8);
is=[];
for i=1:size(TICmat,2)
   if rank([temp,TICmat(:,i)])>60
       is=[is;i];
   end
   
end
    