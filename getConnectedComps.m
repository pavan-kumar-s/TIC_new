%% to get the count of connected components from the Nsnp matrix
load('Recon_nsnp_all_reversible')
Adj_mat = zeros(size(Nsnp,1),size(Nsnp,1));
for row =1:size(Nsnp,1)
    cols = find(Nsnp(row,:));
    Adj_mat(row,sum(abs(Nsnp(:,cols)),2)~=0) = 1;
end
Adj_mat = Adj_mat - diag(diag(Adj_mat));
[bins,binsizes]=conncomp(graph(Adj_mat));
