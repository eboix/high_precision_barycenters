% This is a sample of the input format for the Wasserstein_barycenter.m
% It aims to manipulate the data of images (30 x 50 x 50) of nested ellipses 
load('inputs.mat', 'ellipses');

N = 30;
d = 2;
db=cell(1,1);
w = cell(N,1);
db{1}.stride = zeros(1,N);

for i=1:N
    v  = find(ellipses(i,:,:));
    db{1}.stride(i) = length(v);
end

stride_cumsum = [0,cumsum(db{1}.stride)] ;
db{1}.w = zeros(1,stride_cumsum(end));
db{1}.supp = zeros(d,stride_cumsum(end));

for i=1:N
    A = reshape(ellipses(i,:,:),[50,50]);
    [row,col,v]  = find(A);
    db{1}.supp(:,stride_cumsum(i)+1:stride_cumsum(i+1)) = [row';col'];
    db{1}.w(stride_cumsum(i)+1:stride_cumsum(i+1)) = ones(1,db{1}.stride(i))./db{1}.stride(i);
end

c0 = cell(1,1);
c0{1}.w = ones(1,2500)./2500;

row = reshape(mod(0:2500-1,50)+1,1,[]);
col = kron(1:50,ones(1,50));
c0{1}.supp = [row;col];

options.method='fixed_maaipm'; % {'ibp','gurobi','badmm', 'admm', 'ibp'}
options.ipmouttolog = 1;
options.ipmtol_primal_dual_gap = 1e-3; 
options.largem = 1;
[c, OT]=Wasserstein_Barycenter(db, c0, options);

C = reshape(c{1}.w,50,50)*10000;
image(C)
colorbar