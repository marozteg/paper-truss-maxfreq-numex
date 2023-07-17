clear
close all

gamma = 20;
lambda = 1e-3;
fprintf('given gamma and lambda\n')
fprintf('min v\n')
fprintf('s.t. gamma*K(x)-ff''>=0,\n')
fprintf('     K(x)-lambda*M(x)>=0\n')
fprintf('     sum(x)<=v and\n') 
fprintf('     x>=0.\n')

[coord, elem, dof, adof, E, rho, f] = data_example47();
m = size(elem,1);

draw_truss();
cathetus = coord(elem(:,2),:)-coord(elem(:,1),:); % nel x 3
elem_length = vecnorm(cathetus'); % 1 x nel

e = eye(m);
n = size(f,1);
Kdx = cell(m,1);
Mdx = cell(m,1);
for i = 1:m
    area = 1/elem_length(i)*e(:,i); % area for unitary volume of elem i
    Kdx{i} = sparse(mat_K(E*ones(m,1),area,elem,coord,dof,adof));
    mass = rho*1*e(:,i); % mass for unitary volume of elem i 
    Mdx{i} = sparse(mat_M(mass,elem,coord,dof,adof));
end
K0 = sparse(zeros(n,n));
M0 = sparse(zeros(n,n));
Kfun = @(x) sparse(sum_mat(1,m,K0,Kdx,x));
Mfun = @(x) sparse(sum_mat(1,m,M0,Mdx,x));

c =  [zeros(m,1); 1];
A =  [ones(1,m)  -1];
b =  0;
Volume = 0;           % lower bound for volume*
lb = [zeros(m,1); Volume];
ub = [];
data.ncut = 2;
data.cut{1} = @(x,data) compliance_lmi_cut(x,data);
data.cut{2} = @(x,data) frequency_lmi_cut(x,data);
data.Kdx = Kdx;
data.Mdx = Mdx;
data.K = @(x) Kfun(x);
data.M = @(x) Mfun(x);
data.K0 = K0;
data.M0 = M0;
data.f = f;
data.gamma = gamma;
data.lambda = lambda;
data.TOL = 1e-6;
data.options = cplexoptimset('cplex');
[xsol,A,b,exitflag] = cpas(c,A,b,lb,ub,data,1);

zeroval = 5e-3;
elempos = elem;
elempos(xsol(1:m)<zeroval,:) = [];
xpos = xsol(1:m);
xpos(xsol(1:m)<zeroval) = [];
draw_elem(coord,elempos,xpos);

function cuts = frequency_lmi_cut(x,data)
    nvar = length(x);
    [w,smalleig] = eigs(data.K(x) - data.lambda*data.M(x),1,'smallestreal');
    if smalleig >= -data.TOL
        cuts = [];
        return
    end
    cuts = zeros(1,nvar+1);
    for i = 1:nvar-1
        cuts(i) = - w'*(data.Kdx{i} - data.lambda*data.Mdx{i})*w;
    end
    %cuts(nvar) = 0; % cuts(nvar)*vol
    cuts(nvar+1) = w'*(data.K0 - data.lambda*data.M0)*w;
end

function cuts = compliance_lmi_cut(x,data)
    F = sparse(data.f*(data.f)');
    nvar = length(x);
    gamma = data.gamma;
    [w,smalleig] = eigs(gamma*data.K(x)-F,1,'smallestreal');
    if smalleig >= -data.TOL
        cuts = [];
        return
    end
    cuts = zeros(1,nvar+1);
    for i = 1:nvar-1
        cuts(i) = - data.gamma*w'*data.Kdx{i}*w;
    end
    %cuts(nvar) = 0; % cuts(nvar)*vol
    cuts(nvar+1) = data.gamma*w'*data.K0*w - (data.f'*w)^2;
end