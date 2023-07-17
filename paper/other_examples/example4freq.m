clear
close all

v = 1.20229;
gamma = 1;
fprintf('given v and gamma\n')
fprintf('max lambda\n')
fprintf('s.t. gamma*K(x)-ff''>=0,\n') 
fprintf('     K(x)-lambda*M(x)>=0,\n') 
fprintf('     sum(x)<=v and\n') 
fprintf('     x>=0.\n')

[coord, elem, dof, adof, E, rho, f] = data_example41();
m = size(elem,1);

draw_elem(coord,elem,ones(m,1));
cathetus = coord(elem(:,2),:)-coord(elem(:,1),:); % nel x 3
elem_length = vecnorm(cathetus'); % 1 x nel

e = eye(m);
n = size(f,1);
Kdx = cell(m,1);
Mdx = cell(m,1);
for i = 1:m
    area = 1/elem_length(i)*e(:,i); % area for unitary volume of elem i
    Kdx{i} = mat_K(E*ones(m,1),area,elem,coord,dof,adof);
    mass = rho*1*e(:,i); % mass for unitary volume of elem i 
    Mdx{i} = mat_M(mass,elem,coord,dof,adof);
end
K0 = zeros(n,n);
M0 = zeros(n,n);
Kfun = @(x) sum_mat(1,m,K0,Kdx,x);
Mfun = @(x) sum_mat(1,m,M0,Mdx,x);

c =  [zeros(m,1);-1];
A =  [ones(1,m)   0];
b =  v;
Lambda = 100;       % upper bound for lambda*
lb = [zeros(m,1);    -Inf];
ub = [ones(m,1)*Inf; Lambda];  
data.ncut = 2;
data.cut{1} = @(x,data) compliance_lmi_cut(x,data);
data.cut{2} = @(x,data) frequency_bmi_cut(x,data);
data.Kdx = Kdx;
data.Mdx = Mdx;
data.K = @(x) Kfun(x);
data.M = @(x) Mfun(x);
data.K0 = K0;
data.M0 = M0;
data.f = f;
data.gamma = gamma;
data.vol = v;
data.TOL = 1e-5;
data.options = cplexoptimset('cplex');
[xsol,A,b,exitflag] = cpas(c,A,b,lb,ub,data,0);

zeroval = 1e-3;
elempos = elem;
elempos(xsol(1:m)<zeroval,:) = [];
xpos = xsol(1:m);
xpos(xsol(1:m)<zeroval) = [];
draw_elem(coord,elempos,xpos);


function cuts = frequency_bmi_cut(x,data)
    nvar = length(x);
    lambda = x(nvar);
    [w,smalleig] = eigs(data.K(x)-lambda*data.M(x),1,'smallestreal');
    if smalleig >= -data.TOL
        cuts = [];
        return
    end
    cuts = zeros(1,nvar+1);
    beta = zeros(nvar-1,1);
    beta0 = w'*data.M0*w;
    for i = 1:nvar-1
        beta(i) =  w'*data.Mdx{i}*w;
        cuts(i) = - (w'*data.Kdx{i}*w - lambda*beta(i));
    end
    cuts(nvar) = max(beta)*data.vol + beta0; % cuts(nvar)*lambda
    cuts(nvar+1) = lambda*cuts(nvar) + w'*data.K0*w - lambda*beta0;
end

function cuts = compliance_lmi_cut(x,data)
    f = data.f;
    nvar = length(x);
    [w,smalleig] = eigs(data.gamma*data.K(x)-f*f',1,'smallestreal');
    if smalleig >= -data.TOL
        cuts = [];
        return
    end
    cuts = zeros(1,nvar+1);
    alpha = zeros(nvar,1);
    alpha0 = w'*data.K0*w;
    for i = 1:nvar-1
        alpha(i) = w'*data.Kdx{i}*w;
        cuts(i) = - data.gamma*alpha(i);
    end
    %cuts(nvar) = 0; % cuts(nvar)*lambda
    cuts(nvar+1) = - (w'*f)^2 + data.gamma*alpha0;
end