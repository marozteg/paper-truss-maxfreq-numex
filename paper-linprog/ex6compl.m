clear
close all

load("LP_sol.mat","xLP")
v = sum(xLP);
lambda = 5e-10;
fprintf('given v and lambda\n')
fprintf('min gamma\n')
fprintf('s.t. gamma*K(x)-ff''>=0,\n')
fprintf('     K(x)-lambda*M(x)>=0\n')
fprintf('     sum(x)<=v and\n') 
fprintf('     x>=0.\n')

[coord, elem, dof, adof, E, rho, f] = data_shortcant();
m = size(elem,1);

%draw_truss(coord,elem,ones(m,1),[1 1]);
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

c =  [zeros(m,1); 1];
A =  [ones(1,m)   0];
b =  v;
Gamma = 0;         % lower bound for gamma*
lb = [zeros(m,1); Gamma];
ub = [];
data.ncut = 2;
data.cut{1} = @(x,data) compliance_lmi_cut(x,data);
%data.cut{1} = @(x,data) compliance_bmi_cut(x,data);
data.cut{2} = @(x,data) frequency_lmi_cut(x,data);
data.Kdx = Kdx;
data.Mdx = Mdx;
data.K = @(x) Kfun(x);
data.M = @(x) Mfun(x);
data.K0 = K0;
data.M0 = M0;
data.f = f;
data.vol = v;
data.lambda = lambda;
data.TOL = 1e-6;
data.options = optimoptions('linprog','Display','off');
[xsol,A,b,exitflag] = cpas(c,A,b,lb,ub,data,1);

zeroval = 0;
elempos = elem;
elempos(xsol(1:m)<zeroval,:) = [];
xpos = xsol(1:m);
xpos(xsol(1:m)<zeroval) = [];
draw_truss(coord,elempos,xpos,[1,1],zeroval);

function cuts = frequency_lmi_cut(x,data)
    nvar = length(x);
    %[w,smalleig] = eigs(data.K(x) - data.lambda*data.M(x),1,'smallestreal');
    [W,D] = eig(data.K(x) - data.lambda*data.M(x));
    [sortedeigs,index] = sort(diag(D));
    smalleig = sortedeigs(1);
    w = W(:,index(1));
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
    f = data.f;
    nvar = length(x);
    [w,smalleig] = eigs([x(nvar) -f';-f data.K(x)],1,'smallestreal');
    if smalleig >= -data.TOL
        cuts = [];
        return
    end
    h = w(1);
    w = w(2:end);
    cuts = zeros(1,nvar+1);
    for i = 1:nvar-1
        cuts(i) = - w'*data.Kdx{i}*w;
    end
    cuts(nvar) = - h*h; % cuts(nvar)*gamma
    cuts(nvar+1) = - 2*h*f'*w + w'*data.K0*w;
end

function cuts = compliance_bmi_cut(x,data)
    f = data.f;
    nvar = length(x);
    gamma = x(nvar);
    [w,smalleig] = eigs(gamma*data.K(x)-f*f',1,'smallestreal');
    if smalleig >= -data.TOL
        cuts = [];
        return
    end
    cuts = zeros(1,nvar+1);
    alpha = zeros(nvar-1,1);
    alpha0 = w'*data.K0*w;
    for i = 1:nvar-1
        alpha(i) = w'*data.Kdx{i}*w;
        cuts(i) = - gamma*alpha(i);
    end
    cuts(nvar) = - (max(alpha)*data.vol + alpha0); % cuts(nvar)*gamma
    cuts(nvar+1) = gamma*cuts(nvar) - (w'*f)^2 + gamma*alpha0;
end