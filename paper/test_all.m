% obtain solutions with cpas and penlab
clc
clear
close all
addpath ('/opt/ibm/ILOG/CPLEX_Studio1210/cplex/matlab/x86-64_linux')

[coord, elem, dof, adof, E, rho, f] = data_test(2);
n = length(f);
m = size(elem,1);
draw_truss(coord,elem,ones(m,1),[2 1],0.0);

K0 = zeros(n,n);
Kdx = mat_Kdx_vol(E*ones(m,1),elem,coord,dof,adof);
M0 = zeros(n,n);
Mdx = mat_Mdx_vol(rho*ones(m,1),elem,coord,dof,adof);
Kfun = @(x) sum_mat(1,m,K0,Kdx,x);
Mfun = @(x) sum_mat(1,m,M0,Mdx,x);

x0 = ones(m,1);
gamma = mincompl(Kfun(x0),f);
lambda0 = maxpsdpair(Kfun(x0),Mfun(x0));

data.gamma = gamma;
data.f = f;
data.K = Kfun;
data.K0 = K0;
data.Kdx = Kdx;
data.M = Mfun;
data.M0 = M0;
data.Mdx = Mdx;

% <CPAS>
fprintf('************************************************\n')
fprintf('gamma=%1.1e\n',data.gamma)
fprintf('min v s.t. [gamma -f'';-f K(x)]>=0, sum(x)<v and x>=0\n')
c =  [zeros(m,1); 1];
A =  [ones(1,m)  -1];
b =  0;
Volume = 0;           % lower bound for volume*
lb = [zeros(m,1); Volume];
ub = [];
data.ncut = 1;
data.cut{1} = @(x,data) compliance_lmi_cut(x,data);
data.TOL = 1e-6;
data.options = cplexoptimset('cplex');
[x1,A,b,exitflag] = cpas(c,A,b,lb,ub,data,1);
if exitflag <= 0
    keyboard
end
%draw_truss(coord,elem,x1,[1 1],0.0);
lambda1 = maxpsdpair(Kfun(x1),Mfun(x1),1e-6,0);

fprintf('*********************************************************************\n')
data.lambda = lambda1 + 1e-1;
fprintf('gamma=%1.1e, lambda=%1.1e\n',data.gamma,data.lambda)
fprintf('min v s.t. [gamma -f'';-f A(x)]>=0, sum(x)<v, A(x)-lambda*B(x)>=0 and x>=0\n')
c =  [zeros(m,1); 1];
A =  [ones(1,m)  -1];
b =  0;
Volume = 0;           % lower bound for volume*
lb = [zeros(m,1); Volume];
ub = [];
data.ncut = 2;
data.cut{1} = @(x,data) compliance_lmi_cut(x,data);
data.cut{2} = @(x,data) frequency_lmi_cut(x,data);
data.TOL = 1e-6;
data.options = cplexoptimset('cplex');
[x2,A,b,exitflag] = cpas(c,A,b,lb,ub,data,1);
if exitflag <= 0
    keyboard
end
%draw_truss(coord,elem,x2,[1 1],1e-2);
%lambda2 = maxpsdpair(Kfun(x2),Mfun(x2));

fprintf('*********************************************************************************\n')
data.vol = sum(x2(1:m));
fprintf('gamma=%1.1e, V=%1.1e\n',data.gamma,data.vol)
fprintf('min -lambda s.t. [gamma -f'';-f A(x)]>=0, A(x)-lambda*B(x)>=0, sum(x)<=V and x>=0\n')
c =  [zeros(m,1);-1];
A =  [ones(1,m)   0];
b =  data.vol;
Lambda = 100;       % upper bound for lambda*
lb = [zeros(m,1);    -Inf];
ub = [ones(m,1)*Inf; Lambda];  
data.ncut = 2;
data.cut{1} = @(x,data) compliance_lmi_cut(x,data);
data.cut{2} = @(x,data) frequency_bmi_cut(x,data);
data.TOL = 1e-6;
data.options = cplexoptimset('cplex');
[x3,A,b,exitflag] = cpas(c,A,b,lb,ub,data,1);
if exitflag <= 0
    keyboard
end
%draw_truss(coord,elem,x3,[1 1],1e-2);
%lambda3 = maxpsdpair(Kfun(x3(1:m)),Mfun(x3(1:m)));
% </GCPS >





% <PENLAB>
fprintf('************************************************\n')
fprintf('gamma=%1.1e\n',data.gamma)
fprintf('min v s.t. [gamma -f'';-f A(x)]>=0, sum(x)<=v and x>=0\n')
penm = [];
penm.Nx = m;
% objective function f
penm.objfun  = @(x,Y,ud) deal(ones(m,1)'*x, ud);
penm.objgrad = @(x,Y,ud) deal(ones(m,1), ud);
penm.objhess = @(x,Y,ud) deal(sparse(m,m), ud);
% matrix constraints A_{k}
penm.NALIN = 1; % number of Linear Matrix Ineq constraints (LMI)
penm.lbA = zeros(penm.NALIN,1); % A_i(x) must be positive semidef
penm.mconfun  = @(x,Y,k,ud)     deal(mconfun( x,k,data),ud);
penm.mcongrad = @(x,Y,k,i,ud)   deal(mcongrad(x,k,i,data),ud);
% x>=0
penm.lbx = zeros(m,1); % lower bound of x
p1 = penlab(penm);
p1.opts = struct('outlev',0);
exitflag = p1.solve();
if exitflag ~= 1
    keyboard
end
fprintf('*********************************************************************\n')
fprintf('gamma=%1.1e, lambda=%1.1e\n',data.gamma,data.lambda)
fprintf('min v s.t. [gamma -f'';-f A(x)]>=0, A(x)-lambda*B(x)>=0, sum(x)<=v and x>=0\n')
penm = [];
penm.Nx = m;
% objective function f
penm.objfun  = @(x,Y,ud) deal(ones(m,1)'*x, ud);
penm.objgrad = @(x,Y,ud) deal(ones(m,1), ud);
penm.objhess = @(x,Y,ud) deal(sparse(m,m), ud);
% matrix constraints A_{k}
penm.NALIN = 2; % number of Linear Matrix Ineq constraints (LMI)
penm.lbA = zeros(penm.NALIN,1); % A_i(x) must be positive semidef
penm.mconfun  = @(x,Y,k,ud)     deal(mconfun( x,k,data),ud);
penm.mcongrad = @(x,Y,k,i,ud)   deal(mcongrad(x,k,i,data),ud);
% x>=0
penm.lbx = zeros(m,1); % lower bound of x
p2 = penlab(penm);
p2.opts = struct('outlev',0);
exitflag = p2.solve();
if exitflag ~= 1
    keyboard
end
fprintf('*********************************************************************************\n')
fprintf('gamma=%1.1e, V=%1.1e\n',data.gamma,data.vol)
fprintf('min -lambda s.t. [gamma -f'';-f A(x)]>=0, A(x)-lambda*B(x)>=0, sum(x)<=V and x>=0\n')
penm = [];
penm.Nx = m+1;
% objective function f
penm.objfun  = @(x,Y,ud) deal([zeros(m,1);-1]'*x, ud);
penm.objgrad = @(x,Y,ud) deal([zeros(m,1);-1],    ud);
penm.objhess = @(x,Y,ud) deal(sparse(m+1,m+1),  ud);
% constraint function g
penm.NgLIN = 1; % linear, in this case
penm.confun  = @(x,Y,ud) deal([ones(m,1);0]'*x,ud);
penm.congrad = @(x,Y,ud) deal([ones(m,1);0]',  ud);
penm.ubg = data.vol;    % upper bound of g
penm.lbg = -Inf;        % lower bound of g
% matrix constraints A_{k}
penm.NANLN = 1; % number of Bilinear Matrix Ineq constraints (BMI)
penm.NALIN = 1; % number of Linear Matrix Ineq constraints (LMI)
penm.lbA = zeros(penm.NANLN+penm.NALIN,1); % A_i(x) must be positive semidef
penm.mconfun  = @(x,Y,k,ud)     deal(mconfun3( x,k,data),ud);
penm.mcongrad = @(x,Y,k,i,ud)   deal(mcongrad3(x,k,i,data),ud);
penm.mconhess = @(x,Y,k,i,j,ud) deal(mconhess3(x,k,i,j,data),ud);
% x>=0
penm.lbx = [zeros(m,1);-Inf]; % lower bound of x, -Inf gives no lang mult
bigLambPEN = 100;
penm.ubx = [+inf(m,1); bigLambPEN];
p3 = penlab(penm);
p3.opts = struct('outlev',0);%2);
exitflag = p3.solve();
if exitflag ~= 1
    keyboard
end

% KKT
% Lagrange mult. (OMG!! penlab gives positive semidef. Lag. multipliers.!!!)
lambC1 = -p1.UA{1}; % matrix, [gamma -f';-f A(x(1:m))]>=0
lambX1 = -p1.uxbox; % vector, x(1:m)>=0
probx1 =  p1.x;
lambC2 = -p2.UA{1}; % matrix, [gamma -f';-f A(x(1:m))]>=0
lambF2 = -p2.UA{2}; % matrix, A(x(1:m))-lambda*B(x(1:m))>=0
lambX2 = -p2.uxbox; % vector, x(1:m)>=0
probx2 =  p2.x;
lambC3 = - p3.UA{2}; % matrix, [gamma -f';-f A(x(1:m))]>=0
lambF3 = - p3.UA{1}; % matrix, A(x(1:m))-x(m+1)*B(x(1:m))>=0
lambX3 = - p3.uxbox; % vector, x(1:m)>=0
lambV3 = + p3.uineq; % scalar, sum(x(1:m))-V<=0 (no sign change)
probx3 =  p3.x;
% KKTPvolWOfreq(probx1,lambC1,lambX1,Kfun,Kdx,gamma,f)
% KKTPvol(probx2,lambC2,lambF2,lambX2,Kfun,Mfun,Kdx,Mdx,gamma,f,data.lambda)
% KKTPeig(probx3,lambC3,lambF3,lambV3,lambX3,Kfun,Mfun,Kdx,Mdx,gamma,f,data.vol)

% % comparing solution
zeroval = 1e-3;
% xtopo1 = x1(1:m);
% xtopo2 = x2(1:m);
% xtopo3 = x3(1:m);
% xtopo1(xtopo1<zeroval) = 0;
% xtopo2(xtopo2<zeroval) = 0;
% xtopo3(xtopo3<zeroval) = 0;
% probxtopo1 = probx1;
% probxtopo2 = probx2;
% probxtopo3 = probx3(1:m);
% probxtopo1(probxtopo1<zeroval) = 0;
% probxtopo2(probxtopo2<zeroval) = 0;
% probxtopo3(probxtopo3<zeroval) = 0;
noopt =          [x0(1:m);   sum(x0(1:m))];%   maxpsdpair(Kfun(x0),Mfun(x0))];
PvolwofCPAS =    [x1(1:m);   sum(x1(1:m))];%   maxpsdpair(Kfun(x1),Mfun(x1))];
PvolCPAS =       [x2(1:m);   sum(x2(1:m))];%   maxpsdpair(Kfun(x2),Mfun(x2))];
PeigCPAS =       [x3(1:m);   sum(x3(1:m))];%   maxpsdpair(Kfun(x3),Mfun(x3))];
PvolwofPEN =     [p1.x(1:m); sum(p1.x(1:m))];% maxpsdpair(Kfun(p1.x),Mfun(p1.x))];
PvolPEN  =       [p2.x(1:m); sum(p2.x(1:m))];% maxpsdpair(Kfun(p2.x),Mfun(p2.x))];
PeigPEN  =       [p3.x(1:m); sum(p3.x(1:m))];% maxpsdpair(Kfun(p3.x),Mfun(p3.x))];
table(noopt, PvolwofCPAS, PvolwofPEN, PvolCPAS, PvolPEN, PeigCPAS, PeigPEN)

draw_truss(coord,elem,x0,[1,1],zeroval);

elem1 = elem;
elem1(x1<zeroval,:) = [];
xdraw1 = x1;
xdraw1(x1<zeroval) = [];
draw_truss(coord,elem1,xdraw1,[1,1],zeroval);

elem2 = elem;
elem2(x2<zeroval,:) = [];
xdraw2 = x2;
xdraw2(x2<zeroval) = [];
draw_truss(coord,elem2,xdraw2,[1,1],zeroval);

elem3 = elem;
elem3(x3<zeroval,:) = [];
xdraw3 = x3;
xdraw3(x3<zeroval) = [];
draw_truss(coord,elem3,xdraw3,[1,1],zeroval);



% penlab help functions
function Ak = mconfun(x,k,data)
    if k == 1
        Ak = sparse([data.gamma -data.f';-data.f data.K(x)]);
    elseif k == 2
        Ak = sparse(data.K(x) - data.lambda*data.M(x));
    end
end
function Aki = mcongrad(~,k,i,data)
    m = size(data.f,1);
    if k == 1
        Aki = sparse([0 zeros(1,m); zeros(m,1) data.Kdx{i}]);
    elseif k == 2
        Aki = sparse(data.Kdx{i} - data.lambda*data.Mdx{i});
    end
end
function Ak = mconfun3(x,k,data)
    if k == 1
        Ak = sparse(data.K(x) - x(length(x))*data.M(x));
    elseif k == 2
        Ak = sparse([data.gamma -data.f'; -data.f data.K(x)]);
    end
end
function Aki = mcongrad3(x,k,i,data)
    n = length(x); % [x;lambda]
    if k == 1
        if i <= n - 1
            Aki = sparse(data.Kdx{i} - x(n)*data.Mdx{i}); % deriv of A_{k} wrt x(i)
        elseif i == n
            Aki = sparse(- data.M(x)); % deriv of A_{k} wrt x(n)
        end
    elseif k == 2
        m = length(data.f);
        if i <= n - 1
            Aki = sparse([0 zeros(1,m); zeros(m,1) data.Kdx{i}]); % deriv of A_{k} wrt x(i)
        elseif i == n
            Aki = sparse(m+1,m+1); % deriv of A_{k} wrt x(n)
        end
    end
end
function Akij = mconhess3(x,k,i,j,data)
    n = length(x);
    m = length(data.f);
    if k == 1
        if i <= n - 1
            if j <= n - 1
                Akij = sparse(m,m);
            elseif j == n
                Akij = sparse(- data.Mdx{i});
            end
        elseif i == n
            if j <= n - 1
                Akij = sparse(- data.Mdx{j});
            elseif j == n
                Akij = sparse(m,m);
            end
        end
    elseif k == 2
        Akij = sparse(m+1,m+1);
    end
end
% <\PENLAB>



% cpas help functions
function cuts = frequency_lmi_cut(x,data)
    nvar = length(x);
    %[w,smalleig] = eigs(data.K(x) - data.lambda*data.M(x),1,'smallestreal');
    [allEigVec,allEigVal] = eig(data.K(x) - data.lambda*data.M(x));
    [sortedAllEig,sortedIndex] = sort(diag(allEigVal));
    smalleig = sortedAllEig(1);
    w = allEigVec(:,sortedIndex(1));
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
    gamma = data.gamma;
    % [w,smalleig] = eigs(gamma*data.K(x)-f*f',1,'smallestreal');
    [allEigVec,allEigVal] = eig(gamma*data.K(x)-f*f');
    [sortedAllEig,sortedIndex] = sort(diag(allEigVal));
    smalleig = sortedAllEig(1);
    w = allEigVec(:,sortedIndex(1));
    if smalleig >= -data.TOL
        cuts = [];
        return
    end
    cuts = zeros(1,nvar+1);
    for i = 1:nvar-1
        cuts(i) = - data.gamma*w'*data.Kdx{i}*w;
    end
    %cuts(nvar) = 0; % cuts(nvar)*vol
    cuts(nvar+1) = data.gamma*w'*data.K0*w - (f'*w)^2;
end

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
