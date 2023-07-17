clear
close all

v = 20;
gamma = 9.816573522e+00;
fprintf('given v and gamma\n')
fprintf('max lambda\n')
fprintf('s.t. gamma*K(x)-ff''>=0,\n') 
fprintf('     K(x)-lambda*M(x)>=0,\n') 
fprintf('     sum(x)<=v and\n') 
fprintf('     x>=0.\n')

[coord, elem, dof, adof, E, rho, f] = data_bridge_21();
m = size(elem,1);

draw_truss(coord,elem,ones(m,1),[1,1]);
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
data.options = optimoptions('linprog','Display','off');
[xcpas,A,b,exitflag] = cpas(c,A,b,lb,ub,data,1);






% <PENLAB>
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
prob = penlab(penm);
prob.opts = struct('outlev',2);%0);
exitflag2 = prob.solve();
if exitflag2 ~= 1
    keyboard
end
probx = prob.x;



draw_truss(coord,elem,xcpas,[1,1],1e-3);



% comparing solution
xcpas(abs(xcpas)<1e-4) = 0;
CPAS =  [xcpas(1:m);  sum(xcpas(1:m));  mincompl(Kfun(xcpas),f);  maxpsdpair(Kfun(xcpas),Mfun(xcpas))];
probx(abs(probx)<1e-4) = 0;
PEN  =  [probx(1:m); sum(probx(1:m)); mincompl(Kfun(probx),f); maxpsdpair(Kfun(probx),Mfun(probx))];
table(CPAS, PEN)


% penlab help functions
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