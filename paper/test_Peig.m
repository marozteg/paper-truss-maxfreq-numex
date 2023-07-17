% obtain solutions with cpas and penlab
clc
clear
close all
addpath('/opt/ibm/ILOG/CPLEX_Studio1210/cplex/matlab/x86-64_linux')
addpath('/home/marozteg/MATLAB-Drive/tools/')
addpath('/home/marozteg/MATLAB-Drive/paper/examples')
% to use PENLAB run install.m in PENLAB folder

prob_sizes = 2;
n_arr = zeros(prob_sizes,1);
m_arr = zeros(prob_sizes,1);
time_CPAS_arr = zeros(prob_sizes,1);
time_PEN_arr = zeros(prob_sizes,1);
vol_CPAS_arr = zeros(prob_sizes,1);
vol_PEN_arr = zeros(prob_sizes,1);
compl_CPAS_arr = zeros(prob_sizes,1);
compl_PEN_arr = zeros(prob_sizes,1);
lambda_CPAS_arr = zeros(prob_sizes,1);
lambda_PEN_arr  = zeros(prob_sizes,1);
x_diff_arr = zeros(prob_sizes,1);
x_CPAS_cell = cell(prob_sizes,1);
x_PEN_cell = cell(prob_sizes,1);

%problems = {"data_shortcant","data_shortcant2","data_mbbbeam"};
problems = {"data_hempcant"};%"data_halfwheel"};%,"data_hempcant"};
%problems = {"data_benchmark01"};
nprob = size(problems,2);

for jprob = 1:nprob
    problem = problems{jprob};
    data_problem = str2func(problem);
    
    for prob_size = 1:prob_sizes % prob_size must be 1,2,3,...
        [coord, elem, dof, adof, E, rho, f] = data_problem(prob_size);
    
        n = length(f);      n_arr(prob_size) = n;
        m = size(elem,1);   m_arr(prob_size) = m;
       
        fprintf('%s, prob_size=%d, m=%d, n=%d\n',problem,prob_size,m,n)
        % draw_truss(coord,elem,ones(m,1),[1 2],0.0);
        % continue

        K0 = zeros(n,n);
        Kdx = mat_Kdx_vol(E*ones(m,1),elem,coord,dof,adof);
        M0 = zeros(n,n);
        Mdx = mat_Mdx_vol(rho*ones(m,1),elem,coord,dof,adof);
        Kfun = @(x) sum_mat(1,m,K0,Kdx,x);
        Mfun = @(x) sum_mat(1,m,M0,Mdx,x);
        
        % <DATA>
        x0 = ones(m,1);
     %    gammas = [     1.303490445540362e-01
     % 1.186195714574521e-01
     % 1.089519728804086e-01
     % 1.010546564679102e-01
     % 9.453950426800813e-02];
        data.gamma = mincompl(Kfun(x0),f);%gammas(prob_sizes);%mincompl(Kfun(x0),f);
%         volumes = [     8.533834209779938e+02
%      1.278523748879444e+03
%      1.895387248516710e+03
%      2.706347246383653e+03
%      3.721769834336337e+03
% ];
        data.vol = sum(x0);%volumes(prob_sizes);%sum(x0);
        data.f = f;
        % </DATA>
        
        data.K = Kfun;
        data.K0 = K0;
        data.Kdx = Kdx;
        data.M = Mfun;
        data.M0 = M0;
        data.Mdx = Mdx;
                
        % <CPAS>
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
    
        tic
        [xcpas,A,b,exitflag1] = cpas(c,A,b,lb,ub,data,1);
        time_CPAS_arr(prob_size) = toc;

        x_CPAS_cell{prob_size} = xcpas;
    
        if exitflag1 <= 0
            %keyboard
            vol_CPAS_arr(prob_size) = +Inf;
            compl_CPAS_arr(prob_size) = +Inf;
            lambda_CPAS_arr(prob_size) = +Inf;
        else
            vol_CPAS_arr(prob_size) = sum(xcpas(1:m));
            compl_CPAS_arr(prob_size) = mincompl(Kfun(xcpas),f);
            lambda_CPAS_arr(prob_size) = xcpas(m+1);
        end
        %draw_truss(coord,elem,xcpas,[1 1],1e-2);
        %lambda3 = maxpsdpair(Kfun(xcpas(1:m)),Mfun(xcpas(1:m)));
        % </CAPS >
        
        
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
    
        tic
        exitflag2 = prob.solve();
        time_PEN_arr(prob_size) = toc;

        xpen = prob.x;
        x_PEN_cell{prob_size} = xpen;

        if exitflag2 ~= 1
            %keyboard
            vol_PEN_arr(prob_size) = +Inf;
            compl_PEN_arr(prob_size) = +Inf;
            lambda_PEN_arr(prob_size) = +Inf;
        else
            vol_PEN_arr(prob_size) = sum(xpen(1:m));
            compl_PEN_arr(prob_size) = mincompl(Kfun(xpen),f);
            lambda_PEN_arr(prob_size) = xpen(m+1);
        end
        % </PENLAB>
        
        % comparing solution
        if exitflag1 <= 0 || exitflag2 ~= 1
            x_diff_arr(prob_size) = +Inf;
        else
            x_diff_arr(prob_size) = norm(xcpas-xpen)/norm(xcpas);
        end
        % noop =     [x0(1:m);     sum(x0(1:m));     mincompl(Kfun(x0),f);     maxpsdpair(Kfun(x0),Mfun(x0))];
        % xcpas(abs(xcpas)<1e-4) = 0;
        % CPAS =  [xcpas(1:m);  sum(xcpas(1:m));  mincompl(Kfun(xcpas),f);  maxpsdpair(Kfun(xcpas),Mfun(xcpas))];
        % xpen(abs(xpen)<1e-4) = 0;
        % PEN  =  [xpen(1:m); sum(xpen(1:m)); mincompl(Kfun(xpen),f); maxpsdpair(Kfun(xpen),Mfun(xpen))];
        %table(noop, CPAS, PEN)
        % [noop CPAS PEN]
        
        % if exitflag1 > 0
        %     draw_truss(coord,elem,xcpas,[1,1],1e-6);
        %     title('CPAS')
        % elseif exitflag2 == 1
        %     draw_truss(coord,elem,xpen,[1,1],1e-6);
        %     title('PENLAB')
        % end

    end % prob_size
    m = m_arr; n = n_arr;
    time_CPAS = time_CPAS_arr; time_PEN = time_PEN_arr;
    lambda_CPAS = lambda_CPAS_arr; lambda_PEN = lambda_PEN_arr;
    vol_CPAS = vol_CPAS_arr; vol_PEN = vol_PEN_arr;
    compl_CPAS = compl_CPAS_arr; compl_PEN = compl_PEN_arr;
    x_diff = x_diff_arr;
    save(problem,"m","n","time_CPAS","time_PEN","lambda_CPAS","lambda_PEN","vol_CPAS","vol_PEN","compl_CPAS","compl_PEN","x_diff","x_CPAS_cell","x_PEN_cell")
end


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



% cpas help functions
function cuts = compliance_lmi_cut(x,data) % uses all eigenvectors with negative eigenval
    fprintf('2')
    f = data.f;
    nvar = length(x);
    [allEigVec,allEigVal] = eig(data.gamma*data.K(x)-f*f');
    fprintf('3')
    mineig = min(diag(allEigVal));
    negIdx = find(diag(allEigVal) < min(-data.TOL,mineig/2));
    fprintf('4')
    nneigenval = length(negIdx);
    if nneigenval == 0
        cuts = [];
        return
    end
    w = allEigVec(:,negIdx);
    fprintf('5')
    cuts = zeros(nneigenval, nvar+1);
    for i = 1:nneigenval
        eigenvec = w(:,i);
        for j = 1:nvar-1
            cuts(i,j) = - data.gamma*eigenvec'*data.Kdx{j}*eigenvec;
        end
        %cuts(i,nvar) = 0; % cuts(nvar)*vol
        cuts(i,nvar+1) = data.gamma*eigenvec'*data.K0*eigenvec - (f'*eigenvec)^2;
    end
    fprintf('6')
end

function cuts = frequency_bmi_cut(x,data) % uses all eigenvectors with negative eigenval
    fprintf('7')
    nvar = length(x);
    lambda = x(nvar);
    [allEigVec,allEigVal] = eig(data.K(x) - lambda*data.M(x));
    fprintf('8')
    negIdx = find(diag(allEigVal) < -data.TOL);
    fprintf('9')
    nneigenval = length(negIdx);
    if nneigenval == 0
        cuts = [];
        return
    end
    w = allEigVec(:,negIdx);
    fprintf('A')
    cuts = zeros(nneigenval, nvar+1);
    beta = zeros(nvar-1,1);
    for i = 1:nneigenval
        eigenvec = w(:,i);
        beta0 = eigenvec'*data.M0*eigenvec;
        for j = 1:nvar-1
            beta(j) =  eigenvec'*data.Mdx{j}*eigenvec;
            cuts(i,j) = - (eigenvec'*data.Kdx{j}*eigenvec - lambda*beta(j));
        end
        cuts(i,nvar) = max(beta)*data.vol + beta0; % cuts(nvar)*lambda
        cuts(i,nvar+1) = lambda*cuts(i,nvar) + eigenvec'*data.K0*eigenvec - lambda*beta0;
    end
    fprintf('B')
end

function cuts = compliance_lmi_cut1(x,data) % uses the eignvector of smallest eigenval
    f = data.f;
    nvar = length(x);
    [allEigVec,allEigVal] = eig(data.gamma*data.K(x)-f*f');
    [sortedAllEig,sortedIndex] = sort(diag(allEigVal));
    %[w,smalleig] = eigs(data.gamma*data.K(x)-f*f',1,'smallestreal');
    smalleig = sortedAllEig(1);
    if smalleig >= -data.TOL
        cuts = [];
        return
    end
    w = allEigVec(:,sortedIndex(1));
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

function cuts = frequency_bmi_cut1(x,data)  % uses the eignvector of smallest eigenval
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
