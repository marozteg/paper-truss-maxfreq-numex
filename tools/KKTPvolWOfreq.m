function KKTPvolWOfreq(x,lambC,lambX,A,Adx,gamma,f)
% The first order KKT of min sum(x) s.t. [gamma -f;-f A(x)] pos semidef is
%       ones(n,1) + gradC*svec3(lambC) + eye(n)*lambX == 0,
%       [gamma -f'';-f A(x))]*lambC == 0,
%       [gamma -f'';-f A(x))] positive semidefinite,
%       x>=0,
%       lambC negative semidefinite and
%       lambX negative semidefinite

    n = length(x);
    m = size(f,1);
    z = zeros(m,1);
    gradC = zeros(n,(m+1)*(m+2)/2);
    for i = 1:n
        gradC(i,:) = svec3([0 z';z Adx{i}])';
    end
    gradlag = ones(n,1) + gradC*svec3(lambC) + eye(n)*lambX;
    fprintf('*******************************************\n')
    fprintf('First order KKT conditions for (PvolWOfreq)\n')
    fprintf('\tgradient of lagrangian\n')
    fprintf('\t\tnorm(gradlag)\t\t\t\t\t\t\t%1.1e\n',norm(gradlag))
    fprintf('\tcomplementarity\n')
    fprintf('\t\tnorm([gamma -f'';-f A(x))]*lambC)\t\t%1.1e\n',norm([gamma -f';-f A(x)]*lambC))
    fprintf('\t\tx''*lambX\t\t\t\t\t\t\t\t%1.1e\n',x'*lambX)
    fprintf('\tfeasibility\n')
    fprintf('\t\tCompl constraint\n')
    mineig1 = mineig([gamma -f';-f A(x)]);
    mc = mincompl(A(x),f);
    fprintf('\t\t\t[gamma -f'';-f A(x)]>=0\t\t\t\t%s\n',isTrue(mineig1>=0))
    fprintf('\t\t\tmineig([gamma -f'';-f A(x)])\t\t\t%1.1e\n',mineig1)
    fprintf('\t\t\tmin{c: [c -f'';-f A(x)]>=0}<=gamma\t%s\n',isTrue(mc<=gamma))
    fprintf('\t\t\tmin{c: [c -f'';-f A(x)]>=0}=%1.1e\t%1.1e=gamma\n',mc,gamma)
    fprintf('\t\tTopo constraint\n')
    fprintf('\t\t\tall(x>=0)\t\t\t\t\t\t\t%s\n',isTrue(all(x(1:n)>=0)))
    fprintf('\t\t\tmin(x)\t\t\t\t\t\t\t\t%1.1e\n',min(x(1:n)))
    fprintf('\t\tLagrange multipliers\n')
    fprintf('\t\t\tmaxeig(lambC)\t\t\t\t\t\t%1.1e\n',max(eig(lambC)))
    fprintf('\t\t\tmax(lambX)\t\t\t\t\t\t\t%1.1e\n',max(lambX))
end

function s = isTrue(b)
    if b == 1
        s = 'true';
    else
        s = 'false';
    end
end