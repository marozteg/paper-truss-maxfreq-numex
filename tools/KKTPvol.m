function KKTPvol(x,lambC,lambF,lambX,A,B,Adx,Bdx,gamma,f,lambda)
    n = length(x);
    m = size(lambF,1);
    z = zeros(m,1);
    gradC = zeros(n,(m+1)*(m+2)/2);
    gradF = zeros(n,m*(m+1)/2);
    for i = 1:n
        gradC(i,:) = svec3([0 z';z Adx{i}])';
        gradF(i,:) = svec3(Adx{i} - lambda*Bdx{i})';
    end
    gradlag = ones(n,1) + gradC*svec3(lambC) + gradF*svec3(lambF) + eye(n)*lambX;
    fprintf('*************************************\n')
    fprintf('First order KKT conditions for (Pvol)\n')
    fprintf('\tgradient of lagrangian\n')
    fprintf('\t\tnorm(gradlag)\t\t\t\t\t\t\t%1.1e\n',norm(gradlag))
    fprintf('\tcomplementarity\n')
    fprintf('\t\tnorm([gamma -f'';-f A(x))]*lambC)\t\t%1.1e\n',norm([gamma -f';-f A(x)]*lambC))
    fprintf('\t\tnorm((A(x)-lambda*B(x))*lambF)\t\t\t%1.1e\n',norm((A(x)-lambda*B(x))*lambF))
    fprintf('\t\tx''*lambX\t\t\t\t\t\t\t\t%1.1e\n',x'*lambX)
    fprintf('\t\tmaxeig(lambC)\t\t\t\t\t\t\t%1.1e\n',max(eig(lambC)))
    fprintf('\t\tmaxeig(lambF)\t\t\t\t\t\t\t%1.1e\n',max(eig(lambF)))
    fprintf('\t\tmax(lambX)\t\t\t\t\t\t\t\t%1.1e\n',max(lambX))
    fprintf('\tfeasibility\n')
    fprintf('\t\tCompl constraint\n')
    mineig1 = mineig([gamma -f';-f A(x)]);
    mc = mincompl(A(x),f);
    fprintf('\t\t\t[gamma -f'';-f A(x)]>=0\t\t\t\t%s\n',isTrue(mineig1>=0))
    fprintf('\t\t\tmineig([gamma -f'';-f A(x)])\t\t\t%1.1e\n',mineig1)
    fprintf('\t\t\tmin{c: [c -f'';-f A(x)]>=0}<=gamma\t%s\n',isTrue(mc<=gamma))
    fprintf('\t\t\tmin{c: [c -f'';-f A(x)]>=0}=%1.1e\t%1.1e=gamma\n',mc,gamma)
    fprintf('\t\tFreq constraint\n')
    mineig2 = mineig(A(x) - lambda*B(x));
    mp = maxpsdpair(A(x),B(x));
    fprintf('\t\t\tA(x)-lambda*B(x)>=0\t\t\t\t\t%s\n',isTrue(mineig2 >= 0))
    fprintf('\t\t\tmineig(A(x)-lambda*B(x))\t\t\t%1.1e\n',mineig2)
    fprintf('\t\t\tlambda <= max{c: A(x)-c*B(x)>=0}\t%s\n',isTrue(lambda<=mp))
    fprintf('\t\t\tmax{c: A(x)-c*B(x)}=%1.1e\t\t\t%1.1e=lambda\n',mp,lambda)
    fprintf('\t\tTopo constraint\n')
    fprintf('\t\t\tall(x>=0)\t\t\t\t\t\t\t%s\n',isTrue(all(x(1:n)>=0)))
    fprintf('\t\t\tmin(x)\t\t\t\t\t\t\t\t%1.1e\n',min(x(1:n)))
end

function s = isTrue(b)
    if b == 1
        s = 'true';
    else
        s = 'false';
    end
end