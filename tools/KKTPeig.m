function KKTPeig(x,lambC,lambF,lambV,lambX,A,B,Adx,Bdx,gamma,f,V)
    n = length(x);
    m = size(lambF,1);
    z = zeros(m,1);
    gradC = zeros(n,(m+1)*(m+2)/2);
    gradF = zeros(n,m*(m+1)/2);
    for i = 1:n-1
        gradC(i,:) = svec3([0 z';z Adx{i}])';
        gradF(i,:) = svec3(Adx{i} - x(n)*Bdx{i})';
    end
    gradC(n,:) = svec3(zeros(m+1,m+1))';
    gradF(n,:) = svec3(-B(x(1:n-1)))';
    gradlag = [zeros(n-1,1);-1] + gradC*svec3(lambC) + gradF*svec3(lambF) + eye(n)*lambX + [ones(n-1,1);0]*lambV;
    fprintf('*************************************\n')
    fprintf('First order KKT conditions for (Peig)\n')
    fprintf('\tgradient of lagrangian\n')
    fprintf('\t\tnorm(gradlag)\t\t\t\t\t\t\t%1.1e\n',norm(gradlag))
    fprintf('\tcomplementarity\n')
    fprintf('\t\tnorm([gamma -f'';-f A(x))]*lambC)\t\t%1.1e\n',norm([gamma -f';-f A(x(1:n-1))]*lambC))
    fprintf('\t\tnorm((A(x)-lambda*B(x))*lambF)\t\t\t%1.1e\n',norm((A(x(1:n-1))-x(n)*B(x(1:n-1)))*lambF))
    fprintf('\t\tx''*lambX\t\t\t\t\t\t\t\t%1.1e\n',x'*lambX)
    fprintf('\t\t(sum(x(1:n))-V)*lambV\t\t\t\t\t%1.1e\n',(sum(x(1:n-1))-V)*lambV)
    fprintf('\tfeasibility\n')
    fprintf('\t\tCompl constraint\n')
    mineig1 = mineig([gamma -f';-f A(x(1:n-1))]);
    mc = mincompl(A(x(1:n-1)),f);
    fprintf('\t\t\t[gamma -f'';-f A(x)]>=0\t\t\t\t%s\n',isTrue(mineig1>=0))
    fprintf('\t\t\tmineig([gamma -f'';-f A(x)])\t\t\t%1.1e\n',mineig1)
    fprintf('\t\t\tmin{c: [c -f'';-f A(x)]>=0}<=gamma\t%s\n',isTrue(mc<=gamma))
    fprintf('\t\t\tmin{c: [c -f'';-f A(x)]>=0}=%1.1e\t%1.1e=gamma\n',mc,gamma)    
    fprintf('\t\tTopo constraint\n')
    fprintf('\t\t\tall(x>=0)\t\t\t\t\t\t\t%s\n',isTrue(all(x(1:n-1)>=0)))
    fprintf('\t\t\tmin(x)\t\t\t\t\t\t\t\t%1.1e\n',min(x(1:n-1)))
    fprintf('\t\tVol constraint\n')
    fprintf('\t\t\tsum(x)<=V\t\t\t\t\t\t\t%s\n',isTrue(sum(x(1:n-1))<=V))
    fprintf('\t\t\tsum(x)=%1.1e\t\t\t\t\t\t%1.1e=V\n',sum(x(1:n-1)),V)
end

function s = isTrue(b)
    if b == 1
        s = 'true';
    else
        s = 'false';
    end
end