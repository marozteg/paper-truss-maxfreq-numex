function [x,A,b,exitflag] = cpas(c,A,b,lb,ub,data,verbose)
    [x,~,~,~,~] = linprog(c,A,b,[],[],lb,ub,data.options);
    n = length(x);
    iter = 1;
    if nargin == 6
        verbose = 0;
    end
    while 1
        cuts = [];
        for k = 1:data.ncut
            cuts = [cuts; data.cut{k}(x,data)];
        end
        if isempty(cuts)
            exitflag = 1;
            break
        end
        a2 = sum(cuts(end,1:n)*cuts(end,1:n)');
        gap = min(abs(cuts(end,1:n)*x - cuts(end,n+1))./sqrt(a2));
        if verbose
            fprintf(['k=%d fobj=%1.9e ' ...
                'gap=%1.9e ' ...
                'ncuts=%d\n'],iter,c'*x,gap,size(cuts,1))
        end
        if gap < data.TOL
            exitflag = 2;
            break
        end
        A = [A;cuts(:,1:n)];
        b = [b;cuts(:,n+1)];
        [x,~,exitflag,~,~] = linprog(c,A,b,[],[],lb,ub,data.options);
        if exitflag <= 0
            fprintf('linprog existflag %d.\n',exitflag)
            keyboard
        end
        iter = iter + 1;
    end
end