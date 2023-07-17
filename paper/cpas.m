function [x,A,b,exitflag] = cpas(c,A,b,lb,ub,data,verbose)
    nfix = size(A,1);
    cplex = Cplex('');
    cplex.Model.sense = 'minimize';
    cplex.addCols(c,[],lb,ub);
    cplex.addRows(-Inf*ones(length(b),1), A, b);
    cplex.DisplayFunc = @(s) ''; % suppress display output
    cplex.solve();
    x = cplex.Solution.x; % solves initial relaxation (x = argmin{c'*x : A*x<=b, lb<=x<=ub})
    n = length(x);
    iter = 1;
    if nargin == 6
        verbose = 0;
    end
    elapsedtimeid = tic;
    itertimeid = tic;
    while 1
        fprintf('1')
        cuts = [];
        for k = 1:data.ncut
            cuts = [cuts; data.cut{k}(x,data)];
        end
        fprintf('C')
        if isempty(cuts)
            exitflag = 1;
            elapsedtime = toc(elapsedtimeid);
            break
        end
        aa = cuts(:,1:n);
        bb = cuts(:,n+1);
        sa = sqrt(diag(aa*aa'));
        gap = min(abs(aa*x - bb)./sa);
        itertime = toc(itertimeid);
        fprintf('D')
        if verbose && mod(iter,1)==0
            fprintf('k=%d fobj=%1.9e gap=%1.9e ncuts=%d itertime=%1.2e\n',iter,c'*x,gap,size(cuts,1),itertime)
        end
        A = [A;cuts(:,1:n)]; b = [b;cuts(:,n+1)];
        fprintf('E')
        cplex.addRows(-Inf*ones(length(cuts(:,n+1)),1), cuts(:,1:n), cuts(:,n+1));
        cplex.Start.x = x;
        fprintf('F')
        cplex.solve();
        fprintf('G')
        x = cplex.Solution.x;
        if mod(iter,1)==0 && size(A,1)>900
            fprintf('fobj=%1.9e\n',c'*x)
            s = b-A*x;
            s(1:nfix) = -1;
            [~,idx] = sort(s,'descend');
            numdel = size(A,1) - 600;
            whichdel = idx(1:numdel);
            cplex.delRows(whichdel);
            %cplex.addRows(-Inf,[ones(1,n-1) 0],data.vol); % <--- vol constraint
            A(whichdel,:) = [];
            %A = [A;[ones(1,n-1) 0]]; % <--- vol constraint
            b(whichdel) = [];
            %b = [data.vol;b]; % <--- vol constraint
            cplex.solve();
            x = cplex.Solution.x;
            fprintf('fobj=%1.9e\n',c'*x)
        end
        exitflag = cplex.Solution.status;
        if exitflag > 1 && exitflag~=5
            fprintf('cplex.Solution.status %d.\n',exitflag)
            fprintf('%s.\n',cplex.Solution.statusstring)
            fprintf('Method = %d.\n',cplex.Solution.method)
            keyboard
        end
        iter = iter + 1;
        elapsedtime = toc(elapsedtimeid);
        if elapsedtime > 24*3600%60*5
            fprintf('time exceeded... aborting\n')
            exitflag = -999;
            break
        end
    end
    fprintf('Elapsed time: %d seconds\n',elapsedtime)
end