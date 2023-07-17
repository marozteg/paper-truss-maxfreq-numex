function lambda = maxpsdpair(A,B,epsilon,verbose)
% given A and B positive semidefinite
% compute lambda = argmax{c: A-c*B is positive semidefinite}

    P = orth(full(B));
    [V,D] = eig(P'*A*P,P'*B*P);
    d = sort(diag(D));
    lambda = d(1);

    % posNear0 = 1e-7;
    % if nargin == 2
    %     epsilon = 1e-6;
    %     verbose = 0;
    % elseif nargin == 3
    %     verbose = 0;
    % end
    % smalleigA = eigs(A,1,'smallestreal');
    % smalleigB = eigs(B,1,'smallestreal');
    % if smalleigA < -posNear0 || smalleigB < -posNear0
    %     keyboard
    % end
    % 
    % % looking for lambda such that (A-lambda*B) not positive semidef.
    % lambda0 = 10;
    % maxtry = 100;
    % lambda = lambda0;
    % for i = 1:maxtry
    %     [w,feas] = eigs(A-lambda*B,1,'smallestreal');
    %     found_negative_eig = feas < -epsilon;
    %     if found_negative_eig
    %         break
    %     end
    %     lambda = 10*lambda;
    % end
    % if ~found_negative_eig || isnan(feas)
    %     keyboard % unbound solution? or have to increase maxtry? 
    % end
    % if verbose
    %     fprintf('k=1\tlambda=%1.2e\tfeas=%1.2e\n',lambda,feas)
    % end
    % k = 2;
    % lamb_prev = lambda;
    % count = 0;
    % while 1
    %     if norm(B*w) < posNear0 || isnan(feas)
    %         keyboard
    %     end
    %     lambda = w'*A*w/(w'*B*w);
    %     [w,feas] = eigs(A-lambda*B,1,'smallestreal');
    %     if verbose
    %         fprintf('k=%d\tlambda=%1.2e\tfeas=%1.2e\n',k,lambda,feas)
    %     end
    %     if feas > -epsilon
    %         break
    %     end
    %     if abs(lambda-lamb_prev)/(1+lambda) < epsilon
    %         count = count + 1;
    %         if count == 5
    %             fprintf('maxpsdpair: lambda did''t change in %d iteration.\n',count)
    %             break
    %         end
    %     end
    %     lamb_prev = lambda;
    %     k = k + 1;
    % end
    % 
    % if verbose
    %     % some check
    %     delta = 1e-4;
    %     [~,negeig] = eigs(A-(lambda+delta)*B,1,'smallestreal');
    %     [~,poseig] = eigs(A-(lambda-delta)*B,1,'smallestreal');
    %     if poseig<-posNear0 || negeig>posNear0
    %         keyboard
    %     else % seems good
    %         gap = 2*delta;
    %         for lambda2 = lambda-delta:gap/10:lambda+delta
    %             [~,poseig] = eigs(A-lambda2*B,1,'smallestreal');
    %             if poseig<-posNear0
    %                 if abs(lambda-lambda2)>1e-4
    %                     keyboard
    %                 else
    %                     break
    %                 end
    %             end
    %         end
    %     end
    % end
end