function gamma = mincompl(A,f,epsilon,verbose)
% given A positive semidefinite and f non null
% compute gamma = argmin{gamma: [gamma -f^T; -f A] positive semidefinite}

    P = orth(full(A));
    [V,D] = eig(P'*(f*f')*P,P'*A*P);
    d = sort(diag(D),"descend");
    gamma = d(1);


    % posNear0 = 1e-7;
    % if nargin == 2
    %     epsilon = 1e-6;
    %     verbose = 0;
    % elseif nargin == 3
    %     verbose = 0;
    % end
    % smalleigA = eigs(A,1,'smallestreal');
    % if smalleigA < -posNear0 || norm(f) < posNear0
    %     keyboard
    % end
    % 
    % % looking for gamma such that [gamma -f'; -f A] not positive semidef.
    % gamma0 = -10;
    % maxtry = 100;
    % gamma = gamma0;
    % for i = 1:maxtry
    %     [w,feas] = eigs([gamma -f';-f A],1,'smallestreal');
    %     h = w(1);
    %     w = w(2:end);
    %     found_negative_eig = feas < -epsilon;
    %     if found_negative_eig
    %         break
    %     end
    %     gamma = 10*gamma;
    % end
    % if ~found_negative_eig || isnan(feas)
    %     keyboard % unbound solution? or have to increase maxtry? 
    % end
    % if verbose
    %     fprintf('k=1\tgamma=%1.2e\tfeas=%1.2e\n',gamma,feas)
    % end
    % k = 2;
    % gamma_prev = gamma;
    % count = 0;
    % while 1
    %     if abs(h)<posNear0 || isnan(feas)
    %         keyboard
    %     end
    %     gamma = (2*h*w'*f-w'*A*w)/(h*h);
    %     [w,feas] = eigs([gamma -f';-f A],1,'smallestreal');
    %     h = w(1);
    %     w = w(2:end);
    %     if verbose
    %         fprintf('k=%d\tgamma=%1.2e\tfeas=%1.2e\n',k,gamma,feas)
    %     end
    %     if feas > -epsilon
    %         break
    %     end
    %     if abs(gamma-gamma_prev)/(1+gamma) < epsilon
    %         count = count + 1;
    %         if count == 5
    %             fprintf('mincompl: gamma did''t change in %d iteration.\n',count)
    %             break
    %         end
    %     end
    %     gamma_prev = gamma;
    %     k = k + 1;
    % end
    % 
    % if verbose
    %     % some check
    %     delta = 1e-4;
    %     [~,poseig] = eigs([gamma+delta -f';-f A],1,'smallestreal');
    %     [~,negeig] = eigs([gamma-delta -f';-f A],1,'smallestreal');
    %     if poseig<-posNear0 || negeig>posNear0
    %         keyboard
    %     else % seems good
    %         gap = 2*delta;
    %         for gamma2 = gamma-delta:gap/10:gamma+delta
    %             [~,negeig] = eigs([gamma2 -f';-f A],1,'smallestreal');
    %             if negeig>posNear0
    %                 if abs(gamma-gamma2)>1e-4
    %                     keyboard
    %                 else
    %                     break
    %                 end
    %             end
    %         end
    %     end
    % end
end