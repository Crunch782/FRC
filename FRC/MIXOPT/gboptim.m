% [ JDJ Xopt neval ]= gboptim(func,Xn,dir,method,C,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% func   : Function to optimize
% Jopt   : Value
% Xn     : Vector (can't be a matrix)
% dir    : Direction: 1: max, -1: min
% method : Method: 'grad' , 'conjgrad', 'powit'
% C      : Norm of the hsphere of constraints (0 = no hsphere constraint)
% p      : p-norm constraint.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% JDJ       :   1st column: history of J,
%               2nd column, history of dJ,
%               3rd column, history of dJp,
%               4th column, history of angle between X and dJ.
%               5th column, history of ndir.
% Jopt = JDJ(end,1)
% Xopt    : Optimal vector
% neval   : number of function evaluations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Gradient Based Optimization Handling Hypersphere Constraints.
% Dimitry Foures - 26/01/12
% MAJ            - 25/04/12   (debug)
% MAJ       	 - 23/05/12   (debug)
% MAJ            - 04/06/12   (p-norm constraints, no more projection of the previous gradient at each step)
% MAJ            - 04/10/12   (Infinity norm constraint when p is very large)
% MAJ            - 08/11/12   (Cleaning)
% MAJ            - 09/01/13   (History of convergence in file + debug)
% MAJ            - 04/02/13   (debug)
% MAJ            - 18/02/13   (Possibility of stop signal with global variable STOP_SIGNAL)
% MAJ            - 06/03/13   (No more 0 gradient: keep old one when not computed)
% MAJ            - 20/03/13   (History of convergence on the fly)
% MAJ            - 28/03/13   (debug)
% MAJ            - 02/04/13   Added angle information
% MAJ            - 01/05/13   Added flag possibility and cleaning

function [JDJ, Xopt, neval, conv]= gboptim(func,Xn,dir,method,C,p,tag)

    %%%%%%%%%%% Parameters %%%%%%%%%%%
    NN          =   100;                % Max number of iterations
    normeps     =   1e-2;              % Normalized residual tolerance
    tol         =   eps^2;           % Gradient tolerance ( hard to reach computer precision for hypersphere constraints or fluid problems... )
    K           =   2;              % Gradient line search multiplication factor
    r           =   0.5;            % Multiplication factor when evaluation fails (step too large)
    np          =   1;              % Choose 1 if an extra parameter in func chooses to compute the gradient or not
    LS          =   1;               % Line search or not
    LSI         =   1;              % Line search interpolation
    e0init      =   1e-1;           % First step choice (angle if hsphere constraint)
    proj        =   1;               % Projection of the gradient for the conjgrad direction (seems better with a projected gradient)
    cmethod = 'rotation';           % Type on constraining algorithm  (lagrangian or rotation)
    test        =   0;              % Test with 2D function.
    conv        = 0;                % Converged or not
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Print in file
    filename = append('Data/',tag,'/JDJHistory.txt');
    fopen(filename,'w+');

    % Check that C~=0 if powit
    if strcmp(method,'powit')
        if C == 0
            fprintf('Cannot use power iteration when C = 0, using gradient descent instead.\n');
            method = 'grad';
        end
    end

    % If lagrangian method: gradient and no line search
    if strcmp(cmethod,'lagrangian')
        if C ~= 0
            if strcmp(method,'conjgrad')
                fprintf('Cannot use conjugate gradient when Lagrangian constraining is specified, using gradient descent instead.\n');
                method = 'grad';
            end
            LSI = 0;
        end
    end

    % Check if X = vector
    if min(size(Xn))~=1
        fprintf('X must be a vector !\n')
        return; 
    end

    % Transform vector in the right format
    if size(Xn,2)~=1
        Xn = Xn'; 
    end

    % Check
    if C~=0 && (norm(Xn,p)-C) > 1e-10
        fprintf('Unconsistent normalization: norm(X) = %f, C = %f...\n',norm(Xn),C); 
        return;
    end;

    % Build function
    funcXn = strcat(func,'(Xn)');

    % 0 = initial values
    % old = previous iteration
    % n = current value
    % c = current best value

    % Initialization
    JDJ = [];
    ndir = 0;

    % First evaluation
    n = 0;
    cp = 1;

    [J0, dJ0] = eval(funcXn);
    dJ0p = proj_grad(Xn,dJ0,e0init,C,p,cmethod);
    n = n+1;
    JDJ(n,:) = write_history(J0,dJ0,dJ0p,e0init,ndir,filename);

    dJold = dJ0;
    dJoldp = dJ0p;
    Jold = J0;
    Xold = Xn;
    L = dir*dJold;

    if strcmp(method,'powit') == 0
    
        % First step choice
        e0 = e0init;
    
        % Begin loop
        restart_same_dir = 0;
        res = 1;
        resNorm = 1;
        while((res > tol)&&(n<=NN)&&(resNorm > normeps))
            % Direction counter
            if restart_same_dir == 0
                ndir = ndir +1; 
            end
            restart_same_dir = 0;
        
            % First step size
            e = e0;
            % First descent on L
            fprintf('Fixed gradient descent/ascent ...\n');
            % Current value
            Jc = Jold;
            dJc = dJold;
            dJcp = dJoldp;
            Xc = Xold;
            % Update
            Xn = update_pos(Xc,L,e,C,p,cmethod);
        
        
            % Evaluation
            [Jn, dJn] = eval(funcXn);          
            if cp == 0
                dJn = dJold; 
            end
            dJnp = proj_grad(Xn,dJn,e,C,p,cmethod);
            n = n+1;
            JDJ(n,:) = write_history(Jn,dJn,dJnp,e,ndir,filename);
            
            % Initialization of line search
            if LS == 0
                cp = 1; 
            else
                cp = 0;   
            end
            
            nl = 1; s = 0; Js = 0;
            % Value storage
            s(nl) = e;
            Js(nl) = Jn;
        
            % Line search: search step;
            while dir*(Jn-Jc) > 0
            
                % Linesearch counter
                nl = nl + 1;
                % Update of best position and value
                Xc = Xn;
                Jc = Jn;
                dJc = dJn;
                % If no line search at all, one update only.
                if LS == 0; break; end;
                % Increase step size if success.
                e = K * e;
                % Update
                Xn = update_pos(Xc,L,e,C,p,cmethod);
                % Evaluation
                [Jn, dJn] = eval(funcXn);
                if cp == 0
                    dJn = dJold; 
                end
                dJnp = proj_grad(Xn,dJn,e,C,p,cmethod);
                n = n+1;
                JDJ(n,:) = write_history(Jn,dJn,dJnp,e,ndir,filename);
                % Save the J(s) with s = sum(e)
                s(nl) = s(nl-1)+e;
                Js(nl) = Jn;
                       
            end
        
        
            % Gradient computation necessary now
            cp =1;
        
            % If no step done successfully, go back to previous candidate with a shorter step.
            if nl == 1
                e0 = r * e0;
                Xc = Xold;
                dJc = dJold;
                Jc = Jold;
                restart_same_dir = 1;
                fprintf('First evaluation of line search failed...\nRecuding step size... (e0 = %5.5e)\n',e0);
                if e0 < eps
                    fprintf('\nWe have reached gradient precision for the step size e...\nThe gradient might be wrong.\n'); 
                    break;
                end;
            end
        
            % No need to interpolate, the best is the second point.
            if nl == 2
                Xn = Xc;
                if LS == 1
                    [Jc , dJc] = eval(funcXn);
                    dJcp = proj_grad(Xc,dJc,e,C,p,cmethod);
                    n = n+1; 
                    JDJ(n,:) = write_history(Jc,dJc,dJcp,e,ndir,filename);
                end
                if cp == 0
                    dJc = dJold; 
                end
                dJcp = proj_grad(Xc,dJc,e,C,p,cmethod);       
                n = n+1;
                JDJ(n,:) = write_history(Jc,dJc,dJcp,e,ndir,filename);
            end
        
            % Line search: interpolation step
            if nl > 2
                if LSI == 1
                    fprintf('Line search from e = %5.5e to %5.5e ...\n',min(s),max(s));
                    si = linspace(s(1),s(end),100);
                    coefs = polyfit(s,Js,numel(s)-1); 
                    Jsi = 0*si;
                    for nc = 1:numel(coefs)
                        Jsi = Jsi + coefs(nc)*si.^(numel(coefs)-nc);
                    end
                
                    if dir == -1 
                        [Jext ie] = min(Jsi); 
                    else
                        [Jext ie] = max(Jsi); 
                    end
                    e = si(ie) - s(numel(s)-1);
                
                    % Update
                    Xn = update_pos(Xc,L,e,C,p,cmethod);
                    % Test
                    [Jn , dJn] = eval(funcXn);
                    if cp == 0
                        dJn = dJold; 
                    end
                    dJnp = proj_grad(Xn,dJn,e,C,p,cmethod);
                    n = n+1;
                    JDJ(n,:) = write_history(Jn,dJn,dJnp,e,ndir,filename);
                    if -dir*(Jn-Jc) > 0
                        fprintf('Line search interpolation unsuccessful !\n');
                        Xn = Xc;
                        [Jc , dJc] = eval(funcXn);
                        %%%missing cp=0
                        dJcp = proj_grad(Xc,dJc,e,C,p,cmethod);
                        n = n+1;
                        JDJ(n,:) = write_history(Jc,dJc,dJcp,e,ndir,filename);
                    else
                        fprintf('Line search interpolation successful !\n');
                        Xc = Xn;
                        Jc = Jn;
                        dJc = dJn;
                    end
                else
                    Xn = Xc;
                    [Jc , dJc] = eval(funcXn);
                    if cp == 0
                        dJc = dJold; 
                    end
                    dJcp = proj_grad(Xc,dJc,e,C,p,cmethod);
                    n = n+1;
                    JDJ(n,:) = write_history(Jc,dJc,dJcp,e,ndir,filename);
                end
            end
                
            % Projection of the new gradient
            dJcp = proj_grad(Xc,dJc,e,C,p,cmethod);
        
            % Conjugate gradient
            if strcmp(method,'grad') == 1
                % Steepest descent / ascent
                b = 0;
            elseif  strcmp(method,'conjgrad') == 1
                % Polack-Ribiere
                b = dJcp'*(dJcp - dJoldp)./(dJoldp'*dJoldp);
            end
            b = max(0,b);
        
            % Update direction
            if proj == 1
                L = dir*dJcp + b*L;
            else
                L = dir*dJc + b*L;
            end
        
            % New e0 based on gradient value.
            e0 = min(e0init,norm(dJoldp)/norm(dJcp)*e0);
            if e0 < eps
                fprintf('\nWe have reached gradient precision for the step size e...\nThe gradient might be wrong.\n'); 
                break; 
            end
            
            if LS == 0
                e0 = K * e0; 
            end
            fprintf('New step size: e0 = %5.5e\n',e0);
        
        
            dJcp = proj_grad(Xc,dJc,e,C,p,cmethod);
        
        
            [res, resNorm] = display_infos(Jc,JDJ,n);
        
        
            % Update all
            dJoldp = dJcp;
            dJold = dJc;
            Jold = Jc;
            Xold = Xc;
        
        
        end
    
    elseif strcmp(method,'powit')==1
    
        % Begin loop
        n = 1;
        res = 1;
        resNorm = 1;
        dres = 1;
        while((res > tol)&&(dres > eps)&&(n<=NN)&&(resNorm > normeps))
        
       
            Xn = L*C/norm(L);
        
            [Jn, dJn] = eval(funcXn);
            dJnp = proj_grad(Xn,dJn,1,C,p,'rotation');
            n = n+1;
            JDJ(n,:) = write_history(Jn,dJn,dJnp,1,n,filename);
        
        
            res_old = res;
            [res, resNorm] = display_infos(Jn,JDJ,n);
            dres = abs(res_old - res);
        
            L = dir*dJn;
            Jold = Jn;
        
         
        end
    
    
    end

    % Optimal position
    Xopt  = Xold ;
    neval = n;
    if resNorm > normeps
        conv = 0;
    else
        conv = 1;
    end



end





function Xn = update_pos(Xold,L,e,C,p,cmethod)
    L = L/norm(L);
    if C ~= 0
        % Projection
        L = proj_grad(Xold,L,e,C,p,cmethod);
        if strcmp(cmethod,'rotation')
            Xn = cos(e)*Xold + sin(e)*C*L;
            Xn = C*Xn./norm(Xn,p);  % Normalization if not shpere.
        elseif strcmp(cmethod,'lagrangian')
            Xn = Xold + e * L;
            Xn = C*Xn./norm(Xn);           % Normalization or deviation....
        end
        if abs(norm(Xn) - C) > 1e-12; fprintf('Nornalization problem: |norm(Xn)-C| = %5.5e\n',abs(norm(Xn)-C)); 
            pause; 
        end
    else
        Xn = Xold + e * L;
    end
end

function dJn = proj_grad(Xold,dJn,e,C,p,cmethod)
    if C~= 0
        if strcmp(cmethod,'rotation')
            Ghold = abs(Xold).^(p-1).*sign(Xold);
            if isinf(norm(Ghold)) || norm(Ghold)<1e-14   % if p is too large: max norm maximized
                [~ , i] = max(abs(Xold));
                Ghold = zeros(numel(Xold),1); Ghold(i) = sign(Xold(i));
            end
            dJn = dJn - Ghold'*dJn./(Ghold'*Ghold)*Ghold;
        elseif strcmp(cmethod,'lagrangian')
            dJn = dJn/norm(dJn);
            % Lambda equation...
            lambdam = 1i;
            while isreal(lambdam) == 0
                lambdap = 1/e + norm(dJn)/C*dJn'*Xold/(norm(dJn)*norm(Xold)) + (1/e^2 - (norm(dJn)/C)^2*(1-(dJn'*Xold/(norm(dJn)*norm(Xold)))^2))^(1/2);
                lambdam = 1/e + norm(dJn)/C*dJn'*Xold/(norm(dJn)*norm(Xold)) - (1/e^2 - (norm(dJn)/C)^2*(1-(dJn'*Xold/(norm(dJn)*norm(Xold)))^2))^(1/2);
                e = e./10;
                if e < eps
                    fprintf('Warning: close to computer precision for e. (e = %5.4e).\n',e); 
                end
            end
            dJn = dJn-lambdam*Xold;
        end
    end
end



function JDJ = write_history(J,dJ,dJp,e,ndir,filename)

    ndJ = norm(dJ);
    ndJp = norm(dJp);
    a =  norm(dJp)/norm(dJ);
    e = abs(e);
    JDJ = [J, ndJ, ndJp, a, e, ndir];
    save(filename,'JDJ','-ASCII','-append');

end

function [res, resNorm] = display_infos(Jc,JDJ,n)
        
    J = JDJ(n,1);
    rdJ = JDJ(n,2);
    rdJp = JDJ(n,3);
    a = JDJ(n,4);  %alpha ? la sortie ?cran
    e = JDJ(n,5);
    ndir = JDJ(n,6);
    
    res = rdJp.^2;
    resNorm = a.^2;

    % Display
    fprintf('\n%dth Direction, %dth Function evaluation, Current best: Jc = %5.10e,\n\t J = %5.5e,\n\t dJ = %5.5e,\n\t dJp = %5.5e,\n\t alpha = %5.5e\n\t Step size = %5.5e\n\n',ndir,n,Jc,J,rdJ,rdJp,a,e);

end
        