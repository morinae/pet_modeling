function [x_res, x_r, res, iter, x_nof, y_d, XN, x_kprm] = IRGNM_modified(t, delta_x, Imax, K, m_Biexp, Lam, Mu, C_TOT_meas, C_T_noise, a1, a2, b1, b2, g1, g2, scaleflag)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INPUT                                                            %
    % t ........ temporal grid                (      in seconds     )  %
    % delta_x .. error level                  (  for initial guess  )  %
    % Imax ..... maximal iterations                                    %
    % K ........ metabolic parameters                                  %
    % m_Biexp .. parent plasma fraction parameters (biexpontential)    %
    % Lam ...... factors of parameterization of plasma concentration   %
    % Mu ....... exponents of parameterization of plasma concentration %
    % C_TOT_meas measured samples of whole blood tracer concentration  %
    % C_T_noise  measured samples of tracer concentration in tissue    %
    % a1 ....... regularisation parameter for plasma parameters        %
    % a2 ....... regularisation parameter for plasma parameters        %
    % b1 ....... regularisation parameter for ppf. function parameters %
    % b2 ....... regularisation parameter for ppf. function parameters %
    % g1 ....... regularisation parameter for metabolic parameters     %
    % g2 ....... regularisation parameter for metabolic parameters     %
    % scaleflag  flag if scaling approach is applied to method         %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OUTPUT                                                           %
    % x_res ..... array of relative deviations (parameter types local) %
    % x_r ....... array of relative deviation  (       global        ) %
    % res ....... array of residual            (       global        ) %
    % iter ...... amount of executed iterations                        %
    % x_nof ..... array of relative deviation  (without f parameters ) %
    % y_d ....... noisy measurements                                   %
    % XN ........ history of all iterations                            %
    % x_kprm .... array of relative deviations (metabolic separately ) %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % INITIALISATION %
    T = size(t,2);                  % amount of Time measurements
    t_sec = t/60;                   % time normalisation in minutes
    n=size(K,1);                    % number of considered regions
    k_23 = K(:,2)+K(:,3);           % additive metabolic values k_2+k_3
    q=size(m_Biexp,2);              % number of f parameters
    p = size(Lam,2);                % degree of plasma concentration
    eps = 0.001;                    % parameter domain delimiter
    iter=1;                         % iteration parameter
    scale = 1;                      % scale parameter                     

    x_dagg = [Lam, Mu, m_Biexp, reshape(K',[1,3*n])]'; % ground truth
    K_dagg = K; % ground truth metabolic parameters
    Kcoeff_dagg = K_dagg(:,1).*K_dagg(:,3)./(K_dagg(:,2)+K_dagg(:,3));

    % reshape parameters of plasma concentration
    vMu= reshape(Mu,[1,1,p]);
    vLam=reshape(Lam,[1,1,p]);
    
    % Polyexponential Integral Tensors
    % Int(exp(mu_i s))
    I_1 = (exp(vMu.*t_sec)-1)./vMu;
    
    % Int(exp((k2+k3)(s-t))exp(mu_i s))
    I_2=(exp(vMu.*t_sec)-exp(-k_23.*t_sec))./(vMu+k_23);

    I_1 = I_1.*vLam;      % Overload Tensor I_1
    I_2 = I_2.*vLam;      % Overload Tensor I_2

    I_1 = sum(I_1,3);     % Build up Integral values
    I_2 = sum(I_2,3);     % Build up Integral values

    % Initial Guess of Ground Truth 
    g = randn(size(x_dagg))*sqrt(delta_x/4)+delta_x;   % Gaussian Term
    s = 2*floor(2*rand(size(x_dagg)))-1;               % Bernoulli Term
    x_star = x_dagg.*(1+s.*g);                         % Initial Guess
    
    % Projection to Domain of Operator F
    x_star(p+1:2*p) = min(x_star(p+1:2*p), -eps);
    x_star(2*p+1:2*p+q) = [max(x_star(2*p+1),eps); min(x_star(2*p+2),-eps); min(x_star(2*p+3),-eps)]; 
    x_star(2*p+q+1:end) = max(x_star(2*p+q+1:end),eps);
    
    % Parameters according to Initial Guess
    x_n = x_star;
    Lam = x_n(1:p)';                % factors of plasma concentration
    Mu = x_n(p+1:2*p)';             % exponents of plasma concentration
    m_Biexp = x_n(2*p+1:2*p+q)';    % parameters of function f
    K(:,1)=x_n(2*p+q+1:3:end);      % K_1 metabolic parameters
    K(:,2)=x_n(2*p+q+2:3:end);      % k_2 metabolic parameters
    K(:,3)=x_n(2*p+q+3:3:end);      % k_3 metabolic parameters
    k_23 = K(:,2)+K(:,3);           % Additive values k_2+k_3
    
    % reshape parameters of plasma concentration
    vMu= reshape(Mu,[1,1,p]);
    vLam=reshape(Lam,[1,1,p]);
    
    % Polyexponential Integral Tensors
    % Int(exp(mu_i s))
    I_1 = (exp(vMu.*t_sec)-1)./vMu;

    % Int(exp((k2+k3)(s-t))exp(mu_i s))
    I_2=(exp(vMu.*t_sec)-exp(-k_23.*t_sec))./(vMu+k_23);

    I_1 = I_1.*vLam;      % Overload Tensor I_1
    I_2 = I_2.*vLam;      % Overload Tensor I_2

    I_1 = sum(I_1,3);     % Build up Integral values
    I_2 = sum(I_2,3);     % Build up Integral values

    % Res1 (Restriction to forward model of tracer concentration in tissue)
    Res1 = (K(:,1).*K(:,2)./k_23).*I_2+(K(:,1).*K(:,3)./k_23).*I_1; 
    % Res2 (Restriction to forward model of whole blood concentration)
    Res2 = (Lam*exp(Mu'.*t_sec))./(m_Biexp(1).*exp(m_Biexp(2).*t_sec)+(1-m_Biexp(1)).*exp(m_Biexp(3).*t_sec));
    % Mapped initial guess (First approximation of y^+)
    F = [Res1;Res2];                           
    
    y_d =[C_T_noise; C_TOT_meas]; % Noisy Measurement y_d 
    R = F-y_d;                    % Initial residual element in Image Space
    r = sqrt(trace(R'*R));        % Initial l^2-Residual
    
    % Initial residuals and deviations
    res=[r];                   
    x_tild = x_dagg-x_n;       
    x_res = [sqrt(x_tild(1:2*p)'*x_tild(1:2*p))/(sqrt(x_dagg(1:2*p)'*x_dagg(1:2*p))), sqrt(x_tild(2*p+1:2*p+q)'*x_tild(2*p+1:2*p+q))/(sqrt(x_dagg(2*p+1:2*p+q)'*x_dagg(2*p+1:2*p+q))), sqrt(x_tild(2*p+q+1:end)'*x_tild(2*p+q+1:end))/(sqrt(x_dagg(2*p+q+1:end)'*x_dagg(2*p+q+1:end)))];
    x_r = [sqrt(x_tild'*x_tild)/(sqrt(x_dagg'*x_dagg))];
    x_nof = [sqrt(x_tild([1:2*p,2*p+q+1:end])'*x_tild([1:2*p,2*p+q+1:end]))/sqrt(x_dagg([1:2*p,2*p+q+1:end])'*x_dagg([1:2*p,2*p+q+1:end]))];
    XN = [x_star'];
    K_tild = K_dagg-K;
    x_kprm = [sqrt(sum(K_tild.^2))./sqrt(sum(K_dagg.^2)), sqrt(sum((Kcoeff_dagg-K(:,1).*K(:,3)./k_23).^2))./sqrt(sum(Kcoeff_dagg.^2))];
    
    % ITERATION PROCESS
    while iter<Imax % Discrepancy Principle disabled for objective evaluation                 
        % Regularisation parameters
        alpha = a1*2^(-iter/a2);  % for plasma concentration
        beta = b1*2^(-iter/b2);   % for ppf function f
        gamma = g1*2^(-iter/g2);  % for metabolic parameters  

        % reshape parameters of plasma concentration
        vMu= reshape(Mu,[1,1,p]);
        vLam=reshape(Lam,[1,1,p]);
        % Polyexponential Integral Tensors
        % Int(exp(mu_i s))
        I_1 = (exp(vMu.*t_sec)-1)./vMu;

        % Int(exp((k2+k3)(s-t))exp(mu_i s))
        I_2=(exp(vMu.*t_sec)-exp(-k_23.*t_sec))./(vMu+k_23);
        
        % Int(sexp((k2+k3)(s-t))exp(mu_i s))
        I_3 = t_sec.*exp(vMu.*t_sec)./(vMu+k_23)-(exp(vMu.*t_sec)-exp(-k_23.*t_sec))./((vMu+k_23).^2);

        % Int(sexp(mu_i s))
        I_4 = t_sec.*exp(vMu.*t_sec)./vMu-(exp(vMu.*t_sec)-1)./(vMu.^2);
        
        % Fréchet-Differential w.r.t. lambda parameters
        D_Lam = (K(:,1).*K(:,2)./k_23).*I_2+ (K(:,1).*K(:,3)./k_23).*I_1;
        D_Lam = reshape(permute(D_Lam,[3,2,1]),[p,n*T]);     
        D_Lam = [D_Lam,scale*exp(Mu'.*t_sec)./(m_Biexp(1).*exp(m_Biexp(2).*t_sec)+(1-m_Biexp(1)).*exp(m_Biexp(3).*t_sec))];
        
        % Overload Exponential Integral Tensors
        I_1 = I_1.*vLam;
        I_2 = I_2.*vLam;
        I_3 = I_3.*vLam;
        I_4 = I_4.*vLam;

        % Fréchet-Differential w.r.t. mu parameters
        D_Mu = (K(:,1).*K(:,2)./k_23).*I_3+(K(:,1).*K(:,3)./k_23).*I_4;
        D_Mu = reshape(permute(D_Mu,[3,2,1]),[p,n*T]); 
        D_Mu = [D_Mu, scale*(Lam'.*t_sec.*exp(Mu'.*t_sec))./(m_Biexp(1).*exp(m_Biexp(2).*t_sec)+(1-m_Biexp(1)).*exp(m_Biexp(3).*t_sec));];

        % Fréchet-Differential w.r.t. ppf f parameters
        D_m = [(exp(m_Biexp(2)*t_sec)-exp(m_Biexp(3)*t_sec));
               m_Biexp(1)*t_sec.*exp(m_Biexp(2)*t_sec);
              (1-m_Biexp(1))*t_sec.*exp(m_Biexp(3)*t_sec)];
        D_m = (-scale*(Lam*exp(Mu'.*t_sec))./((m_Biexp(1).*exp(m_Biexp(2).*t_sec)+(1-m_Biexp(1)).*exp(m_Biexp(3).*t_sec)).^2)).*D_m;

        % Build up Integrals w.r.t. plasma concentration
        I_1 = sum(I_1,3);
        I_2 = sum(I_2,3);
        I_3 = sum(I_3,3);

        % Fréchet-Differential w.r.t. K_1^i variables stored in common array
        D_K1 = K(:,2).*I_2./k_23+K(:,3).*I_1./k_23;
        % Fréchet-Differential w.r.t. k_2^i variables stored in common array
        D_k2 = (K(:,1).*K(:,3)./(k_23.^2)).*(I_2-I_1)+(K(:,1).*K(:,2)./k_23).*(I_3-I_2.*t_sec);
        % Fréchet-Differential w.r.t. k_3^i variables stored in common array
        D_k3 = (K(:,1).*K(:,2)./(k_23.^2)).*(I_1-I_2)+(K(:,1).*K(:,2)./k_23).*(I_3-I_2.*t_sec);
        % Fréchet-Differential w.r.t. metabolic parameters
        DK = zeros(n,T,3);
        DK(:,:,1)=D_K1;
        DK(:,:,2)=D_k2;
        DK(:,:,3)=D_k3;
        DK = reshape(permute(DK,[3,1,2]),[3*n,T])';
        BK = kron(speye(n),ones(T,3));
        BK(logical(BK))=DK(:);
        
        % Build up Adjoint Operator of Fréchet-Differential
        DF_adj = [D_Lam;
                  D_Mu;
                  zeros(q,n*T),D_m;
                  full(BK)',zeros(3*n,T)];
              
        % Regularisation matrix
        S = blkdiag(alpha*eye(2*p),beta*eye(q),gamma*eye(3*n));
        
        % Iteratively Regularised Gauss-Newton Step
        x_n = x_n+(DF_adj*DF_adj'+S)\(DF_adj*(reshape((y_d-F)',[1,(n+1)*T])')+S*(x_star-x_n));

        % Refresh variables in new iteration step supplied by IRGNM
        % Metric projection to closed convex domain of operator F
        Lam = x_n(1:p)';
        Mu = min(x_n(p+1:2*p)', -eps);
        m_Biexp = [max(x_n(2*p+1),eps), min(x_n(2*p+2),-eps), min(x_n(2*p+3),-eps)]; 
        K(:,1)=max(x_n(2*p+q+1:3:end),eps);
        K(:,2)=max(x_n(2*p+q+2:3:end),eps);
        K(:,3)=max(x_n(2*p+q+3:3:end),eps);   
        k_23 = K(:,2)+K(:,3);

        x_n = [Lam, Mu, m_Biexp, reshape(K',[1,3*n])]'; % new iterative
        XN = [XN;x_n'];                                 % save iterative
        
        vMu= reshape(Mu,[1,1,p]);   % reshape parameters
        vLam=reshape(Lam,[1,1,p]);  % reshape parameters
        I_1 = (exp(vMu.*t_sec)-1)./vMu;
        I_2=(exp(vMu.*t_sec)-exp(-k_23.*t_sec))./(vMu+k_23);
        I_1 = I_1.*vLam;            % Overload Tensor I_1
        I_2 = I_2.*vLam;            % Overload Tensor I_2
        I_1 = sum(I_1,3);           % Build up Integral values
        I_2 = sum(I_2,3);           % Build up Integral values
        
        % Compute Forward Operator - Restriction w.r.t. C_T
        Res1 = (K(:,1).*K(:,2)./k_23).*I_2+(K(:,1).*K(:,3)./k_23).*I_1;  
        % Compute Forward Operator - Restriction w.r.t. C_TOT
        Res2 = (Lam*exp(Mu'.*t_sec))./(m_Biexp(1).*exp(m_Biexp(2).*t_sec)+(1-m_Biexp(1)).*exp(m_Biexp(3).*t_sec));

        % Build up Forward Operator
        F = [Res1;scale*Res2];
        y_d =[C_T_noise; scale*C_TOT_meas];  
        R = F-y_d;           % Residual element in Image Space w.r.t. y_d
        r=sqrt(trace(R'*R)); % Compute new l^2-Residual
        res = [res;r];       % Store new l^2-Residual in array
        iter = iter+1;       % Increase Iterator
        scale = scaleflag*(iter-1)+1;  % increase scale iff scaleflag
        x_tild = x_dagg-x_n; % Compute new deviation element in Preimage Space
        K_tild = K_dagg-K;   % Compute deviation for metabolic parameters
        % Extend Deviation Arrays
        x_res = [x_res; sqrt(x_tild(1:2*p)'*x_tild(1:2*p))/(sqrt(x_dagg(1:2*p)'*x_dagg(1:2*p))), sqrt(x_tild(2*p+1:2*p+q)'*x_tild(2*p+1:2*p+q))/(sqrt(x_dagg(2*p+1:2*p+q)'*x_dagg(2*p+1:2*p+q))), sqrt(x_tild(2*p+q+1:end)'*x_tild(2*p+q+1:end))/(sqrt(x_dagg(2*p+q+1:end)'*x_dagg(2*p+q+1:end)))];
        x_r = [x_r;sqrt(x_tild'*x_tild)/(sqrt(x_dagg'*x_dagg))];
        x_nof = [x_nof;sqrt(x_tild([1:2*p,2*p+q+1:end])'*x_tild([1:2*p,2*p+q+1:end]))/sqrt(x_dagg([1:2*p,2*p+q+1:end])'*x_dagg([1:2*p,2*p+q+1:end]))];
        x_kprm = [x_kprm;sqrt(sum(K_tild.^2))./sqrt(sum(K_dagg.^2)), sqrt(sum((Kcoeff_dagg-K(:,1).*K(:,3)./k_23).^2))./sqrt(sum(Kcoeff_dagg.^2))];
    end
end