function [C_T, f, C_P, C_TOT]=forward_model(K, m_Biexp, Lam, Mu, t)
    % Input
    % K (n x 3):      1st column K_1 parameters
    %                 2nd column k_2 parameters
    %                 3rd column k_3 parameters
    %                 Rows are given by n different regions
    %
    % m_Biexp (1x3):  1st entry is multiplicative factor
    %                 2nd and 3rd entry are exponential exponents
    %                 biexponential model for ppf
    %
    % Lam (1xp):      polyexponential model for C_P
    %                 consists of multiplicative factors
    %
    % Mu (1xp):       polyexponential model for C_P
    %                 consists of exponential exponents
    %
    % t:              Times in seconds

    T = size(t,2);          % Amount of Time measurements
    t_sec = t/60;           % Time normalisation in minutes
    n = size(K,1);          % Number of considered regions
    k_23 = K(:,2)+K(:,3);   % Additive values k_2+k_3
    q = size(m_Biexp,2);    % number of parameters for function f
    p = size(Lam,2);        % degree of plasma concentration
    
    % reshape C_P parameters
    vMu= reshape(Mu,[1,1,p]);
    vLam=reshape(Lam,[1,1,p]);
    
    % Polyexponential Integral Tensors
    % Int(exp(mu_i s)): dimensions 1xTxp
    I_1 = (exp(vMu.*t_sec)-1)./vMu;

    % Int(exp((k2+k3)(s-t))exp(mu_i s)): dimensions nxTxp
    I_2=(exp(vMu.*t_sec)-exp(-k_23.*t_sec))./(vMu+k_23);

    % Overload Tensors
    I_1 = I_1.*vLam;
    I_2 = I_2.*vLam;

    % Build up Integral values
    I_1 = sum(I_1,3);
    I_2 = sum(I_2,3);

    % Foward Model for tracer concentration in tissue
    C_T = (K(:,1).*K(:,2)./k_23).*I_2+(K(:,1).*K(:,3)./k_23).*I_1;
    % Foward Model for plasma concentration
    C_P = Lam*exp(Mu'.*t_sec);
    % Forward Model for function f
    f = m_Biexp(1)*exp(m_Biexp(2)*t_sec)+(1-m_Biexp(1))*exp(m_Biexp(3)*t_sec);
    % Forward Model for whole blood tracer concentration
    C_TOT = C_P./f;
end