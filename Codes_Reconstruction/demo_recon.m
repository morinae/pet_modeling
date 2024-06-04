% Quick demo example based on dataset of patient C

% simulation of FDG-like AIF and TACS
% AIF fit based on Feng model with parameter fitted from https://doi.org/10.1186/s40658-020-00330-x
% K values from https://journals.sagepub.com/doi/epdf/10.1038/jcbfm.1991.65

K = [0.157, 0.174, 0.118;     % Frontal Cortex  [K_1^1, k_2^1, k_3^1; Control (C)
     0.161, 0.179, 0.096;     % Temporal Cortex  K_1^2, k_2^2, k_3^2; Control (C)
     0.177, 0.159, 0.088;     % Occipital Cortex K_1^3, k_2^3, k_3^3] Control (C)
     0.100, 0.161, 0.047];    % White matter     K_1^3, k_2^3, k_3^3] Control (C)

% Number of regions
n = size(K,1);

% Attenuation model (Bi-exponential: see available models Paper)
m_Biexp = [0.2, -0.2, -0.005];  % Bi-exponential model

% Arterial Concentration (4-exponential model)
Lam = [-(9.5450+0.7331+0.6355), 9.5450, 0.7331, 0.6355];    % lambda array in 4-exponential model
Mu = [-13.4522,-3.2672,-0.15324,-0.01055];            		% mu array in 4-exponential model  

% Frames
frm_dur_s = [5, 5, 5, 5, ... 
             10, 10, 10, 10, ...
             30, 30, 30, 30, ...
             60, 60, ...
             150, 150, 150, ...
             300, 300, 300, 300, 300, 300 ...
             600, 600];
frm_end_time_s = cumsum(frm_dur_s);
frm_start_time_s = [0, frm_end_time_s(1:end-1) ];

% Frame time (middle of frame)
t = 0.5*(frm_start_time_s + frm_end_time_s);

% Frame time in minutes
t_min = t/60;

% Number of time frames
T = size(t_min,2);

% Fractional blood volume
fbv = 0.05;

% Ground truth curves
[C_T, f, C_P, C_TOT_true] = forward_model(K, m_Biexp, Lam, Mu, t);
CT_true_frontal = C_T(1,:);
CT_true_temporal = C_T(2,:);
CT_true_occipital = C_T(3,:);
CT_true_wm = C_T(4,:);

% Number of experiments
ex = 20;

% Artificial noise for whole blood tracer conc. and tracer con. in tissue
delta_y = 1e-3;
fac = 2;
Delta_y = fac*ones(ex,1)*delta_y;
rng(0);         % seed random number generator
C_TOT_noisy = abs(C_TOT_true + randn(ex,size(C_TOT_true,2))*delta_y);
CT_recon_frontal = abs(CT_true_frontal + randn(ex,size(CT_true_frontal,2))*delta_y);
CT_recon_temporal = abs(CT_true_temporal + randn(ex,size(CT_true_temporal,2))*delta_y);
CT_recon_occipital = abs(CT_true_occipital + randn(ex,size(CT_true_occipital,2))*delta_y);
CT_recon_wm = abs(CT_true_wm + randn(ex,size(CT_true_wm,2))*delta_y);

% Hyperparameters
hyper = [100, 7, 200, 7, 4000, 7, 6.8;  % C_tot_true, scale, mod, not red.
         100, 8, 400, 8, 3000, 8, 17.6; % C_tot_noise, no s., mod, not red.
         600, 7, 0,0, 10, 5, 9.2];      % C_tot_true, scale, mod, reduced
c=1;

% Reconstruction
warning('off','all')

patient = "C";
counter = 1e-3/delta_y; % note that in this demo "counter" has no physical meaning as opposed to the .py data generation file where it is obtained differently
check_theoretical_assumptions(K,Mu,Lam,t_min)

% Noiseless setting
% Reduced setup
reconstruction_noiseless(0,K, m_Biexp, Lam, Mu, t, ones(size(C_TOT_noisy,1),1)*C_TOT_true, patient, hyper(3,1), hyper(3,2),hyper(3,3),hyper(3,4),hyper(3,5),hyper(3,6),1)
% Full setup
reconstruction_noiseless(1,K, m_Biexp, Lam, Mu, t, ones(size(C_TOT_noisy,1),1)*C_TOT_true, patient, hyper(1,1), hyper(1,2),hyper(1,3),hyper(1,4),hyper(1,5),hyper(1,6), 1)
  
% Noisy setting
% Reduced setup with noiseless whole blood tissue concentration
[Rep_K, Mean_K, Std_K]=reconstruction_noise(0,Delta_y,c,K, m_Biexp, Lam, Mu, t,ones(size(C_TOT_noisy,1),1)*C_TOT_true, CT_recon_frontal, CT_recon_temporal, CT_recon_occipital, CT_recon_wm, patient, counter, hyper(3,1), hyper(3,2),hyper(3,3),hyper(3,4),hyper(3,5),hyper(3,6), hyper(3,7), 1,0,C_TOT_noisy,fbv)
% Full setup with noiseless whole blood tissue concentration
[Rep_K, Mean_K, Std_K]=reconstruction_noise(1,Delta_y,c,K, m_Biexp, Lam, Mu, t,ones(size(C_TOT_noisy,1),1)*C_TOT_true, CT_recon_frontal, CT_recon_temporal, CT_recon_occipital, CT_recon_wm, patient, counter, hyper(1,1), hyper(1,2),hyper(1,3),hyper(1,4),hyper(1,5),hyper(1,6), hyper(1,7), 1,0,C_TOT_noisy,fbv)
% Full setup with noisy whole blood tissue concentration
[Rep_K, Mean_K, Std_K]=reconstruction_noise(1,Delta_y,c,K, m_Biexp, Lam, Mu, t,C_TOT_noisy, CT_recon_frontal, CT_recon_temporal, CT_recon_occipital, CT_recon_wm, patient, counter, hyper(2,1), hyper(2,2),hyper(2,3),hyper(2,4),hyper(2,5),hyper(2,6), hyper(2,7),0,0,C_TOT_noisy,fbv)
