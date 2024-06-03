%% Choose Patient and Counter-Frequency
% patient = "C";
% counter = .1;

% part of file name based on patient dataset
switch patient
    case "AD"
        str = "../pet_data_sim/data/AD_";
        rstr = "tacs_AD_frm_0";
    case "C"
        str = "../pet_data_sim/data/C_";
        rstr = "tacs_C_frm_0";
    otherwise
        warning('Invalid patient type')
end

% part of file name based on on counter setting
switch counter
    case 10
        estr = "_s0_19_sens_10.0_recons.h5";
        cstr = str + "10/";
    case 1
        estr = "_s0_19_sens_1.0_recons.h5";
        cstr = str + "1/";
    case 0.1
        estr = "_s0_19_sens_0.1_recons.h5";
        cstr = str + "01/";
    otherwise
        warning('Invalid counter frequency')
end
cstr = cstr + rstr;
str = cstr + "00" + estr;

% Regions
% 1st: frontal cortex
% 2nd: temporal cortex
% 3rd: occipital cortex
% 4th: white matter
C_P = h5read(str, "/C_P");                  % load plasma conc. C_P
C_T = h5read(str, "/C_T");                  % load tissue conc. C_T
C_TOT = h5read(str, "/C_TOT");              % load total conc. C_TOT (noise)
K = h5read(str, "/K");                      % load metabolic param. K
Lam = h5read(str, "/Lam");                  % load C_P factors Lambda
Mu = h5read(str, "/Mu");                    % load C_P exponents Mu
fbv = h5read(str, "/fbv");                  % load fractional blood vol.
frm_dur_s = h5read(str, "/frm_dur_s");      % load time frame duration
m_Biexp = h5read(str, "/m_Biexp");          % load ppf parameters
t_min = h5read(str, "/t_min");              % load time grid

eps = 10e-8;                                % computational zero  
T = size(t_min,1);                          % number of time frames
ex = size(h5read(str,"/C_recon_frontal"),1);% number of experiments
n = size(K,2);                              % number of regions
C_recon_frontal = zeros(ex,T);              % initialise recon. frontal
C_recon_occipital = zeros(ex,T);            % initialise recon. occipital
C_recon_temporal = zeros(ex,T);             % initialise recon. temporal
C_recon_wm = zeros(ex,T);                   % initialise recon. white matt.
C_TOT_recon = zeros(ex,T);                  % initialise recon. tot. conc.

C_true_frontal = zeros(1,T);                % initialise true frontal
C_true_occipital = zeros(1,T);              % initialise true occipital
C_true_temporal = zeros(1,T);               % initialise true temporal
C_true_wm = zeros(1,T);                     % initialise true white matt.
C_TOT_true = zeros(1,T);                    % initialise true tot. conc.

C_recon_frontal_std = zeros(ex,T);          % init. recon. frontal std
C_recon_temporal_std = zeros(ex,T);         % init. recon. occipital std
C_recon_occipital_std = zeros(ex,T);        % init. recon. temporal std
C_recon_wm_std = zeros(ex,T);               % init. recon. white mat. std

for i = 0:24
    % iterate over data files and generate corresponding dataset name
    if numel(num2str(i))==1
        str = cstr + "0" + num2str(i)+ estr;
    else
        str = cstr + num2str(i) + estr;
    end
    % check if global parameters coincide
    assert(sum(abs(C_P-h5read(str, "/C_P")))<=eps);
    assert(sum(abs(C_T-h5read(str, "/C_T")),'all')<=eps);
    assert(sum(abs(C_TOT-h5read(str, "/C_TOT")))<=eps);
    assert(sum(abs(K-h5read(str, "/K")),'all')<=eps);
    assert(sum(abs(Lam-h5read(str, "/Lam")))<=eps);
    assert(sum(abs(Mu-h5read(str, "/Mu")))<=eps);
    assert(sum(abs(fbv-h5read(str, "/fbv")))<=eps);
    assert(sum(abs(frm_dur_s-h5read(str, "/frm_dur_s")))<=eps);
    assert(sum(abs(m_Biexp-h5read(str, "/m_Biexp")))<=eps);
    assert(sum(abs(t_min-h5read(str, "/t_min")))<=eps);
    % load reconstructions
    C_recon_frontal(:,i+1)=h5read(str, "/C_recon_frontal");
    C_recon_occipital(:,i+1)=h5read(str, "/C_recon_occipital");
    C_recon_temporal(:,i+1)=h5read(str, "/C_recon_temporal");
    C_recon_wm(:,i+1)=h5read(str, "/C_recon_wm");
    C_TOT_recon(:,i+1) = h5read(str, "/C_TOT_recon");
    % load ground truth values
    C_true_frontal(i+1)=h5read(str, "/C_true_frontal");
    C_true_occipital(i+1)=h5read(str, "/C_true_occipital");
    C_true_temporal(i+1)=h5read(str, "/C_true_temporal");
    C_true_wm(i+1)=h5read(str, "/C_true_wm");
    C_TOT_true(i+1) = h5read(str, "/C_TOT_true");
    % load standard deviations
    C_recon_frontal_std(:,i+1) = h5read(str, "/C_recon_frontal_std");
    C_recon_temporal_std(:,i+1) = h5read(str, "/C_recon_temporal_std");
    C_recon_occipital_std(:,i+1) = h5read(str, "/C_recon_occipital_std");
    C_recon_wm_std(:,i+1) = h5read(str, "/C_recon_wm_std");
end
% tranpose parameter arrays for parameter identification algorithm
C_P = C_P';
C_T = C_T';
C_TOT = C_TOT';
K = K';
Lam = Lam';
Mu = Mu';
frm_dur_s = frm_dur_s';
m_Biexp = m_Biexp';
t_min = t_min';

% Quantative noise estimation by mean of variation coefficients
Delta_y = (mean(C_recon_wm_std./C_recon_wm,2)+mean(C_recon_frontal_std./C_recon_frontal,2)+mean(C_recon_temporal_std./C_recon_temporal,2)+mean(C_recon_occipital_std./C_recon_occipital,2))/4;
Delta_y = Delta_y/400; % normalization factor
% Bias correction factors
beta = zeros(n,T);  % bias correction for tissue time activity curve
beta(1,:) = C_true_frontal./mean(C_recon_frontal);
beta(2,:) = C_true_temporal./mean(C_recon_temporal);
beta(3,:) = C_true_occipital./mean(C_recon_occipital);
beta(4,:) = C_true_wm./mean(C_recon_wm);
% bias correction for whole blood concentration
gamma = C_TOT_true./mean(C_TOT_recon);  

% Preprocessing of noisy C_PET data to noisy C_T under bias correction
CT_recon_frontal = (beta(1,:).*C_recon_frontal-fbv*C_TOT)/(1-fbv);
CT_recon_temporal = (beta(2,:).*C_recon_temporal-fbv*C_TOT)/(1-fbv);
CT_recon_occipital = (beta(3,:).*C_recon_occipital-fbv*C_TOT)/(1-fbv);
CT_recon_wm = (beta(4,:).*C_recon_wm-fbv*C_TOT)/(1-fbv);
% Preprocessing of noisy C_TOT data under bias correction
C_TOT_noisy = gamma.*C_TOT_recon;
