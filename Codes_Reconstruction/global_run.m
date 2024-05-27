% patients = ["AD", "C"];
patients = ["C"];
counters = [.1,1,10,0]; % 0 ... noiseless
no_pat = size(patients,2);
no_count = size(counters,2);

% Hyperparameters for reference deviation delta_x = 0.3, pat AD, counter 1
% columns = [a1, a2, b1, b2, g1, g2, tau]
hyper = [100, 7, 200, 7, 4000, 7, 6.8;  % C_tot_true, scale, mod, not red.
         100, 8, 400, 8, 3000, 8, 17.6; % C_tot_noise, no s., mod, not red.
         600, 7, 0,0, 10, 5, 9.2];      % C_tot_true, scale, mod, reduced
c=0;

% Generation of plot in Paper
warning('off','all')
diary reconstruction_info
 for np = 1:no_pat
     patient = patients(np);
     for nc = 1:no_count
         counter = counters(nc);
         if counter ~=0
             run("read_data.m")
             check_theoretical_assumptions(K,Mu,Lam,t_min)
             t = t_min*60;
         end
         setup = [0,1];
         if (counter==1)
             c=1;
         else
             c=0;
         end
         if counter == 0
             reconstruction_noiseless(1,K, m_Biexp, Lam, Mu, t, ones(size(C_TOT_noisy,1),1)*C_TOT_true, patient, hyper(1,1), hyper(1,2),hyper(1,3),hyper(1,4),hyper(1,5),hyper(1,6), 1)
             reconstruction_noiseless(1,K, m_Biexp, Lam, Mu, t, ones(size(C_TOT_noisy,1),1)*C_TOT_true, patient, hyper(2,1), hyper(2,2),hyper(2,3),hyper(2,4),hyper(2,5),hyper(2,6),0)
             reconstruction_noiseless(0,K, m_Biexp, Lam, Mu, t, ones(size(C_TOT_noisy,1),1)*C_TOT_true, patient, hyper(3,1), hyper(3,2),hyper(3,3),hyper(3,4),hyper(3,5),hyper(3,6),1)
         else
             [Rep_K, Mean_K, Std_K]=reconstruction_noise(1,Delta_y,c,K, m_Biexp, Lam, Mu, t,ones(size(C_TOT_noisy,1),1)*C_TOT_true, CT_recon_frontal, CT_recon_temporal, CT_recon_occipital, CT_recon_wm, patient, counter, hyper(1,1), hyper(1,2),hyper(1,3),hyper(1,4),hyper(1,5),hyper(1,6), hyper(1,7), 1,0,C_TOT_noisy,fbv)
             [Rep_K, Mean_K, Std_K]=reconstruction_noise(1,Delta_y,c,K, m_Biexp, Lam, Mu, t,C_TOT_noisy, CT_recon_frontal, CT_recon_temporal, CT_recon_occipital, CT_recon_wm, patient, counter, hyper(2,1), hyper(2,2),hyper(2,3),hyper(2,4),hyper(2,5),hyper(2,6), hyper(2,7),0,0,C_TOT_noisy,fbv)
             [Rep_K, Mean_K, Std_K]=reconstruction_noise(0,Delta_y,c,K, m_Biexp, Lam, Mu, t,ones(size(C_TOT_noisy,1),1)*C_TOT_true, CT_recon_frontal, CT_recon_temporal, CT_recon_occipital, CT_recon_wm, patient, counter, hyper(3,1), hyper(3,2),hyper(3,3),hyper(3,4),hyper(3,5),hyper(3,6), hyper(3,7), 1,0,C_TOT_noisy,fbv)
         end
     end
 end
 diary off