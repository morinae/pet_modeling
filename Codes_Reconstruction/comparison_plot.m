function comparison_plot(setup, XNN, K, m_Biexp, Lam, Mu, YDD, counter, delta_x, t, C_TOT_N, patient, ctotmode, fbv)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INPUT                                                            %
    % setup .... 1 (Full setup) and 0 (Reduced setup)                  %
    % XNN ...... reconstructed parameters for different delta_x        %
    % K ........ ground truth metabolic parameters                     %
    % m_Biexp .. ground truth parameters of function f                 %
    % Lam ...... ground truth multiplicative parameters of C_P         %
    % Mu ....... ground truth exponential parameters of C_P            %
    % YDD ...... representative reconstructions of C_T                 %
    % counter .. count setting                                         %
    % delta_x .. initial error levels                                  %
    % t ........ times                                                 %
    % C_TOT_N .. noisy realizations of C_TOT for representative plots  %
    % patient .. Patient                                               %
    % ctotmode . regularity (noiseless/noisy) of C_TOT                 %
    % fbv ...... fractional blood volume                               %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    d=2; % plot w.r.t. delta_x(2) = 0.3
    if setup
        figure('Name', "Comparison - full at counter = "+num2str(counter) + " and delta_x = "+num2str(delta_x(d))+" for patient "+patient+" and C_TOT " +ctotmode)
    else
        figure('Name', "Comparison - reduced at counter = "+num2str(counter) + " and delta_x = "+num2str(delta_x(d))+" for patient "+patient+" and C_TOT " +ctotmode)
    end
    % ground truth forward model
    [C_T, f, C_P, C_TOT]=forward_model(K, m_Biexp, Lam, Mu,t); 
    C_PET = (1-fbv)*C_T+fbv*C_TOT;  % ground truth tissue time activity                              
    T = size(t,2);                  % number of times              
    t_sec = t/60;                   % renormalize to minutes
    n=size(K,1);                    % number of regions
    q=size(m_Biexp,2);              % number of parameters of function f
    p = size(Lam,2);                % degree of plasma concentration
    
    for i = d
        % reconstruction corresponding to representative plots
        K_rec = zeros(n,3);         % reconstructed metabolic parameters
        Lam_rec = XNN(1:p,i)';      % reconstructed mult. parameters of C_P
        Mu_rec = XNN(p+1:2*p,i)';   % reconstructed exp. parameters of C_P
        m_rec = XNN(2*p+1:2*p+q,i)';% reconstructed param. of function f
        if ~setup                   % if reduced setup no function f
            q=0;
        end
        K_rec(:,1) = XNN(2*p+q+1:3:end,i); % reconstructed K_1 parameters
        K_rec(:,2) = XNN(2*p+q+2:3:end,i); % reconstructed k_2 parameters
        K_rec(:,3) = XNN(2*p+q+3:3:end,i); % reconstructed k_3 parameters
        C_TOT_noisy = C_TOT_N(i,:);        % noisy realization of C_TOT
    
        % Forward model reconstruction of parameters of representative plot
        [C_T_rec, f_rec, C_P_rec, C_TOT_rec]=forward_model(K_rec, m_rec, Lam_rec, Mu_rec,t);
        C_PET_rec = (1-fbv)*C_T_rec+fbv*C_TOT_rec;          % reconstructed time activity curve
        PET_noise = (1-fbv)*YDD(:, 4*i-3:4*i)'+fbv*C_TOT_N; % noisy time activity curve
        
        if setup            % for full setup also f and C_TOT
            subplot(2,2,1)  % regarding function f
            semilogx(t_sec, f_rec,'LineWidth',2,'Color',[0.8500 0.3250 0.0980])
            hold on
            semilogx(t_sec,f,'LineWidth',2,'LineStyle','--','Color',[0 0.4470 0.7410])
            title(['Parent plasma fraction ', '$$f$$'],'fontweight','bold','interpreter','latex','Fontsize',20)
            xlabel('Time (min)','Fontsize',12)
            ylim([0.5,1.05])
            legend('Reconstruction','Ground truth', 'Fontsize',12)
            grid on
            ax = gca;
            ax.FontSize = 16; 

            subplot(2,2,2)  % regarding whole blood concentration C_TOT
            semilogx(t_sec, C_TOT_rec,'LineWidth',2,'Color',[0.8500 0.3250 0.0980])
            hold on
            semilogx(t_sec,C_TOT,'LineWidth',2,'LineStyle','--','Color',[0 0.4470 0.7410])
                if ctotmode == "noisy"
                    semilogx(t_sec, C_TOT_noisy,'r*','Color','#CD0000','Marker','*','MarkerSize',8)
                    legend('Reconstruction','Ground truth','Measurements', 'Fontsize',12)
                else
                    legend('Reconstruction','Ground truth', 'Fontsize',12)
                end
            title(['Total arterial blood tracer concentration ', '$$C_{WB}$$'],'fontweight','bold','interpreter','latex','Fontsize',20)
            xlabel('Time (min)','Fontsize',12)
            grid on 
            ax = gca;
            ax.FontSize = 16; 
        end

        subplot(2,2,3)  % regarding plasma concentration
        semilogx(t_sec, C_P_rec,'LineWidth',2,'Color',[0.8500 0.3250 0.0980])
        hold on
        semilogx(t_sec, C_P,'LineWidth',2,'LineStyle','--','Color',[0 0.4470 0.7410])
        title(['Plasma concentration ', '$$C_P$$'],'fontweight','bold','interpreter','latex','Fontsize',20)
        xlabel('Time (min)','Fontsize',12)
        legend('Reconstruction','Ground truth', 'Fontsize',12)
        grid on
        ax = gca;
        ax.FontSize = 16; 
        
        subplot(2,2,4) % regarding time activity curve
        semilogx(t_sec, C_PET_rec(1,:),'LineWidth',2,'Color',[0.8500 0.3250 0.0980])
        hold on
        semilogx(t_sec,C_PET(1,:),'LineWidth',2,'Color',[0.6350 0.0780 0.1840],'LineStyle','--')
        semilogx(t_sec, PET_noise(1,:),'r*','Color',[205,0,0]/255,'Marker','*','MarkerSize',8)

        semilogx(t_sec, C_PET_rec(2,:),'LineWidth',2,'Color',[154,205,50]/255)
        semilogx(t_sec,C_PET(2,:),'LineWidth',2,'Color',[0,139,0]/255,'LineStyle','--')
        semilogx(t_sec, PET_noise(2,:),'r*','Color','#006400','Marker','*','MarkerSize',8)

        semilogx(t_sec, C_PET_rec(3,:),'LineWidth',2,'Color',[0,205,205]/255)
        semilogx(t_sec,C_PET(3,:),'LineWidth',2,'Color',[0,0,139]/255,'LineStyle','--')
        semilogx(t_sec, PET_noise(3,:),'r*','Color',[0,0,205]/255,'Marker','*','MarkerSize',8)

        semilogx(t_sec, C_PET_rec(4,:),'LineWidth',2,'Color',[205,41,144]/255)
        semilogx(t_sec,C_PET(4,:),'LineWidth',2,'Color',[139,0,139]/255,'LineStyle','--')
        semilogx(t_sec, PET_noise(4,:),'r*','Color',[148,0,211]/255,'Marker','*','MarkerSize',8)
        title(['Tissue time activity curve ', '$$C_{PET}$$'],'fontweight','bold','interpreter','latex','Fontsize',20)
        xlabel('Time (min)','Fontsize',12)
        legend('recon.: frontal', 'gr. truth: frontal', 'meas.: frontal','recon.: temporal','gr. truth: temporal',  'meas.: temporal', ...
             'recon.: occipital', 'gr. truth: occipital', 'meas.: occipital','recon.: white matter','gr. truth: white matter',  'meas.: white matter', ...
             'NumColumns',2,'Location','NorthWest', 'Fontsize',12);
        grid on
        ax = gca;
        ax.FontSize = 16; 
        
        fig = gcf ;
        set(fig,'PaperType','A4')
    end
end