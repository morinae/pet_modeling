function [Rep_K, Mean_K, Std_K]=reconstruction_noise(setup, Delta_y, comp, K, m_Biexp, Lam, Mu, t, C_TOT, CT_frontal, CT_temporal, CT_occipital, CT_wm, patient, counter, a1, a2, b1, b2, g1, g2, tau, scaleflag, kparamflag, C_TOT_noisy, fbv)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INPUT                                                              %
    % setup ...... 1 (full setup with f parameters) and 0 (reduced w.o.) %
    % Delta_y .... noise estimation of samples                           %
    % comp ....... if 1 (give comparison plot as in Figure 2 in Paper)   %
    % K .......... metabolic parameters                                  %
    % m_Biexp .... parent plasma fraction parameters (biexpontential)    %
    % Lam ........ factors of parameterization of plasma concentration   %
    % Mu ......... exponents of parameterization of plasma concentration %
    % t .......... temporal grid                (      in seconds     )  %
    % C_TOT ...... samples of whole blood tracer concentration           %
    % CT_frontal . noisy frontal tissue concentration                    %
    % CT_temporal  noisy temporal tissue concentration                   %
    % CT_occipital noisy occipital tissue concentration                  %
    % CT_wm ...... noisy white matter tissue concentration               %
    % patient .... Patient                                               %
    % counter .... count setting                                         %
    % a1 ......... regularisation parameter for plasma parameters        %
    % a2 ......... regularisation parameter for plasma parameters        %
    % b1 ......... regularisation parameter for ppf. function parameters %
    % b2 ......... regularisation parameter for ppf. function parameters %
    % g1 ......... regularisation parameter for metabolic parameters     %
    % g2 ......... regularisation parameter for metabolic parameters     %
    % tau ........ Discrepancy principle parameter (for stopping crit.)  %
    % scaleflag .. flag if scaling approach is applied to method         %
    % kparamflag . flag if metabolic deviations separately for K1,k2,k3  %
    % C_TOT_noisy  the noisy realizations of C_TOT (can be same as C_TOT)%
    % fbv ........ fractional blood volume                               %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OUTPUT                                                             %
    % Rep_K ...... metabolic reconstruction (for representative plot)    %
    % Mean_K ..... mean over metabolic reconstructions                   %
    % Std_K ...... standard deviation over metabolic reconstructions     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    rng(0);         % seed random number generator
    if setup
        stp = " Not Reduced -";
    else
        stp = " Reduced -";
    end
    % regularity of C_TOT
    if max(abs(C_TOT-mean(C_TOT)))<1e-8
        ctotmode = "noiseless";
    else
        ctotmode = "noisy";
    end
    % scaling
    if scaleflag
        scl = " yes";
    else
        scl = "no";
    end
    % figure name
    figure('Name', "Mode:"+stp+"Patient: "+patient+" - Counter: " + counter + " - C_T: noisy - C_TOT: "+ctotmode+" - Scale: "+scl)

    x_low = 0.00004;            % Deviation plot lower bound for ordinate
    x_upp = 1.5;                % Deviation plot upper bound for ordinate
    y_low = 1e-4/4;             % Residual plot lower bound for ordinate
    y_upp = 5;                  % Residual plot upper bound for ordinate
    Imax = 200;                 % Maximal amount of iterations (IRGNM)
    delta_x = [.4, .3, .2, .1]; % Array of initial error levels
    Opt = [];                   % Optimal Reductions for each error level
    Pos = [];                   % Respective iterations where optima attained
    RHO = [-1];                 % Reductions over multiple experiments for each error level (delimiters in array by value -1)
    RED=[];                     % Optimal reductions based on discrepancy principle
    POS = [];                   % Respective iterations where discrepandy principle attained
    RHO_DIV = zeros(1,4);       % Number of insufficient reductions in experiments for each error level
    exm = size(C_TOT,1);        % Number of experiments per error level 
    XNN = [];                   % Array to save reconstructed parameters (by discrepancy principle) for initial error levels
    YDD = [];                   % Array to save sample of corresponding noisy target
    C_TOT_N = [];               % Array to save sample of corresponding whole blood concentration
    for i = 1:size(delta_x,2)
        if setup
            X_RES = zeros(Imax,3,exm);  % Array to save separated relative deviations for all experiments at fixed error level
        else
            X_RES = zeros(Imax,2,exm);  % Array to save separated relative deviations for all experiments at fixed error level
        end
        X_R = [];                   % Array to save relative deviations for all experiments at fixed error level
        RES = [];                   % Array to save residuals for all experiments at fixed error level
        RHO_D = [];                 % Reduction of each experiment at fixed error level
        X_N = [];                   % Array to save reconstructed parameters for all experiments
        Y_D = [];                   % Array to save sample of noisy target for all experiments
        X_KPRM = [];                % Array to save sample of separated metabolic deviations
        C_TOT_noisy_i = C_TOT_noisy;% same noisy realizations of whole blood concentration for all initial error levels
        for s = 1:exm
            C_T = [CT_frontal(s,:);CT_temporal(s,:);CT_occipital(s,:);CT_wm(s,:)]; % sample of tissue tracer concentration
            C_TOT_s = C_TOT(s,:);   % sample of whole blood concentration
            delta = Delta_y(s);     % corresponding noise level
            % IRGNM
            if setup
                [x_res,x_r, res, iter,~, y_d,XN,x_kprm] = IRGNM_modified(t, delta_x(i),Imax, K, m_Biexp, Lam, Mu, C_TOT_s, C_T, a1, a2, b1, b2, g1, g2, scaleflag); 
            else 
                [x_res, x_r, res, iter, ~,~,y_d,XN,~] = IRGNM_reduced(t, delta_x(i),Imax, K, m_Biexp, Lam, Mu, C_TOT_s, C_T, a1, a2, g1, g2);
            end
            [~, pr]=max(res<=delta*tau); % position where disrepancy principle applies
            xn = XN(pr,:)';              % regarding reconstructed parameters
            yd = y_d(1:4,:)';            % regarding sample of noisy target
            Y_D = [Y_D,yd];              % save sample of noisy target
            if kparamflag
                X_KPRM = [X_KPRM,x_kprm];% save separate reconst. of metabolic parameters
            end
            vr = x_r(pr);                % regarding relative deviation
            red = floor(10000-10000*vr/x_r(1))/100; % regarding reduction (discrepancy principle)
            [v, p] = min(x_r);           % smallest relative deviation 
            opti=floor(10000-10000*v/x_r(1))/100;   % regarding reduction
            X_N = [X_N,xn];              % save reconstructed parameters
            X_RES(:,:,s) = x_res;        % save separated relative deviations
            X_R = [X_R,x_r];             % save relative deviations
            RES = [RES,res];             % save residuals
            RHO_D = [RHO_D, red];        % save discrepancy principle reductions
            if red<=0                    % if discrepancy principle reduction insufficient
                RHO_DIV(i)=RHO_DIV(i)+1; % increase number of occurred divergences
            end
        end
        bnd = 0;                         % lower bound for further classification (here no further classification)
        X_RES = X_RES(:,:,RHO_D>bnd);    % extract regarding separated relative deviations where no divergence occurs
        X_R = X_R(:,RHO_D>bnd);          % extract regarding relative deviations where no divergence occurs
        RES = RES(:,RHO_D>bnd);          % extract regarding residuals where no divergence occurs
        X_N = X_N(:, RHO_D>bnd);         % extract regarding reconstructed parameters where no divergence occurs
        C_TOT_noisy_i = C_TOT_noisy_i(RHO_D>bnd,:); % extract regarding whole blood concentrations where no divergence occurs
        msk = kron((RHO_D>bnd), ones(1,4));         % mask for extraction
        Y_D = Y_D(:,logical(msk));       % extract regarding sample of noisy target where no divergence occurs
        if kparamflag
            X_KPRM = X_KPRM(:,logical(msk));        % extract regarding sample of separated metabolic rec. where no div.
        end
        RHO_D = RHO_D(RHO_D>bnd);        % extract regarding optimal reductions where no divergence occurs
        RHO = [RHO,RHO_D,-1];            % extend reduction array with delimiter -1
        [~,position] = min(abs(RHO_D-median(RHO_D)));   % extract number of experiment with median optimal reduction
        x_res  = X_RES(:,:,position);    % extract regarding separated relative deviation
        x_r = X_R(:,position);           % extract regarding relative deviation
        res = RES(:,position);           % extract regarding residual
        x_n = X_N(:, position);          % extract regarding reconstructed parameters
        K_N = X_N(end-11:end,:);         % extract regarding reconstructed (representative) metabolic parameters
        M_K = mean(K_N,2);               % calculate mean of reconstructed metabolic parameters
        S_K = std(K_N,0,2);              % calculate std of reconstructed metabolic parameters
        % reshape
        Rep_K=[x_n(end-11:3:end),x_n(end-10:3:end),x_n(end-9:3:end)];
        Mean_K=[M_K(end-11:3:end),M_K(end-10:3:end),M_K(end-9:3:end)];
        Std_K=[S_K(end-11:3:end),S_K(end-10:3:end),S_K(end-9:3:end)];
        
        C_TOT_noisy_i = C_TOT_noisy_i(position,:);      % extract regarding sample of noisy whole blood concentration
        C_TOT_N = [C_TOT_N;C_TOT_noisy_i];              % save for all initial error levels
        y_d = Y_D(:,4*position-3:4*position);           % extract regarding sample of noisy target
        if kparamflag
            x_kprm = X_KPRM(:,4*position-3:4*position); % extract regarding sample of separate metabolic parameters
        end
        XNN = [XNN,x_n];                 % save median reconstructed parameters
        delta = Delta_y(s);              % estimated noise level
        YDD = [YDD, y_d];                % save median sample of noisy target
        [~, pr]=max(res<=delta*tau);     % position where disrepancy principle applies
        vr = x_r(pr);                    % regarding relative deviation
        red = floor(10000-10000*vr/x_r(1))/100; % regarding reduction (discrepancy principle)
        Red = num2str(red,'%05.2f');
        [v, p] = min(x_r);                      % smallest relative deviation 
        opti=floor(10000-10000*v/x_r(1))/100;   % regarding reduction
        RED = [RED,red];                        % save median reduction

        subplot(2,size(delta_x,2),i)    % plot resdiual-evolution
        semilogy(res,'LineWidth',1.6)   % logarithmic display
        hold on
        xBox = [1,Imax,Imax,1,1];       % plot discrepancy area
        yBox = [y_low, y_low, delta*tau, delta*tau, y_low];
        patch(xBox, yBox, 'black', 'FaceColor', 'red', 'FaceAlpha', 0.05);
        semilogy(1:size(res,1),delta*tau*ones(size(res)),'LineWidth',1.2)
        ylim([y_low y_upp])
        yticks([1e-4, 1e-3, 1e-2, 1e-1, 1e0])
        grid on
        ax = gca;
        ax.FontSize = 15; 
        axis square
        if (i==1)
            ylabel({'\textbf{Residual Value}','$$\Vert F\left(x^\delta_k\right)-y^\delta\Vert_Y$$',' '},'fontweight','bold','interpreter','latex','Fontsize',24)
        end
        xlabel('Iteration','Fontsize',18)
        xlim([1 iter])
        title({['$$\delta_x = ~ $$', num2str(delta_x(i))],' '},'fontweight','bold','interpreter','latex','Fontsize',24)

        subplot(2,size(delta_x,2),i+size(delta_x,2))    % plot deviation-evolutions
        semilogy(x_res(:,1),'LineWidth',2.6,'Color',[0.9290 0.6940 0.1250])
        hold on
        xBox = [pr,Imax,Imax,pr,pr]; % discrepancy visualisation for deviations
        yBox = [x_low,x_low,x_upp,x_upp,x_low];
        patch(xBox, yBox, 'black', 'FaceColor', 'red', 'FaceAlpha', 0.05,'HandleVisibility','off');
        semilogy([pr,pr], [x_low, x_upp],'LineWidth',1.2,'HandleVisibility','off') 
        if setup
            semilogy(x_res(:,2),'LineWidth',1.6,'Color',[0.4660 0.6740 0.1880])
            semilogy(x_res(:,3),'LineWidth',1.6,'Color',[0.4940 0.1840 0.5560])
            semilogy(x_r,'LineWidth',2,'Color',[0 0.4470 0.7410])
            if kparamflag
                semilogy(x_kprm(:,1),'LineWidth',1,'Color',"#4DBEEE")
                semilogy(x_kprm(:,2),'LineWidth',1,'Color',"#FF3E96")
                semilogy(x_kprm(:,3),'LineWidth',1,'Color',"#458B00")
                semilogy(x_kprm(:,4),'LineWidth',1,'Color',"#D95319")
                legend('C_{P}','f','K','All','K_1','k_2','k_3','K_1k_3/(k_2+k_3)','FontSize',12,'Location','Southwest','NumColumns',2)
            else
                legend('C_{P}','f','K','All','FontSize',12,'Location','Southwest')
            end
        else
            semilogy(x_res(:,2),'LineWidth',1.6,'Color',[0.4940 0.1840 0.5560])
            semilogy(x_r,'LineWidth',1.3,'Color',[0 0.4470 0.7410])
            legend('C_{P}','K','All','','','FontSize',12,'Location','Southwest')
        end
        
        ax = gca;
        ax.FontSize = 15; 
        axis square
        if (i==1)
            ylabel({'\textbf{Relative Deviation}','$$\Vert x^\dagger\Vert_X^{-1}\Vert x_k^\delta-x^\dagger\Vert_X$$',' '},'fontweight','bold','interpreter','latex','Fontsize',24)
        end
        title({['$$~ ~\rho_d = ~$$',Red,'$$\% ~ : ~ k = ~$$', num2str(pr)], ['$$\rho_{opt} = ~$$',num2str(opti,'%05.2f'),'$$\% ~ : ~ k = ~$$', num2str(p)]},'fontweight','light','interpreter','latex','Fontsize',18)
        xlabel('Iteration','FontSize',18)
        xlim([1 iter])
        ylim([x_low x_upp])
        yticks([1e-4, 1e-3, 1e-2, 1e-1, 1e0])
        grid on

        Opt = [Opt, opti];      % median optimal reductions
        Pos = [Pos, p];         % regarding iteration positions
        POS = [POS, pr];        % positions where discrepancy principle applies
    end

    set(gcf, 'PaperPositionMode', 'auto');
    print -depsc figure_size_control 
    Image = getframe(gcf);
    if setup
        str = append('irgnm_dy_count',num2str(counter),'.png');
    else
        str = append('noise_dy_count',num2str(counter),'.png');
    end

    disp(["Mode:"+stp+"Patient: "+patient+" - Counter: " + counter + " - C_T: noisy - C_TOT: "+ctotmode+" - Scale: "+scl])
    disp('Optimal reductions for error levels:    '+string(num2str(Opt)))
    disp('Optimal iterations for error levels:    '+string(num2str(Pos)))
    disp('Discrepancy reductions for error levels:    '+string(num2str(RED)))
    disp('Discrepancy iterations for error levels:    '+string(num2str(POS)))
    disp('Amount of divergent experiments (discrepancy) for error levels:    '+string(num2str(RHO_DIV)))
    disp(' ')
    disp(' ')
    
    if comp
        comparison_plot(setup,XNN, K, m_Biexp, Lam, Mu, YDD, counter, delta_x,t, C_TOT_N,patient,ctotmode,fbv)
    end
end