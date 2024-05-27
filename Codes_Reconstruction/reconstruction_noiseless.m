function reconstruction_noiseless(setup, K, m_Biexp, Lam, Mu, t, C_TOT, patient, a1, a2, b1, b2, g1, g2, scaleflag)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INPUT                                                            %
    % setup .... 1 (full setup with f parameters) and 0 (reduced w.o.) %
    % K ........ metabolic parameters                                  %
    % m_Biexp .. parent plasma fraction parameters (biexpontential)    %
    % Lam ...... factors of parameterization of plasma concentration   %
    % Mu ....... exponents of parameterization of plasma concentration %
    % t ........ temporal grid                (      in seconds     )  %
    % C_TOT .... samples of whole blood tracer concentration           %
    % a1 ....... regularisation parameter for plasma parameters        %
    % a2 ....... regularisation parameter for plasma parameters        %
    % b1 ....... regularisation parameter for ppf. function parameters %
    % b2 ....... regularisation parameter for ppf. function parameters %
    % g1 ....... regularisation parameter for metabolic parameters     %
    % g2 ....... regularisation parameter for metabolic parameters     %
    % scaleflag  flag if scaling approach is applied to method         %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    rng(0);     % seed random number generator
    % setup
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
    figure('Name', "Mode:"+stp+"Patient: "+patient+" - Counter: oo" + " - C_T: noiseless - C_TOT: "+ctotmode+" - Scale: "+scl)

    x_low = 0.000000002;            % Deviation plot lower bound for ordinate
    x_upp = 1.5;                    % Deviation plot upper bound for ordinate
    y_low = 0.000000002;            % Residual plot lower bound for ordinate
    y_upp = 5;                      % Residual plot upper bound for ordinate
    Imax = 300;                     % Maximal amount of iterations (IRGNM)
    delta_x = [.4, .3, .2, .1];  % Array of initial error levels
    Opt = [];                       % Optimal Reductions for each error level
    Pos = [];                       % Respective iterations where optima attained
    RHO = [-1];                     % Reductions over multiple experiments for each error level (delimiters in array by value -1)
    RED = [];                       % In noiseless case RED = Opt
    RHO_DIV = zeros(1,4);           % Number of insufficient reductions in experiments for each error level
    exm = size(C_TOT,1);            % Number of experiments per error level
    [C_T, ~,~,~]=forward_model(K, m_Biexp, Lam, Mu, t); % forward model ground truth parameters
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
        
        for s = 1:exm               % for each realization
            C_TOT_s = C_TOT(s,:);   % corresponding whole blood concentration realization
            % IRGNM
            if setup
                [x_res,x_r, res, iter,~, y_d,XN,~] = IRGNM_modified(t, delta_x(i),Imax,K, m_Biexp, Lam, Mu, C_TOT_s, C_T, a1, a2, b1, b2, g1, g2, scaleflag); 
            else
                [x_res, x_r, res, iter, ~,~,y_d,XN,~] = IRGNM_reduced(t, delta_x(i),Imax,K, m_Biexp, Lam, Mu, C_TOT_s, C_T,  a1, a2, g1, g2);
            end
            [v, pos] = min(x_r);    % smallest relative deviation of single experiment
            opti=floor(10000-10000*v/x_r(1))/100;   % respective reduction
            xn = XN(pos,:)';        % regarding reconstructed parameters
            yd = y_d(3,:)';         % regarding sample of noisy target
            Y_D = [Y_D,yd];         % save sample of noisy target
            X_N = [X_N,xn];         % save reconstructed parameters
            X_RES(:,:,s) = x_res;   % save separated relative deviations
            X_R = [X_R,x_r];        % save relative deviations
            RES = [RES,res];        % save residuals
            RHO_D = [RHO_D, opti];  % save optimal reduction
            if opti<=0              % if no convergence becomes apparent
                RHO_DIV(i)=RHO_DIV(i)+1;    % increase number of occurred divergences
            end
        end
        bnd = 0;                        % lower bound for further classification (here no further classification)
        X_RES = X_RES(:,:,RHO_D>bnd);   % extract regarding separated relative deviations where no divergence occurs
        X_R = X_R(:,RHO_D>bnd);         % extract regarding relative deviations where no divergence occurs
        RES = RES(:,RHO_D>bnd);         % extract regarding residuals where no divergence occurs
        X_N = X_N(:, RHO_D>bnd);        % extract regarding reconstructed parameters where no divergence occurs
        Y_D = Y_D(:, RHO_D>bnd);        % extract regarding sample of noisy target where no divergence occurs
        RHO_D = RHO_D(RHO_D>bnd);       % extract regarding optimal reductions where no divergence occurs
        RHO = [RHO,RHO_D,-1];           % extend reduction array with delimiter -1
        [~,position] = min(abs(RHO_D-median(RHO_D)));   % extract number of experiment with median optimal reduction
        x_res  = X_RES(:,:,position);   % extract regarding separated relative deviation
        x_r = X_R(:,position);          % extract regarding relative deviation
        res = RES(:,position);          % extract regarding residual
        x_n = X_N(:, position);         % extract regarding reconstructed parameters
        y_d = Y_D(:,position);          % extract regarding sample of noisy target
        [v, p] = min(x_r);              % get iteration position and value of median optimal reduction
        opti=vpa(floor(10000000-10000000*vpa(v/x_r(1),5))/100000);  % calculate reduction
        RED = [RED,opti];               % save median reduction

        subplot(2,size(delta_x,2),i)    % plot residual-evolution
        semilogy(res,'LineWidth',1.6)   % logarithmic display
        ylim([y_low y_upp])             % ordinate bounds
        yticks([1e-7 1e-5 1e-3 1e-1])
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
        if setup
            semilogy(x_res(:,1),'LineWidth',2,'Color',[0.9290 0.6940 0.1250])
            hold on
            semilogy(x_res(:,2),'LineWidth',1.6,'Color',[0.4660 0.6740 0.1880])
            semilogy(x_res(:,3),'LineWidth',1.6,'Color',[0.4940 0.1840 0.5560])
            semilogy(x_r,'LineWidth',2,'Color',[0 0.4470 0.7410])
            legend('C_{P}','f','K','All','FontSize',12,'Location','Southwest')
        else
            semilogy(x_res(:,1),'LineWidth',2,'Color',[0.9290 0.6940 0.1250])
            hold on
            semilogy(x_res(:,2),'LineWidth',1.6,'Color',[0.4940 0.1840 0.5560])
            semilogy(x_r,'LineWidth',1.3,'Color',[0 0.4470 0.7410])
            legend('C_{P}','K','All','FontSize',12,'Location','Southwest')
        end
            
        yticks([1e-7 1e-5 1e-3 1e-1])
        ax = gca;
        ax.FontSize = 15;
        axis square
        if (i==1)
            ylabel({'\textbf{Relative Deviation}','$$\Vert x^\dagger\Vert_X^{-1}\Vert x_k^\delta-x^\dagger\Vert_X$$',' '},'fontweight','bold','interpreter','latex','Fontsize',24)
        end
        title(['$$\rho_{opt} = ~$$',num2str(double(opti),7),'$$\% ~ : ~ k = ~$$', num2str(p)],'fontweight','light','interpreter','latex','Fontsize',16)
        xlabel('Iteration','FontSize',18)
        xlim([1 iter])
        ylim([x_low x_upp])
        grid on

        Opt = [Opt, opti];  % median optimal reductions
        Pos = [Pos, p];     % regarding iteration positions
    end

    set(gcf, 'PaperPositionMode', 'auto');
    print -depsc figure_size_control
    Image = getframe(gcf);
    
    disp(["Mode:"+stp+"Patient: "+patient+" - Counter: oo" + " - C_T: noiseless - C_TOT: "+ctotmode+" - Scale: "+scl])
    disp(['Optimal reductions for error levels:    ',num2str(double(Opt))])
    disp('Optimal iterations for error levels:    '+string(num2str(Pos)))
    disp('Amount of divergent experiments for error levels:    '+string(num2str(RHO_DIV)))
    disp(' ')
    disp(' ')
end