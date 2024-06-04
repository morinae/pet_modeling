% simulation of FDG-like AIF and TACS
% AIF fit based on Feng model with parameter fitted from https://doi.org/10.1186/s40658-020-00330-x
% K values from https://journals.sagepub.com/doi/epdf/10.1038/jcbfm.1991.65

flagfig = 0; % flag for plotting ground truth curves

subjects = ["AD", "C"];

for i = 1:2

    subject = subjects(i);
    
    if subject == "C"
        K = [0.157, 0.174, 0.118;     % Frontal Cortex  [K_1^1, k_2^1, k_3^1; Control (C)
             0.161, 0.179, 0.096;     % Temporal Cortex  K_1^2, k_2^2, k_3^2; Control (C)
             0.177, 0.159, 0.088;     % Occipital Cortex K_1^3, k_2^3, k_3^3] Control (C)
             0.100, 0.161, 0.047];    % White matter     K_1^3, k_2^3, k_3^3] Control (C)
    elseif subject == "AD"
        K = [0.127, 0.172, 0.055;     % Frontal Cortex  [K_1^1, k_2^1, k_3^1; AD
             0.126, 0.211, 0.080;     % Temporal Cortex  K_1^2, k_2^2, k_3^2; AD
             0.182, 0.264, 0.079;     % Occipital Cortex K_1^3, k_2^3, k_3^3] AD
             0.102, 0.149, 0.013];    % White matter     K_1^3, k_2^3, k_3^3] AD
    else
        error("Invalid subject")
    end
    
    % Attenuation model (Bi-exponential: see available models Paper)
    m_Biexp = [0.2, -0.2, -0.005];  % Bi-exponential model
    
    % Arterial Concentration (Tri-exponential model)
    %Lam = [-(9.8154+0.8734),9.8154, 0.8734];     % lambda array in Tri-exponential model
    %Mu = [-13.0317,-2.8424,-0.01952];            % mu array in Tri-exponential model  
    
    %% Arterial Concentration (4-exponential model)
    Lam = [-(9.5450+0.7331+0.6355), 9.5450, 0.7331, 0.6355];     % lambda array in 4-exponential model
    Mu = [-13.4522,-3.2672,-0.15324,-0.01055];            % mu array in 4-exponential model  
    
    frm_dur_s = [5, 5, 5, 5, ... 
                 10, 10, 10, 10, ...
                 30, 30, 30, 30, ...
                 60, 60, ...
                 150, 150, 150, ...
                 300, 300, 300, 300, 300, 300 ...
                 600, 600];
    frm_end_time_s = cumsum(frm_dur_s);
    frm_start_time_s = [0, frm_end_time_s(1:end-1) ];
    
    % frame time (middle of frame)
    t = 0.5*(frm_start_time_s + frm_end_time_s);
    
    [C_T, f, C_P, C_TOT] = forward_model(K, m_Biexp, Lam, Mu, t);
    
    save("tacs_" + subject + ".h5","t","C_TOT","C_P","C_T","f","frm_dur_s", "K", "m_Biexp", "Lam", "Mu", "-v7.3");
    if flagfig
        fig = figure;
        fig.Position = [100 100 1200 400];
    
        subplot(1,3,1);
        plot(t/60,C_T(1,:), '-o'); 
        hold on; 
        plot(t/60,C_T(2,:), '-o');
        hold on; 
        plot(t/60,C_T(3,:), '-o'); 
        hold on; 
        plot(t/60,C_T(4,:), '-o'); 
        hold off;
        legend('C_T frontal','C_T temporal','C_T occipital','C_T white matter','Location','northwest')
        
        subplot(1,3,2);
        plot(t/60,C_TOT, '-o');
        hold on;
        plot(t/60,C_P, '-o');
        hold off;
        legend('C_{TOT}', 'C_{P}')
    
        subplot(1,3,3);
        plot(t/60,f, '-o');
        hold off;
        legend('f')
    end
end
