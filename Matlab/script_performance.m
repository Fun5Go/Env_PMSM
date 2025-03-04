clear vars
close all

% Load initial parameters for PMSM
init_params;
Te_ref = [60,120];

%% -------- SAVE DATA INTO CSV FILE -------- %%
% ----- PMSM Data -----%
Id      = Mot.mag_t.imd;
Iq      = Mot.mag_t.imq;
Ld      = 1e3*(Mot.mag_t.Lmidd + Mot.mag_t.Lmidq);
Lq      = 1e3*(Mot.mag_t.Lmiqq + Mot.mag_t.Lmidq);
Ldd     = 1e3*Mot.mag_t.Lmidd;
Lqq     = 1e3*Mot.mag_t.Lmiqq;
Phi_d   = Mot.mag_t.Psid;
Phi_q   = Mot.mag_t.Psiq;

% Reshape data
IQ = kron(Iq,ones(1,size(Id,2)))';
ID = kron(ones(1,size(Iq,2)),Id)';
LD = reshape(Ld,[],1);
LQ = reshape(Lq,[],1);
PHI_D = reshape(Phi_d,[],1);
PHI_Q = reshape(Phi_q,[],1);

% Reshape data 2
IQ_2 = kron(ones(1,size(Id,2)),Iq)';
ID_2 = kron(Id,ones(1,size(Iq,2)))';
LD_2 = reshape(Ld',[],1);
LQ_2 = reshape(Lq',[],1);
PHI_D_2 = reshape(Phi_d',[],1);
PHI_Q_2 = reshape(Phi_q',[],1);

% Ld
csv_data = [ID IQ LD];
writematrix(csv_data,'speedam_Ld_PMSM.csv','Delimiter','space')
plot3d_csv('speedam_Ld_PMSM.csv',size(Id,2));
% Lq
csv_data = [ID IQ LQ];
writematrix(csv_data,'speedam_Lq_PMSM.csv','Delimiter','space')
plot3d_csv('speedam_Lq_PMSM.csv',size(Id,2));
% Phi_d
csv_data = [ID IQ PHI_D];
writematrix(csv_data,'speedam_Phi_d_PMSM.csv','Delimiter','space')
plot3d_csv('speedam_Phi_d_PMSM.csv',size(Id,2));
% Phi_q
csv_data = [ID IQ PHI_Q];
writematrix(csv_data,'speedam_Phi_q_PMSM.csv','Delimiter','space')
plot3d_csv('speedam_Phi_q_PMSM.csv',size(Id,2));
% Ld_2
csv_data = [ID_2 IQ_2 LD_2];
writematrix(csv_data,'speedam_Ld_PMSM_2.csv','Delimiter','space')
plot3d_csv('speedam_Ld_PMSM_2.csv',size(Iq,2));
% Lq_2
csv_data = [ID_2 IQ_2 LQ_2];
writematrix(csv_data,'speedam_Lq_PMSM_2.csv','Delimiter','space')
plot3d_csv('speedam_Lq_PMSM_2.csv',size(Iq,2));
% Phi_d_2
csv_data = [ID_2 IQ_2 PHI_D_2];
writematrix(csv_data,'speedam_Phi_d_PMSM_2.csv','Delimiter','space')
plot3d_csv('speedam_Phi_d_PMSM_2.csv',size(Iq,2));
% Phi_q_2
csv_data = [ID_2 IQ_2 PHI_Q_2];
writematrix(csv_data,'speedam_Phi_q_PMSM_2.csv','Delimiter','space')
plot3d_csv('speedam_Phi_q_PMSM_2.csv',size(Iq,2));


% Plot lambda vs I
figure;
Iq_0 = find(Iq == 0,1);
for i = 1:size(Phi_d,1)
    lambda = abs(Phi_d(i,Iq_0:end) + 1i*Phi_q(i,Iq_0:end));
    plot(Iq(Iq_0:end), lambda, 'DisplayName','\lambda'); hold on
end
xlabel("I_d [A]"); ylabel("\lambda [Wb]"); title("\lambda vs I_q");

%% -------- SIMULATION ADAPTIVE PI -------- %%
fprintf("ADAPTIVE PI\n");
for Tref = Te_ref
    % Simulate
    PI_sim = sim('PMSM_PI_MTPA');
    
    PI = {};
    % Store data from simulation
    PI.Te.data      = PI_sim.Te.data;       % Te [Nm]
    PI.Ld.data      = PI_sim.Ld.data*1e3;   % Ld [mH]
    PI.Lq.data      = PI_sim.Lq.data*1e3;   % Lq [mH]
    PI.Id.data      = PI_sim.Id.data;       % Id [A]
    PI.Iq.data      = PI_sim.Iq.data;       % Iq [A]
    PI.Idref.data   = PI_sim.Idref.data;    % Idref [A]
    PI.Iqref.data   = PI_sim.Iqref.data;    % Iqref [A]
    PI.Phi_d.data   = PI_sim.Phi_d.data;    % Phi_d [Wb]   
    PI.Phi_q.data   = PI_sim.Phi_q.data;    % Phi_q [Wb]
    PI.Phi_m.data   = PI_sim.Phi_m.data;    % Phi_m [Wb]
    PI.noise.data   = PI_sim.noise.data;    % Phi_m [Wb]
    % Store time from simulation in [ms]
    PI.Te.time_ms   = PI_sim.Te.time*1e3;
    PI.Ld.time_ms   = PI_sim.Ld.time*1e3;
    PI.Lq.time_ms   = PI_sim.Lq.time*1e3;
    PI.Id.time_ms   = PI_sim.Id.time*1e3;
    PI.Iq.time_ms   = PI_sim.Iq.time*1e3;
    PI.Idref.time_ms= PI_sim.Idref.time*1e3;
    PI.Iqref.time_ms= PI_sim.Iqref.time*1e3;
    PI.Phi_d.time_ms= PI_sim.Phi_d.time*1e3;
    PI.Phi_q.time_ms= PI_sim.Phi_q.time*1e3;
    PI.Phi_m.time_ms= PI_sim.Phi_m.time*1e3;
    PI.noise.time_ms= PI_sim.noise.time*1e3;
    
    %% -------- PERFORMANCE PARAMETERS -------- %%
    ss_percentage   = 0.98;
    ss_criterion    = Tref*(1 - ss_percentage);
    index_ss        = find(abs(PI.Te.data - Tref) > ss_criterion,1,"last");
    % Transient
    PI.set_time     = PI.Te.time_ms(index_ss);
    PI.undershoot   = min(PI.Te.data(1:index_ss));
    PI.overshoot    = max(PI.Te.data(1:index_ss));
    % Steady-state
    PI.ss_data      = PI.Te.data(index_ss+100:end);
    PI.ss_mean      = mean(PI.ss_data);
    %PI.ss_ripple    = sum(abs(PI.ss_data - PI.ss_mean))/(size(PI.ss_data,1));
    PI.ss_ripple    = 1/Tref*(max(PI.ss_data) - min(PI.ss_data));
    PI.ss_ripple_percentage = 100*PI.ss_ripple;
    % Display relevant parameters
    fprintf("Te^* = %d,\nSettling time = %.2f [ms]\nSS Mean = %.2f [Nm]\n" + ...
        "SS Ripple = %.2f [Nm]\nSS Ripple/Reference = %.2f %% \n", ...
        Tref, PI.set_time, PI.ss_mean, PI.ss_ripple, PI.ss_ripple_percentage);
    
    %% -------- PLOTS -------- %%
    figure;
    plot(PI.Te.time_ms, Tref*ones(size(PI.Te.time_ms)),'DisplayName','T_{ref}'); hold on
    plot(PI.Te.time_ms, PI.Te.data,'DisplayName','T_e');
    xlabel("Time [ms]"); ylabel("Torque [Nm]"); title("Torque vs time");
    % Id vs time
    figure; 
    plot(PI.Id.time_ms, PI.Id.data, 'DisplayName','I_{d}'); hold on
    plot(PI.Idref.time_ms, PI.Idref.data, 'DisplayName','I_{d}');
    xlabel("Time [ms]"); ylabel("I_d [A]"); title("I_d vs time");
    % Iq vs time
    figure; 
    plot(PI.Iq.time_ms, PI.Iq.data, 'DisplayName','I_{q}'); hold on
    plot(PI.Iqref.time_ms, PI.Iqref.data, 'DisplayName','I_{qref}');
    xlabel("Time [ms]"); ylabel("I_q [A]"); title("I_q vs time");
    
   
    % % index
    % for i=0:size(Ld,2)
    %     fprintf("[%d %d]\n", i*25, i*25+24);
    % end
    
    % ----- Simulation Data -----%
    % Te
    csv_data = [PI.Te.time_ms PI.Te.data Tref*ones(size(PI.Te.time_ms))];
    writematrix(csv_data,strcat('speedam_Te_',num2str(Tref,'%02d'),'_adaptive_PI_MTPA.csv'))
    % Id
    csv_data = [PI.Id.time_ms PI.Id.data];
    writematrix(csv_data,strcat('speedam_Id_',num2str(Tref,'%02d'),'_adaptive_PI_MTPA.csv'))
    % Idref
    csv_data = [PI.Idref.time_ms PI.Idref.data];
    writematrix(csv_data,strcat('speedam_Idref_',num2str(Tref,'%02d'),'_adaptive_PI_MTPA.csv'))
    % Iq
    csv_data = [PI.Iq.time_ms PI.Iq.data];
    writematrix(csv_data,strcat('speedam_Iq_',num2str(Tref,'%02d'),'_adaptive_PI_MTPA.csv'))
    % Idref
    csv_data = [PI.Iqref.time_ms PI.Iqref.data];
    writematrix(csv_data,strcat('speedam_Iqref_',num2str(Tref,'%02d'),'_adaptive_PI_MTPA.csv'))
    % Ld Lq
    csv_data = [PI.Ld.time_ms PI.Ld.data PI.Lq.time_ms PI.Lq.data];
    writematrix(csv_data,strcat('speedam_LdLq_',num2str(Tref,'%02d'),'_adaptive_PI_MTPA.csv'))
    % Noise
    csv_data = [PI.noise.time_ms PI.noise.data];
    writematrix(csv_data,strcat('speedam_Te_noise.csv'))
    
end


%% -------- SIMULATION FIXED PI-------- %%
fprintf("FIXED PI\n");
Te_ref = [60,120];
for Tref = Te_ref
    if Tref == 60
        Kp_Id = 4.746;
        Ki_Id = 5.274e3;
        Kp_Iq = 8.861;
        Ki_Iq = 9.846e3;
    elseif Tref == 120
        Kp_Id = 2.957;
        Ki_Id = 3.285e3;
        Kp_Iq = 2.004e3;
        Ki_Iq = 9.846e3;
    else
        fprint("Unknown reference torque value")
    end
    % Simulate
    PI_sim = sim('PMSM_PI_fixed_MTPA');
    
    PI = {};
    % Store data from simulation
    PI.Te.data      = PI_sim.Te.data;       % Te [Nm]
    PI.Ld.data      = PI_sim.Ld.data*1e3;   % Ld [mH]
    PI.Lq.data      = PI_sim.Lq.data*1e3;   % Lq [mH]
    PI.Id.data      = PI_sim.Id.data;       % Id [A]
    PI.Iq.data      = PI_sim.Iq.data;       % Iq [A]
    PI.Idref.data   = PI_sim.Idref.data;    % Idref [A]
    PI.Iqref.data   = PI_sim.Iqref.data;    % Iqref [A]
    PI.Phi_d.data   = PI_sim.Phi_d.data;    % Phi_d [Wb]   
    PI.Phi_q.data   = PI_sim.Phi_q.data;    % Phi_q [Wb]
    PI.Phi_m.data   = PI_sim.Phi_m.data;    % Phi_m [Wb]
    % Store time from simulation in [ms]
    PI.Te.time_ms   = PI_sim.Te.time*1e3;
    PI.Ld.time_ms   = PI_sim.Ld.time*1e3;
    PI.Lq.time_ms   = PI_sim.Lq.time*1e3;
    PI.Id.time_ms   = PI_sim.Id.time*1e3;
    PI.Iq.time_ms   = PI_sim.Iq.time*1e3;
    PI.Idref.time_ms= PI_sim.Idref.time*1e3;
    PI.Iqref.time_ms= PI_sim.Iqref.time*1e3;
    PI.Phi_d.time_ms= PI_sim.Phi_d.time*1e3;
    PI.Phi_q.time_ms= PI_sim.Phi_q.time*1e3;
    PI.Phi_m.time_ms= PI_sim.Phi_m.time*1e3;
    
    %% -------- PERFORMANCE PARAMETERS -------- %%
    ss_percentage   = 0.98;
    ss_criterion    = Tref*(1 - ss_percentage);
    index_ss        = find(abs(PI.Te.data - Tref) > ss_criterion,1,"last");
    if index_ss == size(PI.Te.data,1)
        index_ss = find(PI.Te.time_ms == 5.0);
    end
    % Transient
    PI.set_time     = PI.Te.time_ms(index_ss);
    PI.undershoot   = min(PI.Te.data(1:index_ss));
    PI.overshoot    = max(PI.Te.data(1:index_ss));
    % Steady-state
    PI.ss_data      = PI.Te.data(index_ss+100:end);
    PI.ss_mean      = mean(PI.ss_data);
    %PI.ss_ripple    = sum(abs(PI.ss_data - PI.ss_mean))/(size(PI.ss_data,1));
    PI.ss_ripple    = 1/Tref*(max(PI.ss_data) - min(PI.ss_data));
    PI.ss_ripple_percentage = 100*PI.ss_ripple;
    % Display relevant parameters
    fprintf("Te^* = %d,\nSettling time = %.2f [ms]\nSS Mean = %.2f [Nm]\n" + ...
        "SS Ripple = %.2f [Nm]\nSS Ripple/Reference = %.2f %% \n", ...
        Tref, PI.set_time, PI.ss_mean, PI.ss_ripple, PI.ss_ripple_percentage);
    
    %% -------- PLOTS -------- %%
    figure;
    plot(PI.Te.time_ms, Tref*ones(size(PI.Te.time_ms)),'DisplayName','T_{ref}'); hold on
    plot(PI.Te.time_ms, PI.Te.data,'DisplayName','T_e');
    xlabel("Time [ms]"); ylabel("Torque [Nm]"); title("Torque vs time");
    % Id vs time
    figure; 
    plot(PI.Id.time_ms, PI.Id.data, 'DisplayName','I_{d}'); hold on
    plot(PI.Idref.time_ms, PI.Idref.data, 'DisplayName','I_{d}');
    xlabel("Time [ms]"); ylabel("I_d [A]"); title("I_d vs time");
    % Iq vs time
    figure; 
    plot(PI.Iq.time_ms, PI.Iq.data, 'DisplayName','I_{q}'); hold on
    plot(PI.Iqref.time_ms, PI.Iqref.data, 'DisplayName','I_{qref}');
    xlabel("Time [ms]"); ylabel("I_q [A]"); title("I_q vs time");
    
    % % index
    % for i=0:size(Ld,2)
    %     fprintf("[%d %d]\n", i*25, i*25+24);
    % end
    
    % ----- Simulation Data -----%
    % Te
    csv_data = [PI.Te.time_ms PI.Te.data Tref*ones(size(PI.Te.time_ms))];
    writematrix(csv_data,strcat('speedam_Te_',num2str(Tref,'%02d'),'_fixed_PI_MTPA.csv'))
    % Id
    csv_data = [PI.Id.time_ms PI.Id.data];
    writematrix(csv_data,strcat('speedam_Id_',num2str(Tref,'%02d'),'_fixed_PI_MTPA.csv'))
    % Idref
    csv_data = [PI.Idref.time_ms PI.Idref.data];
    writematrix(csv_data,strcat('speedam_Idref_',num2str(Tref,'%02d'),'_fixed_PI_MTPA.csv'))
    % Iq
    csv_data = [PI.Iq.time_ms PI.Iq.data];
    writematrix(csv_data,strcat('speedam_Iq_',num2str(Tref,'%02d'),'_fixed_PI_MTPA.csv'))
    % Idref
    csv_data = [PI.Iqref.time_ms PI.Iqref.data];
    writematrix(csv_data,strcat('speedam_Iqref_',num2str(Tref,'%02d'),'_fixed_PI_MTPA.csv'))
    % Ld Lq
    csv_data = [PI.Ld.time_ms PI.Ld.data PI.Lq.time_ms PI.Lq.data];
    writematrix(csv_data,strcat('speedam_LdLq_',num2str(Tref,'%02d'),'_fixed_PI_MTPA.csv'))
    
end

% figure;colormap('turbo');surf(X,Y,Ldd);colorbar;
% figure;colormap('jet');surf(X,Y,Ldd);colorbar;
% figure;colormap('parula');surf(X,Y,Ldd);colorbar;

% Function to rewrite csv file for 3dplot
function plot3d_csv(filename,rows)
    file  = fopen(filename,'r+');
    lines = textscan(file,'%s','delimiter','\n');
    fclose(file);
    file  = fopen(filename,'w+');
    for i = 1:numel(lines{1})
        fprintf(file,'%s\n',lines{1}{i});
        if mod(i,rows) == 0
            fprintf(file,'\n');
        end
    end
    fclose(file);
end


