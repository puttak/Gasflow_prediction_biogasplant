function[gasflow_complete]=ADM_Sim_ML_remote(ch, pr, li, dil)

%File for running an implementation of anaerobic digester simulation with 
% the ADM1 model. 

% clear all 
% close all
rng(0,'twister');
% retrive input data and parameter values 
l = IndataADM1_v3;

y0=[0.009;... % S_su 
    0.0009;... % S_aa 
    0.0009;... % S_fa
%      0.0009;... % S_va 
    0.5;... % S_va 
%      0.0009;... % S_bu 
    0.5;... % S_bu 
%      0.0009;... % S_pro 
    0.5;... % S_pro
%      0.0009;... % S_ac 
    1;... % S_ac 
    2.3594e-9;...%S_h2 
    2.3594e-6;...%S_ch4 
    0.039;... %S_IC 
    0.13023;... %S_IN 
    0.009;... %S_I 
    0.30870;... %X_c 
    0.02795;... %X_ch 
    0.10260;... %X_pr 
    0.02948; ... %X_li 
    0.42016;... %X_su 
    1.17917;... %X_aa 
    0.24303;... %X_fa 
%     0.43192;... %X_c4 
%     0.13730;... %X_pro 
%     0.76056;... %X_ac 
    0.283;... %X_c4 
    0.136;... %X_pro 
    0.9;... %X_ac 
    0.31702;... %X_h2 
    25.61739;... %X_I 
    0.04;... %S_cat 
    0.02;... %S_an 
    0.0116;... %S_vam 
    0.01322;... %S_bum
    0.01574;... %S_prom 
    0.19724;... %S_acm 
    0.14278;... %S_hco3m 
    0.00409;... %S_nh3 
    1.023e-5;... %h2 
    1.62125;... %ch4 
    0.01411;... %co2 
    0;... %S_lac 
    0;... %X_lac_f 
    0;... %X_lac_o 
    0]; %S_ca

t=linspace(0,300,100);

t_red = py_teil(t);

% t_red_1 = t(1,2:5:(end/2))';
% t_red_2 = t(1,(end/2):10:end)';
%Solve the ODE-system 
options = odeset('RelTol', 1e-10, 'AbsTol', 1e-10); 

high = 4000; % Anzahl Variationen
n_size = high; 

% high = 5; % Anzahl Variationen
%% für X
y0_size = 3*high;
y1 = repmat(y0,1,y0_size);

%%
% ch(1:high) = 1.*rand(1,n_size); % -
% pr(1:high) = rand(1,n_size).*(1 - (ch)); % -
% li(1:high) = 1 - (ch + pr); % -
% 
% pr(high+1:(2*high)) = 1.*rand(1,n_size); % -
% ch(high+1:(2*high)) = rand(1,n_size).*(1 - (pr(high+1:(2*high)))); % -
% li(high+1:(2*high)) = 1 - (ch(high+1:(2*high)) + pr(high+1:(2*high))); % -
% 
% li((2*high)+1:(3*high)) = 1.*rand(1,n_size); % -
% pr((2*high)+1:(3*high)) = rand(1,n_size).*(1 - (li((2*high)+1:(3*high)))); % -
% ch((2*high)+1:(3*high)) = 1 - (li((2*high)+1:(3*high)) + pr((2*high)+1:(3*high))); % -
%%

% ch = 0.8147;
% pr = 0.0791;
% li = 0.1062;
% dil= 0.1338;

% figure;
% scatter(1:3*high, ch);
% 
% figure;
% scatter(1:3*high, pr);
% 
% figure;
% scatter(1:3*high, li);
%% für ingredients
sI = 0; % -
xI = 0; % -


%% für Dilution
dil_size = 3*high;
a1 = 0.01;
b1 = 0.2;

%%
% gasflow = zeros(size(y1,2),1);
% for iv = 1:size(y0,1)
    
%     X_var = zeros(size(y1,2),1);
%     y1 = y1_orig;

%     number = iv;
%     a = y1(number);
%     b = 100.*y1(number);
%     y1(number,:) = (b-a).*rand(1,y0_size) + a;

%%
% dil = (b1-a1).*rand(1,dil_size) + a1;
%%

% ix = [33:35];
% a2 = y0(ix)/10;
% b2 = y0(ix)*100;
% y1(ix,:) = (b2-a2).*rand(1,dil_size) + a2;

% dil = repmat(0.05,1,high);

% end
%%

% for ii = 1:(y0_size)
% parfor ii = 1:size(dil,2)
ii = 1;
%% optimierte Substratzusammensetzung   
%     in = [0; 0; ch; pr; li; dil(1,ii)];
    in = [0; 0; ch(1,ii); pr(1,ii); li(1,ii); dil(1,ii)];

    [T,solution] = ode15s(@(t,y) ADM1_fun_v2_ODE_1_remote(t,y,l,in), t, y1(:,ii), options);

    S_gas_h2 = solution(:,33);
    S_gas_ch4 = solution(:,34);
    S_gas_co2 = solution(:,35);

    %% gas pressure
    P_gas_h2 = S_gas_h2*l.R*l.T_op/16;
    P_gas_ch4 = S_gas_ch4*l.R*l.T_op/64;
    P_gas_co2 = S_gas_co2*l.R*l.T_op;
    
    %% gas flow
    P_gas = P_gas_h2+P_gas_ch4+P_gas_co2+l.p_gas_h2o;
    
    %% total gas pressure
    q_gas = l.k_p*(P_gas-l.P_atm).*P_gas/l.P_atm;

    %% total gas flow %calculate gas component gas flow

    q_gas_ch4 = P_gas_ch4./P_gas.*q_gas;

    gasfl = (q_gas_ch4);
    disp(gasfl(end,:));
    
    gasflow_zwischen_1(:,ii) = gasfl(2:5:end/2,:);
    gasflow_zwischen_2(:,ii) = gasfl(end/2:10:end,:);

% end

gasflow_points = [gasflow_zwischen_1; gasflow_zwischen_2];
% t_red = [t_red_1;t_red_2];

t_full = t';
% gasflow = [gasflow_points, t_red];

gasflow_complete = [gasfl, T]; 


% in_li  = li';
% in_pr  = pr';
% in_ch  = ch';
% D      = dil';
% % gases  = y1(33:35,:)'; 
% 
% save(strcat('gas_sub_dilhoch_points',num2str(size(gasflow_points,2))),'t_red',...
%     'gasflow_points','D','in_ch','in_pr','in_li');

end