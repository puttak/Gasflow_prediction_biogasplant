function ADM1_dt=ADM1_fun_v2_ODE_1(t,y,l,in)

l.f_sI_xc = in(1);
l.f_xI_xc = in(2);
l.f_ch_xc = in(3);
l.f_pr_xc = in(4);
l.f_li_xc = in(5);
l.D = in(6);

% if t<5
%     l.D=0;
% else
%     l.D=l.D;
% end

S_su=y(1); 
S_aa=y(2); 
S_fa=y(3); 
S_va=y(4); 
S_bu=y(5);
S_pro=y(6); 
S_ac=y(7); 
S_h2=y(8);
S_ch4=y(9);
S_IC=y(10); 
S_IN=y(11); 
S_I=y(12);

X_c=y(13);

X_ch=y(14); 
X_pr=y(15);
X_li=y(16);
X_su=y(17); 
X_aa=y(18); 
X_fa=y(19); 
X_c4=y(20); 
X_pro=y(21);
X_ac=y(22); 
X_h2=y(23); 
X_I=y(24); 
S_cat=y(25); 
S_an=y(26); 
S_vam=y(27); 
S_bum=y(28); 
S_prom=y(29); 
S_acm=y(30); 
S_hco3m=y(31); 
S_nh3=y(32); 
S_gas_h2=y(33); 
S_gas_ch4=y(34); 
S_gas_co2=y(35); 
%.................................................. 
%Extra input to ADM1 model(added more substrates) 
%.................................................. 
S_lac=y(36); 
X_lac_f=y(37); 
X_lac_o=y(38); 
S_ca=y(39);
%==========================================================================
%calculation of gas pressure 
%========================================================================== 
P_gas_h2=S_gas_h2*l.R*l.T_op/16;
P_gas_ch4=S_gas_ch4*l.R*l.T_op/64;
P_gas_co2=S_gas_co2*l.R*l.T_op; 
P_gas=P_gas_h2+P_gas_ch4+P_gas_co2+l.p_gas_h2o; 
q_gas=l.k_p*(P_gas-l.P_atm)*P_gas/l.P_atm;
%========================================================================== 
%calculation of potential Sh 
%========================================================================== 
S_nh4=S_IN-S_nh3; 
S_co2=S_IC-S_hco3m;
Theta=S_cat+S_nh4-S_hco3m-(S_acm/64)-(S_prom/112)-(S_bum/160)...
    -(S_vam/208)-S_an+S_ca;
%added S_ca to orginal equation 
Sh= -Theta/2+sqrt(Theta^2+4*l.K_w)/2; 
if(Sh<=0) 
    Sh=1e-12; 
end
pH=-log10(Sh);

%pH inhibition 
if(pH<l.pH_UL_aa) 
    I_pH_aa=exp(-3*((pH-l.pH_UL_aa)/(l.pH_UL_aa-l.pH_LL_aa))^2); 
else
    I_pH_aa=1; 
end

if(pH<l.pH_UL_ac) 
    I_pH_ac=exp(-3*((pH-l.pH_UL_ac)/(l.pH_UL_ac-l.pH_LL_ac))^2); 
else
    I_pH_ac=1; 
end

if(pH<l.pH_UL_h2) 
    I_pH_h2=exp(-3*((pH-l.pH_UL_h2)/(l.pH_UL_h2-l.pH_LL_h2))^2); 
else
    I_pH_h2=1; 
end

%========================================================================== 
%Process inhibition 
%========================================================================== 
I_IN_lim=1/(1+l.K_S_IN/S_IN); 
I_h2_fa=1/(1+S_h2/l.K_I_h2_fa); 
I_h2_c4=1/(1+S_h2/l.K_I_h2_c4); 
I_h2_pro=1/(1+S_h2/l.K_I_h2_pro); 
I_nh3=1/(1+S_nh3/l.K_I_nh3); 
I_5=I_pH_aa*I_IN_lim; 
I_6=I_5; 
I_7=I_pH_aa*I_IN_lim*I_h2_fa; 
I_8=I_pH_aa*I_IN_lim*I_h2_c4; 
I_9=I_8;
I_10=I_pH_aa*I_IN_lim*I_h2_pro; 
I_11=I_pH_ac*I_IN_lim*I_nh3; 
I_12=I_pH_h2*I_IN_lim; 
%.................................................. 
%Extra input to ADM1 model
%.................................................. 
I_h2_lac_o=1/(1+S_h2/l.K_I_h2_lac_o); 
I_20=I_IN_lim; 
I_21=I_h2_lac_o*I_IN_lim;

%========================================================================== 
%Biochemical processrates 
%========================================================================== 
rho_1=l.k_dis*X_c; 
rho_2=l.k_hyd_ch*X_ch; 
rho_3=l.k_hyd_pr*X_pr; 
rho_4=l.k_hyd_li*X_li; 
rho_5=l.k_m_su*(S_su/(l.K_S_su+S_su))*X_su*I_5; 
rho_6=l.k_m_aa*(S_aa/(l.K_S_aa+S_aa))*X_aa*I_6; 
rho_7=l.k_m_fa*(S_fa/(l.K_S_fa+S_fa))*X_fa*I_7; 
rho_8=l.k_m_c4*(S_va/(l.K_S_c4+S_va))*X_c4*(S_va/(S_bu+S_va+1e-6))*I_8;
rho_9=l.k_m_c4*(S_bu/(l.K_S_c4+S_bu))*X_c4*(S_bu/(S_va+S_bu+1e-6))*I_9;
rho_10=l.k_m_pro*(S_pro/(l.K_S_pro+S_pro))*X_pro*I_10; 
rho_11=l.k_m_ac*(S_ac/(l.K_S_ac+S_ac))*X_ac*I_11; 
rho_12=l.k_m_h2*(S_h2/(l.K_S_h2+S_h2))*X_h2*I_12;
rho_13=l.k_dec_X_su*X_su; 
rho_14=l.k_dec_X_aa*X_aa; 
rho_15=l.k_dec_X_fa*X_fa; 
rho_16=l.k_dec_X_c4*X_c4; 
rho_17=l.k_dec_X_pro*X_pro; 
rho_18=l.k_dec_X_ac*X_ac; 
rho_19=l.k_dec_X_h2*X_h2; 
%.................................................. 
%Extra input to ADM1 model
%.................................................. 
rho_20=l.k_m_lac_f*(S_lac/(l.K_S_lac_f+S_lac))*X_lac_f*I_20; 
rho_21=l.k_m_lac_o*(S_lac/(l.K_S_lac_o+S_lac))*X_lac_o*I_21;
rho_22=l.k_dec_X_lac_f*X_lac_f;
rho_23=l.k_dec_X_lac_o*X_lac_o;

%========================================================================== 
%Acid-base rates 
%========================================================================== 
rho_A_4=l.k_A_B_va*(S_vam*(l.K_a_va+Sh)-l.K_a_va*S_va); 
rho_A_5=l.k_A_B_bu*(S_bum*(l.K_a_bu+Sh)-l.K_a_bu*S_bu); 
rho_A_6=l.k_A_B_pro*(S_prom*(l.K_a_pro+Sh)-l.K_a_pro*S_pro); 
rho_A_7=l.k_A_B_ac*(S_acm*(l.K_a_ac+Sh)-l.K_a_ac*S_ac); 
rho_A_10=l.k_A_B_co2*(S_hco3m*(l.K_a_co2+Sh)-l.K_a_co2*S_IC);
rho_A_11=l.k_A_B_IN*(S_nh3*(l.K_a_IN+Sh)-l.K_a_IN*S_IN);
%========================================================================== 
%Gas transfer rates
%========================================================================== 
rho_T_8=l.k_L_a_h2*(S_h2-16*l.K_H_h2*P_gas_h2); 
rho_T_9=l.k_L_a_ch4*(S_ch4-64*l.K_H_ch4*P_gas_ch4); 
rho_T_10=l.k_L_a_co2*(S_co2-l.K_H_co2*P_gas_co2);
%========================================================================== 
%Precipitation rates
%========================================================================== 
%to avoid complex numbers prevent it from being negative

if(S_ca<1e-6||S_hco3m<1e-6) 
    rho_P_24=0; 
else
    rho_P_24=l.K_r_caco3*((S_ca^.5*S_hco3m^.5-l.K_S_p_caco3^.5)^2); % korrigiert DW 28.04.18
end
%========================================================================== 
%Processes(carbon) 
%========================================================================== 
s_1=-l.C_xc+l.f_sI_xc*l.C_sI+l.f_ch_xc*l.C_ch+l.f_pr_xc*l.C_pr+l.f_li_xc*l.C_li+l.f_xI_xc*l.C_xI; 
s_2=-l.C_ch+l.C_su; 
s_3=-l.C_pr+l.C_aa; 
s_4=-l.C_li+(1-l.f_fa_li)*l.C_su+l.f_fa_li*l.C_fa;
s_5=-l.C_su+(1-l.Y_su)*(l.f_bu_su*l.C_bu+l.f_pro_su*l.C_pro+l.f_ac_su*l.C_ac)+l.Y_su*l.C_bac; 
s_6=-l.C_aa+(1-l.Y_aa)*(l.f_va_aa*l.C_va+l.f_bu_aa*l.C_bu+l.f_pro_aa*l.C_pro+l.f_ac_aa*l.C_ac)+... 
    l.Y_aa*l.C_bac; 

%% abgeändert (mit den % Verhältnissen erweitert)
s_7=-l.C_fa+(1-l.Y_fa)*l.f_ac_fa*l.C_ac+l.Y_fa*l.C_bac;  % l_f_ac_fa eingeführt 
s_8=-l.C_va+(1-l.Y_c4)*l.f_pro_va*l.C_pro+(1-l.Y_c4)*l.f_ac_va*l.C_ac+l.Y_c4*l.C_bac; % l.f_pro_va, l.f_ac_va eingeführt
s_9=-l.C_bu+(1-l.Y_c4)*l.f_ac_bu*l.C_ac+l.Y_c4*l.C_bac; % l.f_ac_bu eingeführt
s_10=-l.C_pro+(1-l.Y_pro)*l.f_ac_pro*l.C_ac+l.Y_pro*l.C_bac; % l.f_ac_pro eingeführt
s_11=-l.C_ac+(1-l.Y_ac)*l.C_ch4+l.Y_ac*l.C_bac; 
s_12=(1-l.Y_h2)*l.C_ch4+l.Y_h2*l.C_bac; 
s_13=-l.C_bac+l.C_xc; 
%.................................................. 
%Extra input to ADM1 model
%.................................................. 
s_20=-l.C_lac+(1-l.Y_lac_f)*0.785*l.C_pro+(1-l.Y_lac_f)*0.215*l.C_ac+l.Y_lac_f*l.C_bac; 
s_21=-l.C_lac+(1-l.Y_lac_o)*2/3*l.C_ac+l.Y_lac_o*l.C_bac;
%========================================================================== 
%Differential equations 
%========================================================================== 
%.......................................................................... 
%Differential equations 1-4,soluble matter 
%.......................................................................... StateNo. 
% d_S_su_dt=l.q_in/l.V_liq*(l.S_su_in-S_su)+rho_2+(1-l.f_fa_li)*rho_4-rho_5;  %1
% d_S_aa_dt=l.q_in/l.V_liq*(l.S_aa_in-S_aa)+rho_3-rho_6;                      %2 
% d_S_fa_dt=l.q_in/l.V_liq*(l.S_fa_in-S_fa)+l.f_fa_li*rho_4-rho_7;            %3 
% d_S_va_dt=l.q_in/l.V_liq*(l.S_va_in-S_va)+(1-l.Y_aa)*l.f_va_aa*rho_6-rho_8; %4 

d_S_su_dt=l.D*(l.S_su_in-S_su)+rho_2+(1-l.f_fa_li)*rho_4-rho_5;  %1
d_S_aa_dt=l.D*(l.S_aa_in-S_aa)+rho_3-rho_6;                      %2 
d_S_fa_dt=l.D*(l.S_fa_in-S_fa)+l.f_fa_li*rho_4-rho_7;            %3 
d_S_va_dt=l.D*(l.S_va_in-S_va)+(1-l.Y_aa)*l.f_va_aa*rho_6-rho_8; %4 

%.......................................................................... 
%Differential equations 5-8, soluble matter 
%.......................................................................... StateNo. 
% d_S_bu_dt=l.q_in/l.V_liq*(l.S_bu_in-S_bu)+(1-l.Y_su)*l.f_bu_su*rho_5+...    %5 
%     (1-l.Y_aa)*l.f_bu_aa*rho_6-rho_9;                                       %? 
% d_S_pro_dt=l.q_in/l.V_liq *(l.S_pro_in-S_pro)+(1-l.Y_su)*l.f_pro_su*rho_5+...%6(addedlac) 
%     (1 -l.Y_aa)*l.f_pro_aa*rho_6+(1-l.Y_c4)*0.54*rho_8-rho_10+...            %? 
%     (1 -l.Y_lac_f)*0.785*rho_20;                                             %? 
% d_S_ac_dt=l.q_in/l.V_liq *(l.S_ac_in-S_ac)+(1-l.Y_su)*l.f_ac_su*rho_5+...   %7(addedlac) 
%     (1 -l.Y_aa)*l.f_ac_aa*rho_6+(1-l.Y_fa)*0.7*rho_7+...                    %? 
%     (1 -l.Y_c4)*0.31*rho_8+(1-l.Y_c4)*0.8*rho_9+...                         %? 
%     (1 -l.Y_pro)*0.57*rho_10-rho_11+(1-l.Y_lac_f)*0.215*rho_20+...          %? 
%     (1 -l.Y_lac_o)*2/3*rho_21;                                              %? 
% d_S_h2_dt=l.q_in/l.V_liq *(l.S_h2_in-S_h2)+(1-l.Y_su)*l.f_h2_su*rho_5+...   %8(addedlac) 
%     (1 -l.Y_aa)*l.f_h2_aa*rho_6+(1-l.Y_fa)*0.3*rho_7+...                    %? 
%     (1 -l.Y_c4)*0.15*rho_8+(1-l.Y_c4)*0.2*rho_9+...                         %? 
%     (1 -l.Y_pro)*0.43*rho_10-rho_12-rho_T_8+(1-l.Y_lac_o)*1/3*rho_21;       %? 

% %% D eingeführt
% d_S_bu_dt=l.D*(l.S_bu_in-S_bu)+(1-l.Y_su)*l.f_bu_su*rho_5+...    %5 
%     (1-l.Y_aa)*l.f_bu_aa*rho_6-rho_9;                                       %? 
% d_S_pro_dt=l.D *(l.S_pro_in-S_pro)+(1-l.Y_su)*l.f_pro_su*rho_5+...%6(addedlac) 
%     (1 -l.Y_aa)*l.f_pro_aa*rho_6+(1-l.Y_c4)*0.54*rho_8-rho_10+...            %? 
%     (1 -l.Y_lac_f)*0.785*rho_20;                                             %? 
% d_S_ac_dt=l.D *(l.S_ac_in-S_ac)+(1-l.Y_su)*l.f_ac_su*rho_5+...   %7(addedlac) 
%     (1 -l.Y_aa)*l.f_ac_aa*rho_6+(1-l.Y_fa)*0.7*rho_7+...                    %? 
%     (1 -l.Y_c4)*0.31*rho_8+(1-l.Y_c4)*0.8*rho_9+...                         %? 
%     (1 -l.Y_pro)*0.57*rho_10-rho_11+(1-l.Y_lac_f)*0.215*rho_20+...          %? 
%     (1 -l.Y_lac_o)*2/3*rho_21;                                              %? 
% d_S_h2_dt=l.D *(l.S_h2_in-S_h2)+(1-l.Y_su)*l.f_h2_su*rho_5+...   %8(addedlac) 
%     (1 -l.Y_aa)*l.f_h2_aa*rho_6+(1-l.Y_fa)*0.3*rho_7+...                    %? 
%     (1 -l.Y_c4)*0.15*rho_8+(1-l.Y_c4)*0.2*rho_9+...                         %? 
%     (1 -l.Y_pro)*0.43*rho_10-rho_12-rho_T_8+(1-l.Y_lac_o)*1/3*rho_21;       %? 

%% D eingeführt, + weitere prozentuale Beziehungen
d_S_bu_dt=l.D*(l.S_bu_in-S_bu)+(1-l.Y_su)*l.f_bu_su*rho_5+...               %5 
    (1-l.Y_aa)*l.f_bu_aa*rho_6-rho_9;                                       
d_S_pro_dt=l.D *(l.S_pro_in-S_pro)+(1-l.Y_su)*l.f_pro_su*rho_5+...          %6(addedlac) 
    (1 -l.Y_aa)*l.f_pro_aa*rho_6+(1-l.Y_c4)*l.f_pro_va*rho_8-rho_10+...     
    (1 -l.Y_lac_f)*0.785*rho_20;                                            
d_S_ac_dt=l.D *(l.S_ac_in-S_ac)+(1-l.Y_su)*l.f_ac_su*rho_5+...              %7(addedlac) 
    (1 -l.Y_aa)*l.f_ac_aa*rho_6+(1-l.Y_fa)*l.f_ac_fa*rho_7+...                    
    (1 -l.Y_c4)*l.f_ac_va*rho_8+(1-l.Y_c4)*l.f_ac_bu*rho_9+...                         
    (1 -l.Y_pro)*l.f_ac_pro*rho_10-rho_11+(1-l.Y_lac_f)*0.215*rho_20+...          
    (1 -l.Y_lac_o)*2/3*rho_21;                                              
d_S_h2_dt=l.D *(l.S_h2_in-S_h2)+(1-l.Y_su)*l.f_h2_su*rho_5+...              %8(addedlac) 
    (1 -l.Y_aa)*l.f_h2_aa*rho_6+(1-l.Y_fa)*l.f_h2_fa*rho_7+...                     
    (1 -l.Y_c4)*l.f_h2_va*rho_8+(1-l.Y_c4)*l.f_h2_bu*rho_9+...                         
    (1 -l.Y_pro)*l.f_h2_pro*rho_10-rho_12-rho_T_8+(1-l.Y_lac_o)*1/3*rho_21;       

% .......................................................................... 
%Differential equations 9-12, soluble matter
%..........................................................................StateNo. 
% d_S_ch4_dt=l.q_in/l.V_liq*(l.S_ch4_in-S_ch4)+(1-l.Y_ac)*rho_11+...          %9
%     (1-l.Y_h2)*rho_12-rho_T_9;                                              %? 
% d_S_IC_dt=l.q_in/l.V_liq *(l.S_IC_in-S_IC)-(s_1*rho_1+s_2*rho_2+...         %10(addedlac) 
%     s_3 *rho_3+s_4*rho_4+s_5*rho_5+s_6*rho_6+s_7*rho_7+...                  %? 
%     s_8 *rho_8+s_9*rho_9+s_10*rho_10+s_11*rho_11+s_12*rho_12+...            %?
%     s_20 *rho_20+s_21*rho_21+s_13*(rho_13+rho_14+rho_15+rho_16...           %? 
%     +rho_17+rho_18+rho_19+rho_22+rho_23))-rho_T_10;                         %? 
% d_S_IN_dt=l.q_in/l.V_liq *(l.S_IN_in-S_IN)-l.Y_su*l.N_bac*rho_5+...         %11(addedlac) 
%     (l.N_aa -l.Y_aa*l.N_bac)*rho_6-l.Y_fa*l.N_bac*rho_7-...                 %? 
%     l.Y_c4 *l.N_bac*rho_8-l.Y_c4*l.N_bac*rho_9-...                          %? 
%     l.Y_pro *l.N_bac*rho_10-l.Y_ac*l.N_bac*rho_11-...                       %?
%     l.Y_h2 *l.N_bac*rho_12-l.Y_lac_f*l.N_bac*rho_20-...                     %? 
%     l.Y_lac_o *l.N_bac*rho_21+(l.N_bac-l.N_xc)*...                          %? 
%     (rho_13+rho_14+rho_15+rho_16+rho_17+rho_18+...                          %?
%     rho_19+rho_22+rho_23)+...                                               %?
% (l.N_xc-l.f_xI_xc*l.N_I-l.f_sI_xc*l.N_I-l.f_pr_xc*l.N_aa)*rho_1;            %? 
% d_S_I_dt=l.q_in/l.V_liq *(l.S_I_in-S_I)+l.f_sI_xc*rho_1;                    %12

d_S_ch4_dt=l.D*(l.S_ch4_in-S_ch4)+(1-l.Y_ac)*rho_11+...          %9
    (1-l.Y_h2)*rho_12-rho_T_9;                                              %? 
d_S_IC_dt=l.D *(l.S_IC_in-S_IC)-(s_1*rho_1+s_2*rho_2+...         %10(addedlac) 
    s_3 *rho_3+s_4*rho_4+s_5*rho_5+s_6*rho_6+s_7*rho_7+...                  %? 
    s_8 *rho_8+s_9*rho_9+s_10*rho_10+s_11*rho_11+s_12*rho_12+...            %?
    s_20 *rho_20+s_21*rho_21+s_13*(rho_13+rho_14+rho_15+rho_16...           %? 
    +rho_17+rho_18+rho_19+rho_22+rho_23))-rho_T_10;                         %? 
% d_S_IC_dt=l.D *(l.S_IC_in-S_IC)-(s_1*rho_1+s_2*rho_2+...         %10(addedlac) 
%     s_3 *rho_3+s_4*rho_4+s_5*rho_5+s_6*rho_6+s_7*rho_7+...                  %? 
%     s_8 *rho_8+s_9*rho_9+s_10*rho_10+s_11*rho_11+s_12*rho_12+...            %?
%     s_13*(rho_13+rho_14+rho_15+rho_16...           %? 
%     +rho_17+rho_18+rho_19+rho_22+rho_23))-rho_T_10; 
d_S_IN_dt=l.D *(l.S_IN_in-S_IN)-l.Y_su*l.N_bac*rho_5+...         %11(addedlac) 
    (l.N_aa -l.Y_aa*l.N_bac)*rho_6-l.Y_fa*l.N_bac*rho_7-...                 %? 
    l.Y_c4 *l.N_bac*rho_8-l.Y_c4*l.N_bac*rho_9-...                          %? 
    l.Y_pro *l.N_bac*rho_10-l.Y_ac*l.N_bac*rho_11-...                       %?
    l.Y_h2 *l.N_bac*rho_12-l.Y_lac_f*l.N_bac*rho_20-...                     %? 
    l.Y_lac_o *l.N_bac*rho_21+(l.N_bac-l.N_xc)*...                          %? 
    (rho_13+rho_14+rho_15+rho_16+rho_17+rho_18+...                          %?
    rho_19+rho_22+rho_23)+...                                               %?
(l.N_xc-l.f_xI_xc*l.N_I-l.f_sI_xc*l.N_I-l.f_pr_xc*l.N_aa)*rho_1;            %? 
d_S_I_dt=l.D *(l.S_I_in-S_I)+l.f_sI_xc*rho_1;                    %12

%.......................................................................... 
%Differential equations 13-16, particulate matter 
%..........................................................................StateNo. 
% d_X_c_dt=l.q_in/l.V_liq*(l.X_xc_in-X_c)-rho_1+...                           %13(addedlac) 
%     rho_13+rho_14+rho_15+rho_16+rho_17+rho_18+rho_19+rho_22+rho_23;         %? 
% d_X_ch_dt=l.q_in/l.V_liq *(l.X_ch_in-X_ch)+l.f_ch_xc*rho_1-rho_2;           %14
% d_X_pr_dt=l.q_in/l.V_liq *(l.X_pr_in-X_pr)+l.f_pr_xc*rho_1-rho_3;           %15
% d_X_li_dt=l.q_in/l.V_liq *(l.X_li_in-X_li)+l.f_li_xc*rho_1-rho_4;           %16

d_X_c_dt=l.D*(l.X_xc_in-X_c)-rho_1+...                           %13(addedlac) 
    rho_13+rho_14+rho_15+rho_16+rho_17+rho_18+rho_19+rho_22+rho_23;         %? 
d_X_ch_dt=l.D *(l.X_ch_in-X_ch)+l.f_ch_xc*rho_1-rho_2;           %14
d_X_pr_dt=l.D *(l.X_pr_in-X_pr)+l.f_pr_xc*rho_1-rho_3;           %15
d_X_li_dt=l.D *(l.X_li_in-X_li)+l.f_li_xc*rho_1-rho_4;           %16

% ..........................................................................
%Differential equations 17-20, particulate matter 
%..........................................................................StateNo. 
% d_X_su_dt=l.q_in/l.V_liq*(l.X_su_in-X_su)+l.Y_su*rho_5-rho_13;              %17
% d_X_aa_dt=l.q_in/l.V_liq*(l.X_aa_in-X_aa)+l.Y_aa*rho_6-rho_14;              %18 
% d_X_fa_dt=l.q_in/l.V_liq*(l.X_fa_in-X_fa)+l.Y_fa*rho_7-rho_15;              %19 
% d_X_c4_dt=l.q_in/l.V_liq*(l.X_c4_in-X_c4)+l.Y_c4*rho_8+l.Y_c4*rho_9-rho_16; %20 

d_X_su_dt=l.D*(l.X_su_in-X_su)+l.Y_su*rho_5-rho_13;              %17
d_X_aa_dt=l.D*(l.X_aa_in-X_aa)+l.Y_aa*rho_6-rho_14;              %18 
d_X_fa_dt=l.D*(l.X_fa_in-X_fa)+l.Y_fa*rho_7-rho_15;              %19 
d_X_c4_dt=l.D*(l.X_c4_in-X_c4)+l.Y_c4*rho_8+l.Y_c4*rho_9-rho_16; %20 

%.......................................................................... 
%Differential equations 21-24, particulate matter 
%..........................................................................StateNo. 
% d_X_pro_dt=l.q_in/l.V_liq*(l.X_pro_in-X_pro)+l.Y_pro*rho_10-rho_17;         %21 
% d_X_ac_dt=l.q_in/l.V_liq*(l.X_ac_in-X_ac)+l.Y_ac*rho_11-rho_18;             %22
% d_X_h2_dt=l.q_in/l.V_liq*(l.X_h2_in-X_h2)+l.Y_h2*rho_12-rho_19;             %23
% d_X_I_dt=l.q_in/l.V_liq*(l.X_I_in-X_I)+l.f_xI_xc*rho_1;                     %24

d_X_pro_dt=l.D*(l.X_pro_in-X_pro)+l.Y_pro*rho_10-rho_17;         %21 
d_X_ac_dt=l.D*(l.X_ac_in-X_ac)+l.Y_ac*rho_11-rho_18;             %22
d_X_h2_dt=l.D*(l.X_h2_in-X_h2)+l.Y_h2*rho_12-rho_19;             %23
d_X_I_dt=l.D*(l.X_I_in-X_I)+l.f_xI_xc*rho_1;  

%.......................................................................... 
%Differential equations 25-26, cations and anions 
%..........................................................................StateNo. 
% d_S_cat_dt=l.q_in/l.V_liq*(l.S_cat_in-S_cat);                               %25
% d_S_an_dt=l.q_in/l.V_liq*(l.S_an_in-S_an);                                  %26

d_S_cat_dt=l.D*(l.S_cat_in-S_cat);                               %25
d_S_an_dt=l.D*(l.S_an_in-S_an);                                  %26

%..........................................................................
%Differential equations 27-32, ion states
%..........................................................................StateNo. 
d_S_vam_dt=-rho_A_4;                                                        %27
d_S_bum_dt=-rho_A_5;                                                        %28
d_S_prom_dt=-rho_A_6;                                                       %29
d_S_acm_dt=-rho_A_7;                                                        %30 
d_S_hco3m_dt=-rho_A_10;                                                     %31 
d_S_nh3_dt=-rho_A_11;                                                       %32 
%.......................................................................... 
%Differential equations 33-35, gas phase equations 
%..........................................................................StateNo. 
d_S_gas_h2_dt=-S_gas_h2*q_gas/l.V_gas+rho_T_8*l.V_liq/l.V_gas;              %33
d_S_gas_ch4_dt=-S_gas_ch4*q_gas/l.V_gas+rho_T_9*l.V_liq/l.V_gas;            %34 
d_S_gas_co2_dt=-S_gas_co2*q_gas/l.V_gas+rho_T_10*l.V_liq/l.V_gas;           %35 
%.......................................................................... 
%Differential equations 36-38, new equations 
%..........................................................................StateNo. 
% d_S_lac_dt=l.q_in/l.V_liq*(l.S_lac_in-S_lac)-rho_20-rho_21;                 %36
% d_X_lac_f_dt=l.q_in/l.V_liq*(l.X_lac_f_in-X_lac_f)+l.Y_lac_f*rho_20-rho_22; %37 
% d_X_lac_o_dt=l.q_in/l.V_liq*(l.X_lac_o_in-X_lac_o)+l.Y_lac_o*rho_21-rho_23; %38 
% d_S_ca_dt=l.q_in/l.V_liq*(l.S_ca_in-S_ca)-rho_P_24;                         %39

d_S_lac_dt=l.D*(l.S_lac_in-S_lac)-rho_20-rho_21;                 %36
d_X_lac_f_dt=l.D*(l.X_lac_f_in-X_lac_f)+l.Y_lac_f*rho_20-rho_22; %37 
d_X_lac_o_dt=l.D*(l.X_lac_o_in-X_lac_o)+l.Y_lac_o*rho_21-rho_23; %38 
d_S_ca_dt=l.D*(l.S_ca_in-S_ca)-rho_P_24;                         %39

%output vector from function 
ADM1_dt=[d_S_su_dt;d_S_aa_dt;d_S_fa_dt;d_S_va_dt;d_S_bu_dt;... 
    d_S_pro_dt;d_S_ac_dt;d_S_h2_dt;d_S_ch4_dt;d_S_IC_dt;d_S_IN_dt;... 
    d_S_I_dt;d_X_c_dt;d_X_ch_dt;d_X_pr_dt;d_X_li_dt;d_X_su_dt;... 
    d_X_aa_dt;d_X_fa_dt;d_X_c4_dt;d_X_pro_dt;d_X_ac_dt;d_X_h2_dt;... 
    d_X_I_dt;d_S_cat_dt;d_S_an_dt;d_S_vam_dt;d_S_bum_dt;d_S_prom_dt;... 
    d_S_acm_dt;d_S_hco3m_dt;d_S_nh3_dt;d_S_gas_h2_dt;d_S_gas_ch4_dt;... 
    d_S_gas_co2_dt;d_S_lac_dt;d_X_lac_f_dt;d_X_lac_o_dt;d_S_ca_dt]; 
end

