function []=plot_resids_ME_indSeis_201711()

%
close all
system('rm vals_residME_indSeis.mat')

% plot logicals
plot_bias_phi_tau=1
plot_within_event=1
plot_within_event_vs30=1
plot_between_event=1
% plot_within_event_together=0

% read data files
if exist('vals_residME_indSeis.mat','file')
  disp('Loading existing mat-file, vals_residME_indSeis.mat ...')
  load vals_residME_indSeis.mat
else
  disp('NGAW')
  [dist_ngaw2,mag_ngaw2,depth_ngaw2,per_arr_ngaw2,resid_c_ci95l_ngaw2,resid_c_ngaw2,resid_c_ci95u_ngaw2,resid_phi_ci95l_ngaw2,resid_phi_ngaw2,resid_phi_ci95u_ngaw2,resid_tau_ci95l_ngaw2,resid_tau_ngaw2,resid_tau_ci95u_ngaw2,resid_intra_ngaw2,resid_inter_ngaw2,vs30_ngaw2]=read_data_NGAW2();
  disp('NGAE1')
  [dist_ngaep01,mag_ngaep01,depth_ngaep01,per_arr_ngaep01,resid_c_ci95l_ngaep01,resid_c_ngaep01,resid_c_ci95u_ngaep01,resid_phi_ci95l_ngaep01,resid_phi_ngaep01,resid_phi_ci95u_ngaep01,resid_tau_ci95l_ngaep01,resid_tau_ngaep01,resid_tau_ci95u_ngaep01,resid_intra_ngaep01,resid_inter_ngaep01,vs30_ngaep01]=read_data_NGAEp01();
  disp('NGAE2')
  [dist_ngaep03,mag_ngaep03,depth_ngaep03,per_arr_ngaep03,resid_c_ci95l_ngaep03,resid_c_ngaep03,resid_c_ci95u_ngaep03,resid_phi_ci95l_ngaep03,resid_phi_ngaep03,resid_phi_ci95u_ngaep03,resid_tau_ci95l_ngaep03,resid_tau_ngaep03,resid_tau_ci95u_ngaep03,resid_intra_ngaep03,resid_inter_ngaep03,vs30_ngaep03]=read_data_NGAEp03();
  disp('CEUS')
  [dist_ceus,mag_ceus,depth_ceus,per_arr_ceus,resid_c_ci95l_ceus,resid_c_ceus,resid_c_ci95u_ceus,resid_phi_ci95l_ceus,resid_phi_ceus,resid_phi_ci95u_ceus,resid_tau_ci95l_ceus,resid_tau_ceus,resid_tau_ci95u_ceus,resid_intra_ceus,resid_inter_ceus,vs30_ceus]=read_data_CEUS();
  save vals_residME_indSeis
end

% read ME results
% bias terms
%ff=csvread('resid_WasM7_ASK14_All_lme.csv',1,0);
% no evid column (e.g., hypoC_ax13_seed103_sv1_rv60)
% ASK14: phi, tau
% s1
  phi_ASK14_smM=[.578 .555 .548 .527 .505 .457 .429];
  tau_ASK14_smM=[.47 .47 .47 .47 .47 .47 .47];
% s2
  phi_ASK14=[.64 .65 .64 .63 .63 .63 .63];
  tau_ASK14=[.36 .36 .36 .36 .36 .36 .36];
  dist_ASK14=logspace(-2,2,50);
  for ii=1:length(dist_ASK14)
    s5_ASK14_J_1p5=0.80;
    s6_ASK14_J_1p5=0.66;
    s5_ASK14_J_2=0.80;
    s6_ASK14_J_2=0.62;
    s5_ASK14_J_3=0.80;
    s6_ASK14_J_3=0.55;
    s5_ASK14_J_4=0.76;
    s6_ASK14_J_4=0.52;
    s5_ASK14_J_5=0.72;
    s6_ASK14_J_5=0.5;
    s5_ASK14_J_7p5=0.67;
    s6_ASK14_J_7p5=0.5;
    s5_ASK14_J_10=0.64;
    s6_ASK14_J_10=0.5;
    if (dist_ASK14(ii)<30)
      phi_ASK14_J_1p5(ii)=s5_ASK14_J_1p5;
      phi_ASK14_J_2(ii)=s5_ASK14_J_2;
      phi_ASK14_J_3(ii)=s5_ASK14_J_3;
      phi_ASK14_J_4(ii)=s5_ASK14_J_4;
      phi_ASK14_J_5(ii)=s5_ASK14_J_5;
      phi_ASK14_J_7p5(ii)=s5_ASK14_J_7p5;
      phi_ASK14_J_10(ii)=s5_ASK14_J_10;
    elseif (dist_ASK14(ii)>=30) & (dist_ASK14(ii)<80)
      phi_ASK14_J_1p5(ii)=s5_ASK14_J_1p5+(s6_ASK14_J_1p5-s5_ASK14_J_1p5)*(dist_ASK14(ii)-30)/50;
      phi_ASK14_J_2(ii)=s5_ASK14_J_2+(s6_ASK14_J_2-s5_ASK14_J_2)*(dist_ASK14(ii)-30)/50;
      phi_ASK14_J_3(ii)=s5_ASK14_J_3+(s6_ASK14_J_3-s5_ASK14_J_3)*(dist_ASK14(ii)-30)/50;
      phi_ASK14_J_4(ii)=s5_ASK14_J_4+(s6_ASK14_J_4-s5_ASK14_J_4)*(dist_ASK14(ii)-30)/50;
      phi_ASK14_J_5(ii)=s5_ASK14_J_5+(s6_ASK14_J_5-s5_ASK14_J_5)*(dist_ASK14(ii)-30)/50;
      phi_ASK14_J_7p5(ii)=s5_ASK14_J_7p5+(s6_ASK14_J_7p5-s5_ASK14_J_7p5)*(dist_ASK14(ii)-30)/50;
      phi_ASK14_J_10(ii)=s5_ASK14_J_10+(s6_ASK14_J_10-s5_ASK14_J_10)*(dist_ASK14(ii)-30)/50;
    else
      phi_ASK14_J_1p5(ii)=s6_ASK14_J_1p5;
      phi_ASK14_J_2(ii)=s6_ASK14_J_2;
      phi_ASK14_J_3(ii)=s6_ASK14_J_3;
      phi_ASK14_J_4(ii)=s6_ASK14_J_4;
      phi_ASK14_J_5(ii)=s6_ASK14_J_5;
      phi_ASK14_J_7p5(ii)=s6_ASK14_J_7p5;
      phi_ASK14_J_10(ii)=s6_ASK14_J_10;
    end

  end

% PLOTTING
if plot_bias_phi_tau
% FIG1
figure
%bias
subplot(1,3,1)
plot(per_arr_ngaw2,resid_c_ci95l_ngaw2,'b-'), hold on 
plot(per_arr_ngaw2,resid_c_ci95u_ngaw2,'b-'), 
plot(per_arr_ngaep01,resid_c_ci95l_ngaep01,'r-'), 
plot(per_arr_ngaep01,resid_c_ci95u_ngaep01,'r-'), 
plot(per_arr_ngaep03,resid_c_ci95l_ngaep03,'r--'), 
plot(per_arr_ngaep03,resid_c_ci95u_ngaep03,'r--'), 
plot(per_arr_ceus,resid_c_ci95l_ceus,'c-'), 
plot(per_arr_ceus,resid_c_ci95u_ceus,'c-'), 
set(gca,'XScale','log')
title('bias')
xlabel('T (s)')
%set(gca,'YLim',[-0.2 0.8])
axis square
subplot(1,3,2)
plot(per_arr_ngaw2,resid_tau_ci95l_ngaw2,'b-'), hold on 
plot(per_arr_ngaw2,resid_tau_ci95u_ngaw2,'b-'), 
plot(per_arr_ngaep01,resid_tau_ci95l_ngaep01,'r-'), 
plot(per_arr_ngaep01,resid_tau_ci95u_ngaep01,'r-'), 
plot(per_arr_ngaep03,resid_tau_ci95l_ngaep03,'r--'), 
plot(per_arr_ngaep03,resid_tau_ci95u_ngaep03,'r--'), 
plot(per_arr_ceus,resid_tau_ci95l_ceus,'c-'), 
plot(per_arr_ceus,resid_tau_ci95u_ceus,'c-'), 
set(gca,'XScale','log')
title('\tau')
xlabel('T (s)')
%set(gca,'YLim',[-0.2 0.8])
axis square
subplot(1,3,3)
plot(per_arr_ngaw2,resid_phi_ci95l_ngaw2,'b-'), hold on 
plot(per_arr_ngaw2,resid_phi_ci95u_ngaw2,'b-'), 
plot(per_arr_ngaep01,resid_phi_ci95l_ngaep01,'r-'), 
plot(per_arr_ngaep01,resid_phi_ci95u_ngaep01,'r-'), 
plot(per_arr_ngaep03,resid_phi_ci95l_ngaep03,'r--'), 
plot(per_arr_ngaep03,resid_phi_ci95u_ngaep03,'r--'), 
plot(per_arr_ceus,resid_phi_ci95l_ceus,'c-'), 
plot(per_arr_ceus,resid_phi_ci95u_ceus,'c-'), 
set(gca,'XScale','log')
title('\phi')
xlabel('T (s)')
%set(gca,'YLim',[-0.2 0.8])
axis square
end
%..............................

% inter-event
if plot_between_event
figure
subplot(2,2,1)
plot1inter(mag_ngaw2,resid_inter_ngaw2)
title('NGAW2, between event')
ylabel('resid, inter')
xlabel('M')
subplot(2,2,2)
plot1inter(mag_ngaep01,resid_inter_ngaep01)
title('NGAE,p01')
ylabel('resid, inter')
xlabel('M')
subplot(2,2,3)
plot1inter(mag_ngaep03,resid_inter_ngaep03)
title('NGAE,p03')
ylabel('resid, inter')
subplot(2,2,4)
plot1inter(mag_ceus,resid_inter_ceus)
title('CEUS')
ylabel('resid, inter')
ylabel('resid, inter')
xlabel('M')
end
%..............................

%..............................
% intra-event
if plot_within_event
[dist_bins_0p1_ngaw2,cnt_distBins_0p1_ngaw2, mean_distBins_0p1_ngaw2,var_distBins_0p1_ngaw2,mean_distBins_sm_0p1_ngaw2,var_distBins_sm_0p1_ngaw2]=bin_logData_by_dist(dist_ngaw2.t0p1,resid_intra_ngaw2.t0p1);
[dist_bins_0p2_ngaw2,cnt_distBins_0p2_ngaw2, mean_distBins_0p2_ngaw2,var_distBins_0p2_ngaw2,mean_distBins_sm_0p2_ngaw2,var_distBins_sm_0p2_ngaw2]=bin_logData_by_dist(dist_ngaw2.t0p2,resid_intra_ngaw2.t0p2);
[dist_bins_0p3_ngaw2,cnt_distBins_0p3_ngaw2, mean_distBins_0p3_ngaw2,var_distBins_0p3_ngaw2,mean_distBins_sm_0p3_ngaw2,var_distBins_sm_0p3_ngaw2]=bin_logData_by_dist(dist_ngaw2.t0p3,resid_intra_ngaw2.t0p3);
[dist_bins_0p5_ngaw2,cnt_distBins_0p5_ngaw2, mean_distBins_0p5_ngaw2,var_distBins_0p5_ngaw2,mean_distBins_sm_0p5_ngaw2,var_distBins_sm_0p5_ngaw2]=bin_logData_by_dist(dist_ngaw2.t0p5,resid_intra_ngaw2.t0p5);
[dist_bins_1p0_ngaw2,cnt_distBins_1p0_ngaw2, mean_distBins_1p0_ngaw2,var_distBins_1p0_ngaw2,mean_distBins_sm_1p0_ngaw2,var_distBins_sm_1p0_ngaw2]=bin_logData_by_dist(dist_ngaw2.t1p0,resid_intra_ngaw2.t1p0);
[dist_bins_2p0_ngaw2,cnt_distBins_2p0_ngaw2, mean_distBins_2p0_ngaw2,var_distBins_2p0_ngaw2,mean_distBins_sm_2p0_ngaw2,var_distBins_sm_2p0_ngaw2]=bin_logData_by_dist(dist_ngaw2.t2p0,resid_intra_ngaw2.t2p0);
[dist_bins_3p0_ngaw2,cnt_distBins_3p0_ngaw2, mean_distBins_3p0_ngaw2,var_distBins_3p0_ngaw2,mean_distBins_sm_3p0_ngaw2,var_distBins_sm_3p0_ngaw2]=bin_logData_by_dist(dist_ngaw2.t3p0,resid_intra_ngaw2.t3p0);
% NGAE,p01
[dist_bins_0p1_ngaep01,cnt_distBins_0p1_ngaep01, mean_distBins_0p1_ngaep01,var_distBins_0p1_ngaep01,mean_distBins_sm_0p1_ngaep01,var_distBins_sm_0p1_ngaep01]=bin_logData_by_dist(dist_ngaep01.t0p1,resid_intra_ngaep01.t0p1);
[dist_bins_0p2_ngaep01,cnt_distBins_0p2_ngaep01, mean_distBins_0p2_ngaep01,var_distBins_0p2_ngaep01,mean_distBins_sm_0p2_ngaep01,var_distBins_sm_0p2_ngaep01]=bin_logData_by_dist(dist_ngaep01.t0p2,resid_intra_ngaep01.t0p2);
[dist_bins_0p3_ngaep01,cnt_distBins_0p3_ngaep01, mean_distBins_0p3_ngaep01,var_distBins_0p3_ngaep01,mean_distBins_sm_0p3_ngaep01,var_distBins_sm_0p3_ngaep01]=bin_logData_by_dist(dist_ngaep01.t0p3,resid_intra_ngaep01.t0p3);
[dist_bins_0p5_ngaep01,cnt_distBins_0p5_ngaep01, mean_distBins_0p5_ngaep01,var_distBins_0p5_ngaep01,mean_distBins_sm_0p5_ngaep01,var_distBins_sm_0p5_ngaep01]=bin_logData_by_dist(dist_ngaep01.t0p5,resid_intra_ngaep01.t0p5);
[dist_bins_1p0_ngaep01,cnt_distBins_1p0_ngaep01, mean_distBins_1p0_ngaep01,var_distBins_1p0_ngaep01,mean_distBins_sm_1p0_ngaep01,var_distBins_sm_1p0_ngaep01]=bin_logData_by_dist(dist_ngaep01.t1p0,resid_intra_ngaep01.t1p0);
[dist_bins_2p0_ngaep01,cnt_distBins_2p0_ngaep01, mean_distBins_2p0_ngaep01,var_distBins_2p0_ngaep01,mean_distBins_sm_2p0_ngaep01,var_distBins_sm_2p0_ngaep01]=bin_logData_by_dist(dist_ngaep01.t2p0,resid_intra_ngaep01.t2p0);
[dist_bins_3p0_ngaep01,cnt_distBins_3p0_ngaep01, mean_distBins_3p0_ngaep01,var_distBins_3p0_ngaep01,mean_distBins_sm_3p0_ngaep01,var_distBins_sm_3p0_ngaep01]=bin_logData_by_dist(dist_ngaep01.t3p0,resid_intra_ngaep01.t3p0);
% NGAE,p03
[dist_bins_0p1_ngaep03,cnt_distBins_0p1_ngaep03, mean_distBins_0p1_ngaep03,var_distBins_0p1_ngaep03,mean_distBins_sm_0p1_ngaep03,var_distBins_sm_0p1_ngaep03]=bin_logData_by_dist(dist_ngaep03.t0p1,resid_intra_ngaep03.t0p1);
[dist_bins_0p2_ngaep03,cnt_distBins_0p2_ngaep03, mean_distBins_0p2_ngaep03,var_distBins_0p2_ngaep03,mean_distBins_sm_0p2_ngaep03,var_distBins_sm_0p2_ngaep03]=bin_logData_by_dist(dist_ngaep03.t0p2,resid_intra_ngaep03.t0p2);
[dist_bins_0p3_ngaep03,cnt_distBins_0p3_ngaep03, mean_distBins_0p3_ngaep03,var_distBins_0p3_ngaep03,mean_distBins_sm_0p3_ngaep03,var_distBins_sm_0p3_ngaep03]=bin_logData_by_dist(dist_ngaep03.t0p3,resid_intra_ngaep03.t0p3);
[dist_bins_0p5_ngaep03,cnt_distBins_0p5_ngaep03, mean_distBins_0p5_ngaep03,var_distBins_0p5_ngaep03,mean_distBins_sm_0p5_ngaep03,var_distBins_sm_0p5_ngaep03]=bin_logData_by_dist(dist_ngaep03.t0p5,resid_intra_ngaep03.t0p5);
[dist_bins_1p0_ngaep03,cnt_distBins_1p0_ngaep03, mean_distBins_1p0_ngaep03,var_distBins_1p0_ngaep03,mean_distBins_sm_1p0_ngaep03,var_distBins_sm_1p0_ngaep03]=bin_logData_by_dist(dist_ngaep03.t1p0,resid_intra_ngaep03.t1p0);
[dist_bins_2p0_ngaep03,cnt_distBins_2p0_ngaep03, mean_distBins_2p0_ngaep03,var_distBins_2p0_ngaep03,mean_distBins_sm_2p0_ngaep03,var_distBins_sm_2p0_ngaep03]=bin_logData_by_dist(dist_ngaep03.t2p0,resid_intra_ngaep03.t2p0);
[dist_bins_3p0_ngaep03,cnt_distBins_3p0_ngaep03, mean_distBins_3p0_ngaep03,var_distBins_3p0_ngaep03,mean_distBins_sm_3p0_ngaep03,var_distBins_sm_3p0_ngaep03]=bin_logData_by_dist(dist_ngaep03.t3p0,resid_intra_ngaep03.t3p0);
% CEUS
[dist_bins_0p1_ceus,cnt_distBins_0p1_ceus, mean_distBins_0p1_ceus,var_distBins_0p1_ceus,mean_distBins_sm_0p1_ceus,var_distBins_sm_0p1_ceus]=bin_logData_by_dist(dist_ceus.t0p1,resid_intra_ceus.t0p1);
[dist_bins_0p2_ceus,cnt_distBins_0p2_ceus, mean_distBins_0p2_ceus,var_distBins_0p2_ceus,mean_distBins_sm_0p2_ceus,var_distBins_sm_0p2_ceus]=bin_logData_by_dist(dist_ceus.t0p2,resid_intra_ceus.t0p2);
[dist_bins_0p3_ceus,cnt_distBins_0p3_ceus, mean_distBins_0p3_ceus,var_distBins_0p3_ceus,mean_distBins_sm_0p3_ceus,var_distBins_sm_0p3_ceus]=bin_logData_by_dist(dist_ceus.t0p3,resid_intra_ceus.t0p3);
[dist_bins_0p5_ceus,cnt_distBins_0p5_ceus, mean_distBins_0p5_ceus,var_distBins_0p5_ceus,mean_distBins_sm_0p5_ceus,var_distBins_sm_0p5_ceus]=bin_logData_by_dist(dist_ceus.t0p5,resid_intra_ceus.t0p5);
[dist_bins_1p0_ceus,cnt_distBins_1p0_ceus, mean_distBins_1p0_ceus,var_distBins_1p0_ceus,mean_distBins_sm_1p0_ceus,var_distBins_sm_1p0_ceus]=bin_logData_by_dist(dist_ceus.t1p0,resid_intra_ceus.t1p0);
[dist_bins_2p0_ceus,cnt_distBins_2p0_ceus, mean_distBins_2p0_ceus,var_distBins_2p0_ceus,mean_distBins_sm_2p0_ceus,var_distBins_sm_2p0_ceus]=bin_logData_by_dist(dist_ceus.t2p0,resid_intra_ceus.t2p0);
%
if plot_within_event_vs30
figure
subplot(2,3,1)
plot(vs30_ngaw2.t0p1,resid_intra_ngaw2.t0p1,'bs'), hold on
%set(gca,'XScale','log')
%set(gca,'XLim',[3 300]),
title('intra-event, NGA-W2 0.1 s')
xlabel('Vs30 (m/s)')
ylabel('mean(intra)')
subplot(2,3,2)
plot(vs30_ngaw2.t0p2,resid_intra_ngaw2.t0p2,'bs'), hold on
%set(gca,'XScale','log')
%set(gca,'XLim',[3 300]),
title('0.2 s')
xlabel('Vs30 (m/s)')
ylabel('mean(intra)')
subplot(2,3,3)
plot(vs30_ngaw2.t0p3,resid_intra_ngaw2.t0p3,'bs'), hold on
%set(gca,'XScale','log')
%set(gca,'XLim',[3 300]),
title('0.3 s')
xlabel('Vs30 (m/s)')
ylabel('mean(intra)')
subplot(2,3,4)
plot(vs30_ngaw2.t1p0,resid_intra_ngaw2.t1p0,'bs'), hold on
%set(gca,'XScale','log')
%set(gca,'XLim',[3 300]),
title('1.0 s')
xlabel('Vs30 (m/s)')
ylabel('mean(intra)')
subplot(2,3,5)
plot(vs30_ngaw2.t2p0,resid_intra_ngaw2.t2p0,'bs'), hold on
%set(gca,'XScale','log')
%set(gca,'XLim',[3 300]),
title('2.0 s')
%xlabel('R (km)')
xlabel('Vs30 (m/s)')
ylabel('mean(intra)')
subplot(2,3,6)
plot(vs30_ngaw2.t3p0,resid_intra_ngaw2.t3p0,'bs'), hold on
title('3.0 s')
xlabel('Vs30 (m/s)')
ylabel('mean(intra)')
%set(gca,'XScale','log')
%set(gca,'XLim',[3 300]),
end % end within event vs30

%
figure
subplot(2,3,1)
plot(dist_ngaw2.t0p1,resid_intra_ngaw2.t0p1,'bs'), hold on
plot(dist_bins_0p1_ngaw2,mean_distBins_0p1_ngaw2,'k-')
set(gca,'XScale','log')
set(gca,'XLim',[3 300]),
title('intra-event, NGA-W2 0.1 s')
xlabel('R (km)')
ylabel('mean(intra)')
subplot(2,3,2)
plot(dist_ngaw2.t0p2,resid_intra_ngaw2.t0p2,'bs'), hold on
plot(dist_bins_0p2_ngaw2,mean_distBins_0p2_ngaw2,'k-')
set(gca,'XScale','log')
set(gca,'XLim',[3 300]),
title('0.2 s')
xlabel('R (km)')
ylabel('mean(intra)')
subplot(2,3,3)
plot(dist_ngaw2.t0p3,resid_intra_ngaw2.t0p3,'bs'), hold on
plot(dist_bins_0p3_ngaw2,mean_distBins_0p3_ngaw2,'k-')
set(gca,'XScale','log')
set(gca,'XLim',[3 300]),
title('0.3 s')
xlabel('R (km)')
ylabel('mean(intra)')
subplot(2,3,4)
plot(dist_ngaw2.t1p0,resid_intra_ngaw2.t1p0,'bs'), hold on
plot(dist_bins_1p0_ngaw2,mean_distBins_1p0_ngaw2,'k-')
set(gca,'XScale','log')
set(gca,'XLim',[3 300]),
title('1.0 s')
xlabel('R (km)')
ylabel('mean(intra)')
subplot(2,3,5)
plot(dist_ngaw2.t2p0,resid_intra_ngaw2.t2p0,'bs'), hold on
plot(dist_bins_2p0_ngaw2,mean_distBins_2p0_ngaw2,'k-')
set(gca,'XScale','log')
set(gca,'XLim',[3 300]),
title('2.0 s')
xlabel('R (km)')
ylabel('mean(intra)')
subplot(2,3,6)
plot(dist_ngaw2.t3p0,resid_intra_ngaw2.t3p0,'bs'), hold on
plot(dist_bins_3p0_ngaw2,mean_distBins_3p0_ngaw2,'k-')
title('3.0 s')
xlabel('R (km)')
ylabel('mean(intra)')
set(gca,'XScale','log')
set(gca,'XLim',[3 300]),
% combine intra-event resid
figure
subplot(2,3,1)
plot(dist_bins_0p1_ngaw2,mean_distBins_0p1_ngaw2,'b-'), hold on
plot(dist_bins_0p1_ngaep01,mean_distBins_0p1_ngaep01,'c-')
plot(dist_bins_0p1_ngaep03,mean_distBins_0p1_ngaep03,'c--')
plot(dist_bins_0p1_ceus,mean_distBins_0p1_ceus,'m-')
legend('NGA-W2','NGAE,0.1','NGAE,0.3','CEUS')
set(gca,'XScale','log')
set(gca,'XLim',[3 300]),
title('intra, all, 0.1 s')
xlabel('R (km)')
ylabel('mean(intra)')
subplot(2,3,2)
plot(dist_bins_0p2_ngaw2,mean_distBins_0p2_ngaw2,'b-'), hold on
plot(dist_bins_0p2_ngaep01,mean_distBins_0p2_ngaep01,'c-')
plot(dist_bins_0p2_ngaep03,mean_distBins_0p2_ngaep03,'c--')
plot(dist_bins_0p2_ceus,mean_distBins_0p2_ceus,'m-')
set(gca,'XScale','log')
set(gca,'XLim',[3 300]),
xlabel('R (km)')
title('0.2 s')
ylabel('mean(intra)')
subplot(2,3,3)
plot(dist_bins_0p3_ngaw2,mean_distBins_0p3_ngaw2,'b-'), hold on
plot(dist_bins_0p3_ngaep01,mean_distBins_0p3_ngaep01,'c-')
plot(dist_bins_0p3_ngaep03,mean_distBins_0p3_ngaep03,'c--')
plot(dist_bins_0p3_ceus,mean_distBins_0p3_ceus,'m-')
set(gca,'XScale','log')
set(gca,'XLim',[3 300]),
xlabel('R (km)')
title('0.3 s')
ylabel('mean(intra)')
subplot(2,3,4)
plot(dist_bins_1p0_ngaw2,mean_distBins_1p0_ngaw2,'b-'), hold on
plot(dist_bins_1p0_ngaep01,mean_distBins_1p0_ngaep01,'c-')
plot(dist_bins_1p0_ngaep03,mean_distBins_1p0_ngaep03,'c--')
plot(dist_bins_1p0_ceus,mean_distBins_1p0_ceus,'m-')
set(gca,'XScale','log')
set(gca,'XLim',[3 300]),
xlabel('R (km)')
ylabel('mean(intra)')
title('1.0 s')
subplot(2,3,5)
plot(dist_bins_2p0_ngaw2,mean_distBins_2p0_ngaw2,'b-'), hold on
plot(dist_bins_2p0_ngaep01,mean_distBins_2p0_ngaep01,'c-')
plot(dist_bins_2p0_ngaep03,mean_distBins_2p0_ngaep03,'c--')
plot(dist_bins_2p0_ceus,mean_distBins_2p0_ceus,'m-')
set(gca,'XScale','log')
set(gca,'XLim',[3 300]),
title('2.0 s')
xlabel('R (km)')
ylabel('mean(intra)')
subplot(2,3,6)
plot(dist_bins_3p0_ngaw2,mean_distBins_3p0_ngaw2,'b-'), hold on
plot(dist_bins_3p0_ngaep01,mean_distBins_3p0_ngaep01,'c-')
plot(dist_bins_3p0_ngaep03,mean_distBins_3p0_ngaep03,'c--')
set(gca,'XScale','log')
set(gca,'XLim',[3 300]),
title('3.0 s')
xlabel('R (km)')
ylabel('mean(intra)')

% combining intra event NGA
figure
subplot(2,2,1),
plot(dist_bins_0p1_ngaw2,mean_distBins_0p1_ngaw2), hold on
plot(dist_bins_0p2_ngaw2,mean_distBins_0p2_ngaw2),
plot(dist_bins_0p5_ngaw2,mean_distBins_0p5_ngaw2),
plot(dist_bins_1p0_ngaw2,mean_distBins_1p0_ngaw2),
plot(dist_bins_2p0_ngaw2,mean_distBins_2p0_ngaw2),
plot(dist_bins_3p0_ngaw2,mean_distBins_3p0_ngaw2),
%legend('0.1','0.2','0.5','1.0','2.0','3.0')
xlabel('R (km)')
ylabel('mean(intra)')
title('NGAW2')
set(gca,'XScale','log')
set(gca,'XLim',[5 300]),
subplot(2,2,2),
plot(dist_bins_0p1_ngaep01,mean_distBins_0p1_ngaep01), hold on
plot(dist_bins_0p2_ngaep01,mean_distBins_0p2_ngaep01),
plot(dist_bins_0p5_ngaep01,mean_distBins_0p5_ngaep01),
plot(dist_bins_1p0_ngaep01,mean_distBins_1p0_ngaep01),
plot(dist_bins_2p0_ngaep01,mean_distBins_2p0_ngaep01),
plot(dist_bins_3p0_ngaep01,mean_distBins_3p0_ngaep01),
legend('0.1','0.2','0.5','1.0','2.0','3.0')
set(gca,'XScale','log')
set(gca,'XLim',[5 300]),
xlabel('R (km)')
ylabel('mean(intra)')
title('NGAE,0.1')
subplot(2,2,3),
plot(dist_bins_0p1_ngaep03,mean_distBins_0p1_ngaep03), hold on
plot(dist_bins_0p2_ngaep03,mean_distBins_0p2_ngaep03),
plot(dist_bins_0p5_ngaep03,mean_distBins_0p5_ngaep03),
plot(dist_bins_1p0_ngaep03,mean_distBins_1p0_ngaep03),
plot(dist_bins_2p0_ngaep03,mean_distBins_2p0_ngaep03),
plot(dist_bins_3p0_ngaep03,mean_distBins_3p0_ngaep03),
set(gca,'XScale','log')
set(gca,'XLim',[5 300]),
xlabel('R (km)')
ylabel('mean(intra)')
title('NGAE,0.3')
subplot(2,2,4),
plot(dist_bins_0p1_ceus,mean_distBins_0p1_ceus), hold on
plot(dist_bins_0p2_ceus,mean_distBins_0p2_ceus),
plot(dist_bins_0p5_ceus,mean_distBins_0p5_ceus),
plot(dist_bins_1p0_ceus,mean_distBins_1p0_ceus),
plot(dist_bins_2p0_ceus,mean_distBins_2p0_ceus),
set(gca,'XScale','log')
set(gca,'XLim',[5 300]),
xlabel('R (km)')
ylabel('mean(intra)')
title('CEUS,0.1')

%plot1intra(dist,resid_intra_ngaw2,1.5)
%subplot(2,3,2)
%plot1intra(rrup,intra_2,2)
%subplot(2,3,3)
%plot1intra(rrup,intra_3,3)
%subplot(2,3,4)
%plot1intra(rrup,intra_5,5)
%subplot(2,3,5)
%plot1intra(rrup,intra_7p5,7.5)
%subplot(2,3,6)
%plot1intra(rrup,intra_10,10)
end
%..............................

% intra-event, together
%if plot_within_event_together
%figure
%plot1intra_meanStd(rrup,intra_1p5,1.5)
%plot1intra_meanStd(rrup,intra_2,2)
%plot1intra_meanStd(rrup,intra_3,3)
%plot1intra_meanStd(rrup,intra_5,5)
%plot1intra_meanStd(rrup,intra_7p5,7.5)
%plot1intra_meanStd(rrup,intra_10,10)
%xlabel('Rrup (km)')
%title('intra, combined')
%ylabel('\delta{W_{es}}')
%end
%..............................

end
%-----------------------------------------------------

%-----------------------------------------------------
function [dist,mag,depth,per_arr,resid_c_ci95l,resid_c,resid_c_ci95u,resid_phi_ci95l,resid_phi,resid_phi_ci95u,resid_tau_ci95l,resid_tau,resid_tau_ci95u,resid_intra,resid_inter,vs30]=read_data_CEUS()

%
[mag_0p1,depth_0p1,dist_0p1,resid_c_ci95l_0p1,resid_c_0p1,resid_c_ci95u_0p1,resid_phi_ci95l_0p1,resid_phi_0p1,resid_phi_ci95u_0p1,resid_tau_ci95l_0p1,resid_tau_0p1,resid_tau_ci95u_0p1,resid_intra_0p1,resid_inter_0p1,vs30_0p1]=read_data_1csvCEUS('induced_lme_20171103/CEUS_resid_0p1_resid_0p1_lme.csv');
[mag_0p2,depth_0p2,dist_0p2,resid_c_ci95l_0p2,resid_c_0p2,resid_c_ci95u_0p2,resid_phi_ci95l_0p2,resid_phi_0p2,resid_phi_ci95u_0p2,resid_tau_ci95l_0p2,resid_tau_0p2,resid_tau_ci95u_0p2,resid_intra_0p2,resid_inter_0p2,vs30_0p2]=read_data_1csvCEUS('induced_lme_20171103/CEUS_resid_0p2_resid_0p2_lme.csv');
[mag_0p3,depth_0p3,dist_0p3,resid_c_ci95l_0p3,resid_c_0p3,resid_c_ci95u_0p3,resid_phi_ci95l_0p3,resid_phi_0p3,resid_phi_ci95u_0p3,resid_tau_ci95l_0p3,resid_tau_0p3,resid_tau_ci95u_0p3,resid_intra_0p3,resid_inter_0p3,vs30_0p3]=read_data_1csvCEUS('induced_lme_20171103/CEUS_resid_0p3_resid_0p3_lme.csv');
[mag_0p5,depth_0p5,dist_0p5,resid_c_ci95l_0p5,resid_c_0p5,resid_c_ci95u_0p5,resid_phi_ci95l_0p5,resid_phi_0p5,resid_phi_ci95u_0p5,resid_tau_ci95l_0p5,resid_tau_0p5,resid_tau_ci95u_0p5,resid_intra_0p5,resid_inter_0p5,vs30_0p5]=read_data_1csvCEUS('induced_lme_20171103/CEUS_resid_0p5_resid_0p5_lme.csv');
[mag_1p0,depth_1p0,dist_1p0,resid_c_ci95l_1p0,resid_c_1p0,resid_c_ci95u_1p0,resid_phi_ci95l_1p0,resid_phi_1p0,resid_phi_ci95u_1p0,resid_tau_ci95l_1p0,resid_tau_1p0,resid_tau_ci95u_1p0,resid_intra_1p0,resid_inter_1p0,vs30_1p0]=read_data_1csvCEUS('induced_lme_20171103/CEUS_resid_1p0_resid_1p0_lme.csv');
[mag_2p0,depth_2p0,dist_2p0,resid_c_ci95l_2p0,resid_c_2p0,resid_c_ci95u_2p0,resid_phi_ci95l_2p0,resid_phi_2p0,resid_phi_ci95u_2p0,resid_tau_ci95l_2p0,resid_tau_2p0,resid_tau_ci95u_2p0,resid_intra_2p0,resid_inter_2p0,vs30_2p0]=read_data_1csvCEUS('induced_lme_20171103/CEUS_resid_2p0_resid_2p0_lme.csv');
%
per_arr=[0.1 0.2 0.3 0.5 1.0 2.0];
resid_c_ci95l=[resid_c_ci95l_0p1 resid_c_ci95l_0p2 resid_c_ci95l_0p3 resid_c_ci95l_0p5 resid_c_ci95l_1p0 resid_c_ci95l_2p0];
resid_c=[resid_c_0p1 resid_c_0p2 resid_c_0p3 resid_c_0p5 resid_c_1p0 resid_c_2p0];
resid_c_ci95u=[resid_c_ci95u_0p1 resid_c_ci95u_0p2 resid_c_ci95u_0p3 resid_c_ci95u_0p5 resid_c_ci95u_1p0 resid_c_ci95u_2p0];
resid_phi_ci95l=[resid_phi_ci95l_0p1 resid_phi_ci95l_0p2 resid_phi_ci95l_0p3 resid_phi_ci95l_0p5 resid_phi_ci95l_1p0 resid_phi_ci95l_2p0 ];
resid_phi=[resid_phi_0p1 resid_phi_0p2 resid_phi_0p3 resid_phi_0p5 resid_phi_1p0 resid_phi_2p0 ];
resid_phi_ci95u=[resid_phi_ci95u_0p1 resid_phi_ci95u_0p2 resid_phi_ci95u_0p3 resid_phi_ci95u_0p5 resid_phi_ci95u_1p0 resid_phi_ci95u_2p0 ];
resid_tau_ci95l=[resid_tau_ci95l_0p1 resid_tau_ci95l_0p2 resid_tau_ci95l_0p3 resid_tau_ci95l_0p5 resid_tau_ci95l_1p0 resid_tau_ci95l_2p0 ];
resid_tau=[resid_tau_0p1 resid_tau_0p2 resid_tau_0p3 resid_tau_0p5 resid_tau_1p0 resid_tau_2p0 ];
resid_tau_ci95u=[resid_tau_ci95u_0p1 resid_tau_ci95u_0p2 resid_tau_ci95u_0p3 resid_tau_ci95u_0p5 resid_tau_ci95u_1p0 resid_tau_ci95u_2p0 ];
%
dist.t0p1=dist_0p1;
dist.t0p2=dist_0p2;
dist.t0p3=dist_0p3;
dist.t0p5=dist_0p5;
dist.t1p0=dist_1p0;
dist.t2p0=dist_2p0;
%
mag.t0p1=mag_0p1;
mag.t0p2=mag_0p2;
mag.t0p3=mag_0p3;
mag.t0p5=mag_0p5;
mag.t1p0=mag_1p0;
mag.t2p0=mag_2p0;
%
vs30.t0p1=vs30_0p1;
vs30.t0p2=vs30_0p2;
vs30.t0p3=vs30_0p3;
vs30.t0p5=vs30_0p5;
vs30.t1p0=vs30_1p0;
vs30.t2p0=vs30_2p0;
%
depth.t0p1=depth_0p1;
depth.t0p2=depth_0p2;
depth.t0p3=depth_0p3;
depth.t0p5=depth_0p5;
depth.t1p0=depth_1p0;
depth.t2p0=depth_2p0;
%
resid_intra.t0p1=resid_intra_0p1;
resid_intra.t0p2=resid_intra_0p2;
resid_intra.t0p3=resid_intra_0p3;
resid_intra.t0p5=resid_intra_0p5;
resid_intra.t1p0=resid_intra_1p0;
resid_intra.t2p0=resid_intra_2p0;
%
resid_inter.t0p1=resid_inter_0p1;
resid_inter.t0p2=resid_inter_0p2;
resid_inter.t0p3=resid_inter_0p3;
resid_inter.t0p5=resid_inter_0p5;
resid_inter.t1p0=resid_inter_1p0;
resid_inter.t2p0=resid_inter_2p0;


end
%-----------------------------------------------------


%-----------------------------------------------------
function [dist,mag,depth,per_arr,resid_c_ci95l,resid_c,resid_c_ci95u,resid_phi_ci95l,resid_phi,resid_phi_ci95u,resid_tau_ci95l,resid_tau,resid_tau_ci95u,resid_intra,resid_inter,vs30]=read_data_NGAEp03()

%
[mag_0p1,depth_0p1,dist_0p1,resid_c_ci95l_0p1,resid_c_0p1,resid_c_ci95u_0p1,resid_phi_ci95l_0p1,resid_phi_0p1,resid_phi_ci95u_0p1,resid_tau_ci95l_0p1,resid_tau_0p1,resid_tau_ci95u_0p1,resid_intra_0p1,resid_inter_0p1,vs30_0p1]=read_data_1csv('induced_lme_20171103/NGAE_p03_resid_0p1_resid_0p1_lme.csv');
[mag_0p2,depth_0p2,dist_0p2,resid_c_ci95l_0p2,resid_c_0p2,resid_c_ci95u_0p2,resid_phi_ci95l_0p2,resid_phi_0p2,resid_phi_ci95u_0p2,resid_tau_ci95l_0p2,resid_tau_0p2,resid_tau_ci95u_0p2,resid_intra_0p2,resid_inter_0p2,vs30_0p2]=read_data_1csv('induced_lme_20171103/NGAE_p03_resid_0p2_resid_0p2_lme.csv');
[mag_0p3,depth_0p3,dist_0p3,resid_c_ci95l_0p3,resid_c_0p3,resid_c_ci95u_0p3,resid_phi_ci95l_0p3,resid_phi_0p3,resid_phi_ci95u_0p3,resid_tau_ci95l_0p3,resid_tau_0p3,resid_tau_ci95u_0p3,resid_intra_0p3,resid_inter_0p3,vs30_0p3]=read_data_1csv('induced_lme_20171103/NGAE_p03_resid_0p3_resid_0p3_lme.csv');
[mag_0p4,depth_0p4,dist_0p4,resid_c_ci95l_0p4,resid_c_0p4,resid_c_ci95u_0p4,resid_phi_ci95l_0p4,resid_phi_0p4,resid_phi_ci95u_0p4,resid_tau_ci95l_0p4,resid_tau_0p4,resid_tau_ci95u_0p4,resid_intra_0p4,resid_inter_0p4,vs30_0p4]=read_data_1csv('induced_lme_20171103/NGAE_p03_resid_0p4_resid_0p4_lme.csv');
[mag_0p5,depth_0p5,dist_0p5,resid_c_ci95l_0p5,resid_c_0p5,resid_c_ci95u_0p5,resid_phi_ci95l_0p5,resid_phi_0p5,resid_phi_ci95u_0p5,resid_tau_ci95l_0p5,resid_tau_0p5,resid_tau_ci95u_0p5,resid_intra_0p5,resid_inter_0p5,vs30_0p5]=read_data_1csv('induced_lme_20171103/NGAE_p03_resid_0p5_resid_0p5_lme.csv');
[mag_1p0,depth_1p0,dist_1p0,resid_c_ci95l_1p0,resid_c_1p0,resid_c_ci95u_1p0,resid_phi_ci95l_1p0,resid_phi_1p0,resid_phi_ci95u_1p0,resid_tau_ci95l_1p0,resid_tau_1p0,resid_tau_ci95u_1p0,resid_intra_1p0,resid_inter_1p0,vs30_1p0]=read_data_1csv('induced_lme_20171103/NGAE_p03_resid_1p0_resid_1p0_lme.csv');
[mag_2p0,depth_2p0,dist_2p0,resid_c_ci95l_2p0,resid_c_2p0,resid_c_ci95u_2p0,resid_phi_ci95l_2p0,resid_phi_2p0,resid_phi_ci95u_2p0,resid_tau_ci95l_2p0,resid_tau_2p0,resid_tau_ci95u_2p0,resid_intra_2p0,resid_inter_2p0,vs30_2p0]=read_data_1csv('induced_lme_20171103/NGAE_p03_resid_2p0_resid_2p0_lme.csv');
[mag_3p0,depth_3p0,dist_3p0,resid_c_ci95l_3p0,resid_c_3p0,resid_c_ci95u_3p0,resid_phi_ci95l_3p0,resid_phi_3p0,resid_phi_ci95u_3p0,resid_tau_ci95l_3p0,resid_tau_3p0,resid_tau_ci95u_3p0,resid_intra_3p0,resid_inter_3p0,vs30_3p0]=read_data_1csv('induced_lme_20171103/NGAE_p03_resid_3p0_resid_3p0_lme.csv');
[mag_4p0,depth_4p0,dist_4p0,resid_c_ci95l_4p0,resid_c_4p0,resid_c_ci95u_4p0,resid_phi_ci95l_4p0,resid_phi_4p0,resid_phi_ci95u_4p0,resid_tau_ci95l_4p0,resid_tau_4p0,resid_tau_ci95u_4p0,resid_intra_4p0,resid_inter_4p0,vs30_4p0]=read_data_1csv('induced_lme_20171103/NGAE_p03_resid_4p0_resid_4p0_lme.csv');
[mag_5p0,depth_5p0,dist_5p0,resid_c_ci95l_5p0,resid_c_5p0,resid_c_ci95u_5p0,resid_phi_ci95l_5p0,resid_phi_5p0,resid_phi_ci95u_5p0,resid_tau_ci95l_5p0,resid_tau_5p0,resid_tau_ci95u_5p0,resid_intra_5p0,resid_inter_5p0,vs30_5p0]=read_data_1csv('induced_lme_20171103/NGAE_p03_resid_5p0_resid_5p0_lme.csv');
[mag_10p0,depth_10p0,dist_10p0,resid_c_ci95l_10p0,resid_c_10p0,resid_c_ci95u_10p0,resid_phi_ci95l_10p0,resid_phi_10p0,resid_phi_ci95u_10p0,resid_tau_ci95l_10p0,resid_tau_10p0,resid_tau_ci95u_10p0,resid_intra_10p0,resid_inter_10p0,vs30_10p0]=read_data_1csv('induced_lme_20171103/NGAE_p03_resid_10p0_resid_10p0_lme.csv');
%
per_arr=[0.1 0.2 0.3 0.4 0.5 1.0 2.0 3.0 4.0 5.0 10.0];
resid_c_ci95l=[resid_c_ci95l_0p1 resid_c_ci95l_0p2 resid_c_ci95l_0p3 resid_c_ci95l_0p4 resid_c_ci95l_0p5 resid_c_ci95l_1p0 resid_c_ci95l_2p0 resid_c_ci95l_3p0 resid_c_ci95l_4p0 resid_c_ci95l_5p0 resid_c_ci95l_10p0];
resid_c=[resid_c_0p1 resid_c_0p2 resid_c_0p3 resid_c_0p4 resid_c_0p5 resid_c_1p0 resid_c_2p0 resid_c_3p0 resid_c_4p0 resid_c_5p0 resid_c_10p0];
resid_c_ci95u=[resid_c_ci95u_0p1 resid_c_ci95u_0p2 resid_c_ci95u_0p3 resid_c_ci95u_0p4 resid_c_ci95u_0p5 resid_c_ci95u_1p0 resid_c_ci95u_2p0 resid_c_ci95u_3p0 resid_c_ci95u_4p0 resid_c_ci95u_5p0 resid_c_ci95u_10p0];
resid_phi_ci95l=[resid_phi_ci95l_0p1 resid_phi_ci95l_0p2 resid_phi_ci95l_0p3 resid_phi_ci95l_0p4 resid_phi_ci95l_0p5 resid_phi_ci95l_1p0 resid_phi_ci95l_2p0 resid_phi_ci95l_3p0 resid_phi_ci95l_4p0 resid_phi_ci95l_5p0 resid_phi_ci95l_10p0];
resid_phi=[resid_phi_0p1 resid_phi_0p2 resid_phi_0p3 resid_phi_0p4 resid_phi_0p5 resid_phi_1p0 resid_phi_2p0 resid_phi_3p0 resid_phi_4p0 resid_phi_5p0 resid_phi_10p0];
resid_phi_ci95u=[resid_phi_ci95u_0p1 resid_phi_ci95u_0p2 resid_phi_ci95u_0p3 resid_phi_ci95u_0p4 resid_phi_ci95u_0p5 resid_phi_ci95u_1p0 resid_phi_ci95u_2p0 resid_phi_ci95u_3p0 resid_phi_ci95u_4p0 resid_phi_ci95u_5p0 resid_phi_ci95u_10p0];
resid_tau_ci95l=[resid_tau_ci95l_0p1 resid_tau_ci95l_0p2 resid_tau_ci95l_0p3 resid_tau_ci95l_0p4 resid_tau_ci95l_0p5 resid_tau_ci95l_1p0 resid_tau_ci95l_2p0 resid_tau_ci95l_3p0 resid_tau_ci95l_4p0 resid_tau_ci95l_5p0 resid_tau_ci95l_10p0];
resid_tau=[resid_tau_0p1 resid_tau_0p2 resid_tau_0p3 resid_tau_0p4 resid_tau_0p5 resid_tau_1p0 resid_tau_2p0 resid_tau_3p0 resid_tau_4p0 resid_tau_5p0 resid_tau_10p0];
resid_tau_ci95u=[resid_tau_ci95u_0p1 resid_tau_ci95u_0p2 resid_tau_ci95u_0p3 resid_tau_ci95u_0p4 resid_tau_ci95u_0p5 resid_tau_ci95u_1p0 resid_tau_ci95u_2p0 resid_tau_ci95u_3p0 resid_tau_ci95u_4p0 resid_tau_ci95u_5p0 resid_tau_ci95u_10p0];
%
dist.t0p1=dist_0p1;
dist.t0p2=dist_0p2;
dist.t0p3=dist_0p3;
dist.t0p4=dist_0p4;
dist.t0p5=dist_0p5;
dist.t1p0=dist_1p0;
dist.t2p0=dist_2p0;
dist.t3p0=dist_3p0;
dist.t4p0=dist_4p0;
dist.t5p0=dist_5p0;
dist.t10p0=dist_10p0;
%
mag.t0p1=mag_0p1;
mag.t0p2=mag_0p2;
mag.t0p3=mag_0p3;
mag.t0p4=mag_0p4;
mag.t0p5=mag_0p5;
mag.t1p0=mag_1p0;
mag.t2p0=mag_2p0;
mag.t3p0=mag_3p0;
mag.t4p0=mag_4p0;
mag.t5p0=mag_5p0;
mag.t10p0=mag_10p0;
%
vs30.t0p1=vs30_0p1;
vs30.t0p2=vs30_0p2;
vs30.t0p3=vs30_0p3;
vs30.t0p4=vs30_0p4;
vs30.t0p5=vs30_0p5;
vs30.t1p0=vs30_1p0;
vs30.t2p0=vs30_2p0;
vs30.t3p0=vs30_3p0;
vs30.t4p0=vs30_4p0;
vs30.t5p0=vs30_5p0;
vs30.t10p0=vs30_10p0;
%
depth.t0p1=depth_0p1;
depth.t0p2=depth_0p2;
depth.t0p3=depth_0p3;
depth.t0p4=depth_0p4;
depth.t0p5=depth_0p5;
depth.t1p0=depth_1p0;
depth.t2p0=depth_2p0;
depth.t3p0=depth_3p0;
depth.t4p0=depth_4p0;
depth.t5p0=depth_5p0;
depth.t10p0=depth_10p0;
%
resid_intra.t0p1=resid_intra_0p1;
resid_intra.t0p2=resid_intra_0p2;
resid_intra.t0p3=resid_intra_0p3;
resid_intra.t0p4=resid_intra_0p4;
resid_intra.t0p5=resid_intra_0p5;
resid_intra.t1p0=resid_intra_1p0;
resid_intra.t2p0=resid_intra_2p0;
resid_intra.t3p0=resid_intra_3p0;
resid_intra.t4p0=resid_intra_4p0;
resid_intra.t5p0=resid_intra_5p0;
resid_intra.t10p0=resid_intra_10p0;
%
resid_inter.t0p1=resid_inter_0p1;
resid_inter.t0p2=resid_inter_0p2;
resid_inter.t0p3=resid_inter_0p3;
resid_inter.t0p4=resid_inter_0p4;
resid_inter.t0p5=resid_inter_0p5;
resid_inter.t1p0=resid_inter_1p0;
resid_inter.t2p0=resid_inter_2p0;
resid_inter.t3p0=resid_inter_3p0;
resid_inter.t4p0=resid_inter_4p0;
resid_inter.t5p0=resid_inter_5p0;
resid_inter.t10p0=resid_inter_10p0;


end
%-----------------------------------------------------

%-----------------------------------------------------
function [dist,mag,depth,per_arr,resid_c_ci95l,resid_c,resid_c_ci95u,resid_phi_ci95l,resid_phi,resid_phi_ci95u,resid_tau_ci95l,resid_tau,resid_tau_ci95u,resid_intra,resid_inter,vs30]=read_data_NGAEp01()

%
[mag_0p1,depth_0p1,dist_0p1,resid_c_ci95l_0p1,resid_c_0p1,resid_c_ci95u_0p1,resid_phi_ci95l_0p1,resid_phi_0p1,resid_phi_ci95u_0p1,resid_tau_ci95l_0p1,resid_tau_0p1,resid_tau_ci95u_0p1,resid_intra_0p1,resid_inter_0p1,vs30_0p1]=read_data_1csv('induced_lme_20171103/NGAE_p01_resid_0p1_resid_0p1_lme.csv');
[mag_0p2,depth_0p2,dist_0p2,resid_c_ci95l_0p2,resid_c_0p2,resid_c_ci95u_0p2,resid_phi_ci95l_0p2,resid_phi_0p2,resid_phi_ci95u_0p2,resid_tau_ci95l_0p2,resid_tau_0p2,resid_tau_ci95u_0p2,resid_intra_0p2,resid_inter_0p2,vs30_0p2]=read_data_1csv('induced_lme_20171103/NGAE_p01_resid_0p2_resid_0p2_lme.csv');
[mag_0p3,depth_0p3,dist_0p3,resid_c_ci95l_0p3,resid_c_0p3,resid_c_ci95u_0p3,resid_phi_ci95l_0p3,resid_phi_0p3,resid_phi_ci95u_0p3,resid_tau_ci95l_0p3,resid_tau_0p3,resid_tau_ci95u_0p3,resid_intra_0p3,resid_inter_0p3,vs30_0p3]=read_data_1csv('induced_lme_20171103/NGAE_p01_resid_0p3_resid_0p3_lme.csv');
[mag_0p4,depth_0p4,dist_0p4,resid_c_ci95l_0p4,resid_c_0p4,resid_c_ci95u_0p4,resid_phi_ci95l_0p4,resid_phi_0p4,resid_phi_ci95u_0p4,resid_tau_ci95l_0p4,resid_tau_0p4,resid_tau_ci95u_0p4,resid_intra_0p4,resid_inter_0p4,vs30_0p4]=read_data_1csv('induced_lme_20171103/NGAE_p01_resid_0p4_resid_0p4_lme.csv');
[mag_0p5,depth_0p5,dist_0p5,resid_c_ci95l_0p5,resid_c_0p5,resid_c_ci95u_0p5,resid_phi_ci95l_0p5,resid_phi_0p5,resid_phi_ci95u_0p5,resid_tau_ci95l_0p5,resid_tau_0p5,resid_tau_ci95u_0p5,resid_intra_0p5,resid_inter_0p5,vs30_0p5]=read_data_1csv('induced_lme_20171103/NGAE_p01_resid_0p5_resid_0p5_lme.csv');
[mag_1p0,depth_1p0,dist_1p0,resid_c_ci95l_1p0,resid_c_1p0,resid_c_ci95u_1p0,resid_phi_ci95l_1p0,resid_phi_1p0,resid_phi_ci95u_1p0,resid_tau_ci95l_1p0,resid_tau_1p0,resid_tau_ci95u_1p0,resid_intra_1p0,resid_inter_1p0,vs30_1p0]=read_data_1csv('induced_lme_20171103/NGAE_p01_resid_1p0_resid_1p0_lme.csv');
[mag_2p0,depth_2p0,dist_2p0,resid_c_ci95l_2p0,resid_c_2p0,resid_c_ci95u_2p0,resid_phi_ci95l_2p0,resid_phi_2p0,resid_phi_ci95u_2p0,resid_tau_ci95l_2p0,resid_tau_2p0,resid_tau_ci95u_2p0,resid_intra_2p0,resid_inter_2p0,vs30_2p0]=read_data_1csv('induced_lme_20171103/NGAE_p01_resid_2p0_resid_2p0_lme.csv');
[mag_3p0,depth_3p0,dist_3p0,resid_c_ci95l_3p0,resid_c_3p0,resid_c_ci95u_3p0,resid_phi_ci95l_3p0,resid_phi_3p0,resid_phi_ci95u_3p0,resid_tau_ci95l_3p0,resid_tau_3p0,resid_tau_ci95u_3p0,resid_intra_3p0,resid_inter_3p0,vs30_3p0]=read_data_1csv('induced_lme_20171103/NGAE_p01_resid_3p0_resid_3p0_lme.csv');
[mag_4p0,depth_4p0,dist_4p0,resid_c_ci95l_4p0,resid_c_4p0,resid_c_ci95u_4p0,resid_phi_ci95l_4p0,resid_phi_4p0,resid_phi_ci95u_4p0,resid_tau_ci95l_4p0,resid_tau_4p0,resid_tau_ci95u_4p0,resid_intra_4p0,resid_inter_4p0,vs30_4p0]=read_data_1csv('induced_lme_20171103/NGAE_p01_resid_4p0_resid_4p0_lme.csv');
[mag_5p0,depth_5p0,dist_5p0,resid_c_ci95l_5p0,resid_c_5p0,resid_c_ci95u_5p0,resid_phi_ci95l_5p0,resid_phi_5p0,resid_phi_ci95u_5p0,resid_tau_ci95l_5p0,resid_tau_5p0,resid_tau_ci95u_5p0,resid_intra_5p0,resid_inter_5p0,vs30_5p0]=read_data_1csv('induced_lme_20171103/NGAE_p01_resid_5p0_resid_5p0_lme.csv');
[mag_10p0,depth_10p0,dist_10p0,resid_c_ci95l_10p0,resid_c_10p0,resid_c_ci95u_10p0,resid_phi_ci95l_10p0,resid_phi_10p0,resid_phi_ci95u_10p0,resid_tau_ci95l_10p0,resid_tau_10p0,resid_tau_ci95u_10p0,resid_intra_10p0,resid_inter_10p0,vs30_10p0]=read_data_1csv('induced_lme_20171103/NGAE_p01_resid_10p0_resid_10p0_lme.csv');
%
per_arr=[0.1 0.2 0.3 0.4 0.5 1.0 2.0 3.0 4.0 5.0 10.0];
resid_c_ci95l=[resid_c_ci95l_0p1 resid_c_ci95l_0p2 resid_c_ci95l_0p3 resid_c_ci95l_0p4 resid_c_ci95l_0p5 resid_c_ci95l_1p0 resid_c_ci95l_2p0 resid_c_ci95l_3p0 resid_c_ci95l_4p0 resid_c_ci95l_5p0 resid_c_ci95l_10p0];
resid_c=[resid_c_0p1 resid_c_0p2 resid_c_0p3 resid_c_0p4 resid_c_0p5 resid_c_1p0 resid_c_2p0 resid_c_3p0 resid_c_4p0 resid_c_5p0 resid_c_10p0];
resid_c_ci95u=[resid_c_ci95u_0p1 resid_c_ci95u_0p2 resid_c_ci95u_0p3 resid_c_ci95u_0p4 resid_c_ci95u_0p5 resid_c_ci95u_1p0 resid_c_ci95u_2p0 resid_c_ci95u_3p0 resid_c_ci95u_4p0 resid_c_ci95u_5p0 resid_c_ci95u_10p0];
resid_phi_ci95l=[resid_phi_ci95l_0p1 resid_phi_ci95l_0p2 resid_phi_ci95l_0p3 resid_phi_ci95l_0p4 resid_phi_ci95l_0p5 resid_phi_ci95l_1p0 resid_phi_ci95l_2p0 resid_phi_ci95l_3p0 resid_phi_ci95l_4p0 resid_phi_ci95l_5p0 resid_phi_ci95l_10p0];
resid_phi=[resid_phi_0p1 resid_phi_0p2 resid_phi_0p3 resid_phi_0p4 resid_phi_0p5 resid_phi_1p0 resid_phi_2p0 resid_phi_3p0 resid_phi_4p0 resid_phi_5p0 resid_phi_10p0];
resid_phi_ci95u=[resid_phi_ci95u_0p1 resid_phi_ci95u_0p2 resid_phi_ci95u_0p3 resid_phi_ci95u_0p4 resid_phi_ci95u_0p5 resid_phi_ci95u_1p0 resid_phi_ci95u_2p0 resid_phi_ci95u_3p0 resid_phi_ci95u_4p0 resid_phi_ci95u_5p0 resid_phi_ci95u_10p0];
resid_tau_ci95l=[resid_tau_ci95l_0p1 resid_tau_ci95l_0p2 resid_tau_ci95l_0p3 resid_tau_ci95l_0p4 resid_tau_ci95l_0p5 resid_tau_ci95l_1p0 resid_tau_ci95l_2p0 resid_tau_ci95l_3p0 resid_tau_ci95l_4p0 resid_tau_ci95l_5p0 resid_tau_ci95l_10p0];
resid_tau=[resid_tau_0p1 resid_tau_0p2 resid_tau_0p3 resid_tau_0p4 resid_tau_0p5 resid_tau_1p0 resid_tau_2p0 resid_tau_3p0 resid_tau_4p0 resid_tau_5p0 resid_tau_10p0];
resid_tau_ci95u=[resid_tau_ci95u_0p1 resid_tau_ci95u_0p2 resid_tau_ci95u_0p3 resid_tau_ci95u_0p4 resid_tau_ci95u_0p5 resid_tau_ci95u_1p0 resid_tau_ci95u_2p0 resid_tau_ci95u_3p0 resid_tau_ci95u_4p0 resid_tau_ci95u_5p0 resid_tau_ci95u_10p0];
%
dist.t0p1=dist_0p1;
dist.t0p2=dist_0p2;
dist.t0p3=dist_0p3;
dist.t0p4=dist_0p4;
dist.t0p5=dist_0p5;
dist.t1p0=dist_1p0;
dist.t2p0=dist_2p0;
dist.t3p0=dist_3p0;
dist.t4p0=dist_4p0;
dist.t5p0=dist_5p0;
dist.t10p0=dist_10p0;
%
mag.t0p1=mag_0p1;
mag.t0p2=mag_0p2;
mag.t0p3=mag_0p3;
mag.t0p4=mag_0p4;
mag.t0p5=mag_0p5;
mag.t1p0=mag_1p0;
mag.t2p0=mag_2p0;
mag.t3p0=mag_3p0;
mag.t4p0=mag_4p0;
mag.t5p0=mag_5p0;
mag.t10p0=mag_10p0;
%
vs30.t0p1=vs30_0p1;
vs30.t0p2=vs30_0p2;
vs30.t0p3=vs30_0p3;
vs30.t0p4=vs30_0p4;
vs30.t0p5=vs30_0p5;
vs30.t1p0=vs30_1p0;
vs30.t2p0=vs30_2p0;
vs30.t3p0=vs30_3p0;
vs30.t4p0=vs30_4p0;
vs30.t5p0=vs30_5p0;
vs30.t10p0=vs30_10p0;
%
depth.t0p1=depth_0p1;
depth.t0p2=depth_0p2;
depth.t0p3=depth_0p3;
depth.t0p4=depth_0p4;
depth.t0p5=depth_0p5;
depth.t1p0=depth_1p0;
depth.t2p0=depth_2p0;
depth.t3p0=depth_3p0;
depth.t4p0=depth_4p0;
depth.t5p0=depth_5p0;
depth.t10p0=depth_10p0;
%
resid_intra.t0p1=resid_intra_0p1;
resid_intra.t0p2=resid_intra_0p2;
resid_intra.t0p3=resid_intra_0p3;
resid_intra.t0p4=resid_intra_0p4;
resid_intra.t0p5=resid_intra_0p5;
resid_intra.t1p0=resid_intra_1p0;
resid_intra.t2p0=resid_intra_2p0;
resid_intra.t3p0=resid_intra_3p0;
resid_intra.t4p0=resid_intra_4p0;
resid_intra.t5p0=resid_intra_5p0;
resid_intra.t10p0=resid_intra_10p0;
%
resid_inter.t0p1=resid_inter_0p1;
resid_inter.t0p2=resid_inter_0p2;
resid_inter.t0p3=resid_inter_0p3;
resid_inter.t0p4=resid_inter_0p4;
resid_inter.t0p5=resid_inter_0p5;
resid_inter.t1p0=resid_inter_1p0;
resid_inter.t2p0=resid_inter_2p0;
resid_inter.t3p0=resid_inter_3p0;
resid_inter.t4p0=resid_inter_4p0;
resid_inter.t5p0=resid_inter_5p0;
resid_inter.t10p0=resid_inter_10p0;


end
%-----------------------------------------------------


%-----------------------------------------------------
function [dist,mag,depth,per_arr,resid_c_ci95l,resid_c,resid_c_ci95u,resid_phi_ci95l,resid_phi,resid_phi_ci95u,resid_tau_ci95l,resid_tau,resid_tau_ci95u,resid_intra,resid_inter,vs30]=read_data_NGAW2()

%
[mag_0p1,depth_0p1,dist_0p1,resid_c_ci95l_0p1,resid_c_0p1,resid_c_ci95u_0p1,resid_phi_ci95l_0p1,resid_phi_0p1,resid_phi_ci95u_0p1,resid_tau_ci95l_0p1,resid_tau_0p1,resid_tau_ci95u_0p1,resid_intra_0p1,resid_inter_0p1,vs30_0p1]=read_data_1csv('induced_lme_20171103/NGAW2_resid_0p1_resid_0p1_lme.csv');
[mag_0p2,depth_0p2,dist_0p2,resid_c_ci95l_0p2,resid_c_0p2,resid_c_ci95u_0p2,resid_phi_ci95l_0p2,resid_phi_0p2,resid_phi_ci95u_0p2,resid_tau_ci95l_0p2,resid_tau_0p2,resid_tau_ci95u_0p2,resid_intra_0p2,resid_inter_0p2,vs30_0p2]=read_data_1csv('induced_lme_20171103/NGAW2_resid_0p2_resid_0p2_lme.csv');
[mag_0p3,depth_0p3,dist_0p3,resid_c_ci95l_0p3,resid_c_0p3,resid_c_ci95u_0p3,resid_phi_ci95l_0p3,resid_phi_0p3,resid_phi_ci95u_0p3,resid_tau_ci95l_0p3,resid_tau_0p3,resid_tau_ci95u_0p3,resid_intra_0p3,resid_inter_0p3,vs30_0p3]=read_data_1csv('induced_lme_20171103/NGAW2_resid_0p3_resid_0p3_lme.csv');
[mag_0p4,depth_0p4,dist_0p4,resid_c_ci95l_0p4,resid_c_0p4,resid_c_ci95u_0p4,resid_phi_ci95l_0p4,resid_phi_0p4,resid_phi_ci95u_0p4,resid_tau_ci95l_0p4,resid_tau_0p4,resid_tau_ci95u_0p4,resid_intra_0p4,resid_inter_0p4,vs30_0p4]=read_data_1csv('induced_lme_20171103/NGAW2_resid_0p4_resid_0p4_lme.csv');
[mag_0p5,depth_0p5,dist_0p5,resid_c_ci95l_0p5,resid_c_0p5,resid_c_ci95u_0p5,resid_phi_ci95l_0p5,resid_phi_0p5,resid_phi_ci95u_0p5,resid_tau_ci95l_0p5,resid_tau_0p5,resid_tau_ci95u_0p5,resid_intra_0p5,resid_inter_0p5,vs30_0p5]=read_data_1csv('induced_lme_20171103/NGAW2_resid_0p5_resid_0p5_lme.csv');
[mag_1p0,depth_1p0,dist_1p0,resid_c_ci95l_1p0,resid_c_1p0,resid_c_ci95u_1p0,resid_phi_ci95l_1p0,resid_phi_1p0,resid_phi_ci95u_1p0,resid_tau_ci95l_1p0,resid_tau_1p0,resid_tau_ci95u_1p0,resid_intra_1p0,resid_inter_1p0,vs30_1p0]=read_data_1csv('induced_lme_20171103/NGAW2_resid_1p0_resid_1p0_lme.csv');
[mag_2p0,depth_2p0,dist_2p0,resid_c_ci95l_2p0,resid_c_2p0,resid_c_ci95u_2p0,resid_phi_ci95l_2p0,resid_phi_2p0,resid_phi_ci95u_2p0,resid_tau_ci95l_2p0,resid_tau_2p0,resid_tau_ci95u_2p0,resid_intra_2p0,resid_inter_2p0,vs30_2p0]=read_data_1csv('induced_lme_20171103/NGAW2_resid_2p0_resid_2p0_lme.csv');
[mag_3p0,depth_3p0,dist_3p0,resid_c_ci95l_3p0,resid_c_3p0,resid_c_ci95u_3p0,resid_phi_ci95l_3p0,resid_phi_3p0,resid_phi_ci95u_3p0,resid_tau_ci95l_3p0,resid_tau_3p0,resid_tau_ci95u_3p0,resid_intra_3p0,resid_inter_3p0,vs30_3p0]=read_data_1csv('induced_lme_20171103/NGAW2_resid_3p0_resid_3p0_lme.csv');
[mag_4p0,depth_4p0,dist_4p0,resid_c_ci95l_4p0,resid_c_4p0,resid_c_ci95u_4p0,resid_phi_ci95l_4p0,resid_phi_4p0,resid_phi_ci95u_4p0,resid_tau_ci95l_4p0,resid_tau_4p0,resid_tau_ci95u_4p0,resid_intra_4p0,resid_inter_4p0,vs30_4p0]=read_data_1csv('induced_lme_20171103/NGAW2_resid_4p0_resid_4p0_lme.csv');
[mag_5p0,depth_5p0,dist_5p0,resid_c_ci95l_5p0,resid_c_5p0,resid_c_ci95u_5p0,resid_phi_ci95l_5p0,resid_phi_5p0,resid_phi_ci95u_5p0,resid_tau_ci95l_5p0,resid_tau_5p0,resid_tau_ci95u_5p0,resid_intra_5p0,resid_inter_5p0,vs30_5p0]=read_data_1csv('induced_lme_20171103/NGAW2_resid_5p0_resid_5p0_lme.csv');
[mag_10p0,depth_10p0,dist_10p0,resid_c_ci95l_10p0,resid_c_10p0,resid_c_ci95u_10p0,resid_phi_ci95l_10p0,resid_phi_10p0,resid_phi_ci95u_10p0,resid_tau_ci95l_10p0,resid_tau_10p0,resid_tau_ci95u_10p0,resid_intra_10p0,resid_inter_10p0,vs30_10p0]=read_data_1csv('induced_lme_20171103/NGAW2_resid_10p0_resid_10p0_lme.csv');
%
per_arr=[0.1 0.2 0.3 0.4 0.5 1.0 2.0 3.0 4.0 5.0 10.0];
resid_c_ci95l=[resid_c_ci95l_0p1 resid_c_ci95l_0p2 resid_c_ci95l_0p3 resid_c_ci95l_0p4 resid_c_ci95l_0p5 resid_c_ci95l_1p0 resid_c_ci95l_2p0 resid_c_ci95l_3p0 resid_c_ci95l_4p0 resid_c_ci95l_5p0 resid_c_ci95l_10p0];
resid_c=[resid_c_0p1 resid_c_0p2 resid_c_0p3 resid_c_0p4 resid_c_0p5 resid_c_1p0 resid_c_2p0 resid_c_3p0 resid_c_4p0 resid_c_5p0 resid_c_10p0];
resid_c_ci95u=[resid_c_ci95u_0p1 resid_c_ci95u_0p2 resid_c_ci95u_0p3 resid_c_ci95u_0p4 resid_c_ci95u_0p5 resid_c_ci95u_1p0 resid_c_ci95u_2p0 resid_c_ci95u_3p0 resid_c_ci95u_4p0 resid_c_ci95u_5p0 resid_c_ci95u_10p0];
resid_phi_ci95l=[resid_phi_ci95l_0p1 resid_phi_ci95l_0p2 resid_phi_ci95l_0p3 resid_phi_ci95l_0p4 resid_phi_ci95l_0p5 resid_phi_ci95l_1p0 resid_phi_ci95l_2p0 resid_phi_ci95l_3p0 resid_phi_ci95l_4p0 resid_phi_ci95l_5p0 resid_phi_ci95l_10p0];
resid_phi=[resid_phi_0p1 resid_phi_0p2 resid_phi_0p3 resid_phi_0p4 resid_phi_0p5 resid_phi_1p0 resid_phi_2p0 resid_phi_3p0 resid_phi_4p0 resid_phi_5p0 resid_phi_10p0];
resid_phi_ci95u=[resid_phi_ci95u_0p1 resid_phi_ci95u_0p2 resid_phi_ci95u_0p3 resid_phi_ci95u_0p4 resid_phi_ci95u_0p5 resid_phi_ci95u_1p0 resid_phi_ci95u_2p0 resid_phi_ci95u_3p0 resid_phi_ci95u_4p0 resid_phi_ci95u_5p0 resid_phi_ci95u_10p0];
resid_tau_ci95l=[resid_tau_ci95l_0p1 resid_tau_ci95l_0p2 resid_tau_ci95l_0p3 resid_tau_ci95l_0p4 resid_tau_ci95l_0p5 resid_tau_ci95l_1p0 resid_tau_ci95l_2p0 resid_tau_ci95l_3p0 resid_tau_ci95l_4p0 resid_tau_ci95l_5p0 resid_tau_ci95l_10p0];
resid_tau=[resid_tau_0p1 resid_tau_0p2 resid_tau_0p3 resid_tau_0p4 resid_tau_0p5 resid_tau_1p0 resid_tau_2p0 resid_tau_3p0 resid_tau_4p0 resid_tau_5p0 resid_tau_10p0];
resid_tau_ci95u=[resid_tau_ci95u_0p1 resid_tau_ci95u_0p2 resid_tau_ci95u_0p3 resid_tau_ci95u_0p4 resid_tau_ci95u_0p5 resid_tau_ci95u_1p0 resid_tau_ci95u_2p0 resid_tau_ci95u_3p0 resid_tau_ci95u_4p0 resid_tau_ci95u_5p0 resid_tau_ci95u_10p0];
%
dist.t0p1=dist_0p1;
dist.t0p2=dist_0p2;
dist.t0p3=dist_0p3;
dist.t0p4=dist_0p4;
dist.t0p5=dist_0p5;
dist.t1p0=dist_1p0;
dist.t2p0=dist_2p0;
dist.t3p0=dist_3p0;
dist.t4p0=dist_4p0;
dist.t5p0=dist_5p0;
dist.t10p0=dist_10p0;
%
mag.t0p1=mag_0p1;
mag.t0p2=mag_0p2;
mag.t0p3=mag_0p3;
mag.t0p4=mag_0p4;
mag.t0p5=mag_0p5;
mag.t1p0=mag_1p0;
mag.t2p0=mag_2p0;
mag.t3p0=mag_3p0;
mag.t4p0=mag_4p0;
mag.t5p0=mag_5p0;
mag.t10p0=mag_10p0;
%
vs30.t0p1=vs30_0p1;
vs30.t0p2=vs30_0p2;
vs30.t0p3=vs30_0p3;
vs30.t0p4=vs30_0p4;
vs30.t0p5=vs30_0p5;
vs30.t1p0=vs30_1p0;
vs30.t2p0=vs30_2p0;
vs30.t3p0=vs30_3p0;
vs30.t4p0=vs30_4p0;
vs30.t5p0=vs30_5p0;
vs30.t10p0=vs30_10p0;
%
depth.t0p1=depth_0p1;
depth.t0p2=depth_0p2;
depth.t0p3=depth_0p3;
depth.t0p4=depth_0p4;
depth.t0p5=depth_0p5;
depth.t1p0=depth_1p0;
depth.t2p0=depth_2p0;
depth.t3p0=depth_3p0;
depth.t4p0=depth_4p0;
depth.t5p0=depth_5p0;
depth.t10p0=depth_10p0;
%
resid_intra.t0p1=resid_intra_0p1;
resid_intra.t0p2=resid_intra_0p2;
resid_intra.t0p3=resid_intra_0p3;
resid_intra.t0p4=resid_intra_0p4;
resid_intra.t0p5=resid_intra_0p5;
resid_intra.t1p0=resid_intra_1p0;
resid_intra.t2p0=resid_intra_2p0;
resid_intra.t3p0=resid_intra_3p0;
resid_intra.t4p0=resid_intra_4p0;
resid_intra.t5p0=resid_intra_5p0;
resid_intra.t10p0=resid_intra_10p0;
%
resid_inter.t0p1=resid_inter_0p1;
resid_inter.t0p2=resid_inter_0p2;
resid_inter.t0p3=resid_inter_0p3;
resid_inter.t0p4=resid_inter_0p4;
resid_inter.t0p5=resid_inter_0p5;
resid_inter.t1p0=resid_inter_1p0;
resid_inter.t2p0=resid_inter_2p0;
resid_inter.t3p0=resid_inter_3p0;
resid_inter.t4p0=resid_inter_4p0;
resid_inter.t5p0=resid_inter_5p0;
resid_inter.t10p0=resid_inter_10p0;

end
%-----------------------------------------------------

%-----------------------------------------------------
function[mag,depth,dist,resid_c_ci95l,resid_c,resid_c_ci95u,resid_phi_ci95l,resid_phi,resid_phi_ci95u,resid_tau_ci95l,resid_tau,resid_tau_ci95u,resid_intra,resid_inter,vs30]=read_data_1csvCEUS(csvfile);

%
ff=csvread(csvfile,1,0);
gmN=ff(:,1);
mag=ff(:,2);
depth=ff(:,3);
dist=ff(:,4);
originTime=ff(:,5);
resid0p1=ff(:,6);
resid0p2=ff(:,7);
resid0p3=ff(:,8);
resid0p5=ff(:,9);
resid1p0=ff(:,10);
resid2p0=ff(:,11);
vs30=ff(:,12);
evid=ff(:,13);
resid_c_ci95l=ff(:,14);
resid_c=ff(:,15);
resid_c_ci95u=ff(:,16);
resid_phi_ci95l=ff(:,17);
resid_phi=ff(:,18);
resid_phi_ci95u=ff(:,19);
resid_tau_ci95l=ff(:,20);
resid_tau=ff(:,21);
resid_tau_ci95u=ff(:,22);
resid_intra=ff(:,23);
resid_inter=ff(:,24);
% means
resid_c_ci95l=mean(resid_c_ci95l);
resid_c=mean(resid_c);
resid_c_ci95u=mean(resid_c_ci95u);
resid_phi_ci95l=mean(resid_phi_ci95l);
resid_phi=mean(resid_phi);
resid_phi_ci95u=mean(resid_phi_ci95u);
resid_tau_ci95l=mean(resid_tau_ci95l);
resid_tau=mean(resid_tau);
resid_tau_ci95u=mean(resid_tau_ci95u);
%

end
%-----------------------------------------------------

%-----------------------------------------------------
function[mag,depth,dist,resid_c_ci95l,resid_c,resid_c_ci95u,resid_phi_ci95l,resid_phi,resid_phi_ci95u,resid_tau_ci95l,resid_tau,resid_tau_ci95u,resid_intra,resid_inter,vs30]=read_data_1csv(csvfile);

%
ff=csvread(csvfile,1,0);
gmN=ff(:,1);
mag=ff(:,2);
depth=ff(:,3);
dist=ff(:,4);
originTime=ff(:,5);
resid0p1=ff(:,6);
resid0p2=ff(:,7);
resid0p3=ff(:,8);
resid0p4=ff(:,9);
resid0p5=ff(:,10);
resid1p0=ff(:,11);
resid2p0=ff(:,12);
resid3p0=ff(:,13);
resid4p0=ff(:,14);
resid5p0=ff(:,15);
resid10p0=ff(:,16);
vs30=ff(:,17);
evid=ff(:,18);
resid_c_ci95l=ff(:,19);
resid_c=ff(:,20);
resid_c_ci95u=ff(:,21);
resid_phi_ci95l=ff(:,22);
resid_phi=ff(:,23);
resid_phi_ci95u=ff(:,24);
resid_tau_ci95l=ff(:,25);
resid_tau=ff(:,26);
resid_tau_ci95u=ff(:,27);
resid_intra=ff(:,28);
resid_inter=ff(:,29);
% means
resid_c_ci95l=mean(resid_c_ci95l);
resid_c=mean(resid_c);
resid_c_ci95u=mean(resid_c_ci95u);
resid_phi_ci95l=mean(resid_phi_ci95l);
resid_phi=mean(resid_phi);
resid_phi_ci95u=mean(resid_phi_ci95u);
resid_tau_ci95l=mean(resid_tau_ci95l);
resid_tau=mean(resid_tau);
resid_tau_ci95u=mean(resid_tau_ci95u);
%

end
%-----------------------------------------------------


%-----------------------------------------------------
function [dist_bins, m_1p5s, m_2s, m_3s, m_4s, m_5s, m_7p5s, m_10s]=plot_withinEvent_rupPar1(rupPar,holdV)

%
per_arr=[1.5 2 3 4 5 7.5 10];
fout=sprintf('arr_intraVar_%s.tmp',rupPar);

% 
ff=load(fout);
rrup=ff(:,1);
intra_1p5s=ff(:,2);
intra_2s=ff(:,3);
intra_3s=ff(:,4);
intra_4s=ff(:,5);
intra_5s=ff(:,6);
intra_7p5s=ff(:,7);
intra_10s=ff(:,8);

% 
[dist_bins_1p5s,cnt_distBins_1p5s, mean_distBins_1p5s,var_distBins_1p5s,mean_distBins_sm_1p5s,var_distBins_sm_1p5s]=bin_logData_by_dist(rrup,intra_1p5s);
[dist_bins_2s,cnt_distBins_2s, mean_distBins_2s,var_distBins_2s,mean_distBins_sm_2s,var_distBins_sm_2s]=bin_logData_by_dist(rrup,intra_2s);
[dist_bins_3s,cnt_distBins_3s, mean_distBins_3s,var_distBins_3s,mean_distBins_sm_3s,var_distBins_sm_3s]=bin_logData_by_dist(rrup,intra_3s);
[dist_bins_4s,cnt_distBins_4s, mean_distBins_4s,var_distBins_4s,mean_distBins_sm_4s,var_distBins_sm_4s]=bin_logData_by_dist(rrup,intra_4s);
[dist_bins_5s,cnt_distBins_5s, mean_distBins_5s,var_distBins_5s,mean_distBins_sm_5s,var_distBins_sm_5s]=bin_logData_by_dist(rrup,intra_5s);
[dist_bins_7p5s,cnt_distBins_7p5s, mean_distBins_7p5s,var_distBins_7p5s,mean_distBins_sm_7p5s,var_distBins_sm_7p5s]=bin_logData_by_dist(rrup,intra_7p5s);
[dist_bins_10s,cnt_distBins_10s, mean_distBins_10s,var_distBins_10s,mean_distBins_sm_10s,var_distBins_sm_10s]=bin_logData_by_dist(rrup,intra_10s);
%
dist_bins=dist_bins_1p5s;
use_sm_val=0
if use_sm_val
  m_1p5s=mean_distBins_sm_1p5s;
  m_2s=mean_distBins_sm_2s;
  m_3s=mean_distBins_sm_3s;
  m_4s=mean_distBins_sm_4s;
  m_5s=mean_distBins_sm_5s;
  m_7p5s=mean_distBins_sm_7p5s;
  m_10s=mean_distBins_sm_10s;
else
  m_1p5s=mean_distBins_1p5s;
  m_2s=mean_distBins_2s;
  m_3s=mean_distBins_3s;
  m_4s=mean_distBins_4s;
  m_5s=mean_distBins_5s;
  m_7p5s=mean_distBins_7p5s;
  m_10s=mean_distBins_10s;
end
subplot(2,3,1)
plot(dist_bins,m_1p5s)
set(gca,'XScale','log'),
if holdV
  hold on,
  titlepr=sprintf('T=1.5s, %s',rupPar);
%  title('T=1.5s')
  title(titlepr)
end
subplot(2,3,2)
plot(dist_bins,m_2s)
set(gca,'XScale','log'),
if holdV
  hold on,
  title('T=2s')
end
subplot(2,3,3)
plot(dist_bins,m_3s)
set(gca,'XScale','log'),
if holdV
  hold on,
  title('T=3s')
end
subplot(2,3,4)
plot(dist_bins,m_5s)
set(gca,'XScale','log'),
if holdV
  hold on,
  title('T=5s')
end
subplot(2,3,5)
plot(dist_bins,m_7p5s)
set(gca,'XScale','log'),
if holdV
  hold on,
  title('T=7.5s')
end
subplot(2,3,6)
plot(dist_bins,m_10s)
set(gca,'XScale','log'),
if holdV
  hold on,
  title('T=10s')
end



end
%-----------------------------------------------------

%-----------------------------------------------------
function []=plot_interEvent_rupPar2()

%
figure
[per_arr_sv1,mean_arr_sv1]=plot_interEvent_rupPar1('sv1','b','s'), hold on
[per_arr_sv2,mean_arr_sv2]=plot_interEvent_rupPar1('sv2','b','d'),
[per_arr_rv60,mean_arr_rv60]=plot_interEvent_rupPar1('rv60','r','s'), 
[per_arr_rv80,mean_arr_rv80]=plot_interEvent_rupPar1('rv80','r','d'),
[per_arr_hypoS,mean_arr_hypoS]=plot_interEvent_rupPar1('hypoS','c','s'), 
[per_arr_hypoN,mean_arr_hypoN]=plot_interEvent_rupPar1('hypoN','c','d'),
[per_arr_hypoC,mean_arr_hypoC]=plot_interEvent_rupPar1('hypoC','c','x'),
[per_arr_ax13,mean_arr_ax13]=plot_interEvent_rupPar1('ax13','m','s'), 
[per_arr_ax3,mean_arr_ax3]=plot_interEvent_rupPar1('ax3','m','d'),
[per_arr_seed103,mean_arr_seed103]=plot_interEvent_rupPar1('seed103','g','s'), 
[per_arr_seed103r,mean_arr_seed103r]=plot_interEvent_rupPar1('seed103r','g','d'), 
[per_arr_seed56,mean_arr_seed56]=plot_interEvent_rupPar1('seed56','g','o'), 
[per_arr_seed83,mean_arr_seed83]=plot_interEvent_rupPar1('seed83','g','x'), 
close
%
per_arr=per_arr_sv1;
mean_arr_sv=[mean_arr_sv1; mean_arr_sv2];
mean_arr_rv=[mean_arr_rv60; mean_arr_rv80];
mean_arr_ax=[mean_arr_ax13; mean_arr_ax3];
mean_arr_hypo=[mean_arr_hypoS; mean_arr_hypoN; mean_arr_hypoC];
mean_arr_seed=[mean_arr_seed103; mean_arr_seed103r; mean_arr_seed56; mean_arr_seed83];
%
figure
plot_interEvent_mean2(per_arr,mean_arr_sv,'b'), hold on
plot_interEvent_mean2(per_arr,mean_arr_rv,'r'), 
plot_interEvent_mean2(per_arr,mean_arr_hypo,'c'), 
plot_interEvent_mean2(per_arr,mean_arr_ax,'m'), 
plot_interEvent_mean2(per_arr,mean_arr_seed,'g'), 
xlabel('T (s)')
title('Effect of rupture parameter on between-event residuals')
set(gca,'XScale','log')
legend('Slip velocity','Rupture speed','Hypocenter','Correlation length','Random field')

end
%-----------------------------------------------------



%-----------------------------------------------------
function []=plot_interEvent_rupPar()

%
figure
subplot(1,2,1)
[per_arr_sv1,mean_arr_sv1]=plot_interEvent_rupPar1('sv1','b','s'), hold on
[per_arr_sv2,mean_arr_sv2]=plot_interEvent_rupPar1('sv2','b','d'),
[per_arr_rv60,mean_arr_rv60]=plot_interEvent_rupPar1('rv60','r','s'), 
[per_arr_rv80,mean_arr_rv80]=plot_interEvent_rupPar1('rv80','r','d'),
[per_arr_hypoS,mean_arr_hypoS]=plot_interEvent_rupPar1('hypoS','c','s'), 
[per_arr_hypoN,mean_arr_hypoN]=plot_interEvent_rupPar1('hypoN','c','d'),
[per_arr_hypoC,mean_arr_hypoC]=plot_interEvent_rupPar1('hypoC','c','x'),
[per_arr_ax13,mean_arr_ax13]=plot_interEvent_rupPar1('ax13','m','s'), 
[per_arr_ax3,mean_arr_ax3]=plot_interEvent_rupPar1('ax3','m','d'),
[per_arr_seed103,mean_arr_seed103]=plot_interEvent_rupPar1('seed103','g','s'), 
[per_arr_seed103r,mean_arr_seed103r]=plot_interEvent_rupPar1('seed103r','g','d'), 
[per_arr_seed56,mean_arr_seed56]=plot_interEvent_rupPar1('seed56','g','o'), 
[per_arr_seed83,mean_arr_seed83]=plot_interEvent_rupPar1('seed83','g','x'), 
set(gca,'XScale','log')
%
per_arr=per_arr_sv1;
mean_arr_sv=[mean_arr_sv1; mean_arr_sv2];
mean_arr_rv=[mean_arr_rv60; mean_arr_rv80];
mean_arr_ax=[mean_arr_ax13; mean_arr_ax3];
mean_arr_hypo=[mean_arr_hypoS; mean_arr_hypoN; mean_arr_hypoC];
mean_arr_seed=[mean_arr_seed103; mean_arr_seed103r; mean_arr_seed56; mean_arr_seed83];
subplot(1,2,2),
plot_interEvent_mean(per_arr,mean_arr_sv,'b'), hold on
plot_interEvent_mean(per_arr,mean_arr_rv,'r'), 
plot_interEvent_mean(per_arr,mean_arr_hypo,'c'), 
plot_interEvent_mean(per_arr,mean_arr_ax,'m'), 
plot_interEvent_mean(per_arr,mean_arr_seed,'g'), 
set(gca,'XScale','log')

end
%-----------------------------------------------------

%-----------------------------------------------------
function []=plot_interEvent_mean2(perV,meanV_arr,plotColor)

%
plotComb=sprintf('%s',plotColor);

% 
maxV=max(meanV_arr,[],1);
minV=min(meanV_arr,[],1);
diffV=maxV-minV;

%
plot(perV,diffV,plotComb,'LineWidth',2),

end
%-----------------------------------------------------


%-----------------------------------------------------
function []=plot_interEvent_mean(perV,meanV_arr,plotColor)

%
plotComb=sprintf('%ss',plotColor);

% 
maxV=max(meanV_arr,[],1);
minV=min(meanV_arr,[],1);
diffV=maxV-minV;

%
plot(perV,diffV,plotComb),

end
%-----------------------------------------------------


%-----------------------------------------------------
function [per_arr,mean_arr]=plot_interEvent_rupPar1(rupPar,plotColor,plotSymbol)

%
fout=sprintf('arr_interVar_%s.tmp',rupPar);
per_arr=[1.5 2 3 4 5 7.5 10];

%
plotComb=sprintf('%s%s',plotColor,plotSymbol);

%
ff=load(fout);
inter_1p5s=ff(:,1);
inter_2s=ff(:,2);
inter_3s=ff(:,3);
inter_4s=ff(:,4);
inter_5s=ff(:,5);
inter_7p5s=ff(:,6);
inter_10s=ff(:,7);
mean_arr=[mean(inter_1p5s) mean(inter_2s) mean(inter_3s) mean(inter_4s) mean(inter_5s) mean(inter_7p5s) mean(inter_10s)];

%
plot(per_arr,mean_arr,plotComb,'MarkerFaceColor',plotColor),


end
%-----------------------------------------------------

%-----------------------------------------------------
function []=plot1intra_meanStd(rrup,intraV,per)

% distance binning for means/std 
[dist_bins,cnt_distBins, mean_distBins,var_distBins,mean_distBins_sm,var_distBins_sm]=bin_logData_by_dist(rrup,intraV);

%
subplot(1,2,1),
plot(dist_bins,mean_distBins_sm), 
if per<2
  hold on
end
if per>9
  set(gca,'XScale','log')
  title('mean-intra-event variability')
end
%
subplot(1,2,2),
plot(dist_bins,sqrt(var_distBins_sm)), 
if per<2
  hold on
end
if per>9
  legend('1.5','2','3','5','7.5','10')
  set(gca,'XScale','log')
  title('SD-intra-event variability')
end


end
%-----------------------------------------------------



%-----------------------------------------------------
function []=plot1intra(rrup,intraV,per)

% distance binning for means/std 
[dist_bins,cnt_distBins, mean_distBins,var_distBins,mean_distBins_sm,var_distBins_sm]=bin_logData_by_dist(rrup,intraV);

% plotting
plot(rrup,intraV,'bs'), hold on
% means/std
plot(dist_bins,mean_distBins,'k-'), 
plot(dist_bins,mean_distBins+sqrt(var_distBins),'k--'), 
plot(dist_bins,mean_distBins-sqrt(var_distBins),'k--'), 
plot(dist_bins,mean_distBins_sm,'r-'), 
plot(dist_bins,mean_distBins_sm+sqrt(var_distBins_sm),'r--'), 
plot(dist_bins,mean_distBins_sm-sqrt(var_distBins_sm),'r--'), 

set(gca,'XScale','log'),



end
%-----------------------------------------------------


%-----------------------------------------------------
function []=plot1inter(mag,resid)

  [mag_bins_0p1,resid_bins_0p1]=bin_mag_resid(mag.t0p1,resid.t0p1);
  [mag_bins_0p2,resid_bins_0p2]=bin_mag_resid(mag.t0p2,resid.t0p2);
  [mag_bins_0p3,resid_bins_0p3]=bin_mag_resid(mag.t0p3,resid.t0p3);
  [mag_bins_0p5,resid_bins_0p5]=bin_mag_resid(mag.t0p5,resid.t0p5);
  [mag_bins_1p0,resid_bins_1p0]=bin_mag_resid(mag.t1p0,resid.t1p0);
  [mag_bins_2p0,resid_bins_2p0]=bin_mag_resid(mag.t2p0,resid.t2p0);

%
%  figure
  plot(mag_bins_0p1,resid_bins_0p1),hold on
  plot(mag_bins_0p2,resid_bins_0p2),
  plot(mag_bins_0p3,resid_bins_0p3),
  plot(mag_bins_0p5,resid_bins_0p5),
  plot(mag_bins_1p0,resid_bins_1p0),
  plot(mag_bins_2p0,resid_bins_2p0),
  legend('0.1','0.2','0.3','0.5','1.0','2.0')

end
%-----------------------------------------------------


%-----------------------------------------------------
function [mag_bins,resid_bins]=bin_mag_resid(mag,resid)

%
  magInc=0.25;
  mag_bins_orig=4:magInc:6;
  mag_bins=mag_bins_orig-0.05;
  for ii=1:length(mag_bins)-1
    magmin=mag_bins(ii);
    magmax=mag_bins(ii+1);
    eleF=find(mag>=magmin&mag<magmax);
    if (eleF>0) 
      resid_bins(ii)=mean(resid(eleF));
      cntBins(ii)=length(resid(eleF));
    else
      resid_bins(ii)=-999;
    end
  end

%
  eleF2=find(resid_bins>-999);
  mag_bins=mag_bins(eleF2);
  mag_bins=mag_bins+0.05;
  cntBins=cntBins(eleF2);
  resid_bins=resid_bins(eleF2);
%save aaa

%  figure
%  plot(mag_bins,cntBins,'bs-')

end
%-----------------------------------------------------


%-----------------------------------------------------
function [c_array]=get_bias_array(ff,ele1)

%
c_array=[mean(ff(:,ele1)) mean(ff(:,ele1+11)) mean(ff(:,ele1+11*2)) mean(ff(:,ele1+11*3)) mean(ff(:,ele1+11*4)) mean(ff(:,ele1+11*5)) mean(ff(:,ele1+11*6))];

end
%-----------------------------------------------------
