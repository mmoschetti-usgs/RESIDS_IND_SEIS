clear all
close all

load aaa3

[vs30_bins_0p1_ngaw2,cnt_vs30Bins_0p1_ngaw2, mean_vs30Bins_0p1_ngaw2,var_vs30Bins_0p1_ngaw2,mean_vs30Bins_sm_0p1_ngaw2,var_vs30Bins_sm_0p1_ngaw2]=bin_data_vs30(vs30_ngaw2.t0p1,resid_intra_ngaw2.t0p1);

figure
plot(vs30_ngaw2.t0p1,resid_intra_ngaw2.t0p1,'bs'), hold on
plot(vs30_bins_0p1_ngaw2,mean_vs30Bins_0p1_ngaw2,'k-'),
plot(vs30_bins_0p1_ngaw2,mean_vs30Bins_sm_0p1_ngaw2,'k--'),
%set(gca,'XScale','log')
%set(gca,'XLim',[3 300]),
title('intra-event, NGA-W2 0.1 s')
xlabel('Vs30 (m/s)')
ylabel('mean(intra)')

