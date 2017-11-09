function [vs30_bins,cnt_vs30Bins, mean_vs30Bins,var_vs30Bins,mean_vs30Bins_sm,var_vs30Bins_sm]=bin_data_vs30(vs30_arr_1d,data_arr_1d)

%
plotVals=0;

% vs30 bins
nbins=20
nbins=15
vs30_bins=logspace(log10(120),log10(1860),nbins);

save eee
% loop over all data, calculate mean/var/cnt for all distance bins
length(vs30_bins)
cnt=1;
for ii=1:length(vs30_bins)-1
%  ii2=ii
  cnt=1;
  tmp_arr=-999;
  vs30_binC(ii)=0.5*(vs30_bins(ii)+vs30_bins(ii+1));
  for jj=1:length(data_arr_1d)
    if (vs30_arr_1d(jj)>=vs30_bins(ii)) & (vs30_arr_1d(jj) < vs30_bins(ii+1)) 
      tmp_arr(cnt)=data_arr_1d(jj);
      cnt=cnt+1;
    end
  end
  if ( tmp_arr(1) > -999 ) 
    mean_vs30Bins(ii)=mean(tmp_arr);
    var_vs30Bins(ii)=var(tmp_arr);
    cnt_vs30Bins(ii)=length(tmp_arr);
  else
    mean_vs30Bins(ii)=0;
    var_vs30Bins(ii)=0;
    cnt_vs30Bins(ii)=0;
  end
  clear tmp_arr
%
%
%  save all_PGV_params mean_bins_ln mean_bins sigma_bins_ln var_bins_ln var_bins Z_p1 Z_m1 Z_p2 Z_m2 vs30_bins
%  save all_3s_params mean_bins_ln mean_bins sigma_bins_ln var_bins_ln var_bins Z_p1 Z_m1 Z_p2 Z_m2 vs30_bins
%  clear tmp_arr2
end
%vs30_bins=vs30_bins(1:length(vs30_bins)-1);
save eee2

% remove data points with no data counts 
%nthresh=3
%nthresh=10
nthresh=5;
vs30_bins=vs30_binC(cnt_vs30Bins>nthresh);
mean_vs30Bins=mean_vs30Bins(cnt_vs30Bins>nthresh);
var_vs30Bins=var_vs30Bins(cnt_vs30Bins>nthresh);
cnt_vs30Bins=cnt_vs30Bins(cnt_vs30Bins>nthresh);

% smoothed values
supsmuL=1;
if supsmuL
  mean_vs30Bins_sm=supsmu(vs30_bins,mean_vs30Bins);
  var_vs30Bins_sm=supsmu(vs30_bins,var_vs30Bins);
else
% fixed smoothing
%save dsupsmu
  spanv=1e-3
  mean_vs30Bins_sm=supsmu(vs30_bins,mean_vs30Bins,'span',spanv);
  var_vs30Bins_sm=supsmu(vs30_bins,var_vs30Bins,'span',spanv);
end

%
if plotVals
%  vs30_bins
  figure
  plot(vs30_bins,mean_vs30Bins,'k-'), hold on
  plot(vs30_bins,mean_vs30Bins+sqrt(var_vs30Bins),'k-')
  plot(vs30_bins,mean_vs30Bins-sqrt(var_vs30Bins),'k-')
  plot(vs30_bins,mean_vs30Bins_sm,'r-'), hold on
  plot(vs30_bins,mean_vs30Bins_sm+sqrt(var_vs30Bins_sm),'r-')
  plot(vs30_bins,mean_vs30Bins_sm-sqrt(var_vs30Bins_sm),'r-')
  set(gca,'XScale','log'),
end


end
%----------------------------------------------------------------------
