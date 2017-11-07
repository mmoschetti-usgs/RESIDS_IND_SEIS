function [dist_bins,cnt_distBins, mean_distBins,var_distBins,mean_distBins_sm,var_distBins_sm]=bin_logData_by_dist(dist_arr_1d,data_arr_1d)

%
plotVals=0;

% distance bins
%dist_incr=.75
%dist_bins=0:dist_incr:83-dist_incr;
%dist_bins=logspace(0.477,2.5,15);
dist_bins=logspace(0.699,2.5,10);
dist_bins=logspace(0.699,2.5,11);
%dist_bins=[0.01 dist_bins];

%dist_bins

save eee
% loop over all data, calculate mean/var/cnt for all distance bins
length(dist_bins)
for ii=1:length(dist_bins)-1
%  ii2=ii
  cnt=1;
  tmp_arr=-999;
  for jj=1:length(data_arr_1d)
    if (dist_arr_1d(jj)>=dist_bins(ii)) & (dist_arr_1d(jj) < dist_bins(ii+1)) 
      tmp_arr(cnt)=data_arr_1d(jj);
      cnt=cnt+1;
    end
  end
  if ( tmp_arr(1) > -999 ) 
    mean_distBins(ii)=mean(tmp_arr);
    var_distBins(ii)=var(tmp_arr);
    cnt_distBins(ii)=length(tmp_arr);
  else
    mean_distBins(ii)=0;
    var_distBins(ii)=0;
    cnt_distBins(ii)=0;
  end
  clear tmp_arr
%
%
%  save all_PGV_params mean_bins_ln mean_bins sigma_bins_ln var_bins_ln var_bins Z_p1 Z_m1 Z_p2 Z_m2 dist_bins
%  save all_3s_params mean_bins_ln mean_bins sigma_bins_ln var_bins_ln var_bins Z_p1 Z_m1 Z_p2 Z_m2 dist_bins
%  clear tmp_arr2
end
dist_bins=dist_bins(1:length(dist_bins)-1);
save eee2

% remove data points with no data counts 
%nthresh=3
nthresh=10
nthresh=5
dist_bins=dist_bins(cnt_distBins>nthresh);
mean_distBins=mean_distBins(cnt_distBins>nthresh);
var_distBins=var_distBins(cnt_distBins>nthresh);
cnt_distBins=cnt_distBins(cnt_distBins>nthresh);

% smoothed values
supsmuL=1;
if supsmuL
  mean_distBins_sm=supsmu(dist_bins,mean_distBins);
  var_distBins_sm=supsmu(dist_bins,var_distBins);
else
% fixed smoothing
%save dsupsmu
  spanv=1e-3
  mean_distBins_sm=supsmu(dist_bins,mean_distBins,'span',spanv);
  var_distBins_sm=supsmu(dist_bins,var_distBins,'span',spanv);
end

%
if plotVals
%  dist_bins
  figure
  plot(dist_bins,mean_distBins,'k-'), hold on
  plot(dist_bins,mean_distBins+sqrt(var_distBins),'k-')
  plot(dist_bins,mean_distBins-sqrt(var_distBins),'k-')
  plot(dist_bins,mean_distBins_sm,'r-'), hold on
  plot(dist_bins,mean_distBins_sm+sqrt(var_distBins_sm),'r-')
  plot(dist_bins,mean_distBins_sm-sqrt(var_distBins_sm),'r-')
  set(gca,'XScale','log'),
end


end
%----------------------------------------------------------------------
