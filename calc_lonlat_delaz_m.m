function []=calc_lonlat_delaz_m(lon1_deg,lat1_deg,lon2_deg,lat2_deg)

%
close all

%
coef=pi/180;
lon1=lon1_deg*coef;
lon2=lon2_deg*coef;
lat1=lat1_deg*coef;
lat2=lat2_deg*coef;

%
cdelt=cos(lat1)*sin(lat2)*cos(lon2-lon1)+sin(lat1)*sin(lat2);
x=cos(lat1)*sin(lat2)-cos(lat2)*sin(lat1)*cos(lon2-lon1);
y=cos(lat2)*sin(lon2-lon1);
%
sdelt=sqrt(x^2+y^2);
delta=atan(sdelt/cdelt)/coef;

%
az=atan(y/x)/coef

%
x=cos(lat2)*sin(lat1)-cos(lat1)*sin(lat2)*cos(lon2-lon1);
y=-sin(lon2-lon1)*cos(lat1);

%
baz=atan(y/x)/coef
delta=delta*111.11

figure
plot(lon1_deg,lat1_deg,'bs'), hold on
plot(lon2_deg,lat2_deg,'rs')

