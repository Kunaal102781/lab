% HATA MODEL
%tx height vs path loss
clc
clear all;
close all;
hre=input('enter the rx height:');
hte=30:200;
d=input('enter the distance btw tx and rx:');
f=input('enter the frequency:');
% MEDIUM-SMALL CITY (SUB URBAN AREA)
ahre_1=(1.1*log10(f)-0.7)*hre-(1.56*log10(f)-0.8);
L50_1=69.55+26.16*log10(f)-13.82*log10(hte)+(44.9-6.55*log10(hte))*log10(d)-ahre_1;
L50_dB_1=L50_1-(2*log10(f/28)*log10(f/28))-5.4;
% MEDIUM-SMALL CITY (RURAL/OPEN AREA)
ahre_2=(1.1*log10(f)-0.7)*hre-(1.56*log10(f)-0.8);
L50_2=69.55+26.16*log10(f)-13.82*log10(hte)+(44.9-6.55*log10(hte))*log10(d)-ahre_2;
L50_dB_2=L50_2-(4.78*log10(f)*log10(f))+18.33*log10(f)-40.98;
% LARGE CITY (URBAN AND SUBURBAN AREA)
if f<300000000
ahre_3=8.29*(log10(1.54*hre)).*(log10(1.54*hre))-1.1;
else
ahre_3=3.2*(log10(11.75*hre)).*(log10(11.75*hre))-4.97;
end
L50_3=69.55+26.16*log10(f)-13.82*log10(hte)+(44.9-6.55*log10(hte))*log10(d)-ahre_3;
L50_dB_3=L50_3-(2*log10(f/28)*log10(f/28))-5.4;
% tx height vs path loss
figure
plot(hte,L50_dB_1,'r',hte,L50_dB_2,'g',hte,L50_3,'b',hte,L50_dB_3,'m');
legend('L_5_0(dB)(medium-small city (suburban area))','L_5_0(dB)(medium-smallcity (rural/open area))','L_5_0(dB)(large city (rural area))','L_5_0(dB)(large city (suburban area))');
xlabel('tx height in m');
ylabel('path loss in dB');
title('HATA MODEL');
grid on;