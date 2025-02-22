%SUB URBAN
clc
clear all;
close all;
hte=input('enter the tx height');
hre=input('enter the rx height');
f=input('enter the frequency');
ame=43;
garea=9;
lam=3e8/f
eirp=60;
d=input('enter the distance between the antennas');
ghte=20*log10(hte/200)
ghre=20*log10(hre/3)
lf=-10*log10((lam^2)/((4*3.14*d)^2))
l50=lf+ame-ghte-ghre-garea
pr=(eirp-l50+0)
%plot distance vs path loss
figure(1)
ghte=20*log10(hte/200)
ghre=20*log10(hre/3)
d1=1000:100:50000
lf=-10*log10((lam.^2)./((4*3.14*d1).^2))
l50=lf+ame-ghte-ghre-garea
subplot(2,3,1)
plot(d1,l50)
xlabel('distance in m')
ylabel('path loss')
%rx height vs path loss
ghte=20*log10(hte/200)
hre1=0:10
ghre1=20*log10(hre1/3)
lf=-10*log10((lam.^2)./((4*3.14*d).^2))
l50=lf+ame-ghte-ghre1-garea
subplot(2,3,2)
plot(hre1,l50)
xlabel('height in m')
ylabel('path loss')
%distance vs power rxved
ghte=20*log10(hte/200)
ghre=20*log10(hre/3)
d2=1000:100:50000
lf=-10*log10((lam^2)./((4*3.14*d2).^2))
l50=lf+ame-ghte-ghre-garea
pr=(eirp-l50+0)
subplot(2,3,3)
plot(d2,pr)
xlabel('distance')
ylabel('power rx')
% rx height vs power rx
ghte=20*log10(hte/200)
hre2=0:10
ghre2=20*log10(hre2/3)
lf=-10*log10((lam^2)./((4*3.14*d).^2))
l50=lf+ame-ghte-ghre2-garea
pr=(eirp-l50+0)
subplot(2,3,4)
plot(hre2,pr)
xlabel('rx height')
ylabel('power rx')

% tx height vs power rx
hte1=0:100
ghte1=20*log10(hte1/200)
ghre=20*log10(hre/3)
lf=-10*log10((lam^2)./((4*3.14*d).^2))
l50=lf+ame-ghte1-ghre-garea
pr=(eirp-l50+0)
subplot(2,3,5)
plot(hte1,pr)
xlabel('rx height')
ylabel('power rx')
% tx height vs path loss
hte2=0:100
ghte2=20*log10(hte2/200)
ghre=20*log10(hre/3)
lf=-10*log10((lam^2)./((4*3.14*d).^2))
l50=lf+ame-ghte2-ghre-garea
pr=(eirp-l50+0)
subplot(2,3,6)
plot(hte2,l50)
xlabel('rx height')
ylabel('power rx')