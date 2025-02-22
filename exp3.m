clc; clear all; close all;
d0=100;
d=[0 10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170 180 190 200 210 220 230 240 250];
n=[2.2 1.8 3 2.4 2.6 2.0 2.1 1.8 1.6 3.0 2.1 3.3];
f=[914 914 1500 900 1900 1300 4000 1300 1300 900 4000 1300];
s=[8.7 5.2 7 9.6 14.1 3 7 6 5.8 7 9.7 6.8];
for i=1:12
v(i)=3e8/f(i);
end
for i=1:12
PL_d0(i)=20*log10(v(i)/(4*pi*d0));
X(i)=s(i).* randn(size(PL_d0(i)))
end
for i=1:12
for j=1:26
PL(i,j)=PL_d0(i)+10*n(i)*log(d(j)/d0)+X(i);
end
end
d=0:10:250;
plot(d,PL);
xlabel('Distance(m)');
ylabel('path loss(db)');
title('Distance vs path loss');
legend('Retail store','Grocery store','Office-hard partition','Office-soft partition','Office-soft partition','Textile/chemical');