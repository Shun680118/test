clc
clear all
close all

% load the data
startdate = '01/01/1994';
enddate = '01/01/2022';
f = fred
JPN = fetch(f,'JPNRGDPEXP',startdate,enddate)
KOR = fetch(f,'NGDPRSAXDCKRQ',startdate,enddate)
jpn = log(JPN.Data(:,2));
kor = log(KOR.Data(:,2));
q = JPN.Data(:,1);

T = size(jpn,1);
U = size(kor,1);

% Hodrick-Prescott filter
lam = 1600;
A = zeros(T,T);

% unusual rows
A(1,1)= lam+1; A(1,2)= -2*lam; A(1,3)= lam;
A(2,1)= -2*lam; A(2,2)= 5*lam+1; A(2,3)= -4*lam; A(2,4)= lam;

A(T-1,T)= -2*lam; A(T-1,T-1)= 5*lam+1; A(T-1,T-2)= -4*lam; A(T-1,T-3)= lam;
A(T,T)= lam+1; A(T,T-1)= -2*lam; A(T,T-2)= lam;

% generic rows
for i=3:T-2
    A(i,i-2) = lam; A(i,i-1) = -4*lam; A(i,i) = 6*lam+1;
    A(i,i+1) = -4*lam; A(i,i+2) = lam;
end

% Hodrick-Prescott filter
lam = 1600;
B = zeros(U,U);

% unusual rows
B(1,1)= lam+1; B(1,2)= -2*lam; B(1,3)= lam;
B(2,1)= -2*lam; B(2,2)= 5*lam+1; B(2,3)= -4*lam; B(2,4)= lam;

B(T-1,T)= -2*lam; B(T-1,T-1)= 5*lam+1; B(T-1,T-2)= -4*lam; B(T-1,T-3)= lam;
B(T,T)= lam+1; B(T,T-1)= -2*lam; B(T,T-2)= lam;

% generic rows
for i=3:T-2
    B(i,i-2) = lam; B(i,i-1) = -4*lam; B(i,i) = 6*lam+1;
    B(i,i+1) = -4*lam; B(i,i+2) = lam;
end

tauJGDP = A\jpn;
tauAGDP = B\kor;

% detrended GDP
ytilde = jpn-tauJGDP;
ctilde = kor-tauAGDP;

% plot detrended GDP
dates = 1994:1/4:2022.1/4; zerovec = zeros(size(jpn));
figure
title('Detrended log(real GDP) 1994Q1-2022Q1'); hold on
plot(q, ytilde,'r', q, ctilde,'b')
datetick('x', 'yyyy-qq')
legend('JAPAN', 'KOR', 'Location', 'best')

% compute sd(y), sd(c), rho(y), rho(c), corr(y,c) (from detrended series)
ysd = std(ytilde);
csd = std(ctilde);
corryc = corrcoef(ytilde(1:T),ctilde(1:T)); corryc = corryc(1,2);

disp(['Standard deviation of detrended log real JAPANGDP: ', num2str(ysd),'.']); disp(' ')
disp(['Standard deviation of detrended log real KORGDP: ', num2str(csd),'.']); disp(' ')
disp(['Contemporaneous correlation between detrended log real JAPAN and KOR: ', num2str(corryc),'.']);