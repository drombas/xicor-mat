close all;clc;clearvars;

% x =  linspace(-10,10,200);
% y = x.^2 + randn(1,200);

n_iter = 1000;
N = 40;

x = randn(n_iter,N);
y = randn(n_iter,N);
for n=1:n_iter
xi(n) = xicor(x(n,:),y(n,:));
end

e = sqrt(n)*xi;
x = max(abs(e))*[-1:.01:1];
y = normpdf(x,0,sqrt(2/5));

histogram(e,'Normalization','pdf');
hold on
plot(x,y,'r');