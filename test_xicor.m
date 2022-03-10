clc;clearvars;close all;
%% Basic test: uncorrelated
x = rand(1,100);
y = rand(1,100);
xi = xicor(x,y);
my_assert(abs(xi) < 0.3);

%% Basic test: totally dependent
x = randn(1,100);
y = x.^2;
xi = xicor(x,y);
my_assert(xi > 0.9);

%% Break tests
testCase = matlab.unittest.TestCase.forInteractiveUse;
verifyError(testCase,@() xicor('0'),"err1:MoreInputsRequired")
verifyError(testCase,@() xicor('0',2),"err2:TypeError")
verifyError(testCase,@() xicor([1 2],[2 4 5]),"err3:IncorrectLength")
verifyError(testCase,@() xicor([1 2 3],[2 4 5],'symmetric',2),"err2:TypeError")

%% With - without ties comparison
n = 4000;
x = 1:n;
y = cos(x./30) + 4*randn(1,n);
[xi, p] = xicor(x,y);

idx = randi(n);
x2 = [x x(idx)];
y2 = [y y(idx)];
[xi2, p2] = xicor(x2,y2);

my_assert(abs(xi - xi2) < 0.1);

%% Theoretical p_value vs. null
n = 4000;
n_boot = 10000;

x = 1:n;
y = cos(x./30) + 5*randn(1,n);

[xi, p] = xicor(x,y);

x_boot = randn(n_boot, n);
y_boot = randn(n_boot, n);
xi_boot = nan(1,n_boot);
for i_boot=1:n_boot
   xi_boot(i_boot) = xicor(x_boot(i_boot,:),y_boot(i_boot,:));
end

p_boot = sum(xi_boot > xi)/n_boot;

subplot(121);plot(x,y);
subplot(122);
h = histogram(sqrt(n)*xi_boot);hold on;
plot(xi*sqrt(n)*[1 1],[0 max(h.Values)],'r','LineWidth',1.5);

my_assert(abs(p - p_boot) < 0.1);
%% Basic assert function
function my_assert(cond)
if cond
    disp('Verification passed.');
else
    warning('Verification failed');
end
end
