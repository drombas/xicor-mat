clc;clearvars;close all;
%% Basic test: uncorrelated
x = rand(1,100);
y = rand(1,100);
xi = xicor(x,y);
if abs(xi) < 0.1
    disp('Verification passed.');
else
    warning('Verification failed');
end
%% Basic test: totally dependent
x = randn(1,100);
y = x.^2;
xi = xicor(x,y);
assert(xi > 0.9);
%% Break test inputs
testCase = matlab.unittest.TestCase.forInteractiveUse;
verifyError(testCase,@() xicor('0'),"err1:MoreInputsRequired")
verifyError(testCase,@() xicor('0',2),"err2:TypeError")
verifyError(testCase,@() xicor([1 2],[2 4 5]),"err3:IncorrectLength")
verifyError(testCase,@() xicor([1 2],[2 4 5],'symmetric',2),"err2:TypeError")


