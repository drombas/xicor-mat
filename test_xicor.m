%% Basic test: uncorrelated
x = rand(1,100);
y = rand(1,100);
xi = xicor(x,y);
assert(abs(xi) < 0.1);
%% Basic test: totally dependent
x = randn(1,100);
y = x.^2;
xi = xicor(x,y);
assert(xi > 0.9);
%% Break test:
% testCase = matlab.unittest.TestCase.forInteractiveUse;
% verifyError(testCase,@() xicor('0'),"add5:InputMustBeNumeric")