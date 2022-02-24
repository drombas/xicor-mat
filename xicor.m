function [xi,p] = xicor(x,y,method)
%XICOR Computes Chaterjee's xi correlation between x and y variables
%
%   [xi,p] = xicor(x,y,method)
%   Returns the xi-correlationa as well as the corresponding p-value
%   for the pair of variables x and y.
%
%   Input arguments:
%  
%   'x'              Independent variable. Numeric 1D array.
%            
%   'y'              Dependent variable. Numeric 1D array.
%  
%   'method'         Method to be used to compute the correlation. Options
%                    Options are
%  
%   Output arguments:
%  
%   'xi'             Computed xi-correlation. In range [0,1]
%
%   'p'              Estimated p-value.
%
%
%   
%   Notes
%   -----
%   The xi-correlation is not simmetric. 
%
%
%   References
%   ----------
%   [1]  Sourav Chatterjee, A New Coefficient of Correlation, Journal of 
%   the American Statistical Association, 116:536, 2009-2022, 2021.
%   DOI: 10.1080/01621459.2020.1758115
%
%   [2] Zhexiao Lin* and Fang Han†, On boosting the power of Chatterjee’s
%   rank correlation, arXiv, 2021. https://arxiv.org/abs/2108.06828
%
%
%   Example
%   ---------      
%   % Compute the xi-correlation between two variables
%
%     x = linspace(-10,10,50);
%     y = x.^2 + randn(1,50); 
%     [xi, p] = xicor(x,y);
%     
%  
%   David Romero-Bascones, dromero@mondragon.edu
%   Biomedical Engineering Department, Mondragon Unibertsitatea, 2022

if nargin == 2
    method = 'original'; 
end

n = length(x);

if n ~= length(y)
    warning('x and y have different number of samples');
end

switch method
    case 'original'
        % Compute y ranks
        [~, si] = sort(y, 'ascend');
        r = 1:n;
        r(si) = r;
        
        % Reorder based on x
        [~, si] = sort(x, 'ascend');
        r = r(si);
        
        % Compute correlation
        xi = 1 - 3*sum(abs(diff(r)))/(n^2 - 1);
        
        % Compute p-values (only valid for large n)        
        p = 1 - normcdf(sqrt(n)*xi,0,sqrt(2/5));
        
    case 'modified'
        error("Not implemented yet. Use 'original'");
        
    otherwise
        error("Not supported method. Use 'original' or 'modified'");
end
