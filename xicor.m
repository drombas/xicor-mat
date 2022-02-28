function [xi,p] = xicor(x,y,varargin)
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
%
%   Name-value arguments:
%  
%   'method'         Method to be used to compute the correlation. Options
%                    Options are
%
%   'symmetric'      If true xi is computed as (r(x,y)+r(y,x))/2. Default is 
%                    false.
%  
%   Output arguments:
%  
%   'xi'             Computed xi-correlation.
%
%   'p'              Estimated p-value.
%
%
%   
%   Notes
%   -----
%   The xi-correlation is not symmetric. 
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
    symmetric = false; 
end

n = length(x);

if n ~= length(y)
    error('x and y have different number of samples');
end

% Check for NaN values
is_nan = isnan(x) | isnan(y);

if sum(is_nan) > 0
    warning('NaN values encountered');           
    x = x(~is_nan);
    y = y(~is_nan);
end

% Reorder based on x
[~, si] = sort(x, 'ascend');
y = y(si);

% Compute y ranks
[~, si] = sort(y, 'ascend');
r = 1:n;
r(si) = r;

% If no Y ties compute it directly
if length(unique(y)) == n
    xi = 1 - 3*sum(abs(diff(r)))/(n^2 - 1);
else
    % Get r (yj<=yi) and l (yj>=yi)
    l = n - r + 1;
    
    y_unique = unique(y);
    idx_tie = find(groupcounts(y)>1);
        
    for i=1:length(idx_tie)
        tie_mask = (y == y_unique(idx_tie));                
        r(tie_mask) = max(r(tie_mask))*ones(1,sum(tie_mask));    
        l(tie_mask) = max(l(tie_mask))*ones(1,sum(tie_mask));
    end    
    
    % Compute correlation
    xi = 1 - n*sum(abs(diff(r)))/(2*sum(l * (n - l)));
end

% If only one output return xi
if nargout == 1
    return
end

if nargout == 2
    % Compute p-values (only valid for large n)
    p = 1 - normcdf(sqrt(n)*xi,0,sqrt(2/5));        
end
