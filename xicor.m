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

% Initial checks
if nargin == 1
    error('err1:MoreInputsRequired','xicor requires at least two inputs');
end

parser = inputParser;
addRequired(parser,'x');
addRequired(parser,'y');
addOptional(parser,'symmetric',false)
addOptional(parser,'p_value_method','theoretical')
addOptional(parser,'n_perm',100)

parse(parser,x,y,varargin{:})

x = parser.Results.x;
y = parser.Results.y;
symmetric = parser.Results.symmetric;
p_value_method = parser.Results.p_value_method;
n_perm = parser.Results.n_perm;

if ~isnumeric(x) || ~isnumeric(y)
    error('err2:TypeError','x and y are must be numeric');
end

n = length(x);

if n ~= length(y)
    error('err3:IncorrectLength','x and y must have the same length');
end

if ~islogical(symmetric)
    error('err2:TypeError','symmetric must be true or false');
end

% Check for NaN values
is_nan = isnan(x) | isnan(y);

if sum(is_nan) == n
    warning('No points remaining after excluding NaN.');
    xi = nan;
    return
elseif sum(is_nan) > 0
    warning('NaN values encountered.');           
    x = x(~is_nan);
    y = y(~is_nan);
    n = length(x);
end

if n < 10
    warning(['Running xicor with only ', num2str(n),...
        ' points. This might produce unstable results.']);
end

[xi, r, l] = compute_xi(x, y);

if symmetric
    xi = (xi + compute_xi(y,x))/2;
end

% If only one output return xi
if nargout <= 1
    return
end

if strcmp(p_value_method, 'numeric') && symmetric==true
    p_value_method = 'permutation';
end

% Compute p-values (only valid for large n)
switch p_value_method
    case 'theoretical'
        if length(unique(y)) == n
            p = 1 - normcdf(sqrt(n)*xi,0,sqrt(2/5));                
        else
            u = sort(r);
            v = cumsum(u);
            i = 1:n;
            
            a = 1/n^4*sum((2*n -2*i +1) .* u.^2);
            b = 1/n^5*sum((v + (n - i) .* u).^2);
            c = 1/n^3*sum((2*n -2*i +1) .* u);
            d = 1/n^3*sum(l .* (n - l));
            
            tau = sqrt((a - 2*b + c^2)/d^2);
            
            p = 1 - normcdf(sqrt(n)*xi,0,tau);
        end
    case 'permutation'
        xi_perm = nan(1,n_perm);
        
        if symmetric
            for i_perm=1:n_perm
                xi_perm(i_perm) = compute_xi(x(rand_perm(n)),y);
            end
        else
            for i_perm=1:n_perm
                x_perm = x(randperm(n));
                xi_perm(i_perm) = (compute_xi(x_perm,y) + ...
                                   compute_xi(y,x_perm))/2;
            end            
        end
        
        p = sum(xi_perm > xi)/n_perm;
    otherwise
        error("Wrong p_value_method. Use 'numeric' or 'permutation'");        
end

function [xi, r, l] = compute_xi(x,y)

n = length(x);

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
    r = nan;
    l = nan;
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
    xi = 1 - n*sum(abs(diff(r)))/(2*sum((n - l) .* l));
end
