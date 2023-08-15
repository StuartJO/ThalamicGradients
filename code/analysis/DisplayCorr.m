function DisplayCorr(X, Y)
% A function to display Pearson's correlation coefficient and related statistics
% 
% Inputs:
%   X - First variable (column vector)
%   Y - Second variable (column vector)
%
% Outputs:
%   None (displays correlation information in the command window)
%
% Usage example:
%   DisplayCorr(X, Y);
%
%   This function calculates and displays Pearson's correlation coefficient,
%   p-value, and confidence interval for two input variables X and Y.

if size(X,2) ~= 1
    error('X must be a column vector');
end

if size(Y,2) ~= 1
    error('X must be a column vector');
end
% Check if X and Y are column vectors of the same length
if ~isvector(X) || ~isvector(Y) || length(X) ~= length(Y)
    error('X and Y must be column vectors of the same length.');
end

% Calculate correlation coefficients, p-value, and confidence intervals
[RHO, P, CIL, CIU] = corrcoef(X, Y, 'Rows', 'complete');

% Calculate degrees of freedom
DF = sum(~isnan(sum([X, Y], 2))) - 2;

% Display correlation information
disp(['Pearson''s r(', num2str(DF), ') = ', num2str(RHO(1, 2)), ...
    ', p = ', num2str(P(1, 2)), ', CI = [', num2str(CIL(1, 2)), ...
    ', ', num2str(CIU(1, 2)), '], two-tailed'])
end