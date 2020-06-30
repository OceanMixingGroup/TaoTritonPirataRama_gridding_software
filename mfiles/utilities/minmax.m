function y = minmax(x)

% Finds the minimum and maximum of a variable that has multiple dimensions.
% For example, for a 4 dimensional variable Y:
% minmax(Y) = [min(min(min(min(Y)))) max(max(max(max(Y))))]
%
% see also: MEAN, MAX, MIN, MEANALL, MAXALL

y = [min(min(min(min(x)))) max(max(max(max(x))))];