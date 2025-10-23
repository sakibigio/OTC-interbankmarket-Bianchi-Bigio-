function [RF_a] = averagerate_leontief(theta_0, delta_r, matchtech)
%AVERAGERATE_LEONTIEF Compute average OTC rate for Leontief matching
%   Computes RF_a = CHIp_0 / Gammas
%
%   Inputs:
%       theta_0   : Initial market tightness
%       delta_r   : Interest rate spread (r^w - r^m)
%       matchtech : Matching technology parameters
%
%   Outputs:
%       RF_a      : Average OTC rate

% Get initial liquidity yield
[~, CHIp_0] = analytic_leontief(0, theta_0, delta_r, matchtech);

% Get matching probability  
[~, Gammas] = probs_leontief(theta_0, matchtech);

% Compute average rate
RF_a = CHIp_0 / Gammas;

end
