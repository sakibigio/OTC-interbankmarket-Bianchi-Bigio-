function [RF_a] = averagerate_cobbdouglas(theta_0, delta_r, matchtech)
%AVERAGERATE_COBBDOUGLAS Compute average OTC rate for Cobb-Douglas matching
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
[~, CHIp_0] = analytic_cobbdouglas(0, theta_0, delta_r, matchtech);

% Get matching probability
[~, Gammas] = probs_cobbdouglas(theta_0, matchtech);

% Compute average rate
RF_a = CHIp_0 / Gammas;

end
