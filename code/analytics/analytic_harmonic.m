function [CHIm_t,CHIp_t,r_f_t,gammad_t,gammas_t,theta_t]=analyticharmonic(t,theta_0,delta_r,matchtech)
% Analytic solution for harmonic mean matching function (p = -1)
% Inputs:
%   t         : time points (can be vector)
%   theta_0   : initial market tightness
%   delta_r   : spread r^w - r^m (penalty rate minus reserve rate)
%   matchtech : structure with lambda (matching efficiency) and eta (bargaining power)

lambda = matchtech.lambda; % Matching Efficiency
eta = matchtech.eta;       % Bargaining Power of Borrower

% Compute constant K from initial condition
K = (theta_0 - 1)^2 / theta_0;

% Compute theta_t using the closed-form solution
% theta_t = (2 + K*e^(2*lambda*t) + sqrt((2 + K*e^(2*lambda*t))^2 - 4)) / 2
A = theta_0 / (1 - theta_0)^2;     % >0 since theta_0>0, theta_0~=1
a = A .* exp(-lambda * t);
s = sqrt(1 + 4 .* a);
if theta_0 < 1
    % stable "minus" branch (returns theta_0 at t=0)
    theta_t = (2 .* a) ./ (1 + 2 .* a + s);
elseif theta_0 > 1
    % "plus" branch (returns theta_0 at t=0)
    theta_t = (1 + 2 .* a + s) ./ (2 .* a);
else
    theta_t = theta_0 == 1;
end

a = A .* exp(-lambda);
s = sqrt(1 + 4 .* a);
if theta_0 < 1
    % stable "minus" branch (returns theta_0 at t=0)
    theta_1 = (2 .* a) ./ (1 + 2 .* a + s);
elseif theta_0 > 1
    % "plus" branch (returns theta_0 at t=0)
    theta_1 = (1 + 2 .* a + s) ./ (2 .* a);
else
    theta_1 = theta_0 == 1;
end

% Matching rates for harmonic mean: gamma(theta) = 2*theta/(1+theta)
gammas_t = lambda * theta_t ./ (1 + theta_t);      % gamma(theta)
gammad_t = lambda ./ (1 + theta_t);                % gamma(1/theta)

% % Matching probabilities (closed-form for harmonic case)
if theta_0 ~= 1
     % Case: theta_0 >= 1 (deficit market)
     PSIp = (theta_1 - theta_0) / (theta_1 - 1);
     PSId = (theta_1 - theta_0) / ((theta_0 - 1) * theta_1);
else
     % Case: theta_0 < 1 (surplus market)
     % Need to use negative sign in solution and adjust formulas
     PSIp = 1-exp(-lambda/2);
     PSId = 1-exp(-lambda/2);
end

% Liquidity yield slopes
CHIp_t = delta_r * ((theta_1 ./ theta_t).^eta) .* ...
         ((theta_t.^eta .* theta_1.^(1-eta) - theta_t) ./ (theta_1 - 1));
     
CHIm_t = delta_r * ((theta_1 ./ theta_t).^eta) .* ...
         ((theta_t.^eta .* theta_1.^(1-eta) - 1) ./ (theta_1 - 1));

% Average OTC rate spread (r^f - r^m)
r_f_t = delta_r * ((theta_1 ./ theta_t).^eta) .* ...
        ((theta_t.^eta .* theta_1.^(1-eta) - 1) ./ (theta_1 - 1) - ...
         eta * (theta_t - 1) ./ (theta_1 - 1));

end