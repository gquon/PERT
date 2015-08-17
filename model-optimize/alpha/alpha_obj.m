function [myalhood, deriv_alhood] = alpha_obj(ww, sum_log_theta, D)

ww = ww';
a = exp(ww) + 1;

myalhood = D * (gammaln(sum(a)) - sum(gammaln(a))) + ((a - 1)*sum_log_theta');

deriv_alhood = (D * (psi(sum(a)) - psi(a)) +  sum_log_theta) .* exp(ww);

deriv_alhood = deriv_alhood';

myalhood = -myalhood;

deriv_alhood = -deriv_alhood;
