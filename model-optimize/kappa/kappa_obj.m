function [likelihood, gradient] = kappa_obj(log_kappa, docs, model)


kappa = exp(log_kappa);
G=size(model.log_prob_w,2);

likelihood = sum((kappa-1).*log(model.rho) - (kappa.*model.rho)) - G*gammaln(kappa) + (G*kappa*log(kappa));

gradient = sum(log(model.rho) - model.rho) - G*psi(kappa) + G*log(kappa) + G;

gradient = gradient * exp(log_kappa);

likelihood = -likelihood;

gradient = -gradient;
