function likelihood = compute_likelihood(cultureddata, model)

D = size(cultureddata,2);
G = size(model.log_prob_w,2);

likelihood = sum((model.kappa-1).*log(model.rho) - (model.kappa.*model.rho)) - G*gammaln(model.kappa) + (G*model.kappa*log(model.kappa));

likelihood = likelihood + D*(gammaln(sum(model.alpha)) - sum(gammaln(model.alpha))) + sum((model.alpha-1) * log(model.theta)');

for dd=1:D
    log_ptgt = logsum(repmat( log(model.theta(dd,:)), G, 1)' + model.log_prob_w, 1);
    likelihood = likelihood + (log_ptgt * cultureddata(:,dd));
end
