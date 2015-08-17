function [likelihood, gradient] = rho_obj(ww, cultureddata, model)

rho=exp(ww);
logbeta = model.log_ref_profiles + repmat(ww',size(model.log_ref_profiles,1),1);
logbeta = logbeta -repmat(logsum(logbeta,2),1,size(logbeta,2));

G = size(model.log_ref_profiles,2);
K = size(model.theta,2);
D = size(cultureddata,2);

likelihood = sum((model.kappa-1).*log(rho) - (model.kappa.*rho));

for dd=1:D
    log_P_t_given_theta = logsum(   repmat( log(model.theta(dd,:)), G, 1)' + logbeta, 1);
    likelihood = likelihood + (log_P_t_given_theta * cultureddata(:,dd));
end

gradient = (model.kappa-1)./rho - model.kappa + sum(cultureddata./repmat(rho,1,D),2);

%pi_k_g = \frac{\gamma_{k,g}}{\sum_{g'}\gamma_{k,g'}rho_{g'}}
pi_k_g = exp(    model.log_ref_profiles - repmat(logsum(model.log_ref_profiles + repmat(ww',K,1),2),1,G)   );
MM_d_g = cultureddata' ./ (model.theta * pi_k_g);

for dd=1:D
    V_g_g = pi_k_g' * diag(model.theta(dd,:)) * pi_k_g;
    gradient=gradient - (V_g_g * (MM_d_g(dd,:)'));
end

gradient = gradient .* rho;

likelihood = -likelihood;

gradient = -gradient;
