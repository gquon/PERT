function likelihood = theta_compute_likelihood(theta, cultureddata, doc_index, model)

likelihood = (model.alpha-1) * log(theta');
W = size(model.log_prob_w,2);
K = size(model.log_prob_w,1);

doccounts = cultureddata(:,doc_index);

log_P_t_given_theta = logsum(   repmat( log(theta), W, 1)' + model.log_prob_w    , 1);
likelihood = likelihood + (log_P_t_given_theta * doccounts);
