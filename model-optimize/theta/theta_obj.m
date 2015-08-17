function [myalhood, deriv_alhood] = theta_obj(ww, cultureddata, doc_index, model)

expww = exp(ww');
theta = expww ./ sum(expww);

myalhood = (model.alpha-1) * log(theta');
W = size(model.log_prob_w,2);
K = size(model.log_prob_w,1);

doccounts = cultureddata(:,doc_index);

log_P_t_given_theta = logsum(   repmat( log(theta), W, 1)' + model.log_prob_w    , 1);
myalhood = myalhood + (log_P_t_given_theta * doccounts);


dLdtheta = ((model.alpha-1)./theta) +   sum(  exp(model.log_prob_w - repmat(log_P_t_given_theta, K,1))' .* repmat(doccounts,1,K)         );


bb = expww ./sum(expww);
deriv_alhood = real(exp(logsum(log(dLdtheta) + log(bb),2) + log(bb)));
deriv_alhood = -deriv_alhood + (dLdtheta.*bb);

deriv_alhood = deriv_alhood';

deriv_alhood(1) = 0;

myalhood = -myalhood;
deriv_alhood = -deriv_alhood;
