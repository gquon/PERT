function likelihood = alpha_compute_likelihood(alpha, sum_log_theta, D)

likelihood = D * (gammaln(sum(alpha)) - sum(gammaln(alpha))) + ((alpha - 1)*sum_log_theta');
