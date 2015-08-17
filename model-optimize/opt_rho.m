function model = opt_rho(cultureddata, model, NUM_ITERATIONS_RMINIMIZE, iter, NUM_GRID_SEARCH_ITERATIONS)

init_xx = log(model.rho);

[xx tmp num_iter] = rminimize(init_xx, 'rho_obj', NUM_ITERATIONS_RMINIMIZE,cultureddata, model);

model.rho = exp(xx);

logbeta= model.log_ref_profiles+ repmat(xx',size(model.log_ref_profiles,1),1);

model.log_prob_w = logbeta -repmat(logsum(logbeta,2),1,size(logbeta,2));

