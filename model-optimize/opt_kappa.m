function model = opt_kappa(cultureddata, model, NUM_ITERATIONS_RMINIMIZE, iter, NUM_GRID_SEARCH_ITERATIONS)

kappa = model.kappa;

init_xx = log(kappa);

[xx tmp num_iter] = rminimize(init_xx, 'kappa_obj', NUM_ITERATIONS_RMINIMIZE,cultureddata, model);
    
model.kappa = exp(xx');
