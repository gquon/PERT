function model = opt_theta(cultureddata, model, NUM_ITERATIONS_RMINIMIZE, iter, NUM_GRID_SEARCH_ITERATIONS)

D = size(model.theta,1);

if ~isfield(model, 'theta_weights');
        model.theta_weights = log(model.theta);
end

for dd=1:D
	init_xx = model.theta_weights(dd,:)';
	[xx tmp num_iter] = rminimize(init_xx, 'theta_obj', NUM_ITERATIONS_RMINIMIZE,cultureddata, dd, model);
	model.theta(dd,:) = exp(xx)'/sum(exp(xx));
	model.theta_weights(dd,:) = xx';
end

