function [model likelihood] = optmodel(cultureddata, model)

likelihood_old = -Inf;

NUM_GRID_SEARCH_ITERATIONS=5;
NUM_ITERATIONS=35;

for iter=1:NUM_ITERATIONS
	if (iter <= NUM_GRID_SEARCH_ITERATIONS)
		NUM_ITERATIONS_RMINIMIZE=5;
	else
		NUM_ITERATIONS_RMINIMIZE=20;
	end
	
	disp('--- optimizing rho ...');
	model = opt_rho(cultureddata, model, NUM_ITERATIONS_RMINIMIZE, iter, NUM_GRID_SEARCH_ITERATIONS);
       
	disp('--- optimizing theta...');
	model = opt_theta(cultureddata, model, NUM_ITERATIONS_RMINIMIZE, iter, NUM_GRID_SEARCH_ITERATIONS);

	disp('--- optimizing alpha...');
	model = opt_alpha(cultureddata, model, NUM_ITERATIONS_RMINIMIZE, iter, NUM_GRID_SEARCH_ITERATIONS);

	disp('--- optimizing kappa...');
	model = opt_kappa(cultureddata, model, NUM_ITERATIONS_RMINIMIZE, iter, NUM_GRID_SEARCH_ITERATIONS);
 
	likelihood = compute_likelihood(cultureddata, model);
	change_ll = likelihood - likelihood_old;

	disp(['iter: ' num2str(iter) '/' num2str(NUM_ITERATIONS) ', likelihood: ' num2str(likelihood) ', change: ' num2str(change_ll)]);
    
	likelihood_old = likelihood;

end	
