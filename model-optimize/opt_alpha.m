function model = opt_alpha(cultureddata, model, NUM_ITERATIONS_RMINIMIZE, iter, NUM_GRID_SEARCH_ITERATIONS)

alpha = reshape(model.alpha, length(model.alpha), 1);

init_ww = log(alpha-1);

sum_log_theta = sum(log(model.theta),1);

ww = rminimize(init_ww, 'alpha_obj', NUM_ITERATIONS_RMINIMIZE, sum_log_theta, size(cultureddata,2));

alpha = exp(ww)+1;

model.alpha = reshape(alpha, 1, length(alpha));

if (iter <= NUM_GRID_SEARCH_ITERATIONS)
    likelihood = alpha_compute_likelihood(model.alpha, sum_log_theta, size(cultureddata,2));
    
    for cancer_alpha=10:10:100
        newalpha = ones(size(model.alpha'));  newalpha(end) = cancer_alpha;
        
        [xx tmp num_iter] = rminimize(log(newalpha-1), 'alpha_obj', NUM_ITERATIONS_RMINIMIZE, sum_log_theta, size(cultureddata,2));
        newalpha = exp(xx')+1;
        
        newll = alpha_compute_likelihood(newalpha, sum_log_theta, size(cultureddata,2));
        
        if (newll > likelihood)
            disp('alpha grid search found better solution.');
            model.alpha = newalpha;
            likelihood = newll;
        end
    end
    
end

