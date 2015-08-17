function model = new_model(data, INITIAL_ALPHA,SourcePanelRates)

model.alpha = INITIAL_ALPHA;
model.log_prob_w = log(SourcePanelRates); %log SourcePanelRates
model.log_ref_profiles = log(SourcePanelRates);
model.kappa = 1;

D = size(data,2);
K = length(model.alpha);

model.theta = rand(D,K);
model.theta = model.theta ./ repmat(sum(model.theta,2), 1, K);

model.rho = ones(size(data,1),1);
