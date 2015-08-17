function [model likelihood] = directoptlearn(cultureddata, referenceProfilesTranspose)
%function [model likelihood] = directoptlearn(cultureddata, referenceProfilesTranspose)
%
% INPUT:
%  cultureddata: a GxD matrix representing gene expression profiles of heterogeneous (mixed) tumor samples, where G is the number of genes, D is the number of tumor samples
%  referenceProfilesTranspose: a KxG matrix, representing the expression profiles of the healthy contaminants of the sample.  So for lung cancer samples for example, this would be a set of expression profiles of healthy lung tissues (not necessarily patient-matched).
% SiteOfOriginPanelRatesTranspose: a MxG matrix, representing the expression profiles whose convex combination form the prior over the purified cancer profile learned.  Most of the time, is just the same matrix as referenceProfilesTranspose.
%
% OUTPUT:
%
% likelihood: likelihood of the final model
%
% newmodel: a structure with the following important fields:
%
%  theta: a Dx(K) matrix, giving the fractional composition of each tumor sample.  Each row represents a tumor sample that was part of the input, and the first K columns correspond to the fractional composition with respect to the Source Panel contaminants.  The last column represents the fractional composition of the pure cancer cells.
%  omega: a Mx1 vector describing the convex combination weights learned by ISOLATE over the SiteOfOriginPanelRatesTranspose matrix, that when applied to the Site of Origin Panel, forms the prior over the purified cancer profile



NTOPICS=size(referenceProfilesTranspose,1);

INITIAL_ALPHA = rand(1,NTOPICS) + 1;

referenceProfilesTranspose = referenceProfilesTranspose ./ repmat(sum(referenceProfilesTranspose,2),1, size(referenceProfilesTranspose,2));

disp('-----------------');
disp('Initializing...');

INIT_MODEL = new_model(cultureddata, INITIAL_ALPHA, referenceProfilesTranspose);
[model likelihood] = optmodel(cultureddata, INIT_MODEL);
