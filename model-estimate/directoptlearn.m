function [model likelihood] = directoptlearn(HeterogeneousProfiles, ReferenceProfiles)
%function [newmodel likelihood] = directoptlearn(HeterogeneousProfiles, ReferenceProfiles)
%
% INPUT:
%  HeterogeneousProfiles: a GxD matrix representing gene expression profiles of heterogeneous (mixed) samples, where G is the number of genes, D is the number of heterogeneous samples
%  ReferenceProfiles: a KxG matrix, representing the expression profiles of the K reference cell populations.  Note the dimensions (KxG)
%
% OUTPUT:
%
% likelihood: likelihood of the final model
%
% newmodel: a structure with the following important fields:
%
%  theta: a DxK matrix, giving the fractional composition of each heterogeneous sample.  Each row represents a heterogeneous sample that was part of the input, and the first K columns correspond to the fractional composition with respect to reference population K.
%  omega: a Mx1 vector describing the convex combination weights learned by NNMLnp over the reference profiles, that when applied to the reference profiles, forms the prior over the new population.
%  rho: the scaling factors that, when applied to each of the reference profiles (then re-normalized -- see paper for details), obtains the estimated expression profiles of the constituent populations.


NTOPICS=size(ReferenceProfiles,1);

INITIAL_ALPHA = rand(1,NTOPICS) + 1;

ReferenceProfiles = ReferenceProfiles ./ repmat(sum(ReferenceProfiles,2),1, size(ReferenceProfiles,2));

disp('-----------------');
disp('Initializing...');

INIT_MODEL = new_model(HeterogeneousProfiles, INITIAL_ALPHA, ReferenceProfiles);
[model likelihood] = optmodel(HeterogeneousProfiles, INIT_MODEL);
