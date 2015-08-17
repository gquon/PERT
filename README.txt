Please read the following paper first to understand the notation:

PERT: a method for expression deconvolution of human blood samples from varied microenvironmental and developmental conditions
Wenlian Qiao, Gerald Quon, Elizabeth Csaszar, Mei Yu, Quaid Morris, Peter W. Zandstra

1.	Initialization of PERT parameters
The hidden variable \alpha of PERT was initialized such that each entry of the vector was set to one more than a random number drawn uniformly from [1,2]. The parameter \kappa was initialized to 1. The hidden variables \theta_d were initialized such that each component was assigned a random number drawn uniformly from [0,1], then all numbers re-scaled to sum to one. The model parameter \rho_g was all initialized to one. 

2.	Deconvolution protocol
The NNML, NNMLnp and PERT were written in Octave. Two inputs are required for each model: one matrix contains the heterogeneous profiles with genes as rows and samples as columns, and the other matrix is the reference profiles with genes as rows and samples as columns. Elements of the output vector \theta are fractions of mixed profiles attributed to each reference profiles. In order to get the fraction of reference population i that has N_i replicates of gene expression profiles, deconvolved \theta’s for the Ni replicates are summed. 

To run PERT:
>> load(‘gene_expression_data.mat’)
>> addpath(genpath(‘directory_of_PERT_model’))
>> [pert, loglikelihood] = directoptlearn(MixedProfiles, transpose(ReferenceProfiles));
