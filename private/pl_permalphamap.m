function [alphamap] = pl_permalphavalue(data, varargin)
%
% Returns the location-specific statistical map values that correspond to an 'alpha' threshold. The statistical map is created
% from a statistic applied on 'data'

% This function is part of the permutationlab software:
% Author: Dimitrios Pantazis
% The code is provided "as is," without any guarantees or warranties, and is available for unrestricted use.


%% parse inputs 

numpermutation      = pl_inputparser(varargin,'numpermutation',1000,@(x) isscalar(x) && x>0 && x == round(x) );
alpha               = pl_inputparser(varargin,'alpha',0.05,@(x) isscalar(x) && (x > 0) && (x < 1));
statistic           = pl_inputparser(varargin,'statistic','tstat',{'tstat','mean'});
verbose             = pl_inputparser(varargin,'verbose',true);
slicesize           = pl_inputparser(varargin,'slicesize',50000);


%% assign statistic function

switch statistic
    case 'mean'
        statisticfun = @pl_cellmean;
    case 'tstat'
        statisticfun = @pl_celltstat;
end


%% initialize variables

numobservation = length(data); 
dn = size(data{1});
rng('shuffle'); %seed the random number generator based on the current time


%% calculate alpha map

alphamap = zeros(size(data{1}));
n = prod(dn);
[sndx,ssize,nslices] = pl_sliceblocks(n,slicesize);
for s = 1:nslices
    
    %verbose output
    if verbose & ~rem(s,verbose)
        disp(['Slice: ' num2str(s) ' out of ' num2str(nslices)]);
    end
    
    dataslice = pl_cellslice(data,sndx{s});
    statmap = statisticfun(dataslice); %first permutation sample is original data
    statmapperm = zeros(ssize(s),numpermutation,'single'); 
    statmapperm(:,1) = statmap'; %original sample is included in the permutations 
    for i = 2:numpermutation
        permndx = sign(randn(numobservation,1,'single')); %create permutation index
        statmapperm(:,i) = statisticfun(dataslice,permndx); %create permutation sample
    end
    if exist('prctile') %if image processing toolbox installed
        alphamap(sndx{s}) = prctile(statmapperm',100*(1-alpha));
    else
        statmapperm = sort(statmapperm,2);
        alphamap(sndx{s}) = statmapperm(:,ceil(numpermutation*(1-alpha)));
    end

end

