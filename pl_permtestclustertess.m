function stat = pl_permtestclustertess(data,vertices,faces,varargin)
%
% Performs cluster-size statistical inference on 'data' for difference against 0. 
% For a difference against any mean value m, use input data - m.
%
% Data from each observation are randomly multiplied by +-1 to create permutation samples. 
% Statistical maps are then created by either averaging across the permuted observations (statistic = 'mean') or 
% computing a t-statistic across the permuted observations (statistic = 'tstat'). 
% All permutation statistical maps are converted to pvalues and then the clusteralpha threshold (in pvalue units)
% is applied to each permutation sample to identify suprathreshold clusters. The empirical
% distribution of the size of maximum suprathreshold cluster is used to assign statistical
% significance to the clusters (also in pvalue units) of the original data.
%   
% stat = pl_permclustertest(data) uses default values
%
% stat = pl_permclustertest(data,'param1',val1,'param2',val2,...) specifies one or
% more of the following name/value pairs
%
%       Parameter           Value
%       'alpha'             A value between 0-1 (default 0.05), specifying the significance level.
%       'statistic'         'tstat' (default), 'mean'
%                           Defines the test statistic, either a t-statistic or the mean  
%       'numpermutation'    An integer value, the number of permutations (default 1000)
%       'tail'              'right' (default),'left' 
%                           Determines the alternative hypothesis mean > 0, mean < 0
%       'clusteralpha'      A value between 0-1 (default 0.05), specifying the cluster defining threshold
%       'clusterstatistic'  'maxsize', 'maxsum', 'wcm', determines the
%                           cluster statistic which can be equal to the number of voxels (maxsize),
%                           the sum of voxel values (maxsum), 
%                           or a weighted cluster mass (wcm, see 'Hayasaka & Nichols (NeuroImage, 2004)
%       'theta'             A scalar used when 'wcm' is requested (default 0.5). Determines whether the 
%                           cluster statistic is dominated by the cluster size (theta -> 0)
%                           or the largest voxel intensity in the cluster (theta -> 1)
%       'verbose'           Display iteration progress  
%
% stat  .statmap                    Statistical map (variable1 x variable2 x ...)
%       .criticalmap                Critical map (binary map, 1 significant, 0 non-significant
%       .clustermassdist            Distribution of the maximum cluster statistic
%       .clusters                   Indices of significant clusters
%       .clustermass                Mass of significant clusters
%       .clusteralphamap            Statistic values at cluster defining threshold (e.g. when clusteralpha = 0.05)
%       .tail                       Same as input or default values (for bookkeeping)
%       .alpha                      Same
%       .statistic                  Same
%       .numpermutation             Same
%       .clusteralpha               Same
%       .clusterstatistic           Same


% This function is part of the permutationlab software:
% Author: Dimitrios Pantazis
% The code is provided "as is," without any guarantees or warranties, and is available for unrestricted use.


%% parse inputs 

alpha               = pl_inputparser(varargin,'alpha',0.05,@(x) isscalar(x) && (x > 0) && (x < 1));
statistic           = pl_inputparser(varargin,'statistic','tstat',{'tstat','mean'});
numpermutation      = pl_inputparser(varargin,'numpermutation',1000,@(x) isscalar(x) && x>0 && x == round(x) );
tail                = pl_inputparser(varargin,'tail','both',{'both','right','left'});
clusteralpha        = pl_inputparser(varargin,'clusteralpha',0.05,@(x) isscalar(x) && (x > 0) && (x < 1));
clusterstatistic    = pl_inputparser(varargin,'clusterstatistic','maxsize',{'maxsum','maxsize','wcm'});
theta               = pl_inputparser(varargin,'theta',0.5,@(x) isscalar(x) && (x > 0) && (x < 1));
verbose             = pl_inputparser(varargin,'verbose',false);


%% convert input data to cell arrays (if not already)

if ~iscell(data) %convert data to cell (faster)
    data = pl_mat2cell(data);
end


%% assign statistic function

switch statistic
    case 'mean'
        statisticfun = @pl_cellmean;
    case 'tstat'
        statisticfun = @pl_celltstat;
end


%% initialize variables

numobservation = length(data); 
rng('shuffle'); %seed the random number generator based on the current time


%% warning for processing speed

if verbose & ~isvector(data{1}) & ~exist('bwconncomp') & ~exist('spm_bwlabel')
    if ndims(data{1})<=3
        warning('Installing a) the image processing toolbox, or b) the spm software can significantly increase speed for multidimensional data!');
    else
        warning('Installing the image processing toolbox can significantly increase speed for multidimensional data!');
    end
end


%% compute statistic map values at clusteralpha level (symmetric for right/left tail since it is constructed with permutations)

if verbose
    disp(['Computing statistic values at cluster-defining threshold (clusteralpha = ' num2str(clusteralpha) ')...']);
end
slicesize = 50000; %slice data to reduce memory demand
switch tail
    case {'right' 'left'}
        [clusteralphamap] = pl_permalphamap(data, 'numpermutation',numpermutation,'alpha',clusteralpha,'statistic',statistic,'verbose',verbose,'slicesize',slicesize);
    case 'both' %if both, use half the clusteralpha threshold to account for the two tails
        [clusteralphamap] = pl_permalphamap(data, 'numpermutation',numpermutation,'alpha',clusteralpha/2,'statistic',statistic,'verbose',verbose,'slicesize',slicesize);
end

%% compute clustermass of original data

statmap = statisticfun(data); %original statistic map
[clustermass,clusters,numcluster] = compute_clustermass(statmap,vertices,faces,clusterstatistic,tail,clusteralphamap,theta); %clustermass of statmap
if numcluster == 0 %if no clusters found, return
    stat.statmap = statmap;
    stat.criticalmap = zeros(size(stat.statmap),'single');
    stat.clustermassdist = []; %clustermass max distribution
    stat.clusters = [];
    stat.clustermass = 0;
    stat.clusteralphamap = clusteralphamap;
    stat.tail = tail;
    stat.alpha = alpha;
    stat.statistic = statistic;
    stat.numpermutation = numpermutation;
    stat.clusteralpha = clusteralpha;
    stat.clusterstatistic = clusterstatistic;
    if strcmp(clusterstatistic,'wcm')
        stat.theta = theta;
    end
    return
end



%% compute clustermass max distribution using permutations

if verbose
    disp('Computing clustermass distribution using permutations...');
end
clustermassdist = zeros(numpermutation,1,'single');
clustermassdist(1) = max(clustermass); %first permutation is original data
for i = 2:numpermutation %leave one permutation for original data
    %verbose output
    if verbose & ~rem(i,verbose)
        disp(['Permutation: ' num2str(i) ' out of ' num2str(numpermutation)]);
    end
    permndx = sign(randn(numobservation,1,'single')); %create permutation index
    statmapperm = statisticfun(data,permndx); %create permutation sample
    [clustermassperm] = compute_clustermass(statmapperm,vertices,faces,clusterstatistic,tail,clusteralphamap,theta);
    clustermassdist(i) = max(clustermassperm); %distribution of maximum clustermass
end


%% compute clustermass threshold and significant clusters

clustermassdist = sort(clustermassdist); 
clustermassth = clustermassdist(ceil(numpermutation*(1-alpha)));
criticalclusters = clustermass>clustermassth;
clustermass = clustermass(criticalclusters); %keep only critical clusters
clusters = clusters(criticalclusters); %keep only critical clusters
criticalndx = cat(1,clusters{:}); %find critical locations
clusterpvalues = []; %initialize in case no critical clusters and clustermass is zero
for i = 1:length(clustermass)
    clusterpvalues(i) = sum( clustermass(i)<=clustermassdist )/numpermutation; %pvalues of each cluster
end

%% parse output

%clusters
stat.statmap = statmap; 
stat.criticalmap = zeros(size(stat.statmap),'single'); stat.criticalmap(criticalndx) = 1;
stat.clusters = clusters;
stat.clustermass = clustermass;
stat.clusterpvalues = clusterpvalues;
stat.clustermassdist = clustermassdist; %clustermass max distribution
stat.clusteralphamap = clusteralphamap;

% store input parameters
stat.tail = tail;
stat.alpha = alpha;
stat.statistic = statistic;
stat.numpermutation = numpermutation;
stat.clusteralpha = clusteralpha;
stat.clusterstatistic = clusterstatistic;
if strcmp(clusterstatistic,'wcm')
    stat.theta = theta;
end






function [clustermass,clusters,numcluster] = compute_clustermass(statmap,vertices,faces,clusterstatistic,tail,clusteralphamap,theta)

switch tail
    
    case 'right'
        
        statmap(statmap<=clusteralphamap) = 0; %thresholded statistic map 
        
        switch clusterstatistic %what type of clustermass to return
            case 'maxsize'
                [clustermass,clusters,numcluster] = pl_conncomptess(statmap,vertices,faces);
            case 'maxsum'
                [~,clusters,numcluster] = pl_conncomptess(statmap,vertices,faces);
                clustermass = zeros(1,numcluster,'single');
                for c = 1:numcluster
                    clustermass(c) = sum(statmap(clusters{c}));
                end
            case 'wcm'
                [~,clusters,numcluster] = pl_conncomptess(statmap,vertices,faces);
                clustermass = zeros(1,numcluster,'single');
                for c = 1:numcluster
                    clustermass(c) = sum( abs(statmap(clusters{c}) - clusteralphamap(clusters{c})).^(theta / (1-theta)) );
                end
        end
    
    
    case 'left'
        
        statmap(statmap>=-clusteralphamap) = 0; %thresholded statistic map 
        statmap = -statmap; %to return positive values of clustermass   
    
        switch clusterstatistic %what type of clustermass to return
            case 'maxsize'
                [clustermass,clusters,numcluster] = pl_conncomptess(statmap,vertices,faces);
            case 'maxsum'
                [~,clusters,numcluster] = pl_conncomptess(statmap,vertices,faces);
                clustermass = zeros(1,numcluster,'single');
                for c = 1:numcluster
                    clustermass(c) = sum(statmap(clusters{c}));
                end
            case 'wcm'
                [~,clusters,numcluster] = pl_conncomptess(statmap,vertices,faces);
                clustermass = zeros(1,numcluster,'single');
                for c = 1:numcluster
                    clustermass(c) = sum( abs(statmap(clusters{c}) - clusteralphamap(clusters{c})).^(theta / (1-theta)) );
                end
        end
    
    case 'both'
        
        statmap_left = statmap; %keep a copy
        statmap(statmap<=clusteralphamap) = 0; %thresholded statistic map
        
        switch clusterstatistic %what type of clustermass to return
            case 'maxsize'
                [clustermass,clusters,numcluster] = pl_conncomptess(statmap,vertices,faces);
            case 'maxsum'
                [~,clusters,numcluster] = pl_conncomptess(statmap,vertices,faces);
                clustermass = zeros(1,numcluster,'single');
                for c = 1:numcluster
                    clustermass(c) = sum(statmap(clusters{c}));
                end
            case 'wcm'
                [~,clusters,numcluster] = pl_conncomptess(statmap,vertices,faces);
                clustermass = zeros(1,numcluster,'single');
                for c = 1:numcluster
                    clustermass(c) = sum( abs(statmap(clusters{c}) - clusteralphamap(clusters{c})).^(theta / (1-theta)) );
                end
        end
        
        statmap_left(statmap_left>=-clusteralphamap) = 0; %thresholded statistic map
        statmap_left = -statmap_left; %to return positive values of clustermass       
        
        switch clusterstatistic %what type of clustermass to return
            case 'maxsize'
                    [clustermass_left,clusters_left,numcluster_left] = pl_conncomptess(statmap_left,vertices,faces);
            case 'maxsum'
                [~,clusters_left,numcluster_left] = pl_conncomptess(statmap_left,vertices,faces);
                clustermass_left = zeros(1,numcluster_left,'single');
                for c = 1:numcluster_left
                    clustermass_left(c) = sum(statmap_left(clusters_left{c}));
                end
            case 'wcm'
                [~,clusters_left,numcluster_left] = pl_conncomptess(statmap_left,vertices,faces);
                clustermass_left = zeros(1,numcluster_left,'single');
                for c = 1:numcluster_left
                    clustermass_left(c) = sum( abs(statmap_left(clusters_left{c}) - clusteralphamap(clusters_left{c})).^(theta / (1-theta)) );
                end
        end        
        
        %combine right and left tails
        clustermass = [clustermass clustermass_left];
        if exist('clusters')
            clusters = [clusters clusters_left];
            numcluster = length(clusters);          
        end

end

if isempty(clustermass)
    clustermass = 0; %never return empty
end


