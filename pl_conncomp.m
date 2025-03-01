function [clustersize,clusters,numcluster] = pl_conncomp(map)
%
% Finds clusters (connected components) in maps of any dimension. The clusters consist of map locations with non-zero elements. 
% Useful for cluster-based statistical procedures 
% 

% This function is part of the permutationlab software:
% Author: Dimitrios Pantazis
% The code is provided "as is," without any guarantees or warranties, and is available for unrestricted use.

if exist('bwconncomp') %if image processing toolbox exists
    
    C = bwconncomp(map); %compute connected components
    numcluster = C.NumObjects;
    clusters = C.PixelIdxList;
    clustersize = cellfun(@length,clusters);
    if isempty(clustersize)
        clustersize = 0;
    end
    
elseif exist('spm_bwlabel') & ~isvector(map) & ndims(map)<=3 %if spm is installed, and 2D, 3D data

    [L,numcluster] = spm_bwlabel(map); %compute connected components
    ndx = find(L);
    clustersize = histcounts(L(ndx),numcluster);
    
    if nargout > 1 %if cluster indices are also needed
        clusters = cell(1,numcluster);
        for c = ndx'
            clusters{L(c)}(end+1,1)=c;
        end
    end   
    
else %use Matlab slower code
    
    [clustersize,clusters,numcluster] = pl_conncomp_matlab(map,'maximal');    
    
end




function [clustersize,clusters,numcluster] = pl_conncomp_matlab(map,conn)


%initialize
N = ndims(map); %map number of dimensions
dn = size(map); %dimensions
nvoxels = prod(size(map)); %number of voxels
mapmask = map~=0; %find non-zero voxels, convert to logical (possible NaNs exist)
ndx = find(mapmask); %index of suprathreshold voxels
ndx(isnan(map(ndx)))=[]; %remove NaNs
nsupvoxels = length(ndx); %number of suprathreshold voxels
i = 0; %no clusters initially
clusters{1} = []; %initially empty cluster

%if no suprathreshold voxels exist
if isempty(ndx)
    numcluster = 0;
    clusters = {[]};
    clustersize = 0;
    return
end

%template of neighbouring coordinates, used to expand a cluster (expndx)
if N>3 %only needed for general N case (>3)
    expndx = rem(floor([0:3^N-1]' * 3.^(-N+1:0)),3) - 1; %build a list of neighbor indexes
    s = sum(abs(expndx),2);
    if strcmp(conn,'maximal')
        expndx(s == 0,:) = []; %remove same point
    elseif strcmp(conn,'minimal')
        expndx(s == 0 | s == N,:) = []; %remove points with no common edge (all dimensions change) or same point
    end
    expndx = expndx';
    connectivity = size(expndx,2);
end


%find clusters
for k = 1:nsupvoxels %for all supratheshold voxels
    
    %index of suprathreshold voxel
    svox_ndx = ndx(k);
    
    %if voxel already belongs to a cluster, continue to next iteration
    if ~mapmask(svox_ndx) %if map is 0
        continue
    end
    
    %create new cluster
    i = i+1;
    clustersize(i) = 1; %new cluster
    clusters{i} = {svox_ndx}; %new cluster voxel
    mapmask(svox_ndx) = 0; %found this voxel
    cluster_exp = svox_ndx; %initialize cluster expand variable
    
    %expand cluster
    while ~isempty(cluster_exp) %while the cluster keeps expanding
        
        if N==1 %fast implementation for N=1
            [x,y] = ind2sub(dn,cluster_exp(:));
            expndx = [x-1 x+1];
            ndxout = expndx<1 | expndx>dn(1); %find out of bound voxels            
            expndx(ndxout)=[]; %remove out of bound voxels
            cluster_exp = sub2ind(dn,expndx); %convert to single indexing scheme

        elseif N==2 %fast implementation for N=2
            [x,y] = ind2sub(dn,cluster_exp(:));
            expndx = [x-1 y-1; x-1 y; x-1 y+1; x y-1; x y+1; x+1 y-1; x+1 y; x+1 y+1];
            ndxout = expndx(:,1)<1 | expndx(:,1)>dn(1) | expndx(:,2)<1 | expndx(:,2)>dn(2); %find out of bound voxels
            expndx(ndxout,:)=[]; %remove out of bound voxels
            cluster_exp = sub2ind(dn,expndx(:,1),expndx(:,2)); %convert to single indexing scheme

        elseif N==3 %fast implementation for N=3
            [x,y,z] = ind2sub(dn,cluster_exp(:));
            expndx = [x-1 y-1 z-1; x-1 y-1 z; x-1 y-1 z+1; x-1 y   z-1; x-1 y   z; x-1 y   z+1; x-1 y+1 z-1; x-1 y+1 z; x-1 y+1 z+1; x y-1 z-1; x y-1 z; x y-1 z+1; x y   z-1; x y   z+1; x y+1 z-1; x y+1 z; x y+1 z+1; x+1 y-1 z-1; x+1 y-1 z; x+1 y-1 z+1; x+1 y   z-1; x+1 y   z; x+1 y   z+1; x+1 y+1 z-1; x+1 y+1 z; x+1 y+1 z+1]; 
            ndxout = expndx(:,1)<1 | expndx(:,1)>dn(1) | expndx(:,2)<1 | expndx(:,2)>dn(2) | expndx(:,3)<1 | expndx(:,3)>dn(3); %find out of bound voxels
            expndx(ndxout,:)=[]; %remove out of bound voxels
            cluster_exp = sub2ind(dn,expndx(:,1),expndx(:,2),expndx(:,3)); %convert to single indexing scheme
        
        else %generic (but slower) N case
            
            %convert cluster_exp to multiple indexing scheme (from single indexing scheme)
            v = cell(N,1);
            [v{:}] = ind2sub(dn,cluster_exp(:));
            
            %apply neighborhood template to voxels v (separately per dimension)
            for n = 1:N
                vexp{n} = bsxfun(@plus,v{n},expndx(n,:)); %apply in each dimension
                vexp{n} = vexp{n}(:); %convert to column
                ndx0{n} = find(~vexp{n})'; %indices of 0 points
                ndxout{n} = find(vexp{n} == dn(n)+1)'; %indices of points exceeding map dimension
            end
            ndx_rm = [ndx0{:} ndxout{:}]; %indices of all out of boundary voxels
            for n = 1:N
                vexp{n}(ndx_rm) = []; %remove out of boundary voxels
            end
            cluster_exp = sub2ind(dn,vexp{:}); %convert to single indexing
            
        end
        
        cluster_exp = cluster_exp(mapmask(cluster_exp)); %remove subthreshold voxels (or already found)
        cluster_exp = unique(cluster_exp); %keep unique voxels
        clustersize(i) = clustersize(i) + numel(cluster_exp); %update cluster size
        clusters{i}{end+1} = cluster_exp; %expand cluster, using multiple cells is fast
        mapmask(cluster_exp) = 0; %remove voxels from map (already found)
        
    end
end

%assemble clusters
numcluster = length(clusters);
for i = 1:length(clusters)
    clusters{i} = cat(1,clusters{i}{:});
end





