function [clustersize,clusters,numcluster] = pl_conncomptess(map,vertices,faces)
%
% Finds clusters (connected components) in maps of tesselation x linear or tesselation x linear x linear dimension. 
% The clusters consist of map locations with non-zero elements. 
% Examples of tesselation are a cortical surface or a sensor surface.
% Useful for cluster-based statistical procedures 
% 

% This function is part of the permutationlab software:
% Author: Dimitrios Pantazis
% The code is provided "as is," without any guarantees or warranties, and is available for unrestricted use.


%initialize
nVerts = size(map,1); %number of tess vertices
nDim2 = size(map,2); %size of linear dimension
nDim3 = size(map,3); %size of linear dimension
nvoxels = prod(size(map)); %number of voxels
mapmask = map~=0; %find non-zero voxels, convert to logical (but NaNs may exist)
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

VertConn = pl_tess_vertconn(vertices, faces);

%find clusters
for k = 1:nsupvoxels %for all supratheshold voxels
    
    %index of suprathreshold voxel
    svox_ndx = ndx(k);
    
    %if voxel already belongs to a cluster, continue to next iteration
    if ~mapmask(svox_ndx) %if map is 0
%        continue
    end
    
    %create new cluster
    i = i+1;
    clustersize(i) = 1; %new cluster
    clusters{i} = {svox_ndx}; %new cluster voxel (seed point)
    mapmask(svox_ndx) = 0; %found this voxel
    cluster_exp = svox_ndx; %initialize cluster expand variable with seed
    
    %expand cluster
    while ~isempty(cluster_exp) %while the cluster keeps expanding
        
        if nDim3==1 %if only one linear dimension
            
            [x,y] = ind2sub([nVerts,nDim2],[cluster_exp]); %convert cluster_exp to multiple indexing scheme (from single indexing scheme)
            %expand over tess dimension
            cluster_exp = [];
            for j = 1:length(x)
                cluster_exp{j} = VertConn{x(j)} + nVerts*(y(j)-1); %notice the use of vertices connectivity
            end
            %expand over dim2 dimension
            expndx = [x y-1; x y+1];
            ndxout = expndx(:,2)<1 | expndx(:,2)>nDim2; %find out of bound voxels
            expndx(ndxout,:)=[]; %remove out of bound voxels
            cluster_exp{end+1} = sub2ind([nVerts,nDim2],expndx(:,1),expndx(:,2)); %convert to single indexing scheme
            cluster_exp = cat(1,cluster_exp{:}); %combine expansion over tess and dim2
            
        else %if two linear dimensions
            
            [x,y,z] = ind2sub([nVerts,nDim2,nDim3],[cluster_exp]); %convert cluster_exp to multiple indexing scheme (from single indexing scheme)
            %expand over tess dimension
            cluster_exp = [];
            for j = 1:length(x)
                cluster_exp{j} = VertConn{x(j)} + nVerts*(y(j)-1) + nVerts*nDim2*(z(j)-1); %notice the use of vertices connectivity
            end
            %expand over dim2 and dim3 linear dimensions: y-1,y,y+1 and z-1,z,z+1 (and combinations)
            expndx = [x y-1 z-1; x y-1 z; x y-1 z+1; x y z-1; x y z+1; x y+1 z-1; x y+1 z; x y+1 z+1];
            ndxout = expndx(:,2)<1 | expndx(:,2)>nDim2 | expndx(:,3)<1 | expndx(:,3)>nDim3; %find out of bound voxels
            expndx(ndxout,:)=[]; %remove out of bound voxels
            cluster_exp{end+1} = sub2ind([nVerts,nDim2,nDim3],expndx(:,1),expndx(:,2),expndx(:,3)); %convert to single indexing scheme
            cluster_exp = cat(1,cluster_exp{:}); %combine expansion over tess and dim2 and dim3

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

%sort clusters to decending order (from largerst to smallest)
[~,I] = sort(clustersize,'descend');
clusters = clusters(I);
clustersize = clustersize(I);





function [VertConn,C] = pl_tess_vertconn(vertices, faces)
% TESS_VERTCONN: Computes vertices connectivity.
% 
% INPUT:
%     - Vertices   : Mx3 double matrix
%     - Faces      : Nx3 double matrix
% OUTPUT:
%    - VertConn: Connectivity structure
%    - C: Connectivity sparse matrix with dimension nVertices x nVertices. 
%                It has 1 at (i,j) when vertices i and j are connected.



% Check matrices orientation
if (size(vertices, 2) ~= 3) || (size(faces, 2) ~= 3)
    error('Faces and Vertices must have 3 columns (X,Y,Z).');
end

% Disable the stupid warnings in old Matlab versions
warning('off', 'MATLAB:conversionToLogical');

% Build connectivity matric
rowno = double([faces(:,1); faces(:,1); faces(:,2); faces(:,2); faces(:,3); faces(:,3)]);
colno = double([faces(:,2); faces(:,3); faces(:,1); faces(:,3); faces(:,1); faces(:,2)]);
data = ones(size(rowno));
n = size(vertices,1);
C = logical(sparse(rowno, colno, data, n, n));

%find nonzero elements and initialize
[rows,cols,vals] = find(C);
d = find(diff([cols' cols(end)+1]));
VertConn = cell(n,1); % one empty cell per vertex

%if there are vertices without any connections, do a slow loop
if length(d)~=n
    for i = 1:length(cols)
        VertConn{cols(i)} = [VertConn{cols(i)}; rows(i)];
    end
    return
end

%fast calculation of vertices connectivity (if no isolated vertices)
VertConn{1} = rows(1:d(1));
for i = 2:n
    VertConn{i} = rows((d(i-1)+1):d(i));
end

