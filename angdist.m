function angmat = angdist(eigvecs,modes,ploton)
%ANGDIST Computes pairwise angular distances between eigenvectors in a list
% Inputs: eigvecs - eigenvectors as columns with additional dimensions (Nchans x Nchans x Nmw/Nsubs x Niter...)

if ~exist('modes','var') || isempty(modes)
    modes = 1;
end

n_eigvecs = size(eigvecs,3);
eigvecs = eigvecs(:,modes,:);

angmat = zeros(n_eigvecs,n_eigvecs);

for i = 1:n_eigvecs
    for j = 1:n_eigvecs
        
        if i < j
            angmat(i,j) = subspace(eigvecs(:,:,i),eigvecs(:,:,j));
        end
    
    end
end

angmat = angmat' + angmat;

end

