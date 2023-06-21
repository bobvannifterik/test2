% Bob van Nifterik - MS3 - TU Delft
% image: image NxN gray scale
% order: max order
% p,q: free parameter
% moments: matrix containing [0,max order] moments

% This implementation is based on the implementation of
% https://doi.org/10.48550/arXiv.2103.14799

function [moments] = JFM_all(image,order,p,q)
[theta,rho,image] = create_grid_polar(image);

moment_mat=zeros(order,order);
for n=0:1:order
    for m=-order:1:order

        Rad=JFM_basis(n,rho,p,q);       % get the radial polynomial
        scaling_w = (1/(2*pi));  % To statisfy orthogonality
        
        temp = double(image) .*Rad.*exp(-1j*m * theta);
         
        moment_mat(n+1,m+order+1)=sum(temp(:))*scaling_w*(4/(nnz(Rad)+1));   
    
    end

end

moments = moment_mat;
end



