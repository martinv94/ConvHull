function [h,h0] = GetHyperplane(pnts, dims)
     % Compute the Hessian Normal form of a set of points
     global n_dec_p;
     h = zeros(1,size(pnts,1));
     dis = -ones(size(pnts,2),1);
     C = horzcat(pnts(dims,:).',dis); 
     hess = null(C); %calculates the null space of C 
     h(dims) = round(hess(1:length(dims),1),n_dec_p); %gets values for up to dims 
     h0  = round(hess(end,1),n_dec_p); %gets remaining values 
end
    
