function flag = HPinChull(h,h0,v,array,dims,flag_points)
% checks if hyperplane is in the convex hull
% function with a logical output flag
    global tol;
    flag = 0;
    for i=1:size(array,2)
        % check if h0 matches
        if abs(h0 - array{2,i}{2}) < tol
            % check if h matches
            if all(abs(h(dims) - array{2,i}{1}(dims)) < tol)
                if flag_points
                    ps1 = sort(reshape(v(dims,:),1,[]));
                    ps2 = sort(reshape(array{3,i}(dims,:),1,[]));
                    % check if generating points match
                    if all(abs(ps1 - ps2) < tol)
                        flag = 1;
                        break
                    end
                else
                    flag = 1;
                    break
                end
            end
        end
    end
end
