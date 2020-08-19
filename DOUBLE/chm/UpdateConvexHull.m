function UpdateConvexHull(newPt,dims)
    % Given a new extreme point, compute all possible HP with the new EP.
    % Removes hyperplanes that lie inside the hull
    global chull;
    global ePoints;
    global tol_mem;
    global n_dec_c;
    global hull_ids;
    for i=1:size(chull,2)
        pts = chull{3,i};
        % check if the new EP is already a member of the set of EPs
        if ismembertol(newPt(dims).',pts(dims,:).',tol_mem,'ByRows', true, 'DataScale',1)
            continue   
        end
        if round(chull{2,i}{1}*newPt,n_dec_c) <= round(chull{2,i}{2},n_dec_c) 
             continue
        end
        for j=1:size(pts,2)
            v = pts;
            v(:,j) = newPt;
            [h,h0] = GetHyperplane(v,dims);
            if HPinChull(h,h0,v,chull,dims,1) || HPinChull(-h,-h0,v,chull,dims,1)
               continue
            end      
            eh = round(h*ePoints,n_dec_c);
            % Adding new hyperplanes to chull
            if max(eh) <= round(h0,n_dec_c) 
                 chull = [chull, {hull_ids; {h,h0}; v;1}];
                 hull_ids = hull_ids + 1;
            if min(eh) >= round(h0,n_dec_c)
                chull = [chull, {hull_ids; {-h,-h0}; v;1}];
                hull_ids = hull_ids + 1;
            end
            end
        end
    end
    % Remove HP laying inside the hull
    to_remove = [];
    for j=1:size(chull,2) 
       ec = round(chull{2,j}{1}*ePoints,n_dec_c);
       h0 = round(chull{2,j}{2},n_dec_c);
       if min(ec) < h0 && max(ec) > h0
            to_remove = [to_remove j];
       end
    end
    chull(:,to_remove) = [];
    rm = size(to_remove,2);
end
