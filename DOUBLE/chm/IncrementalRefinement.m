function IncrementalRefinement(polyt)
    % Initial convex hull is refined by maximizing/minimizing the hyperplanes
    % containing the extreme points until all the facets of the projection are terminal.
    global chull; % initial hull, is updated throughout 
    global ePoints;
    global n_dec_p;
    global n_dec_c;
    global tol_mem;
    global tol;
    global bugreport;
    global adjustol;

    [A,b,Aeq,beq,lb,ub,dims] = polyt{:};
    
    counter = 0;
    while sum(cell2mat(chull(end,:))) ~= 0 % # of HPs that are non-terminal
        non_term = find(cell2mat(chull(end,:))~=0);
        if adjustol==1 && counter >= 1
            tol_mem = tol_mem*10;
            tol = tol*10;
            if bugreport == 1
                fprintf("p2 changed to %d.\n",tol)
            end
        end
        % new extreme points are computed and HPs inside the hull are
        % removed
        nt_chull = find(cell2mat(chull(end,:))~=0); 
        all_ids = cell2mat(chull(1,:));
        nt_ids = all_ids(nt_chull);
        for nt_id=nt_ids
            nt = find(cell2mat(chull(1,:))==nt_id);
            %disp(nt)
            %disp(size(chull))
            if cell2mat(chull(end,nt)) ~=0
                h = chull{2,nt}{1};
                h0 = chull{2,nt}{2};
                % maximize HP
                [xopt,~,sol_flag,~,~] = linprog(-h,A,b,Aeq,beq,lb,ub);
                if sol_flag ~= 1
                    error('Error. \nNo feasible solution found for HP in CHULL.')
                end
                hx = round(h*xopt,n_dec_p);
                % If HP is terminal
                if abs(hx - h0) < tol
                    % mark the HP as terminal (=1 if non-terminal)
                    chull{end,nt} = 0;
                else
                % if HP not terminal, compute new EP
                    ep = ExtremePoint(h, hx,-1,polyt);
                    if ~ismembertol(ep(dims).',ePoints(dims,:).',tol_mem,'ByRows', true,'DataScale',1)
                        ePoints = [ePoints, ep];
                        % update CH with new EP
                        UpdateConvexHull(ep,dims);
                    end
                end
            end     
        end
        %disp(size(ePoints))
        % Reset back to the original tolerance if terminal HPs were found
        new_non_term = find(cell2mat(chull(end,:))==1);
        if counter >= 1 && sum(ismember(non_term,new_non_term)) ~= size(new_non_term,2)
            tol_mem = tol_mem/10;
            tol = tol/10;
            counter = 0;
            if bugreport == 1
                fprintf("p2 reset to %d.\n",tol)
            end
        end
        % Adjust the tolerance if no new EPs were found and the same set of
        % HPs are still non-terminal
        if adjustol==1 && sum(ismember(non_term,new_non_term)) == size(new_non_term,2)
            counter = counter + 1;
        end
        if adjustol==0 && sum(ismember(non_term,new_non_term)) == size(new_non_term,2)
            error('Error. \nNo feasible solution found. Consider setting adjustol=1.\n')
        end
    end
    
    % Remove HP laying inside the hull
    	to_remove = [] ;
        for j=1:size(chull,2)
           ec = round(chull{2,j}{1}*ePoints,n_dec_c);
           h0 = round(chull{2,j}{2},n_dec_c);
           if min(ec) < h0  && max(ec) > h0
                to_remove = [to_remove j];
                if bugreport == 1
                    disp('HP removed in Incremental Refinement.\n')
                end
           end
        end
        chull(:,to_remove) = [];
end
