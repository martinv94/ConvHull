function pnts = InitialPoints(polyt)
    global tol_mem;
    global tol;
    global n_dec_p;
    
    [A,b,Aeq,beq,lb,ub,dims] = polyt{:};
    nA = max(size(A,2),size(Aeq,2));
    h = zeros(1,nA); 
    h(dims(1)) = 1;
    % maximum
    [xopt,~,sol_flag] = linprog(-h,A,b,Aeq,beq,lb,ub);
    if sol_flag ~= 1
    	error('Error. \nNo feasible solution found.')
    end
    xopt=round(xopt,n_dec_p); %fixed-point estimate 
    pnts = ExtremePoint(h,h*xopt,-1,polyt);
    % minimum
    [xopt,~,sol_flag] = linprog(h,A,b,Aeq,beq,lb,ub); %POINT2
    if sol_flag ~= 1
    	error('Error. \nNo feasible solution found.')
    end
    xopt=round(xopt,n_dec_p);
    ep = ExtremePoint(h,h*xopt,1,polyt);
    if ~ismembertol(ep(dims).', pnts(dims,:).', tol_mem,'ByRows', true, 'DataScale',1) %0 false, 1 true 
        % if the two points are not the same store them both 
        pnts = [pnts, ep];
    end   
    while size(pnts,2) <= length(dims)
        [h,h0] = GetHyperplane(pnts,dims);
        [xopt,~,sol_flag] = linprog(h,A,b,Aeq,beq,lb,ub);%find new points from hessian normal form POINT3
        if sol_flag ~= 1
            error('Error. \nNo feasible solution found.')
        end
        xopt=round(xopt,n_dec_p);
        hx =h*xopt;
        if abs(hx - h0) > tol %checks if the new point sits on the hyperplane or not ????
            ep = ExtremePoint(h,hx,1,polyt); % if no, then calculate the new extreme point 
            if ~ismembertol(ep(dims).', pnts(dims,:).', tol_mem,'ByRows', true, 'DataScale',1)
                %checks if the extreme point has already been calculated
                %so where would gather all points that are "similar" 
                 pnts = [pnts, ep];
            end
        else
            [xopt,~,sol_flag] = linprog(-h,A,b,Aeq,beq,lb,ub); %checks the other side of the plane if no feasible solution on the other side 
            if sol_flag ~= 1
                error('Error. \nNo feasible solution found.')
            end
            xopt=round(xopt,n_dec_p);
            %Why do we not need to check here if the point sits on the
            %hyperplance or not???? 
            ep = ExtremePoint(h,h*xopt,-1,polyt);
            if ~ismembertol(ep(dims).', pnts(dims,:).', tol_mem,'ByRows', true, 'DataScale',1)
                 pnts = [pnts, ep];
            end
        end
    end
end

