function xopt = ExtremePoint(hyrp,hyp0,opti,polyt)
    % compute an extreme point for a given projection
    global n_dec_p;
    global tol;
    global bugreport;
    global adjustol;
    [A,b,Aeq,beq,lb,ub,dims] = polyt{:};
    nA = max(size(A,2),size(Aeq,2));
    A=[A;hyrp;-hyrp];
    b=[b;hyp0;-hyp0];
    h = zeros(1,nA); 
    h(dims) = ones(1,length(dims));
    [xopt,~,sol_flag,output] = linprog(opti*h,A,b,Aeq,beq,lb,ub);
    % include tolerance for computing the linear program
    if sol_flag ~=1 
        [A,b,Aeq,beq,lb,ub,~] = polyt{:};
        A=[A;round(hyrp,n_dec_p);-round(hyrp,n_dec_p)];
        b=[b;hyp0  + tol;-hyp0 + tol];
        [xopt,~,sol_flag,~,~] = linprog(opti*h,A,b,Aeq,beq,lb,ub);
        adj_tol = 1; %counter for adjusting the tolerance 
    end
    % if tolerance is not enough, then adjust tolerance if necessary
    if adjustol==1
        while sol_flag ~=1 
            % checks if rounding causes there to be no optimal solution 
            new_tol = tol/adj_tol;
            [A,b,Aeq,beq,lb,ub,~] = polyt{:};
            A=[A;round(hyrp,n_dec_p);-round(hyrp,n_dec_p)];
            b=[b;hyp0  + new_tol;-hyp0 + new_tol];
            [xopt,~,sol_flag,~,~] = linprog(opti*h,A,b,Aeq,beq,lb,ub);
            if bugreport==1 && adj_tol > 1
                fprintf("p2 is temporarily adjusted to %d (in ExtremePoint calculation).\n",new_tol)
            end 
            adj_tol = adj_tol*10;
        end
    end
    if sol_flag ~=1 && adjustol==0
        disp(output)
        error('Error. \nNo feasible solution found. Consider setting adjustol=1.\n')
    end
    xopt = round(xopt,n_dec_p);
end
