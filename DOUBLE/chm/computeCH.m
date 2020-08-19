function CH = computeCH(Aeq,lbs,ubs,dims,np,nc,br,at)
    % Computes the convex hull for production envelopes of metabolic network. 
    % Solution is the list of hyperplanes and set of extreme points of the Convex hull. 
    % Inputs are:
    % - Aeq: the Stoichiometric matrix
    % - lbs ubs : and constraints on the fluxes (upper and lower bounds)
    % - dims: indices of dimensions onto which the projection should be computed
    % Optional arguments:
    % - n_dec_p: general precision for arithmetic operations. Default is 12.
    % - n_dec_c: precision for value comparison >=,<=,=,etc. Default is 8.
    % - bugreport: if True then automatic changes in tolerance are shown
    % Outputs are:
    % - CH.eps contains the extreme points
    % - CH.hps contains the hyperplanes

    global chull; %hull parameters
    global ePoints; %computed extreme points
    global lp_count; %counts the no. of linear programs 
    global hull_ids; %used to call HPs as you would from a dictionary
    global tol; %input parameters
    global tol_mem; %input parameter
    global n_dec_c; %input parameter
    global n_dec_p; %input parameter
    global bugreport; %input parameter
    global adjustol; %input parameter

    if(nargin<4)
        error('Error using computeCH. Not enough input arguments.')
    end
    if(nargin == 4)
	% default values
        n_dec_p = 12;
        n_dec_c = 8;
        tol_mem = 1*10^-n_dec_c; %format used by the function uniquetol
        tol = power(10,-8); %tolerance parameter for equality constraints
        bugreport = 0;
        adjustol = 1;
    else
    % input values 
        n_dec_p = np;
        n_dec_c = nc;
        tol_mem = 1*10^-n_dec_c;
        tol = power(10,-n_dec_c);
        bugreport = br;
        adjustol = at;
    end

    num_mets= size(Aeq,1); %number of metabolites
    beq = zeros(num_mets,1);
    A = []; %no inequalities in the optimization problem 
    b = [];
    lp_count = 0; %count the number of LPs done
    hull_ids = 1; %ids of the HPs store in chull
    
    polytope = {A,b,Aeq,beq,lbs,ubs,dims};
    
    % Initial points
    ePoints  = InitialPoints(polytope); 
    % Initial Hull
    chull = InitialHull(ePoints,dims);
    % Hull refinement
    IncrementalRefinement(polytope);
    % return list of HPs
    ht = GetListHPs(dims);
    
    %Return output
    CH ={};
    CH.hps = ht;
    CH.eps = ePoints(dims,:).';
    
    %Display output
    fprintf('%d linear programs computed, %d extreme points found.\n',lp_count,size(CH.eps,1))
end