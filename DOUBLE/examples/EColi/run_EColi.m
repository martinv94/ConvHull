%This script requires gurobi to be set up for matlab 

% Load the Model

% % The cobratoolbox is used to load the SBML model (.xml)
% model = readCbModel('iJ01366.xml');
% dims = [7,36,124,138,208,293];
% Aeq = model.S;
% beq = zeros(size(Aeq,1));
% lbs = model.lb;
% ubs = model.ub;

%Load the stoichiometric matrix and domain of the model (instead of xml)
Aeq = load('ECcp_Aeq.txt'); %stoichiometric matrix
dom = load('ECcp_domain.txt');
dims = [1,3]; %dimensions for the reactions of interest 
lbs = dom(:,1); %lower bounds
ubs = dom(:,2); %upper bounds

% Run double-precision CHM implementation

cd ../../chm
p1 = 12; %numerical accuracy parameter 
p2 = 8; %numerical comparison parameter
tic
CH=computeCH(Aeq,lbs,ubs,dims,p1,p2,0,1);
toc

