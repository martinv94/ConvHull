%This script requires gurobi and the cobratoolbox

% Load the Model
model = readCbModel('ThreeSpcModel.xml');

model = changeObjective(model,'mue_comm');
FBAsolution = optimizeCbModel(model,'max');
FBAsolution.f
model = changeRxnBounds(model, 'mue_comm',FBAsolution.f*0.99, 'l');
model = changeRxnBounds(model, 'mue_comm',FBAsolution.f, 'u');

cd ../../chm

Aeq = model.S;
beq = zeros(size(Aeq,1));
dims = [306, 300, 310, 296, 312, 324]; 
%Form_MM, H2_MB, H2_MM, Ac_MB, CO2_in, CH_in
lbs = model.lb;
ubs = model.ub;
p1 = 12;
p2 = 8;

CH=computeCH(Aeq,lbs,ubs,dims,p1,p2,0,1);
