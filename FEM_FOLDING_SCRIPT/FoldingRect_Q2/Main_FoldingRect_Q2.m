clear all
close all
clc

%% INPUT

% Number of elements in circumferential direction (15,30,45,60,75,90,105,120)
nelc = 60;
[nodaba eldaba cortex_elements subcortex_elements] = readAbaqusFile(nelc);
sizeY_cortex = 2;
pert = 0.005;            % pertubation of center nodes in cortex (normalized wrt cortical thickness)
pert_width = 0.05;        % width of perturbed area (normalized wrt width of plate)

% Compute sizes and numbers
sizeX = max(nodaba(:,2));
sizeY = max(nodaba(:,3));
sizeY_subcortex = sizeY-sizeY_cortex;
nel_cortex = length(cortex_elements);
nel_subcortex = length(subcortex_elements);
nel = nel_cortex+nel_subcortex;
nnode = size(nodaba,1);

% Element direction
elemDir_cortex = [1 0]';
elemDir_subcortex = 'vertical';       % random-horizontal-vertical-radial
rand('seed',0);

% Parameter
param = Parameters();
max_iter = 10;
max_iter_inc = 6;          % if niter smaller, increase time step by 1.25
max_step = 201;

% Growth and time
tEnd = 60;
dt0  = 5;
dtMin = 0.05;
dtMax = dt0;
th_rate = 0.01;      % dtheta/dt
Gs = 0*0.1;          % Gs for strain driven growth
Jcrit = 1.0;           % Jcrit for strain driven growth
lambdaCrit = 1.0;

% OutputVariables (needed for plotting)
outputVar.time=0;
outputVar.movDuration = 5;      % seconds
outputVar.tEnd = tEnd;
outputVar.sizeY = sizeY;
outputVar.sizeX = sizeX;

% Filenames
strDescr = 'Plate';
file_out        = ['ModelsCpp/',strDescr,'.dat'];
rw_option = 1;          % 0 - no file writing/reading;  % 1 - file writing 

% Add paths
addpath('../GeneralFunctions')

fprintf('Building model...')

%% Materials
matProp.E = 1000;
matProp.nu = 0.4;
matProp.alpha = 0;
matProp.tref = 0;
matProp.lambda = 0;
matProp.rho = 1;
materials(1) = Material('MAT_NEOH',matProp);
materials(1).localID = 12;
matProp.E = matProp.E/5;
materials(2) = Material('MAT_NEOH',matProp);
materials(2).localID = 23;

%% Properties
var.th_rate = th_rate; var.Gs = Gs; var.Jcrit=Jcrit; var.lambdaCrit = lambdaCrit;
properties(1) = Property('PSOLID11',materials(1),var);
properties(1).localID = 111;
properties(2) = Property('PSOLID21',materials(2),var);
properties(2).localID = 222;

%% Nodes

% Nodal locations
for id = 1:nnode
    xi = nodaba(id,2);
    yi = nodaba(id,3);
    
    if(abs(xi-sizeX/2)<= pert_width*sizeX/2 && yi>=sizeY_subcortex) 
        yi = yi+ pert* sizeY_cortex; 
    end
    
    nodes(id) = Node([nodaba(id,2), nodaba(id,3)]);
    nodes(id).localID = id;
    nodes(id).dofID = [(id-1)*param.numDofPerNode+1:id*param.numDofPerNode];
end

%% Elements
count = 0;
for id=1:nel_subcortex
    count = count+1;
    nID = eldaba(subcortex_elements(id),2:5);
    prop = properties(2);   % subcortex: strain driven growth
    
    % Get center location of element
    locN = [nodes(nID(1)).loc; nodes(nID(2)).loc; nodes(nID(3)).loc; nodes(nID(4)).loc];
    locAvg = mean(locN);
    xA = locAvg(1); yA = locAvg(2);
    elemDir = getElemDir(xA, yA, elemDir_subcortex);

    elements(count) = Element('CQUAD',nodes(nID),elemDir,prop);
    elements(count).localID = count;
end

for id=1:nel_cortex
    count = count+1;
    nID = eldaba(cortex_elements(id),2:5);
    prop = properties(1);   % cortex: strain driven growth
    elemDir = elemDir_cortex;

    elements(count) = Element('CQUAD',nodes(nID),elemDir,prop);
    elements(count).localID = count;
end

%% SPC

% clamp bottom nodes
count = 0;
clampID = find(abs(nodaba(:,3))<1e-6);
for k=1:length(clampID)
    count = count+1;
    spc(count) = SPC(nodes(clampID(count)),1:param.numDofPerNode,0);
end

% constrain sides in x,y direction
count = 0;
bcID = [find(abs(nodaba(:,2)+sizeX)<1e-6); find(abs(nodaba(:,2)-sizeX)<1e-6)];
for k=1:length(bcID)
    count = count+1;
    spc(end+1) = SPC(nodes(bcID(count)),1:param.numDofPerNode,0);
    %spc(count) = SPC(nodes(bcID(count)),1,0);
end

%% LOAD
loads =[];
fprintf('\t\t Done!\n\n')

fprintf(['Number of nodes:    \t\t',num2str(length(nodes)),'\n']);
fprintf(['Number of elements: \t\t',num2str(length(elements)),'\n']);
fprintf(['Number of dof:      \t\t',num2str(length(nodes)*param.numDofPerNode),'\n\n']);

%% CREATE OUTPUT FOLDER
folders = CreateOutputFolders();

%% COMPUTATIONS
NonLinFem

%an_deflection = loadY*sizeZ^3*12/3/matProp.E/(sizeX^4);
%sim_deflection = max(u(2:3:end));

%%% Plot results
%}
%{
%% Write cpp file
if(rw_option==1)
    loadcases(1).spcSetID = 5;
    loadcases(1).loadSetID = 9;
    loadcases(1).mpcSetID = 0;
    loadcases(1).ploadSetID = 0;
    loadcases(1).tempSetID = 0;



    FEMmodel.nnode = nnode;
    FEMmodel.nel = nel;
    FEMmodel.materials = materials;
    FEMmodel.properties = properties;
    FEMmodel.nodes = nodes;
    FEMmodel.elements = elements;
    FEMmodel.spc = spc;
    FEMmodel.loads = loads;
    FEMmodel.loadcases = loadcases;
    FEMmodel.solution = 'QuasiStatic';
    FEMmodel.tEnd = tEnd;
    FEMmodel.dt0  = dt0;
    FEMmodel.dtMin = dtMin;
    FEMmodel.dtMax = dtMax;

    % Add some paths
    addpath('../../GeneralFunctions');

    WriteCppFileNewFemConcept(FEMmodel,500,file_out,rw_option);
end
%}

%printFunctions.Nodes(nodes);
%printFunctions.Elements(elements);
%printFunctions.SPC(spc);
%printFunctions.LOAD(loads);