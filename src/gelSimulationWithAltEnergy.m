% Gel simulation using the Alternative energy function


% *** HELPER FUNCTIONS ***

% For the simulation:

% poissonPoints
% laguerreTessellation
% skeleton
% gaussianSmoothing
% largestPeriodicComponent
% getStartingNodes
% getEdgeLengths
% getPadding
% intervalIndex
% getUpdatedAngles
% getAngles
% getProbabilityAlt1
% f1newEnergy
% closestCoordsVec
% getNewPosition
% rayTraversal
% deltaNodeI
% getProbabilityAlt2

% For the printing and plotting of statistics:

% printAcceptanceRates
% plotTracesAndTargets
% plotTracesUpdateSteps
% plotSimulatedNetwork
% densityPlot
% cumuPlot
% printLengthStatistics
% plotNodeDisplacements
% plotEdgeLengthsUpdateStep1
% plotEdgeLengthEvolution
% plotSimulatedNetworkSlice

% Toolboxes needed:

% Image Processing Toolbox
% Statistics and Machine Learning Toolbox

tic


% PARAMETER DEFINITIONS %

% GENERAL PARAMETERS %

% Random seed
seed = 1337;

% Number of iterations
nIterationsOuter = 30;
nIterationsInner1 = 10000;
nIterationsInner2 = 10000;

% Max edge length allowed
r = 40;

% Poisson parameter for the number of nodes 
% lambdaNodes = expected number of nodes in the simulation domain
lambdaNodes = 1604; %1341;%308;

% Standard deviation used in the truncated normal distribution in
% the second update step when a node relocation is proposed
sigmaNewPos = 50;


% DOMAIN AND NODE REGION PARAMETERS %

% Domain: [0,nx)x[0,ny)x[0,nz)
nx = 750;%427;%327;
ny = 750;%427;
nz = 125;%427;%128;


% Poisson parameter for the number of voids
% lambdaVoids = expected number of generator spheres (in total)
lambdaVoids = 100;%63;%15;


% Mean and std for the logNormal distribution from which the radii for the 
% generator spheres are sampled
muLogN = 1; 
sigmaLogN = 0.5;


% Gaussian smoothing parameter
sigmaSmooth = 1;%20;

% Threshold for the creation of the node region
threshold = 0.2;%0.02; 


% ENERGY FUNCTION PARAMETERS %

% Limits for the three intervals I1, I2 and I3
% I1 = (0 , l(1)]
% I2 = (l(1) , l(2)]
% I3 = (l(2) , r]
% Note that we must have 0 < l(1) < l(2) < r
l = zeros(2,1);
l(1) = 18;     
l(2) = 36;

% Number of edges in I1, I2 and I3, respectively.
% Should always be initialized to zero. 
nI = zeros(3,1);

% Target values used in the alternative energy function
propTarget = zeros(3,1);
propTarget(1) = 0.56; %0.24;  %0.56;   %2/5;
propTarget(2) = 0.358; %0.53;  %0.358;  %2/5;
propTarget(3) = 0.082; %0.23;  %0.082;  %1/5;

angleTarget = 2*pi/3;
degreeTarget = 2.74;

% Energy function tuning constants
% Values of c2 and c3 taken from Nisslert
c1 = zeros(3,1);
c1(1) = 4;     
c1(2) = 4;
c1(3) = 4;

c2 = 25;
c3 = 8;

% END OF PARAMETER DEFINITIONS %


rng(seed)


% CREATION OF NODE REGION

% Number of generator spheres for the tessellation
N = poissrnd(lambdaVoids); 

% Centre points for the generator spheres in the (interior of the) 
% domain (0,nx)x(0,ny)x(0,nz).
[X,Y,Z] = poissonPoints(nx,ny,nz,N);

% Radii for the generator spheres
R = lognrnd(muLogN,sigmaLogN,size(X));

% Number of regions for the tessellation
nRegionsX = 10;
nRegionsY = 10;
nRegionsZ = 10;

% Tessellation, creation of skeleton, smoothing and thresholding
disp('Tessellation performed...')
T = laguerreTessellation( ...
    nx,ny,nz, ...
    nRegionsX, ...
    nRegionsY, ...
    nRegionsZ, ...
    X,Y,Z, ...
    R);


disp('Skeleton created...')
S = skeleton(T);

disp('Smoothing performed...')
tempD = gaussianSmoothing(S,sigmaSmooth);

disp('Thresholding...')
D = tempD >= threshold;


% Finds the largest periodic connected component of D
[D, largestSize, prop] = largestPeriodicComponent(D);


if prop < 0.6
    error(['The created region is too disconnected (the greatest connected component only make up ', ...
           num2str(round(100*prop,2)), ...
           '% of total ones), please try again.'])
end

disp(['The greatest connected component (', ...
      num2str(round(100*prop,2)), ...
      '% of total ones) is used as Node region.'])


% INITIALIZATION OF NODES AND EDGES

% Initialization of the number of nodes and their positions
nNodes = poissrnd(lambdaNodes);
Nodes = getStartingNodes(D,nNodes);

% Initialization of the number of edges and the adjacency matrix
nEdges = zeros(nNodes,1); % Number of edges connected to each node
Edges = logical(sparse(nNodes, nNodes)); % Adjacency matrix

% intAngles - Stores all the angles for each node
intAngles = cell(1,nNodes);


% INITIALIZATION OF STATISTICS ARRAYS
% To calculate acceptance rates of the first and second update step,
% respectively
% Update step 1:
% Edge-removal in it. j is accepted: acc1(j) = -1
% Edge-removal in it. j is rejected: acc1(j) = -2
% Edge-addition in it. j is accepted: acc1(j) = 1
% Edge-addition in it. j is rejected: acc1(j) = 2
% Update step 2:
% Node relocation in it. j is accepted/rejected: acc2(j) = True/False
acc1 = zeros(nIterationsOuter*nIterationsInner1,1);
acc2 = false(nIterationsOuter*nIterationsInner2,1); 

% To store the lengths of the edges that are proposed to be shifted in the
% first update step
edgeLengths = zeros(nIterationsOuter*nIterationsInner1,1); 

% To store the mean degree, mean angle and mean edge length of the current
% simulated graph after each outer iteration. 
meanStatsOuter = zeros(nIterationsOuter+1,3); 

% To store the mean degree, mean angle and mean edge length of the current
% simulated graph every 100th inner iteration.
meanStatsInner = zeros(nIterationsOuter/100*(nIterationsInner1+nIterationsInner2),3);
count = 1;

% To store the t-values corresponding to the limits of the line segment, as
% well as the new node position proposal, in the second update step.
tCount = 1;
tValues = zeros(nIterationsOuter*nIterationsInner2,3);

% To store edge lengths for histogram plotting
histoCount = 1;
histoNum = min(6,nIterationsOuter);
histoIdx = floor(linspace(1,nIterationsOuter,histoNum));
histoLengths = cell(1,histoNum);


% PRINTING OF INITIAL VALUES
disp([newline 'Initial values'])
meanDegree = sum(nEdges)/nNodes;
disp("Mean degree: " + num2str(meanDegree))
meanAngle = 180/pi*mean([intAngles{:}]);
disp("Mean angle: " + num2str(meanAngle))
meanLength = getEdgeLengths(Nodes,Edges,nx,ny,nz,"mean"); 
disp("Mean edge length: " + num2str(meanLength))
disp("Edge length distribution: [0,0,0]")

% STORING OF STATISTICS
meanStatsOuter(1,:) = [meanDegree,meanAngle,meanLength];



% MCMC
% Outer loop
for n = 1:nIterationsOuter 

    % Creating periodic copies of the nodes
    [padding, realIndex] = getPadding(Nodes,r,0,nx,0,ny,0,nz);
    NodesAll = [Nodes;padding];

    % Creating k-d tree, neighborhoods and distances for all nodes in unit
    % cell
    Mdl = KDTreeSearcher(NodesAll);

    % Preallocations
    intNeighbors = cell(1,nNodes);
    intNeighborsIndex = cell(1,nNodes);

    % Number of nodes with no other nodes in its neighborhood
    noNbs = 0;

    [Idx,distance]=rangesearch(Mdl,Nodes,r);
    distance = distance';

    for i = 1:nNodes
        % The first element is deleted since it corresponds to node i 
        % itself.
        distance{i}(1) = []; 
        intNeighbors{i} = NodesAll(Idx{i}(2:end),:);  % Nodes in node i's neighborhood
        
        % Preallocation
        intNeighborsIndex{i} = zeros(length(distance{i}),1);

        if isempty(intNeighbors{i})
            noNbs = noNbs + 1;
        else
            for j = 1:size(intNeighbors{i},1) 
                [Lia,Locb] = ismember(intNeighbors{i}(j,:),Nodes,"rows");
                % If the neighbour is a "real" node (in unit cell)
                if (Lia == 1)  
                    intNeighborsIndex{i}(j) = Locb;
                % If the neighbour is a "virtual" node (from padding)
                else
                    [Lia, Locb] = ismember(intNeighbors{i}(j,:),padding,"rows");
                    tempIndex = realIndex(Locb);
                    intNeighborsIndex{i}(j) = tempIndex;  
                end
            end
        end
    end


    if noNbs == 1
        disp('1 node without neighbors')
    elseif noNbs > 1
        disp([num2str(noNbs), ' nodes without neighbors'])
    end 


    % First update step: Addition and removal of edges
    for i = 1:nIterationsInner1

        % First node
        randNode = randi(nNodes);
        nNbs = length(intNeighborsIndex{randNode});

        if (nNbs > 0)
            % Second node, in the neighbourhood of the first node
            randNeighborInd = randi(nNbs);
            randNeighbor = intNeighborsIndex{randNode}(randNeighborInd);
            dist = distance{randNode}(randNeighborInd);
            intIdx = intervalIndex(dist,l,r); 
            conn = Edges(randNode,randNeighbor); % Boolean indicating if 
                                                 % there is an edge between
                                                 % the two nodes

            % New angles of first node, if edge added/removed
            updatedAnglesNode = getUpdatedAngles( ...
                randNeighbor, ...
                randNeighborInd, ...
                conn, ...
                Edges(randNode,:), ...
                Nodes(randNode,:), ...
                intNeighborsIndex{randNode}, ...
                intNeighbors{randNode} ...
                );

            % New angles of second node, if edge added/removed
            updatedAnglesNeighbor = getUpdatedAngles( ...
                randNode, ...
                find(intNeighborsIndex{randNeighbor}==randNode), ...
                conn, ...
                Edges(randNeighbor,:), ...
                Nodes(randNeighbor,:), ...
                intNeighborsIndex{randNeighbor}, ...
                intNeighbors{randNeighbor} ...
                );

            totalEdges = nnz(Edges)/2;

            % Acceptance probability
            prob = getProbabilityAlt1( ...
                intIdx, ...
                propTarget, ...
                c1, ...
                totalEdges, ...
                nI, ...
                nEdges(randNode), ...
                nEdges(randNeighbor), ...
                [updatedAnglesNode,updatedAnglesNeighbor], ...
                [intAngles{randNode},intAngles{randNeighbor}], ...
                degreeTarget, ...
                angleTarget, ...
                c2,c3, ...
                conn);


            if (rand < prob)
                % The proposal is accepted!
                % Update of adjacency matrix, number of edges, angles and
                % acceptance rate array.

                if conn
                    nI(intIdx) = nI(intIdx) - 1;
                    acc1((n-1)*nIterationsInner1+i) = -1;  % Edge is removed
                    Edges(randNode,randNeighbor) = 0;
                    Edges(randNeighbor,randNode) = 0;
                    nEdges([randNode,randNeighbor]) = nEdges([randNode,randNeighbor]) - 1; 
                else
                    nI(intIdx) = nI(intIdx) + 1;
                    acc1((n-1)*nIterationsInner1+i) = 1;  % Edge is added
                    Edges(randNode,randNeighbor) = 1;
                    Edges(randNeighbor,randNode) = 1;
                    nEdges([randNode,randNeighbor]) = nEdges([randNode,randNeighbor]) + 1;
                end

                intAngles{randNode} = updatedAnglesNode;
                intAngles{randNeighbor} = updatedAnglesNeighbor;

            else
                % The proposal is rejected!
                % Update of acceptance rate array, rest needs not to be
                % updated. 
                if conn
                    acc1((n-1)*nIterationsInner1+i) = -2; % Edge-removal is rejected
                else
                    acc1((n-1)*nIterationsInner1+i) = 2; % Edge-addition is rejected
                end
            end
        end

        % STORING OF STATISTICS %
        edgeLengths((n-1)*nIterationsInner1+i)=dist;

        if mod(i,100) == 0
            meanDegree = sum(nEdges)/nNodes;
            meanAngle = 180/pi*mean([intAngles{:}]);
            meanLength = getEdgeLengths(Nodes,Edges,nx,ny,nz,"mean");
            meanStatsInner(count,:) = [meanDegree,meanAngle,meanLength];
            count = count+1;
        end

        % % % % % % % % % % % % %

    end



    % Second update step: Node relocations
    connectedNodesIndex = cell(nNodes,1);
    connectedNodesCoords = cell(nNodes,1);

    for j = 1:nNodes
        connectedNodesIndex{j} = find(Edges(j,:));
        connectedNodesCoords{j} = closestCoordsVec( ...
            Nodes(j,:), ...
            Nodes(connectedNodesIndex{j},:), ...
            nx,ny,nz ...
            );
    end

    for i = 1:nIterationsInner2
        % Node to relocate
        randNode = randi(nNodes);

        % Direction unit vector
        u = randn(1,3);
        u = u / norm(u);

        % Getting the new suggested node position
        % newPos: coordinates in unit box
        % q1 and q2 to use in Hastings ratio
        % ts = [-tNeg,tPos,t], i.e. the limits of the line segment
        % and the t-value corresponding to the new node position proposal
        [newPos,q1,q2,ts] = getNewPosition(Nodes(randNode,:),u,D,sigmaNewPos);

        % Getting the coordinates and distances of the connected nodes, wrt
        % periodic boundaries
        [newConnectedNodesCoords,newDists] = closestCoordsVec( ...
            newPos, ...
            Nodes(connectedNodesIndex{randNode},:), ...
            nx,ny,nz ...
            );

        % Boolean indicating if any edge is too long due to the node
        % relocation. If so, the proposal is rejected. 
        tooLongEdge = any(newDists>r);

        % If not, the acceptance probability is calculated
        if ~tooLongEdge
            % Preallocations
            % newAngle: Cell array where all new angles (from randNode and all 
            % its connected nodes) are stored.
            % newCoords: Cell array where the coordinates of randNode are
            % updated, for each of its connected nodes.
            newAngles = cell(1,nEdges(randNode)); 
            newCoords = cell(nEdges(randNode),1);

            for j = 1:nEdges(randNode)
                % The jth connected node of randNode
                tempNodeIndex = connectedNodesIndex{randNode}(j);
                % Coordinates in unit cell
                tempNodeCoords = Nodes(tempNodeIndex,:);
                % The jth connected node's connected nodes
                tempConnNodesIndex = connectedNodesIndex{tempNodeIndex};

                % Updating the coordinates of randNode
                tempRandNodeIndex = tempConnNodesIndex == randNode; 
                newPosCoords = closestCoordsVec(tempNodeCoords,newPos,nx,ny,nz);
                tempConnNodesCoords = connectedNodesCoords{tempNodeIndex};
                tempConnNodesCoords(tempRandNodeIndex',:) = newPosCoords;
                newCoords{j} = tempConnNodesCoords;


                % Getting the new angles of the jth connected node of
                % randNode
                newAngles{j} = getAngles(tempNodeCoords,tempConnNodesCoords);

              
            end
            
            % Adding the angles of randNode after relocation
            randNodeAngles = getAngles(newPos,newConnectedNodesCoords);
            allNewAngles = [newAngles{:},randNodeAngles];
            allOldAngles = [intAngles{[connectedNodesIndex{randNode},randNode]}];

            % Distances before relocation
            oldDists = vecnorm(connectedNodesCoords{randNode}-Nodes(randNode,:),2,2);

            
            deltaI = deltaNodeI(newDists,oldDists,l);
            totalEdges = sum(nEdges)/2;

            % Get acceptance probability
            prob = getProbabilityAlt2( ...
                nI, ...
                deltaI, ...
                c1, ...
                totalEdges, ...
                propTarget, ...
                allOldAngles, ...
                allNewAngles, ...
                angleTarget, ...
                c3, ...
                q1, q2 ...
                );
           

            
            if rand < prob
                % The proposal is accepted!
                % Update of the acceptance rates array and the node
                % coordinates
                acc2((n-1)*nIterationsInner2+i) = 1;
                Nodes(randNode,:) = newPos;

                % Update the connected nodes of randNode:
                connectedNodesCoords{randNode} = newConnectedNodesCoords;

                % Update the connected nodes of the connected nodes of randNode:
                for m = 1:nEdges(randNode)
                    tempNodeIndex = connectedNodesIndex{randNode}(m);
                    connectedNodesCoords{tempNodeIndex} = newCoords{m};
                    intAngles{tempNodeIndex} = newAngles{m};
                end

                intAngles{randNode} = randNodeAngles;
                nI = nI + deltaI;
            end
        end


        % STORING OF STATISTICS %
 
        tValues(tCount,:) = ts;
        tCount = tCount + 1;

        if mod(i,100) == 0
            meanDegree = sum(nEdges)/nNodes;
            meanAngle = 180/pi*mean([intAngles{:}]);
            meanLength = getEdgeLengths(Nodes,Edges,nx,ny,nz,"mean");
            meanStatsInner(count,:) = [meanDegree,meanAngle,meanLength];
            count = count+1;
        end

        % % % % % % % % % % % % %
    end




    % PRINTING OF STATISTICS %
    disp([newline 'Iteration ' num2str(n)])
    meanDegree = sum(nEdges)/nNodes;
    disp("Mean degree: " + num2str(meanDegree))
    meanAngle = 180/pi*mean([intAngles{:}]);
    disp("Mean angle: " + num2str(meanAngle))
    meanLength = getEdgeLengths(Nodes,Edges,nx,ny,nz,"mean");
    disp("Mean edge length: " + num2str(meanLength))
    totalEdges = sum(nEdges)/2;
    props = round(nI/totalEdges,2);
    disp("Edge length distribution: ")
    disp("[ " + num2str(props(1)) + " , " + num2str(props(2)) + " , " + num2str(props(3)) + " ]")

    % % % % % % % % % % % % % %

    % STORING OF STATISTICS %

    if ismember(n,histoIdx)
        histoLengths{histoCount} = getEdgeLengths(Nodes,Edges,nx,ny,nz,"all");
        histoCount = histoCount + 1;
    end

    meanStatsOuter(n+1,:) = [meanDegree,meanAngle,meanLength];

    % % % % % % % % % % % % % % %


end


disp([newline 'Inputs:'])
disp(['Mean degree: ' num2str(degreeTarget)])
disp(['Mean angle: ' num2str(180/pi*angleTarget)])
disp(['Edge length distribution: [ ' num2str(round(propTarget(1),2)) ' , ' num2str(round(propTarget(2),2)) ' , ' num2str(round(propTarget(3),2)) ' ]' newline])

time = toc;

disp("Elapsed time is " + num2str(toc/60) + " minutes.")
disp(newline)

% PRINTING SECTION - UNCOMMENT IF YOU WANT TO SEE THE PLOTS AND THE STATISTICS

printAcceptanceRates(nIterationsOuter,nIterationsInner1,nIterationsInner2,acc1,acc2)

plotTracesAndTargets(nIterationsOuter,meanStatsOuter,degreeTarget,angleTarget,0,1)

plotTracesUpdateSteps(nIterationsOuter,nIterationsInner1,nIterationsInner2,meanStatsInner)

plotSimulatedNetwork(Nodes,Edges,nx,ny,nz,r)

densityPlot(meanStatsOuter);

cumuPlot(meanStatsOuter);

printLengthStatistics(tValues, acc2)

plotNodeDisplacements(tValues,acc2)

plotEdgeLengthsUpdateStep1(edgeLengths,nEdges,acc1)

plotEdgeLengthEvolution(histoNum,histoLengths,histoIdx)

plotSimulatedNetworkSlice(Nodes,Edges,nx,ny,nz,r,200,300,0,ny,0,nz);

