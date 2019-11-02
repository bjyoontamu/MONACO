function [outputAlignment] = MONACO(netIdList, inputFolder, idFlag, outFileName, alpha, alignmentConstructionStrategy)

% MONACO constructs multiple alignment of a set of biological networks
%         outputAlignment = MONACO (netIdList, inputFolder, idFlag, outFileName, alpha, alignmentConstructionStrategy)
% 
% MONACO Input and Output Description
%         Input   netIdList                       ID list of input networks.
%                 inputFolder                     Input folder.
%                 idFlag                          This flag is used to expedite the reading process of input files. 
%                                                 id_flag = 1, if every node id is represented as "networkID+number".
%                                                 id_flag = 0, otherwise.
%                 outFileName                     Output file name including path.
%                 alpha                           A balancing parameter between node-level similarity and topological similarity.
%                                                 alpha < 0.5 works well in most cases.
%                 alignmentConstructionStrategy   0: Many-to-Many mapping
%                                                 1: One-to-one mapping
%                                                 2: Maximum Weighted Bipartite Matching (Pairwise network alignment only)
%         output  outputAlignment                 Alignment results
% 
% Input File Formats
%         Suppose we have three PPI networks 'a', 'b', 'c'.
%         To run MONACO, the files listed below are required.
% 
%         Tab-separated undirected PPI network files: a.net, b.net, and c.net 
%                 Example)  a.net
%                         a1      a2
%                         a3	a1
%                         a4	a2
%                         a2	a3
% 
%         Node-level similarity score file for each network pair: a-b.sim, a-c.sim, b-c.sim
%                 Example) a-b.sim
%                         a1	b1	153
%                         a1	b3	55
%                         a1	b7	49
%                         a2	b3	444
%                         a3	b3	211
%                         a3	b4	122
%                         a4	b5	251
%                         a4	b8	71
%                 * Nodes in the first column must be the nodes of the first network in the file name.
%                   Similarly, nodes in the second column have to be the nodes in the second network of the file name.
% 
% Output Files format:
%         Each line of the output file corresponds to individual cluster
%                 Example) 
%                         a4 b5 c4
%                         a1 b1 b7 c1
%                         a2 b3 c2
%                         a3 b4 c3
% 
% Example:
%        alignment = MONACO({'a', 'b', 'c'}, 'test', 1, 'output.txt', 0.4, 1)
% 
% For more information on the algorithms, please see:
% 
% Hyun-Myung Woo and Byung-Jun Yoon (2019)
% MONACO: accurate biological network alignment through optimal neighborhood matching between focal nodes
% 
% Contact: bjyoon@ece.tamu.edu

[G, S, L, edges, nodes, nodeIndices, Sims] = readNetworks(inputFolder, netIdList, idFlag);

tic;
numberOfIterations = 5;
clusterSize = 10;

[NCS, S_n, T_n] = computeNodeCorrespondenceScore(G, S, numberOfIterations, alpha);

if alignmentConstructionStrategy == 2
    [val, m1, m2] = bipartite_matching(NCS{1, 2});
    timeComplexity = toc; 
    fod = fopen(outFileName,'w');
    for ind = 1:length(m1)
        if idFlag==0
            fprintf(fod,'%s %s \n', cell2mat(nodes{1}(m1(ind))), cell2mat(nodes{2}(m2(ind))));
            alignmentResult{ind} = { cell2mat(nodes{1}(m1(ind))), cell2mat(nodes{2}(m2(ind))) };
        else
            fprintf(fod,'%s %s\n',strcat(netIdList{1},num2str(m1(ind))),strcat(netIdList{2},num2str(m2(ind))));
            alignmentResult{ind} = { strcat(netIdList{1},num2str(m1(ind))), strcat(netIdList{2},num2str(m2(ind))) };
        end
    end
    fclose(fod);
    outputAlignment = alignmentResult;
else
    P = [];
    for l = 1: length(G)
        for m = l + 1: length(G)
            if ~isempty(NCS{l, m})
                P_temp = -1 * ones(length(find(NCS{l, m})), length(G) + 1);
                [c, r, v] = find(NCS{l, m});
                P_temp(:, l) = c;
                P_temp(:, m) = r;
                P_temp(:, end) = v;
                P = [P; P_temp];
            end
        end
    end
    
    NCS2 = P;
    
    if alignmentConstructionStrategy == 1
        AL = injectiveAlignmentConstruction(NCS2, nodeIndices, L);
    else
        AL = alignmentConstruction(NCS2, nodeIndices, L, clusterSize);
    end
    timeComplexity = toc; 

    fod = fopen(outFileName,'w');
    outputAlignment = cell(1, length(AL));
    for i = 1: length(AL)
        outputAlignment{i} = cell(1, size(AL{i}, 1));
        for j = 1: size(AL{i}, 1)
            if idFlag == 0
                fprintf(fod, '%s ', upper(nodes{AL{i}(j, 1)}{AL{i}(j, 2)}));
                outputAlignment{i}{j} = nodes{AL{i}(j, 1)}{AL{i}(j, 2)};
            else
                fprintf(fod, '%s%d ', netIdList{AL{i}(j, 1)}, AL{i}(j, 2));
                outputAlignment{i}{j} = strcat(netIdList{AL{i}(j, 1)}, num2str(AL{i}(j, 2)));
            end
            
        end
        fprintf(fod, '\n');
    end
    fclose(fod);
end

end

function [G, S, L, edges, nodes, nodeIndices, Sims] = readNetworks(inputFolder, netIdList, idFlag)

M = length(netIdList);
nets = cell(1,M);
nodes = cell(1,M);
edges = cell(1,M);
Sims = cell(M,M);
netIdList = lower(netIdList);

for i = 1: M
    fid = fopen(strcat(inputFolder,'/',upper(netIdList{i}),'.net'));
    if idFlag
        nets{i} = textscan(fid, strcat(netIdList{i}, '%d\t', netIdList{i}, '%d\t%s'), 1000000);
    else
        nets{i} = textscan(fid, '%s\t%s\t%s', 1000000);
    end
    fclose(fid);
    if isempty(nets{i}{3}{1})
        nets{i}(3) = [];
    end
end

for i = 1: M
    for j = i + 1: M
        fid = fopen(strcat(inputFolder, '/', netIdList{i}, '-', netIdList{j}, '.sim'));
        if idFlag
            Sims{i, j}= textscan(fid, strcat(netIdList{i}, '%d\t', netIdList{j}, '%d\t%f'), 1000000);
        else
            Sims{i, j} = textscan(fid, '%s\t%s\t%f', 10000000);
        end
        fclose(fid);
    end
end

if idFlag
    nnodes = cell(1,M);
    for i = 1: M
        nnodes{i} = max([nets{i}{1}; nets{i}{2}]);
        for j = i + 1: M
            nnodes{i} = max([nnodes{i}; Sims{i, j}{1}]);
        end
        for j = 1: i-1
            nnodes{i} = max([nnodes{i}; Sims{j, i}{2}]);
        end
    end
    for i = 1: M
        nodes{i} = 1: nnodes{i};
    end
    for i = 1: M
        if length(nets{i}) == 2
            edges{i} = [nets{i}{1}, nets{i}{2}, ones(length(nets{i}{1}), 1)];
        else
            edges{i} = [nets{i}{1}, nets{i}{2}, cell2mat(nets{i, j}{3})];
        end
    end
    for i = 1: M
        for j = i + 1: M
            Sims{i, j} = [Sims{i, j}{1}, Sims{i, j}{2}, Sims{i, j}{3}];
        end
    end
else
    for i = 1: M
        nodes{i} = unique([nets{i}{1}; nets{i}{2}]);
        for j = i + 1: M
            nodes{i} = unique([nodes{i}; Sims{i, j}{1}]);
        end
        
        for j = 1: i - 1
            nodes{i} = unique([nodes{i}; Sims{j, i}{2}]);
        end
        nodes{i} = sort_nat(nodes{i});
        ids = cell(1, M);
        for j = 1: M
            ids{j} = [];
        end
        
        for j = i + 1: M
            [~, tempList] = ismember(Sims{i, j}{1}, nodes{i});
            Sims{i, j}{1} = {tempList};
        end
        
        for j = 1: i - 1
            [~, tempList] = ismember(Sims{j, i}{2}, nodes{i});
            Sims{j, i}{2} = {tempList};
        end
    end
    
    for i = 1: M
        edges{i} = zeros(length(nets{i}{1}), 3);
        [~, edges{i}(:, 1)] = ismember(nets{i}{1}, nodes{i});
        [~, edges{i}(:, 2)] = ismember(nets{i}{2}, nodes{i});
        
        if length(nets{i}) == 2
            edges{i}(:, 3) = 1;
        else
            edges{i}(:, 3) = str2double(nets{i}{3});
        end
    end
    
    for i = 1: M
        for j = i + 1: M
            Sims{i, j}=[cell2mat(Sims{i, j}{1}), cell2mat(Sims{i, j}{2}), (Sims{i, j}{3})];
        end
    end
    
end

L = zeros(1, M);
nodeIndices = cell(1, M);
for i = 1: M
    L(i) = length(nodes{i});
    nodeIndices{i} = 1: length(nodes{i});
end

G = cell(1, M);
for i = 1: M
    ee = double(edges{i});
    Q = sparse(ee(:, 1), ee(:, 2), ee(:, 3), L(i), L(i));
    indices  = sub2ind(size(Q), ee(:, 2), ee(:, 1));
    Q(indices) = ee(:, 3);
    G{i} = Q;
end

S = cell(M, M);
for i = 1: M
    for j = i + 1: M
        S{i,j} = sparse(double(Sims{i, j}(:, 1)),double(Sims{i, j}(:, 2)),double(Sims{i, j}(:, 3)), L(i), L(j));
    end
end
end % ends of function

function [GG,SS,J,LL]=update(L,G,S)
M=length(L);
J = cell(1,M);
for i=1:M
    J{i}=1:L(i);
end
Rj=J;
for i=1:M
    for j=i+1:M
        Oj=find(sum(S{i,j})==0);
        Oi=find(sum(S{i,j}')==0);
        Rj{i}=intersect(Rj{i},Oi);
        Rj{j}=intersect(Rj{j},Oj);
    end
end

GG=cell(1,M);
for i=1:M
    J{i}=1:L(i);
    J{i}(Rj{i})=[];
    GG{i}=sparse(length(J{i}),length(J{i}));
    GG{i}=G{i}(J{i},J{i});
end
SS=cell(M,M);
for i=1:M
    for j=i+1:M
        SS{i,j}=sparse(length(J{i}),length(J{j}));
        SS{i,j}=S{i,j}(J{i},J{j});
    end
end

LL=zeros(1,M);
for i=1:M
    LL(i)=size(GG{i},1);
end

end

function [NCS, S_n, T_n] = computeNodeCorrespondenceScore(G, S, numberOfIterations, alpha)
NCS = cell(length(G), length(G));

for l = 1: length(G)
    for m = l + 1: length(G)
        Snormalized = S{l, m};
        J1 = spdiags(1./sum(Snormalized, 2), 0, size(Snormalized, 1), size(Snormalized, 1));
        J2 = spdiags(1./sum(Snormalized, 1)', 0, size(Snormalized, 2), size(Snormalized, 2));
        Snormalized = (J1 * Snormalized + Snormalized * J2)/2;
        
        T = S{l, m};
        prevT = T;
        for o = 1: numberOfIterations
            Ttemp = calculateTopologySimilarity(G(l), G(m), T, 1);
            J1 = spdiags(1./sum(Ttemp, 2), 0, size(Ttemp, 1), size(Ttemp, 1));
            J2 = spdiags(1./sum(Ttemp, 1)', 0, size(Ttemp, 2), size(Ttemp, 2));
            Ttemp = (J1 * Ttemp + Ttemp * J2)/2;
            T = Ttemp;
            if max(abs(T-prevT)) < 0.001
                break;
            end
            prevT = T;           
        end
        
        NCS{l, m} = alpha*Snormalized + (1-alpha)*T;
        T_n{l, m} = T;
        S_n{l, m} = Snormalized;
    end
end
end

function MMM = alignmentConstruction(C, J, L, n_max)

M = length(J);

[s, y]=sort(C(:, M+1), 'descend');

CC = cell(M, M);
for i = 1: M
    CC{i, i} = sparse(L(i), L(i));
end

for i = 1: M
    for j = i + 1: M
        f = find(C(:, i) > 0 & C(:, j) > 0);
        CC{i, j} = sparse(C(f, i), C(f, j), C(f, M+1), L(i), L(j));
        CC{j, i} = CC{i, j}';
    end
end

MM = cell(1, M);
for i = 1: M
    MM{i} = sparse(1, L(i));
end
kkkk = .8;

k = 1;
g_max = 1;
max_l = ones(1, M) * n_max;

EC = {};
missed = [];
ssss = [];

while (k <= length(y))
    h = C(y(k), 1: M);
    if C(y(k), M + 1) < 0
        k = k + 1;
        continue;
    end
    if (prod(h) ~= 0)
        F = find(h ~= -1);
        f1 = F(1);
        f2 = F(2);
        h1 = h(f1);
        h2 = h(f2);
        c1 = MM{f1}(h1);
        c2 = MM{f2}(h2);
        if ((c1 == 0) & (c2 == 0))
            MM{f1}(h1) = g_max;
            MM{f2}(h2) = g_max;
            g_max = g_max + 1;
            EC = [EC, [f1, h1; f2, h2]];
            
        elseif ((c1 > 0) & (c2 == 0))
            if length(find(EC{c1}(:, 1) == f2)) >= max_l(f2)
                k = k + 1;
                continue;
            end
            
            if ~isempty((find(EC{c1}(:, 1) == f2)))
                others = find((EC{c1}(:, 1) ~= f2));
                own = find((EC{c1}(:, 1) == f2));
                m_1 = zeros(1, length(own));
                m_2 = 0;
                for j = 1: length(others)
                    jj = others(j);
                    for i = 1: length(own)
                        ii = own(i);
                        m_1(i) = m_1(i) + CC{f2, EC{c1}(jj, 1)}(EC{c1}(ii, 2), EC{c1}(jj, 2));
                    end
                    m_2 = m_2 + CC{f2, EC{c1}(jj, 1)}(h2, EC{c1}(jj, 2));
                end
                
                if m_2 < (kkkk * mean(m_1))
                    k = k + 1;
                    continue;
                end
            end
            
            EC{c1} = [EC{c1}; f2, h2];
            MM{f2}(h2) = MM{f1}(h1);
            
        elseif ((c1 == 0) & (c2 > 0))
            if length(find(EC{c2}(:, 1) == f1)) >= max_l(f1)
                k = k + 1;
                continue;
            end
            
            if ~isempty((find(EC{c2}(:, 1) == f1)))
                others = find((EC{c2}(:, 1) ~= f1));
                own = find((EC{c2}(:, 1) == f1));
                m_1 = zeros(1, length(own));
                m_2 = 0;
                
                for j = 1: length(others)
                    jj = others(j);
                    for i = 1: length(own)
                        ii = own(i);
                        m_1(i) = m_1(i) + CC{f1, EC{c2}(jj, 1)}(EC{c2}(ii, 2), EC{c2}(jj, 2));
                    end
                    m_2 = m_2 + CC{f1, EC{c2}(jj, 1)}(h1, EC{c2}(jj, 2));
                end
                
                if m_2 < (kkkk * mean(m_1))
                    k = k + 1;
                    continue;
                end
            end
            
            EC{c2} = [EC{c2}; f1, h1];
            MM{f1}(h1) = MM{f2}(h2);
            
        elseif (c1 ~= c2)
            ec_temp = [EC{c1}; EC{c2}];
            flag = 0;
            [tt1, htt1] = hist(ec_temp(:, 1), 1: M);
            
            if max(tt1) > 1
                k = k + 1;
                continue;
            end
            
            for i = 1: M
                MM{i}(find(MM{i} == c1 | MM{i} == c2)) = g_max;
            end
            
            EC = [EC, [EC{c1}; EC{c2}]];
            EC{c1} = [];
            EC{c2} = [];
            g_max = g_max+1;
        end
    end
    
    k = k + 1;
end

MMM = cell(g_max, 1);
for j = 1: g_max
    for i = 1: M
        ff = find(MM{i} == j);
        if ~isempty(ff)
            MMM{j} = [MMM{j}; ones(length(ff), 1) * i, J{i}(ff)'];
        end
    end
end

cnt = 1;
while cnt <= length(MMM)
    if isempty(MMM{cnt})
        MMM(cnt) = [];
    else
        cnt = cnt + 1;
    end
end
end

function MMM = injectiveAlignmentConstruction(C, J, L)

M = length(J);

[s, y]=sort(C(:, M+1), 'descend');

CC = cell(M, M);
for i = 1: M
    CC{i, i} = sparse(L(i), L(i));
end

for i = 1: M
    for j = i + 1: M
        f = find(C(:, i) > 0 & C(:, j) > 0);
        CC{i, j} = sparse(C(f, i), C(f, j), C(f, M+1), L(i), L(j));
        CC{j, i} = CC{i, j}';
    end
end

MM = cell(1, M);
for i = 1: M
    MM{i} = sparse(1, L(i));
end

k = 1;
g_max = 1;
EC = {};

while (k <= length(y))
    h = C(y(k), 1: M);
    if C(y(k), M + 1) < 0
        k = k + 1;
        continue;
    end
    if (prod(h) ~= 0)
        F = find(h ~= -1);
        f1 = F(1);
        f2 = F(2);
        h1 = h(f1);
        h2 = h(f2);
        c1 = MM{f1}(h1);
        c2 = MM{f2}(h2);
        
        if ((c1 == 0) && (c2 == 0))
            MM{f1}(h1) = g_max;
            MM{f2}(h2) = g_max;
            g_max = g_max + 1;
            EC = [EC, [f1, h1; f2, h2]];
            
        elseif ((c1 > 0) && (c2 == 0))
            if length(find(EC{c1}(:, 1) == f2)) >= 1
                k = k + 1;
                continue;
            end
            EC{c1} = [EC{c1}; f2, h2];
            MM{f2}(h2) = MM{f1}(h1);
            
        elseif ((c1 == 0) & (c2 > 0))
            if length(find(EC{c2}(:, 1) == f1)) >= 1
                k = k + 1;
                continue;
            end
            EC{c2} = [EC{c2}; f1, h1];
            MM{f1}(h1) = MM{f2}(h2);
        end
    end
    
    k = k + 1;
end

MMM = cell(g_max, 1);
for j = 1: g_max
    for i = 1: M
        ff = find(MM{i} == j);
        if ~isempty(ff)
            MMM{j} = [MMM{j}; ones(length(ff), 1) * i, J{i}(ff)'];
        end
    end
end

cnt = 1;
while cnt <= length(MMM)
    if isempty(MMM{cnt})
        MMM(cnt) = [];
    else
        cnt = cnt + 1;
    end
end
end

function [val, m1, m2, mi] = bipartite_matching(varargin)
% BIPARTITE_MATCHING Solve a maximum weight bipartite matching problem
%
% [val m1 m2]=bipartite_matching(A) for a rectangular matrix A
% [val m1 m2 mi]=bipartite_matching(x,ei,ej,n,m) for a matrix stored
% in triplet format.  This call also returns a matching indicator mi so
% that val = x'*mi.
%
% The maximum weight bipartite matching problem tries to pick out elements
% from A such that each row and column get only a single non-zero but the
% sum of all the chosen elements is as large as possible.
%
% This function is slightly atypical for a graph library, because it will
% be primarily used on rectangular inputs.  However, these rectangular
% inputs model bipartite graphs and we take advantage of that stucture in
% this code.  The underlying graph adjency matrix is
%   G = spaugment(A,0);
% where A is the rectangular input to the bipartite_matching function.
%
% Matlab already has the dmperm function that computes a maximum
% cardinality matching between the rows and the columns.  This function
% gives us the maximum weight matching instead.  For unweighted graphs, the
% two functions are equivalent.
%
% Note: If ei and ej contain duplicate edges, the results of this function
% are incorrect.
%
% See also DMPERM
%
% Example:
%   A = rand(10,8); % bipartite matching between random data
%   [val mi mj] = bipartite_matching(A);
%   val

% David F. Gleich and Ying Wang
% Copyright, Stanford University, 2008-2009
% Computational Approaches to Digital Stewardship

% 2008-04-24: Initial coding (copy from Ying Wang matching_sparse_mex.cpp)
% 2008-11-15: Added triplet input/output
% 2009-04-30: Modified for gaimc library
% 2009-05-15: Fixed error with empty inputs and triple added example.

[rp ci ai tripi n m] = bipartite_matching_setup(varargin{:});

if isempty(tripi)
    error(nargoutchk(0,3,nargout,'struct'));
else
    error(nargoutchk(0,4,nargout,'struct'));
end


if ~isempty(tripi) && nargout>3
    [val m1 m2 mi] = bipartite_matching_primal_dual(rp, ci, ai, tripi, n, m);
else
    [val m1 m2] = bipartite_matching_primal_dual(rp, ci, ai, tripi, n, m);
end
end

function [rp ci ai tripi n m]= bipartite_matching_setup(A,ei,ej,n,m)
% convert the input

if nargin == 1
    if isstruct(A)
        [nzi nzj nzv]=csr_to_sparse(A.rp,A.ci,A.ai);
    else
        [nzi nzj nzv]=find(A);
    end
    [n m]=size(A);
    triplet = 0;
elseif nargin >= 3 && nargin <= 5
    nzi = ei;
    nzj = ej;
    nzv = A;
    if ~exist('n','var') || isempty(n), n = max(nzi); end
    if ~exist('m','var') || isempty(m), m = max(nzj); end
    triplet = 1;
else
    error(nargchk(3,5,nargin,'struct'));
end
nedges = length(nzi);

rp = ones(n+1,1); % csr matrix with extra edges
ci = zeros(nedges+n,1);
ai = zeros(nedges+n,1);
if triplet, tripi = zeros(nedges+n,1); % triplet index
else tripi = [];
end

%
% 1. build csr representation with a set of extra edges from vertex i to
% vertex m+i
%
rp(1)=0;
for i=1:nedges
    rp(nzi(i)+1)=rp(nzi(i)+1)+1;
end
rp=cumsum(rp);
for i=1:nedges
    if triplet, tripi(rp(nzi(i))+1)=i; end % triplet index
    ai(rp(nzi(i))+1)=nzv(i);
    ci(rp(nzi(i))+1)=nzj(i);
    rp(nzi(i))=rp(nzi(i))+1;
end
for i=1:n % add the extra edges
    if triplet, tripi(rp(i)+1)=-1; end % triplet index
    ai(rp(i)+1)=0;
    ci(rp(i)+1)=m+i;
    rp(i)=rp(i)+1;
end
% restore the row pointer array
for i=n:-1:1
    rp(i+1)=rp(i);
end
rp(1)=0;
rp=rp+1;

%
% 1a. check for duplicates in the data
%
colind = false(m+n,1);
for i=1:n
    for rpi=rp(i):rp(i+1)-1
        if colind(ci(rpi)), error('bipartite_matching:duplicateEdge',...
                'duplicate edge detected (%i,%i)',i,ci(rpi));
        end
        colind(ci(rpi))=1;
    end
    for rpi=rp(i):rp(i+1)-1, colind(ci(rpi))=0; end % reset indicator
end

end

function [val m1 m2 mi]=bipartite_matching_primal_dual(...
    rp, ci, ai, tripi, n, m)
% BIPARTITE_MATCHING_PRIMAL_DUAL

alpha=zeros(n,1); % variables used for the primal-dual algorithm
beta=zeros(n+m,1);
queue=zeros(n,1);
t=zeros(n+m,1);
match1=zeros(n,1);
match2=zeros(n+m,1);
tmod = zeros(n+m,1);
ntmod=0;

%
% initialize the primal and dual variables
%
for i=1:n
    for rpi=rp(i):rp(i+1)-1
        if ai(rpi) > alpha(i), alpha(i)=ai(rpi); end
    end
end
% dual variables (beta) are initialized to 0 already
% match1 and match2 are both 0, which indicates no matches
i=1;
while i<=n
    % repeat the problem for n stages
    
    % clear t(j)
    for j=1:ntmod, t(tmod(j))=0; end
    ntmod=0;
    
    % add i to the stack
    head=1; tail=1;
    queue(head)=i; % add i to the head of the queue
    while head <= tail && match1(i)==0
        k=queue(head);
        for rpi=rp(k):rp(k+1)-1
            j = ci(rpi);
            if ai(rpi) < alpha(k)+beta(j) - 1e-8, continue; end % skip if tight
            if t(j)==0,
                tail=tail+1; queue(tail)=match2(j);
                t(j)=k;
                ntmod=ntmod+1; tmod(ntmod)=j;
                if match2(j)<1,
                    while j>0,
                        match2(j)=t(j);
                        k=t(j);
                        temp=match1(k);
                        match1(k)=j;
                        j=temp;
                    end
                    break; % we found an alternating path
                end
            end
        end
        head=head+1;
    end
    
    if match1(i) < 1 % still not matched, so update primal, dual and repeat
        theta=inf;
        for j=1:head-1
            t1=queue(j);
            for rpi=rp(t1):rp(t1+1)-1
                t2=ci(rpi);
                if t(t2) == 0 && alpha(t1) + beta(t2) - ai(rpi) < theta
                    theta = alpha(t1) + beta(t2) - ai(rpi);
                end
            end
        end
        
        for j=1:head-1, alpha(queue(j)) = alpha(queue(j)) - theta; end
        
        for j=1:ntmod, beta(tmod(j)) = beta(tmod(j)) + theta; end
        
        continue;
    end
    
    i=i+1; % increment i
end

val=0;
for i=1:n
    for rpi=rp(i):rp(i+1)-1
        if ci(rpi)==match1(i), val=val+ai(rpi); end
    end
end
noute = 0; % count number of output edges
for i=1:n
    if match1(i)<=m, noute=noute+1; end
end
m1=zeros(noute,1); m2=m1; % copy over the 0 array
noute=1;
for i=1:n
    if match1(i)<=m, m1(noute)=i; m2(noute)=match1(i);noute=noute+1; end
end

if nargout>3
    mi= false(length(tripi)-n,1);
    for i=1:n
        for rpi=rp(i):rp(i+1)-1
            if match1(i)<=m && ci(rpi)==match1(i), mi(tripi(rpi))=1; end
        end
    end
end

end

function [cs,index] = sort_nat(c,mode)
%sort_nat: Natural order sort of cell array of strings.
% usage:  [S,INDEX] = sort_nat(C)
%
% where,
%    C is a cell array (vector) of strings to be sorted.
%    S is C, sorted in natural order.
%    INDEX is the sort order such that S = C(INDEX);
%
% Natural order sorting sorts strings containing digits in a way such that
% the numerical value of the digits is taken into account.  It is
% especially useful for sorting file names containing index numbers with
% different numbers of digits.  Often, people will use leading zeros to get
% the right sort order, but with this function you don't have to do that.
% For example, if C = {'file1.txt','file2.txt','file10.txt'}, a normal sort
% will give you
%
%       {'file1.txt'  'file10.txt'  'file2.txt'}
%
% whereas, sort_nat will give you
%
%       {'file1.txt'  'file2.txt'  'file10.txt'}
%
% See also: sort

% Version: 1.4, 22 January 2011
% Author:  Douglas M. Schwarz
% Email:   dmschwarz=ieee*org, dmschwarz=urgrad*rochester*edu
% Real_email = regexprep(Email,{'=','*'},{'@','.'})

% Copyright (c) 2008, Douglas M. Schwarz
% All rights reserved.


% Set default value for mode if necessary.
if nargin < 2
    mode = 'ascend';
end

% Make sure mode is either 'ascend' or 'descend'.
modes = strcmpi(mode,{'ascend','descend'});
is_descend = modes(2);
if ~any(modes)
    error('sort_nat:sortDirection',...
        'sorting direction must be ''ascend'' or ''descend''.')
end

% Replace runs of digits with '0'.
c2 = regexprep(c, '\d+', '0');

% Compute char version of c2 and locations of zeros.
s1 = char(c2);
z = s1 == '0';

% Extract the runs of digits and their start and end indices.
[digruns,first,last] = regexp(c, '\d+', 'match', 'start', 'end');

% Create matrix of numerical values of runs of digits and a matrix of the
% number of digits in each run.
num_str = length(c);
max_len = size(s1,2);
num_val = NaN(num_str,max_len);
num_dig = NaN(num_str,max_len);
for i = 1:num_str
    num_val(i,z(i,:)) = sscanf(sprintf('%s ',digruns{i}{:}),'%f');
    num_dig(i,z(i,:)) = last{i} - first{i} + 1;
end

% Find columns that have at least one non-NaN.  Make sure activecols is a
% 1-by-n vector even if n = 0.
activecols = reshape(find(~all(isnan(num_val))),1,[]);
n = length(activecols);

% Compute which columns in the composite matrix get the numbers.
numcols = activecols + (1:2:2*n);

% Compute which columns in the composite matrix get the number of digits.
ndigcols = numcols + 1;

% Compute which columns in the composite matrix get chars.
charcols = true(1,max_len + 2*n);
charcols(numcols) = false;
charcols(ndigcols) = false;

% Create and fill composite matrix, comp.
comp = zeros(num_str,max_len + 2*n);
comp(:,charcols) = double(s1);
comp(:,numcols) = num_val(:,activecols);
comp(:,ndigcols) = num_dig(:,activecols);

% Sort rows of composite matrix and use index to sort c in ascending or
% descending order, depending on mode.
[unused,index] = sortrows(comp);
if is_descend
    index = index(end:-1:1);
end
index = reshape(index,size(c));
cs = c(index);
end

function [spectralSignatureMatrix]= calculateTopologySimilarity(G1, G2, T, r)
metric = 0; % 0: MWBM, 1: CIQ

G1 = G1{1};
G2 = G2{1};

M = size(G1, 1);
N = size(G2, 2);

% r = 1;

spectralSignatureMatrix = zeros(M, N);
[nonZeroRow, nonZeroColumn] = find(T > 0);

for index = 1: length(nonZeroRow)
    m = nonZeroRow(index);
    n = nonZeroColumn(index);
    % for m = 1: 1: M
    %     for n = 1: 1: N
    % Generate induced subgraph based on the length of shortest path
    neighbors = find(G1(m, :) ~= 0);
    inducedSubGraph1Indices = neighbors;
    for i = 2: 1: r
        inducedSubGraph1IndicesTemp = [];
        for j = 1: 1: length(neighbors)
            inducedSubGraph1IndicesTemp = [inducedSubGraph1IndicesTemp find(G1(neighbors(j), :) ~= 0)];
        end
        neighbors = unique(inducedSubGraph1IndicesTemp);
        inducedSubGraph1Indices = [inducedSubGraph1Indices neighbors];
    end
    inducedSubGraph1Indices = unique(sort(inducedSubGraph1Indices));
    inducedSubGraph1Indices = inducedSubGraph1Indices(inducedSubGraph1Indices ~= m);
    
    % Generate induced subgraph based on the length of shortest path
    %         neighbors = find(G2(nonZeroColumn(index), :) ~= 0);
    neighbors = find(G2(n, :) ~= 0);
    inducedSubGraph2Indices = neighbors;
    for i = 2: 1: r
        inducedSubGraph2IndicesTemp = [];
        for j = 1: 1: length(neighbors)
            inducedSubGraph2IndicesTemp = [inducedSubGraph2IndicesTemp find(G2(neighbors(j), :) ~= 0)];
        end
        neighbors = unique(inducedSubGraph2IndicesTemp);
        inducedSubGraph2Indices = [inducedSubGraph2Indices neighbors];
    end
    %         inducedSubGraph2Indices = unique(sort([inducedSubGraph2Indices nonZeroColumn(index)]));
    inducedSubGraph2Indices = unique(sort(inducedSubGraph2Indices));
    inducedSubGraph2Indices = inducedSubGraph2Indices(inducedSubGraph2Indices ~= n);
    
    inducedSimilarityMatrix = T(inducedSubGraph1Indices, inducedSubGraph2Indices);
    
    if sum(abs(inducedSimilarityMatrix(:))) ~= 0
        [a, b, c] = bipartite_matching(inducedSimilarityMatrix);
 
        if metric == 1
            % For CIQ
            G = {G1, G2};
            
            Amn = [[m; inducedSubGraph1Indices(b)'] [n;inducedSubGraph2Indices(c)']];
            Nmn = length(Amn);
            sumOfE_Ci_Cj = 0;
            sumOfNumeratorOfCIQ = 0;
            for i = 1: Nmn
                for j = i + 1: Nmn
                    E_Ci_Cj = 0;
                    num = 0;
                    den = 0;
                    
                    for k = 1: 2
                        % Calculate E_Ci_Cj
                        E_Ci_Cj = E_Ci_Cj + length(find(G{k}(Amn(i, k), Amn(j, k)) > 0));
                        
                        if and(~isempty(Amn(i, k)), ~isempty(Amn(j, k)))
                            den = den + 1;
                            if ~isempty(find(G{k}(Amn(i, k), Amn(j, k)) > 0))
                                num = num + 1;
                            end
                        end
                    end
                    
                    if (num == 1)
                        cs = 0;
                    else
                        if (num == 0) && (den == 0)
                            cs = 0;
                        else
                            cs = num/den;
                        end
                    end
                    
                    sumOfNumeratorOfCIQ = sumOfNumeratorOfCIQ + (cs * E_Ci_Cj);
                    sumOfE_Ci_Cj = sumOfE_Ci_Cj + E_Ci_Cj;
                end
            end
            if (sumOfNumeratorOfCIQ == 0) && (sumOfE_Ci_Cj == 0)
                CIQ = 0;
            else
                CIQ = sumOfNumeratorOfCIQ/sumOfE_Ci_Cj;
            end
            
            a = CIQ;
            % End CIQ
        end
    else
        a = 0;
    end
    
    spectralSignatureMatrix(m, n) = a;
end
% end

spectralSignatureMatrix = sparse(spectralSignatureMatrix);

end