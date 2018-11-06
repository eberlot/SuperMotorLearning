function colors = colorDendrogram(Z, N, varargin);
%% colors = colorDendrogram(Z, N, varargin);
% Gives sensible colors given linkage matrix Z (=linkage(X)).
% 
% 
% 
% 
% 
% 

colorspace = 'rgb';
fig = 0;
rotation = [];
order = [1,2,3];
weight = 1;
vararginoptions(varargin,{'colorspace','fig','rotation','order','weight'});

% set rand seed
randSeed = 1234;
method = 'twister';
rng(randSeed, method); % init rng with fixed seed for replicability (kmeans inside spectralclustering uses rng)
            

% check input size
Nrow = size(Z,1);

% get membership
pairs = Z(:,1:2);
clusters = zeros(Nrow,3);
for i=1:Nrow
    clusters(i,1:3) = [N+i, Z(i,1), Z(i,2)];
end

% take 4 most separating clusters
bases = {[1 0 0], [0 1 0],[0 0 1]};

% 1
i = Nrow;
fprintf('C%d=C%d+C%d\n',clusters(i,1),clusters(i,2),clusters(i,3));
Curr.cluster = cat(1,clusters(i,2),clusters(i,3));
Curr.pos = [0 0 0; bases{1+mod(i,3)}*Z(i,3)];
Curr.base = [bases{1+mod(i,3)};bases{1+mod(i,3)}];

% 2
i = Nrow-1;
fprintf('C%d=C%d+C%d\n',clusters(i,1),clusters(i,2),clusters(i,3));
Curr.cluster = cat(1,Curr.cluster,clusters(i,2),clusters(i,3));
idx=Curr.cluster==clusters(i,1);
p0=Curr.pos(idx,:);
p1 = [p0+0.5*Z(i,3)*bases{1+mod(i,3)}; p0-0.5*Z(i,3)*bases{1+mod(i,3)}];
Curr.pos = [Curr.pos; p1];
Curr.base = [Curr.base; bases{1+mod(i,3)}; bases{1+mod(i,3)}];

Curr.cluster=Curr.cluster(~idx);
Curr.pos=Curr.pos(~idx,:);
Curr.base = Curr.base(~idx,:);

% 3
i = Nrow-2;
fprintf('C%d=C%d+C%d\n',clusters(i,1),clusters(i,2),clusters(i,3));
Curr.cluster = cat(1,Curr.cluster,clusters(i,2),clusters(i,3));
idx=Curr.cluster==clusters(i,1);
p0=Curr.pos(idx,:);
p1 = [p0+0.5*Z(i,3)*bases{1+mod(i,3)}; p0-0.5*Z(i,3)*bases{1+mod(i,3)}];
Curr.pos = [Curr.pos; p1];
Curr.base = [Curr.base; bases{1+mod(i,3)}; bases{1+mod(i,3)}];

Curr.cluster=Curr.cluster(~idx);
Curr.pos=Curr.pos(~idx,:);
Curr.base = Curr.base(~idx,:);

% Deal with the rest of clusters
for j=3:Nrow-1
    i = Nrow-j;
    fprintf('C%d=C%d+C%d\n',clusters(i,1),clusters(i,2),clusters(i,3));
    idx=find(Curr.cluster==clusters(i,1));
    
    Curr.cluster = cat(1,Curr.cluster,clusters(i,2),clusters(i,3));
    p0=Curr.pos(idx,:);
    b0=Curr.base(idx,:); 
    if isequal(b0,bases{1+mod(i,3)});
        remain=setdiff([1:3],1+mod(i,3));
        b0=bases{remain(randperm(2))};
    else
        b0=bases{1+mod(i,3)};
    end        
    p1 = [p0+weight/2*Z(i,3)*b0; p0-weight/2*Z(i,3)*b0];
    
    Curr.pos = [Curr.pos; p1];
    Curr.base = [Curr.base; b0; b0];
    
    Ncluster=length(Curr.cluster);
    
    Curr.cluster = Curr.cluster(setdiff([1:Ncluster] , idx));
    Curr.pos = Curr.pos(setdiff([1:Ncluster],idx) , :);    
    Curr.base = Curr.base(setdiff([1:Ncluster],idx) , :);
end
% Sort
[~,idx]=sort(Curr.cluster);

% Normalize to 0-1
U=Curr.pos(idx,order);
U=bsxfun(@minus,U,min(U,[],1));
U=bsxfun(@rdivide,U,max(U,[],1));

% Make color vector
switch lower(colorspace)
    case 'hsv'        
        % prevent full saturation and blackout
        U(:,2) = U(:,2) + 0.3; U(:,2) = min([U(:,2),ones(size(U(:,2)))],[],2);
        U(:,3) = U(:,3) + 0.3; U(:,3) = min([U(:,3),ones(size(U(:,3)))],[],2);
        colors=hsv2rgb(U);
    case 'rgb'
        colors=U;    
    otherwise
        error('not implemented')
end

if fig==1
    %myFigure([35,20]*0.7);
    figure
    subplot(1,2,1);
    for r=0:0.5:1
        for g=0:0.5:1
            for b=0:0.5:1
                plot3(r,g,b,'o',...
            'markerfacecolor',[r,g,b],...
            'markeredgecolor',[0 0 0],...
            'markersize',30); hold on; grid on;
            %alpha(0.5);
            end
        end
    end; axis square
    xlabel('e1');ylabel('e2');zlabel('e3')
    
    subplot(1,2,2);
    for i=1:size(U,1)
        plot3(U(i,1),U(i,2),U(i,3),'o',...
            'markerfacecolor',colors(i,:),...
            'markeredgecolor',colors(i,:),...
            'markersize',10); hold on; grid on;
        %alpha(1);
        %arrow3([0,0,0],[U(i,1),U(i,2),U(i,3)]);
        text(U(i,1),U(i,2),U(i,3),sprintf('C%d',i));
    end; axis square
    plot3(0,0,0,'*',...
            'markerfacecolor','k',...
            'markeredgecolor','k',...
            'markersize',10);
    xlabel('e1');ylabel('e2');zlabel('e3')
end

end

function member = getmember(N,clusters,num)
    member = clusters(clusters(:,1)==num,2:3);    
    iscombo = member>N;
    combos=member(iscombo);
    if ~isempty(combos)
        member(iscombo) = [];    
        for c=1:numel(combos)
            member = [member, getmember(N,clusters,combos(c))];
        end
    end
end