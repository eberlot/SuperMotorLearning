function  varargout = sml1_sequencestimuli(what,varargin)

DirData = '/Users/eberlot/Documents/Data/SuperMotorLearning';


% Stimuli
D.sequences = [ 3 2 1 4 5 3 1 2 5; 3 4 5 2 1 3 5 4 1; 1 4 2 4 3 1 4 5 3; 5 2 4 2 3 5 2 1 3; 1 2 5 2 3 5 3 4 5; 5 4 1 4 3 1 3 2 1];
D.chunk = [ 3 2 1; 3 4 5; 1 4 2 ;5 2 4; 1 2 5; 5 4 1; 2 1 3; 4 5 3; 4 3 1 ; 2 3 5];
D.trans = [ 3 2; 2 1; 1 4; 4 1; 1 3; 3 1; 2 5; 3 4; 4 5; 5 2; 5 3; 3 5; 5 4; 4 2; 2 4; 4 3; 2 3; 1 2];
D.transchunk = [ 1 4; 3 1; 5 2; 3 5; 2 4; 4 2; 5 3; 1 3 ];
% Practice stimuli (4 seq) - TO DO!!!!!
% scramble chosen sequences
% compare with actual stimuli, between themselves

% Feature models
F.FirstFing = zeros(6,5);   
F.AllFing = zeros(6,5);
F.TransFing = zeros(6,25);
F.Chunks =  zeros(6,10);
F.TransChunk = zeros(6,25);

switch(what)
    
    case 'featureModels'            % Create feature models - FirstFing, AllFing, TransFing, ChunkFing (add for practice seq)
        
        A = load(fullfile(DirData,'sequenceType.mat'));   
        % First finger, All fingers
        for i = 1:size(D.sequences,1)
            firstfing = D.sequences(:,1);
            F.FirstFing(i,firstfing(i)) = 1;
            for j = 1:5
                placenumb = find(D.sequences(i,:) == j);
                F.AllFing(i,j) = length(placenumb);
            end
        end
        
        
        % Transitions
        Rd.Doubles = A.Sequence(126:150,1:2);
        
        for i = 1:size(D.sequences,1)
            for p = 1:size(Rd.Doubles,1)
                if any(ismember(D.sequences(i,1:8),Rd.Doubles(p,1)))
                    ind = find(ismember(D.sequences(i,1:8),Rd.Doubles(p,1)));
                    
                    if any(ismember(D.sequences(i,ind+1),Rd.Doubles(p,2)))
                        F.TransFing(i,p) = F.TransFing(i,p) + sum(ismember(D.sequences(i,ind+1),Rd.Doubles(p,2)));
                    end
                end
            end
            
        end
        
        % Chunks
        
        for i = 1:size(D.sequences,1)
            for p = 1:size(D.chunk,1)
                if strcmp(num2str(D.sequences(i,1:3)),num2str(D.chunk(p,1:3))) | strcmp(num2str(D.sequences(i,4:6)),num2str(D.chunk(p,1:3))) | strcmp(num2str(D.sequences(i,7:9)),num2str(D.chunk(p,1:3)))
                    F.Chunks(i,p) = 1;
                end
            end
        end
        
        % Transitions between chunks
         
        for i = 1:size(D.sequences,1)
            for p = 1:size(Rd.Doubles,1)
                if strcmp(num2str(D.sequences(i,3:4)),num2str(Rd.Doubles(p,:))) | strcmp(num2str(D.sequences(i,6:7)),num2str(Rd.Doubles(p,:))) 
                    F.TransChunk(i,p) = 1;
                end
            end
        end
        
        
        
        % Save stimuli - D
        save(fullfile(DirData,'Stimuli.mat'),'-struct','D');
        
        % Save feature models
        save(fullfile(DirData,'FeatureModels.mat'),'-struct','F');
        
        keyboard;
    case 'distModels'               % Calculate distances between sequences per feature model, correlations (add for practice seq)
        
        F = load(fullfile(DirData,'FeatureModels.mat'));
        
        Dist = [];
        
        % G matrices
        G_First       = F.FirstFing * F.FirstFing';
        G_All         = F.AllFing * F.AllFing';
        G_TransFing   = F.TransFing * F.TransFing';
        G_Chunks      = F.Chunks * F.Chunks';
        G_TransChunk  = F.TransChunk * F.TransChunk';
        
        S_ind         = [indicatorMatrix('allpairs',[1:6])];      % indicator matrix for sequences
        
        % distances
        distFirst     = diag(S_ind*G_First*S_ind')';
        distAll       = diag(S_ind*G_All*S_ind')';
        distTransFing = diag(S_ind*G_TransFing*S_ind')';
        distChunk     = diag(S_ind*G_Chunks*S_ind')';
        distTransChunk= diag(S_ind*G_TransChunk*S_ind')';

        
        % save G, distances to structure
        D.G_First       = G_First(:);
        D.G_All         = G_All(:);
        D.G_TransFing   = G_TransFing(:);
        D.G_Chunks      = G_Chunks(:);
        D.G_TransChunk  = G_TransChunk(:);
        D.S_ind         = S_ind;
        D.distFirst     = distFirst;
        D.distAll       = distAll;
        D.distTransF    = distTransFing;
        D.distTransC    = distTransChunk;
        D.distChunk     = distChunk;
        
        
        % correlate the models (All/First/Trans/Chunk)
        c.FirstAll      = corr(distFirst', distAll');
        c.FirstTrans    = corr(distFirst', distTransFing');
        c.FirstChunk    = corr(distFirst', distChunk');
        c.AllTransFing  = corr(distAll', distTransFing');
        c.AllChunk      = corr(distAll', distChunk');
        c.TransFChunk   = corr(distTransFing', distChunk');     % very problematic!
        c.TransCChunk   = corr(distTransChunk', distChunk');    % problematic!
        
        Dist = addstruct(Dist,D);
        Dist = addstruct(Dist,c);
        
        save(fullfile(DirData,'Distances.mat'),'-struct','Dist');
    case 'plot_featureModels'       % Plot feature models
        
        FM = load(fullfile(DirData,'FeatureModels.mat'));
        figure(1)
        subplot(2,2,1)
        set(gcf,'Color',[1 1 1])
        
        imagesc(FM.FirstFing)
        colorbar
        caxis([0, 1])
        title('FirstFing')
        
        subplot(2,2,2)
        imagesc(FM.AllFing)
        title('AllFing')
        colorbar
        caxis([0, 1])
        
        subplot(2,2,3)
        imagesc(FM.TransFing)
        title('TransFing')
        colorbar
        caxis([0, 1])
        colorbar
        
        subplot(2,2,4)
        imagesc(FM.ChunksFing)
        title('ChunkFing')
        colorbar
        caxis([0, 1])
        
    otherwise
        error('no such case!')
        
end

