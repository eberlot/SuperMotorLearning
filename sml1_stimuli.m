function varargout=sml1_stimuli(what,varargin);

baseDir = '/Users/eberlot/Documents/Data/SuperMotorLearning/stimuli_struct';
% baseDir = '/Volumes/MotorControl/data/SuperMotorLearning/stimuli_struct';
cd(baseDir)
A = load('sequenceType.mat');   % all possible doubles, triplets
load('SeqAll.mat','trip');      % sequences


F.FirstFing = zeros(6,5);
F.AllFing = zeros(6,5);
F.TransFing = zeros(6,25);
F.Chunks =  zeros(6,52);

switch(what)
    case 'featureModel'
    Seq=varargin{1};
    figAllow=varargin{2};
    clear p;
        switch(Seq)
            case 1    
                SeqChos = [1 5 3 5 3 4 3 1 5; ...
                1 5 2 1 4 5 3 4 2; ...
                3 4 2 1 4 5 1 5 2; ...
                3 1 5 3 4 2 1 4 5; ...
                5 3 4 5 4 2 1 5 3; ...
                5 4 2 3 4 2 3 1 5];
            case 2
                SeqChos = [5 1 3 1 3 2 3 5 1; ...
                5 1 4 5 2 1 3 2 4; ...
                3 2 4 5 2 1 5 1 4; ...
                3 5 1 3 2 4 5 2 1; ...
                1 3 2 1 2 4 5 1 3; ...
                1 2 4 3 2 4 3 5 1];
            case 3    
                SeqChos = [1 4 3 5 3 4 3 4 2; ...
                1 5 2 1 4 5 3 4 2; ...
                3 4 2 1 4 5 1 5 2; ...
                3 1 5 3 4 2 1 4 5; ...
                5 3 4 5 4 2 1 4 3; ...
                5 4 2 3 4 2 3 1 5];
            case 4
                SeqChos = [5 2 3 1 3 2 3 2 4; ...
                5 1 4 5 2 1 3 2 4; ...
                3 2 4 5 2 1 5 1 4; ...
                3 5 1 3 2 4 5 2 1; ...
                1 3 2 1 2 4 5 2 3; ...
                1 2 4 3 2 4 3 5 1];
        end
        
        % First finger, All fingers
        for u = 1:size(SeqChos,1)    % set of 6 seq
        firstfing = SeqChos(:,1);
        F.FirstFing(u,firstfing(u)) = 1;
        for j = 1:5                 % 5 fingers
            placenumb = find(SeqChos(u,:) == j);
            F.AllFing(u,j) = length(placenumb);
        end
        end


        % Transitions
        Rd.Doubles = A.Sequence(126:150,1:2);

        for k = 1:size(SeqChos,1)    % set of 6 seq
            for p = 1:size(Rd.Doubles,1)
                if any(ismember(SeqChos(k,1:8),Rd.Doubles(p,1)))
                    ind = find(ismember(SeqChos(k,1:8),Rd.Doubles(p,1)));    % find matching digits to first finger in doublet
                    
                    if any(ismember(SeqChos(k,ind+1),Rd.Doubles(p,2)))       % compare the second finger in doublet
                        F.TransFing(k,p) = F.TransFing(k,p) + sum(ismember(SeqChos(k,ind+1),Rd.Doubles(p,2)));
                    end
                end
            end
            
        end
        
       % Chunks

       clear p;
         for h = 1:size(SeqChos,1)
            for p = 1:size(trip,1) % all allowed triplets
                if strcmp(num2str(SeqChos(h,1:3)),num2str(trip(p,1:3))) | strcmp(num2str(SeqChos(h,4:6)),num2str(trip(p,1:3))) | strcmp(num2str(SeqChos(h,7:9)),num2str(trip(p,1:3)))
                    F.Chunks(h,p) = 1;
                end
            end
         end
        
        for h = 1:size(SeqChos,1)
            F.whichChunks(h,:) = find(F.Chunks(h,:)==1);
        end
        F.uniqueChunks(1,:) = unique(F.whichChunks);
        
        for ch = 1:length(F.uniqueChunks)
            HowoftenChunk(ch) = sum(F.whichChunks(:) == F.uniqueChunks(ch)); % all allowed triplets
        end

        ChunkAmount = length(F.uniqueChunks);
        WhichChunks = F.whichChunks;
        ChunkHowOft = HowoftenChunk;

        % G matrices
        G.First       = F.FirstFing * F.FirstFing';
        G.All         = F.AllFing * F.AllFing';
        G.TransFing   = F.TransFing * F.TransFing';
        G.Chunks      = F.Chunks * F.Chunks';

        S_ind         = [indicatorMatrix('allpairs',[1:6])];      % indicator matrix for sequences

        % distances
        distFirst     = diag(S_ind*G.First*S_ind')';
        distAll       = diag(S_ind*G.All*S_ind')';
        distTransFing = diag(S_ind*G.TransFing*S_ind')';
        distChunk     = diag(S_ind*G.Chunks*S_ind')';
        
        % save G, distances to structure
        D.G_First       = G.First(:);
        D.G_All         = G.All(:);
        D.G_TransFing   = G.TransFing(:);
        D.G_Chunks      = G.Chunks(:);
        D.distFirst     = distFirst;
        D.distAll       = distAll;
        D.distTrans     = distTransFing;
        D.distChunk     = distChunk;  
        D.distSeq       = ones(1,15); 
        % correlate the models (All/First/Trans/Chunk)
        % correlations with chunk are 0 if chunks in all chosen
        % sequences are different
        A=[D.distFirst' D.distAll' D.distTrans' D.distChunk' D.distSeq']; 
        X  = corrN(A); 
        
        
        switch(figAllow)
            case 1;
                figure(4)
                subplot(2,2,1)
                imagesc(rsa.rdm.squareRDM(D.distFirst))
                title('FirstFing')
                subplot(2,2,2)
                imagesc(rsa.rdm.squareRDM(D.distAll))
                title('AllFing')
                subplot(2,2,3)
                imagesc(rsa.rdm.squareRDM(D.distTrans))
                title('TransFing')
                subplot(2,2,4)
                imagesc(rsa.rdm.squareRDM(D.distChunk))
                title('Chunks')
        end
        
        keyboard;
        % criterion - max(max(triu(X,1))); mean(X(triu(X,1)>0))
    otherwise
        error('no such case!')
end