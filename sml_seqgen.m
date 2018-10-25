function varargout=sml_seqgen(what,varargin);

baseDir = '/Users/eberlot/Documents/Data/SuperMotorLearning/stimuli_struct';
cd (baseDir);
rng('shuffle');

%%%%%%%%%%%%%% SeqGen %%%%%%%%%%%%%%
% Generates sets of 6 sequences with restrictions
% starting with 1/3/5
% max of 9 chunks
% no finger repetitions, trills
% SeqGen.mat - produced structure


%%%%%%%%%%%%% SeqEval %%%%%%%%%%%%%%%%
% Evaluates produced sequence sets
% Creates structure with all correlations between feature models 
% Feature models: first fing, all fing, transitions, chunks, seq
% AllCorrelations.mat - produced structure


%%%%%%%%%%%%% EvalCorr %%%%%%%%%%%%%%%
% Plots correlation values - mean corr, max corr
% Determines best sequence sets
% Saves the chosen sequence set with correlations, distances - SeqFinal.mat


switch (what)
    case 'SeqGen'            % Generate all possible seq with constrictions
        
        persistent TRIP TRILL FingRep TripDif 
        
        NumChunkSet = 100000;    % how many chunk sets to use in generating sequences
        found_ind = 0;
        fail=0;
        
        %all possible triplets
        TRIP = [kron([1:5]',ones(25,1)) kron(kron(ones(5,1),[1:5]'),ones(5,1)) kron(ones(25,1),[1:5]')];
        TripDif = diff(TRIP,[],2);% difference of triplets
        AbsTripDif = abs(TripDif);
        TRILL = AbsTripDif(:,1) == AbsTripDif(:,2);
        
        FingRep = (TripDif(:,1)==0 | TripDif(:,2)==0);
        crit = TRILL + FingRep;     % criterion - no trills, no fing rep;
        
        AllChunksOrd = TRIP(crit==0,:);
        AllChunks = AllChunksOrd(randperm(size(AllChunksOrd,1)),:);     % shuffle the order of chunks around
        numTrip = size(AllChunks,1);
        
        % Choose appropriate chunk sets     
        for seq = 1:NumChunkSet
            clear chunkchos;
            
            % first choose 6 chunks so that 2 start with 1, 3, 5
            chunkchos(1:2,:) = datasample(AllChunks(AllChunks(:,1)==1,:),2,'Replace',false);
            chunkchos(3:4,:) = datasample(AllChunks(AllChunks(:,1)==3,:),2,'Replace',false);
            chunkchos(5:6,:) = datasample(AllChunks(AllChunks(:,1)==5,:),2,'Replace',false);
            
            good1 = false;
            
            % choose additional chunks (different from above) 
            % a total of 9 chunks
            while good1 == false
                chunkchos(7:9,:) = datasample(AllChunks(AllChunks(:,1)==1,:),3);
                if (sum(ismember(chunkchos(1:6,:),chunkchos(7,:),'rows'))==1 |...
                    sum(ismember(chunkchos(1:6,:),chunkchos(8,:),'rows'))==1 | ...
                    sum(ismember(chunkchos(1:6,:),chunkchos(9,:),'rows'))==1)
                else
                    good1 = 1;
                end
            end 
            
            ChunkChos{seq}=chunkchos;
            chunknum = 9;
        end
        
        for seq = 1:NumChunkSet
            
            % determine the first chunk for the 6 sequences
            Seq = zeros(6,9);
            Seq(:,1:3)=ChunkChos{seq}(1:6,:);
            
            
            % add 2 chunks to each sequence (out of set of 9)
            % no chunk repetition, no finger repetition 
            
            % Seq 1
            good2=false;
            while  good2 == false % no two chunks repeated, no finger rep
                s1 = (datasample(ChunkChos{seq}(2:end,:),2,'Replace',false))';    
                Seq(1,4:9) = s1(:);
                if (sum(ismember(diff(Seq(1,:)),0))==0)
                    good2=1;
                end
            end
            
            % Seq 2
            good2=false;
            while  good2 == false % no two chunks repeated, no finger rep
                s2 = (datasample(ChunkChos{seq}([1 3:end],:),2,'Replace',false))';    
                Seq(2,4:9) = s2(:);
                if (sum(ismember(diff(Seq(2,:)),0))==0)
                    good2=1;
                end
            end
            
            % Seq 3
            good2=false;
            while  good2 == false % no two chunks repeated, no finger rep
                s3 = (datasample(ChunkChos{seq}([1:2 4:end],:),2,'Replace',false))';    
                Seq(3,4:9) = s3(:);
                if (sum(ismember(diff(Seq(3,:)),0))==0)
                    good2=1;
                end
            end
            
            % Seq 4
            good2=false;
            while  good2 == false % no two chunks repeated, no finger rep
               s4 = (datasample(ChunkChos{seq}([1:3 5:end],:),2,'Replace',false))';
               Seq(4,4:9) = s4(:);
                if ((sum(ismember(diff(Seq(4,:)),0))==0))
                        good2=1;
                end
            end
            
            % Seq 5
            good2=false;
            while  good2 == false % no two chunks repeated, no finger rep
               s5 = (datasample(ChunkChos{seq}([1:4 6:end],:),2,'Replace',false))';
               Seq(5,4:9) = s5(:);
                if ((sum(ismember(diff(Seq(5,:)),0))==0))
                        good2=1;
                end
            end
            
            % Seq 6
            good2=false;
            while  good2 == false % no two chunks repeated, no finger rep
                s6 = (datasample(ChunkChos{seq}([1:5 7:end],:),2,'Replace',false))';
                Seq(6,4:9) = s6(:);
                if ((sum(ismember(diff(Seq(6,:)),0))==0))
                        good2=1;
                end
            end
            
            % Repetitions of fingers per sequence
            for seq=1:6
                for fing = 1:5
                    RepFing(seq,fing)=sum(ismember(Seq(seq,:),fing));
                end
            end
            
            % Check that all sequences contain all fingers
            i=1:5;
            Fing(i)=sum(ismember(RepFing(:,i),0));
            
            if sum(Fing)>0
                fail=fail+1;
                %fprintf('Failed because of finger rep %d\n',fail);
                %continue
            else
                found_ind = found_ind+1;
                SeqChos{found_ind} = Seq;
                fprintf('Found %d sequence set\n',found_ind);
            end
        end
        
        save(['SeqGen.mat'],'SeqChos');
        save(['AllTrip.mat'],'AllChunksOrd');
   
    case 'SeqEval'           % Eval sequences based on criteria
        
        figure=varargin{1}; % whether to plot figures
        load('SeqGen.mat');  % sequences
        A = load('sequenceType.mat');   % all possible doubles, triplets
        load('AllTrip.mat');      % sequences
        
        trip=AllChunksOrd;  % to keep the naming

        seqnum = size(SeqChos,2);
        
        for s = 1:1000
        %for s = 1:seqnum
            F.FirstFing = zeros(6,5);
            F.AllFing = zeros(6,5);
            F.TransFing = zeros(6,25);
            F.Chunks =  zeros(6,52);
            
            % First finger, All fingers
            for u = 1:size(SeqChos{s},1)    % set of 6 seq
                firstfing = SeqChos{s}(:,1);
                F.FirstFing(u,firstfing(u)) = 1;
                for j = 1:5                 % 5 fingers
                    placenumb = find(SeqChos{s}(u,:) == j);
                    F.AllFing(u,j) = length(placenumb);
                end
            end
         
            % Transitions
            Rd.Doubles = A.Sequence(126:150,1:2);
            
            for k = 1:size(SeqChos{s},1)    % set of 6 seq
                for p = 1:size(Rd.Doubles,1)
                    if any(ismember(SeqChos{s}(k,1:8),Rd.Doubles(p,1)))
                        ind = find(ismember(SeqChos{s}(k,1:8),Rd.Doubles(p,1)));    % find matching digits to first finger in doublet
                        
                        if any(ismember(SeqChos{s}(k,ind+1),Rd.Doubles(p,2)))       % compare the second finger in doublet
                            F.TransFing(k,p) = F.TransFing(k,p) + sum(ismember(SeqChos{s}(k,ind+1),Rd.Doubles(p,2)));
                        end
                    end
                end
                
            end
            
            % Chunks         
            clear p;
            for h = 1:size(SeqChos{s},1)
                for p = 1:size(trip,1) % all allowed triplets
                    if strcmp(num2str(SeqChos{s}(h,1:3)),num2str(trip(p,1:3))) | strcmp(num2str(SeqChos{s}(h,4:6)),num2str(trip(p,1:3))) | strcmp(num2str(SeqChos{s}(h,7:9)),num2str(trip(p,1:3)))
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
            
            D.ChunkAmount = length(F.uniqueChunks);
            WhichChunks = F.whichChunks;
            D.ChunkHowOft = HowoftenChunk;
            
            % build features models
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
            Corr_var=[D.distFirst' D.distAll' D.distTrans' D.distChunk' D.distSeq'];
            AllCorr  = corrN(Corr_var);    % AllCorrelations of interest
            
            D.AllCorr=AllCorr;
            D.MaxCorr = max(max(triu(AllCorr,1))); 
            D.MeanCorr = mean(AllCorr(triu(AllCorr,1)>0));
            
            D=addstruct(D,G);
            
            Corr{s}=D;
            
        switch(figure)
            case 1;
                figure(s)
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
            
        end
        
        save('AllCorrelations.mat','Corr');     % save structure with all correlations
        
    case 'Eval_Corr'         % Plot all correlations, determine set of sequences with lowest corr           
         load('AllCorrelations.mat');   % Corr
         load('SeqGen.mat');                % SeqChos - sequences
         figure=varargin{1};
         SeqField={};
         
         for c = 1:size(Corr,2)
             MeanCorr(c) = Corr{c}.MeanCorr;
             MaxCorr(c) = Corr{c}.MaxCorr;
         end
         
         switch(figure)
             case 1
                 figure(1)
                 subplot(1,2,1);
                 histogram(MeanCorr);
                 title('Mean Correlation');
                 subplot(1,2,2);
                 histogram(MaxCorr);
                 title('Max Correlation');
         end
         
         
         % find sequence set with lowest correlations
         S_BestMax = find(ismember(MaxCorr,min(MaxCorr)));      % sequence with the lowest max correlation
         S_BestMean = find(ismember(MeanCorr,min(MeanCorr)));   % sequence with the lowest mean correlation
         
         % Create structure S containing info about chosen sequence
         S.SeqFinal = SeqChos{S_BestMax};
         S.AllCorr = Corr{S_BestMax}.AllCorr;
         S.MaxCorr = Corr{S_BestMax}.MaxCorr;
         S.MeanCorr = Corr{S_BestMax}.MaxCorr;
         S.ChunkAmount = Corr{S_BestMax}.ChunkAmount;
         S.ChunkHowOft = Corr{S_BestMax}.ChunkHowOft;
         S.distFirst = Corr{S_BestMax}.distFirst;
         S.distAll = Corr{S_BestMax}.distAll;
         S.distTrans = Corr{S_BestMax}.distTrans;
         S.distChunk = Corr{S_BestMax}.distChunk;
         S.distSeq = Corr{S_BestMax}.distSeq;
         
         
         save('SeqFinal.mat','SeqFinal','-struct','S');   % save the final sequence with its correlations, distances etc.
         
         keyboard;
                
    case 'getChunks'         % Separate sequences into chunks - called from SeqEval
        Seq = varargin{1};
        for n=1:size(Seq,1);
            t(n*3-2:n*3,:) = [Seq(n,[1:3]);Seq(n,[4:6]);Seq(n,[7:9])];
            
        end
        varargout = {t};
        
    otherwise
        error('no such case!')
        
end