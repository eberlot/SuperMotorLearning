function varargout=sml_findSeq_run(what,varargin);

baseDir = '/Users/eberlot/Documents/Data/SuperMotorLearning/stimuli_struct';
%baseDir = 'C:\Nicola\SuperMotorLearning\data';
cd (baseDir);

switch (what)
    case 'MakeSequences'            % Make all possible seq; varargin - allow repetition of chunk (0), finger (0), trill (0)
        
        persistent TRIP Tisrun TrepTyp TripDif TRILL
        
        % Same numbering as in sh3
        
        %all possible triplets
        TRIP = [kron([1:5]',ones(25,1)) kron(kron(ones(5,1),[1:5]'),ones(5,1)) kron(ones(25,1),[1:5]')];
        TripDif =diff(TRIP,[],2);% difference of triplets
        
        AbsTripDif = abs(TripDif);
        TRILL = AbsTripDif(:,1) == AbsTripDif(:,2);
        Tisrun = TripDif(:,1).*TripDif(:,2)==1; %is it a run
        TrepTyp = double(TripDif(:,1)==0)+2*double(TripDif(:,2)==0); % if repetition between first two letters than 1
        %if repetition second pair then 2, if all three same number then 3
        
        
        %TRIP
        allowRepChunk = varargin{1}; %repetition of chunks in same sequence
        allowRepFinger = varargin{2};
        allowTrill = varargin{3};
        if allowTrill == 0
            TRIP = TRIP(TRILL ==0 | Tisrun == 1 ,:);
        end
        numTrip = size(TRIP,1);
        trip = TRIP;
        % - ADD that no repetition within chunk or across chunk
        if (allowRepChunk)
            COMB = [kron([1:numTrip]',ones(numTrip*numTrip,1)) kron(kron(ones(numTrip,1),[1:numTrip]'),ones(numTrip,1)) kron(ones(numTrip*numTrip,1),[1:numTrip]')];
        else
            COMB = [kron([1:numTrip]',ones(numTrip*numTrip,1)) kron(kron(ones(numTrip,1),[1:numTrip]'),ones(numTrip,1)) kron(ones(numTrip*numTrip,1),[1:numTrip]')];
            rep=(COMB(:,1)==COMB(:,2) | COMB(:,1)==COMB(:,3) | COMB(:,2)==COMB(:,3));
            COMB=COMB(~rep,:);
            
        end;
        D.Seq = [trip(COMB(:,1),:) trip(COMB(:,2),:) trip(COMB(:,3),:)];
        
        D.Diff = diff(D.Seq,1,2);
        x = any(D.Diff== 0,2);
        D = getrow(D,x~=1);
        
        
        for i = 1:5
            RepFing =  sum(ismember(D.Seq,i),2);
            
            D = getrow(D,RepFing <4);
        end
        

        save('SeqAll.mat','D','trip');
             
    case 'featureModels'            % Select sequences based on criteria imposed, correlations from feature models
        
        load('SeqAll.mat','D','trip');
        A = load('sequenceType.mat');   % all possible doubles, triplets
        
        SeqFinal = {};
        ChunkAmount = [];
        FingRep = [];
        good = false;
        NumbFinalSequenc = 1;
        while good == false;
            for i = 1:30
                good1 = false;
                %make new sequences until the same finger is
                %in the beginning of 2 sequences but not in more than 2
                while good1 == false;
                    seqpick = datasample(D.Seq,3);
                    seqrev =  6-seqpick;
                    seqchos = [seqpick; seqrev];
                    for f = 1:5
                        RepFing(f) =  sum(ismember(seqchos(:,1),f));
                        
                    end
                    if sum(ismember(RepFing,2)) == 3    % 3 repetitions (of the same first fing)
                        good1 = true;
                    end
                    
                end
                SeqChos{i,:} = seqchos;
                
            end
            
            for i = 1:size(SeqChos,1)       % check all the features for all chosen seq; i - set of 6 seq
                
                % Feature models
                F.FirstFing = zeros(6,5);
                F.AllFing = zeros(6,5);
                F.TransFing = zeros(6,25);
                F.Chunks =  zeros(6,52);
                F.TransChunk = zeros(6,25);
                
                % First finger, All fingers
                for u = 1:size(SeqChos{i},1)    % set of 6 seq
                    firstfing = SeqChos{i}(:,1);
                    F.FirstFing(u,firstfing(u)) = 1;
                    for j = 1:5                 % 5 fingers
                        placenumb = find(SeqChos{i}(u,:) == j);
                        F.AllFing(u,j) = length(placenumb);
                    end
                end
                
                
                % Transitions
                Rd.Doubles = A.Sequence(126:150,1:2);
                
                for k = 1:size(SeqChos{i},1)    % set of 6 seq
                    for p = 1:size(Rd.Doubles,1)
                        if any(ismember(SeqChos{i}(k,1:8),Rd.Doubles(p,1)))
                            ind = find(ismember(SeqChos{i}(k,1:8),Rd.Doubles(p,1)));    % find matching digits to first finger in doublet
                            
                            if any(ismember(SeqChos{i}(k,ind+1),Rd.Doubles(p,2)))       % compare the second finger in doublet
                                F.TransFing(k,p) = F.TransFing(k,p) + sum(ismember(SeqChos{i}(k,ind+1),Rd.Doubles(p,2)));
                            end
                        end
                    end
                    
                end
                
                % Chunks
                
                for h = 1:size(SeqChos{i},1)
                    for p = 1:size(trip,1) % all allowed triplets
                        if strcmp(num2str(SeqChos{i}(h,1:3)),num2str(trip(p,1:3))) | strcmp(num2str(SeqChos{i}(h,4:6)),num2str(trip(p,1:3))) | strcmp(num2str(SeqChos{i}(h,7:9)),num2str(trip(p,1:3)))
                            F.Chunks(h,p) = 1;
                        end
                    end
                end
                
                % Transitions between chunks
                
                for m = 1:size(SeqChos{i},1)
                    for p = 1:size(Rd.Doubles,1)    % all possible doubles - 25
                        if strcmp(num2str(SeqChos{i}(m,3:4)),num2str(Rd.Doubles(p,:))) | strcmp(num2str(SeqChos{i}(m,6:7)),num2str(Rd.Doubles(p,:)))
                            F.TransChunk(m,p) = 1;
                        end
                    end
                end
                
                
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
                D.distTrans     = distTransFing;
                D.distChunk     = distChunk;
                
                
                % correlate the models (All/First/Trans/Chunk)
                % correlations with chunk are 0 if chunks in all chosen
                % sequences are different
                c.CorTest(1,i)    = corr(distFirst', distAll'); % First & All
                c.CorTest(2,i)    = corr(distFirst', distTransFing'); %First & Transfing
                %c.CorTest(4,i)    = corr(distFirst', distChunk'); %FirstFinger & Chunk
                c.CorTest(3,i)    = corr(distAll', distTransFing'); % All & TransFing
                %c.CorTest(5,i)    = corr(distAll', distChunk'); %All & Chunk
                %c.CorTest(6,i)    = corr(distTransFing', distChunk'); % Transfing & Chunk
                %c.CorTest(7,i)    = corr(distTransChunk', distChunk'); % TransChunk & Chunk
                
                Dist = addstruct(Dist,D);
                Dist = addstruct(Dist,c);
                
                %how often is each finger pressed
                for o = 1:5
                    AmountFing(i,o)=sum(sum(ismember(SeqChos{i},o)));
                end
                
                %get all the possible chunks of the sequence
                allchunks = sml_findSeq('getChunks',SeqChos{i});
                uniqchunks = unique(allchunks,'rows');
                
                % find if all three chunks start with the same finger
                startFingSame = zeros(size(SeqChos{i},1),1);
                for f = 1:size(SeqChos{i},1)
                    startFingSame(f) = (SeqChos{i}(f,1)==SeqChos{i}(f,4) & SeqChos{i}(f,1)==SeqChos{i}(f,7) & SeqChos{i}(f,4)==SeqChos{i}(f,7));
                end
                
                
                % the three chunks should not start with same finger
                if any(ismember(startFingSame,1))
                    continue
                else
                    
                % none of the 6 sequences should start with the same chunk              
                StartChunk = unique(SeqChos{i}(:,1:3),'rows');
                if size(StartChunk,1) ~=6
                    continue
                else
        
                    % the amount of chunks has to be above 7 but below 12
                    if size(uniqchunks,1) <12 & size(uniqchunks,1) > 7
                        chunkamount = zeros(1,size(uniqchunks,1));
                        
                        %find out how often each chunk is in all the sequences
                        for z = 1:size(uniqchunks,1)
                            if any(ismember(SeqChos{i}(:,1:3),uniqchunks(z,:),'rows')) | ...
                                    any(ismember(SeqChos{i}(:,4:6),uniqchunks(z,:),'rows'))| ...
                                    any(ismember(SeqChos{i}(:,7:9),uniqchunks(z,:),'rows'));
                                
                                chunkamount(z) = sum(ismember(SeqChos{i}(:,1:3),uniqchunks(z,:),'rows')+...
                                    ismember(SeqChos{i}(:,4:6),uniqchunks(z,:),'rows')+...
                                    ismember(SeqChos{i}(:,7:9),uniqchunks(z,:),'rows'));
                            end
                        end
                        
                         %at least one run (chunk) in the sequence set
                        chunkrun = diff(uniqchunks,1,2);
                        
                       if any(ismember(chunkrun(:,1:2),[1 1],'rows')) | ...
                                any(ismember(chunkrun(:,1:2),[-1 -1],'rows'));
                            
                            %the chunk that is represented most in the sequences
                            %can only be in the sequences 4 times or less and the
                            %least represented chunk has to be in the sequences at
                            %least 3 times
                            if max(chunkamount)< 5
                                
                                % each finger has to be at least 6 times in the
                                % sequences but no more than 15 times
                                if  min(AmountFing(i,:)>6) & max(AmountFing(i,:)) < 15
                                    
                                    %all correlations have to be under 0.6 and
                                    %can't be NaNs
                                    if sum((c.CorTest(:,i)<0.6)) == 3 & sum(isnan(c.CorTest(:,i))) == 0
                                        
                                        %save the sequences and the correlations in
                                        %variable
                                        SeqFinal{NumbFinalSequenc} = SeqChos{i,:};
                                        TestCorr(:,NumbFinalSequenc) = c.CorTest(:,i);
                                        ChunkAmount(NumbFinalSequenc) = size(uniqchunks,1);
                                        FingRep(:,NumbFinalSequenc) = AmountFing(i,:)';
                                        
                                        fprintf('Final sequence: %d', NumbFinalSequenc);
                                        
                                        NumbFinalSequenc = NumbFinalSequenc+1;
                                        
                                        
                                        
                                    else
                                        continue;
                                    end             % correlation
                                else
                                    continue;
                                    
                                end                 % finger amount
                            else
                                continue;
                            end                     % chunk amount max. 5
                        else
                            continue;
                        end                         % run in the chunk
                    else
                        continue;
                    end                             % chunk amount - 8-11
                end                                 % start same finger
                end                                 % 6 different starting chunks
            end                                     % for 30 6-sequence combinations
            clear c
            %if we get more than a certain amount of different good
            %sequences we can save them
            if size(SeqFinal,2) >19
                good = 1;
            end
            
        end
        save('SeqFinal.mat','SeqFinal','TestCorr','ChunkAmount','FingRep');
        
    case 'getChunks'                % Separate sequences into chunks - called from featureModels
        Seq = varargin{1};
        for n=1:size(Seq,1);
            t(n*3-2:n*3,:) = [Seq(n,[1:3]);Seq(n,[4:6]);Seq(n,[7:9])];
            
        end
        varargout = {t};
        
    otherwise
        error('no such case!')
        
end


