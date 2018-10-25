function sml_randomSeq(what,varargin)

% Generate random sequences based on stimuli in the experiment
% Seq1/2 - sequence set for group 1/2
% randomSeq - random sequences with 0 known chunks; no rep, no trills
% randomChunk - known chunks arranged in a new way

baseDir = '/Users/eberlot/Documents/Data/SuperMotorLearning/stimuli_struct';
cd (baseDir);
rng('shuffle');
load('AllTrip.mat');

Seq1 = [1 5 3 5 3 4 3 1 5; ...
        1 5 2 1 4 5 3 4 2; ...
        3 4 2 1 4 5 1 5 2; ...
        3 1 5 3 4 2 1 4 5; ...
        5 3 4 5 4 2 1 5 3; ...
        5 4 2 3 4 2 3 1 5];

Seq2 = [5 1 3 1 3 2 3 5 1; ...
        5 1 4 5 2 1 3 2 4; ...
        3 2 4 5 2 1 5 1 4; ...
        3 5 1 3 2 4 5 2 1; ...
        1 3 2 1 2 4 5 1 3; ...
        1 2 4 3 2 4 3 5 1];

Chunks1 = [1 5 3; 5 3 4; 3 1 5; 1 5 2; 1 4 5; 3 4 2; 5 4 2];
Chunks2 = [5 1 3; 1 3 2; 3 5 1; 5 1 4; 5 2 1; 3 2 4; 1 2 4];

RandSet = [Seq1(:); Seq2(:)];

numRandS = 6*6;     % 6 sets needed (2BL, 4 tests)
numRandC = 6*4;     % 4 sessions

switch(what)
    case 'randomSeq'   % create random sequences with 0 chunks in common
       crit1 = ismember(AllChunksOrd,Chunks1,'rows');
       crit2 = ismember(AllChunksOrd,Chunks2,'rows');
       crit = crit1+crit2;
       
       trip = AllChunksOrd(crit==0,:); % different triplets from experiment
       
       for r = 1:numRandS         
           good1=false;
            while  good1 == false % no two chunks repeated, no finger rep
                c = (datasample(trip,3,'Replace',false))';    
                RandomSeq(r,:) = c(:);
                if (sum(ismember(diff(RandomSeq(r,:)),0))==0)
                    good1=1;
                end
            end         
       end

       save('RandomSeq.mat','RandomSeq');
    case 'randomChunk'  % create sequences with known chunks for G1/2
        
        for r = 1:numRandC      % Group 1  
           good1=false;
            while  good1 == false % no two chunks repeated, no finger rep
                c = (datasample(Chunks1,3,'Replace',false))';    
                RandomChunkSeq1(r,:) = c(:);
                if (sum(ismember(diff(RandomChunkSeq1(r,:)),0))==0)
                    if sum(ismember(Seq1,RandomChunkSeq1(r,:),'rows')==0)
                        good1=1;
                    end
                end
            end         
        end
       
        for r = 1:numRandC      % Group 2  
           good1=false;
            while  good1 == false % no two chunks repeated, no finger rep
                c = (datasample(Chunks2,3,'Replace',false))';    
                RandomChunkSeq2(r,:) = c(:);
                if (sum(ismember(diff(RandomChunkSeq2(r,:)),0))==0)
                    if sum(ismember(Seq2,RandomChunkSeq2(r,:),'rows')==0)
                        good1=1;
                    end
                end
            end         
        end

       save('RandomChunkSeq1.mat','RandomChunkSeq1');
       save('RandomChunkSeq2.mat','RandomChunkSeq2');
       
    otherwise
        error('no such case!')

end
