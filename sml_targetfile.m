function varargout=sml_targetfile(what,varargin)
% make target files

%%%%%%%%% USEFUL INFO  %%%%%%%%%%%%

%%%%%% groups %%%%%%
% G1, G2
% both training with the right hand

%%%%%% stimuli %%%%%%%
% chunks - triplets
% 6 sequences of 9 presses

%%%%% seqType %%%%%%%  
% 1 - training 
% 2 - untrained (other group)
% 3 - finger mapping
% 4 - random seq
% 5 - random chunks
% 6 - other hand (trained intrinsic)
% 7 - other hand (trained extrinsic)
% 8 - other hand (random intrinsic)
% 9 - other hand (random extrinsic)

%%%%% seqNumb %%%%%%
% 1-6 trained
% 7-12 other group (untrained)
% 13-17 finger mapping
% 18-53 random sequences
% 54-77 random chunks
% 78-83 trained seq other hand
% 84-113 random seq other hand
% ADD CHUNKS

%%%%%% scan sessions %%%%%%%
% 8 functional runs with train seq (6) and untrain seq (6)
% 2 localisers (finger mapping)
% 3 scanning sessions with metronome, one without (ScanGo)

%%%%%% isScan %%%%%%%%%
% 0 - training, test blocks - behaviour
% 1 - scanning sessions
% 2 - behavioural prescan sessions which include a metronome (trials paced)

%%%%% blockType %%%%%%% 
% which type of the task is done
% 1 - BeforeLearn
% 2 - BeforeLearnMet
% 3 - Scan
% 4 - ScanLoc
% 5 - ChunkLearn
% 6 - SeqLearn
% 7 - PreScanMet
% 8 - PreScanGo
% 9 - ScanGo
% 10 - OtherHandRand
% 11 - SeqTest
% 12 - OtherHandTest
% 13 - BeforeLearnFingMap


%%%%%%%%%%%%%% INITIALISING COMMON STUFF %%%%%%%%%%%%%%%%%%

baseDir = '/Users/eberlot/Documents/RobotCode/SuperMotorLearning/target/NewTarget';
seqDir = '/Users/eberlot/Documents/Data/SuperMotorLearning/stimuli_struct';

cd(baseDir);

rng('default');

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

ChunkAmount = 1:7;
ChunkNumb = [114 115 116 117 118 119 120]; 

Chunks1 = [1 5 3; 5 3 4; 3 1 5; 1 5 2; 1 4 5; 3 4 2; 5 4 2];
Chunks2 = [5 1 3; 1 3 2; 3 5 1; 5 1 4; 5 2 1; 3 2 4; 1 2 4];

SeqChunks1 = [1,2,3; ...
    4,5,6; ...
    6,5,4; ...
    3,6,5; ...
    2,7,1; ...
    7,6,3];

SeqChunks2 = [1,2,3; ...
    4,5,6; ...
    6,5,4; ...
    3,6,5; ...
    2,7,1; ...
    7,6,3];

FingMap = [1 1 1 1 1 1 1 1 1; ...
    2 2 2 2 2 2 2 2 2; ...
    3 3 3 3 3 3 3 3 3; ...
    4 4 4 4 4 4 4 4 4; ...
    5 5 5 5 5 5 5 5 5];

% Used for BL, PreScan, PreScanMet, tests - 36 (6 sets of 6 seq)
load(fullfile(seqDir,'RandomSeq.mat'));

% Used for SeqTests
R1 = load(fullfile(seqDir,'RandomChunkSeq1.mat'));
R2 = load(fullfile(seqDir,'RandomChunkSeq2.mat'));
RandChunkSeq1 = R1.RandomChunkSeq1;
RandChunkSeq2 = R2.RandomChunkSeq2;

C1 = {'£'; '£';'£';'£';'£';'£'};
C2 = {'£'; '£'; '£';'£'; '£';'£';'£';'£';'£';'£';'£'};
C3 = {'***'};
C3Seq = {'*********'};

% Make sure numbers are adapted if changing num of rand seq etc.
SeqNumb = 1:6;
SeqNumbUntrain = 7:12;
SeqNumbFingMap = 13:17;
SeqRandNumb = (1:length(RandomSeq))+17; % 18-53
SeqRandChunkNum = (1:length(RandChunkSeq1))+(length(RandomSeq))+17;  %54-77
SeqNumbTrainOtherHand = (1:6)+17+length(RandomSeq)+length(RandChunkSeq1);   % 78-83 (intr., extr. given same number)
SeqNumbRandOtherHand = 84:113;  % 30 seq in total - extr./intr. given same seq number, different in seqType


OrderFields = {'startTime','blockType','seqNumb','seqType','press1','press2','press3','press4','press5','press6','press7','press8','press9','hand','cueS','cueC','cueP','iti','sounds','isScan','isMetronome','lastTrial','FoSEX','announceTime','feedbackTime','ScanSess'};

SeqG1 = Seq1;
SeqG2 = Seq2;
SeqChunksG1 = SeqChunks1;
SeqChunksG2 = SeqChunks2;
ChunksG1 = Chunks1;
ChunksG2 = Chunks2;

switch what
    case 'ChunkLearn'      
    % 70 trials per block, 10 chunk rep per run
    % each chunk repeated twice in a row, 2nd time *** (recall)

    Trials = 1:70;  % 7 chunks - each repeated 10 times (5 doubles)

    runs = 44;

    for e = 1:runs
        ChunkLearn_G1.cueS = cellstr(repmat('£',length(Trials),1));
        ChunkLearn_G2.cueS = cellstr(repmat('£',length(Trials),1));

        ChunkLearn_G1.iti(1:length(Trials),:) = 500;
        ChunkLearn_G2.iti(1:length(Trials),:) = 500;
        ChunkLearn_G1.blockType(1:length(Trials),:) = 5;
        ChunkLearn_G2.blockType(1:length(Trials),:) = 5;

        ChunkLearn_G1.hand(1:length(Trials),:) = 2; % both groups right hand
        ChunkLearn_G2.hand(1:length(Trials),:) = 2;

        ChunkLearn_G1.announceTime(1:length(Trials),:) = 0;
        ChunkLearn_G2.announceTime(1:length(Trials),:) = 0;
        ChunkLearn_G1.feedbackTime(1:length(Trials),:) = 400;
        ChunkLearn_G2.feedbackTime(1:length(Trials),:) = 400;

        ChunkLearn_G1.sounds(1:length(Trials),:) = 1;
        ChunkLearn_G2.sounds(1:length(Trials),:) = 1;

        %randomly choose order
        X = sample_wor(ChunkAmount,length(ChunkAmount));
        p = 1;
        l = 1;
        for i = 1:(length(Trials)/2)

            % press cue
            ChunkLearn_G1.cueP(p,:) = regexprep(cellstr(num2str(ChunksG1(X(l),:))),'\s','');
            ChunkLearn_G1.cueP(p+1,:) = C3; % stars every second trial

            ChunkLearn_G2.cueP(p,:) = regexprep(cellstr(num2str(ChunksG2(X(l),:))),'\s','');
            ChunkLearn_G2.cueP(p+1,:) = C3;

            ChunkLearn_G1.press1(p:p+1,:) = ChunksG1(X(l),1);
            ChunkLearn_G1.press2(p:p+1,:) = ChunksG1(X(l),2);
            ChunkLearn_G1.press3(p:p+1,:) = ChunksG1(X(l),3);

            ChunkLearn_G2.press1(p:p+1,:) = ChunksG2(X(l),1);
            ChunkLearn_G2.press2(p:p+1,:) = ChunksG2(X(l),2);
            ChunkLearn_G2.press3(p:p+1,:) = ChunksG2(X(l),3);

            ChunkLearn_G1.seqNumb(p:p+1,:) = ChunkNumb(X(l));
            ChunkLearn_G2.seqNumb(p:p+1,:) = ChunkNumb(X(l));

            
            %chunk cue
            ChunkLearn_G1.cueC(p:p+1,:) = C2(X(l),:) ;
            ChunkLearn_G2.cueC(p:p+1,:) = C2(X(l),:) ;

            p = p+2;
            l = l+1;
            if l == 8;
                l = 1;
                X = sample_wor(ChunkAmount,length(ChunkAmount));
            end
        end
        ChunkLearn_G1.press4 = zeros(length(Trials),1);
        ChunkLearn_G1.press5 = zeros(length(Trials),1);
        ChunkLearn_G1.press6 = zeros(length(Trials),1);
        ChunkLearn_G1.press7 = zeros(length(Trials),1);
        ChunkLearn_G1.press8 = zeros(length(Trials),1);
        ChunkLearn_G1.press9 = zeros(length(Trials),1);

        ChunkLearn_G2.press4 = zeros(length(Trials),1);
        ChunkLearn_G2.press5 = zeros(length(Trials),1);
        ChunkLearn_G2.press6 = zeros(length(Trials),1);
        ChunkLearn_G2.press7 = zeros(length(Trials),1);
        ChunkLearn_G2.press8 = zeros(length(Trials),1);
        ChunkLearn_G2.press9 = zeros(length(Trials),1);

        ChunkLearn_G1.seqType = zeros(length(Trials),1);
        ChunkLearn_G2.seqType = zeros(length(Trials),1);

        ChunkLearn_G1.isScan = zeros(length(Trials),1);
        ChunkLearn_G2.isScan = zeros(length(Trials),1);

        ChunkLearn_G1.ScanSess = zeros(length(Trials),1);
        ChunkLearn_G2.ScanSess = zeros(length(Trials),1);

        ChunkLearn_G1.seqType = zeros(length(Trials),1);
        ChunkLearn_G2.seqType = zeros(length(Trials),1);

        ChunkLearn_G1.FoSEX = repmat([1:2]',length(Trials)/2,1);
        ChunkLearn_G2.FoSEX = repmat([1:2]',length(Trials)/2,1);

        ChunkLearn_G1.isMetronome = zeros(length(Trials),1);
        ChunkLearn_G2.isMetronome = zeros(length(Trials),1);

        ChunkLearn_G1.lastTrial(1:length(Trials),1) = 0;
        ChunkLearn_G1.lastTrial(length(Trials),:) = 1;
        ChunkLearn_G2.lastTrial(1:length(Trials),1) = 0;
        ChunkLearn_G2.lastTrial(length(Trials),:) = 1;

        ChunkLearn_G1.startTime = zeros(length(Trials),1);
        ChunkLearn_G2.startTime = zeros(length(Trials),1);

        dirname = fullfile(baseDir,'ChunkLearn');
        dircheck(dirname);

        dsave(fullfile(baseDir,['ChunkLearn/','sml_CLB_G1_B', num2str(e), '.tgt']),orderfields(ChunkLearn_G1,OrderFields));
        dsave(fullfile(baseDir,['ChunkLearn/','sml_CLB_G2_B', num2str(e), '.tgt']),orderfields(ChunkLearn_G2,OrderFields));

    end
        
    case 'SeqLearn'      
    % every sequence repeated twice in a row, 2nd time ********* (recall) - FoSEx
    % each sequence repeated 4x per sequence (2 pairs FoSEx)
    % 24 trials per block

    Trials = 1:24;       

    runs = 174;

    for e = 1:runs

        SequenceLearn_G2.cueC = cellstr(repmat('£',length(Trials),1));
        SequenceLearn_G1.cueC = cellstr(repmat('£',length(Trials),1));

        SequenceLearn_G1.iti(1:length(Trials),:) = 500;
        SequenceLearn_G2.iti(1:length(Trials),:) = 500;

        SequenceLearn_G1.hand(1:length(Trials),:) = 2;
        SequenceLearn_G2.hand(1:length(Trials),:) = 2;

        SequenceLearn_G1.blockType(1:length(Trials),:) = 6;
        SequenceLearn_G2.blockType(1:length(Trials),:) = 6;

        SequenceLearn_G1.announceTime(1:length(Trials),:) = 0;
        SequenceLearn_G2.announceTime(1:length(Trials),:) = 0;
        SequenceLearn_G1.feedbackTime(1:length(Trials),:) = 400;
        SequenceLearn_G2.feedbackTime(1:length(Trials),:) = 400;

        SequenceLearn_G1.sounds(1:length(Trials),:) = 1;
        SequenceLearn_G2.sounds(1:length(Trials),:) = 1;

        X = sample_wor(SeqNumb,length(SeqNumb));
        p = 1;
        l = 1;
        for i = 1:(length(Trials)/2)

            SequenceLearn_G1.cueP(p,:) = regexprep(cellstr(num2str(SeqG1(X(l),:))),'\s','');
            SequenceLearn_G1.cueP(p+1,:) = C3Seq;   % stars for every second trial
            SequenceLearn_G2.cueP(p,:) = regexprep(cellstr(num2str(SeqG2(X(l),:))),'\s','');
            SequenceLearn_G2.cueP(p+1,:) = C3Seq;
            SequenceLearn_G1.cueS(p:p+1,:) = C1(X(l));  % £
            SequenceLearn_G2.cueS(p:p+1,:) = C1(X(l));
            SequenceLearn_G1.cueC(p:p+1,:) = regexprep(cellstr(strjoin(C2(SeqChunksG1(X(l),:))')),'\s',''); % £
            SequenceLearn_G2.cueC(p:p+1,:) = regexprep(cellstr(strjoin(C2(SeqChunksG2(X(l),:))')),'\s','');

            SequenceLearn_G1.press1(p:p+1,:) = SeqG1(X(l),1);
            SequenceLearn_G1.press2(p:p+1,:) = SeqG1(X(l),2);
            SequenceLearn_G1.press3(p:p+1,:) = SeqG1(X(l),3);
            SequenceLearn_G1.press4(p:p+1,:) = SeqG1(X(l),4);
            SequenceLearn_G1.press5(p:p+1,:) = SeqG1(X(l),5);
            SequenceLearn_G1.press6(p:p+1,:) = SeqG1(X(l),6);
            SequenceLearn_G1.press7(p:p+1,:) = SeqG1(X(l),7);
            SequenceLearn_G1.press8(p:p+1,:) = SeqG1(X(l),8);
            SequenceLearn_G1.press9(p:p+1,:) = SeqG1(X(l),9);

            SequenceLearn_G2.press1(p:p+1,:) = SeqG2(X(l),1);
            SequenceLearn_G2.press2(p:p+1,:) = SeqG2(X(l),2);
            SequenceLearn_G2.press3(p:p+1,:) = SeqG2(X(l),3);
            SequenceLearn_G2.press4(p:p+1,:) = SeqG2(X(l),4);
            SequenceLearn_G2.press5(p:p+1,:) = SeqG2(X(l),5);
            SequenceLearn_G2.press6(p:p+1,:) = SeqG2(X(l),6);
            SequenceLearn_G2.press7(p:p+1,:) = SeqG2(X(l),7);
            SequenceLearn_G2.press8(p:p+1,:) = SeqG2(X(l),8);
            SequenceLearn_G2.press9(p:p+1,:) = SeqG2(X(l),9);

            SequenceLearn_G1.seqNumb(p:p+1,:) = SeqNumb(X(l));
            SequenceLearn_G2.seqNumb(p:p+1,:) = SeqNumb(X(l));

            SequenceLearn_G1.isScan = zeros(length(Trials),1);
            SequenceLearn_G2.isScan = zeros(length(Trials),1);

            SequenceLearn_G1.ScanSess = zeros(length(Trials),1);
            SequenceLearn_G2.ScanSess = zeros(length(Trials),1);
           
            % one for trained sequence
            SequenceLearn_G1.seqType = ones(length(Trials),1);
            SequenceLearn_G2.seqType = ones(length(Trials),1);

            SequenceLearn_G1.FoSEX = repmat([1:2]',length(Trials)/2,1);
            SequenceLearn_G2.FoSEX = repmat([1:2]',length(Trials)/2,1);

            SequenceLearn_G1.isMetronome = zeros(length(Trials),1);
            SequenceLearn_G2.isMetronome = zeros(length(Trials),1);

            SequenceLearn_G1.lastTrial(1:length(Trials),1) = 0;
            SequenceLearn_G1.lastTrial(length(Trials),:) = 1;
            SequenceLearn_G2.lastTrial(1:length(Trials),1) = 0;
            SequenceLearn_G2.lastTrial(length(Trials),:) = 1;

            SequenceLearn_G1.startTime = zeros(length(Trials),1);
            SequenceLearn_G2.startTime = zeros(length(Trials),1);

            p = p+2;
            l = l+1;
            if l == length(SeqNumb)+1;
                l = 1;
                X = sample_wor(SeqNumb,length(SeqNumb));
            end
        end

        dirname = fullfile(baseDir,'SeqLearn');
        dircheck(dirname);

        dsave(fullfile(baseDir,['SeqLearn/','sml_SLB_G1_B' ,num2str(e), '.tgt']),orderfields(SequenceLearn_G1,OrderFields));
        dsave(fullfile(baseDir,['SeqLearn/','sml_SLB_G2_B' ,num2str(e), '.tgt']),orderfields(SequenceLearn_G2,OrderFields));
    end
    
    case 'BeforeLearn' 
    % first session - to get familiar with the keyboard/task
    % random sequences - first 12 (2 sets of 6 like for scanning)
    % sequence presented twice in a row 
    % each sequence pair (FoSEx) repeated twice per block - 4x6 = 24 trials per block

    runs = 16;
    trial_run = 24;

    Trials = 1:(runs*trial_run);  % 4 blocks of 24

    SeqG1 = RandomSeq(1:12,:);
    SeqG2 = RandomSeq(1:12,:);

    X = [sample_wor(1:size(SeqG1,1),size(SeqG1,1));...
        sample_wor(1:size(SeqG1,1),size(SeqG1,1));...
        sample_wor(1:size(SeqG1,1),size(SeqG1,1));...
        sample_wor(1:size(SeqG1,1),size(SeqG1,1));...
        sample_wor(1:size(SeqG1,1),size(SeqG1,1));...
        sample_wor(1:size(SeqG1,1),size(SeqG1,1));...
        sample_wor(1:size(SeqG1,1),size(SeqG1,1));...
        sample_wor(1:size(SeqG1,1),size(SeqG1,1));...
        sample_wor(1:size(SeqG1,1),size(SeqG1,1));...
        sample_wor(1:size(SeqG1,1),size(SeqG1,1));...
        sample_wor(1:size(SeqG1,1),size(SeqG1,1));...
        sample_wor(1:size(SeqG1,1),size(SeqG1,1));...
        sample_wor(1:size(SeqG1,1),size(SeqG1,1));...
        sample_wor(1:size(SeqG1,1),size(SeqG1,1));...
        sample_wor(1:size(SeqG1,1),size(SeqG1,1));...
        sample_wor(1:size(SeqG1,1),size(SeqG1,1))];

    BeforeLearn_G1.cueC = cellstr(repmat('£',length(Trials),1));
    BeforeLearn_G2.cueC = cellstr(repmat('£',length(Trials),1));

    BeforeLearn_G1.iti(1:length(Trials),:) = 500;
    BeforeLearn_G2.iti(1:length(Trials),:) = 500;

    BeforeLearn_G1.hand(1:length(Trials),:) = 2;
    BeforeLearn_G2.hand(1:length(Trials),:) = 2;

    BeforeLearn_G1.blockType(1:length(Trials),:) = 1;
    BeforeLearn_G2.blockType(1:length(Trials),:) = 1;

    BeforeLearn_G1.announceTime(1:length(Trials),:) = 0;
    BeforeLearn_G2.announceTime(1:length(Trials),:) = 0;
    BeforeLearn_G1.feedbackTime(1:length(Trials),:) = 400;
    BeforeLearn_G2.feedbackTime(1:length(Trials),:) = 400;

    BeforeLearn_G1.sounds(1:length(Trials),:) = 1;
    BeforeLearn_G2.sounds(1:length(Trials),:) = 1;

    BeforeLearn_G1.isScan(1:length(Trials),:) = 0;
    BeforeLearn_G2.isScan(1:length(Trials),:) = 0;

    BeforeLearn_G1.ScanSess(1:length(Trials),:) = 0;
    BeforeLearn_G2.ScanSess(1:length(Trials),:) = 0;

    BeforeLearn_G1.isMetronome(1:length(Trials),:) = 0;
    BeforeLearn_G2.isMetronome(1:length(Trials),:) = 0;

    p = 1;
    l = 1;

    for i = 1:(length(Trials)/2)

        BeforeLearn_G1.cueP(p:p+1,:) = regexprep(cellstr(num2str(SeqG1(X(l),:))),'\s','');
        BeforeLearn_G2.cueP(p:p+1,:) = regexprep(cellstr(num2str(SeqG2(X(l),:))),'\s','');
        BeforeLearn_G1.cueS(p:p+1,:) = C1(1);  %£
        BeforeLearn_G2.cueS(p:p+1,:) = C1(1);

        BeforeLearn_G1.press1(p:p+1,:) = SeqG1(X(l),1);
        BeforeLearn_G1.press2(p:p+1,:) = SeqG1(X(l),2);
        BeforeLearn_G1.press3(p:p+1,:) = SeqG1(X(l),3);
        BeforeLearn_G1.press4(p:p+1,:) = SeqG1(X(l),4);
        BeforeLearn_G1.press5(p:p+1,:) = SeqG1(X(l),5);
        BeforeLearn_G1.press6(p:p+1,:) = SeqG1(X(l),6);
        BeforeLearn_G1.press7(p:p+1,:) = SeqG1(X(l),7);
        BeforeLearn_G1.press8(p:p+1,:) = SeqG1(X(l),8);
        BeforeLearn_G1.press9(p:p+1,:) = SeqG1(X(l),9);

        BeforeLearn_G2.press1(p:p+1,:) = SeqG2(X(l),1);
        BeforeLearn_G2.press2(p:p+1,:) = SeqG2(X(l),2);
        BeforeLearn_G2.press3(p:p+1,:) = SeqG2(X(l),3);
        BeforeLearn_G2.press4(p:p+1,:) = SeqG2(X(l),4);
        BeforeLearn_G2.press5(p:p+1,:) = SeqG2(X(l),5);
        BeforeLearn_G2.press6(p:p+1,:) = SeqG2(X(l),6);
        BeforeLearn_G2.press7(p:p+1,:) = SeqG2(X(l),7);
        BeforeLearn_G2.press8(p:p+1,:) = SeqG2(X(l),8);
        BeforeLearn_G2.press9(p:p+1,:) = SeqG2(X(l),9);

        BeforeLearn_G1.seqNumb(p:p+1,:) = SeqRandNumb(X(l));
        BeforeLearn_G2.seqNumb(p:p+1,:) = SeqRandNumb(X(l));

        p = p+2;
        l = l+1;

    end

    BeforeLearn_G1.isScan = zeros(length(Trials),1);
    BeforeLearn_G2.isScan = zeros(length(Trials),1);

    BeforeLearn_G1.seqType = ones(length(Trials),1)*4;  %4 - random
    BeforeLearn_G2.seqType = ones(length(Trials),1)*4;  %4 - random
            
    BeforeLearn_G1.FoSEX = repmat([1:2]',length(Trials)/2,1);
    BeforeLearn_G2.FoSEX = repmat([1:2]',length(Trials)/2,1);

    BeforeLearn_G1.isMetronome = zeros(length(Trials),1);
    BeforeLearn_G2.isMetronome = zeros(length(Trials),1);

    BeforeLearn_G1.lastTrial(1:length(Trials),1) = 0;
    BeforeLearn_G1.lastTrial([24,48,72,96,120,144,168,192,216,240,264,288,312,336,360,384],:) = 1;      % every 24th trial - last trial in a block
    BeforeLearn_G2.lastTrial(1:length(Trials),1) = 0;
    BeforeLearn_G2.lastTrial([24,48,72,96,120,144,168,192,216,240,264,288,312,336,360,384],:) = 1;

    BeforeLearn_G1.startTime = zeros(length(Trials),1);
    BeforeLearn_G2.startTime = zeros(length(Trials),1);

    B1 = 1;
    B2 = trial_run;

    dirname = fullfile(baseDir,'BeforeLearn');
    dircheck(dirname);
    
    for i = 1:runs
        dsave(fullfile(baseDir,['BeforeLearn/','sml_BL_G1_B', num2str(i), '.tgt']),orderfields(getrow(BeforeLearn_G1,B1:B2),OrderFields));
        dsave(fullfile(baseDir,['BeforeLearn/','sml_BL_G2_B', num2str(i), '.tgt']),orderfields(getrow(BeforeLearn_G2,B1:B2),OrderFields));
        B1 = B1+trial_run;
        B2 = B2+trial_run;

    end
        
    case 'BeforeLearnMet'
    % to get familiar with the metronome, fmri setup (bed)
    % random sequences - 12 (2 sets of 6)
    % sequences presented twice in a row 
    % each sequence presented 6x per block (3 repetitions) 
    % 12 x 6 = 72 trials per block

    Trials = 1:72;
    runs = 2;

    SeqG1 = RandomSeq(1:12,:);
    SeqG2 = RandomSeq(1:12,:);

    for e = 1:runs
        X = [sample_wor(1:size(SeqG1,1),size(SeqG1,1));sample_wor(1:size(SeqG1,1),size(SeqG1,1));sample_wor(1:size(SeqG1,1),size(SeqG1,1))];

        BLM_G1.cueC = cellstr(repmat('£',length(Trials),1));
        BLM_G2.cueC = cellstr(repmat('£',length(Trials),1));

        BLM_G1.iti(1:length(Trials),:) = 0;
        BLM_G2.iti(1:length(Trials),:) = 0;

        BLM_G1.hand(1:length(Trials),:) = 2;
        BLM_G2.hand(1:length(Trials),:) = 2;

        BLM_G1.blockType(1:length(Trials),:) = 2;
        BLM_G2.blockType(1:length(Trials),:) = 2;

        BLM_G1.announceTime(1:length(Trials),:) = 1000;
        BLM_G2.announceTime(1:length(Trials),:) = 1000;
        BLM_G1.feedbackTime(1:length(Trials),:) = 400;
        BLM_G2.feedbackTime(1:length(Trials),:) = 400;

        BLM_G1.sounds(1:length(Trials),:) = 1;
        BLM_G2.sounds(1:length(Trials),:) = 1;

        BLM_G1.isScan(1:length(Trials),:) = 2;
        BLM_G2.isScan(1:length(Trials),:) = 2;

        BLM_G1.ScanSess(1:length(Trials),:) = 0;
        BLM_G2.ScanSess(1:length(Trials),:) = 0;

        BLM_G1.isMetronome(1:length(Trials),:) = 1;
        BLM_G2.isMetronome(1:length(Trials),:) = 1;

        p = 1;
        l = 1;

        trialLength = 5000;
        startRun = 5000;

        RandBreaks = sort(sample_wor(1:length(Trials),5));
        while (sum(find(RandBreaks == 1)) || sum(rem(RandBreaks,2))<5) % breaks cannot be on first trial, between two consecutive repetitions (only odd trials)
            RandBreaks = sort(sample_wor(1:length(Trials),5));
        end

        for i = 1:length(Trials)

            if any(ismember(RandBreaks,i))
                startRun = startRun+10000;   % pause of 10000 - for baseline
            end
            BLM_G1.startTime(i,:) = startRun;
            BLM_G2.startTime(i,:) = startRun;

            startRun = startRun+trialLength;
        end
        
        for i = 1:(length(Trials)/2)

            BLM_G1.cueP(p:p+1,:) = regexprep(cellstr(num2str(SeqG1(X(l),:))),'\s','');
            BLM_G2.cueP(p:p+1,:) = regexprep(cellstr(num2str(SeqG2(X(l),:))),'\s','');
            BLM_G1.cueS(p:p+1,:) = C1(1);  %£
            BLM_G2.cueS(p:p+1,:) = C1(1);

            BLM_G1.press1(p:p+1,:) = SeqG1(X(l),1);
            BLM_G1.press2(p:p+1,:) = SeqG1(X(l),2);
            BLM_G1.press3(p:p+1,:) = SeqG1(X(l),3);
            BLM_G1.press4(p:p+1,:) = SeqG1(X(l),4);
            BLM_G1.press5(p:p+1,:) = SeqG1(X(l),5);
            BLM_G1.press6(p:p+1,:) = SeqG1(X(l),6);
            BLM_G1.press7(p:p+1,:) = SeqG1(X(l),7);
            BLM_G1.press8(p:p+1,:) = SeqG1(X(l),8);
            BLM_G1.press9(p:p+1,:) = SeqG1(X(l),9);

            BLM_G2.press1(p:p+1,:) = SeqG2(X(l),1);
            BLM_G2.press2(p:p+1,:) = SeqG2(X(l),2);
            BLM_G2.press3(p:p+1,:) = SeqG2(X(l),3);
            BLM_G2.press4(p:p+1,:) = SeqG2(X(l),4);
            BLM_G2.press5(p:p+1,:) = SeqG2(X(l),5);
            BLM_G2.press6(p:p+1,:) = SeqG2(X(l),6);
            BLM_G2.press7(p:p+1,:) = SeqG2(X(l),7);
            BLM_G2.press8(p:p+1,:) = SeqG2(X(l),8);
            BLM_G2.press9(p:p+1,:) = SeqG2(X(l),9);
            
            BLM_G1.seqNumb(p:p+1,:) = SeqRandNumb(X(l));
            BLM_G2.seqNumb(p:p+1,:) = SeqRandNumb(X(l));

            p = p+2;
            l = l+1;
        end

        BLM_G1.seqType = ones(length(Trials),1)*4;  %4 - random
        BLM_G2.seqType = ones(length(Trials),1)*4;  %4 - random

        BLM_G1.FoSEX = repmat([1:2]',length(Trials)/2,1);
        BLM_G2.FoSEX = repmat([1:2]',length(Trials)/2,1);

        BLM_G1.lastTrial(1:length(Trials),1) = 0;
        BLM_G1.lastTrial(length(Trials),:) = 1;
        BLM_G2.lastTrial(1:length(Trials),1) = 0;
        BLM_G2.lastTrial(length(Trials),:) = 1;

        dirname = fullfile(baseDir,'BeforeLearnMet');
        dircheck(dirname);

        dsave(fullfile(baseDir,['BeforeLearnMet/','sml_BLM_G1_B', num2str(e), '.tgt']),orderfields(BLM_G1,OrderFields));
        dsave(fullfile(baseDir,['BeforeLearnMet/','sml_BLM_G2_B', num2str(e), '.tgt']),orderfields(BLM_G2,OrderFields));

    end
        
    case 'PreScanMet'
    % trained and random sequences with metronome (12 in total)
    % before 2nd, 3rd scanning sessions 

    sequencesG1 = [SeqG1;RandomSeq(1:6,:)]; 
    sequencesG2 = [SeqG2;RandomSeq(1:6,:)];

    runs=2;

    seqforRandG1 = [sequencesG1;sequencesG1;sequencesG1];
    seqforRandG2 = [sequencesG2;sequencesG2;sequencesG2];

    seqG1 = kron(seqforRandG1,ones(2,1));
    seqG2 = kron(seqforRandG2,ones(2,1));

    randomnumbG1(:,1) = [1:36]';
    randomnumbG1(:,2) = [1:6,1:6,1:6,7:12,7:12,7:12]';

    randomnumbG1(:,1) = [1:36]';
    randomnumbG2(:,2) = [1:6,1:6,1:6,7:12,7:12,7:12]';

    % Group 1
    PreScanMetG1.cueC = cellstr(repmat('£',length(seqG1),1));
    PreScanMetG1.cueS = cellstr(repmat('£',length(seqG1),1));

    PreScanMetG1.hand(1:length(seqG1),:) = 2;

    PreScanMetG1.iti(1:length(seqG1),:) = 0;
    PreScanMetG1.sounds(1:length(seqG1),:) = 1;

    PreScanMetG1.blockType(1:length(seqG1),:) = 7;

    PreScanMetG1.isScan(1:length(seqG1),:) = 2;
    PreScanMetG1.ScanSess(1:length(seqG2),:) = 0;
    PreScanMetG1.isMetronome(1:length(seqG1),:) = 1;

    PreScanMetG1.FoSEX(1:2:length(seqG1),:) = 1;
    PreScanMetG1.FoSEX(2:2:length(seqG1),:) = 2;

    PreScanMetG1.lastTrial(1:length(seqG1),:) = 0;
    PreScanMetG1.lastTrial(length(seqG1),:) = 1;

    PreScanMetG1.announceTime(1:length(seqG1),:) = 1000;
    PreScanMetG1.feedbackTime(1:length(seqG1),:) = 400;

    % Group 2
    PreScanMetG2.cueC = cellstr(repmat('£',length(seqG2),1));
    PreScanMetG2.cueS = cellstr(repmat('£',length(seqG2),1));

    PreScanMetG2.hand(1:length(seqG2),:) = 2;

    PreScanMetG2.iti(1:length(seqG2),:) = 0;
    PreScanMetG2.sounds(1:length(seqG2),:) = 1;

    PreScanMetG2.blockType(1:length(seqG2),:) = 7;

    PreScanMetG2.isScan(1:length(seqG2),:) = 2;
    PreScanMetG2.ScanSess(1:length(seqG2),:) = 0;
    PreScanMetG2.isMetronome(1:length(seqG2),:) = 1;

    PreScanMetG2.FoSEX(1:2:length(seqG2),:) = 1;
    PreScanMetG2.FoSEX(2:2:length(seqG2),:) = 2;

    PreScanMetG2.lastTrial(1:length(seqG2),:) = 0;
    PreScanMetG2.lastTrial(length(seqG2),:) = 1;

    PreScanMetG2.announceTime(1:length(seqG2),:) = 1000;
    PreScanMetG2.feedbackTime(1:length(seqG2),:) = 400;

    for e = 1:runs

        X1 = [sample_wor(1:size(randomnumbG1,1),size(randomnumbG1,1))];
        X2 = kron(X1,ones(2,1));    % for all trials - with repetitions

        trialLength = 5000;
        startRun = 5000;

        RandBreaks = sort(sample_wor(1:length(X2),5));
        while (sum(find(RandBreaks == 1)) || sum(rem(RandBreaks,2))<5) % breaks cannot be on first trial, between two consecutive repetitions (only odd trials)
            RandBreaks = sort(sample_wor(1:length(X2),5));
        end

        for i = 1:length(X2)

            if any(ismember(RandBreaks,i))
                startRun = startRun+10000;   % pause of 10000 - for baseline
            end
            PreScanMetG1.startTime(i,:) = startRun;
            PreScanMetG2.startTime(i,:) = startRun;

            startRun = startRun+trialLength;
        end

        p = 1;
        l = 1;

        for i = 1:length(X2)
            PreScanMetG1.press1(p,:) = sequencesG1(randomnumbG1(X2(l),2),1);
            PreScanMetG1.press2(p,:) = sequencesG1(randomnumbG1(X2(l),2),2);
            PreScanMetG1.press3(p,:) = sequencesG1(randomnumbG1(X2(l),2),3);
            PreScanMetG1.press4(p,:) = sequencesG1(randomnumbG1(X2(l),2),4);
            PreScanMetG1.press5(p,:) = sequencesG1(randomnumbG1(X2(l),2),5);
            PreScanMetG1.press6(p,:) = sequencesG1(randomnumbG1(X2(l),2),6);
            PreScanMetG1.press7(p,:) = sequencesG1(randomnumbG1(X2(l),2),7);
            PreScanMetG1.press8(p,:) = sequencesG1(randomnumbG1(X2(l),2),8);
            PreScanMetG1.press9(p,:) = sequencesG1(randomnumbG1(X2(l),2),9);

            if randomnumbG1(X2(l),2) <7
                PreScanMetG1.seqNumb(p,:) = SeqNumb(randomnumbG1(X2(l),2));
                PreScanMetG1.seqType(p,:) = 1;  % trained
            else
                PreScanMetG1.seqNumb(p,:) = SeqRandNumb(randomnumbG1(X2(l),2)-6);
                PreScanMetG1.seqType(p,:) = 4; % random
            end 
            
            z(l) = regexprep(cellstr(num2str(sequencesG1(randomnumbG1(X2(l),2),:))),'\s','');
            PreScanMetG1.cueP(p,:) = repmat(z(l),1,1);

            PreScanMetG2.press1(p,:) = sequencesG2(randomnumbG2(X2(l),2),1);
            PreScanMetG2.press2(p,:) = sequencesG2(randomnumbG2(X2(l),2),2);
            PreScanMetG2.press3(p,:) = sequencesG2(randomnumbG2(X2(l),2),3);
            PreScanMetG2.press4(p,:) = sequencesG2(randomnumbG2(X2(l),2),4);
            PreScanMetG2.press5(p,:) = sequencesG2(randomnumbG2(X2(l),2),5);
            PreScanMetG2.press6(p,:) = sequencesG2(randomnumbG2(X2(l),2),6);
            PreScanMetG2.press7(p,:) = sequencesG2(randomnumbG2(X2(l),2),7);
            PreScanMetG2.press8(p,:) = sequencesG2(randomnumbG2(X2(l),2),8);
            PreScanMetG2.press9(p,:) = sequencesG2(randomnumbG2(X2(l),2),9);

            if randomnumbG2(X2(l),2) <7
                %if trained sequencesG2
                PreScanMetG2.seqNumb(p,:) = SeqNumb(randomnumbG2(X2(l),2));
                PreScanMetG2.seqType(p,:) = 1;
            else
                %random
                PreScanMetG2.seqNumb(p,:) = SeqRandNumb(randomnumbG2(X2(l),2)-6);
                PreScanMetG2.seqType(p,:) = 4;
            end

            z(l) = regexprep(cellstr(num2str(sequencesG2(randomnumbG2(X2(l),2),:))),'\s','');
            PreScanMetG2.cueP(p,:) = repmat(z(l),1,1);

            p = p+1;
            l = l+1;

        end

        dirname = fullfile(baseDir,'PreScanMet');
        dircheck(dirname);
        dsave(fullfile(baseDir,['PreScanMet/','sml_PSM_G1_B', num2str(e), '.tgt']),orderfields(PreScanMetG1,OrderFields));
        dsave(fullfile(baseDir,['PreScanMet/','sml_PSM_G2_B', num2str(e), '.tgt']),orderfields(PreScanMetG2,OrderFields));

    end
   
    case 'ScanSeq'  

    % 8 functional runs per session (3 sessions)
    % 12 seq - 6 trained, 6 untrained
    % each sequence repeated twice - FoSEx
    % 3 trial pairs of each seq (6 in total) per run
    % 12 x 6 - 72 trials per run
    % 6.833 min per run - 54 minutes of scanning

    runs = 8;
 
    ScanSess=3;
    
    sequencesG1 = [SeqG1;SeqG2];
    sequencesG2 = [SeqG2;SeqG1];

    seqforRandG1 = [SeqG1;SeqG1;SeqG1;SeqG2;SeqG2;SeqG2];
    seqforRandG2 = [SeqG2;SeqG2;SeqG2;SeqG1;SeqG1;SeqG1];

    seqG1 = kron(seqforRandG1,ones(2,1));
    seqG2 = kron(seqforRandG2,ones(2,1));

    randomnumbG1(:,1) = [1:36]';
    randomnumbG1(:,2) = [1:6,1:6,1:6,7:12,7:12,7:12]';

    randomnumbG1(:,1) = [1:36]';
    randomnumbG2(:,2) = [1:6,1:6,1:6,7:12,7:12,7:12]';

    % Group 1
    ScanG1.cueC = cellstr(repmat('£',length(seqG1),1));
    ScanG1.cueS = cellstr(repmat('£',length(seqG1),1));

    ScanG1.iti(1:length(seqG1),:) = 0;
    ScanG1.sounds(1:length(seqG1),:) = 1;

    ScanG1.blockType(1:length(seqG1),:) = 3;

    ScanG1.isScan(1:length(seqG1),:) = 1;
    ScanG1.isMetronome(1:length(seqG1),:) = 1;

    ScanG1.FoSEX(1:2:length(seqG1),:) = 1;
    ScanG1.FoSEX(2:2:length(seqG1),:) = 2;

    ScanG1.lastTrial(1:length(seqG1),:) = 0;
    ScanG1.lastTrial(length(seqG1),:) = 1;

    ScanG1.announceTime(1:length(seqG1),:) = 1000;
    ScanG1.feedbackTime(1:length(seqG1),:) = 400;


    % Group 2
    ScanG2.cueC = cellstr(repmat('£',length(seqG2),1));
    ScanG2.cueS = cellstr(repmat('£',length(seqG2),1));

    ScanG2.iti(1:length(seqG2),:) = 0;
    ScanG2.sounds(1:length(seqG2),:) = 1;

    ScanG2.blockType(1:length(seqG2),:) = 3;

    ScanG2.isScan(1:length(seqG2),:) = 1;
    ScanG2.isMetronome(1:length(seqG2),:) = 1;

    ScanG2.FoSEX(1:2:length(seqG2),:) = 1;
    ScanG2.FoSEX(2:2:length(seqG2),:) = 2;

    ScanG2.lastTrial(1:length(seqG2),:) = 0;
    ScanG2.lastTrial(length(seqG2),:) = 1;

    ScanG2.announceTime(1:length(seqG2),:) = 1000;
    ScanG2.feedbackTime(1:length(seqG2),:) = 400;

    for sess = 1:ScanSess
    for e = 1:runs

        X1 = [sample_wor(1:size(randomnumbG1,1),size(randomnumbG1,1))];
        X2 = kron(X1,ones(2,1));

        trialLength = 5000;
        startRun = 5000;

        RandBreaks = sort(sample_wor(1:length(X2),5));
        while (sum(find(RandBreaks == 1)) || sum(rem(RandBreaks,2))<5) % breaks cannot be on first trial, between two consecutive repetitions (only odd trials)
            RandBreaks = sort(sample_wor(1:length(X2),5));
        end

        for i = 1:length(X2)

            if any(ismember(RandBreaks,i))
                startRun = startRun+10000;   % pause of 10000 - for baseline
            end
            ScanG1.startTime(i,:) = startRun;
            ScanG2.startTime(i,:) = startRun;

            startRun = startRun+trialLength;
        end

        p = 1;
        l = 1;

        for i = 1:size(X2,1)

            % Group1
            ScanG1.press1(p,:) = sequencesG1(randomnumbG1(X2(l),2),1);
            ScanG1.press2(p,:) = sequencesG1(randomnumbG1(X2(l),2),2);
            ScanG1.press3(p,:) = sequencesG1(randomnumbG1(X2(l),2),3);
            ScanG1.press4(p,:) = sequencesG1(randomnumbG1(X2(l),2),4);
            ScanG1.press5(p,:) = sequencesG1(randomnumbG1(X2(l),2),5);
            ScanG1.press6(p,:) = sequencesG1(randomnumbG1(X2(l),2),6);
            ScanG1.press7(p,:) = sequencesG1(randomnumbG1(X2(l),2),7);
            ScanG1.press8(p,:) = sequencesG1(randomnumbG1(X2(l),2),8);
            ScanG1.press9(p,:) = sequencesG1(randomnumbG1(X2(l),2),9);

            ScanG1.seqNumb(p,:) = randomnumbG1(X2(l),2);
            ScanG1.hand(p,:) = 2;

            if ScanG1.seqNumb(p,:) <7
                %if trained sequencesG1
                ScanG1.seqType(p,:) = 1;
            elseif ScanG1.seqNumb(p,:) >6
                %untrained sequences (seqG2)
                ScanG1.seqType(p,:) = 2;
            end

            ScanG1.ScanSess(p,:) = sess;    

            z(l) = regexprep(cellstr(num2str(sequencesG1(randomnumbG1(X2(l),2),:))),'\s','');
            ScanG1.cueP(p,:) = repmat(z(l),1,1);

            % Group 2
            ScanG2.press1(p,:) = sequencesG2(randomnumbG2(X2(l),2),1);
            ScanG2.press2(p,:) = sequencesG2(randomnumbG2(X2(l),2),2);
            ScanG2.press3(p,:) = sequencesG2(randomnumbG2(X2(l),2),3);
            ScanG2.press4(p,:) = sequencesG2(randomnumbG2(X2(l),2),4);
            ScanG2.press5(p,:) = sequencesG2(randomnumbG2(X2(l),2),5);
            ScanG2.press6(p,:) = sequencesG2(randomnumbG2(X2(l),2),6);
            ScanG2.press7(p,:) = sequencesG2(randomnumbG2(X2(l),2),7);
            ScanG2.press8(p,:) = sequencesG2(randomnumbG2(X2(l),2),8);
            ScanG2.press9(p,:) = sequencesG2(randomnumbG2(X2(l),2),9);

            ScanG2.seqNumb(p,:) = randomnumbG2(X2(l),2);
            ScanG2.hand(p,:) = 2;

            if ScanG2.seqNumb(p,:) <7
                %if trained sequencesG2
                ScanG2.seqType(p,:) = 1;
            elseif ScanG2.seqNumb(p,:) >6
                %untrained sequences (seqG1)
                ScanG2.seqType(p,:) = 2;
            end
            
            ScanG2.ScanSess(p,:) = sess;        
            
            z(l) = regexprep(cellstr(num2str(sequencesG2(randomnumbG2(X2(l),2),:))),'\s','');
            ScanG2.cueP(p,:) = repmat(z(l),1,1);

            p = p+1;
            l = l+1;

        end

        dirname = fullfile(baseDir,'Scan');
        dircheck(dirname);

        dsave(fullfile(baseDir,['Scan/','sml_Scan_G1_B', num2str((sess-1)*8+e), '.tgt']),orderfields(ScanG1,OrderFields));
        dsave(fullfile(baseDir,['Scan/','sml_Scan_G2_B', num2str((sess-1)*8+e), '.tgt']),orderfields(ScanG2,OrderFields));

    end
    end

    case 'ScanLoc'  

    % finger mapping localiser
    % 2 functional runs
    % 4 repetitions of each finger per run
    % 2.167 min per run - 4.3 in total

    ScanSess=4;
    
    seq = repmat(FingMap,4,1);
    runs = 2;   % 2 localisers

    randomnumb(:,1) = [1:20]';
    randomnumb(:,2) = repmat([1:5]',4,1);

    Scan.cueC = cellstr(repmat('£',length(seq),1));
    Scan.cueS = cellstr(repmat('£',length(seq),1));

    Scan.iti(1:length(seq),:) = 0;
    Scan.sounds(1:length(seq),:) = 1;

    Scan.blockType(1:length(seq),:) = 4;

    Scan.isScan(1:length(seq),:) = 1;
    Scan.isMetronome(1:length(seq),:) = 1;

    Scan.FoSEX(1:length(seq),:) = 1;

    Scan.lastTrial(1:length(seq),:) = 0;
    Scan.lastTrial(length(seq),:) = 1;

    Scan.announceTime(1:length(seq),:) = 1000;
    Scan.feedbackTime(1:length(seq),:) = 400;

    for sess=1:ScanSess
        for e = 1:runs
            
            X2 = [sample_wor(1:size(FingMap,1),size(FingMap,1));...
                sample_wor(1:size(FingMap,1),size(FingMap,1));...
                sample_wor(1:size(FingMap,1),size(FingMap,1));...
                sample_wor(1:size(FingMap,1),size(FingMap,1))];
            
            trialLength = 5000;
            startRun = 5000;
            
            RandBreaks = sort(sample_wor(1:length(X2),3));
            while (find(RandBreaks == 1)) % breaks cannot be on the first trial
                RandBreaks = sort(sample_wor(1:length(X2),3));
            end
            
            for i = 1:length(X2)
                
                if any(ismember(RandBreaks,i))
                    startRun = startRun+10000;   % pause of 10000 - for estimating baseline
                end
                Scan.startTime(i,:) = startRun;
                
                startRun = startRun+trialLength;
            end
            
            p = 1;
            l = 1;
            
            for i = 1:size(X2,1)
                
                Scan.press1(p,:) = seq(randomnumb(X2(l),2),1);
                Scan.press2(p,:) = seq(randomnumb(X2(l),2),2);
                Scan.press3(p,:) = seq(randomnumb(X2(l),2),3);
                Scan.press4(p,:) = seq(randomnumb(X2(l),2),4);
                Scan.press5(p,:) = seq(randomnumb(X2(l),2),5);
                Scan.press6(p,:) = seq(randomnumb(X2(l),2),6);
                Scan.press7(p,:) = seq(randomnumb(X2(l),2),7);
                Scan.press8(p,:) = seq(randomnumb(X2(l),2),8);
                Scan.press9(p,:) = seq(randomnumb(X2(l),2),9);
                
                Scan.seqNumb(p,:) = SeqNumbFingMap(randomnumb(X2(l)));
                Scan.hand(p,:) = 2;
                Scan.seqType(p,:) = 3;  % finger mapping
                Scan.ScanSess(p,:) = sess;   
                
                z(l) = regexprep(cellstr(num2str(seq(randomnumb(X2(l),2),:))),'\s','');
                
                Scan.cueP(p,:) = repmat(z(l),1,1);
                
                p = p+1;
                l = l+1;
                
            end

        dirname = fullfile(baseDir,'Scan');
        dircheck(dirname);
        dsave(fullfile(baseDir,['Scan/','sml_Scan_FingMap', num2str((sess-1)*2+e),'.tgt']),orderfields(Scan,OrderFields));
        end
    end

    case 'PreScanGo'

    % Practice for ScanGo session (no metronome, just GO signal)    

    sequencesG1 = [SeqG1;RandomSeq(1:6,:)];
    sequencesG2 = [SeqG2;RandomSeq(1:6,:)];

    runs=1;

    seqforRandG1 = [sequencesG1;sequencesG1;sequencesG1];
    seqforRandG2 = [sequencesG2;sequencesG2;sequencesG2];

    seqG1 = kron(seqforRandG1,ones(2,1));
    seqG2 = kron(seqforRandG2,ones(2,1));

    randomnumbG1(:,1) = [1:36]';
    randomnumbG1(:,2) = [1:6,1:6,1:6,7:12,7:12,7:12]';

    randomnumbG1(:,1) = [1:36]';
    randomnumbG2(:,2) = [1:6,1:6,1:6,7:12,7:12,7:12]';

    % Group 1
    PreScanGoG1.cueC = cellstr(repmat('£',length(seqG1),1));
    PreScanGoG1.cueS = cellstr(repmat('£',length(seqG1),1));

    PreScanGoG1.hand(1:length(seqG1),:) = 2;

    PreScanGoG1.iti(1:length(seqG1),:) = 0;
    PreScanGoG1.sounds(1:length(seqG1),:) = 1;

    PreScanGoG1.blockType(1:length(seqG1),:) = 8;

    PreScanGoG1.isScan(1:length(seqG1),:) = 2;
    PreScanGoG1.ScanSess(1:length(seqG2),:) = 0;
    PreScanGoG1.isMetronome(1:length(seqG1),:) = 2;

    PreScanGoG1.FoSEX(1:2:length(seqG1),:) = 1;
    PreScanGoG1.FoSEX(2:2:length(seqG1),:) = 2;

    PreScanGoG1.lastTrial(1:length(seqG1),:) = 0;
    PreScanGoG1.lastTrial(length(seqG1),:) = 1;

    PreScanGoG1.announceTime(1:length(seqG1),:) = 1000;
    PreScanGoG1.feedbackTime(1:length(seqG1),:) = 400;

    % Group 2
    PreScanGoG2.cueC = cellstr(repmat('£',length(seqG2),1));
    PreScanGoG2.cueS = cellstr(repmat('£',length(seqG2),1));

    PreScanGoG2.hand(1:length(seqG2),:) = 2;

    PreScanGoG2.iti(1:length(seqG2),:) = 0;
    PreScanGoG2.sounds(1:length(seqG2),:) = 1;

    PreScanGoG2.blockType(1:length(seqG2),:) = 8;

    PreScanGoG2.isScan(1:length(seqG2),:) = 2;
    PreScanGoG2.ScanSess(1:length(seqG2),:) = 0;
    PreScanGoG2.isMetronome(1:length(seqG2),:) = 2;

    PreScanGoG2.FoSEX(1:2:length(seqG2),:) = 1;
    PreScanGoG2.FoSEX(2:2:length(seqG2),:) = 2;

    PreScanGoG2.lastTrial(1:length(seqG2),:) = 0;
    PreScanGoG2.lastTrial(length(seqG2),:) = 1;

    PreScanGoG2.announceTime(1:length(seqG2),:) = 1000;
    PreScanGoG2.feedbackTime(1:length(seqG2),:) = 400;

    for e = 1:runs

        X1 = [sample_wor(1:size(randomnumbG1,1),size(randomnumbG1,1))];
        X2 = kron(X1,ones(2,1));    % for all trials - with repetitions

        trialLength = 5000;
        startRun = 5000;

        RandBreaks = sort(sample_wor(1:length(X2),5));
        while (sum(find(RandBreaks == 1)) || sum(rem(RandBreaks,2))<5) % breaks cannot be on first trial, between two consecutive repetitions (only odd trials)
            RandBreaks = sort(sample_wor(1:length(X2),5));
        end


        for i = 1:length(X2)

            if any(ismember(RandBreaks,i))
                startRun = startRun+10000;   % pause of 8000 - for baseline
            end
            PreScanGoG1.startTime(i,:) = startRun;
            PreScanGoG2.startTime(i,:) = startRun;

            startRun = startRun+trialLength;
        end

        p = 1;
        l = 1;

        for i = 1:length(X2)
            PreScanGoG1.press1(p,:) = sequencesG1(randomnumbG1(X2(l),2),1);
            PreScanGoG1.press2(p,:) = sequencesG1(randomnumbG1(X2(l),2),2);
            PreScanGoG1.press3(p,:) = sequencesG1(randomnumbG1(X2(l),2),3);
            PreScanGoG1.press4(p,:) = sequencesG1(randomnumbG1(X2(l),2),4);
            PreScanGoG1.press5(p,:) = sequencesG1(randomnumbG1(X2(l),2),5);
            PreScanGoG1.press6(p,:) = sequencesG1(randomnumbG1(X2(l),2),6);
            PreScanGoG1.press7(p,:) = sequencesG1(randomnumbG1(X2(l),2),7);
            PreScanGoG1.press8(p,:) = sequencesG1(randomnumbG1(X2(l),2),8);
            PreScanGoG1.press9(p,:) = sequencesG1(randomnumbG1(X2(l),2),9);

            if randomnumbG1(X2(l),2) <7
                PreScanGoG1.seqNumb(p,:) = SeqNumb(randomnumbG1(X2(l),2));
                PreScanGoG1.seqType(p,:) = 1;  % trained
            else
                PreScanGoG1.seqNumb(p,:) = SeqRandNumb(randomnumbG1(X2(l),2)-6);
                PreScanGoG1.seqType(p,:) = 4; % random
            end

            z(l) = regexprep(cellstr(num2str(sequencesG1(randomnumbG1(X2(l),2),:))),'\s','');
            PreScanGoG1.cueP(p,:) = repmat(z(l),1,1);

            PreScanGoG2.press1(p,:) = sequencesG2(randomnumbG2(X2(l),2),1);
            PreScanGoG2.press2(p,:) = sequencesG2(randomnumbG2(X2(l),2),2);
            PreScanGoG2.press3(p,:) = sequencesG2(randomnumbG2(X2(l),2),3);
            PreScanGoG2.press4(p,:) = sequencesG2(randomnumbG2(X2(l),2),4);
            PreScanGoG2.press5(p,:) = sequencesG2(randomnumbG2(X2(l),2),5);
            PreScanGoG2.press6(p,:) = sequencesG2(randomnumbG2(X2(l),2),6);
            PreScanGoG2.press7(p,:) = sequencesG2(randomnumbG2(X2(l),2),7);
            PreScanGoG2.press8(p,:) = sequencesG2(randomnumbG2(X2(l),2),8);
            PreScanGoG2.press9(p,:) = sequencesG2(randomnumbG2(X2(l),2),9);

            if randomnumbG2(X2(l),2) <7
                %if trained sequencesG2
                PreScanGoG2.seqNumb(p,:) = SeqNumb(randomnumbG2(X2(l),2));
                PreScanGoG2.seqType(p,:) = 1;
            else
                %random
                PreScanGoG2.seqNumb(p,:) = SeqRandNumb(randomnumbG2(X2(l),2)-6);
                PreScanGoG2.seqType(p,:) = 4;
            end
            
            z(l) = regexprep(cellstr(num2str(sequencesG2(randomnumbG2(X2(l),2),:))),'\s','');
            PreScanGoG2.cueP(p,:) = repmat(z(l),1,1);

            p = p+1;
            l = l+1;

        end

        dirname = fullfile(baseDir,'PreScanGo');
        dircheck(dirname);
        dsave(fullfile(baseDir,['PreScanGo/','sml_PSGo_G1_B', num2str(e), '.tgt']),orderfields(PreScanGoG1,OrderFields));
        dsave(fullfile(baseDir,['PreScanGo/','sml_PSGo_G2_B', num2str(e), '.tgt']),orderfields(PreScanGoG2,OrderFields));

    end
 
    case 'ScanGo'

    % NO METRONOME - perform as fast as possible
    % 8 functional runs
    % 12 seq - 6 trained, 6 untrained
    % each sequence repeated twice - FoSEx
    % 3 trial pairs of each seq (6 in total) per run
    % 12 x 6 - 72 trials per run
    % 6.67 min per run - 54 minutes of scanning
    
    runs = 8;

    sequencesG1 = [SeqG1;SeqG2];
    sequencesG2 = [SeqG2;SeqG1];

    seqforRandG1 = [SeqG1;SeqG1;SeqG1;SeqG2;SeqG2;SeqG2];
    seqforRandG2 = [SeqG2;SeqG2;SeqG2;SeqG1;SeqG1;SeqG1];

    seqG1 = kron(seqforRandG1,ones(2,1));
    seqG2 = kron(seqforRandG2,ones(2,1));

    randomnumbG1(:,1) = [1:36]';
    randomnumbG1(:,2) = [1:6,1:6,1:6,7:12,7:12,7:12]';

    randomnumbG1(:,1) = [1:36]';
    randomnumbG2(:,2) = [1:6,1:6,1:6,7:12,7:12,7:12]';

    % Group 1
    ScanGoG1.cueC = cellstr(repmat('£',length(seqG1),1));
    ScanGoG1.cueS = cellstr(repmat('£',length(seqG1),1));

    ScanGoG1.iti(1:length(seqG1),:) = 0;
    ScanGoG1.sounds(1:length(seqG1),:) = 1;

    ScanGoG1.blockType(1:length(seqG1),:) = 9;

    ScanGoG1.isScan(1:length(seqG1),:) = 1;
    ScanGoG1.isMetronome(1:length(seqG1),:) = 2;

    ScanGoG1.FoSEX(1:2:length(seqG1),:) = 1;
    ScanGoG1.FoSEX(2:2:length(seqG1),:) = 2;

    ScanGoG1.lastTrial(1:length(seqG1),:) = 0;
    ScanGoG1.lastTrial(length(seqG1),:) = 1;

    ScanGoG1.announceTime(1:length(seqG1),:) = 1000;
    ScanGoG1.feedbackTime(1:length(seqG1),:) = 400;

    ScanGoG1.ScanSess(1:length(seqG1),:) = 4;

    % Group 2
    ScanGoG2.cueC = cellstr(repmat('£',length(seqG2),1));
    ScanGoG2.cueS = cellstr(repmat('£',length(seqG2),1));

    ScanGoG2.iti(1:length(seqG2),:) = 0;
    ScanGoG2.sounds(1:length(seqG2),:) = 1;

    ScanGoG2.blockType(1:length(seqG2),:) = 9;

    ScanGoG2.isScan(1:length(seqG2),:) = 1;
    ScanGoG2.isMetronome(1:length(seqG2),:) = 2;

    ScanGoG2.FoSEX(1:2:length(seqG2),:) = 1;
    ScanGoG2.FoSEX(2:2:length(seqG2),:) = 2;

    ScanGoG2.lastTrial(1:length(seqG2),:) = 0;
    ScanGoG2.lastTrial(length(seqG2),:) = 1;

    ScanGoG2.announceTime(1:length(seqG2),:) = 1000;
    ScanGoG2.feedbackTime(1:length(seqG2),:) = 400;

    ScanGoG2.ScanSess(1:length(seqG2),:) = 4;

    for e = 1:runs

        X1 = [sample_wor(1:size(randomnumbG1,1),size(randomnumbG1,1))];
        X2 = kron(X1,ones(2,1));

        trialLength = 5000;
        startRun = 5000;

        RandBreaks = sort(sample_wor(1:length(X2),5));
        while (sum(find(RandBreaks == 1)) || sum(rem(RandBreaks,2))<5) % breaks cannot be on first trial, between two consecutive repetitions (only odd trials)
            RandBreaks = sort(sample_wor(1:length(X2),5));
        end

        for i = 1:length(X2)

            if any(ismember(RandBreaks,i))
                startRun = startRun+10000;   % pause of 8000 - for baseline
            end
            ScanGoG1.startTime(i,:) = startRun;
            ScanGoG2.startTime(i,:) = startRun;

            startRun = startRun+trialLength;
        end

        p = 1;
        l = 1;

        for i = 1:size(X2,1)

            % Group1
            ScanGoG1.press1(p,:) = sequencesG1(randomnumbG1(X2(l),2),1);
            ScanGoG1.press2(p,:) = sequencesG1(randomnumbG1(X2(l),2),2);
            ScanGoG1.press3(p,:) = sequencesG1(randomnumbG1(X2(l),2),3);
            ScanGoG1.press4(p,:) = sequencesG1(randomnumbG1(X2(l),2),4);
            ScanGoG1.press5(p,:) = sequencesG1(randomnumbG1(X2(l),2),5);
            ScanGoG1.press6(p,:) = sequencesG1(randomnumbG1(X2(l),2),6);
            ScanGoG1.press7(p,:) = sequencesG1(randomnumbG1(X2(l),2),7);
            ScanGoG1.press8(p,:) = sequencesG1(randomnumbG1(X2(l),2),8);
            ScanGoG1.press9(p,:) = sequencesG1(randomnumbG1(X2(l),2),9);

            ScanGoG1.seqNumb(p,:) = randomnumbG1(X2(l),2);
            ScanGoG1.hand(p,:) = 2;

            if ScanGoG1.seqNumb(p,:) <7
                %if trained sequencesG1
                ScanGoG1.seqType(p,:) = 1;
            elseif ScanGoG1.seqNumb(p,:) >6
                %untrained sequences (seqG2)
                ScanGoG1.seqType(p,:) = 2;
            end

            z(l) = regexprep(cellstr(num2str(sequencesG1(randomnumbG1(X2(l),2),:))),'\s','');
            ScanGoG1.cueP(p,:) = repmat(z(l),1,1);

            % Group 2
            ScanGoG2.press1(p,:) = sequencesG2(randomnumbG2(X2(l),2),1);
            ScanGoG2.press2(p,:) = sequencesG2(randomnumbG2(X2(l),2),2);
            ScanGoG2.press3(p,:) = sequencesG2(randomnumbG2(X2(l),2),3);
            ScanGoG2.press4(p,:) = sequencesG2(randomnumbG2(X2(l),2),4);
            ScanGoG2.press5(p,:) = sequencesG2(randomnumbG2(X2(l),2),5);
            ScanGoG2.press6(p,:) = sequencesG2(randomnumbG2(X2(l),2),6);
            ScanGoG2.press7(p,:) = sequencesG2(randomnumbG2(X2(l),2),7);
            ScanGoG2.press8(p,:) = sequencesG2(randomnumbG2(X2(l),2),8);
            ScanGoG2.press9(p,:) = sequencesG2(randomnumbG2(X2(l),2),9);

            ScanGoG2.seqNumb(p,:) = randomnumbG2(X2(l),2);
            ScanGoG2.hand(p,:) = 2;

            if ScanGoG2.seqNumb(p,:) <7
                %if trained sequencesG2
                ScanGoG2.seqType(p,:) = 1;
            elseif ScanGoG2.seqNumb(p,:) >6
                %untrained sequences (seqG1)
                ScanGoG2.seqType(p,:) = 2;
            end
            
            z(l) = regexprep(cellstr(num2str(sequencesG2(randomnumbG2(X2(l),2),:))),'\s','');
            ScanGoG2.cueP(p,:) = repmat(z(l),1,1);

            p = p+1;
            l = l+1;

        end

        dirname = fullfile(baseDir,'Scan');
        dircheck(dirname);
        dsave(fullfile(baseDir,['Scan/','sml_ScanGo_G1_B', num2str(e), '.tgt']),orderfields(ScanGoG1,OrderFields));
        dsave(fullfile(baseDir,['Scan/','sml_ScanGo_G2_B', num2str(e), '.tgt']),orderfields(ScanGoG2,OrderFields));

    end

    case 'OtherHandRand'

    % Test for the transfer to the other hand - baseline (before learning)
    % Intrinsic and extrinsic versions of random sequences
    % Each sequence repeated twice
    % 2 extrinsic, 2 intrinsic sequences per block 
    % 1 FoSEx pair per block (extr., intr.)

    Trials = 1:24;       

    runs = 2;

    for e = 1:runs

        OtherHandRand_G2.cueC = cellstr(repmat('£',length(Trials),1));
        OtherHandRand_G1.cueC = cellstr(repmat('£',length(Trials),1));

        OtherHandRand_G1.iti(1:length(Trials),:) = 500;
        OtherHandRand_G2.iti(1:length(Trials),:) = 500;

        OtherHandRand_G1.hand(1:length(Trials),:) = 1;
        OtherHandRand_G2.hand(1:length(Trials),:) = 1;

        OtherHandRand_G1.blockType(1:length(Trials),:) = 10;
        OtherHandRand_G2.blockType(1:length(Trials),:) = 10;
        
        OtherHandRand_G1.announceTime(1:length(Trials),:) = 0;
        OtherHandRand_G2.announceTime(1:length(Trials),:) = 0;
        OtherHandRand_G1.feedbackTime(1:length(Trials),:) = 400;
        OtherHandRand_G2.feedbackTime(1:length(Trials),:) = 400;

        OtherHandRand_G1.sounds(1:length(Trials),:) = 1;
        OtherHandRand_G2.sounds(1:length(Trials),:) = 1;

        OtherHandRand_G1.cueS = cellstr(repmat('£',length(Trials),1));  % £
        OtherHandRand_G2.cueS = cellstr(repmat('£',length(Trials),1));
        OtherHandRand_G1.cueC = cellstr(repmat('£',length(Trials),1)); % £
        OtherHandRand_G2.cueC = cellstr(repmat('£',length(Trials),1));

        OtherHandRand_G1.isScan = zeros(length(Trials),1);
        OtherHandRand_G2.isScan = zeros(length(Trials),1);

        OtherHandRand_G1.ScanSess = zeros(length(Trials),1);
        OtherHandRand_G2.ScanSess = zeros(length(Trials),1);

        OtherHandRand_G1.FoSEX(1:2:length(Trials),:) = 1;
        OtherHandRand_G1.FoSEX(2:2:length(Trials),:) = 2;

        OtherHandRand_G2.FoSEX(1:2:length(Trials),:) = 1;
        OtherHandRand_G2.FoSEX(2:2:length(Trials),:) = 2;

        OtherHandRand_G1.isMetronome(1:length(Trials),:) = 0;
        OtherHandRand_G2.isMetronome(1:length(Trials),:) = 0;

        OtherHandRand_G1.lastTrial(1:length(Trials),1) = 0;
        OtherHandRand_G1.lastTrial(length(Trials),:) = 1;
        OtherHandRand_G2.lastTrial(1:length(Trials),1) = 0;
        OtherHandRand_G2.lastTrial(length(Trials),:) = 1;

        OtherHandRand_G1.startTime = zeros(length(Trials),1);
        OtherHandRand_G2.startTime = zeros(length(Trials),1);

        SeqG1 = [RandomSeq(1:6,:);6-RandomSeq(1:6,:)]; % intrinsic(1:6), extrinsic(7:12)
        SeqG2 = [RandomSeq(1:6,:);6-RandomSeq(1:6,:)];

        X = [sample_wor(1:12,12);sample_wor(1:12,12)];

        p = 1;
        l = 1;
        for i = 1:length(Trials)/2

            OtherHandRand_G1.cueP(p:p+1,:) = regexprep(cellstr(num2str(6-SeqG1(X(l),:))),'\s','');  % 6 - Seq - left hand
            OtherHandRand_G2.cueP(p:p+1,:) = regexprep(cellstr(num2str(6-SeqG2(X(l),:))),'\s','');

            OtherHandRand_G1.press1(p:p+1,:) = SeqG1(X(l),1);
            OtherHandRand_G1.press2(p:p+1,:) = SeqG1(X(l),2);
            OtherHandRand_G1.press3(p:p+1,:) = SeqG1(X(l),3);
            OtherHandRand_G1.press4(p:p+1,:) = SeqG1(X(l),4);
            OtherHandRand_G1.press5(p:p+1,:) = SeqG1(X(l),5);
            OtherHandRand_G1.press6(p:p+1,:) = SeqG1(X(l),6);
            OtherHandRand_G1.press7(p:p+1,:) = SeqG1(X(l),7);
            OtherHandRand_G1.press8(p:p+1,:) = SeqG1(X(l),8);
            OtherHandRand_G1.press9(p:p+1,:) = SeqG1(X(l),9);

            OtherHandRand_G2.press1(p:p+1,:) = SeqG2(X(l),1);
            OtherHandRand_G2.press2(p:p+1,:) = SeqG2(X(l),2);
            OtherHandRand_G2.press3(p:p+1,:) = SeqG2(X(l),3);
            OtherHandRand_G2.press4(p:p+1,:) = SeqG2(X(l),4);
            OtherHandRand_G2.press5(p:p+1,:) = SeqG2(X(l),5);
            OtherHandRand_G2.press6(p:p+1,:) = SeqG2(X(l),6);
            OtherHandRand_G2.press7(p:p+1,:) = SeqG2(X(l),7);
            OtherHandRand_G2.press8(p:p+1,:) = SeqG2(X(l),8);
            OtherHandRand_G2.press9(p:p+1,:) = SeqG2(X(l),9);

            if X(l) <7  % intrinsic
                OtherHandRand_G1.seqNumb(p:p+1,:) = SeqNumbRandOtherHand(X(l));
                OtherHandRand_G1.seqType(p:p+1,:) = 8;
                OtherHandRand_G2.seqNumb(p:p+1,:) = SeqNumbRandOtherHand(X(l));
                OtherHandRand_G2.seqType(p:p+1,:) = 8;
            else        % extrinsic
                OtherHandRand_G1.seqNumb(p:p+1,:) = SeqNumbRandOtherHand(X(l)-6);
                OtherHandRand_G1.seqType(p:p+1,:) = 9;
                OtherHandRand_G2.seqNumb(p:p+1,:) = SeqNumbRandOtherHand(X(l)-6);
                OtherHandRand_G2.seqType(p:p+1,:) = 9;
            end

            p = p+2;
            l = l+1;

        end

        dirname = fullfile(baseDir,'OtherHandRand');
        dircheck(dirname);
        dsave(fullfile(baseDir,['OtherHandRand/','sml_OHR_G1_B' ,num2str(e), '.tgt']),orderfields(OtherHandRand_G1,OrderFields));
        dsave(fullfile(baseDir,['OtherHandRand/','sml_OHR_G2_B' ,num2str(e), '.tgt']),orderfields(OtherHandRand_G2,OrderFields));

    end
        
    case 'SeqTest'              % varargin - test_sess(1-4) (e.g. 'SeqTest',1)
        % Testing performance on trained seq, generalisation
        % - Trained sequences
        % - Random sequences (no chunk in common)
        % - Same chunks in sequences, order randomised
        % 4 sessions, always 3 runs
        % Each sequence presented twice in a row
        % per session - 4 rep per sequence (2 doubles FoSEx)
        
        test_sess = varargin{1};    % which test - session - different stimuli
        
        runs = 3;
        trial_run = 24;
        
        Trials = 1:(runs*trial_run);  % 3 blocks of 24 - 72
        
        switch (test_sess)
            case 1
                SeqG1 = [SeqG1; RandomSeq(13:18,:); RandChunkSeq1(1:6,:)];
                SeqG2 = [SeqG2; RandomSeq(13:18,:); RandChunkSeq2(1:6,:)];
                RandomSeqNumb = SeqRandNumb(13:18);
                RandomChunkNumb = SeqRandChunkNum(1:6);
            case 2
                SeqG1 = [SeqG1; RandomSeq(19:24,:); RandChunkSeq1(7:12,:)];
                SeqG2 = [SeqG2; RandomSeq(19:24,:); RandChunkSeq2(7:12,:)];
                RandomSeqNumb = SeqRandNumb(19:24);
                RandomChunkNumb = SeqRandChunkNum(7:12);
            case 3
                SeqG1 = [SeqG1; RandomSeq(25:30,:); RandChunkSeq1(13:18,:)];
                SeqG2 = [SeqG2; RandomSeq(25:30,:); RandChunkSeq2(13:18,:)];
                RandomSeqNumb = SeqRandNumb(25:30);
                RandomChunkNumb = SeqRandChunkNum(13:18);
            case 4
                SeqG1 = [SeqG1; RandomSeq(31:36,:); RandChunkSeq1(19:24,:)];
                SeqG2 = [SeqG2; RandomSeq(31:36,:); RandChunkSeq2(19:24,:)];
                RandomSeqNumb = SeqRandNumb(31:36);
                RandomChunkNumb = SeqRandChunkNum(19:24);
        end

        
        X = [sample_wor(1:size(SeqG1,1),size(SeqG1,1));sample_wor(1:size(SeqG1,1),size(SeqG1,1))];
        
        SeqTest_G1.cueC = cellstr(repmat('£',length(Trials),1));
        SeqTest_G2.cueC = cellstr(repmat('£',length(Trials),1));
        
        SeqTest_G1.iti(1:length(Trials),:) = 500;
        SeqTest_G2.iti(1:length(Trials),:) = 500;
        
        SeqTest_G1.hand(1:length(Trials),:) = 2;
        SeqTest_G2.hand(1:length(Trials),:) = 2;
        
        SeqTest_G1.blockType(1:length(Trials),:) = 11;
        SeqTest_G2.blockType(1:length(Trials),:) = 11;
        
        SeqTest_G1.announceTime(1:length(Trials),:) = 0;
        SeqTest_G2.announceTime(1:length(Trials),:) = 0;
        SeqTest_G1.feedbackTime(1:length(Trials),:) = 400;
        SeqTest_G2.feedbackTime(1:length(Trials),:) = 400;
        
        SeqTest_G1.sounds(1:length(Trials),:) = 1;
        SeqTest_G2.sounds(1:length(Trials),:) = 1;
        
        SeqTest_G1.isScan(1:length(Trials),:) = 0;
        SeqTest_G2.isScan(1:length(Trials),:) = 0;
        
        SeqTest_G1.ScanSess(1:length(Trials),:) = 0;
        SeqTest_G2.ScanSess(1:length(Trials),:) = 0;
        
        SeqTest_G1.isMetronome(1:length(Trials),:) = 0;
        SeqTest_G2.isMetronome(1:length(Trials),:) = 0;
        
        p = 1;
        l = 1;
        
        for i = 1:(length(Trials)/2)
            
            SeqTest_G1.cueP(p:p+1,:) = regexprep(cellstr(num2str(SeqG1(X(l),:))),'\s','');
            SeqTest_G2.cueP(p:p+1,:) = regexprep(cellstr(num2str(SeqG2(X(l),:))),'\s','');
            SeqTest_G1.cueS(p:p+1,:) = C1(1);  %£
            SeqTest_G2.cueS(p:p+1,:) = C1(1);
            
            SeqTest_G1.press1(p:p+1,:) = SeqG1(X(l),1);
            SeqTest_G1.press2(p:p+1,:) = SeqG1(X(l),2);
            SeqTest_G1.press3(p:p+1,:) = SeqG1(X(l),3);
            SeqTest_G1.press4(p:p+1,:) = SeqG1(X(l),4);
            SeqTest_G1.press5(p:p+1,:) = SeqG1(X(l),5);
            SeqTest_G1.press6(p:p+1,:) = SeqG1(X(l),6);
            SeqTest_G1.press7(p:p+1,:) = SeqG1(X(l),7);
            SeqTest_G1.press8(p:p+1,:) = SeqG1(X(l),8);
            SeqTest_G1.press9(p:p+1,:) = SeqG1(X(l),9);
            
            SeqTest_G2.press1(p:p+1,:) = SeqG2(X(l),1);
            SeqTest_G2.press2(p:p+1,:) = SeqG2(X(l),2);
            SeqTest_G2.press3(p:p+1,:) = SeqG2(X(l),3);
            SeqTest_G2.press4(p:p+1,:) = SeqG2(X(l),4);
            SeqTest_G2.press5(p:p+1,:) = SeqG2(X(l),5);
            SeqTest_G2.press6(p:p+1,:) = SeqG2(X(l),6);
            SeqTest_G2.press7(p:p+1,:) = SeqG2(X(l),7);
            SeqTest_G2.press8(p:p+1,:) = SeqG2(X(l),8);
            SeqTest_G2.press9(p:p+1,:) = SeqG2(X(l),9);
            
             if X(l) <7
                SeqTest_G1.seqNumb(p:p+1,:) = SeqNumb(X(l));
                SeqTest_G1.seqType(p:p+1,:) = 1;  % trained
                SeqTest_G2.seqNumb(p:p+1,:) = SeqNumb(X(l));
                SeqTest_G2.seqType(p:p+1,:) = 1;  % trained
             elseif (X(l)>6 & X(l)<13)
                SeqTest_G1.seqNumb(p:p+1,:) = RandomSeqNumb(X(l)-6);
                SeqTest_G1.seqType(p:p+1,:) = 4;  % random seq
                SeqTest_G2.seqNumb(p:p+1,:) = RandomSeqNumb(X(l)-6);
                SeqTest_G2.seqType(p:p+1,:) = 4;  % random seq 
             else
                SeqTest_G1.seqNumb(p:p+1,:) = RandomChunkNumb(X(l)-12);
                SeqTest_G1.seqType(p:p+1,:) = 5;  % random chunks
                SeqTest_G2.seqNumb(p:p+1,:) = RandomChunkNumb(X(l)-12);
                SeqTest_G2.seqType(p:p+1,:) = 5;  % random chunks
             end
            
            p = p+2;
            l = l+1;
            
        end
        
        SeqTest_G1.isScan = zeros(length(Trials),1);
        SeqTest_G2.isScan = zeros(length(Trials),1);

        SeqTest_G1.FoSEX = repmat([1:2]',length(Trials)/2,1);
        SeqTest_G2.FoSEX = repmat([1:2]',length(Trials)/2,1);
        
        SeqTest_G1.lastTrial(1:length(Trials),1) = 0;
        SeqTest_G1.lastTrial([24,48,72],:) = 1;      % every 24th trial - last trial in a block
        SeqTest_G2.lastTrial(1:length(Trials),1) = 0;
        SeqTest_G2.lastTrial([24,48,72],:) = 1;
        
        SeqTest_G1.startTime = zeros(length(Trials),1);
        SeqTest_G2.startTime = zeros(length(Trials),1);
        
        B1 = 1;
        B2 = trial_run;
        
        dirname = fullfile(baseDir,'SeqTest');
        dircheck(dirname);
        
        for i = 1:runs
            dsave(fullfile(baseDir,['SeqTest/','sml_STB_G1_B', num2str((test_sess-1)*runs+i), '.tgt']),orderfields(getrow(SeqTest_G1,B1:B2),OrderFields));
            dsave(fullfile(baseDir,['SeqTest/','sml_STB_G2_B', num2str((test_sess-1)*runs+i), '.tgt']),orderfields(getrow(SeqTest_G2,B1:B2),OrderFields));
            B1 = B1+trial_run;
            B2 = B2+trial_run;
            
        end
          
    case 'OtherHandTest'        % varargin - test_sess(1-4) (e.g. 'OtherHandTest',1)
    
    % Tests of generalisation to the other hand
    % - Trained sequences (intr., extr.)
    % - Random sequences (intr., extr.)
    % 4 sessions, always 4 runs
    % Each sequence presented twice in a row
    % per session - 4 rep per sequence (2 doubles FoSEx)

    test_sess = varargin{1};    % which test - session - different stimuli
        
        runs = 4;
        trial_run = 24;
        
        Trials = 1:(runs*trial_run);  % 4 blocks of 24 - 96
        
        switch (test_sess)
            case 1
                SeqG1 = [SeqG1; 6-SeqG1; RandomSeq(13:18,:); 6-RandomSeq(13:18,:)];
                SeqG2 = [SeqG2; 6-SeqG2; RandomSeq(13:18,:); 6-RandomSeq(13:18,:)];
                RandomSeqNumb = SeqNumbRandOtherHand(7:12);
            case 2
                SeqG1 = [SeqG1; 6-SeqG1; RandomSeq(19:24,:); 6-RandomSeq(19:24,:)];
                SeqG2 = [SeqG2; 6-SeqG2; RandomSeq(19:24,:); 6-RandomSeq(19:24,:)];
                RandomSeqNumb = SeqNumbRandOtherHand(13:18);
            case 3
                SeqG1 = [SeqG1; 6-SeqG1; RandomSeq(25:30,:); 6-RandomSeq(25:30,:)];
                SeqG2 = [SeqG2; 6-SeqG2; RandomSeq(25:30,:); 6-RandomSeq(25:30,:)];
                RandomSeqNumb = SeqNumbRandOtherHand(19:24);
            case 4
                SeqG1 = [SeqG1; 6-SeqG1; RandomSeq(31:36,:); 6-RandomSeq(31:36,:)];
                SeqG2 = [SeqG2; 6-SeqG2; RandomSeq(31:36,:); 6-RandomSeq(31:36,:)];
                RandomSeqNumb = SeqNumbRandOtherHand(25:30);
        end

        
        X = [sample_wor(1:size(SeqG1,1),size(SeqG1,1));sample_wor(1:size(SeqG1,1),size(SeqG1,1))];
        
        OtherHandTest_G1.cueC = cellstr(repmat('£',length(Trials),1));
        OtherHandTest_G2.cueC = cellstr(repmat('£',length(Trials),1));
        
        OtherHandTest_G1.iti(1:length(Trials),:) = 500;
        OtherHandTest_G2.iti(1:length(Trials),:) = 500;
        
        OtherHandTest_G1.hand(1:length(Trials),:) = 1;
        OtherHandTest_G2.hand(1:length(Trials),:) = 1;
        
        OtherHandTest_G1.blockType(1:length(Trials),:) = 12;
        OtherHandTest_G2.blockType(1:length(Trials),:) = 12;
        
        OtherHandTest_G1.announceTime(1:length(Trials),:) = 0;
        OtherHandTest_G2.announceTime(1:length(Trials),:) = 0;
        OtherHandTest_G1.feedbackTime(1:length(Trials),:) = 400;
        OtherHandTest_G2.feedbackTime(1:length(Trials),:) = 400;
        
        OtherHandTest_G1.sounds(1:length(Trials),:) = 1;
        OtherHandTest_G2.sounds(1:length(Trials),:) = 1;
        
        OtherHandTest_G1.isScan(1:length(Trials),:) = 0;
        OtherHandTest_G2.isScan(1:length(Trials),:) = 0;
        
        OtherHandTest_G1.ScanSess(1:length(Trials),:) = 0;
        OtherHandTest_G2.ScanSess(1:length(Trials),:) = 0;
        
        OtherHandTest_G1.isMetronome(1:length(Trials),:) = 0;
        OtherHandTest_G2.isMetronome(1:length(Trials),:) = 0;
        
        p = 1;
        l = 1;
        
        for i = 1:(length(Trials)/2)
            
            OtherHandTest_G1.cueP(p:p+1,:) = regexprep(cellstr(num2str(6-SeqG1(X(l),:))),'\s',''); 
            OtherHandTest_G2.cueP(p:p+1,:) = regexprep(cellstr(num2str(6-SeqG2(X(l),:))),'\s',''); 
            OtherHandTest_G1.cueS(p:p+1,:) = C1(1);  %£
            OtherHandTest_G2.cueS(p:p+1,:) = C1(1);
            
            OtherHandTest_G1.press1(p:p+1,:) = SeqG1(X(l),1);
            OtherHandTest_G1.press2(p:p+1,:) = SeqG1(X(l),2);
            OtherHandTest_G1.press3(p:p+1,:) = SeqG1(X(l),3);
            OtherHandTest_G1.press4(p:p+1,:) = SeqG1(X(l),4);
            OtherHandTest_G1.press5(p:p+1,:) = SeqG1(X(l),5);
            OtherHandTest_G1.press6(p:p+1,:) = SeqG1(X(l),6);
            OtherHandTest_G1.press7(p:p+1,:) = SeqG1(X(l),7);
            OtherHandTest_G1.press8(p:p+1,:) = SeqG1(X(l),8);
            OtherHandTest_G1.press9(p:p+1,:) = SeqG1(X(l),9);
            
            OtherHandTest_G2.press1(p:p+1,:) = SeqG2(X(l),1);
            OtherHandTest_G2.press2(p:p+1,:) = SeqG2(X(l),2);
            OtherHandTest_G2.press3(p:p+1,:) = SeqG2(X(l),3);
            OtherHandTest_G2.press4(p:p+1,:) = SeqG2(X(l),4);
            OtherHandTest_G2.press5(p:p+1,:) = SeqG2(X(l),5);
            OtherHandTest_G2.press6(p:p+1,:) = SeqG2(X(l),6);
            OtherHandTest_G2.press7(p:p+1,:) = SeqG2(X(l),7);
            OtherHandTest_G2.press8(p:p+1,:) = SeqG2(X(l),8);
            OtherHandTest_G2.press9(p:p+1,:) = SeqG2(X(l),9);
            
             if X(l) <7                     % intrinsic train
                OtherHandTest_G1.seqNumb(p:p+1,:) = SeqNumbTrainOtherHand(X(l));
                OtherHandTest_G1.seqType(p:p+1,:) = 6;  
                OtherHandTest_G2.seqNumb(p:p+1,:) = SeqNumbTrainOtherHand(X(l));
                OtherHandTest_G2.seqType(p:p+1,:) = 6;  
             elseif (X(l)>6 & X(l)<13)      % extrinsic train
                OtherHandTest_G1.seqNumb(p:p+1,:) = SeqNumbTrainOtherHand(X(l)-6);
                OtherHandTest_G1.seqType(p:p+1,:) = 7; 
                OtherHandTest_G2.seqNumb(p:p+1,:) = SeqNumbTrainOtherHand(X(l)-6);
                OtherHandTest_G2.seqType(p:p+1,:) = 7;  
             elseif (X(l)>12 & X(l)<19)     % intrinsic random
                OtherHandTest_G1.seqNumb(p:p+1,:) = RandomSeqNumb(X(l)-12);
                OtherHandTest_G1.seqType(p:p+1,:) = 8;  % random chunks
                OtherHandTest_G2.seqNumb(p:p+1,:) = RandomSeqNumb(X(l)-12);
                OtherHandTest_G2.seqType(p:p+1,:) = 8;  % random chunks
             else                           % extrinsic random
                OtherHandTest_G1.seqNumb(p:p+1,:) = RandomSeqNumb(X(l)-18);
                OtherHandTest_G1.seqType(p:p+1,:) = 9;  % random chunks
                OtherHandTest_G2.seqNumb(p:p+1,:) = RandomSeqNumb(X(l)-18);
                OtherHandTest_G2.seqType(p:p+1,:) = 9;  % random chunks
             end
            
            p = p+2;
            l = l+1;
            
        end
        
        OtherHandTest_G1.isScan = zeros(length(Trials),1);
        OtherHandTest_G2.isScan = zeros(length(Trials),1);

        OtherHandTest_G1.FoSEX = repmat([1:2]',length(Trials)/2,1);
        OtherHandTest_G2.FoSEX = repmat([1:2]',length(Trials)/2,1);
        
        OtherHandTest_G1.lastTrial(1:length(Trials),1) = 0;
        OtherHandTest_G1.lastTrial([24,48,72,96],:) = 1;      % every 24th trial - last trial in a block
        OtherHandTest_G2.lastTrial(1:length(Trials),1) = 0;
        OtherHandTest_G2.lastTrial([24,48,72,96],:) = 1;
        
        OtherHandTest_G1.startTime = zeros(length(Trials),1);
        OtherHandTest_G2.startTime = zeros(length(Trials),1);
        
        B1 = 1;
        B2 = trial_run;
        
        dirname = fullfile(baseDir,'OtherHandTest');
        dircheck(dirname);
        
        for i = 1:runs
            dsave(fullfile(baseDir,['OtherHandTest/','sml_OHT_G1_B', num2str((test_sess-1)*runs+i), '.tgt']),orderfields(getrow(OtherHandTest_G1,B1:B2),OrderFields));
            dsave(fullfile(baseDir,['OtherHandTest/','sml_OHT_G2_B', num2str((test_sess-1)*runs+i), '.tgt']),orderfields(getrow(OtherHandTest_G2,B1:B2),OrderFields));
            B1 = B1+trial_run;
            B2 = B2+trial_run;
            
        end

    case 'BeforeLearnFingMap'  

    % finger mapping localiser - before scanner
    % 1 functional run
    % 4 repetitions of each finger per run
    % 2.167 min per run
    
    seq = repmat(FingMap,4,1);
    runs = 1;
    
    randomnumb(:,1) = [1:20]';
    randomnumb(:,2) = repmat([1:5]',4,1);
    
    Scan.cueC = cellstr(repmat('£',length(seq),1));
    Scan.cueS = cellstr(repmat('£',length(seq),1));
    
    Scan.iti(1:length(seq),:) = 0;
    Scan.sounds(1:length(seq),:) = 1;
    
    Scan.blockType(1:length(seq),:) = 13;
    
    Scan.isScan(1:length(seq),:) = 2;
    Scan.isMetronome(1:length(seq),:) = 1;
    
    Scan.FoSEX(1:length(seq),:) = 1;
    
    Scan.lastTrial(1:length(seq),:) = 0;
    Scan.lastTrial(length(seq),:) = 1;
    
    Scan.announceTime(1:length(seq),:) = 1000;
    Scan.feedbackTime(1:length(seq),:) = 400;  
    
    for e = 1:runs
        
        X2 = [sample_wor(1:size(FingMap,1),size(FingMap,1));...
            sample_wor(1:size(FingMap,1),size(FingMap,1));...
            sample_wor(1:size(FingMap,1),size(FingMap,1));...
            sample_wor(1:size(FingMap,1),size(FingMap,1))];
        
        trialLength = 5000;
        startRun = 5000;
        
        RandBreaks = sort(sample_wor(1:length(X2),3));
        while (find(RandBreaks == 1)) % breaks cannot be on the first trial
            RandBreaks = sort(sample_wor(1:length(X2),3));
        end
        
        for i = 1:length(X2)
            
            if any(ismember(RandBreaks,i))
                startRun = startRun+10000;   % pause of 10000 - for estimating baseline
            end
            Scan.startTime(i,:) = startRun;
            
            startRun = startRun+trialLength;
        end
        
        p = 1;
        l = 1;
        
        for i = 1:size(X2,1)
            
            Scan.press1(p,:) = seq(randomnumb(X2(l),2),1);
            Scan.press2(p,:) = seq(randomnumb(X2(l),2),2);
            Scan.press3(p,:) = seq(randomnumb(X2(l),2),3);
            Scan.press4(p,:) = seq(randomnumb(X2(l),2),4);
            Scan.press5(p,:) = seq(randomnumb(X2(l),2),5);
            Scan.press6(p,:) = seq(randomnumb(X2(l),2),6);
            Scan.press7(p,:) = seq(randomnumb(X2(l),2),7);
            Scan.press8(p,:) = seq(randomnumb(X2(l),2),8);
            Scan.press9(p,:) = seq(randomnumb(X2(l),2),9);
            
            Scan.seqNumb(p,:) = SeqNumbFingMap(randomnumb(X2(l)));
            Scan.hand(p,:) = 2;
            Scan.seqType(p,:) = 3;  % finger mapping
            Scan.ScanSess(p,:) = 0;
            
            z(l) = regexprep(cellstr(num2str(seq(randomnumb(X2(l),2),:))),'\s','');
            
            Scan.cueP(p,:) = repmat(z(l),1,1);
            
            p = p+1;
            l = l+1;
            
        end
        
        dirname = fullfile(baseDir,'BeforeLearnFingMap');
        dircheck(dirname);
        dsave(fullfile(baseDir,['BeforeLearnFingMap/','sml_BLFM_B', num2str(e),'.tgt']),orderfields(Scan,OrderFields));
    end

    otherwise
        error('no such case!')
end


% Local functions
    function dircheck(dir)
        if ~exist(dir,'dir');
            warning('%s doesn''t exist. Creating one now. You''re welcome! \n',dir);
            mkdir(dir);
        end
    end

end





