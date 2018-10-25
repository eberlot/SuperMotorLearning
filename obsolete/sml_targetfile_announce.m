function varargout=sml_targetfile(what,varargin)
% make target files

%%%%%% groups %%%%%%
% G1 - training with right hand (hand=2)
% G2 - training with left hand  (hand=1)

%%%%%% stimuli %%%%%%%
% chunks - triplets, sequences - 9 presses
% training - 6 sequences
% seqType - 1 training, 2 , 3 , 4 random

%%%%%% scan sessions %%%%%%%
% train sequences
% test seq - 4 (train from the other group)
% train other hand (4)

switch what
    case 'targetfile'
        basedir = '/Users/eberlot/Documents/MATLAB/Projects/SuperMotorLearning/targetFiles'; % change!
        %basedir = 'C:\Nicola\SequenceHierarchical\data\sh4\analysis';
        cd(basedir);
        
        rng('shuffle');
        
        Seq1 = [3 2 4 5 3 2 5 4 3; ...
            5 2 1 4 1 2 3 2 4; ...
            5 4 3 1 2 3 4 1 2; ...
            3 4 2 1 3 4 1 2 3; ...
            1 4 5 2 5 4 3 4 2; ...
            1 2 3 5 4 3 2 5 4];
        
        Seq2 = [1 4 3 5 4 2 1 2 4; ...
            5 2 3 1 2 4 5 4 2; ...
            1 2 4 3 5 2 4 5 2; ...
            5 4 2 3 1 4 2 1 4; ...
            3 1 4 5 2 3 4 3 2; ...
            3 5 2 1 4 3 2 3 4];
        
        
        ChunkAmount = 1:10;
        SequenceLength = 9;
        ChunkNumb = [101 102 103 104 105 106 107 108 109 110];
        
        
        Chunks1 = [ 3 2 4; 5 3 2; 5 4 3; 5 2 1; 4 1 2; 1 2 3; 3 4 2; 1 3 4; 1 4 5; 2 5 4];
        Chunks2 = [ 1 2 4; 1 4 3; 3 1 4; 4 3 2; 5 2 3; 3 5 2; 4 5 2; 5 4 2; 2 3 4; 2 1 4];
        
        
        SeqChunks1 = [1,2,3; ...
            4,5,1; ...
            3,6,5; ...
            7,8,6; ...
            9,10,7; ...
            6,3,10];
        
        SeqChunks2 = [2,8,1; ...
            5,1,8; ...
            1,6,7; ...
            8,3,10; ...
            3,5,4; ...
            6,2,9];
        
        UnSeq1 = [1 4 3 5 4 2 1 2 4; ...
            5 2 3 1 2 4 5 4 3; ...
            1 2 4 3 5 2 4 5 2; ...
            5 4 2 3 1 4 2 1 4];
        
        UnSeq2 = [5 2 1 4 1 2 3 2 4; ...
            1 4 5 2 5 4 3 4 2; ...
            5 4 3 1 2 3 4 1 2; ...
            1 2 3 5 4 3 2 5 4];
        
        OtherHand1 = [1 3 2 5 3 4 5 2 4; ...
            5 3 4 1 3 2 1 4 2; ...
            1 3 2 4 3 2 3 5 1; ...
            5 3 4 2 3 4 3 1 5];
        
        OtherHand2 = [5 2 1 4 1 2 3 2 4; ...
            1 4 5 2 5 4 3 4 2; ...
            5 4 3 1 2 3 4 1 2; ...
            1 2 3 5 4 3 2 5 4];
        
        FingMap = [ 1 1 1 1 1 1 1 1 1; ...
            2 2 2 2 2 2 2 2 2; ...
            3 3 3 3 3 3 3 3 3; ...
            4 4 4 4 4 4 4 4 4; ...
            5 5 5 5 5 5 5 5 5];
        
        RandomSeq = [4 1 2 4 2 3 5 2 1; ...
            4 3 1 3 2 4 5 2 4; ...
            3 5 1 4 3 5 2 5 2; ...
            3 1 2 3 2 3 5 4 2; ...
            1 5 1 3 1 5 3 5 2; ...
            3 1 4 3 5 2 5 1 4; ...
            1 3 1 4 2 3 2 5 1; ...
            2 5 4 5 3 1 2 4 5; ...
            4 1 4 5 2 3 2 1 3; ...
            4 1 2 4 1 5 3 1 4; ...
            3 1 5 4 2 4 5 3 2; ...
            3 5 1 3 2 1 5 1 3; ...
            4 5 4 1 2 5 4 1 2; ...
            5 1 5 1 4 1 2 3 2; ...
            3 4 2 5 4 2 3 4 2; ...
            4 1 4 2 5 2 3 5 2; ...
            1 2 3 5 2 3 4 1 3; ...
            2 4 1 5 1 5 4 1 2; ...
            3 4 3 1 4 3 5 4 2; ...
            5 3 5 3 4 2 5 4 1; ...
            4 2 1 3 2 1 4 3 1; ...
            2 1 5 3 1 4 2 1 5; ...
            5 1 5 4 3 5 2 3 2; ...
            1 3 2 1 5 2 3 2 3; ...
            5 3 1 5 3 2 1 5 4; ...
            4 5 3 1 4 3 4 2 5; ...
            3 1 5 4 1 2 1 4 3; ...
            4 3 1 5 3 4 5 4 2];
        
        
        C1 = {'£'; '£';'£';'£';'£';'£'};
        C2 = {'£'; '£'; '£';'£'; '£';'£';'£';'£';'£';'£'};
        C3 = {'***'};
        C3Seq = {'*********'};
        
        SequenceNumb = 1:6;
        
        OrderFields = {'startTime','FT','seqNumb','seqType','press1','press2','press3','press4','press5','press6','press7','press8','press9','hand','cueS','cueC','cueP','iti','sounds','isScan','isMetronome','lastTrial','FoSEX','announceTime','feedbackTime','feedbackType'};
        
        GroupHand = 1;      % or 2 - define here!!! always change!!!
        
        switch GroupHand
            case 1
                SeqG1 = Seq1;
                SeqG2 = Seq2;
                UnSeqG1 = UnSeq1;
                UnSeqG2 = UnSeq2;
                SeqChunksG1 = SeqChunks1;
                SeqChunksG2 = SeqChunks2;
                ChunksG1 = Chunks1;
                ChunksG2 = Chunks2;
                OtherHandG1 = OtherHand1;
                OtherHandG2 = OtherHand2;
            case 2
                SeqG1 = Seq2;
                SeqG2 = Seq1;
                UnSeqG1 = UnSeq2;
                UnSeqG2 = UnSeq1;
                SeqChunksG1 = SeqChunks2;
                SeqChunksG2 = SeqChunks1;
                ChunksG1 = Chunks2;
                ChunksG2 = Chunks1;
                OtherHandG1 = OtherHand2;
                OtherHandG2 = OtherHand1;
        end
        
        %% ChunkTrain Block
        Trials = 1:110;
        
        for e = 1:34 %64 - 34(Eva test)
            ChunkTrain_G1.cueS = cellstr(repmat('£',length(Trials),1));
            ChunkTrain_G2.cueS = cellstr(repmat('£',length(Trials),1));
            
            ChunkTrain_G1.iti(1:length(Trials),:) = 500;
            ChunkTrain_G2.iti(1:length(Trials),:) = 500;
            ChunkTrain_G1.FT(1:length(Trials),:) = 1;
            ChunkTrain_G2.FT(1:length(Trials),:) = 1;
            ChunkTrain_G1.hand(1:length(Trials),:) = 2;
            ChunkTrain_G2.hand(1:length(Trials),:) = 1;
            
            ChunkTrain_G1.announceTime(1:length(Trials),:) = 0;
            ChunkTrain_G2.announceTime(1:length(Trials),:) = 0;
            ChunkTrain_G1.feedbackTime(1:length(Trials),:) = 350;
            ChunkTrain_G2.feedbackTime(1:length(Trials),:) = 350;
            ChunkTrain_G1.feedbackType(1:length(Trials),:) = 1;
            ChunkTrain_G2.feedbackType(1:length(Trials),:) = 1;
            
            ChunkTrain_G1.sounds(1:length(Trials),:) = 1;
            ChunkTrain_G2.sounds(1:length(Trials),:) = 1;
            
            %randomly choose order
            X = sample_wor(ChunkAmount,length(ChunkAmount));
            p = 1;
            l = 1;
            for i = 1:(length(Trials)/2)
                
                % press cue
                ChunkTrain_G1.cueP(p+1,:) = C3;
                ChunkTrain_G1.cueP(p,:) = regexprep(cellstr(num2str(ChunksG1(X(l),:))),'\s','');
                
                ChunkTrain_G2.cueP(p+1,:) = C3;
                ChunkTrain_G2.cueP(p,:) = regexprep(cellstr(num2str(6-ChunksG2(X(l),:))),'\s','');
                
                
                ChunkTrain_G1.press1(p:p+1,:) = ChunksG1(X(l),1);
                ChunkTrain_G1.press2(p:p+1,:) = ChunksG1(X(l),2);
                ChunkTrain_G1.press3(p:p+1,:) = ChunksG1(X(l),3);
                
                ChunkTrain_G2.press1(p:p+1,:) = ChunksG2(X(l),1);
                ChunkTrain_G2.press2(p:p+1,:) = ChunksG2(X(l),2);
                ChunkTrain_G2.press3(p:p+1,:) = ChunksG2(X(l),3);
                
                ChunkTrain_G1.seqNumb(p:p+1,:) = ChunkNumb(X(l));
                ChunkTrain_G2.seqNumb(p:p+1,:) = ChunkNumb(X(l));
                
                
                %chunk cue
                ChunkTrain_G1.cueC(p:p+1,:) = C2(X(l),:) ;
                ChunkTrain_G2.cueC(p:p+1,:) = C2(X(l),:) ;
                
                p = p+2;
                l = l+1;
                if l == 11;
                    l = 1;
                    X = sample_wor(ChunkAmount,length(ChunkAmount));
                end
            end
            ChunkTrain_G1.press4 = zeros(length(Trials),1);
            ChunkTrain_G1.press5 = zeros(length(Trials),1);
            ChunkTrain_G1.press6 = zeros(length(Trials),1);
            ChunkTrain_G1.press7 = zeros(length(Trials),1);
            ChunkTrain_G1.press8 = zeros(length(Trials),1);
            ChunkTrain_G1.press9 = zeros(length(Trials),1);
            
            
            ChunkTrain_G2.press4 = zeros(length(Trials),1);
            ChunkTrain_G2.press5 = zeros(length(Trials),1);
            ChunkTrain_G2.press6 = zeros(length(Trials),1);
            ChunkTrain_G2.press7 = zeros(length(Trials),1);
            ChunkTrain_G2.press8 = zeros(length(Trials),1);
            ChunkTrain_G2.press9 = zeros(length(Trials),1);
            
            
            ChunkTrain_G1.seqType = zeros(length(Trials),1);
            ChunkTrain_G2.seqType = zeros(length(Trials),1);
            
            
            ChunkTrain_G1.isScan = zeros(length(Trials),1);
            ChunkTrain_G2.isScan = zeros(length(Trials),1);
            
            ChunkTrain_G1.seqType = zeros(length(Trials),1);
            ChunkTrain_G2.seqType = zeros(length(Trials),1);
            
            ChunkTrain_G1.FoSEX = repmat([1:2]',length(Trials)/2,1);
            ChunkTrain_G2.FoSEX = repmat([1:2]',length(Trials)/2,1);
            
            ChunkTrain_G1.isMetronome = zeros(length(Trials),1);
            ChunkTrain_G2.isMetronome = zeros(length(Trials),1);
            
            ChunkTrain_G1.lastTrial = zeros(length(Trials),1);
            ChunkTrain_G2.lastTrial = zeros(length(Trials),1);
            
            ChunkTrain_G1.startTime = zeros(length(Trials),1);
            ChunkTrain_G2.startTime = zeros(length(Trials),1);
            
            
            switch GroupHand
                
                case 1
                    
                    dsave(fullfile(basedir,['ChunkLearn/','sml_CLB_G1_R_B', num2str(e), '.tgt']),orderfields(ChunkTrain_G1,OrderFields));
                    dsave(fullfile(basedir,['ChunkLearn/','sml_CLB_G2_L_B', num2str(e), '.tgt']),orderfields(ChunkTrain_G2,OrderFields));
                    
                case 2
                    
                    dsave(fullfile(basedir,['ChunkLearn/','sml_CLB_G2_R_B', num2str(e), '.tgt']),orderfields(ChunkTrain_G1,OrderFields));
                    dsave(fullfile(basedir,['ChunkLearn/','sml_CLB_G1_L_B', num2str(e), '.tgt']),orderfields(ChunkTrain_G2,OrderFields));
                    
            end
            clear X
            clear ChunkTrain_G1
            clear ChunkTrain_G2
        end
        
        
        %% Sequence Learning Block
        
        Trials = 1:24;       % 2 pairs of repetition per sequence per block (FoSEx - 4 in total)
        
        for e = 1:54 % 128 - 54(Eva test)
            
            SequenceLearn_G2.cueC = cellstr(repmat('£',length(Trials),1));
            SequenceLearn_G1.cueC = cellstr(repmat('£',length(Trials),1));
            
            SequenceLearn_G1.iti(1:length(Trials),:) = 500;
            SequenceLearn_G2.iti(1:length(Trials),:) = 500;
            
            SequenceLearn_G1.hand(1:length(Trials),:) = 2;
            SequenceLearn_G2.hand(1:length(Trials),:) = 1;
            
            SequenceLearn_G1.FT(1:length(Trials),:) = 2;
            SequenceLearn_G2.FT(1:length(Trials),:) = 2;
            
            SequenceLearn_G1.announceTime(1:length(Trials),:) = 0;
            SequenceLearn_G2.announceTime(1:length(Trials),:) = 0;
            SequenceLearn_G1.feedbackTime(1:length(Trials),:) = 350;
            SequenceLearn_G2.feedbackTime(1:length(Trials),:) = 350;
            SequenceLearn_G1.feedbackType(1:length(Trials),:) = 1;
            SequenceLearn_G2.feedbackType(1:length(Trials),:) = 1;
            
            SequenceLearn_G1.sounds(1:length(Trials),:) = 1;
            SequenceLearn_G2.sounds(1:length(Trials),:) = 1;
            
            X = sample_wor(SequenceNumb,length(SequenceNumb));
            p = 1;
            l = 1;
            for i = 1:(length(Trials)/2)
                
                SequenceLearn_G1.cueP(p+1,:) = C3Seq;
                SequenceLearn_G1.cueP(p,:) = regexprep(cellstr(num2str(SeqG1(X(l),:))),'\s','');
                SequenceLearn_G2.cueP(p+1,:) = C3Seq;
                SequenceLearn_G2.cueP(p,:) = regexprep(cellstr(num2str(SeqG2(X(l),:))),'\s','');
                SequenceLearn_G1.cueS(p:p+1,:) = C1(X(l));
                SequenceLearn_G2.cueS(p:p+1,:) = C1(X(l));
                SequenceLearn_G1.cueC(p,:) = regexprep(cellstr(strjoin(C2(SeqChunksG1(X(l),:))')),'\s','');
                SequenceLearn_G2.cueC(p,:) = regexprep(cellstr(strjoin(C2(SeqChunksG2(X(l),:))')),'\s','');
                
                
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
                
                SequenceLearn_G1.seqNumb(p:p+1,:) = SequenceNumb(X(l));
                SequenceLearn_G2.seqNumb(p:p+1,:) = SequenceNumb(X(l));
                
                SequenceLearn_G1.isScan = zeros(length(Trials),1);
                SequenceLearn_G2.isScan = zeros(length(Trials),1);
                
                % one for trained sequence
                SequenceLearn_G1.seqType = ones(length(Trials),1);
                SequenceLearn_G2.seqType = ones(length(Trials),1);
                
                SequenceLearn_G1.FoSEX = repmat([1:2]',length(Trials)/2,1);
                SequenceLearn_G2.FoSEX = repmat([1:2]',length(Trials)/2,1);
                
                SequenceLearn_G1.isMetronome = zeros(length(Trials),1);
                SequenceLearn_G2.isMetronome = zeros(length(Trials),1);
                
                SequenceLearn_G1.lastTrial = zeros(length(Trials),1);
                SequenceLearn_G2.lastTrial = zeros(length(Trials),1);
                
                SequenceLearn_G1.startTime = zeros(length(Trials),1);
                SequenceLearn_G2.startTime = zeros(length(Trials),1);
                
                
                p = p+2;
                l = l+1;
                if l == length(SequenceNumb)+1;
                    l = 1;
                    X = sample_wor(SequenceNumb,length(SequenceNumb));
                end
            end
            switch GroupHand
                case 1
                    dsave(fullfile(basedir,['SeqLearn/','sml_SLB_G1_R_B' ,num2str(e), '.tgt']),orderfields(SequenceLearn_G1,OrderFields));
                    dsave(fullfile(basedir,['SeqLearn/','sml_SLB_G2_L_B' ,num2str(e), '.tgt']),orderfields(SequenceLearn_G2,OrderFields));
                case 2
                    dsave(fullfile(basedir,['SeqLearn/','sml_SLB_G2_R_B' ,num2str(e), '.tgt']),orderfields(SequenceLearn_G1,OrderFields));
                    dsave(fullfile(basedir,['SeqLearn/','sml_SLB_G1_L_B' ,num2str(e), '.tgt']),orderfields(SequenceLearn_G2,OrderFields));
                    
            end
            
            clear x
            clear SequenceLearn_G1
            clear SequenceLearn_G2
        end
        
        %% Sequence recall block
        for t = 1:9 %20 - 9(Eva test)
            Trials = 1:24;
            
            SequenceRecall_G1.iti(1:length(Trials),:) = 500;
            SequenceRecall_G2.iti(1:length(Trials),:) = 500;
            
            SequenceRecall_G1.hand(1:length(Trials),:) = 2;
            SequenceRecall_G2.hand(1:length(Trials),:) = 1;
            
            SequenceRecall_G1.FT(1:length(Trials),:) = 7;
            SequenceRecall_G2.FT(1:length(Trials),:) = 7;
            
            SequenceRecall_G1.announceTime(1:length(Trials),:) = 0;
            SequenceRecall_G2.announceTime(1:length(Trials),:) = 0;
            SequenceRecall_G1.feedbackTime(1:length(Trials),:) = 350;
            SequenceRecall_G2.feedbackTime(1:length(Trials),:) = 350;
            SequenceRecall_G1.feedbackType(1:length(Trials),:) = 1;
            SequenceRecall_G2.feedbackType(1:length(Trials),:) = 1;
            
            SequenceRecall_G1.sounds(1:length(Trials),:) = 1;
            SequenceRecall_G2.sounds(1:length(Trials),:) = 1;
            
            SequenceRecall_G1.isScan = zeros(length(Trials),1);
            SequenceRecall_G2.isScan = zeros(length(Trials),1);
            
            % one for trained sequence
            SequenceRecall_G1.seqType = ones(length(Trials),1);
            SequenceRecall_G2.seqType = ones(length(Trials),1);
            
            SequenceRecall_G1.FoSEX = repmat([1:2]',length(Trials)/2,1);
            SequenceRecall_G2.FoSEX = repmat([1:2]',length(Trials)/2,1);
            
            SequenceRecall_G1.isMetronome = zeros(length(Trials),1);
            SequenceRecall_G2.isMetronome = zeros(length(Trials),1);
            
            SequenceRecall_G1.lastTrial = zeros(length(Trials),1);
            SequenceRecall_G2.lastTrial = zeros(length(Trials),1);
            
            SequenceRecall_G1.startTime = zeros(length(Trials),1);
            SequenceRecall_G2.startTime = zeros(length(Trials),1);
            
            SequenceRecall_G1.cueC = cellstr(repmat('£',length(Trials),1));
            SequenceRecall_G2.cueC = cellstr(repmat('£',length(Trials),1));
            
            
            SequenceRecall_G1.iti(1:length(Trials),:) = 500;
            SequenceRecall_G2.iti(1:length(Trials),:) = 500;
            
            SequenceRecall_G1.hand(1:length(Trials),:) = 2;
            SequenceRecall_G2.hand(1:length(Trials),:) = 1;
            
            SequenceRecall_G1.sounds(1:length(Trials),:) = 1;
            SequenceRecall_G2.sounds(1:length(Trials),:) = 1;
            
            X = sample_wor(SequenceNumb,length(SequenceNumb));
            p = 1;
            l = 1;
            for i = 1:(length(Trials))
                
                SequenceRecall_G1.cueS(p,:) = C1(X(l));
                SequenceRecall_G2.cueS(p,:) = C1(X(l));
                
                SequenceRecall_G1.cueP(p,:) = regexprep(cellstr(num2str(SeqG1(X(l),:))),'\s','');
                SequenceRecall_G2.cueP(p,:) = regexprep(cellstr(num2str(SeqG1(X(l),:))),'\s','');
                
                SequenceRecall_G1.press1(p,:) = SeqG1(X(l),1);
                SequenceRecall_G1.press2(p,:) = SeqG1(X(l),2);
                SequenceRecall_G1.press3(p,:) = SeqG1(X(l),3);
                SequenceRecall_G1.press4(p,:) = SeqG1(X(l),4);
                SequenceRecall_G1.press5(p,:) = SeqG1(X(l),5);
                SequenceRecall_G1.press6(p,:) = SeqG1(X(l),6);
                SequenceRecall_G1.press7(p,:) = SeqG1(X(l),7);
                SequenceRecall_G1.press8(p,:) = SeqG1(X(l),8);
                SequenceRecall_G1.press9(p,:) = SeqG1(X(l),9);
                
                SequenceRecall_G2.press1(p,:) = SeqG2(X(l),1);
                SequenceRecall_G2.press2(p,:) = SeqG2(X(l),2);
                SequenceRecall_G2.press3(p,:) = SeqG2(X(l),3);
                SequenceRecall_G2.press4(p,:) = SeqG2(X(l),4);
                SequenceRecall_G2.press5(p,:) = SeqG2(X(l),5);
                SequenceRecall_G2.press6(p,:) = SeqG2(X(l),6);
                SequenceRecall_G2.press7(p,:) = SeqG2(X(l),7);
                SequenceRecall_G2.press8(p,:) = SeqG2(X(l),8);
                SequenceRecall_G2.press9(p,:) = SeqG2(X(l),9);
                
                SequenceRecall_G1.seqNumb(p,:) = SequenceNumb(X(l));
                SequenceRecall_G2.seqNumb(p,:) = SequenceNumb(X(l));
                
                p = p+1;
                l = l+1;
                if l == length(SequenceNumb)+1;
                    l = 1;
                    X = sample_wor(SequenceNumb,length(SequenceNumb));
                end
            end
            
            switch GroupHand
                case 1
                    dsave(fullfile(basedir,['SeqRecall/','sml_SRB_G1_R_B', num2str(t), '.tgt']),orderfields(SequenceRecall_G1,OrderFields));
                    dsave(fullfile(basedir,['SeqRecall/','sml_SRB_G2_L_B', num2str(t), '.tgt']),orderfields(SequenceRecall_G2,OrderFields));
                case 2
                    dsave(fullfile(basedir,['SeqRecall/','sml_SRB_G2_R_B', num2str(t), '.tgt']),orderfields(SequenceRecall_G1,OrderFields));
                    dsave(fullfile(basedir,['SeqRecall/','sml_SRB_G1_L_B', num2str(t), '.tgt']),orderfields(SequenceRecall_G2,OrderFields));
            end
            
            clear X
            clear SequenceRecall_G1
            clear SequenceRecall_G2
        end
        
        %% Test Block Before Learning
        
        % so 50% left and 50% right
        % twice each - per block 36 trials
        % need 30 extra sequences each sequence repeated 6 times
        
        sequencesG1 = [SeqG1;UnSeqG1;RandomSeq;SeqG2(5:end,:)];
        sequencesG2 = [SeqG2;UnSeqG2;RandomSeq;SeqG1(5:end,:)];
        
        X = [sample_wor(1:size(sequencesG1,1),size(sequencesG1,1));sample_wor(1:size(sequencesG1,1),size(sequencesG1,1));sample_wor(1:size(sequencesG1,1),size(sequencesG1,1))];
        X2  = [[X;X],[ones(size(X,1),1);2*ones(size(X,1),1)]];
        
        Y = [sample_wor(1:size(X2,1),size(X2,1))];
        
        SequenceTestG1.cueC = cellstr(repmat('£',length(X2),1));
        
        SequenceTestG1.cueS = cellstr(repmat('£',length(X2),1));
        
        SequenceTestG1.iti(1:length(X2),:) = 500;
        
        SequenceTestG1.FT(1:length(X2),:) = 4;
        
        SequenceTestG1.announceTime(1:length(Trials),:) = 0;
        SequenceTestG2.announceTime(1:length(Trials),:) = 0;
        SequenceTestG1.feedbackTime(1:length(Trials),:) = 350;
        SequenceTestG2.feedbackTime(1:length(Trials),:) = 350;
        SequenceTestG1.feedbackType(1:length(Trials),:) = 1;
        SequenceTestG2.feedbackType(1:length(Trials),:) = 1;
        
        SequenceTestG1.sounds(1:length(X2),:) = 1;
        
        SequenceTestG1.isScan(1:length(X2),:) = 0;
        
        SequenceTestG1.isMetronome(1:length(X2),:) = 0;
        
        SequenceTestG1.FoSEX(1:length(X2),:) = 1;
        SequenceTestG1.lastTrial = zeros(length(X2),1);
        SequenceTestG1.lastTrial = zeros(length(X2),1);
        
        SequenceTestG1.startTime = zeros(length(X2),1);
        
        
        SequenceTestG2.cueC = cellstr(repmat('£',length(X2),1));
        
        SequenceTestG2.cueS = cellstr(repmat('£',length(X2),1));
        
        SequenceTestG2.iti(1:length(X2),:) = 500;
        
        SequenceTestG2.sounds(1:length(X2),:) = 1;
        
        SequenceTestG2.FT(1:length(X2),:) = 4;
        
        SequenceTestG2.isScan(1:length(X2),:) = 0;
        
        SequenceTestG2.isMetronome(1:length(X2),:) = 0;
        
        SequenceTestG2.FoSEX(1:length(X2),:) = 1;
        
        
        SequenceTestG2.lastTrial = zeros(length(X2),1);
        SequenceTestG2.lastTrial = zeros(length(X2),1);
        
        SequenceTestG2.startTime = zeros(length(X2),1);
        
        p = 1;
        l = 1;
        
        for i = 1:size(X2,1)
            
            
            SequenceTestG1.hand(p,:) = X2(Y(l),2);
            
            SequenceTestG1.press1(p,:) = sequencesG1(X2(Y(l),1),1);
            SequenceTestG1.press2(p,:) = sequencesG1(X2(Y(l),1),2);
            SequenceTestG1.press3(p,:) = sequencesG1(X2(Y(l),1),3);
            SequenceTestG1.press4(p,:) = sequencesG1(X2(Y(l),1),4);
            SequenceTestG1.press5(p,:) = sequencesG1(X2(Y(l),1),5);
            SequenceTestG1.press6(p,:) = sequencesG1(X2(Y(l),1),6);
            SequenceTestG1.press7(p,:) = sequencesG1(X2(Y(l),1),7);
            SequenceTestG1.press8(p,:) = sequencesG1(X2(Y(l),1),8);
            SequenceTestG1.press9(p,:) = sequencesG1(X2(Y(l),1),9);
            SequenceTestG1.seqNumb(p,:) = X2(Y(l),1);
            
            if SequenceTestG1.seqNumb(p,:) <7
                %if trained sequencesG1
                SequenceTestG1.seqType(p,:) = 1;
            elseif SequenceTestG1.seqNumb(p,:) >6 && SequenceTestG1.seqNumb(p,:) < 11
                %if other group
                SequenceTestG1.seqType(p,:) = 2;
            else % random other sequencesG1
                SequenceTestG1.seqType(p,:) = 7;
            end
            if SequenceTestG1.hand(p,:) == 2
                
                z(l) = regexprep(cellstr(num2str(sequencesG1(X2(Y(l),1),:))),'\s','');
                
                SequenceTestG1.cueP(p,:) = repmat(z(l),1,1);
                
            elseif SequenceTestG1.hand(p,:) == 1
                
                z(l) = regexprep(cellstr(num2str(6-sequencesG1(X2(Y(l),1),:))),'\s','');
                
                SequenceTestG1.cueP(p,:) = repmat(z(l),1,1);
                
            end
            
            
            SequenceTestG2.hand(p,:) = X2(Y(l),2);
            
            SequenceTestG2.press1(p,:) = sequencesG2(X2(Y(l),1),1);
            SequenceTestG2.press2(p,:) = sequencesG2(X2(Y(l),1),2);
            SequenceTestG2.press3(p,:) = sequencesG2(X2(Y(l),1),3);
            SequenceTestG2.press4(p,:) = sequencesG2(X2(Y(l),1),4);
            SequenceTestG2.press5(p,:) = sequencesG2(X2(Y(l),1),5);
            SequenceTestG2.press6(p,:) = sequencesG2(X2(Y(l),1),6);
            SequenceTestG2.press7(p,:) = sequencesG2(X2(Y(l),1),7);
            SequenceTestG2.press8(p,:) = sequencesG2(X2(Y(l),1),8);
            SequenceTestG2.press9(p,:) = sequencesG2(X2(Y(l),1),9);
            SequenceTestG2.seqNumb(p,:) = X2(Y(l),1);
            
            if SequenceTestG2.seqNumb(p,:) <7
                %if trained sequencesG1
                SequenceTestG2.seqType(p,:) =1;
            elseif SequenceTestG2.seqNumb(p,:) >6 && SequenceTestG2.seqNumb(p,:) < 11
                %if other group
                SequenceTestG2.seqType(p,:) = 2;
            else % random other sequencesG1
                SequenceTestG2.seqType(p,:) = 4;
            end
            if SequenceTestG2.hand(p,:) == 2
                
                z(l) = regexprep(cellstr(num2str(sequencesG1(X2(Y(l),1),:))),'\s','');
                
                SequenceTestG2.cueP(p,:) = repmat(z(l),1,1);
                
            elseif SequenceTestG2.hand(p,:) == 1
                
                z(l) = regexprep(cellstr(num2str(6-sequencesG1(X2(Y(l),1),:))),'\s','');
                
                SequenceTestG2.cueP(p,:) = repmat(z(l),1,1);
                
            end
            
            p = p+1;
            l = l+1;
            
        end
        B1 = 1;
        B2 = 30;
        
        for i = 1:6 % change to 6
            switch GroupHand
                case 1
                    
                    dsave(fullfile(basedir,['BeforeLearn/','sml_BL_G1_R_B', num2str(i), '.tgt']),orderfields(getrow(SequenceTestG1,B1:B2),OrderFields));
                    dsave(fullfile(basedir,['BeforeLearn/','sml_BL_G2_L_B', num2str(i), '.tgt']),orderfields(getrow(SequenceTestG2,B1:B2),OrderFields));
                    B1 = B1+30;
                    B2 = B2+30;
                case 2
                    dsave(fullfile(basedir,['BeforeLearn/','sml_BL_G2_R_B', num2str(i), '.tgt']),orderfields(getrow(SequenceTestG1,B1:B2),OrderFields));
                    dsave(fullfile(basedir,['BeforeLearn/','sml_BL_G1_L_B', num2str(i), '.tgt']),orderfields(getrow(SequenceTestG2,B1:B2),OrderFields));
                    B1 = B1+30;
                    B2 = B2+30;
                    
            end
        end
        
        %% Before Learning Scan/Metronome
        % so 50% left and 50% right
        % twice each - per block 36 trials
        % add other two sequences at the end of  randomseq
        % need 30 extra sequences each sequence repeated 6 times
        
        sequencesG1 = [SeqG1;UnSeqG1;RandomSeq;SeqG2(5:end,:) ];
        sequencesG2 = [SeqG2;UnSeqG2;RandomSeq;SeqG1(5:end,:)];
        
        X = [sample_wor(1:size(sequencesG1,1),size(sequencesG1,1));sample_wor(1:size(sequencesG1,1),size(sequencesG1,1));sample_wor(1:size(sequencesG1,1),size(sequencesG1,1))];
        X2  = [[X;X],[ones(size(X,1),1);2*ones(size(X,1),1)]];
        
        Y = [sample_wor(1:size(X2,1),size(X2,1))];
        
        BLMG1.cueC = cellstr(repmat('£',length(X2),1));
        
        BLMG1.cueS = cellstr(repmat('£',length(X2),1));
        
        BLMG1.iti(1:length(X2),:) = 0;
        
        BLMG1.sounds(1:length(X2),:) = 1;
        
        BLMG1.FT(1:length(X2),:) = 4;
        
        BLMG1.announceTime(1:2:length(Trials),:) = 2000;
        BLMG2.announceTime(1:2:length(Trials),:) = 2000;
        BLMG1.announceTime(2:2:length(Trials),:) = 0;
        BLMG2.announceTime(2:2:length(Trials),:) = 0;
        BLMG1.feedbackTime(1:length(Trials),:) = 350;
        BLMG2.feedbackTime(1:length(Trials),:) = 350;
        BLMG1.feedbackType(1:2:length(Trials),:) = 1;
        BLMG1.feedbackType(2:2:length(Trials),:) = 2;
        BLMG2.feedbackType(1:2:length(Trials),:) = 1;
        BLMG2.feedbackType(2:2:length(Trials),:) = 2;
        
        BLMG1.isScan(1:length(X2),:) = 1;
        
        BLMG1.isMetronome(1:length(X2),:) = 1;
        
        BLMG1.FoSEX(1:length(X2),:) = 1;
        
        BLMG1.lastTrial(1:length(X2),:) = 0;
        
        BLMG1.lastTrial(length(X2),:) = 1;
        
        
        BLMG2.cueC = cellstr(repmat('£',length(X2),1));
        
        BLMG2.cueS = cellstr(repmat('£',length(X2),1));
        
        BLMG2.iti(1:length(X2),:) = 500;
        
        BLMG2.sounds(1:length(X2),:) = 1;
        
        BLMG2.isScan(1:length(X2),:) = 1;
        
        BLMG2.isMetronome(1:length(X2),:) = 1;
        
        BLMG2.FoSEX(1:length(X2),:) = 1;
        
        BLMG2.lastTrial(1:length(X2),:) = 0;
        
        BLMG2.lastTrial(length(X2),:) = 1;
        
        BLMG2.FT(1:length(X2),:) = 4;
        
        
        p = 1;
        l = 1;
        
        for i = 1:size(X2,1)
            
            
            BLMG1.hand(p,:) = X2(Y(l),2);
            
            BLMG1.press1(p,:) = sequencesG1(X2(Y(l),1),1);
            BLMG1.press2(p,:) = sequencesG1(X2(Y(l),1),2);
            BLMG1.press3(p,:) = sequencesG1(X2(Y(l),1),3);
            BLMG1.press4(p,:) = sequencesG1(X2(Y(l),1),4);
            BLMG1.press5(p,:) = sequencesG1(X2(Y(l),1),5);
            BLMG1.press6(p,:) = sequencesG1(X2(Y(l),1),6);
            BLMG1.press7(p,:) = sequencesG1(X2(Y(l),1),7);
            BLMG1.press8(p,:) = sequencesG1(X2(Y(l),1),8);
            BLMG1.press9(p,:) = sequencesG1(X2(Y(l),1),9);
            BLMG1.seqNumb(p,:) = X2(Y(l),1);
            
            if BLMG1.seqNumb(p,:) <7
                %if trained sequencesG1
                BLMG1.seqType(p,:) = 1;
            elseif BLMG1.seqNumb(p,:) >6 && BLMG1.seqNumb(p,:) < 11;
                %if other group
                BLMG1.seqType(p,:) = 2;
            else % random other sequencesG1
                BLMG1.seqType(p,:) = 4;
            end
            if BLMG1.hand(p,:) == 2
                
                z(l) = regexprep(cellstr(num2str(sequencesG1(X2(Y(l),1),:))),'\s','');
                
                BLMG1.cueP(p,:) = repmat(z(l),1,1);
                
            elseif BLMG1.hand(p,:) == 1
                
                z(l) = regexprep(cellstr(num2str(6-sequencesG1(X2(Y(l),1),:))),'\s','');
                
                BLMG1.cueP(p,:) = repmat(z(l),1,1);
                
            end
            
            trialLength = 8700;
            startRun = 5000;
            
            RandBreaks = sort(sample_wor(1:length(X2),5));
            while find(RandBreaks ==1)
                RandBreaks = sort(sample_wor(1:length(X2),5));
            end
            
            for i = 1:length(X2)
                
                if any(ismember(RandBreaks,i))
                    startRun = startRun+trialLength;
                end
                BLMG1.startTime(i,:) = startRun;
                BLMG1.startTime(i,:) = startRun;
                
                startRun = startRun+trialLength;
            end
            
            
            BLMG2.hand(p,:) = X2(Y(l),2);
            
            BLMG2.press1(p,:) = sequencesG2(X2(Y(l),1),1);
            BLMG2.press2(p,:) = sequencesG2(X2(Y(l),1),2);
            BLMG2.press3(p,:) = sequencesG2(X2(Y(l),1),3);
            BLMG2.press4(p,:) = sequencesG2(X2(Y(l),1),4);
            BLMG2.press5(p,:) = sequencesG2(X2(Y(l),1),5);
            BLMG2.press6(p,:) = sequencesG2(X2(Y(l),1),6);
            BLMG2.press7(p,:) = sequencesG2(X2(Y(l),1),7);
            BLMG2.press8(p,:) = sequencesG2(X2(Y(l),1),8);
            BLMG2.press9(p,:) = sequencesG2(X2(Y(l),1),9);
            BLMG2.seqNumb(p,:) = X2(Y(l),1);
            
            if BLMG2.seqNumb(p,:) <7
                %if trained sequencesG1
                BLMG2.seqType(p,:) =1;
            elseif BLMG2.seqNumb(p,:) >6 && BLMG2.seqNumb(p,:) < 11;
                %if other group
                BLMG2.seqType(p,:) = 2;
            else % random other sequencesG1
                BLMG2.seqType(p,:) = 4;
            end
            if BLMG2.hand(p,:) == 2
                
                z(l) = regexprep(cellstr(num2str(sequencesG1(X2(Y(l),1),:))),'\s','');
                
                BLMG2.cueP(p,:) = repmat(z(l),1,1);
                
            elseif BLMG2.hand(p,:) == 1
                
                z(l) = regexprep(cellstr(num2str(6-sequencesG1(X2(Y(l),1),:))),'\s','');
                
                BLMG2.cueP(p,:) = repmat(z(l),1,1);
                
            end
            
            trialLength = 8700;
            startRun = 5000;
            
            RandBreaks = sort(sample_wor(1:length(X2),5));
            while find(RandBreaks ==1)
                RandBreaks = sort(sample_wor(1:length(X2),5));
            end
            
            block_trials=0;
            for block = 1:4
                for i = 1:60;  %number of trials per block (altogether 240 - length X2)
                    
                    if any(ismember(RandBreaks,i+block_trials))
                        startRun = startRun+trialLength;
                    end
                    
                    BLMG1.startTime(i+block_trials,:) = startRun;
                    BLMG2.startTime(i+block_trials,:) = startRun;
                    
                    startRun = startRun+trialLength;
                end
                startRun = 5000;
                block_trials = block_trials+60;
            end
            p = p+1;
            l = l+1;
            
        end
        B1 = 1;
        B2 = 60;
        for i = 1:4
            
            BLMG1B = getrow(BLMG1,B1:B2);
            BLMG1B.lastTrial(length(BLMG1B.lastTrial)) = 1;
            
            BLMG2B = getrow(BLMG2,B1:B2);
            BLMG2B.lastTrial(length(BLMG1B.lastTrial)) = 1;
            switch GroupHand
                case 1
                    dsave(fullfile(basedir,['BeforeLearnMetronome/','sml_BLM_G1_R_B', num2str(i), '.tgt']),orderfields(BLMG1B,OrderFields));
                    dsave(fullfile(basedir,['BeforeLearnMetronome/','sml_BLM_G2_L_B', num2str(i), '.tgt']),orderfields(BLMG2B,OrderFields));
                    B1 = B1+60;
                    B2 = B2+60;
                    
                case 2
                    dsave(fullfile(basedir,['BeforeLearnMetronome/','sml_BLM_G2_R_B', num2str(i), '.tgt']),orderfields(BLMG1B,OrderFields));
                    dsave(fullfile(basedir,['BeforeLearnMetronome/','sml_BLM_G1_L_B', num2str(i), '.tgt']),orderfields(BLMG2B,OrderFields));
                    B1 = B1+60;
                    B2 = B2+60;
                    
            end
        end
        
        %% Test Block Mondays
        
        % 5 Blocks trained sequencesG1 + intermixed untrained sequencesG1 + 4 left hand same
        % rate
        
        % 5 Blocks with Metronome - same rate left and right intermixed
        
        blocks = 5;
        sequencesG1 = [SeqG1;UnSeqG1;OtherHandG1];
        sequencesG2 = [SeqG2;UnSeqG2;OtherHandG2];
        
        
        for e = 1:blocks
            
            
            X2 = [sample_wor(1:size(sequencesG1,1),size(sequencesG1,1)); sample_wor(1:size(sequencesG1,1),size(sequencesG1,1))];
            
            Y = [sample_wor(1:size(X2,1),size(X2,1))];
            
            PreScanTestG1.cueC = cellstr(repmat('£',length(X2),1));
            
            PreScanTestG1.cueS = cellstr(repmat('£',length(X2),1));
            
            PreScanTestG1.iti(1:length(X2),:) = 500;
            
            PreScanTestG1.sounds(1:length(X2),:) = 1;
            
            PreScanTestG1.FT(1:length(X2),:) = 5;
            
            PreScanTestG1.isScan(1:length(X2),:) = 0;
            
            PreScanTestG1.isMetronome(1:length(X2),:) = 0;
            
            PreScanTestG1.FoSEX =  ones(length(X2),1);
            
            PreScanTestG1.startTime = zeros(length(X2),1);
            PreScanTestG1.startTime = zeros(length(X2),1);
            
            PreScanTestG1.lastTrial(1:length(X2),:) = 0;
            
            
            PreScanTestG2.cueC = cellstr(repmat('£',length(X2),1));
            
            PreScanTestG2.cueS = cellstr(repmat('£',length(X2),1));
            
            PreScanTestG2.iti(1:length(X2),:) = 500;
            
            PreScanTestG2.sounds(1:length(X2),:) = 1;
            
            PreScanTestG2.FT(1:length(X2),:) = 5;
            
            PreScanTestG2.isScan(1:length(X2),:) = 0;
            
            PreScanTestG2.isMetronome(1:length(X2),:) = 0;
            
            PreScanTestG2.FoSEX = ones(length(X2),1);
            
            PreScanTestG2.lastTrial(1:length(X2),:) = 0;
            
            
            
            PreScanTestG2.startTime = zeros(length(X2),1);
            PreScanTestG2.startTime = zeros(length(X2),1);
            
            
            p = 1;
            l = 1;
            
            for i = 1:size(X2,1)
                
                
                PreScanTestG1.press1(p,:) = sequencesG1(X2(Y(l),1),1);
                PreScanTestG1.press2(p,:) = sequencesG1(X2(Y(l),1),2);
                PreScanTestG1.press3(p,:) = sequencesG1(X2(Y(l),1),3);
                PreScanTestG1.press4(p,:) = sequencesG1(X2(Y(l),1),4);
                PreScanTestG1.press5(p,:) = sequencesG1(X2(Y(l),1),5);
                PreScanTestG1.press6(p,:) = sequencesG1(X2(Y(l),1),6);
                PreScanTestG1.press7(p,:) = sequencesG1(X2(Y(l),1),7);
                PreScanTestG1.press8(p,:) = sequencesG1(X2(Y(l),1),8);
                PreScanTestG1.press9(p,:) = sequencesG1(X2(Y(l),1),9);
                PreScanTestG1.seqNumb(p,:) = X2(l);
                
                if PreScanTestG1.seqNumb(p,:) <7
                    %if trained sequencesG1
                    PreScanTestG1.seqType(p,:) = 1;
                    PreScanTestG1.hand(p,:) = 2;
                    
                elseif PreScanTestG1.seqNumb(p,:) >6 && PreScanTestG1.seqNumb(p,:) < 11
                    PreScanTestG1.hand(p,:) = 2;
                    
                    %if other group
                    PreScanTestG1.seqType(p,:) = 2;
                else % other hand sequences
                    PreScanTestG1.seqType(p,:) = 3;
                    PreScanTestG1.hand(p,:) = 1;
                    
                end
                
                
                if PreScanTestG1.hand(p,:) == 2
                    
                    z(l) = regexprep(cellstr(num2str(sequencesG1(X2(Y(l),1),:))),'\s','');
                    
                    PreScanTestG1.cueP(p,:) = repmat(z(l),1,1);
                    
                elseif PreScanTestG1.hand(p,:) == 1
                    
                    z(l) = regexprep(cellstr(num2str(6-sequencesG1(X2(Y(l),1),:))),'\s','');
                    
                    PreScanTestG1.cueP(p,:) = repmat(z(l),1,1);
                    
                end
                PreScanTestG2.press1(p,:) = sequencesG2(X2(Y(l),1),1);
                PreScanTestG2.press2(p,:) = sequencesG2(X2(Y(l),1),2);
                PreScanTestG2.press3(p,:) = sequencesG2(X2(Y(l),1),3);
                PreScanTestG2.press4(p,:) = sequencesG2(X2(Y(l),1),4);
                PreScanTestG2.press5(p,:) = sequencesG2(X2(Y(l),1),5);
                PreScanTestG2.press6(p,:) = sequencesG2(X2(Y(l),1),6);
                PreScanTestG2.press7(p,:) = sequencesG2(X2(Y(l),1),7);
                PreScanTestG2.press8(p,:) = sequencesG2(X2(Y(l),1),8);
                PreScanTestG2.press9(p,:) = sequencesG2(X2(Y(l),1),9);
                PreScanTestG2.seqNumb(p,:) = X2(l);
                
                
                if PreScanTestG2.seqNumb(p,:) <7
                    %if trained sequencesG1
                    PreScanTestG2.seqType(p,:) = 1;
                    PreScanTestG2.hand(p,:) = 1;
                    
                elseif PreScanTestG2.seqNumb(p,:) >6 && PreScanTestG2.seqNumb(p,:) < 11
                    PreScanTestG2.hand(p,:) = 1;
                    
                    %if other group
                    PreScanTestG2.seqType(p,:) = 2;
                else % other hand sequences
                    PreScanTestG2.seqType(p,:) = 3;
                    PreScanTestG2.hand(p,:) = 2;
                    
                end
                
                if PreScanTestG2.hand(p,:) == 2
                    
                    z(l) = regexprep(cellstr(num2str(sequencesG2(X2(Y(l),1),:))),'\s','');
                    
                    PreScanTestG2.cueP(p,:) = repmat(z(l),1,1);
                    
                elseif PreScanTestG2.hand(p,:) == 1
                    
                    z(l) = regexprep(cellstr(num2str(6-sequencesG2(X2(Y(l),1),:))),'\s','');
                    
                    PreScanTestG2.cueP(p,:) = repmat(z(l),1,1);
                    
                end
                p = p+1;
                l = l+1;
                
            end
            
            switch GroupHand
                case 1
                    dsave(fullfile(basedir,['PreScanTest/','sml_PST_G1_R_B', num2str(e), '.tgt']),orderfields(PreScanTestG1,OrderFields));
                    dsave(fullfile(basedir,['PreScanTest/','sml_PST_G2_L_B', num2str(e), '.tgt']),orderfields(PreScanTestG2,OrderFields));
                case 2
                    dsave(fullfile(basedir,['PreScanTest/','sml_PST_G2_R_B', num2str(e), '.tgt']),orderfields(PreScanTestG1,OrderFields));
                    dsave(fullfile(basedir,['PreScanTest/','sml_PST_G1_L_B', num2str(e), '.tgt']),orderfields(PreScanTestG2,OrderFields));
            end
        end
        %% pre scan met
        sequencesG1 = [SeqG1;RandomSeq([5,7,14,20],:)];
        sequencesG2 = [SeqG2;RandomSeq([5,7,14,20],:)];
        
        blocks=5;
        
        seqforRandG1 = [SeqG1;SeqG1;SeqG1;UnSeqG1;UnSeqG1;UnSeqG1];
        seqforRandG2 = [SeqG2;SeqG2;SeqG2;UnSeqG2;UnSeqG2;UnSeqG2];
        
        seqG1 = kron(seqforRandG1,ones(2,1));
        seqG2 = kron(seqforRandG2,ones(2,1));
        
        randomnumbG1(:,1) = [1:30]';
        randomnumbG1(:,2) = [1:6,1:6,1:6,7:10,7:10,7:10]';
        
        randomnumbG1(:,1) = [1:30]';
        randomnumbG2(:,2) = [1:6,1:6,1:6,7:10,7:10,7:10]';
        
                % Group 1
        PreScanMetG1.cueC = cellstr(repmat('£',length(seqG1),1));       
        PreScanMetG1.cueS = cellstr(repmat('£',length(seqG1),1));
        
        PreScanMetG1.iti(1:length(seqG1),:) = 0;        
        PreScanMetG1.sounds(1:length(seqG1),:) = 1;
        
        PreScanMetG1.FT(1:length(seqG1),:) = 6;
        
        PreScanMetG1.isScan(1:length(seqG1),:) = 1;
        PreScanMetG1.isMetronome(1:length(seqG1),:) = 1;
        
        PreScanMetG1.FoSEX(1:length(seqG1),:) = 1;
        
        PreScanMetG1.lastTrial(1:length(seqG1),:) = 0;       
        PreScanMetG1.lastTrial(length(seqG1),:) = 1;
        
        PreScanMetG1.announceTime(1:2:length(seqG1),:) = 2000;
        PreScanMetG1.announceTime(2:2:length(seqG1),:) = 0;
        PreScanMetG1.feedbackTime(1:length(seqG1),:) = 400;
        PreScanMetG1.feedbackType(1:2:length(seqG1),:) = 1;
        PreScanMetG1.feedbackType(2:2:length(seqG1),:) = 2;
        
        % Group 2
        PreScanMetG2.cueC = cellstr(repmat('£',length(seqG2),1));    
        PreScanMetG2.cueS = cellstr(repmat('£',length(seqG2),1));
        
        PreScanMetG2.iti(1:length(seqG2),:) = 0; 
        PreScanMetG2.sounds(1:length(seqG2),:) = 1;
        
        PreScanMetG2.FT(1:length(seqG2),:) = 6;
        
        PreScanMetG2.isScan(1:length(seqG2),:) = 1;
        PreScanMetG2.isMetronome(1:length(seqG2),:) = 1;
        
        PreScanMetG2.FoSEX(1:length(seqG2),:) = 1;
        
        PreScanMetG2.lastTrial(1:length(seqG2),:) = 0;
        PreScanMetG2.lastTrial(length(seqG2),:) = 1;
        
        PreScanMetG2.announceTime(1:2:length(seqG2),:) = 2000; 
        PreScanMetG2.announceTime(2:2:length(seqG2),:) = 0;
        PreScanMetG2.feedbackTime(1:length(seqG2),:) = 400; 
        PreScanMetG2.feedbackType(1:2:length(seqG2),:) = 1;
        PreScanMetG2.feedbackType(2:2:length(seqG2),:) = 2;
        
        for e = 1:blocks
            
            X1 = [sample_wor(1:size(randomnumbG1,1),size(randomnumbG1,1))];  
            X2 = kron(X1,ones(2,1));
            
            trialLength = repmat([7000;5000],length(seqG1)/2,1);  % first trial 7000, second 5000
            startRun = 5000;
            
            RandBreaks = sort(sample_wor(1:length(X2),5));
            while (sum(find(RandBreaks == 1)) || sum(rem(RandBreaks,2))<5) % breaks cannot be on first trial, between two consecutive repetitions (only odd trials)
                RandBreaks = sort(sample_wor(1:length(X2),5));
            end
            
            for i = 1:length(X2)
                
                if any(ismember(RandBreaks,i))
                    startRun = startRun+8000;   % pause of 8000 - for baseline
                end
                PreScanMetG1.startTime(i,:) = startRun;
                PreScanMetG2.startTime(i,:) = startRun;
                
                startRun = startRun+trialLength(i);
            end
               
            p = 1;
            l = 1;
            
            for i = 1:size(X2,1)
                
                
                PreScanMetG1.press1(p,:) = sequencesG1(randomnumbG1(X2(l),2),1);
                PreScanMetG1.press2(p,:) = sequencesG1(randomnumbG1(X2(l),2),2);
                PreScanMetG1.press3(p,:) = sequencesG1(randomnumbG1(X2(l),2),3);
                PreScanMetG1.press4(p,:) = sequencesG1(randomnumbG1(X2(l),2),4);
                PreScanMetG1.press5(p,:) = sequencesG1(randomnumbG1(X2(l),2),5);
                PreScanMetG1.press6(p,:) = sequencesG1(randomnumbG1(X2(l),2),6);
                PreScanMetG1.press7(p,:) = sequencesG1(randomnumbG1(X2(l),2),7);
                PreScanMetG1.press8(p,:) = sequencesG1(randomnumbG1(X2(l),2),8);
                PreScanMetG1.press9(p,:) = sequencesG1(randomnumbG1(X2(l),2),9);
                
                PreScanMetG1.seqNumb(p,:) = randomnumbG1(X2(l),2);
                PreScanMetG1.hand(p,:) = 2;
                
                if PreScanMetG1.seqNumb(p,:) <7
                    %if trained sequencesG1
                    PreScanMetG1.seqType(p,:) = 1;                  
                elseif PreScanMetG1.seqNumb(p,:) >6 
                    %untrained sequences
                    PreScanMetG1.seqType(p,:) = 3;    
                end
                
                if PreScanMetG1.hand(p,:) == 2
                    
                    z(l) = regexprep(cellstr(num2str(sequencesG1(randomnumbG1(X2(l),2),:))),'\s','');
                    
                    PreScanMetG1.cueP(p,:) = repmat(z(l),1,1);
                    
                elseif PreScanMetG1.hand(p,:) == 1
                    
                    z(l) = regexprep(cellstr(num2str(6-sequencesG1(randomnumbG1(X2(l),2),:))),'\s','');
                    
                    PreScanMetG1.cueP(p,:) = repmat(z(l),1,1);
                    
                end
                
                PreScanMetG2.press1(p,:) = sequencesG2(randomnumbG2(X2(l),2),1);
                PreScanMetG2.press2(p,:) = sequencesG2(randomnumbG2(X2(l),2),2);
                PreScanMetG2.press3(p,:) = sequencesG2(randomnumbG2(X2(l),2),3);
                PreScanMetG2.press4(p,:) = sequencesG2(randomnumbG2(X2(l),2),4);
                PreScanMetG2.press5(p,:) = sequencesG2(randomnumbG2(X2(l),2),5);
                PreScanMetG2.press6(p,:) = sequencesG2(randomnumbG2(X2(l),2),6);
                PreScanMetG2.press7(p,:) = sequencesG2(randomnumbG2(X2(l),2),7);
                PreScanMetG2.press8(p,:) = sequencesG2(randomnumbG2(X2(l),2),8);
                PreScanMetG2.press9(p,:) = sequencesG2(randomnumbG2(X2(l),2),9);
                
                PreScanMetG2.seqNumb(p,:) = randomnumbG2(X2(l),2);
                PreScanMetG2.hand(p,:) = 1;
                
                if PreScanMetG2.seqNumb(p,:) <7
                    %if trained sequencesG2
                    PreScanMetG2.seqType(p,:) = 1;
                elseif PreScanMetG2.seqNumb(p,:) >6 
                    %untrained
                    PreScanMetG2.seqType(p,:) = 3;  
                end
                
                if PreScanMetG2.hand(p,:) == 2
                    
                    z(l) = regexprep(cellstr(num2str(sequencesG2(randomnumbG2(X2(l),2),:))),'\s','');
                    
                    PreScanMetG2.cueP(p,:) = repmat(z(l),1,1);
                    
                elseif PreScanMetG2.hand(p,:) == 1
                    
                    z(l) = regexprep(cellstr(num2str(6-sequencesG2(randomnumbG2(X2(l),2),:))),'\s','');
                    
                    PreScanMetG2.cueP(p,:) = repmat(z(l),1,1);
                    
                end
                
                p = p+1;
                l = l+1;                  
            
                
                
            end
            switch GroupHand
                case 1
                    dsave(fullfile(basedir,['PreScanMetronome/','sml_PSM_G1_R_B', num2str(e), '.tgt']),orderfields(PreScanMetG1,OrderFields));
                    dsave(fullfile(basedir,['PreScanMetronome/','sml_PSM_G2_L_B', num2str(e), '.tgt']),orderfields(PreScanMetG2,OrderFields));
                case 2
                    dsave(fullfile(basedir,['PreScanMetronome/','sml_PSM_G2_R_B', num2str(e), '.tgt']),orderfields(PreScanMetG1,OrderFields));
                    dsave(fullfile(basedir,['PreScanMetronome/','sml_PSM_G1_L_B', num2str(e), '.tgt']),orderfields(PreScanMetG2,OrderFields));
            end
        end
        
        %% Scan Days
        
        % 10 functional runs
        % 10 seq - each 2x2 times (repetition) - 6 trained + 4 untrained
        % metronome 
        % 8.9333min per run - 100min pure scanning - 2 hours
        runs = 8;

        sequencesG1 = [SeqG1;UnSeqG1];
        sequencesG2 = [SeqG2;UnSeqG2];
        
        seqforRandG1 = [SeqG1;SeqG1;SeqG1;UnSeqG1;UnSeqG1;UnSeqG1];
        seqforRandG2 = [SeqG2;SeqG2;SeqG2;UnSeqG2;UnSeqG2;UnSeqG2];
        
        seqG1 = kron(seqforRandG1,ones(2,1));
        seqG2 = kron(seqforRandG2,ones(2,1));
        
        randomnumbG1(:,1) = [1:30]';
        randomnumbG1(:,2) = [1:6,1:6,1:6,7:10,7:10,7:10]';
        
        randomnumbG1(:,1) = [1:30]';
        randomnumbG2(:,2) = [1:6,1:6,1:6,7:10,7:10,7:10]';
        
        % Group 1
        ScanG1.cueC = cellstr(repmat('£',length(seqG1),1));       
        ScanG1.cueS = cellstr(repmat('£',length(seqG1),1));
        
        ScanG1.iti(1:length(seqG1),:) = 0;        
        ScanG1.sounds(1:length(seqG1),:) = 1;
        
        ScanG1.FT(1:length(seqG1),:) = 6;
        
        ScanG1.isScan(1:length(seqG1),:) = 1;
        ScanG1.isMetronome(1:length(seqG1),:) = 1;
        
        ScanG1.FoSEX(1:length(seqG1),:) = 1;
        
        ScanG1.lastTrial(1:length(seqG1),:) = 0;       
        ScanG1.lastTrial(length(seqG1),:) = 1;
        
        ScanG1.announceTime(1:2:length(seqG1),:) = 2000;
        ScanG1.announceTime(2:2:length(seqG1),:) = 0;
        ScanG1.feedbackTime(1:length(seqG1),:) = 400;
        ScanG1.feedbackType(1:2:length(seqG1),:) = 1;
        ScanG1.feedbackType(2:2:length(seqG1),:) = 2;
        
        % Group 2
        ScanG2.cueC = cellstr(repmat('£',length(seqG2),1));    
        ScanG2.cueS = cellstr(repmat('£',length(seqG2),1));
        
        ScanG2.iti(1:length(seqG2),:) = 0; 
        ScanG2.sounds(1:length(seqG2),:) = 1;
        
        ScanG2.FT(1:length(seqG2),:) = 6;
        
        ScanG2.isScan(1:length(seqG2),:) = 1;
        ScanG2.isMetronome(1:length(seqG2),:) = 1;
        
        ScanG2.FoSEX(1:length(seqG2),:) = 1;
        
        ScanG2.lastTrial(1:length(seqG2),:) = 0;
        ScanG2.lastTrial(length(seqG2),:) = 1;
        
        ScanG2.announceTime(1:2:length(seqG2),:) = 2000; 
        ScanG2.announceTime(2:2:length(seqG2),:) = 0;
        ScanG2.feedbackTime(1:length(seqG2),:) = 400; 
        ScanG2.feedbackType(1:2:length(seqG2),:) = 1;
        ScanG2.feedbackType(2:2:length(seqG2),:) = 2;
        
        for e = 1:runs
                 
            X1 = [sample_wor(1:size(randomnumbG1,1),size(randomnumbG1,1))];  
            X2 = kron(X1,ones(2,1));
            
            trialLength = repmat([7000;5000],length(seqG1)/2,1);  % first trial 7000, second 5000
            startRun = 5000;
            
            RandBreaks = sort(sample_wor(1:length(X2),5));
            while (sum(find(RandBreaks == 1)) || sum(rem(RandBreaks,2))<5) % breaks cannot be on first trial, between two consecutive repetitions (only odd trials)
                RandBreaks = sort(sample_wor(1:length(X2),5));
            end
            
            for i = 1:length(X2)
                
                if any(ismember(RandBreaks,i))
                    startRun = startRun+8000;   % pause of 8000 - for baseline
                end
                ScanG1.startTime(i,:) = startRun;
                ScanG2.startTime(i,:) = startRun;
                
                startRun = startRun+trialLength(i);
            end
               
            p = 1;
            l = 1;
            
            for i = 1:size(X2,1)
                
                
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
                    %untrained sequences
                    ScanG1.seqType(p,:) = 2;    
                end
                
                if ScanG1.hand(p,:) == 2
                    
                    z(l) = regexprep(cellstr(num2str(sequencesG1(randomnumbG1(X2(l),2),:))),'\s','');
                    
                    ScanG1.cueP(p,:) = repmat(z(l),1,1);
                    
                elseif ScanG1.hand(p,:) == 1
                    
                    z(l) = regexprep(cellstr(num2str(6-sequencesG1(randomnumbG1(X2(l),2),:))),'\s','');
                    
                    ScanG1.cueP(p,:) = repmat(z(l),1,1);
                    
                end
                
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
                ScanG2.hand(p,:) = 1;
                
                if ScanG2.seqNumb(p,:) <7
                    %if trained sequencesG2
                    ScanG2.seqType(p,:) = 1;
                elseif ScanG2.seqNumb(p,:) >6 
                    %untrained
                    ScanG2.seqType(p,:) = 2;  
                end
                
                if ScanG2.hand(p,:) == 2
                    
                    z(l) = regexprep(cellstr(num2str(sequencesG2(randomnumbG2(X2(l),2),:))),'\s','');
                    
                    ScanG2.cueP(p,:) = repmat(z(l),1,1);
                    
                elseif ScanG2.hand(p,:) == 1
                    
                    z(l) = regexprep(cellstr(num2str(6-sequencesG2(randomnumbG2(X2(l),2),:))),'\s','');
                    
                    ScanG2.cueP(p,:) = repmat(z(l),1,1);
                    
                end
                
                p = p+1;
                l = l+1;
                
            end
            
            switch GroupHand
                case 1
                    dsave(fullfile(basedir,['Scan/','sml_Scan_G1_R_B', num2str(e), '.tgt']),orderfields(ScanG1,OrderFields));
                    dsave(fullfile(basedir,['Scan/','sml_Scan_G2_L_B', num2str(e), '.tgt']),orderfields(ScanG2,OrderFields));
                case 2
                    dsave(fullfile(basedir,['Scan/','sml_Scan_G2_R_B', num2str(e), '.tgt']),orderfields(ScanG1,OrderFields));
                    dsave(fullfile(basedir,['Scan/','sml_Scan_G1_L_B', num2str(e), '.tgt']),orderfields(ScanG2,OrderFields));
            end
        end
        
        %% Scan localiser
        
        %sequences - FingMap
        seq = repmat(FingMap,4,1);
        runs = 2;   % 2 localisers
        
        randomnumb(:,1) = [1:20]';
        randomnumb(:,2) = repmat([1:5]',4,1);
        
        Scan.cueC = cellstr(repmat('£',length(seq),1));       
        Scan.cueS = cellstr(repmat('£',length(seq),1));
        
        Scan.iti(1:length(seq),:) = 0;        
        Scan.sounds(1:length(seq),:) = 1;
        
        Scan.FT(1:length(seq),:) = 6;
        
        Scan.isScan(1:length(seq),:) = 1;
        Scan.isMetronome(1:length(seq),:) = 1;
        
        Scan.FoSEX(1:length(seq),:) = 1;
        
        Scan.lastTrial(1:length(seq),:) = 0;       
        Scan.lastTrial(length(seq),:) = 1;
        
        Scan.announceTime(1:length(seq),:) = 2000;
        Scan.feedbackTime(1:length(seq),:) = 400;
        Scan.feedbackType(1:length(seq),:) = 2;
        
        for e = 1:runs
        
        X2 = [sample_wor(1:size(FingMap,1),size(FingMap,1));...
             sample_wor(1:size(FingMap,1),size(FingMap,1));...
             sample_wor(1:size(FingMap,1),size(FingMap,1));...
             sample_wor(1:size(FingMap,1),size(FingMap,1))];
       
        trialLength = 7000;  
        startRun = 5000;
        
        RandBreaks = sort(sample_wor(1:length(X2),5));
        while (find(RandBreaks == 1)) % breaks cannot be on the first trial
            RandBreaks = sort(sample_wor(1:length(X2),5));
        end
        
        for i = 1:length(X2)
            
            if any(ismember(RandBreaks,i))
                startRun = startRun+14000;   % pause of 14000 - for baseline
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
            
            Scan.seqNumb(p,:) = randomnumb(X2(l),2)+10; % 11-15
            Scan.hand(p,:) = 2;
            Scan.seqType(p,:) = 3;  % finger mapping
            
            z(l) = regexprep(cellstr(num2str(seq(randomnumb(X2(l),2),:))),'\s','');
            
            Scan.cueP(p,:) = repmat(z(l),1,1);
            
            p = p+1;
            l = l+1;
            
        end
        dsave(fullfile(basedir,['Scan/','sml_Scan_FingMap', num2str(e),'.tgt']),orderfields(Scan,OrderFields));
        
        end
 
        %% PreScan Go
        
        sequencesG1 = [SeqG1;OtherHandG1];
        sequencesG2 = [SeqG2;OtherHandG2];
        
        for e = 1:blocks
            
            
            X2 = [sample_wor(1:size(sequencesG1,1),size(sequencesG1,1));...
                sample_wor(1:size(sequencesG1,1),size(sequencesG1,1));...
                sample_wor(1:size(sequencesG1,1),size(sequencesG1,1));...
                sample_wor(1:size(sequencesG1,1),size(sequencesG1,1))];
            
            Y = [sample_wor(1:size(X2,1),size(X2,1))];
            
            
            PreScanGoG1.cueC = cellstr(repmat('£',length(X2),1));
            
            PreScanGoG1.cueS = cellstr(repmat('£',length(X2),1));
            
            PreScanGoG1.iti(1:length(X2),:) = 500;
            
            PreScanGoG1.sounds(1:length(X2),:) = 1;
            
            PreScanGoG1.announceTime(1:2:length(Trials),:) = 2000;
            PreScanGoG2.announceTime(1:2:length(Trials),:) = 2000;
            PreScanGoG1.announceTime(2:2:length(Trials),:) = 0;
            PreScanGoG2.announceTime(2:2:length(Trials),:) = 0;
            PreScanGoG1.feedbackTime(1:length(Trials),:) = 350;
            PreScanGoG2.feedbackTime(1:length(Trials),:) = 350;
            PreScanGoG1.feedbackType(1:2:length(Trials),:) = 1;
            PreScanGoG1.feedbackType(2:2:length(Trials),:) = 2;
            PreScanGoG2.feedbackType(1:2:length(Trials),:) = 1;
            PreScanGoG2.feedbackType(2:2:length(Trials),:) = 2;
            
            PreScanGoG1.FT(1:length(X2),:) = 5;
            
            PreScanGoG1.isScan(1:length(X2),:) = 1;
            
            PreScanGoG1.isMetronome(1:length(X2),:) = 2;
            
            PreScanGoG1.FoSEX =  ones(length(X2),1);
            
            PreScanGoG1.lastTrial(1:length(X2),:) = 0;
            
            PreScanGoG1.lastTrial(length(X2),:) = 1;
            
            
            PreScanGoG2.cueC = cellstr(repmat('£',length(X2),1));
            
            PreScanGoG2.cueS = cellstr(repmat('£',length(X2),1));
            
            PreScanGoG2.iti(1:length(X2),:) = 500;
            
            PreScanGoG2.sounds(1:length(X2),:) = 1;
            
            PreScanGoG2.FT(1:length(X2),:) = 5;
            
            PreScanGoG2.isScan(1:length(X2),:) = 1;
            
            PreScanGoG2.isMetronome(1:length(X2),:) = 2;
            
            PreScanGoG2.FoSEX =  ones(length(X2),1);
            
            PreScanGoG2.lastTrial(1:length(X2),:) = 0;
            
            PreScanGoG2.lastTrial(length(X2),:) = 1;
            
            
            p = 1;
            l = 1;
            
            for i = 1:size(X2,1)
                
                
                
                PreScanGoG1.press1(p,:) = sequencesG1(X2(Y(l),1),1);
                PreScanGoG1.press2(p,:) = sequencesG1(X2(Y(l),1),2);
                PreScanGoG1.press3(p,:) = sequencesG1(X2(Y(l),1),3);
                PreScanGoG1.press4(p,:) = sequencesG1(X2(Y(l),1),4);
                PreScanGoG1.press5(p,:) = sequencesG1(X2(Y(l),1),5);
                PreScanGoG1.press6(p,:) = sequencesG1(X2(Y(l),1),6);
                PreScanGoG1.press7(p,:) = sequencesG1(X2(Y(l),1),7);
                PreScanGoG1.press8(p,:) = sequencesG1(X2(Y(l),1),8);
                PreScanGoG1.press9(p,:) = sequencesG1(X2(Y(l),1),9);
                PreScanGoG1.seqNumb(p,:) = X2(l);
                
                if PreScanGoG1.seqNumb(p,:) <7
                    %if trained sequencesG1
                    PreScanGoG1.seqType(p,:) = 1;
                    PreScanGoG1.hand(p,:) = 1;
                elseif PreScanGoG1.seqNumb(p,:) >6
                    PreScanGoG1.hand(p,:) = 2;
                    PreScanGoG1.seqType(p,:) = 3;
                end
                
                if PreScanGoG1.hand(p,:) == 2
                    
                    z(l) = regexprep(cellstr(num2str(sequencesG1(X2(Y(l),1),:))),'\s','');
                    
                    PreScanGoG1.cueP(p,:) = repmat(z(l),1,1);
                    
                elseif PreScanGoG1.hand(p,:) == 1
                    
                    z(l) = regexprep(cellstr(num2str(6-sequencesG1(X2(Y(l),1),:))),'\s','');
                    
                    PreScanGoG1.cueP(p,:) = repmat(z(l),1,1);
                    
                end
                
                trialLength = 8700;
                startRun = 5000;
                
                RandBreaks = sort(sample_wor(1:length(X2),5));
                while find(RandBreaks ==1)
                    RandBreaks = sort(sample_wor(1:length(X2),5));
                end
                
                for i = 1:length(X2)
                    
                    if any(ismember(RandBreaks,i))
                        startRun = startRun+trialLength;
                    end
                    PreScanGoG1.startTime(i,:) = startRun;
                    PreScanGoG1.startTime(i,:) = startRun;
                    
                    startRun = startRun+trialLength;
                end
                
                
                
                PreScanGoG2.press1(p,:) = sequencesG2(X2(Y(l),1),1);
                PreScanGoG2.press2(p,:) = sequencesG2(X2(Y(l),1),2);
                PreScanGoG2.press3(p,:) = sequencesG2(X2(Y(l),1),3);
                PreScanGoG2.press4(p,:) = sequencesG2(X2(Y(l),1),4);
                PreScanGoG2.press5(p,:) = sequencesG2(X2(Y(l),1),5);
                PreScanGoG2.press6(p,:) = sequencesG2(X2(Y(l),1),6);
                PreScanGoG2.press7(p,:) = sequencesG2(X2(Y(l),1),7);
                PreScanGoG2.press8(p,:) = sequencesG2(X2(Y(l),1),8);
                PreScanGoG2.press9(p,:) = sequencesG2(X2(Y(l),1),9);
                PreScanGoG2.seqNumb(p,:) = X2(l);
                
                
                if PreScanGoG2.seqNumb(p,:) <7
                    %if trained sequencesG1
                    PreScanGoG2.seqType(p,:) = 1;
                    PreScanGoG2.hand(p,:) = 1;
                    
                elseif PreScanGoG2.seqNumb(p,:) >6
                    PreScanGoG2.hand(p,:) = 2;
                    PreScanGoG2.seqType(p,:) = 3;
                end
                
                if PreScanGoG2.hand(p,:) == 2
                    
                    z(l) = regexprep(cellstr(num2str(sequencesG2(X2(Y(l),1),:))),'\s','');
                    
                    PreScanGoG2.cueP(p,:) = repmat(z(l),1,1);
                    
                elseif PreScanGoG2.hand(p,:) == 1
                    
                    z(l) = regexprep(cellstr(num2str(6-sequencesG2(X2(Y(l),1),:))),'\s','');
                    
                    PreScanGoG2.cueP(p,:) = repmat(z(l),1,1);
                    
                end
                
                trialLength = 8700;
                startRun = 5000;
                
                RandBreaks = sort(sample_wor(1:length(X2),5));
                while find(RandBreaks ==1)
                    RandBreaks = sort(sample_wor(1:length(X2),5));
                end
                
                for i = 1:length(X2)
                    
                    if any(ismember(RandBreaks,i))
                        startRun = startRun+trialLength;
                    end
                    PreScanGoG2.startTime(i,:) = startRun;
                    PreScanGoG2.startTime(i,:) = startRun;
                    
                    startRun = startRun+trialLength;
                end
                
                p = p+1;
                l = l+1;
                
            end
            switch GroupHand
                case 1
                    dsave(fullfile(basedir,['PreScanGo/','sml_PSMGO_G1_R_B', num2str(e), '.tgt']),orderfields(PreScanGoG1,OrderFields));
                    dsave(fullfile(basedir,['PreScanGo/','sml_PSMGO_G2_L_B', num2str(e), '.tgt']),orderfields(PreScanGoG2,OrderFields));
                case 2
                    dsave(fullfile(basedir,['PreScanGo/','sml_PSMGO_G2_R_B', num2str(e), '.tgt']),orderfields(PreScanGoG1,OrderFields));
                    dsave(fullfile(basedir,['PreScanGo/','sml_PSMGO_G1_L_B', num2str(e), '.tgt']),orderfields(PreScanGoG2,OrderFields));
            end
        end
        %% Scan Go last scan session
        
        % 8 runs
        % 14 seq * 3 + 5fing + 7pauses= 62+5blank Trials + 1min pause
        % metronome
        % 8.9333min per run - 100min pure scanning - 2 hours
        runs = 8;
        
        sequencesG1 = [SeqG1;UnSeqG1;OtherHandG1;FingMap];
        sequencesG2 = [SeqG2;UnSeqG2;OtherHandG2;FingMap];
        
        seqforRandG1 = [SeqG1;SeqG1;SeqG1;UnSeqG1;UnSeqG1;UnSeqG1;OtherHandG1;OtherHandG1;OtherHandG1;FingMap];
        seqforRandG2 = [SeqG2;SeqG2;SeqG2;UnSeqG2;UnSeqG2;UnSeqG2;OtherHandG2;OtherHandG2;OtherHandG2;FingMap];
        
        randomnumbG1(:,1) = [1:47]';
        randomnumbG1(:,2) = [1:6,1:6,1:6,7:10,7:10,7:10,11:14,11:14,11:14,15:19]';
        
        randomnumbG1(:,1) = [1:47]';
        randomnumbG2(:,2) = [1:6,1:6,1:6,7:10,7:10,7:10,11:14,11:14,11:14,15:19]';
        
        
        ScanGoG1.cueC = cellstr(repmat('£',length(seqforRandG1),1));
        
        ScanGoG1.cueS = cellstr(repmat('£',length(seqforRandG1),1));
        
        ScanGoG1.iti(1:length(seqforRandG1),:) = 0;
        
        ScanGoG1.sounds(1:length(seqforRandG1),:) = 1;
        
        ScanGoG1.announceTime(1:2:length(Trials),:) = 2000;
        ScanGoG2.announceTime(1:2:length(Trials),:) = 2000;
        ScanGoG1.announceTime(2:2:length(Trials),:) = 0;
        ScanGoG2.announceTime(2:2:length(Trials),:) = 0;
        ScanGoG1.feedbackTime(1:length(Trials),:) = 350;
        ScanGoG2.feedbackTime(1:length(Trials),:) = 350;
        ScanGoG1.feedbackType(1:2:length(Trials),:) = 1;
        ScanGoG1.feedbackType(2:2:length(Trials),:) = 2;
        ScanGoG2.feedbackType(1:2:length(Trials),:) = 1;
        ScanGoG2.feedbackType(2:2:length(Trials),:) = 2;
        
        ScanGoG1.FT(1:length(seqforRandG1),:) = 6;
        
        ScanGoG1.isScan(1:length(seqforRandG1),:) = 1;
        
        ScanGoG1.isMetronome(1:length(seqforRandG1),:) = 2;
        
        ScanGoG1.FoSEX(1:length(seqforRandG1),:) = 1;
        
        ScanGoG1.lastTrial(1:length(seqforRandG1),:) = 0;
        
        ScanGoG1.lastTrial(length(seqforRandG1),:) = 1;
        
        %
        ScanGoG2.cueC = cellstr(repmat('£',length(seqforRandG1),1));
        
        ScanGoG2.cueS = cellstr(repmat('£',length(seqforRandG1),1));
        
        ScanGoG2.iti(1:length(seqforRandG1),:) = 0;
        
        ScanGoG2.sounds(1:length(seqforRandG1),:) = 1;
        
        ScanGoG2.FT(1:length(seqforRandG1),:) = 6;
        
        ScanGoG2.isScan(1:length(seqforRandG1),:) = 1;
        
        ScanGoG2.isMetronome(1:length(seqforRandG1),:) = 2;
        
        ScanGoG2.FoSEX(1:length(seqforRandG1),:) = 1;
        
        ScanGoG2.lastTrial(1:length(seqforRandG1),:) = 0;
        
        ScanGoG2.lastTrial(length(seqforRandG1),:) = 1;
        
        for e = 1:runs
            
            
            X2 = [sample_wor(1:size(randomnumbG1,1),size(randomnumbG1,1))];
            
            
            trialLength = 8700;
            startRun = 5000;
            
            RandBreaks = sort(sample_wor(1:length(X2),5));
            while find(RandBreaks ==1)
                RandBreaks = sort(sample_wor(1:length(X2),5));
            end
            
            for i = 1:length(X2)
                
                if any(ismember(RandBreaks,i))
                    startRun = startRun+12000;  % 12000 for baseline estimation - pause
                end
                ScanGoG1.startTime(i,:) = startRun;
                ScanGoG2.startTime(i,:) = startRun;
                
                startRun = startRun+trialLength;
            end
            
            
            p = 1;
            l = 1;
            
            for i = 1:size(X2,1)
                
                
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
                
                if ScanGoG1.seqNumb(p,:) <7
                    %if trained sequencesG1
                    ScanGoG1.seqType(p,:) = 1;
                    ScanGoG1.hand(p,:) = 2;
                elseif ScanGoG1.seqNumb(p,:) >6 &ScanGoG1.seqNumb(p,:) <11
                    ScanGoG1.hand(p,:) = 2;
                    ScanGoG1.seqType(p,:) = 2;
                elseif ScanGoG1.seqNumb(p,:) >10 & ScanGoG1.seqNumb(p,:)<15
                    ScanGoG1.hand(p,:) = 1;
                    ScanGoG1.seqType(p,:) = 3;
                    
                elseif ScanGoG1.seqNumb(p,:) >14
                    ScanGoG1.hand(p,:) = 2;
                    ScanGoG1.seqType(p,:) = 5;
                end
                
                if ScanGoG1.hand(p,:) == 2
                    
                    z(l) = regexprep(cellstr(num2str(sequencesG1(randomnumbG1(X2(l),2),:))),'\s','');
                    
                    ScanGoG1.cueP(p,:) = repmat(z(l),1,1);
                    
                elseif ScanGoG1.hand(p,:) == 1
                    
                    z(l) = regexprep(cellstr(num2str(6-sequencesG1(randomnumbG1(X2(l),2),:))),'\s','');
                    
                    ScanGoG1.cueP(p,:) = repmat(z(l),1,1);
                    
                end
                
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
                
                if ScanGoG2.seqNumb(p,:) <7
                    %if trained sequencesG2
                    ScanGoG2.seqType(p,:) = 1;
                    ScanGoG2.hand(p,:) = 1;
                elseif ScanGoG2.seqNumb(p,:) >6 &ScanGoG2.seqNumb(p,:) <11
                    ScanGoG2.hand(p,:) = 1;
                    ScanGoG2.seqType(p,:) = 2;
                elseif ScanGoG2.seqNumb(p,:) >10 & ScanGoG2.seqNumb(p,:)<15
                    ScanGoG2.hand(p,:) = 2;
                    ScanGoG2.seqType(p,:) = 3;
                    
                elseif ScanGoG2.seqNumb(p,:) >14
                    ScanGoG2.hand(p,:) = 1;
                    ScanGoG2.seqType(p,:) = 5;
                end
                
                if ScanGoG2.hand(p,:) == 2
                    
                    z(l) = regexprep(cellstr(num2str(sequencesG2(randomnumbG2(X2(l),2),:))),'\s','');
                    
                    ScanGoG2.cueP(p,:) = repmat(z(l),1,1);
                    
                elseif ScanGoG2.hand(p,:) == 1
                    
                    z(l) = regexprep(cellstr(num2str(6-sequencesG2(randomnumbG2(X2(l),2),:))),'\s','');
                    
                    ScanGoG2.cueP(p,:) = repmat(z(l),1,1);
                    
                end
                
                p = p+1;
                l = l+1;
                
            end
            
            switch GroupHand
                case 1
                    dsave(fullfile(basedir,['ScanGo/','sml_ScanGo_G1_R_B', num2str(e), '.tgt']),orderfields(ScanGoG1,OrderFields));
                    dsave(fullfile(basedir,['ScanGo/','sml_ScanGo_G2_L_B', num2str(e), '.tgt']),orderfields(ScanGoG2,OrderFields));
                case 2
                    dsave(fullfile(basedir,['ScanGo/','sml_ScanGo_G2_R_B', num2str(e), '.tgt']),orderfields(ScanGoG1,OrderFields));
                    dsave(fullfile(basedir,['ScanGo/','sml_ScanGo_G1_L_B', num2str(e), '.tgt']),orderfields(ScanGoG2,OrderFields));
                    
            end
        end
end

end





