function sml1_bids(what,varargin)
%function sml1_bids(what,varargin)
% makes current folder into bids files / filenames
%bidsDir ='/Volumes/MotorControl/data/longtermseq_bids';
bidsDir ='/Volumes/Eva_Data/longtermseq_bids';
origDir ='/Volumes/MotorControl/data/SuperMotorLearning';
locDir  = '/Users/eberlot/Documents/Data/SuperMotorLearning'; % local dir
eventDir = '/Volumes/Eva_Data/events'; % temporary for the events.tsv

sn_new = 1:26;
sn_old = [5:9,11:31];
runID = [1:3 5:7 9:10;
         11:13 15:17 19:20;
         21:23 25:27 29:30;
         31:33 35:37 39:40];  
taskName = 'motorseq';
switch what
    case 'CREATE:folder_struct'
        % here just create the needed folder structure
        for s=sn_new
            subDir = fullfile(bidsDir,sprintf('sub-%02d',s));
            mkdir(subDir);
            for ss=1:4 % per session
                sesDir = fullfile(subDir,sprintf('ses-%02d',ss));
                mkdir(sesDir);
                funcDir = fullfile(sesDir,'func');
                fmapDir = fullfile(sesDir,'fmap');
                mkdir(funcDir); mkdir(fmapDir);
                if ss==1 % add anatomy
                    anatDir = fullfile(sesDir,'anat');
                    mkdir(anatDir);
                end
                fprintf('Done subject-%d, sess-%d.\n',s,ss);
            end
        end
    case 'COPY:anatomical'
        % copy the anatomical raw file
        % also rename it
        for s=sn_new
            INfile = fullfile(origDir,'anatomicals',sprintf('s%02d',sn_old(s)),sprintf('s%02d_anatomical_raw.nii',sn_old(s)));
            OUTfile = fullfile(bidsDir,sprintf('sub-%02d',s),'ses-01','anat',sprintf('sub-%02d_ses-01_T1w.nii',s));
            copyfile(INfile,OUTfile,'f');
            gzip(OUTfile);
            delete(OUTfile);
            fprintf('Done subject-%d.\n',s);
        end 
    case 'COPY:fmap'
        % copy fieldmap files
        sn = 1:26;
        sess = 1:4;
        vararginoptions(varargin,{'sn','sess'});
        for s=sn
            for ss=sess
                % phase file
                INfile = fullfile(origDir,'fieldmaps',sprintf('s%02d',sn_old(s)),sprintf('sess%d',ss),sprintf('s%02d_phase.nii',sn_old(s)));
                OUTfile = fullfile(bidsDir,sprintf('sub-%02d',s),sprintf('ses-%02d',ss),'fmap',sprintf('sub-%02d_ses-%02d_phasediff.nii',s,ss));
                copyfile(INfile,OUTfile,'f');
                gzip(OUTfile);
                delete(OUTfile);
                % magnitude file
                INfile = fullfile(origDir,'fieldmaps',sprintf('s%02d',sn_old(s)),sprintf('sess%d',ss),sprintf('s%02d_magnitude.nii',sn_old(s)));
                OUTfile = fullfile(bidsDir,sprintf('sub-%02d',s),sprintf('ses-%02d',ss),'fmap',sprintf('sub-%02d_ses-%02d_magnitude1.nii',s,ss));
                copyfile(INfile,OUTfile,'f');
                gzip(OUTfile);
                delete(OUTfile);
                fprintf('Done subject-%d, sess-%d.\n',s,ss);
            end
        end
    case 'COPY:func'
        % copy functional data
        sn = 1:26;
        sess = 1:4;
        vararginoptions(varargin,{'sn','sess'});
        fprintf('Copying functional data:\n');
        for s=sn
            fprintf('Subject-%d:\n',s);
            for ss=sess
                fprintf('-session %d:',ss);
                tElapsed = tic;
                for r=1:8
                    INfile = fullfile(origDir,'imaging_data_raw',sprintf('s%02d',sn_old(s)),sprintf('sess%d',ss),sprintf('s%02d_run_%02d.nii',sn_old(s),runID(ss,r)));
                    OUTfile = fullfile(bidsDir,sprintf('sub-%02d',s),sprintf('ses-%02d',ss),'func',sprintf('sub-%02d_ses-%02d_task-%s_run-%02d_bold.nii',s,ss,taskName,r));
                    copyfile(INfile,OUTfile,'f');
                    fprintf('%d.',r);
                end
                fprintf('...done\n');
                toc(tElapsed);
            end
        end
        
    case 'GUNZIP:fmap'
        fprintf('Gunzipping fieldmap data:\n');
         for s=sn_new
            for ss=1:4
                % phase file
                INfile = fullfile(bidsDir,sprintf('sub-%02d',s),sprintf('ses-%02d',ss),'fmap',sprintf('sub-%02d_ses-%02d_phasediff.nii',s,ss));
                gzip(INfile);
                delete(INfile);
                % magnitude file
                INfile = fullfile(bidsDir,sprintf('sub-%02d',s),sprintf('ses-%02d',ss),'fmap',sprintf('sub-%02d_ses-%02d_magnitude1.nii',s,ss));
                gzip(INfile);
                delete(INfile);
                fprintf('Done subject-%d, sess-%d.\n',s,ss);
            end
        end
    case 'GUNZIP:func'
        sn = 1:26;
        sess = 1:4;
        vararginoptions(varargin,{'sn','sess'});
        % copy functional data
        fprintf('Gunzipping functional data:\n');
        for s=sn
            fprintf('Subject-%d:\n',s);
            for ss=sess
                fprintf('-session %d:',ss);
                tElapsed = tic;
                for r=1:8
                    INfile = fullfile(bidsDir,sprintf('sub-%02d',s),sprintf('ses-%02d',ss),'func',sprintf('sub-%02d_ses-%02d_task-%s_run-%02d_bold.nii',s,ss,taskName,r));
                    gzip(INfile);
                    delete(INfile);
                    fprintf('%d.',r);
                end
                fprintf('...done\n');
                toc(tElapsed);
            end
        end
        
    case 'CREATE:events'
        % create events tsv
        sn = 1:26;
        sess = 1:4;
        vararginoptions(varargin,{'sn','sess'});
        format shortg;
        for s=sn
            D = dload(fullfile(locDir,'behavioral_data','data',sprintf('sml1_s%02d.dat',sn_old(s))));
            for ss=sess
                L = getrow(D,D.ScanSess==ss);    % only blocks of that scan session
                if ss == 4
                    Rr = getrow(L,L.blockType==9); %blockType==9 - func imaging run without metronome
                else
                    Rr = getrow(L,L.blockType==3); %blockType==3 - funct imaging run with metornome
                end   
                uniqrun=unique(Rr.BN);
                for r=1:8 % each run
                    R = getrow(Rr,Rr.BN==uniqrun(r)); 
                    % necessary fields: onset, duration
                    % others: is_error, points, movement_time
                    onset = R.startTimeReal./1000 - 4; % 4 dummies
                    duration = ones(size(onset))*2.5;
                    is_error = R.numError~=0;
                    points = ~is_error.*3;
                    movement_time = R.MT;
                    movement_time(movement_time==0 & is_error==0)=randi([2400,2600],1,1); % insert for other MT
                    movement_time = movement_time./1000;
                    movement_time = round(movement_time,2);
                    sequence_num = R.seqNumb;
                    sequence_type = R.seqType;
                    repetition = R.FoSEx;
                    t = table(onset,duration,is_error,points,movement_time,sequence_num,sequence_type,repetition);
                    writetable(t,fullfile(eventDir,sprintf('sub-%02d',s),sprintf('sub-%02d_ses-%02d_task-%s_run-%02d_events.tsv',s,ss,taskName,r)),'FileType','text','Delimiter','\t');
                end
                fprintf('Done ses-%d sub-%02d.\n',ss,s);
            end
        end
    case 'MOVE:events'
        % move the events tsv
        sn = 1:26;
        sess = 1:4;
        vararginoptions(varargin,{'sn','sess'});
        for s=sn
            for ss=sess
                for r=1:8
                    movefile(fullfile(eventDir,sprintf('sub-%02d',s),sprintf('sub-%02d_ses-%02d_task-%s_run-%02d_events.tsv',s,ss,taskName,r)),...
                        fullfile(bidsDir,sprintf('sub-%02d',s),sprintf('ses-%02d',ss),'func',sprintf('sub-%02d_ses-%02d_task-%s_run-%02d_events.tsv',s,ss,taskName,r)));
                end
                fprintf('Done ses-%d sub-%02d.\n',ss,s);
            end
        end
        
    case 'COPY:json_fmap'
        % here just copy and rename the json files for
        sn = 1:26;
        sess = 1:4;
        for s=sn
            for ss=sess
                fileName = sprintf('sub-%02d_ses-%02d_magnitude1.json',s,ss);
                copyfile(fullfile(bidsDir,'magnitude1.json'),fullfile(bidsDir,sprintf('sub-%02d',s),sprintf('ses-%02d',ss),'fmap',fileName));
                fileName = sprintf('sub-%02d_ses-%02d_phasediff.json',s,ss);
                copyfile(fullfile(bidsDir,'phase.json'),fullfile(bidsDir,sprintf('sub-%02d',s),sprintf('ses-%02d',ss),'fmap',fileName));
                fprintf('Done ses-%d sub-%02d.\n',ss,s);
            end
        end
        
    case 'ANAT:rename_defaced'
        sn = 2:26;
        vararginoptions(varargin,{'sn'});
        for s=sn
            fileORIG = fullfile(bidsDir,sprintf('sub-%02d',s),'ses-01','anat',sprintf('sub-%02d_ses-01_T1w.nii.gz',s)); 
            fileDEFAC = fullfile(bidsDir,sprintf('sub-%02d',s),'ses-01','anat',sprintf('sub-%02d_ses-01_T1w_defaced.nii.gz',s));
            delete(fileORIG);
            movefile(fileDEFAC,fileORIG);
        end
        
    case 'FIND:nii_file'
        % find nifti file
        sn = 1:26;
        sess = 1:4;
        for s=sn
            for ss=sess
                cd(fullfile(bidsDir,sprintf('sub-%02d',s),sprintf('ses-%02d',ss),'func'));
                d=dir('*.nii');
                if size(d,1)~=0
                    keyboard;
                end
                cd(fullfile(bidsDir,sprintf('sub-%02d',s),sprintf('ses-%02d',ss),'fmap'));
                d=dir('*.nii');
                if size(d,1)~=0
                    keyboard;
                end
                if ss==1
                    cd(fullfile(bidsDir,sprintf('sub-%02d',s),sprintf('ses-%02d',ss),'anat'));
                    d=dir('*.nii');
                    if size(d,1)~=0
                        keyboard;
                    end
                end
            end
        end
        
    case 'RENAME:func'
        % remove the 'new' part of the name
        sn = 25;
        sess = 1:4;
        vararginoptions(varargin,{'sn','sess'});
        for s=sn
            for ss=sess
                for r=1:8
                    movefile(fullfile(bidsDir,sprintf('sub-%02d',s),sprintf('ses-%02d',ss),'func',sprintf('sub-%02d_ses-%02d_task-motorseq_run-%02d_bold_new.nii.gz',s,ss,r)),...
                        fullfile(bidsDir,sprintf('sub-%02d',s),sprintf('ses-%02d',ss),'func',sprintf('sub-%02d_ses-%02d_task-motorseq_run-%02d_bold.nii.gz',s,ss,r)));
                    % delete the json file
                    delete(fullfile(bidsDir,sprintf('sub-%02d',s),sprintf('ses-%02d',ss),'func',sprintf('sub-%02d_ses-%02d_task-motorseq_run-%02d_bold_new.json',s,ss,r)));
                end
                fprintf('Done sub-%02d, ses-%02d\n',s,ss);
            end
        end
        
    case 'run_job'
        sml1_bids('COPY:fmap');
        sml1_bids('COPY:func');
        
        
    case 'test_RH' % extra
        % here create a file for right hand tests
        behDir =[locDir '/behavioral_data'];

        B=load(fullfile(behDir,'analyze','testsRH.mat'));
        B = getrow(B,B.SN~=4 & ismember(B.seqType,[1,4]));
        B.seqType(B.seqType>1) = 2;
        B.day(B.day==19)=20;
        B.day(B.day==11)=13;

        % test plotting
        [MT,idx,~] = pivottable([B.day B.SN],[B.seqType],B.MT,'nanmedian','subset',B.isError==0);
        K.MT = MT(:);
        K.seqType = [ones(size(MT,1),1);ones(size(MT,1),1)*2];
        K.day = [idx(:,1);idx(:,1)];
        K.sn  = [idx(:,2);idx(:,2)];
        figure
        plt.box(K.day,K.MT,'split',K.seqType); title('behavior');
        
        % continue
        for s=1:length(sn_old)
            B.SN(B.SN==sn_old(s))=sn_new(s);
        end
        % so seqNumb is continuous
        B.seqNumb(B.seqNumb>29) = B.seqNumb(B.seqNumb>29)-23;
       
        % structure with relevant fields
        SN = B.SN;
        %day = B.day;
        test_day = B.day;
        test_day(test_day==3)=1;
        test_day(test_day==8)=2;
        test_day(test_day==13)=3;
        test_day(test_day==20)=4;
        seq_type = B.seqType;
        seq_num = B.seqNumb;
        digit_cue = B.cueP;
        points = B.points;
        is_error = B.isError;
        movement_time = B.MT;
        repetition = B.FoSEx;
        % get into tsv
        t = table(SN,test_day,seq_type,seq_num,points,is_error,movement_time,repetition,digit_cue);
        writetable(t,'behavioralData.tsv','FileType','text','Delimiter','\t');
        keyboard;
    otherwise 
        fprintf('Wrong case!\n');
end
