function varargout = sml_connect(what,varargin)
% ------------------------- Directories -----------------------------------
baseDir         ='/Users/eberlot/Documents/Data/SuperMotorLearning';
%baseDir         ='/Volumes/MotorControl/data/SuperMotorLearning';
betaDir         =[baseDir '/betas'];
connectDir      =[baseDir '/connectivity_new'];
codeDir         ='/Users/eberlot/Documents/MATLAB/projects/SuperMotorLearning';
wbDir           =[baseDir '/surfaceWB'];
% other info
hemI            = {'L','R'};     
hemName         = {'CortexLeft','CortexRight'};   
% ------------------------- Subject things --------------------------------
% The variables in this section must be updated for every new subject.
subj_name  = {'s01','s02','s03','s04','s05','s06','s07','s08','s09','s10','s11','s12','s13','s14','s15','s16','s17','s18',...
                's19','s20','s21','s22','s23','s24','s25','s26','s27','s28','s29','s30','s31'};  
sn          = [5:9,11:31]; % final selection of subjects
sessN       = 1:4; % all sessions unless specified otherwise 
style.file(fullfile(codeDir,'sml_style.m'));
style.use('default');

switch what
    case 'RDM:intersubject'
        % calculate consistency of RDMs (per tessel) across subjects
        % upper / lower noise ceiling
        sessN = 1:4;
        hemi=1:2;
        betaChoice='multi';
        nTessel = 162;
        seqType = 'all'; % all / trained / untrained
        tesselSelect = 0; % whether or not to select tessels; otherwise all of them
        vararginoptions(varargin,{'sessN','sn','hemi','betaChoice','distType','seqType','nTessel','tesselSelect'});
        
        for h=hemi
            if tesselSelect
                %tessels2{h} = sml_connect('TESSEL:select','hemi',hemi(h),'betaChoice',betaChoice,'sessN',sessN,'nTessel',nTessel);
                tessels{h} = sml_connect('TESSEL:select','hemi',hemi(h),'betaChoice',betaChoice,'sessN',sessN,'nTessel',nTessel);
            else
                BB=load(fullfile(betaDir,'s05',sprintf('betas_tesselsWB_%d_s05_sess1',nTessel))); % used for 642
                tessels{h} = (BB.regType(BB.regSide==h))';
            end
        end
        switch seqType
            case 'all'
                metric = 'RDM';
            case 'trained'
                metric = 'RDM_train';
            case 'untrained'
                metric = 'RDM_untrain';
        end
        RR = [];
        for ss=sessN
            % load in stats file with RDMs
            TT=load(fullfile(betaDir,'group',sprintf('stats_tesselsWB_%d_%sPW_sess%d',nTessel,betaChoice,ss)));
            % only select tessels with data present
            for h=hemi
                %T = getrow(TT,ismember(TT.region,tessels{h}));
                T = getrow(TT,ismember(TT.regType,tessels{h}));
                for i = tessels{h}
                    t = getrow(T,T.regType==i&T.regSide==h);
                    %t = getrow(T,T.region==i);
                    [noise_high,noise_low] = sml_connect('CALC:noiseCeiling',t.(metric),t.SN);
                    R.noise_high = mean(noise_high);
                    R.noise_low = mean(noise_low);
                    R.tessel = i;
                    R.regSide = h;
                    R.regType = unique(t.regType);
                    R.region = unique(t.region);
                    R.sessN = ss;
                    RR = addstruct(RR,R);
                    fprintf('Done sess-%d: tessel-%d\n',ss,i);
                end
            end
            fprintf('Done sess-%d.\n\n',ss);
        end
        save(fullfile(connectDir,'RDMs',sprintf('rdm_tessel-%d_intersubj-consistency_seq-%s',nTessel,seqType)),'-struct','RR');
    case 'RDM:intersess'
        % calculate consistency of RDMs (per tessel) across sessions
        sessN = 1:4;
        hemi=1:2;
        betaChoice='multi';
        seqType = 'all'; % all / trained / untrained
        nTessel = 162;
        tesselSelect = 0;
        vararginoptions(varargin,{'sessN','sn','hemi','betaChoice','distType','seqType','nTessel','tesselSelect'});
        
        switch seqType
            case 'all'
                metric = 'RDM';
            case 'trained'
                metric = 'RDM_train';
            case 'untrained'
                metric = 'RDM_untrain';
        end
        RR = [];
        
        for ss=1:numel(sessN)
            % load in stats file with RDMs
            T{ss}=load(fullfile(betaDir,'group',sprintf('stats_tesselsWB_%d_%sPW_sess%d',nTessel,betaChoice,sessN(ss))));
        end
        for h=hemi
            fprintf('%s:\n',hemName{h})
            if tesselSelect
                tessels = sml_connect('TESSEL:select','hemi',hemi(h),'betaChoice',betaChoice,'sessN',sessN,'nTessel',nTessel);
            else
                %tessels = sml_connect('TESSEL:select_old','hemi',hemi(h),'betaChoice',betaChoice,'sessN',sessN,'nTessel',nTessel);
                BB=load(fullfile(betaDir,'s05',sprintf('betas_tesselsWB_%d_s05_sess1',nTessel))); % used for 642
                tessels = (BB.regType(BB.regSide==h))';
            end
            ind = indicatorMatrix('allpairs',(sessN));
            for ses=1:size(ind,1)
                sess1 = sessN(ind(ses,:)==1);
                sess2 = sessN(ind(ses,:)==-1);
                for s=sn
                    for i = tessels
                        t1 = getrow(T{sess1},T{sess1}.regType==i & T{sess1}.SN==s & T{sess1}.regSide==h);
                        t2 = getrow(T{sess2},T{sess2}.regType==i & T{sess2}.SN==s & T{sess2}.regSide==h);
                        rdm = [t1.(metric);t2.(metric)];
                        R.corr_sess = rsa_vectorizeRDM(corr(rdm'));
                        R.sn = s;
                        R.sess1 = sess1;
                        R.sess2 = sess2;
                        R.sessTr = ses; % transition
                        R.tessel = i;
                        R.regSide = h;
                        R.regType = unique(t1.regType);
                        R.region = unique(t1.region);
                        RR = addstruct(RR,R);
                    end
                    fprintf('Done %s - sess %d-%d\n',subj_name{s},sess1,sess2);
                end
            end
            clear tessels;
        end
        save(fullfile(connectDir,'RDMs',sprintf('rdm_tessel-%d_intersess-consistency_seq-%s',nTessel,seqType)),'-struct','RR');
    case 'RDM:tessels_surface'
        % project the RDM correlations to surface
        var = 'intersubj'; % intersubj or intersess
        seqType = {'all','trained','untrained'};
        nTessel = 362;
        hemi = 1:2;
        vararginoptions(varargin,{'var','nTessel','hemi'});
        
        regStart = 0;
        for h=hemi
            col = 1;
            G = gifti(fullfile(wbDir,'connect_FS_LR_164',sprintf('Icosahedron-%d.164k.%s.label.gii',nTessel,hemI{h})));
            nReg = numel(unique(G.cdata))-1; % exclude 0 - medial wall
            for st = seqType
                T = load(fullfile(connectDir,'RDMs',sprintf('rdm_tessel-%d_%s-consistency_seq-%s',nTessel,var,st{:})));
                allReg = (unique(T.regType(T.regSide==h)))';
                if strcmp(var,'intersubj')
                    sess = unique(T.sessN); sessM = 'sessN'; metric = {'noise_high','noise_low'};
                elseif strcmp(var,'intersess');
                    sess = unique(T.sessTr); sessM = 'sessTr'; metric = {'corr_sess'};
                end
                for str=sess'
                    if strcmp(var,'intersess')
                        ses1 = unique(T.sess1(T.sessTr==str));
                        ses2 = unique(T.sess2(T.sessTr==str));
                    end
                    for m=1:size(metric,2)
                        data(:,col) = zeros(size(G.cdata,1),1);
                        for r=allReg%nReg %regStart:(regStart+nReg)
                            t = getrow(T,T.regType==r & T.regSide==h & T.(sessM)==str);
                            idx = find(G.cdata(:,1)==r);
                            data(idx,col) = mean(t.(metric{m}));
                            clear idx;
                        end
                        if strcmp(var,'intersubj')
                            column_name{col} = sprintf('%s-corr_%s_rdm-%s_%s-%d',var,metric{m},st{:},sessM,str);
                        else
                            column_name{col} = sprintf('%s-corr_%s_rdm-%s_%s-%d-%d',var,metric{m},st{:},sessM,ses1,ses2);
                        end
                        col = col+1;
                    end
                end
            end
            outfile         = fullfile(wbDir,'connect_FS_LR_164',sprintf('%s.tesselsWB_%d_RDM-%s-consistency.func.gii',hemI{h},nTessel,var));
            G               = surf_makeFuncGifti(data,'anatomicalStruct',hemName{h},'columnNames',column_name);
            save(G,outfile);
            regStart = regStart+nReg;
        end
    case 'RDM:histogram'
        % here plot the consistency (across subjects / sessions) as a
        % histogram
        var = 'intersubj'; % intersubj or intersess
        seqType = {'all','trained','untrained'};
        nTessel = 362;
        hemi = 1:2;
        vararginoptions(varargin,{'var','nTessel','hemi'});
        
        for st=seqType
            figure
            T = load(fullfile(connectDir,'RDMs',sprintf('rdm_tessel-%d_%s-consistency_seq-%s',nTessel,var,st{:})));
            if strcmp(var,'intersubj')
                sess = unique(T.sessN); sessM = 'sessN'; metric = {'noise_high','noise_low'};
            elseif strcmp(var,'intersess');
                sess = unique(T.sessTr); sessM = 'sessTr'; metric = {'corr_sess'};
            end
            for str=sess'
                if strcmp(var,'intersess')
                    ses1 = unique(T.sess1(T.sessTr==str));
                    ses2 = unique(T.sess2(T.sessTr==str));
                end
                for m=1:size(metric,2)
                    t = getrow(T,ismember(T.hemi,hemi) & T.(sessM)==str);
                    subplot(1,numel(sess),find(str==sess));
                    histogram(t.(metric{m})); hold on;
                    drawline(mean(t.(metric{m})),'dir','vert');
                    if strcmp(var,'intersess')
                        title(sprintf('RDM %s - %d-%d',st{1},ses1,ses2));
                    else
                        title(sprintf('RDM %s - %d',st{1},str));
                    end
                end
            end
        end
        
    case 'ALPHA:calculate'
        % calculate alpha - first level - connectivity matrix
        % per session
        sessN = 1:4;
        hemi=1:2;
        betaChoice='multi';
        distType = 'cosine';
        nTessel = 362;
        seqType = 'all'; % all / trained / untrained
        tesselSelect = 0;
        vararginoptions(varargin,{'sessN','sn','hemi','betaChoice','distType','seqType','nTessel','tesselSelect'});

        if length(hemi)==2
            T=load(fullfile(betaDir,'group',sprintf('stats_tesselsWB_%d_%sPW_sess1',nTessel,betaChoice)));
            for i=1:2
                if tesselSelect
                    tess1 = sml_connect('TESSEL:select','hemi',i,'betaChoice',betaChoice,'sessN',sessN,'nTessel',nTessel);
                    tess{i} = unique(T.region(ismember(T.regType,tess1)&T.regSide==i));
                else
                    %tessels = sml_connect('TESSEL:select_old','hemi',hemi(h),'betaChoice',betaChoice,'sessN',sessN,'nTessel',nTessel);
                    BB=load(fullfile(betaDir,'s05',sprintf('betas_tesselsWB_%d_s05_sess1',nTessel))); % used for 642
                    tess{i} = (BB.regType(BB.regSide==i))';
                end
                
            end
            tessels = [tess{1},tess{2}];
            hemiType = 'CortexBoth';
        else
            tessels = sml_connect('TESSEL:select','hemi',hemi,'betaChoice',betaChoice,'sessN',sessN,'nTessel',nTessel);
            hemiType = hemName{hemi};
        end

        switch seqType
            case 'all'
                metric = 'RDM';
            case 'trained'
                metric = 'RDM_train';
            case 'untrained'
                metric = 'RDM_untrain';
        end             
        AA = [];
        for ss=sessN
            % load in stats file with RDMs
            T=load(fullfile(betaDir,'group',sprintf('stats_tesselsWB_%d_%sPW_sess%d',nTessel,betaChoice,ss)));
            % only select tessels with data present
            T = getrow(T,ismember(T.region,tessels) & ismember(T.regSide,hemi));
            for s=sn
                t = getrow(T,T.SN==s);
                A = sml_connect('CALC:distance',t.(metric),distType,1);
                A.similarity = sml_connect('CALC:gaussian_similarity',A.dist);
                A.tessel1 = tessels(A.l1);
                A.tessel2 = tessels(A.l2);
                A.sn = s;
                A.sessN = ss;
                AA = addstruct(AA,A);
                fprintf('Done sess-%d: %s\n',ss,subj_name{s});
            end
            fprintf('Done sess-%d, RDM: %s.\n\n',ss,seqType);
        end
        save(fullfile(connectDir,'alpha',sprintf('alpha_%s_nTessel-%d_seq-%s_%s',hemiType,nTessel,seqType,distType)),'-struct','AA');        
    case 'ALPHA:intersubject'
        % intersubject correlation (low / high noise ceiling)
        distType = 'cosine'; % cosine, correlation, cosineW
        seqType = 'all'; % all / trained / untrained
        hemiType = 'CortexLeft'; % CortexLeft, CortexRight,CortexBoth
        nTessel = 362;
        sessN = 1:4;
        vararginoptions(varargin,{'sessN','distType','seqType','nTessel','hemiType'});
        A = load(fullfile(connectDir,'alpha',sprintf('alpha_%s_nTessel-%d_seq-%s_%s',hemiType,nTessel,seqType,distType)));
        figure
        for ss=sessN
            a = getrow(A,A.sessN==ss);
            [noise_high,noise_low]=sml_connect('CALC:noiseCeiling',a.dist,a.sn);
            subplot(1,numel(sessN),find(ss==sessN))
            histogram(noise_high); hold on; histogram(noise_low); drawline(mean(noise_low),'dir','vert');
            xlabel('intersubject correlation'); title(sprintf('%s - %s: session-%d',hemiType,distType,ss));
        end
    case 'ALPHA:intersess'
        % intersession correlation
        distType = 'cosine';
        seqType = 'all'; % all / trained / untrained
        nTessel = 362;
        hemiType = 'CortexLeft';
        vararginoptions(varargin,{'sessN','distType','seqType','sn','nTessel','hemiType'});
        A = load(fullfile(connectDir,'alpha',sprintf('alpha_%s_nTessel-%d_seq-%s_%s',hemiType,nTessel,seqType,distType)));

        TT = [];
        allSess = unique(A.sessN);
        ind = indicatorMatrix('allpairs',(allSess'));
        figure
        for i=1:size(ind,1)
            sess1 = allSess(ind(i,:)==1);
            sess2 = allSess(ind(i,:)==-1);
            for s=sn
                a = getrow(A,A.sn==s & ismember(A.sessN,[sess1,sess2]));
                T.corr_sess = rsa_vectorizeRDM(corr(a.dist'));
                T.sn = s;
                T.sess1 = sess1;
                T.sess2 = sess2;
                TT = addstruct(TT,T);
            end
            subplot(1,size(ind,1),i)
            histogram(TT.corr_sess(TT.sess1==sess1 & TT.sess2==sess2));
            hold on; drawline(median(TT.corr_sess(TT.sess1==sess1 & TT.sess2==sess2)),'color',[1 0 0]);
            xlabel('Correlation');
            title(sprintf('%s - %s: interses: %d-%d',hemiType,distType,sess1,sess2));
        end  
    case 'ALPHA:cross_sess_calculate'
        % calculate alpha across sessions
        sessN = 1:4;
        hemi = 1;
        betaChoice='multi';
        distType = 'cosineW';
        crossType = 'splithalf'; % all or splithalf
        tesselSelect = 0;
        seqType = 'all'; % all / trained / untrained
        nTessel = 642;
        vararginoptions(varargin,{'sessN','sn','hemi','betaChoice','distType','seqType','nTessel','crossType','tesselSelect'});
        
        if length(hemi)==2
            T_select=load(fullfile(betaDir,'group',sprintf('stats_tesselsWB_%d_%sPW_sess1',nTessel,betaChoice)));
            for i=1:2
                if tesselSelect
                   %tess1 = sml_connect('TESSEL:select','hemi',i,'betaChoice',betaChoice,'sessN',sessN,'nTessel',nTessel);
                    tess1 = sml_connect('TESSEL:select_old','hemi',i,'betaChoice',betaChoice,'sessN',sessN,'nTessel',nTessel);
                    %tess{i} = unique(T_select.region(ismember(T_select.regType,tess1)&T_select.regSide==i));
                    tess{i} = tess1';
                else
                    tess{i} = unique(T_select.region(T_select.regSide==i));
                end
            end
            tessels = [tess{1};tess{2}];
            hemiType = 'CortexBoth';
        else
            if tesselSelect
                %tessels = sml_connect('TESSEL:select','hemi',hemi,'betaChoice',betaChoice,'sessN',sessN,'nTessel',nTessel);
                tessels = sml_connect('TESSEL:select_old','hemi',hemi,'betaChoice',betaChoice,'sessN',sessN,'nTessel',nTessel);
            else
                T_select=load(fullfile(betaDir,'group',sprintf('stats_tesselsWB_%d_%sPW_sess1',nTessel,betaChoice)));
                tessels = unique(T_select.region(T_select.regSide==hemi));
            end
            hemiType = hemName{hemi};
        end
        switch seqType
            case 'all'
                metric = 'RDM';
            case 'trained'
                metric = 'RDM_train';
            case 'untrained'
                metric = 'RDM_untrain';
        end
        switch crossType
            case 'alldata'
                metric_raw = metric;
            case 'splithalf'
                switch seqType
                    case 'all'
                        metric_raw = {'partA_all','partB_all'};
                    case 'trained'
                        metric_raw = {'partA_train','partB_train'};
                    case 'untrained'
                        metric_raw = {'partA_untrain','partB_untrain'};
                end
        end
        TT = []; TT2 = []; AA = [];
        for ss=sessN
            if strcmp(crossType,'alldata')
                T=load(fullfile(betaDir,'group',sprintf('stats_tesselsWB_%d_%sPW_sess%d',nTessel,betaChoice,ss)));
                T.sessN=ones(size(T.SN))*ss;
                T2 = T;
            elseif strcmp(crossType,'splithalf')
                V=load(fullfile(betaDir,'group',sprintf('betas_partition_tesselsWB_%d_%sPW_sess%d',nTessel,betaChoice,ss)));
                T.(metric) = [V.(metric_raw{1});V.(metric_raw{2})];
                T.sessN     = [V.sessN; V.sessN];
                T.SN        = [V.sn; V.sn];
                T.region    = [V.roi; V.roi];
                T2.SN       = V.sn;
                T2.region   = V.roi;
                T2.regSide  = V.regSide;
                T2.regType  = V.regType;
                T2.sessN    = V.sessN;
            elseif strcmp(crossType,'doubleCrossval')
                T=load(fullfile(betaDir,'group',sprintf('betas_tesselsWB_%d_sess%d',nTessel,ss)));
                T.sessN     = ones(size(T.SN))*ss;
                T2 = T;
            end
            TT = addstruct(TT,T);
            TT2 = addstruct(TT2,T2);
        end
        % only select tessels with data present
        TT = getrow(TT,ismember(TT.region,tessels));
        TT2 = getrow(TT2,ismember(TT2.region,tessels));
        for s=sn
            t   = getrow(TT,TT.SN==s);
            t2  = getrow(TT2,TT2.SN==s);
            if any(strcmp(crossType,{'alldata','splithalf'}))
                A = sml_connect('CALC:distance',t.(metric),distType);
            elseif strcmp(crossType,'doubleCrossval')
                condVec = repmat(1:12,1,8)';
                nPart = 8;
                switch seqType
                    case 'all'
                        nCond = numel(unique(condVec));
                        data = cell(size(t.SN,1),1);
                        for rr=1:size(t.SN,1)
                            data{rr,:} = t.betaW{rr}(1:size(condVec,1),:);
                        end
                        data = sml_connect('HOUSEKEEPING:removeRunMean',data,nPart,nCond);
                        tElapsed = tic;
                        lCKA = doubleCrossval_lcka_cv(data,nPart,nCond);
                        toc(tElapsed);
                        A.dist = rsa_vectorizeRDM(lCKA.ccv);
                    case 'trained'
                        condVec = condVec(condVec<7,:);
                        nCond = numel(unique(condVec));
                        for rr=1:size(t.SN,1)
                            data{rr,:} = t.betaW{rr}(condVec>0,:);
                        end
                        data = sml_connect('HOUSEKEEPING:removeRunMean',data,nPart,nCond);
                        lCKA = doubleCrossval_lcka_multiReg(data,nPart,nCond);
                        A.dist = rsa_vectorizeRDM(lCKA.ccv);
                    case 'untrained'
                        condVec = condVec(condVec>6,:);
                        nCond = numel(unique(condVec));
                        for rr=1:size(t.SN,1)
                            data{rr,:} = t.betaW{rr}(condVec>0,:);
                        end
                        data = sml_connect('HOUSEKEEPING:removeRunMean',data,nPart,nCond);
                        lCKA = doubleCrossval_lcka_multiReg(data,nPart,nCond);
                        A.dist = rsa_vectorizeRDM(lCKA.ccv);
                end
            end
            A.similarity    = sml_connect('CALC:gaussian_similarity',A.dist);
            if any(strcmp(crossType,{'alldata','doubleCrossval'}))
                A.confidence = ones(1,size(squareform(A.dist),1));
                A.tessel        = (t.region)';
                A.regSide       = (t.regSide)';
                A.regType       = (t.regType)';
                A.sessN         = (t.sessN)';
            elseif strcmp(crossType,'splithalf')
                nDist = size(squareform(A.dist),1);
                % take only off-diagonal
                D = squareform(A.dist);
                S = squareform(A.similarity);
                ll1 = squareform(A.l1);
                ll2 = squareform(A.l2);
                A.dist = rsa_vectorizeRDM(D(nDist/2+1:end,1:nDist/2));
                A.similarity = rsa_vectorizeRDM(S(nDist/2+1:end,1:nDist/2));
                A.l1 = rsa_vectorizeRDM(ll1(nDist/2+1:end,1:nDist/2));
                A.l2 = rsa_vectorizeRDM(ll2(nDist/2+1:end,1:nDist/2));
                A.confidence = (diag(squareform(A.dist)))';
                A.tessel        = (t2.region)';
                A.regSide       = (t2.regSide)';
                A.regType       = (t2.regType)';
                A.sessN         = (t2.sessN)';
            end
            A.sn            = s;
            AA = addstruct(AA,A);
            fprintf('\nDone %s\n',subj_name{s});
        end
        save(fullfile(connectDir,'alpha',sprintf('alpha_crossSess-%s_%s_nTessel-%d_seq-%s_%s',crossType,hemiType,nTessel,seqType,distType)),'-struct','AA');       
        fprintf('Done %s alpha - %s RDM on %d tessels, %s.\n\n\n',distType,seqType,nTessel,hemiType);
    
    case 'CONNECT:group_louvain' % depreciated
        % make group connectivity analysis (average alpha across subjects)
        % tesselation (162 / 362)
        % first try with predefined gamma and omega
        % otherwise estimate across many of them
        % use iterated_genlouvain
        sessN = 1:3; % options 1-3, 1-4, 3-4
        seqType = 'all';
        distType = 'cosine';
        gamma = 1.1; % modularity
        omega = 0.05; % determines the strength of the inter-slice connection
        
        vararginoptions(varargin,{'seqType','distType','sessN','gamma','omega'});
        
        W = load(fullfile(connectDir,sprintf('alpha_seq-%s_%s',seqType,distType)));
        for ss=sessN
            D = squareform(mean(W.dist(W.sessN==ss,:),1));
            % make connectivity matrix A from distance matrix D
            thres = quantile(rsa_vectorizeRDM(D),0.05);
            A{ss}     = exp(-D.^2 ./ (2*thres^2));
        end
        N=length(A{1});
        T=length(A);
        
        B = multiord(A,gamma,omega); % multiord.m computes multilayer
        % modularity matrix B with homogeneous ordinal interlayer coupling w and
        % a Newman-Girvan null model on each layer (see more detail in
        % documentation of multiord.m)  
        PP = @(S) postprocess_ordinal_multilayer(S,T); % define postprocessing
        % function handle that increases multilayer modularity without changing
        % intralayer partitions in an ordered multilayer networks (see
        % postprocess_ordinal_multilayer.m for more detail)
        [S,Q,n_it] = iterated_genlouvain(B,10000,1,1,'moverandw',[], PP);
        S = reshape(S, N, T);
    case 'CONNECT:group'
        nTessel     = 642; % 162 or 362
        seqType     = 'all'; % all, trained, untrained
        crossType   = 'doubleCrossval'; % alldata or splithalf
        distType    = 'cosineW'; % distance type
        maxClust    = 10;
        hemiType    = 'CortexBoth'; % CortexBoth, CortexLeft
        vararginoptions(varargin,{'sn','sessN','nTessel','sessType','seqType','distType','crossType','maxClust','seqType','hemiType'});
        
        T=load(fullfile(connectDir,'alpha',sprintf('alpha_crossSess-%s_%s_nTessel-%d_seq-%s_%s',crossType,hemiType,nTessel,seqType,distType)));       
        
        % make a group mean
        A = squareform(nanmean(T.similarity));
        [C,L,U]=SpectralClustering(A,maxClust,3);
       % Clust = kmeans(T1.A,maxClust);
       % T1.community = community_louvain(T1.A);
        T1.cluster   = C;
        T1.laplace   = L;
        T1.eigenv    = U;
        T1.A        = squareform(mean(T.similarity)); % similarity matrix
        % here now add which regions / sessions nodes correspond to
        t1 = getrow(T,T.sn==5);
        T1.tessel   = t1.tessel';
        T1.regSide  = t1.regSide';
        T1.regType  = t1.regType';
        T1.sessN    = t1.sessN';
        T1.confidence = (mean(T.confidence))';
        
        save(fullfile(connectDir,'clusters',sprintf('%s_nTessel-%d_seq-%s_%s_cluster-%d_crossSess-%s',hemiType,nTessel,seqType,distType,maxClust,crossType)),'-struct','T1');
        fprintf('Done %s %d clusters - %s RDM on %d tessels, %s, %s.\n\n\n',distType,maxClust,seqType,nTessel,crossType,hemiType);
    case 'CONNECT:cluster_dendrogram'
        % here plot the clusters and dendrogram
        nTessel     = 642; % 162 or 362
        seqType     = 'all'; % all, trained, untrained
        crossType   = 'doubleCrossval'; % all or splithalf
        distType    = 'cosineW'; % distance type
        hemiType    = 'CortexBoth';
        nCluster    = 10;
        figOn       = 0;
        vararginoptions(varargin,{'seqType','distType','crossType','nTessel','seqType','nCluster','hemiType','figOn'});
     
        T = load(fullfile(connectDir,'clusters',sprintf('%s_nTessel-%d_seq-%s_%s_cluster-%d_crossSess-%s',hemiType,nTessel,seqType,distType,nCluster,crossType)));
       % Evaluate similarity graph (W)
        [Csort,idx]=sort(T.cluster);
        A = full(T.A);
        Asort=A(idx,:);
        Asort=Asort(:,idx);
        Ncluster = numel(unique(T.cluster));
        
        % Laplacian eigenvector
        Lsort = T.laplace(idx,:); Lsort=Lsort(:,idx);
        Usort = T.eigenv(idx,:);
        
        % Merge similarities according to clusters
        Asort(logical(eye(size(Asort)))) = NaN;
        for i=1:Ncluster
            for j=1:Ncluster
                idxi=Csort==i;
                idxj=Csort==j;
                rWsrot(i,j) = nanmean(vec(Asort(idxi,idxj)));
                rLsort(i,j) = nanmean(vec(Lsort(idxi,idxj)));
            end
            rUsort(i,:) = nanmean(Usort(idxi,:),1);
        end
        % get dendogram linkage
        Z = linkage(real(rUsort),'ward','euclidean');
        % check if the colour assignment already exists
        colFile = fullfile(connectDir,'clusters',sprintf('myColors_%s_nTessel-%d_seq-%s_%s_cluster-%d_crossSess-%s.txt',hemiType,nTessel,seqType,distType,nCluster,crossType));
        if exist(colFile)
            Col=dlmread(colFile); % read in colours
        else
            Col = sml_connect('CONNECT:assign_color',Z,size(rUsort,1));
            dlmwrite(fullfile(connectDir,'clusters',sprintf('myColors_%s_nTessel-%d_seq-%s_%s_cluster-%d_crossSess-%s.txt',hemiType,nTessel,seqType,distType,nCluster,crossType)),Col);
        end
        if figOn % visually inspect the sorted connectivity and dendrogram 
            figure
            subplot(121)
            imagesc(Asort); colorbar; title(sprintf('%s %s %s RDM - %d clusters, %s',hemiType,distType,seqType,nCluster,crossType));
            subplot(122)
            [h,~,per] = dendrogram(Z,'Orientation','right');
            hold on;
            xlim=get(gca,'xlim');
            for i=1:numel(per)
                plot(xlim(1),i,'o','markersize',8,...
                    'markeredgecolor','k',...
                    'markerfacecolor',Col(per(i),:));
                hold on;
            end
            for i=1:numel(h)
                set(h(i),'color',[0 0 0],'linewidth',1);
            end
        end
    case 'CONNECT:assign_color'
        % assign the colors for plotting, make a dendrogram
        Z = varargin{1};
        sizeZ = varargin{2};
       
        Col = colorDendrogram(Z,sizeZ,'colorspace','rgb','order',[1,2,3],'fig',0,'weight',1);
        varargout{1} = Col;
    case 'CONNECT:surface'
        % here project to surface
        nTessel     = 642; % 162 or 362
        seqType     = 'all'; % all, trained, untrained
        crossType   = 'splithalf'; % all or splithalf
        distType    = 'cosineW'; % distance type
        hemiType    = 'CortexBoth';
        nCluster    = 10;
        vararginoptions(varargin,{'sn','sessN','parcelType','seqType','distType','crossType','nTessel','seqType','nCluster','hemiType'});
        
        switch hemiType
            case 'CortexBoth'
                hemi = [1:2];
            case 'CortexLeft'
                hemi = 1;
            case 'CortexRight'
                hemi = 2;
        end
        T = load(fullfile(connectDir,'clusters',sprintf('%s_nTessel-%d_seq-%s_%s_cluster-%d_crossSess-%s',hemiType,nTessel,seqType,distType,nCluster,crossType)));
        colFile = fullfile(connectDir,'clusters',sprintf('myColors_%s_nTessel-%d_seq-%s_%s_cluster-%d_crossSess-%s.txt',hemiType,nTessel,seqType,distType,nCluster,crossType));
        if exist(colFile)
            Col=dlmread(colFile); % read in colours
        else
            error('The corresponding colour file has not been created. Run "CONNECT:plot_cluster_dendrogram"!\n');
        end
        labelRGBA = [Col ones(size(Col,1),1)];
        labelRGBA = [zeros(1,4);labelRGBA]; % add for 0s
        for h=hemi
            G = gifti(fullfile(wbDir,'connect_FS_LR_164',sprintf('Icosahedron-%d.164k.%s.label.gii',nTessel,hemI{h})));
            data = zeros(size(G.cdata,1),numel(sessN));
            %nReg = numel(unique(G.cdata))-1; % exclude 0 - medial wall
            for ss=sessN
                S = getrow(T,T.regSide==h & T.sessN==ss);
                clustS = unique(S.cluster);
                for i=clustS'
                    SS=getrow(S,S.cluster==i);
                    % data from that hemi, tessel, per cluster
                    indx = ismember(G.cdata,SS.regType);    % indicate the right tessels
                    data(indx,ss) = i;
                    labelName{i} = sprintf('cluster-%d',i);
                end
                columnName{ss} = sprintf('sess-%d',ss);
            end
            G=surf_makeLabelGifti(data,'anatomicalStruct',hemName{h},'columnName',columnName,'labelRGBA',labelRGBA);
            outfile = fullfile(wbDir,'connect_FS_LR_164',sprintf('%s.clusters-%d_%s_nTessel-%d_seq-%s_%s_crossSess-%s.label.gii',hemI{h},nCluster,hemiType,nTessel,seqType,distType,crossType));
            % here save
            save(G,outfile);
        end
              
    case 'PCM:runModels'
        % here run first finger, all fingers, sequences, sequence type
        % per extracted cluster
        % option 1) ignore which session it belongs to
        % option 2) do it per session
        
    case 'TESSEL:select_old'
        % here select tessels where all subjects have data
        sessN=1;
        hemi=1;
        betaChoice='multi';
        nTessel = 162;
        vararginoptions(varargin,{'sessN','betaChoice','hemi','nTessel'});
        
        for ss=1:3       
            indx=[];
            T=load(fullfile(betaDir,'group',sprintf('stats_tesselsWB_%d_%sPW_sess%d',nTessel,betaChoice,ss)));
            T = getrow(T,ismember(T.regSide,hemi));
            for r=unique(T.region)'
                T1 = getrow(T,T.region==r & T.regSide==hemi);
                if ~sum(any(isnan(T1.RDM))) % make sure all subjects have data there
                    %[~,p]=ttest(T1.dist_all,[],2,'onesample');
                    %if p<.01 % maybe still not the best way to select tessels
                        indx = [indx,r];
                    %end
                end
            end
            if ss==sessN(1)
                choice=indx;
            else
                choice=unique([indx,choice]);
            end
        end
        varargout={choice};
    case 'TESSEL:select'
        % here select based on the distance mask
        hemi=1;
        nTessel = 362;
        vararginoptions(varargin,{'hemi','nTessel','betaChoice','sessN'});
        % load in distances and icosahedron
        I = gifti(fullfile(wbDir,'FS_LR_164',sprintf('Icosahedron-%d.164k.%s.label.gii',nTessel,hemI{hemi})));
        G = gifti(fullfile(wbDir,'FS_LR_164',sprintf('s%s.summary_dist.func.gii',hemI{hemi})));
        %maskDist = G.cdata(:,1)>.005; % all vertices where distances (all) in sess1 > .004
        %maskDist = G.cdata(:,17)>3.4; % all vertices where distances (all) in sess1 significant with p<.001
        maskDist = G.cdata(:,1)>.001 & G.cdata(:,17)>3.4;
        tessels = unique(I.cdata)';
        tessels = tessels(tessels~=0); % exclude 0 - medial wall
        choice = [];
        for i=tessels
            numAll = sum(I.cdata==i);
            distPres = sum(maskDist(I.cdata==i)); % how many times significant distance (dist presence)
            if distPres>(numAll*.4)
                choice = [choice,i];
            end
        end
        varargout={double(choice)};
    case 'CALC:distance'
        % calculate correlation / cosine distances
            rdm = varargin{1};
            distType = varargin{2};
            % calculate distance metric from the input
            % input: N x D matrix (N - number of RDMs; D - distance pairs)
            % dist types: 'correlation' or 'cosine'
            % output: structure D, with fields:
            %   - D: pairwise distances of RDMs
            %   - l1: indicator which rdm / layer taken as first
            %   - l2: indicator which rdm ; layer taken as second
            switch (distType)
                case 'correlation'
                    % additional step for correlation - first remove the mean
                    rdm  = bsxfun(@minus,rdm,mean(rdm,2));
            end
            nRDM = size(rdm,1);
            rdm  = normalizeX(rdm);
            tmpR  = rdm*rdm'; % correlation across RDMs
            D.dist = 1-rsa_vectorizeRDM(tmpR); % distances
            N=size(squareform(D.dist),1);
            A=squareform(pdist([1:N]')); % pairs
            [i2,i1]=find(tril(A,-1));
            D.l1 = i1';
            D.l2 = i2';
           % this is equivalent to indicatormatrix
%            ind=indicatorMatrix('allpairs',(1:nRDM));
%            % determine each element of the pair
%            [~,D.l1]=find(ind==1); D.l1=D.l1';
%            i=ismember(ind,-1);
%            D.l2 = sum(cumprod(i==0,2),2)+1; D.l2=D.l2';
            %D.dist = rsa_vectorizeRDM(tmpR)'; % distances
            if strcmp(distType,'cosineW')
                if size(rdm,2)==66
                    nCond = 12;
                elseif size(rdm,2)==15
                    nCond = 6;
                else
                    error('Wrong size of inputs!\n');
                end
                varD = rsa_varianceLDC(zeros(nCond),indicatorMatrix('allpairs',1:nCond),eye(nCond),8,nRDM); % Get the variance
                tmpR = sml_connect('CALC:distance_covWeighted',rdm,varD);
                D.dist = 1-rsa_vectorizeRDM(tmpR); % distances
            end
            varargout{1}=D;     
    case 'CALC:distance_covWeighted'
        % calculated weighted cosine similarity measure
        A = varargin{1}; % N x q vector (distance)
        Sig = varargin{2}; % covariance structure - q x q (under H0)
        % N*N connectivity matrix of weighted inner product of RDMs (cosineW)
        %        C = indicatorMatrix('allpairs',1:size(A,1));
        %idx1 = repmat((1:size(A,1))',1,size(A,1));
        %idx2 = repmat(1:size(A,1),size(A,1),1);
        wCos = zeros(size(A,1));
        [V,L]=eig(Sig);
        l=diag(L);
        sq = V*bsxfun(@rdivide,V',sqrt(l)); % Slightly faster than sq = V*diag(1./sqrt(l))*V';
        
        for i1=1:size(A,1)
            for i2=i1:size(A,1)
                %ind1 = idx1(i1);
                %ind2 = idx2(i2);
                w1 = A(i1,:)*sq;
                w2 = A(i2,:)*sq;
                w1=bsxfun(@rdivide,w1,sqrt(sum(w1.^2,2)));
                w2=bsxfun(@rdivide,w2,sqrt(sum(w2.^2,2)));
                r=w1*w2';
                %wCos(ind1,ind2) = r; wCos(ind2,ind1) = r;
                wCos(i1,i2) = r; wCos(i2,i1) = r;
            end
        end
%         for i=1:size(C,1)
%             ind1 = C(i,:)==1;
%             ind2 = C(i,:)==-1;
%             w1 = A(ind1,:)*sq;
%             w2 = A(ind2,:)*sq;
%             w1=bsxfun(@rdivide,w1,sqrt(sum(w1.^2,2)));
%             w2=bsxfun(@rdivide,w2,sqrt(sum(w2.^2,2)));
%             r=w1*w2';
%             wCos(ind1,ind2) = r; wCos(ind2,ind1) = r;
%         end
        varargout{1} = wCos;
    case 'CALC:gaussian_similarity'
        % here calculate gaussian similarity (from distances)
        % create a matrix of similarity across regions
        A       = varargin{1};
        thres   = quantile(A,0.05);
        W       = rsa_squareRDM(exp(-A.^2 ./ (2*thres^2)));
        varargout{1} = rsa_vectorizeRDM(W);
    case 'CALC:noiseCeiling'
        % here calculate the upper and lower noise ceilings
        data = varargin{1}; % assumes data comes in dimensions subj x var
        sn = varargin{2};
        noise_high = zeros(length(sn),1);
        noise_low = zeros(length(sn),1);
        for s=1:length(sn)
            noise_high(s) = corr(data(s,:)',mean(data(rem(sn,2)==rem(sn(s),2),:),1)');
            noise_low(s) = corr(data(s,:)',mean(data(~ismember(sn,sn(s)) & rem(sn,2)==rem(sn(s),2),:),1)');
        end
        varargout{1}=noise_high;
        varargout{2}=noise_low;
    case 'CALC:double_crossval'
        % here execute doubleCrossval (lcka) for all pairs of regions
        data = varargin{1};
        nPart = varargin{2};
        nCond = varargin{3};
        nReg = size(data,1);
        C = indicatorMatrix('allpairs',1:nReg);
        TT = [];
        fprintf('Done:\n');
        tElapsed = tic;
        for i=1:size(C,1)
            l1 = C(i,:)==1;
            l2 = C(i,:)==-1;
            T = doubleCrossval_lcka([data(l1) data(l2)],nPart,nCond);
            T.l1 = find(l1);
            T.l2 = find(l2);
            TT = addstruct(TT,T);
            fprintf('%d/%d.\n',i,size(C,1));
        end
        toc(tElapsed);
        %D{1} = rsa_squareRDM(TT.ncv');
        %D{2} = rsa_squareRDM(TT.cv');
        %D{3} = rsa_squareRDM(TT.ccv');
        D.dist = TT.ccv';
        D.l1 = TT.l1';
        D.l2 = TT.l2';
        varargout{1}=D;
    
    case 'HOUSEKEEPING:removeRunMean'
        act = varargin{1};
        nPart = varargin{2};
        nCond = varargin{3};
        partVec = kron((1:nPart)',ones(nCond,1));
        numLayer = size(act,1);
        actN = cell(numLayer,1);
        for i=1:numLayer
            actN{i} = [];
            for j=1:max(partVec)
                actN{i} = [actN{i}; bsxfun(@minus,act{i}(partVec==j,:),mean(act{i}(partVec==j,:),1))];
            end
        end
        varargout{1}=actN;  
        
    case 'run_all'
        sT = {'trained','all'};
        nClust = [10,12,14];
        for st=1:2
            sml_connect('ALPHA:cross_sess_calculate','crossType','doubleCrossval','tesselSelect',1,'seqType',sT{st});
            sml_connect('ALPHA:cross_sess_calculate','crossType','doubleCrossval','tesselSelect',1,'hemi',1:2,'seqType',sT{st});
            for n=nClust
                sml_connect('CONNECT:group','seqType',sT{st},'nTessel',642,'hemiType','CortexBoth','maxClust',n,'crossType','doubleCrossval');
                sml_connect('CONNECT:group','seqType',sT{st},'nTessel',642,'hemiType','CortexLeft','maxClust',n,'crossType','doubleCrossval');
                sml_connect('CONNECT:cluster_dendrogram','seqType',sT{st},'nTessel',642,'hemiType','CortexBoth','nCluster',n,'crossType','doubleCrossval');
                sml_connect('CONNECT:cluster_dendrogram','seqType',sT{st},'nTessel',642,'hemiType','CortexLeft','nCluster',n,'crossType','doubleCrossval');
            end
            fprintf('***** Done all for %s *****\n\n\n',sT{st});
        end
        fprintf('***** Done all *****\n\n\n');
    case 'run_tessel'
        sml_connect('ALPHA:cross_sess_calculate','distType','cosineW','nTessel',642,'hemi',1:2);
        sml_connect('ALPHA:cross_sess_calculate','distType','cosineW','nTessel',642,'hemi',1,'seqType','trained');
        sml_connect('ALPHA:cross_sess_calculate','distType','cosineW','nTessel',642,'hemi',1:2,'seqType','trained');

    otherwise
        error('Wrong case!\n');
        
end
end