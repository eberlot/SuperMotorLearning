function varargout = sml1_pcm(what, varargin);
% SuperMoterLearning PCM modelling
% on real data.
baseDir         = '/Volumes/G_Thunderbolt/Yokoi_Research/data/SequenceLearning/sh_eva'; %

addpath('/Users/atsushiyokoi/Dropbox/Matlab/matlab/imaging/mva/pcm_develop/');

figDir = fullfile(baseDir, 'cluster_sess4','figures'); 

% Subjects with odd number practiced sequences 1-6
% Subjects with even number practiced sequences 7-12
% Prewhitened beta 1-6 are always for trained sequences 
% (7-12 are always for untrained)
% Always with digit cues

runEffect   = 'fixed';
singledigit = 'naturalstats'; % single digit structure for pcm. default is to use sh2 data (or 'naturalstats')
freqtype    = 'detailed'; % use training info for transition frequency or not ('detailed' uses training info)
fitAlgorithm= 'NR';
atlasA      = {'x'};
hem = {'lh','rh'};
hemName = {'LeftHem','RightHem'};
hemispheres = {'LeftHem','RightHem'};
interceptstrs = {'no_intercept','with_intercept'};
naturalstatstrs = {'no_naturalstats','with_naturalstats'};
flipstrs = {'','_flipgroup'};
prewBeta = fullfile(baseDir, 'cluster_sess4', 'clusterData.mat');


switch (what)
    case 'show_predictedRDM' % construct model and plot them
        group=1;
        useNaturalStat = 0;
        vararginoptions(varargin, {'group','useNaturalStat'});
        
        S = load(fullfile(baseDir, 'cluster_sess4', 'seqInfo.mat'));
        
        % rename
        S.modelName{1} = 'trained/untrained';
        % add two finger transition
        S.modelName{end+1} = 'TwoFingTransition';
        % modify feature
        
        % first finger
        if useNaturalStat
            sh1Dir         = '/Volumes/G_Thunderbolt/Yokoi_Research/data/SequenceLearning/sh1'; % external HDD 2
            load(fullfile(sh1Dir,'analyze','naturalstatisticmodel.mat'));
            singledigit = NatStats.G_cent;
            singledigit = (singledigit+singledigit')/2;
        else            
            singledigit = eye(5);
        end
        
        switch group
            case 1
            case 2
                S.M{3} = S.M{3}([7:12,1:6],:);
                S.M{4} = S.M{4}([7:12,1:6],:);
                S.Seq = S.Seq([7:12,1:6],:);
        end
        
        % chunk
        S.M{end} = blockdiag(S.M{end}(1:6,:), S.M{end}(7:12,:));
        % transition
        S.M{end+1} = TwoFingTransition(S.Seq);
        
        % calc RDM and show
        myFigure([40,15],'name', 'Predicted RDMs');
        con=indicatorMatrix('allpairs',[1:12]);
        for m=1:length(S.modelName)
            if m==3||m==4
                O=singledigit;
            else
                O=eye(size(S.M{m},2));
            end
            G(:,:,m) = S.M{m}*O*S.M{m}';
            D(:,:,m)=squareform(diag(con*G(:,:,m)*con'));
            subplot(2,length(S.M)/2,m);
            imagesc(D(:,:,m)); axis square;
            title(S.modelName{m});
        end;
        
        varargout = {S, G, D};
    case 'show_clusterRDM' % plot cluster-averaged empirical RDMs
        sess = 4;
        T=load(fullfile(baseDir, sprintf('cluster_sess%d', sess), 'tessel_RDMs.mat'));
        C=load(fullfile(baseDir, sprintf('cluster_sess%d', sess), 'clusterResults.mat'));
        
        % average RDM within each subject using C info
        %T=getrow(T, T.regSide==1); % discard right hemisphere for now
        T.cluster=zeros(size(T.roi));
        for i=1:length(C.roi)
            T.cluster(T.roi==C.roi(i)) = C.cluster(i);
        end
        T=tapply(T, {'sn','cluster',}, {'RDM_train','nanmean(x,1)','name','RDM_train'},...
            {'RDM_untrain','nanmean(x,1)','name','RDM_untrain'},...
            {'RDM_all','nanmean(x,1)','name','RDM_all'});
        
        % calculate mean distance between trained and untrained sequence
        I = blockdiag(true(6),true(6));
        for i=1:length(T.sn)
            D=rsa_squareRDM(T.RDM_all(i,:));
            T.meandist_tut(i,1) = nanmean(D(~I));
        end
        
        % show cluster-average RDM for both groups together
        for group=[1,0]
            Tg=getrow(T, mod(T.sn,2)==(group));
            figure('name',sprintf('cluster RDMs (G%d)',group+2*(group<1)));
            for c=unique(Tg.cluster)'
                subplot(2,5,c);
                D = nanmean(Tg.RDM_all(Tg.cluster==c,:),1);
                %D = nanmean(Tg.RDM_train(Tg.cluster==c,:),1);
                %D = nanmean(Tg.RDM_untrain(Tg.cluster==c,:),1);
                D = rsa_squareRDM(D);
                imagesc(D); axis square;
                title(sprintf('C%d',c));
            end
        end
        
        varargout = {};
    case 'pcm_makeDataDesign' % get data for trained sequences and save
        % ============================== %
        % load data (12 conditions x 8 runs (+intercepts))
        Data = load(prewBeta);
        Data.group(mod(Data.sn,2)==1,1)=1; % odd number (trained:1-6)
        Data.group(mod(Data.sn,2)==0,1)=2; % even number (trained:7-12)
        for i=1:length(Data.sn)
            Data.beta{i}(end-7:end,:) = []; % discard intercept
            idx = [];
            switch Data.group(i)
                case 1
                    for j=1:6 % trained sequences are always condition 1-6
                        idx=[idx, [j:12:96]];
                    end
                case 2
                    for j=1:6 % trained sequences are always condition 1-6
                        idx=[idx, [j:12:96]];
                    end
                otherwise
                    error();
            end
            idx=sort(idx);
            Data.beta_train{i,1}=Data.beta{i}(idx,:);
            Data.beta_untrain{i,1}=Data.beta{i}(setdiff([1:96],idx),:);
        end
        Data.beta_all=Data.beta; % all conditions
        Data=rmfield(Data, 'beta'); % save memory
        
        % create design
        condition = repmat([1:6]',8,1);
        %condition = [condition; zeros(8,1)];
        partition = kron([1:8]',ones(6,1));
        %partition = [partition; [1:8]'];
        Design.condVec=condition;
        Design.partVec=partition;
        
        save(fullfile(baseDir, 'cluster_sess4', 'prewhBeta.mat'), 'Data', 'Design', '-v7.3');
    case 'pcm_fitGroupCrossval_train' % use only trained sequences
        sess = 4;
        intercept = 0;
        useNaturalStat = 0;
        isSave = 0;
        flipgroup = 0;
        vararginoptions(varargin, {'sess','intercept','useNaturalStat', 'isSave','flipgroup'});
        
        interceptstr = interceptstrs{intercept+1};
        nstatstr = naturalstatstrs{useNaturalStat+1};
        flipstr = flipstrs{flipgroup+1};
        
        % ============================== %
        % load prewhitened beta and design
        load(fullfile(baseDir, 'cluster_sess4', 'prewhBeta.mat'));
        
        % ============================== %
        % load model
        [S,G] = sml1_pcm('show_predictedRDM', 'useNaturalStat', useNaturalStat);
        
        % run pcm at each cluster
        for g=[1,2]
            % Make model component for both groups
            c=0;M=[];MF = [];
            for m=[3,4,6,5,2]; % firstFing, allFing, fingeTran, Chunk, Seq
                c=c+1;
                M{c}.type = 'component';
                if flipgroup
                    if g==2
                        M{c}.Gc = G(1:6,1:6,m);
                    elseif g==1
                        M{c}.Gc = G(7:12,7:12,m);
                    end
                else
                    if g==1
                        M{c}.Gc = G(1:6,1:6,m);
                    elseif g==2
                        M{c}.Gc = G(7:12,7:12,m);
                    end
                end
                M{c}.name = S.modelName{m};
                M{c}.numGparams = 1;
            end
            
            if intercept~=0 % Add intercept if required            
                M{end+1}.type = 'component';
                M{end}.Gc = ones(6); % explains variance (2nd moment) but not for distance
                M{end}.name = 'intercept';
                M{end}.numGparams = 1;
                
                % Normalise the model components to make them comparable size
                for i=1:numel(M)
                    M{i}.Gc=(M{i}.Gc+M{i}.Gc')/2;
                    M{i}.Gc=M{i}.Gc./mean(diag(M{i}.Gc));
                end;
                
                % Make the whole model family
                [MF,compI]=pcm_constructModelFamily(M, 'alwaysInclude', numel(M)); % intercept should always be included                
            else
                % Normalise the model components to make them comparable size
                for i=1:numel(M)
                    M{i}.Gc=(M{i}.Gc+M{i}.Gc')/2;
                    M{i}.Gc=M{i}.Gc./mean(diag(M{i}.Gc));
                end;
                
                % Make the whole model family
                [MF,compI]=pcm_constructModelFamily(M);
            end
            % Give initial parameters
            for i=1:numel(MF)
                MF{i}.theta0 = log(ones(MF{i}.numGparams,1));
            end
            
            for c=unique(Data.cluster)';
                D=getrow(Data, Data.cluster==c&Data.group==g);
                Y = D.beta_train;
                
                [T{c,g}, Tcv{c,g}, thetacv{c,g}, compidx] = ...
                    pcm_clusterROI_component(Y, Design.partVec, Design.condVec, MF, compI,...
                    'fitAlgorithm',fitAlgorithm,'runEffect',runEffect);
            end
        end
        
        % save&show result
        fname = sprintf('pcm_cluster_train_%s_%s%s.mat', interceptstr,nstatstr,flipstr);
        if isSave
            save(fullfile(baseDir, 'cluster_sess4', fname), 'T','Tcv','thetacv','compidx','M','MF')
        end
        sml1_pcm('pcm_show_cluster_result', Tcv, M, compI, fname, isSave);
        varargout = {};
    case 'pcm_fitGroupCrossval_untrain' % use only untrained sequences
        intercept = 0;
        useNaturalStat = 0;
        isSave = 0;
        flipgroup = 0;
        vararginoptions(varargin, {'sess','intercept','useNaturalStat', 'isSave','flipgroup'});
        
        interceptstr = interceptstrs{intercept+1};
        nstatstr = naturalstatstrs{useNaturalStat+1};
        flipstr = flipstrs{flipgroup+1};
        
        % ============================== %
        % load prewhitened beta and design
        load(fullfile(baseDir, 'cluster_sess4', 'prewhBeta.mat'));
        
        % ============================== %
        % load model
        [S,G] = sml1_pcm('show_predictedRDM', 'useNaturalStat', useNaturalStat);
        
        % run pcm at each cluster
        for g=[1,2]
            % Make model component for both groups
            c=0;M=[];
            for m=[3,4,6,5,2]; % firstFing, allFing, fingeTran, Chunk, Seq
                c=c+1;
                M{c}.type = 'component';
                if flipgroup==0
                    if g==2
                        M{c}.Gc = G(1:6,1:6,m);
                    elseif g==1
                        M{c}.Gc = G(7:12,7:12,m);
                    end
                else
                    if g==1
                        M{c}.Gc = G(1:6,1:6,m);
                    elseif g==2
                        M{c}.Gc = G(7:12,7:12,m);
                    end
                end
                M{c}.name = S.modelName{m};
                M{c}.numGparams = 1;
            end
            
            if intercept~=0 % Add intercept if required
                M{end+1}.type = 'component';
                M{end}.Gc = ones(6);
                M{end}.name = 'intercept';
                M{end}.numGparams = 1;
                
                % Normalise the model components to make them comparable size
                for i=1:numel(M)
                    M{i}.Gc=(M{i}.Gc+M{i}.Gc')/2;
                    M{i}.Gc=M{i}.Gc./mean(diag(M{i}.Gc));
                end;
                
                % Make the whole model family
                [MF,compI]=pcm_constructModelFamily(M, 'alwaysInclude', numel(M));
            else
                % Normalise the model components to make them comparable size
                for i=1:numel(M)
                    M{i}.Gc=(M{i}.Gc+M{i}.Gc')/2;
                    M{i}.Gc=M{i}.Gc./mean(diag(M{i}.Gc));
                end;
                
                % Make the whole model family
                [MF,compI]=pcm_constructModelFamily(M);
            end
            % give initial parameters
            for i=1:numel(MF)
                MF{i}.theta0 = log(ones(MF{i}.numGparams,1));
            end
            
            for c=unique(Data.cluster)';
                D=getrow(Data, Data.cluster==c&Data.group==g);
                Y = D.beta_untrain;
                
                [T{c,g}, Tcv{c,g}, thetacv{c,g}, compidx] = ...
                    pcm_clusterROI_component(Y, Design.partVec, Design.condVec, MF, compI,...
                    'fitAlgorithm',fitAlgorithm,'runEffect',runEffect);
            end
        end
        
        % save&show result
        fname = sprintf('pcm_cluster_untrain_%s_%s%s.mat',interceptstr,nstatstr,flipstr);
        if isSave
            save(fullfile(baseDir, 'cluster_sess4', fname), 'T','Tcv','thetacv','compidx','M','MF')
        end
        sml1_pcm('pcm_show_cluster_result', Tcv, M, compI, fname, isSave);
        varargout = {};
    
    case 'pcm_fitGroupCrossval_train_allgroup' % use only trained sequences, but all subjects
        intercept = 0;
        useNaturalStat = 0;
        isSave = 0;
        flipgroup=0;
        vararginoptions(varargin, {'sess','intercept','useNaturalStat','isSave','flipgroup'});
        
        interceptstr = interceptstrs{intercept+1};
        nstatstr = naturalstatstrs{useNaturalStat+1};
        flipstr = flipstrs{flipgroup+1};
        
        % ============================== %
        % load model
        [S,G] = sml1_pcm('show_predictedRDM','useNaturalStat', useNaturalStat);
        Gc=G;
        % Make model component for both groups
        c=0; M=[];
        for m=[3,4,5,6];% firstFing, allFing, chunk, fingeTran
            c=c+1;
            if flipgroup
                M{c}.Gc = blockdiag(Gc(7:12,7:12,m), Gc(1:6,1:6,m));
            else
                M{c}.Gc = blockdiag(Gc(1:6,1:6,m), Gc(7:12,7:12,m));
            end
            M{c}.name = S.modelName{m};
        end
        M{end+1}.type = 'component';
        M{end}.Gc = eye(12);
        M{end}.name = 'sequence';
        
        if intercept==1
            M{end+1}.type = 'component';
            M{end}.Gc = blockdiag(ones(6),ones(6));
            M{end}.name = 'intercept';
            
            % Normalise the model components to make them comparable size
            for i=1:numel(M)
                M{i}.type = 'component';
                M{i}.numGparams = 1;
                M{i}.Gc=(M{i}.Gc+M{i}.Gc')/2;
                M{i}.Gc=M{i}.Gc./mean(diag(M{i}.Gc));
            end;
            
            % Make the whole model family
            [MF,compI]=pcm_constructModelFamily(M, 'alwaysInclude', numel(M));
        else
            % Normalise the model components to make them comparable size
            for i=1:numel(M)
                M{i}.type = 'component';
                M{i}.numGparams = 1;
                M{i}.Gc=(M{i}.Gc+M{i}.Gc')/2;
                M{i}.Gc=M{i}.Gc./mean(diag(M{i}.Gc));
            end;
            
            % Make the whole model family
            [MF,compI]=pcm_constructModelFamily(M);
        end
        
        % Give initial parameters
        for i=1:numel(MF)
            MF{i}.theta0 = log(ones(MF{i}.numGparams,1));
        end
        
        % ====
        % run pcm
        % ============================== %        
        % load prewhitened beta and design
        load(fullfile(baseDir, 'cluster_sess4', 'prewhBeta.mat'));                
        
        % Define new Z
        D = getrow(Data, Data.cluster==1);
        Z1 = indicatorMatrix('identity', repmat([1:6]',8,1));
        Z0 = zeros(size(Z1));
        for s=1:length(D.cluster)
            switch D.group(s)
                case 1
                    Z{s} = [Z1, Z0];
                case 2
                    Z{s} = [Z0, Z1];
            end
            partVec{s} = kron([1:8]',ones(6,1));
        end
        
        for c=unique(Data.cluster)';
            D=getrow(Data, Data.cluster==c);
            Y = D.beta_train;
                        
            [T{c,1}, Tcv{c,1}, thetacv{c,1}, compidx] = ...
                pcm_clusterROI_component(Y, partVec, Z, MF, compI,...
                'fitAlgorithm',fitAlgorithm,'runEffect',runEffect);
        end
        
        % save&show result
        fname = sprintf('pcm_cluster_train_all_%s_%s%s.mat', interceptstr, nstatstr,flipstr);
        save(fullfile(baseDir, 'cluster_sess4', fname), 'T','Tcv','thetacv','compidx','M','MF')
        sml1_pcm('pcm_show_cluster_result', Tcv, M, compI, fname, isSave);
        
        varargout = {};    
    case 'pcm_fitGroupCrossval_untrain_allgroup' % use only trained sequences, but all subjects
        intercept = 0;
        useNaturalStat = 0;
        flipgroup = 0;
        vararginoptions(varargin, {'sess','intercept','useNaturalStat','flipgroup'});
        
        interceptstr = interceptstrs{intercept+1};
        nstatstr = naturalstatstrs{useNaturalStat+1};
        flipstr = flipstrs{flipgroup+1};
        
        % ============================== %
        % load model
        [S,G] = sml1_pcm('show_predictedRDM');
        Gc=G;
        % Make model component for both groups
        c=0; M=[];
        for m=[3,4,5,6];% firstFing, allFing, chunk, fingeTran
            c=c+1;
            if flipgroup
                M{c}.Gc = blockdiag(Gc(1:6,1:6,m),Gc(7:12,7:12,m));
            else
                M{c}.Gc = blockdiag(Gc(7:12,7:12,m),Gc(1:6,1:6,m));
            end
            M{c}.name = S.modelName{m};
        end
        M{end+1}.type = 'component';
        M{end}.Gc = eye(12);
        M{end}.name = 'sequence';
        
        % Normalise the model components to make them comparable size
        con = indicatorMatrix('allpairs', 1:12);
        for i=1:numel(M)
            M{i}.type = 'component';
            M{i}.numGparams = 1;
            M{i}.Gc=(M{i}.Gc+M{i}.Gc')/2;
            M{i}.Gc=M{i}.Gc./mean(diag(M{i}.Gc));            
        end;
        
        % Make the whole model family
        [MF,compI]=pcm_constructModelFamily(M);
        for i=1:numel(MF)
            MF{i}.theta0 = log(ones(MF{i}.numGparams,1));
        end
        
        % ====
        % run pcm
        % ============================== %
        % load prewhitened beta and design
        load(fullfile(baseDir, 'cluster_sess4', 'prewhBeta.mat'));                
        for c=unique(Data.cluster)';
            D=getrow(Data, Data.cluster==c);
            Z1 = indicatorMatrix('identity', repmat([1:6]',8,1));
            Z0 = zeros(size(Z1));
            for s=1:length(D.cluster)
                switch D.group(s)
                    case 1
                        Z{s} = [Z1, Z0];
                    case 2
                        Z{s} = [Z0, Z1];
                end
                partVec{s} = kron([1:8]',ones(6,1));
            end
            [T{c,1}, Tcv{c,1}, thetacv{c,1}, compidx] = ...
                pcm_clusterROI_component(D.beta_untrain, partVec, Z, MF, compI,...
                'fitAlgorithm',fitAlgorithm,'runEffect',runEffect);
        end
        % save&show result
        fname = sprintf('pcm_cluster_untrain_all_%s_%s%s.mat', interceptstr, nstatstr,flipstr);
        if isSave
            save(fullfile(baseDir, fname), 'T','Tcv','thetacv','compidx','M','MF')
        end        
        sml1_pcm('pcm_show_cluster_result', Tcv, M, compI,fname,isSave);
        
        varargout = {};
    
    case 'pcm_show_cluster_result' % plot component logBF for each cluster
        Tcv = varargin{1};
        M = varargin{2};
        compI = varargin{3};
        figfname = varargin{4};        
        isSave = varargin{5};
        
        if ischar(Tcv)&&isempty(M)&&isempty(compI)
            fname=fullfile(baseDir, 'cluster_sess4', sprintf('pcm_cluster_%s.mat', Tcv));
            load(fname);
            compI = compidx;
            [~,figfname,ext] = fileparts(fname);
        end
        
        % model name
        for m=1:numel(M)
            modelname{m} = M{m}.name;
        end
        modelname{end+1} = 'noiseceiling';
        
        % plot
        for g=1:size(Tcv,2)
            h{g}=myFigure([45,21],'name', 'PCM result');
            T = Tcv(:,g);
            for c=1:numel(T);
                like = T{c,1}.likelihood;
                
                like = bsxfun(@minus, like, like(:,1));
                modelLike = like(:,1:end-1);
                noiseCeiling = like(:,end);
                [posterior, logBF] = pcm_componentPosterior(modelLike, compI);
                
                subplot(3,4,c);
                %[x,y,e]=barplot([], [logBF, noiseCeiling]); title(sprintf('C=%d',c));
                [x,y,e]=barplot([], [logBF]);
                %myboxplot(kron([1:size(logBF,2)]', ones(size(logBF,1),1)), reshape(logBF, numel(logBF),1), 'plotall',1);
                drawline(0,'dir','horz');
                title(sprintf('C=%d',c));
                xlabel('Models'); ylabel('component logBF');
                if c==1; set(gca, 'xticklabel', modelname); end;
                set(gca,'view', [90,90]);
                set(gca,'ylim',[-20,30]);
            end
            
            % Save figure
            [~,figfname,ext] = fileparts(figfname);
            figName = fullfile(figDir,sprintf('%s_%d',figfname,g));
            mySaveFig(h{g}, [figName,'.eps'], isSave, '-dpsc2');
            mySaveFig(h{g}, figName, isSave, '-dpng','-r200');
        end        
        varargout = {h};
        
    case '*pcm_fitGroupCrossval_allseq' % use all sequences, but separate fit for two groups
        sess = 4;
        group=[1,2];
        vararginoptions(varargin, {'group'});
        
        % ============================== %
        % load model
        [S,G] = sml1_pcm('show_predictedRDM');
        Gc=G;
        % Make model component for both groups
        c=0; M=[];
        for m=[3,4,6];% firstFing, allFing, fingeTran
            c=c+1;
            M{c}.Gc = Gc(:,:,m);
            M{c}.name = S.modelName{m};
        end
        M{end+1}.type = 'component';
        M{end}.Gc = blockdiag(Gc(1:6,1:6,5), zeros(6));
        M{end}.name = 'chunk_train';
        M{end+1}.type = 'component';
        M{end}.Gc = blockdiag(zeros(6), Gc(7:12,7:12,5));
        M{end}.name = 'chunk_untrain';
        M{end+1}.type = 'component';
        M{end}.Gc = blockdiag(Gc(1:6,1:6,2), zeros(6));
        M{end}.name = 'seq_train';
        M{end+1}.type = 'component';
        M{end}.Gc = blockdiag(zeros(6), Gc(7:12,7:12,2));
        M{end}.name = 'seq_untrain';
        % add train/untrain difference as intercept
        M{end+1}.type='component';
        M{end}.Gc = blockdiag(ones(6),ones(6));
        M{end}.name = 'train-untrain';
        
        % Normalise the model components to make them comparable size
        con = indicatorMatrix('allpairs', 1:12);
        for i=1:numel(M)
            M{i}.type = 'component';
            M{i}.numGparams = 1;
            M{i}.Gc=(M{i}.Gc+M{i}.Gc')/2;
            M{i}.Gc=M{i}.Gc./mean(diag(M{i}.Gc));            
        end;
        
        % Make the whole model family
        [MF,compI]=pcm_constructModelFamily(M, 'alwaysInclude',8);
        % remove model containing only either of train or untrain
        idxc=compI(:,4)+compI(:,5)~=1;
        idxs=compI(:,6)+compI(:,7)~=1;
        idxcs=idxc&idxs;
        MF=MF(idxcs);
        compI=compI(idxcs,:);        
        for i=1:numel(MF)
            MF{i}.theta0 = log(ones(MF{i}.numGparams,1));
        end
        % ====
        % run pcm
        % ============================== %
        % load prewhitened beta and design
        load(fullfile(baseDir, 'cluster_sess4', 'prewhBeta.mat'));
        for g=group;
            Design.partVec = kron([1:8]',ones(12,1));
            switch g
                case 1
                    Design.condVec = repmat([1:12]',8,1);
                case 2
                    Design.condVec = repmat([7:12,1:6]',8,1);
            end
            
            for c=unique(Data.cluster)';
                D=getrow(Data, Data.cluster==c&Data.group==g);
                
                [T{c,g}, Tcv{c,g}, thetacv{c,g}, compidx] = ...
                    pcm_clusterROI_component(D.beta_all, Design.partVec, Design.condVec, MF, compI,...
                    'fitAlgorithm',fitAlgorithm,'runEffect',runEffect);
            end
        end
        
        % save&show result
        %save(fullfile(baseDir, 'cluster_sess4', 'pcm_cluster_train.mat'), 'T','Tcv','thetacv','compidx','M','MF')
        save(fullfile(baseDir, 'cluster_sess4', 'pcm_cluster_all.mat'), 'T','Tcv','thetacv','compidx','M','MF')
        for g=[1,2]
            sml1_pcm('pcm_show_cluster_result', Tcv(:,g), M, compI);
        end

        varargout = {};
    case '*pcm_fitGroupCrossval_allseq_allgroup' % use all sequences, but separate fit for two groups
        sess = 4;
        group=[1,2];
        vararginoptions(varargin, {'group'});
        
        % ============================== %
        % load model
        [S,G] = sml1_pcm('show_predictedRDM');
        Gc=G;
        % Make model component for both groups
        c=0; M=[];
        for m=[3,4,6];% firstFing, allFing, fingeTran
            c=c+1;
            M{c}.Gc = Gc(:,:,m);
            M{c}.name = S.modelName{m};
        end
        M{end+1}.type = 'component';
        M{end}.Gc = blockdiag(Gc(1:6,1:6,5), zeros(6));
        M{end}.name = 'chunk_train';
        M{end+1}.type = 'component';
        M{end}.Gc = blockdiag(zeros(6), Gc(7:12,7:12,5));
        M{end}.name = 'chunk_untrain';
        M{end+1}.type = 'component';
        M{end}.Gc = blockdiag(Gc(1:6,1:6,2), zeros(6));
        M{end}.name = 'seq_train';
        M{end+1}.type = 'component';
        M{end}.Gc = blockdiag(zeros(6), Gc(7:12,7:12,2));
        M{end}.name = 'seq_untrain';
        % add train/untrain difference as intercept
        M{end+1}.type='component';
        M{end}.Gc = blockdiag(ones(6),ones(6));
        M{end}.name = 'train-untrain';
        
        % Normalise the model components to make them comparable size
        con = indicatorMatrix('allpairs', 1:12);
        for i=1:numel(M)
            M{i}.type = 'component';
            M{i}.numGparams = 1;
            M{i}.Gc=(M{i}.Gc+M{i}.Gc')/2;
            M{i}.Gc=M{i}.Gc./mean(diag(M{i}.Gc));            
        end;
        
        % Make the whole model family
        [MF,compI]=pcm_constructModelFamily(M, 'alwaysInclude',8);
        % remove model containing only either of train or untrain
        idxc=compI(:,4)+compI(:,5)~=1;
        idxs=compI(:,6)+compI(:,7)~=1;
        idxcs=idxc&idxs;
        MF=MF(idxcs);
        compI=compI(idxcs,:);        
        for i=1:numel(MF)
            MF{i}.theta0 = log(ones(MF{i}.numGparams,1));
        end
        % ====
        % run pcm
        % ============================== %
        % load prewhitened beta and design
        load(fullfile(baseDir, 'cluster_sess4', 'prewhBeta.mat'));                
        for c=unique(Data.cluster)';
            D=getrow(Data, Data.cluster==c);
            for s=1:length(D.cluster)
                switch D.group(s)
                    case 1
                        Design.condVec = repmat([1:12]',8,1);
                    case 2
                        Design.condVec = repmat([7:12,1:6]',8,1);
                end
                Z{s} = indicatorMatrix('identity', Design.condVec);
                partVec{s} = kron([1:8]',ones(12,1));
            end
            [T{c,1}, Tcv{c,1}, thetacv{c,1}, compidx] = ...
                pcm_clusterROI_component(D.beta_all, partVec, Z, MF, compI,...
                'fitAlgorithm',fitAlgorithm,'runEffect',runEffect);
        end
        % save&show result
        %save(fullfile(baseDir, 'cluster_sess4', 'pcm_cluster_train.mat'), 'T','Tcv','thetacv','compidx','M','MF')
        save(fullfile(baseDir, 'cluster_sess4', 'pcm_cluster_all_all.mat'), 'T','Tcv','thetacv','compidx','M','MF')
        sml1_pcm('pcm_show_cluster_result', Tcv, M, compI);
        
        varargout = {};
    
end
end
function F = TwoFingTransition(Seq)
[Nseq,Npress] = size(Seq);
F=zeros(Nseq,25);
for s=1:Nseq
    for p=1:Npress-1;
        two(s,p) = 5*(Seq(s,p) - 1) + Seq(s,p+1);
        F(s, two(s,p)) = F(s, two(s,p)) + 1;
    end
    
end

end
function varargout = pcm_clusterROI_component(Yprewh, partitionVec, conditionVec, MF, CompIdx, varargin)
%% function out = sh1_pcm_clusterROI_component(Yprewh, partitionVec, conditionVec, MF, CompIdx, varargin)
%
% Inputs:
%	Yprewh: A 1xNsubj cell array of prewhitened beta estimates
%
%	partitionVec: A 1xNsubj cell array containing partition vectors.
%
%	conditionVec: A 1xNsubj cell array containing either condition vectors
%						 or design matrices (Z).
%
%	MF : Model family.
%
%   CompIdx: Component indices.
%
%
% Optional Inputs:
%	runEffect: The way how the run-effect is treated in fitting. Default is 'fixed'.
%
%	verbose: If set to 1, the program prints status in the command line. Default is 0.
%
%	fitAlgorithm: Algorithm for fitting. Default is 'NR'.
%
%	prior: Prior probability for calculating component posterior. Default is 0.5;
%
% a-yokoi (2018)

runEffect   = 'fixed'; % run effect
verbose     = 1;
fitAlgorithm= 'NR';
prior = 0.5;
calcCeiling = 1;
sn = [];
numIter = 1e3;
estEffVox=1;
pcm_vararginoptions(varargin,{'runEffect', 'verbose','fitAlgorithm', 'prior','calcCeiling', 'sn','numIter','estEffVox'});

%-----------------------------------------------------------------%
% Deal with missng subjects
if ~isempty(sn)&&islogical(sn);
    Yprewh=Yprewh(sn);
    conditionVec=conditionVec(sn);
    partitionVec=partitionVec(sn);
end

if calcCeiling
    MF{end+1}.type = 'freedirect';
    MF{end}.numGparams = 0;
    MF{end}.theta0 = [];
    MF{end}.name = 'noise_ceiling';
end

% Run group-fit
[T, theta] = pcm_fitModelGroup(Yprewh, MF, partitionVec, conditionVec,...
    'runEffect', runEffect, ...
    'fitAlgorithm', 'NR', ...
    'verbose', verbose, 'MaxIteration', numIter);

% Do the crossvalidated group-fit
[Tcv, theta_cv] = pcm_fitModelGroupCrossval(Yprewh, MF, partitionVec, conditionVec,...
    'runEffect', runEffect, ...
    'fitAlgorithm', fitAlgorithm, ...
    'groupFit', theta, 'MaxIteration', numIter);

varargout = {T,Tcv, theta_cv, CompIdx};
end