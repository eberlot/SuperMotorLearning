function varargout=sml_analyze (what, varargin)

prefix = 'sml1_';
%baseDir = '/Users/eberlot/Documents/Data/SuperMotorLearning';
baseDir = '/Volumes/MotorControl/data/SuperMotorLearning';
behDir  = fullfile(baseDir,'behavioral_data');
anaDir  = fullfile(behDir,'analyze');
memoryDir = fullfile(behDir,'memory');
subj_name  = {'s01','s02','s03','s04','s05','s06','s07','s08','s09','s10','s11','s12','s13','s14','s15','s16','s17','s18',...
                's19','s20','s21','s22','s23','s24','s25','s26','s27','s28','s29','s30','s31'};
%cd(behDir)

% -------------------------- For plotting ---------------------------------
stySeq=style.custom({'red','blue'},'markersize',12);
c1=[39 38 124]/255;
c2=[140 140 185]/255;
c3=[249 191 193]/255;
%styOtherHand= style.custom({'red','blue','green'},'markersize',10,'errorbars','shade');

gray=[80 80 80]/255;
lightgray=[160 160 160]/255;
black=[0 0 0]/255;
blue=[0 0 1];
red=[1 0 0];
blue=[49,130,189]/255;
lightblue=[158,202,225]/255;
red=[222,45,38]/255;
lightred=[252,146,114]/255;
styRSbeh = style.custom({red,lightred,blue,lightblue},'markersize',12);
stySeq2 = style.custom({red,blue,lightred,lightblue},'markersize',12);
%styOtherHand= style.custom({c1,c2,c3},'markersize',10,'errorbars','shade');
%styOtherHand= style.custom({gray,lightgray,black},'markersize',10,'errorbars','shade','linestyle',{'-','--','-.'});
styOtherHand= style.custom({blue,red,black},'markersize',10,'errorbars','shade','linestyle',{'-','--','-.'});
styTest = style.custom({c1,c2});
styTrio = style.custom({red,c1,black,c3,c2,gray});
switch what
    case 'save_subj'  % create .mat structure for each subject from .dat file
        sn=[4:9,11:31];
        vararginoptions(varargin,{'sn'});
        
        for i=sn
            datafilename = fullfile(behDir,['data/sml1_' subj_name{i} '.dat']);
            outfilename  = fullfile(behDir,['analyze/sml1_' subj_name{i} '.mat']);
            D=dload(datafilename);
            save(outfilename,'-struct','D');
            fprintf('%d.',i);
        end;
    case 'all_subj'   % go through trial structure
        sn=1;
        vararginoptions(varargin,{'sn'});
        
        for i=sn
            sml_subj(subj_name{i});
            fprintf('%d.',i);
        end;
    case 'addDay'
       vararginoptions(varargin,{'sn'});
       
       for s = sn   
            D=load(fullfile(anaDir,['sml1_',subj_name{s},'.mat']));
            dayIndx=dlmread(fullfile(anaDir,[subj_name{s},'_BNday.txt']));
            
            Tt=[];
            for d = 1:(length(dayIndx)-1)
                idx=[(dayIndx(d)+1):dayIndx(d+1)];
                A = getrow(D,ismember(D.BN,idx));
                T.day = ones(length(A.BN),1)*d;
                Tt=addstruct(Tt,T);
            end
            D.day = Tt.day;
            
            save(fullfile(anaDir,['sml1_',subj_name{s},'.mat']),'-struct','D');
        end
    
    case 'make_alldat'      % load all subjects dat files, concatenate them - one structure for all subjects 
        sn=[4:9,11:31];
        T=[]; 
        for i=sn
            D=load(fullfile(behDir,'analyze',[prefix subj_name{i} '.mat'])); 
            D.SN=ones(size(D.BN))*i; 
            T=addstruct(T,D); 
        end; 
        save(fullfile(anaDir,'alldata.mat'),'-struct','T'); 
    
    case 'PLOT_learning_group'            % plot learning curve (MT) for seqLearn blocks
        vararginoptions(varargin,{'sn'});
        
        D=load(fullfile(behDir,'analyze','alldata.mat'));
     
        sty= style.custom({'magenta'},'markersize',12,'errorbars','shade');
        figure
        %plt.line([D.day>6 D.day],D.MT,'subset',D.blockType==6,'style',sty);
        plt.line(D.day,D.MT,'subset',D.blockType==6,'style',sty);
        ylabel('Movement time (msec)'); xlabel('Days');
        
        figure
        plt.line(D.BN,D.MT,'style',styOtherHand,'subset',D.blockType==6);
        ylabel('Movement time (msec)'); xlabel('Block number');
        
        figure
        plt.line(D.day,D.MT,'style',styOtherHand,'subset',D.blockType==6,'split',D.FoSEx,'leg',{'1st','2nd'});
        ylabel('Movement time (msec)'); xlabel('Days of training');
    case 'PLOT_learning_sn'
        sn=1; 
        vararginoptions(varargin,{'sn'});
        figure
        for i=1:length(subj_name)
            load(fullfile(behDir,'analyze',['sml1_' subj_name{i} '.mat']));
            lineplot(D.BN,D.MT,'subset',D.blockType==6);
            xlabel('Block number'); ylabel('Movement time');
            hold on
        end
    case 'FIT_learningCurve'
        sn=[4:9,11:31];
        var='BN'; % use blocknumber or day to fit
        fig=1;
        vararginoptions(varargin,{'sn','var','fig'});
        D=load(fullfile(behDir,'analyze','alldata.mat'));
        % function to fit - exponential
        F = @(x,xdata)x(1)*exp(-x(2)*xdata)+x(3);
        SS=[];
        for s=sn
            T = getrow(D,D.blockType==6&D.SN==s);
            T.dayIndx = T.day;
            dayIndx=unique(T.day);
            for i=1:length(dayIndx)
                T.dayIndx(T.dayIndx==dayIndx(i))=i;
            end
            % get the median
            M=tapply(T,{var},{'MT','median'},'subset',T.isError~=1);
            
            if strcmp(var,'BN')
                M.BN=(1:size(M.MT))';
            end
            ydata = M.MT;
            xdata = M.(var);
            x0 = [mean(M.MT(M.(var)==1)) 0.5 mean(M.MT(M.(var)==i))];
            xunc = lsqcurvefit(F,x0,xdata,ydata);
            % now calculate residuals
            Y_pred = F(xunc,M.(var));
            res = M.MT - Y_pred;
            EE=[];
            for i = 1:1000  % start the bootstrapping loop
                % create data (fit + residuals - sample with replacement)
                Y_new = Y_pred + sample_wr(res,1,length(res));
                xunc_boot = lsqcurvefit(F,x0,xdata,Y_new);
                E.x0 = xunc_boot(1);
                E.x1 = xunc_boot(2);
                E.x2 = xunc_boot(3);
                E.init = i;
                EE=addstruct(EE,E);
            end
            S.x0 = xunc(1);
            S.x1 = xunc(2);
            S.x2 = xunc(3);
            S.x0_mean = mean(EE.x0);
            S.x1_mean = mean(EE.x1);
            S.x2_mean = mean(EE.x2);
            S.x0_std  = std(EE.x0);
            S.x1_std  = std(EE.x1);
            S.x2_std  = std(EE.x2);
            x0_SEM    = S.x0_std/sqrt(length(EE.x0));
            x1_SEM    = S.x1_std/sqrt(length(EE.x1));
            x2_SEM    = S.x2_std/sqrt(length(EE.x2));
            ts = tinv([0.025  0.975],length(EE.x0)-1);      % T-Score
            S.x0_CI_upp = S.x0_mean + ts(2)*x0_SEM;  
            S.x0_CI_low = S.x0_mean + ts(1)*x0_SEM;  
            S.x1_CI_upp = S.x1_mean + ts(2)*x1_SEM;  
            S.x1_CI_low = S.x1_mean + ts(1)*x1_SEM; 
            S.x2_CI_upp = S.x2_mean + ts(2)*x2_SEM;  
            S.x2_CI_low = S.x2_mean + ts(1)*x2_SEM; 
            S.sn      = s;
            SS=addstruct(SS,S);
            % refit the curve
            % mean: x0:1167 x1:0.098; std: x0: 149.50, x1: 0.179
            if fig
                figure
                plt.line(M.(var),M.MT);
                hold on;
                plot(xdata,F(xunc,M.(var)),'-k','LineWidth',2);
                hold on;
                plot(xdata,F([mean(EE.x0) mean(EE.x1) mean(EE.x2)],M.(var)),'-r','LineWidth',2);
                xlabel(var);
                ylabel('Movement time');
                title(subj_name{s});
            end
        end
        save(fullfile(anaDir,'exponentialFits'),'-struct','SS');
        % to check the overlap (note: s04 VERY different)
        figure
        plt.hist(SS.x2_mean,'subset',SS.sn~=4);
        hold on
        for i=2:27
            drawline(SS.x2_CI_upp(i),'dir','vert');
            drawline(SS.x2_CI_low(i),'dir','vert','color',[1 0 0]);
        end
    case 'HIST_errors_points_group' 
        sn = [5:9,11:31];
        vararginoptions(varargin,{'sn'});
        D=load(fullfile(behDir,'analyze','alldata.mat'));
     
        figure
        subplot(1,2,1) 
        histplot(D.points(D.blockType==6),'subset',ismember(D.SN,sn));
        xlabel('Points'); ylabel('Proportion');
        subplot(1,2,2)
        histplot(D.isError(D.blockType==6),'subset',ismember(D.SN,sn));
        xlabel('Error'); ylabel('Proportion');       
    case 'HIST_errors_points_sn'
        sn=1; 
        vararginoptions(varargin,{'sn'});
        
        load(fullfile(behDir,'analyze',['sml1_' subj_name{sn} '.mat']));

        figure
        subplot(1,2,1) 
        histogram(D.points(D.blockType==6),'Normalization','probability');
        xlabel('Points'); ylabel('Proportion');
        subplot(1,2,2)
        histogram(D.isError(D.blockType==6),'Normalization','probability');
        xlabel('Error'); ylabel('Proportion');

    case 'CL_MT'            % plot learning curve (MT) for chunkLearn blocks
        figure
        for i=1:length(subj_name)
            load(fullfile(behDir,'analyze',['sml1_' subj_name{i} '.mat']));
            subplot(1,length(subj_name),i);
            lineplot(D.BN,D.MT,'subset',D.blockType==5);
            xlabel('Block number'); ylabel('Movement time');
            hold on
        end
    case 'CL_error_points'
        figure
        for i=1:length(subj_name)
            load(fullfile(behDir,'analyze',['sml1_' subj_name{i} '.mat']));
            subplot(2,length(subj_name),i);
            histogram(D.points(D.blockType==6),'Normalization','probability')
            xlabel('Points'); ylabel('Proportion');
            subplot(2,length(subj_name),i+3)
            histogram(D.isError(D.blockType==6),'Normalization','probability')
            xlabel('Error'); ylabel('Proportion');
            hold on
        end
        
    case 'PLOT_FoSEx_sn'
        vararginoptions(varargin,{'sn'});
        
        for i=1:length(sn)
            D=load(fullfile(behDir,'analyze',['sml1_' subj_name{sn(i)} '.mat']));
            R = getrow(D,D.blockType==11);
            figure(1)
            subplot(1,length(sn),i)
            myboxplot(R.day,R.MT,'split',R.FoSEx,'style_tukey','leg',{'first','second'},'leglocation','northeast');
            xlabel('Day'); ylabel('Movement time');
            hold on;
        end
        
    case 'PLOT_testFoSEx_MT'
        D=load(fullfile(anaDir,'testsRH.mat'));
        
%         sty = style.custom({'red','blue'},'markersize',12);
%         figure
%         subplot(1,3,1)
%         plt.line(D.testDay,D.MT,'split',D.FoSEx,'subset',D.seqType==1,'leg',{'1st','2nd'},'style',sty);
%         title('Trained sequences'); xlabel('Test day'); ylabel('Movement time');
%         
%         subplot(1,3,2)
%         plt.line(D.testDay,D.MT,'split',D.FoSEx,'subset',D.seqType==4,'leg',{'1st','2nd'},'style',sty);
%         title('Random sequences'); ylabel('');
%         
%         subplot(1,3,3)
%         plt.line(D.testDay,D.MT,'split',D.FoSEx,'subset',D.seqType==5,'leg',{'1st','2nd'},'style',sty);
%         title('Chunked sequences'); ylabel('');
%         
%         plt.match('y');
%         
        D = getrow(D,D.seqType~=4);
        D.seqType(D.seqType==4)=2;
         % reduced structure
        [fa ra ca] = pivottable([D.testDay D.SN],[D.FoSEx D.seqType],D.MT,'mean');
        N.MT = [fa(:,1);fa(:,2);fa(:,3);fa(:,4)];
        N.FoSEx = [ones(size(fa(:,1)));ones(size(fa(:,1)));ones(size(fa(:,1)))*2;ones(size(fa(:,1)))*2];
        N.seqType = [ones(size(fa(:,1)));ones(size(fa(:,1)))*2;ones(size(fa(:,1)));ones(size(fa(:,1)))*2];
        N.sn = [ra(:,2);ra(:,2);ra(:,2);ra(:,2)];
        N.testDay = [ra(:,1);ra(:,1);ra(:,1);ra(:,1)];
        
        figure
        plt.dot(N.testDay,N.MT,'split',[N.seqType N.FoSEx],'leg',{'1st train','2nd train','1st untrain','2nd untrain'},'leglocation','northeast','style',styRSbeh);
        xlabel('Test day');
        ylabel('Movement time');
        %
        figure
        plt.box(N.testDay,N.MT,'split',[N.seqType N.FoSEx],'plotall',0,'leg',{'1st train','2nd train','1st untrain','2nd untrain'},'leglocation','northeast','style',styRSbeh);
        xlabel('Test day');
        ylabel('Movement time');
        %      
        % stats: ANOVA (seqType x repetition)
        for t=unique(N.testDay)'
            fprintf('Anova day %d\n',t);
            anovaMixed(N.MT,N.sn,'within',[N.seqType N.FoSEx],{'seqType','FoSEx'},'subset',N.testDay==t);
        end
        keyboard;
        % calculate the percentage per person
        P = getrow(N,N.FoSEx==1);
        P.MT = N.MT(N.FoSEx==2)./N.MT(N.FoSEx==1);
        figure
        plt.dot(P.testDay, P.MT,'split',P.seqType,'leg',{'trained','untrained'},'style',stySeq);
        drawline(1,'dir','horz');
        xlabel('Test day');
        ylabel('Movement time (ratio: 2nd/1st)');
        % stats on ratios
        keyboard;
        for t=unique(P.testDay)'
            fprintf('Day %d\n',t);
            fprintf('One sample trained\n');
            ttestDirect(P.MT-1,[P.sn],2,'onesample','subset',P.testDay==t & P.seqType==1);           
            fprintf('One sample random\n');
            ttestDirect(P.MT-1,[P.sn],2,'onesample','subset',P.testDay==t & P.seqType==2);
            fprintf('Paired\n');
            ttestDirect(P.MT,[P.seqType P.sn],2,'paired','subset',P.testDay==t);
        end
        
    case 'PLOT_testFoSEx_ER'
        D=load(fullfile(anaDir,'testsRH.mat'));
        
        sty= style.custom({'red','blue'},'markersize',12);
        figure
        subplot(1,3,1)
        plt.line(D.testDay,D.isError,'split',D.FoSEx,'subset',D.seqType==1,'leg',{'1st','2nd'},'style',sty);
        title('Trained sequences'); xlabel('Test day'); ylabel('Movement time');
        
        subplot(1,3,2)
        plt.line(D.testDay,D.isError,'split',D.FoSEx,'subset',D.seqType==4,'leg',{'1st','2nd'},'style',sty);
        title('Random sequences'); ylabel('');
        
        subplot(1,3,3)
        plt.line(D.testDay,D.isError,'split',D.FoSEx,'subset',D.seqType==5,'leg',{'1st','2nd'},'style',sty);
        title('Chunked sequences'); ylabel('');
        
        plt.match('y');
    case 'FORM_testRH_psc'
        D=load(fullfile(anaDir,'testsRH.mat'));
        
        D=getrow(D,D.SN~=2 & D.seqType~=5); % remove the first test day - only random
        D.seqType(D.seqType==4)=2;
        % form a new structure - percentage of decrease
        T=tapply(D,{'testDay','SN','seqType','FoSEx'},{'MT','median'},'subset',D.isError~=1);
       
        
        % create a new structure
        N1 = getrow(T,T.FoSEx==1);
        N2 = getrow(T,T.FoSEx==2);
        N = N1;
        N.MT = N2.MT./N1.MT;
        T.MT_psc=ones(size(T.MT));
        
         for s=unique(T.SN)' % subject-specific normalisation
             MT_init_sn=mean(T.MT(T.testDay==1&T.SN==s));
             for t=2:4  % test days
                T.MT_psc(T.testDay==t & T.SN==s)=T.MT(T.testDay==t & T.SN==s)/MT_init_sn;             
             end
         end
         
         figure
         plt.line(T.testDay,T.MT_psc,'split',T.seqIdx,'style',styOtherHand,'leg',{'Intrinsic','Extrinsic','Random'});
         xlabel('Test days'); ylabel('Relative movement time');
         
         save(fullfile(anaDir,'tests_OtherHand_psc.mat'),'-struct','T');

    case 'FORM_testRH_group'
        sn=[4,5,7:9,11:31];
        vararginoptions(varargin,{'sn'});
        D=load(fullfile(anaDir,'alldata_wChunks.mat'));
        
        TT=[];  
        for i=1:length(sn)
            T=getrow(D,D.SN==sn(i)&D.blockType==11);
            
            testDay=unique(T.day);
            for t = 1:length(testDay) % test days
                T.testDay(T.day==testDay(t),:)=t;
            end
            TT=addstruct(TT,T);        
        end
        
         save(fullfile(anaDir,'testsRH.mat'),'-struct','TT');
    case 'PLOT_testRH'
        D=load(fullfile(anaDir,'testsRH.mat'));
        
        sty= style.custom({'magenta','green'},'markersize',12);
        figure
        plt.box(D.testDay,D.MT,'split',D.seqType,'subset',D.seqType==1 | D.seqType==4,'leg',{'trained','random'},'style',sty,'plotall',0);
        figure
        plt.box(D.testDay,D.MT,'split',D.seqType,'leg',{'trained','random','chunked'},'plotall',0);
        %title('Trained sequences'); 
        xlabel('Test day'); ylabel('Movement time');
    case 'STATS_testRH'
        D=load(fullfile(anaDir,'testsRH.mat'));
        
        var=D.MT; % or D.isError
        % only trained vs. random, not the chunked ones
        anovaMixed(var,D.SN,'within',[D.testDay D.seqType],{'testDay','seqType'},'subset',D.seqType~=4);
        
        % follow-up t-tests per testing day

        % test day 1 
        ttestDirect(var,[D.seqType D.SN],2,'paired','subset',D.testDay==1 & D.seqType~=4);         % trained vs. random
        % test day 2 
        ttestDirect(var,[D.seqType D.SN],2,'paired','subset',D.testDay==2 & D.seqType~=4);
        % test day 3 
        ttestDirect(var,[D.seqType D.SN],2,'paired','subset',D.testDay==3 & D.seqType~=4);    
        % test day 4 
        ttestDirect(var,[D.seqType D.SN],2,'paired','subset',D.testDay==4 & D.seqType~=4);
        
    case 'PLOT_FoSEx'
        D=load(fullfile(anaDir,'testsRH.mat'));
        
        figure
        subplot(311)
        plt.line(D.testDay,D.MT,'split',D.FoSEx,'subset',D.seqType==1,'style',styTest);
        ylabel('MT');
        subplot(312)
        plt.line(D.testDay,D.MT,'split',D.FoSEx,'subset',D.seqType==4,'style',styTest);
        ylabel('MT');
        subplot(313)
        plt.line(D.testDay,D.MT,'split',D.FoSEx,'subset',D.seqType==5,'style',styTest);
        ylabel('MT');
        xlabel('Test Day');
    case 'STATS_FoSEx'
        D=load(fullfile(anaDir,'testsRH.mat'));
        
        var=D.MT; % or D.isError
        % only trained vs. random, not the chunked ones
        anovaMixed(var,D.SN,'within',[D.testDay D.FoSEx],{'testDay','FoSEx'},'subset',D.seqType==1);
        anovaMixed(var,D.SN,'within',[D.testDay D.FoSEx],{'testDay','FoSEx'},'subset',D.seqType==4);
        anovaMixed(var,D.SN,'within',[D.testDay D.FoSEx],{'testDay','FoSEx'},'subset',D.seqType==5);

        % follow-up t-tests per testing day

        % test day 1 
        ttestDirect(var,[D.FoSEx D.SN],2,'paired','subset',D.testDay==1 & D.seqType==4);         % trained vs. random
        % test day 2 
        ttestDirect(var,[D.FoSEx D.SN],2,'paired','subset',D.testDay==2 & D.seqType==4);
        % test day 3 
        ttestDirect(var,[D.FoSEx D.SN],2,'paired','subset',D.testDay==3 & D.seqType==4);    
        % test day 4 
        ttestDirect(var,[D.FoSEx D.SN],2,'paired','subset',D.testDay==4 & D.seqType==4);

    case 'FORM_OtherHand'
        sn=[4:9,11:31];
        vararginoptions(varargin,{'sn'});
        D=load(fullfile(anaDir,'alldata.mat'));
        
        TT=[];     
        for i=1:length(sn)
            T=getrow(D,D.SN==sn(i) & (D.blockType==10 | D.blockType==12));      
            % determine days
            BN=unique(T.BN);
            BN_indx=[0 find(diff(BN)>1)' length(BN)];
            T.testDay=zeros(size(T.day));  
            % identify days based on separate in numbers of consequent BN
            testDay=unique(T.day);
            for t = 1:length(BN_indx)-1 % test days
                indx=[BN(BN_indx(t)+1):BN(BN_indx(t+1))];
                T.testDay(ismember(T.BN,indx'))=t;
            end
            TT=addstruct(TT,T);
        end  
        TT.OH_seqIdx=TT.seqType;
        TT.OH_seqIdx(TT.seqType==6)=1;  % trained intrinsic
        TT.OH_seqIdx(TT.seqType==7)=2;  % trained extrinsic
        TT.OH_seqIdx(TT.seqType==8 | TT.seqType==9)=3; %random
        
        save(fullfile(anaDir,'tests_OtherHand.mat'),'-struct','TT');
    case 'PLOT_OtherHand'
        D=load(fullfile(anaDir,'tests_OtherHand.mat'));
        
        figure
        plt.line(D.testDay,D.isError,'split',D.OH_seqIdx,'leg',{'train intr','train extr','rand'},'style',styOtherHand);
        title('Other hand transfer'); xlabel('Test day'); ylabel('Error rate');
        
        figure
        plt.line(D.testDay,D.MT,'split',D.OH_seqIdx,'leg',{'train intr','train extr','rand'},'style',styOtherHand);
        title('Other hand transfer'); xlabel('Test day'); ylabel('Movement time');
        
        figure
        plt.line(D.testDay,D.MT,'split',D.OH_seqIdx,'leg',{'train intr','train extr','rand'},'style',styOtherHand,'subset',D.testDay>1 & D.isError~=1);
        title('Other hand transfer'); xlabel('Test day'); ylabel('Movement time');
    case 'STATS_OtherHand'
        D=load(fullfile(anaDir,'tests_OtherHand.mat'));
        
        var=D.MT; % or D.isError
        anovaMixed(var,D.SN,'within',[D.testDay D.OH_seqIdx],{'testDay','seqType'},'subset',D.testDay>1& D.isError~=1);
        % extrinsic vs. intrinsic
        anovaMixed(var,D.SN,'within',[D.testDay D.OH_seqIdx],{'testDay','seqType'},'subset',D.testDay>1 & D.OH_seqIdx<3 & D.isError~=1);
        % intrinsic vs. random
        anovaMixed(var,D.SN,'within',[D.testDay D.OH_seqIdx],{'testDay','seqType'},'subset',D.testDay>1 & D.OH_seqIdx~=2 & D.isError~=1);
        % extrinsic vs. random
        anovaMixed(var,D.SN,'within',[D.testDay D.OH_seqIdx],{'testDay','seqType'},'subset',D.testDay>1 & D.OH_seqIdx~=1 & D.isError~=1);
        % extrinsic vs. intrinsic test days 2 and 4    
        anovaMixed(var,D.SN,'within',[D.testDay D.OH_seqIdx],{'testDay','seqType'},'subset',D.testDay>1 & D.testDay~=3 & D.OH_seqIdx<3 & D.isError~=1);

        % follow-up ANOVAs & t-tests per testing day
        keyboard;
        % test day 1 pre-test -> only random
        % test day 2 (before learning)
        anovaMixed(var,D.SN,'within',D.OH_seqIdx,{'seqType'},'subset',D.testDay==2&D.isError~=1);
        
        % test day 3 (after 1 week)
        anovaMixed(D.MT,D.SN,'within',D.OH_seqIdx,{'seqType'},'subset',D.testDay==3);
        ttestDirect(D.MT,[D.OH_seqIdx D.SN],2,'paired','subset',D.testDay==3 & D.OH_seqIdx~=2&D.isError~=1);         % intr vs. random
        ttestDirect(D.MT,[D.OH_seqIdx D.SN],2,'paired','subset',D.testDay==3 & D.OH_seqIdx~=1&D.isError~=1);         % extr vs. random
        ttestDirect(D.MT,[D.OH_seqIdx D.SN],2,'paired','subset',D.testDay==3 & D.OH_seqIdx~=3&D.isError~=1);         % intr vs. extr
        
        % test day 4 (after 3 weeks)
        anovaMixed(var,D.SN,'within',D.OH_seqIdx,{'seqType'},'subset',D.testDay==4);
        ttestDirect(var,[D.OH_seqIdx D.SN],2,'paired','subset',D.testDay==4 & D.OH_seqIdx~=2&D.isError~=1);         % intr vs. random
        ttestDirect(var,[D.OH_seqIdx D.SN],2,'paired','subset',D.testDay==4 & D.OH_seqIdx~=1&D.isError~=1);         % extr vs. random
        ttestDirect(var,[D.OH_seqIdx D.SN],2,'paired','subset',D.testDay==4 & D.OH_seqIdx~=3&D.isError~=1);         % intr vs. extr
        
        % test day 5 (after 5 weeks)
        anovaMixed(var,D.SN,'within',D.OH_seqIdx,{'seqType'},'subset',D.testDay==5);
        ttestDirect(var,[D.OH_seqIdx D.SN],2,'paired','subset',D.testDay==5 & D.OH_seqIdx~=2&D.isError~=1);         % intr vs. random
        ttestDirect(var,[D.OH_seqIdx D.SN],2,'paired','subset',D.testDay==5 & D.OH_seqIdx~=1&D.isError~=1);         % extr vs. random
        ttestDirect(var,[D.OH_seqIdx D.SN],2,'paired','subset',D.testDay==5 & D.OH_seqIdx~=3&D.isError~=1);         % intr vs. extr
    
    case 'PLOT_OtherHand_interIndivid'
        D=load(fullfile(anaDir,'tests_OtherHand.mat'));
        
        figure
        for i=2:5 % each test-day
            [MT,subj]   = pivottable(D.SN,D.OH_seqIdx,D.MT,'nanmean','subset',D.testDay==i);
            subplot(1,4,i-1)
            % first remove performance on random sequences
            MT(:,1) = MT(:,3) - MT(:,1);
            MT(:,2) = MT(:,3) - MT(:,2);
            plt.scatter(MT(:,1),MT(:,2),'label',(1:27));
            hold on;
            ax = xlim;
            ay = ylim;
            xlim([min([ax(1) ay(1)]), max([ax(2) ay(2)])]);
            ylim([min([ax(1) ay(1)]), max([ax(2) ay(2)])]);
            xlabel('Intrinsic');
            ylabel('Extrinsic');
            title(sprintf('Test day %d',i-1));
            %refline(1,0);
            hold on;
            %axis equal;
            drawline(0,'dir','horz');
            drawline(0,'dir','vert');
            ind = MT(:,1)./MT(:,2);
            cIE=corr(ind,MT(:,3))
        end
    case 'PLOT_OtherHand_diff'
        D=load(fullfile(anaDir,'tests_OtherHand.mat'));
        TT=[];
        figure
        for i=2:5 % each test-day
            [MT,subj]   = pivottable(D.SN,D.OH_seqIdx,D.MT,'nanmean','subset',D.testDay==i);
            subplot(1,4,i-1)
            % first remove performance on random sequences
            %MT(:,1) = MT(:,3) - MT(:,1);
            %MT(:,2) = MT(:,3) - MT(:,2);
            %MT(:,1) = -(MT(:,1) - MT(:,3));
            %MT(:,1) = -(MT(:,2) - MT(:,3));
            diffTime = MT(:,1)-MT(:,2);
            plt.scatter([1:27]',diffTime,'label',(1:27));
           % plot([1:27],diffTime,'o');
            hold on;
            ylabel('Extrinsic - intrinsic');
            title(sprintf('Test day %d',i-1));
            %refline(1,0);
            hold on;
            %axis equal;
            drawline(0,'dir','horz');
            drawline(0,'dir','vert');
            ind = MT(:,1)./MT(:,2);
            cIE=corr(ind,MT(:,3))
            T.diffTime = diffTime;
            T.SN       = [1:27]';
            T.testDay  = ones(size(diffTime))*i;
            TT=addstruct(TT,T);
        end

        keyboard;
        T2=getrow(TT,TT.SN~=25); % remove outlier
        figure
        subplot(131)
        plt.scatter(T2.diffTime(T2.testDay==2),T2.diffTime(T2.testDay==3),'label',([1:24,26,27]));
        title(sprintf('corr %d',corr(T2.diffTime(T2.testDay==2),T2.diffTime(T2.testDay==3))));
        subplot(132)
        plt.scatter(T2.diffTime(T2.testDay==3),T2.diffTime(T2.testDay==4),'label',([1:24,26,27]));
        title(sprintf('corr %d',corr(T2.diffTime(T2.testDay==3),T2.diffTime(T2.testDay==4))));
        subplot(133)
        plt.scatter(T2.diffTime(T2.testDay==4),T2.diffTime(T2.testDay==5),'label',([1:24,26,27]));
        title(sprintf('corr %d',corr(T2.diffTime(T2.testDay==4),T2.diffTime(T2.testDay==5))));
        keyboard;
        
        
    case 'FORM_OtherHand_psc'
        D=load(fullfile(anaDir,'tests_OtherHand.mat'));
        
        D=getrow(D,D.testDay>1); % remove the first test day - only random
        D.testDay=D.testDay-1; % relabel test days - start with 1
        D.seqIdx=D.seqType;
        D.seqIdx(D.seqType==6)=1;  % trained intrinsic
        D.seqIdx(D.seqType==7)=2;  % trained extrinsic
        D.seqIdx(D.seqType==8 | D.seqType==9)=3; %random
        % form a new structure - percentage of decrease
        T=tapply(D,{'testDay','SN','seqIdx'},{'MT','mean'},'subset',D.isError~=1);
        
        % percentage of MT
        T.MT_psc=ones(size(T.MT));
        
         for s=unique(T.SN)' % subject-specific normalisation
             MT_init_sn=mean(T.MT(T.testDay==1&T.SN==s));
             for t=2:4  % test days
                T.MT_psc(T.testDay==t & T.SN==s)=T.MT(T.testDay==t & T.SN==s)/MT_init_sn;             
             end
         end
         
         figure
         plt.line(T.testDay,T.MT_psc,'split',T.seqIdx,'style',styOtherHand,'leg',{'Intrinsic','Extrinsic','Random'});
         xlabel('Test days'); ylabel('Relative movement time');
         
         save(fullfile(anaDir,'tests_OtherHand_psc.mat'),'-struct','T');
    case 'FORM_OtherHand_subtract'
        D=load(fullfile(anaDir,'tests_OtherHand.mat'));
        
        D=getrow(D,D.testDay>1); % remove the first test day - only random
        %D=getrow(D,D.testDay>1 & D.SN~=2); % remove the first test day - only random
        D.testDay=D.testDay-1; % relabel test days - start with 1
        D.seqIdx=D.seqType;
        D.seqIdx(D.seqType==6)=1;  % trained intrinsic
        D.seqIdx(D.seqType==7)=2;  % trained extrinsic
        D.seqIdx(D.seqType==8 | D.seqType==9)=3; %random
        % form a new structure - percentage of decrease
        T=tapply(D,{'testDay','SN','seqIdx'},{'MT','mean'},'subset',D.isError~=1);
        %T=tapply(D,{'testDay','SN','seqIdx'},{'MT','mean'}); 
       
        % percentage of MT
        T.MT_subtract=zeros(size(T.MT));
        
        for s=unique(T.SN)' % subject-specific normalisation
            for t=2:4  % test days
                for st=1:3 % seqIndex
                    MT_before_sn=mean(T.MT(T.testDay==t-1&T.SN==s&T.seqIdx==st));
                    T.MT_subtract(T.testDay==t & T.SN==s & T.seqIdx==st)=MT_before_sn-T.MT(T.testDay==t & T.SN==s & T.seqIdx==st);
                end
            end
        end
         
        figure
        plt.bar(T.testDay,T.MT_subtract,'split',T.seqIdx,'subset',T.testDay>1,'style',styOtherHand,'leg',{'Intrinsic','Extrinsic','Random'});
        drawline(0,'dir','horz');
        xlabel('Test days'); ylabel('Between subsequent session learning (msec)'); 
        
        
        anovaMixed(T.MT_subtract,T.SN,'within',[T.seqIdx T.testDay],{'seqType','day'},'subset',T.testDay>1);
        % day 2
        ttestDirect(T.MT_subtract,[T.seqIdx T.SN],2,'paired','subset',T.testDay==2 & T.seqIdx~=2);         % intr vs. random
        ttestDirect(T.MT_subtract,[T.seqIdx T.SN],2,'paired','subset',T.testDay==2 & T.seqIdx~=1);         % extr vs. random
        ttestDirect(T.MT_subtract,[T.seqIdx T.SN],2,'paired','subset',T.testDay==2 & T.seqIdx~=3);         % intr vs. extr
        % day 3
        ttestDirect(T.MT_subtract,[T.seqIdx T.SN],2,'paired','subset',T.testDay==3 & T.seqIdx~=2);         % intr vs. random
        ttestDirect(T.MT_subtract,[T.seqIdx T.SN],2,'paired','subset',T.testDay==3 & T.seqIdx~=1);         % extr vs. random
        ttestDirect(T.MT_subtract,[T.seqIdx T.SN],2,'paired','subset',T.testDay==3 & T.seqIdx~=3);         % intr vs. extr
        % day 4
        ttestDirect(T.MT_subtract,[T.seqIdx T.SN],2,'paired','subset',T.testDay==4 & T.seqIdx~=2);         % intr vs. random
        ttestDirect(T.MT_subtract,[T.seqIdx T.SN],2,'paired','subset',T.testDay==4 & T.seqIdx~=1);         % extr vs. random
        ttestDirect(T.MT_subtract,[T.seqIdx T.SN],2,'paired','subset',T.testDay==4 & T.seqIdx~=3);         % intr vs. extr
         
    case 'SAVE_memoryTests'
        vararginoptions(varargin,{'sn','sessN'});
       
        MM=[];
        for s=sn
            for ss=sessN
                file=fopen(fullfile(memoryDir,sprintf('%s_sess%d.txt',subj_name{s},ss)));
                datacell=textscan(file,'%f%f%f%f','HeaderLines',1,'CollectOutput',1);
                S_file=datacell{1};
                M.sn=[s;s];
                M.sessN=[ss;ss];
                group=rem(s,2);
                if group==0
                    group=2;
                end
                M.group=[group;group];
                M.seqType=[1;2];
                M.recall=S_file(:,1);
                M.recog_d=S_file(:,2);
                M.recog_bias=S_file(:,3);
                MM=addstruct(MM,M);
            end
        end
        keyboard;
        save(fullfile(memoryDir,'MemoryTests.mat'),'-struct','MM');
    case 'PLOT_memory'
        M=load(fullfile(memoryDir,'MemoryTests.mat'));
        
        figure
        subplot(1,3,1)
        lineplot([M.sessN>3 M.sessN],M.recall,'split',M.seqType,'style_thickline','plotfcn','nanmean');
        title('Recall'); xlabel('Session'); ylabel('Proportion of correct recall');
        subplot(1,3,2)
        lineplot([M.sessN>3 M.sessN],M.recog_d,'split',M.seqType,'style_thickline','plotfcn','nanmean');
        title('Recognition - d prime');  ylabel('D prime');
        subplot(1,3,3)
        lineplot([M.sessN>3 M.sessN],M.recog_bias,'split',M.seqType,'style_thickline','plotfcn','nanmean','leg',{'trained','untrained'});
        title('Recognition - bias');  ylabel('Bias');  
   
    case 'CALC_chunks'
        sn=[4:31];
        vararginoptions(varargin,{'sn'});
        D=load(fullfile(anaDir,'alldata.mat'));
        TT=[];
        for i=1:8
            t=eval(sprintf('D.pressTime%d',i))-eval(sprintf('D.pressTime%d',i-1));
            var = sprintf('IPI_%d',i);
            D.(var)=t;
        end
        % calculate within / between-chunk
        D.withinChunk1 = nanmean([D.IPI_1 D.IPI_2],2);
        D.withinChunk2 = nanmean([D.IPI_4 D.IPI_5],2);
        D.withinChunk3 = nanmean([D.IPI_7 D.IPI_8],2);
        D.betweenChunk1 = D.IPI_3;
        D.betweenChunk2 = D.IPI_6;
        D.withinChunk = nanmean([D.withinChunk1 D.withinChunk2 D.withinChunk3],2);
        D.betweenChunk = nanmean([D.betweenChunk1 D.betweenChunk2],2);
        % save the new structure
        save(fullfile(anaDir,'alldata_wChunks'),'-struct','D');
    case 'ADD_scanDay'
        sn=[4:9,11:31];
        D=load(fullfile(anaDir,'alldata_wChunks'));
        TT=[];
        for s=sn
            T=getrow(D,D.SN==s);
            % add field
            T.scanDay = zeros(size(T.SN));
            % find the actual days
            t = unique(T.day(ismember(T.blockType,[3,4,9])));
            % replace with 4 scanDays
            for i=1:4
            %    T.scanDay(T.day==t(i))=i;
                T.scanDay(T.day==t(i)&ismember(T.blockType,[3,4,9]))=i;
            end
            TT=addstruct(TT,T);
        end
        % save the new structure with scanDays added
        save(fullfile(anaDir,'alldata_wChunks'),'-struct','TT');
    case 'PLOT_chunks_learning'
        D=load(fullfile(anaDir,'alldata_wChunks'));
        D = getrow(D,D.blockType==6);
        T.chunkTime = [D.withinChunk; D.betweenChunk];
        T.chunkType = [ones(size(D.withinChunk));ones(size(D.withinChunk))*2];
        T.day       = [D.day; D.day];
        T.SN        = [D.SN; D.SN];
        T.seqType   = [D.seqType; D.seqType];
        T.MT        = [D.MT; D.MT];
        [fa ra ca] = pivottable([T.day T.SN],[T.chunkType T.seqType],T.chunkTime,'mean');
        % create a new structure
        C.chunkTime = fa(:);
        C.day       = repmat(ra(:,1),2,1);
        C.SN        = repmat(ra(:,2),2,1);
        C.chunkType = [ones(size(fa,1),1);ones(size(fa,1),1)*2];
                
        figure
        subplot(211)
        plt.line(C.day,C.chunkTime,'split',C.chunkType,'leg',{'withinChunk','betweenChunk'},'style',styTest);
        xlabel('Day');
        ylabel('Movement time');
        
        % normalise by the movement time fo the day
        [fa ra ca] = pivottable([T.day T.SN],[T.chunkType T.seqType],T.chunkTime./T.MT,'mean');
        C.chunkTime_norm = fa(:);
        subplot(212)
        plt.line(C.day,C.chunkTime_norm,'split',C.chunkType,'leg',{'withinChunk','betweenChunk'},'style',styTest);
        xlabel('Day');
        ylabel('Nomralized movement time');
        
    case 'PLOT_chunks_testDays'
        D=load(fullfile(anaDir,'testsRH'));
        
        T.chunkTime = [D.withinChunk; D.betweenChunk];
        T.chunkType = [ones(size(D.withinChunk));ones(size(D.withinChunk))*2];
        T.testDay   = [D.testDay; D.testDay];
        T.SN        = [D.SN; D.SN];
        T.seqType   = [D.seqType; D.seqType];
        [fa ra ca] = pivottable([T.testDay T.SN],[T.chunkType T.seqType],T.chunkTime,'mean');
        % create a new structure
        C.chunkTime = fa(:);
        C.testDay   = repmat(ra(:,1),6,1);
        C.SN        = repmat(ra(:,2),6,1);
        C.chunkType = [ones(size(fa,1)*3,1);ones(size(fa,1)*3,1)*2];
        C.seqType   = repmat([ones(size(fa,1),1);ones(size(fa,1),1)*3;ones(size(fa,1),1)*2],2,1);
        % do stats
        for i=unique(C.testDay)'
            fprintf('Day - %d\n',i);
            fprintf('All sequences:\n');
            anovaMixed(C.chunkTime,C.SN,'within',[C.chunkType C.seqType],{'chunkType','seqType'},'subset',C.testDay==i);
            fprintf('Trained vs. chunked:\n');
            anovaMixed(C.chunkTime,C.SN,'within',[C.chunkType C.seqType],{'chunkType','seqType'},'subset',C.seqType~=3 & C.testDay==i);
            fprintf('Trained vs. random:\n');
            anovaMixed(C.chunkTime,C.SN,'within',[C.chunkType C.seqType],{'chunkType','seqType'},'subset',C.seqType~=2 & C.testDay==i);
            fprintf('Chunked vs. random:\n');
            anovaMixed(C.chunkTime,C.SN,'within',[C.chunkType C.seqType],{'chunkType','seqType'},'subset',C.seqType~=1 & C.testDay==i);
            % follow-up t-tests per testing day
            for c=1:2  % per chunkType (effect of seqType)
                fprintf('Trained vs. chunked - chunkType %d\n',c);
                ttestDirect(C.chunkTime,[C.seqType C.SN],2,'paired','subset',C.testDay==i & C.seqType~=3 & C.chunkType==c);  
                fprintf('Trained vs. random - chunkType %d\n',c);
                ttestDirect(C.chunkTime,[C.seqType C.SN],2,'paired','subset',C.testDay==i & C.seqType~=2 & C.chunkType==c);  
                fprintf('Chunked vs. random - chunkType %d\n',c);
                ttestDirect(C.chunkTime,[C.seqType C.SN],2,'paired','subset',C.testDay==i & C.seqType~=1 & C.chunkType==c);  
            end
            for s=1:3 % per seqType (effect of chunkType)
                fprintf('SeqType %d: effect of chunkType\n',s);
                ttestDirect(C.chunkTime,[C.chunkType C.SN],2,'paired','subset',C.testDay==i & C.seqType==s);  
            end
        end
        
        figure
        plt.box(C.testDay,C.chunkTime,'split',[C.chunkType C.seqType],'leg',{'w-trained','w-chunked','w-random','b-trained','b-chunked','b-random'},'style',styTrio,'plotall',0);
        xlabel('Test day');
        ylabel('Chunk movement time');
    case 'PLOT_chunks_scanning'
        D = load(fullfile(anaDir,'alldata_wChunks'));
        D = getrow(D,D.scanDay>0 & D.seqType~=3);
        
        T.chunkTime = [D.withinChunk; D.betweenChunk];
        T.chunkType = [ones(size(D.withinChunk));ones(size(D.withinChunk))*2];
        T.scanDay   = [D.scanDay; D.scanDay];
        T.SN        = [D.SN; D.SN];
        T.seqType   = [D.seqType; D.seqType];
        [fa ra ca] = pivottable([T.scanDay T.SN],[T.chunkType T.seqType],T.chunkTime,'mean');
        % create a new structure
        C.chunkTime = fa(:);
        C.scanDay   = repmat(ra(:,1),4,1);
        C.SN        = repmat(ra(:,2),4,1);
        C.chunkType = [ones(size(fa,1)*2,1);ones(size(fa,1)*2,1)*2];
        C.seqType   = repmat([ones(size(fa,1),1);ones(size(fa,1),1)*2],2,1);
        figure
        plt.box(C.scanDay,C.chunkTime,'split',[C.chunkType C.seqType],'leg',{'w-trained','w-untrained','b-trained','b-untrained'},'style',stySeq2);
        xlabel('Scan day');
        ylabel('Chunk movement time');
        keyboard;
        % stats
        for i=unique(C.scanDay)'
            fprintf('Scan day %d\n',i);
             anovaMixed(C.chunkTime,C.SN,'within',[C.chunkType C.seqType],{'chunkType','seqType'},'subset',C.scanDay==i & ~isnan(C.chunkTime));
        end
        
    case 'BEH_ipi'
        vararginoptions(varargin,{'sn'});
        
        AllSubj = [];
        for s = sn
            
            D = load(fullfile(anaDir,['sml1_',subj_name{s},'.mat']));
            
            C = getrow(D,D.blockType==5);   % chunk training
            S = getrow(D,D.blockType==6);   % seq training
            R = getrow(D,D.blockType==11);  % test blocks
            BL = getrow(D,D.blockType==1);  % before learning
            
            indx=intersect(unique(C.day),unique(S.day)); % only take days where ChunkLearn and SeqLearn both happen
            
            C=getrow(C,ismember(C.day,indx));
            S=getrow(S,ismember(S.day,indx));
            
            chunkIPI = mean([C.pressTime1 - C.pressTime0, C.pressTime2 - C.pressTime1],2);
            seqchunkIPI = mean([S.pressTime1-S.pressTime0,S.pressTime2-S.pressTime1,S.pressTime4-S.pressTime3,S.pressTime5-S.pressTime4,S.pressTime7-S.pressTime6,S.pressTime8-S.pressTime7],2);
            seqbetchunkIPI = mean([S.pressTime3-S.pressTime2,S.pressTime6-S.pressTime5],2);
            
            rand_seqWithin = mean([R.pressTime1-R.pressTime0,R.pressTime2-R.pressTime1,R.pressTime4-R.pressTime3,R.pressTime5-R.pressTime4,R.pressTime7-R.pressTime6,R.pressTime8-R.pressTime7],2);
            rand_seqBetween = mean([R.pressTime3-R.pressTime2,R.pressTime6-R.pressTime5],2);
            
            BL_seqWithin = mean([BL.pressTime1-BL.pressTime0,BL.pressTime2-BL.pressTime1,BL.pressTime4-BL.pressTime3,BL.pressTime5-BL.pressTime4,BL.pressTime7-BL.pressTime6,BL.pressTime8-BL.pressTime7],2);
            BL_seqBetween = mean([BL.pressTime3-BL.pressTime2,BL.pressTime6-BL.pressTime5],2);
            
            
            figure(1) % chunk IPIs and trained seq
            subplot(2,length(sn),s)
            lineplot([C.day;S.day],[chunkIPI;seqchunkIPI],'split',[ones(size(C.day));ones(size(S.day))*2],'style_thickline','errorbars','shade','leg',{'chunks','sequences'},'leglocation','northeast');

            subplot(2,length(sn),length(sn)+s)
            lineplot([S.day;S.day],[seqchunkIPI;seqbetchunkIPI],'split',[ones(size(S.day));ones(size(S.day))*2],'style_thickline','errorbars','shade','leg',{'within','between'},'leglocation','northeast');

            figure(2) % trained chunks, seq
            subplot(1,length(sn),s)
            lineplot([C.day;S.day;S.day],[chunkIPI;seqchunkIPI;seqbetchunkIPI],'split',[ones(size(C.day));ones(size(S.day))*2;ones(size(S.day))*3],'style_thickline','errorbars','shade','leg',{'chunks','seq-within','seq-between'},'leglocation','northeast');
            
            figure(3) % random seq
            subplot(1,length(sn),s)
            lineplot([R.day;R.day],[rand_seqWithin;rand_seqBetween],'split', [[R.seqType;R.seqType] [ones(size(R.day));ones(size(R.day))*2]],'style_thickline2x3','leg',{'trained W','trained B','random W','random B','random chunk W', 'random chunk B'});
            
            figure(4)
            subplot(1,length(sn),s)
            lineplot([BL.BN;BL.BN],[BL_seqWithin;BL_seqBetween],'split',[ones(size(BL.day));ones(size(BL.day))*2],'style_thickline','errorbars','shade','leg',{'within','between'},'leglocation','northeast');

            
        end

    case 'CALC_dprime'
        T = load(fullfile(anaDir,'recogScores'));
        
        d1 = zeros(length(T.sn),1);
        d2 = zeros(length(T.sn),1);
        c1 = d1;
        c2 = d2;
        for s=1:length(T.sn)
            [d1(s) c1(s)] = dprime(T.hit1(s),T.falseAlarm1(s),6);
            [d2(s) c2(s)] = dprime(T.hit2(s),T.falseAlarm2(s),6);      
        end
        T.dprime1 = d1;
        T.dprime2 = d2;
        T.bias1 = c1;
        T.bias2 = c2;
        save(fullfile(anaDir,'recogScore'),'-struct','T');
   
    otherwise 
        fprintf('No such case');
        
end