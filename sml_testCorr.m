function varargout = sml_testCorr(what,varargin)
% code for estimating across-session correlation

corrLim = [0 1]; % correlation limits
nModel  = 30; % number of models to estimate between corrLim (determines resolution of corr sampling)
baseDir = cd; % same directory as for the code

switch what
    case 'runAll'
        reg = 'PMd'; % PMd or PMv;
        figOn = 1; % plot the evidence
        vararginoptions(varargin,{'reg','figOn'});
        load(fullfile(baseDir,sprintf('PCM_%s_data',reg)));
        M = sml_testCorr('constructModel','corrLim',corrLim,'nModel',nModel);
        P = pcm_fitModels(Data,M,partVec,condVec,'random','NR');
        if figOn
            sml_testCorr('plot:likelihood','Struct',P,'metric','bayesEst_cross'); % here bayesEst or bayesEst_cross
        end
        keyboard;
    case 'constructModel'
        corrLim = [0 1]; % defaults
        nModel  = 10; 
        vararginoptions(varargin,{'corrLim','nModel'});
        % specific correlation models with additional modelling of session
        corrS = linspace(corrLim(1),corrLim(2),nModel); % correlation models to assess
        % Model 1: Null model
        M{1}.type       = 'feature';
        M{1}.numGparams = 1;
        M{1}.name       = 'null';
        M{1}.Ac(:,1:12,1)  = zeros(12);
        
        % Model 2: First vs. second session
        M{2}.type       = 'feature';
        M{2}.numGparams = 2;
        M{2}.name       = 'Session';
        M{2}.Ac(:,1,1)  = [ones(6,1);zeros(6,1)];
        M{2}.Ac(:,2,2)  = [zeros(6,1);ones(6,1)];
        
        % Other models: specific correlation
        for c=1:length(corrS)
            [th2,th3]       = sml_testCorr('determine_thetaCorr','r',corrS(c));
            % Build the models
            M{c+2}.type         = 'feature';
            M{c+2}.numGparams   = 14;
            M{c+2}.name         = sprintf('SpecCorr_%1.3f',corrS(c));
            M{c+2}.Ac(:,1,1)    = [ones(6,1);zeros(6,1)];
            M{c+2}.Ac(:,2,2)    = [zeros(6,1);ones(6,1)];
            % for sequence-specific modelling- one parameter per session
            for i=1:6
                A=zeros(6);
                A(i,i)=1;
                M{c+2}.Ac(:,3:8,i+2)     = [A;zeros(6)];      % Unique sess1 sequence patterns
                M{c+2}.Ac(:,9:14,8+i)    = [zeros(6);A*th2];  % Unique sess2 sequence pattterns
                M{c+2}.Ac(:,3:8,8+i)     = [zeros(6);A*th3];  % Correlation sess1-sess2
            end;
        end
        varargout{1}=M;
    case 'determine_thetaCorr'
        vararginoptions(varargin,{'r'}); % which correlation to provide thetas for
        % determine thetas such that correlation is of specific value
        if r==1
            th2=0;
            th3=1;
        else
            th2=1;
            th3=sqrt((r^2)/(1-r^2)); % equation for determining thetas
        end
        varargout{1}=th2;
        varargout{2}=th3;
    case 'plot:likelihood'
        vararginoptions(varargin,{'Struct','metric'});
        T = Struct; % pcm structure
        % subtract the session effect
        T.(metric) = bsxfun(@minus,T.(metric),T.(metric)(:,2));
        T.(metric) = T.(metric)(:,2:end);
        % subtract the mean evidence across all corr models
        T.metric2 = bsxfun(@minus,T.(metric),mean(T.(metric),2)); 
        % concatenate
        D.model     = kron((1:size(T.(metric),2))',ones(size(T.(metric),1),1));      
        D.(metric)  = T.(metric)(:);
        D.metric2   = T.metric2(:);
        figure
        subplot(121)
        lineplot(D.model,D.(metric)); title(metric);
        subplot(122)
        lineplot(D.model,D.metric2); title('mean corrModel evidence subtracted');       
end
        
end
% local function for fitting PCM
function T = pcm_fitModels(Data,M,partVec,condVec,runEffect,algorithm)
    % --------------------------------------
    % Run group fit in 3 stages 
    % 1) first minimize (to provide decent starting values)
    % 2) then NR 
    % 3) then minimize again (if any furhter improvement)
     [~,theta_hat,~,~]  = pcm_fitModelGroup(Data,M,partVec,condVec,'runEffect',runEffect,'fitScale',1,'fitAlgorithm','minimize');
     fprintf('Group fit with minimize algorithm done.\n');
     [~,theta_hat,~,~]  = pcm_fitModelGroup(Data,M,partVec,condVec,'runEffect',runEffect,'fitScale',1,'fitAlgorithm','NR','theta0',theta_hat);
     fprintf('Group fit with NR algorithm done.\n');
     [T,theta_hat,~,~]  = pcm_fitModelGroup(Data,M,partVec,condVec,'runEffect',runEffect,'fitScale',1,'fitAlgorithm','minimize','theta0',theta_hat);
     fprintf('Group fit with minimize algorithm done.\n');
    % Run crossvalidated group fit
    [Tcross,theta_hat]       = pcm_fitModelGroupCrossval(Data,M,partVec,condVec,'runEffect',runEffect,'groupFit',theta_hat,'fitScale',1,'fitAlgorithm',algorithm);
    T.theta = theta_hat;
    fprintf('Crossvalidated fit with %s algorithm done.\n',algorithm);
    T.cross_likelihood  = Tcross.likelihood;
    T.bayesEst          = bsxfun(@minus,T.likelihood,T.likelihood(:,1));
    T.bayesEst_cross    = bsxfun(@minus,T.cross_likelihood,T.cross_likelihood(:,1));
end