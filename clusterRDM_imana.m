function varargout = clusterRDM_imana(what, varargin);
%%
% Access to:
%   - glm directory
%   - surface caret directory
%       - pre-calculated distance result
%   - region of interest directory
%   -
% Does:
%   - clustering using RDM (or second-moment matrix)
%   - some visualization
%
%
% ayokoi (Jul 2018)


%% Data set
datasets = {'sh1'; 'sh2'; 'bmw1'};
glmtypes = [3,2,10];
subjects = {[2:6,8,10:13 14:15]; [1:9]; [1:7]};
subj_names = {
    {'p01','t01','t02','s03','s04','s07',... % s05
    's09','s10','s14','s15', 's17','s19','s22','s23','s25'};
    {'s01','s02','s03','s04','s05','s06','s07','s08','s09'};
    {'s02','s04','s05','s06','s07','s09','s10'};
    };
distanceforms = {
    's%s_dist_raw_glm3.metric';
    's%s_alldistance.metric';
    '%s_clusterRDM_D.metric'; % only defined on reduced searchlight, hence under _reduced folder
    };
distloci = {'','','_reduced'};
groupdistnames = {
    '%s.summary_meanDist_smooth.metric'; % col: 1,2 (ssqrt-ed)
    '%s.DistSummary.metric'; % col: 1-4, 9-12 (separate threshold for single/multi finger active) (ssqrt-ed)
    '%s.dist_summary.metric'; % col: 1-4, 7-10 (overall mean) (ssqrt-ed)
    };
groupcols = {[1],[1,2],[1:4]};
maskthresholds = {[0.03]; [0.045,0.025]; [0.05,0.05,0.03,0.03]};
maskuse = {'sh1','sh1','bmw1'};
%% fresurfer and caret
atlasA      = {'x'};
atlasname   = {'fsaverage_sym'};
hemispheres = {'LeftHem','RightHem'};
hemName = hemispheres;
hemtags = {'lh','rh'};
hem = hemtags;
%% Directory
dataDir = '/Volumes/G_Thunderbolt/Yokoi_Research/data/';
baseDir = fullfile(dataDir, 'RepresentationalClustering'); dircheck(baseDir);
baseDir = fullfile(baseDir, 'clusterRDM1'); dircheck(baseDir);
analyzeDir = fullfile(baseDir, 'analyze'); dircheck(analyzeDir);
figDir = fullfile(baseDir, 'figures'); dircheck(figDir);
surfaceDir = fullfile(baseDir, 'surfaceCaret', 'fsaverage_sym'); dircheck(analyzeDir);
dircheck(fullfile(surfaceDir, hemName{1}));
dircheck(fullfile(surfaceDir, hemName{2}));
% get directory for each peoject
for p=1:numel(datasets)
    str = sprintf('%s = getDirectory(datasets{p}, dataDir);', datasets{p});
    eval(str);
end

%% Main
switch (what)
    case 'SL_make_distmask'        % Create distance-thresholding mask defined on full-surface from .metric file of mean distance
        pj      = varargin{1}; % project
        smoothparam = [25];
        maxcsize = 200;
        vararginoptions(varargin(2:end),{'smoothparam','maxcsize'});
        
        % loop over projects
        for p=pj
            project = datasets{p};
            sn = subjects{p};
            subj_name = subj_names{p};
            
            % Get project-wise directory structure
            str=sprintf('Dir=%s;', datasets{p});
            eval(str);
            caretDir = Dir.caretDir;
            
            % Load group mean distance data, smooth it and save to project directory
            groupdist= groupdistnames{p};
            col = groupcols{p};
            for h=[1,2]
                caretGroupDir = fullfile(caretDir, 'fsaverage_sym', hemName{h});
                summaryfile = fullfile(caretGroupDir, sprintf(groupdist, hem{h}));
                M = caret_load(summaryfile);
                dist = M.data(:,col);
                
                %coord = fullfile(caretGroupDir, sprintf('%s.INFLATED.coord', hem{h}));
                coord = fullfile(caretGroupDir, sprintf('%s.WHITE.coord', hem{h}));
                topo = fullfile(caretGroupDir, sprintf('%s.CLOSED.topo', hem{h}));
                surface = caret_getsurface(coord, topo);
                surface.Edges=surface.Tiles;
                
                idx = false(length(M.index),1);
                for c=1:length(col)
                    % save distance data
                    C = M;
                    C.data = dist(:,c);
                    C.num_cols = length(col(c));
                    C.column_name = M.column_name(col(c));
                    C.column_color_mapping = M.column_color_mapping(col(c),:);
                    fname = fullfile(surfaceDir, hemName{h}, sprintf('%s.%s_dist_%d.metric', hem{h}, project, c));
                    caret_savemetric(fname, C);
                    
                    % smooth
                    fname_smooth = caret_smooth(fname, 'coord', coord, 'topo', topo, 'iterations',smoothparam);
                    
                    % apply threshold
                    D = caret_load(fname_smooth{1});
                    idx = idx + D.data >= maskthresholds{p}(c);
                end
                % find clusters of mask to remove small scums
                cidx = caret_clusters(surface, idx);
                cs = unique(cidx);
                for i=2:length(cs)
                    curridx = cidx==cs(i);
                    if sum(curridx)<maxcsize % if cluster is small
                        idx(curridx) = 0;
                    end
                end
                
                % save as .metric file
                C = D;
                C.data = double(idx);
                C.num_cols = 1;
                C.column_name = {'distancemask'};
                fname = fullfile(surfaceDir, hemName{h}, sprintf('%s.%s_distmask.metric', hem{h}, project));
                caret_savemetric(fname, C);
                
                % smooth again
                fname_smooth = caret_smooth(fname, 'coord', coord, 'topo', topo, 'iterations',smoothparam);
                S = caret_load(fname_smooth{1});
                
                % save as .mat file
                mask(:,h) = double(S.data>0.1);idx;
                
                % update  .metric file
                C.data = double(mask(:,h));
                C.num_cols = 1;
                C.column_name = {'distancemask'};
                fname = fullfile(surfaceDir, hemName{h}, sprintf('%s.%s_distmask2.metric', hem{h}, project));
                caret_savemetric(fname, C);
            end
            fname = fullfile(analyzeDir, sprintf('%s_distmask.mat', project));
            save(fname, 'mask');
        end
        
        % make joint mask
        jmask = zeros(length(idx),2);
        for p=1:numel(datasets)
            project=datasets{p};
            fname = fullfile(analyzeDir, sprintf('%s_distmask.mat', project));
            load(fname);
            jmask = jmask+logical(mask);
        end
        jmask=jmask>0;
        for h=1:2
            fname = fullfile(surfaceDir, hemName{h}, sprintf('%s.joint_distmask.metric', hem{h}));
            C.data = double(jmask(:,h));
            caret_savemetric(fname,C);
            
            % smooth
            caretGroupDir = fullfile(caretDir, 'fsaverage_sym', hemName{h});
            %coord = fullfile(caretGroupDir, sprintf('%s.INFLATED.coord', hem{h}));
            coord = fullfile(caretGroupDir, sprintf('%s.WHITE.coord', hem{h}));
            topo = fullfile(caretGroupDir, sprintf('%s.CLOSED.topo', hem{h}));
            fname_smooth = caret_smooth(fname, 'coord', coord, 'topo', topo, 'iterations',smoothparam);
            S = caret_load(fname_smooth{1});
            jmask(:,h) = S.data>0;
        end
        fname = fullfile(analyzeDir, sprintf('joint_distmask.mat'));
        save(fname,'jmask');
        
        varargout = {};
    case 'SL_make_distmask_nii'        % TODO: make this fit into current function Make .img mask based on distance
        sn          = varargin{1};
        threshold=0.03;
        atlas       = 1;
        glm         = 3;
        vararginoptions({varargin{2:end}},{'atlas','glm','threshold'});
        
        % load group mean distance
        for h=1:2
            groupDir = fullfile(caretDir,'fsaverage_sym',hemName{h});
            cd(groupDir);
            
            meanDist    = caret_load([hem{h} '.summary_meanDist_smooth.metric']);
            indx{h}          = meanDist.data(:,1)>threshold;
        end
        
        for s=sn
            if exist(fullfile(regDir,subj_name{s}))
                cd(fullfile(regDir,subj_name{s}));
            else
                mkdir(fullfile(regDir,subj_name{s}));
                cd(fullfile(regDir,subj_name{s}));
            end
            glmDir = fullfile(baseDir,glmName{glm});
            
            for h=1:2
                caretSubjDir = fullfile(caretDir,['x' subj_name{s}]);
                
                % Nodes whose mean distance is above threshold
                numNodes    = sum(indx{h});
                
                R{h}.type = 'surf_nodes';
                R{h}.white=fullfile(caretSubjDir,hemName{h},[hem{h} '.WHITE.coord']);
                R{h}.pial=fullfile(caretSubjDir,hemName{h},[hem{h} '.PIAL.coord']);
                R{h}.topo=fullfile(caretSubjDir,hemName{h},[hem{h} '.CLOSED.topo']);
                R{h}.linedef=[5,0,1];
                R{h}.image=fullfile(glmDir,subj_name{s},'mask.img');
                R{h}.name=[subj_name{s} '_meanDistThredholded_' hem{h}];
                R{h}.location=find(indx{h});
            end;
            R = region_calcregions(R);
            
            % Apply resultant [x y z] to mask.img to make new thresholded
            % individual mask
            V = spm_vol(R{1}.image);
            V2 = V;
            V.mask = spm_read_vols(V);
            V2.mask = zeros(V.dim);
            for h=1:2
                V2.mask(R{h}.linvoxidxs) = 1;
            end
            figure(1);
            for i=1:V.dim(3)
                subplot(4,8,i); imagesc(V.mask(:,:,i)');
            end
            figure(2);
            for i=1:V.dim(3)
                subplot(4,8,i); imagesc(V2.mask(:,:,i)');
            end
            % Write volume
            V2.fname = fullfile(glmDir,subj_name{s},'distmask.img');
            V2.descrip = 'spm mask thresholded by mean distance';
            V2 = spm_write_vol(V2, V2.mask);
            
            fprintf('\n subject no. %d done.\n',s);
        end
        
        varargout = {};
    case 'SL_define_nodes'         % Subsample "fsaverage_sym" surface to make new reduced surface for group pca searchlight
        pj = varargin{1};
        separation = 7;
        fig     = 1;
        shownum = 0;
        mask = '%s_distmask.mat';        
        vararginoptions(varargin(2:end),{'separation','fig','shownum','mask'});
        
        % fsaverage_sym is the same acros projects
        caretDir = eval(sprintf('%s.caretDir',datasets{1}));
        for p=pj;%1:numel(datasets)
            project = datasets{p};
            
            % load mask if specified
            if ~isempty(mask)
                D=load(fullfile(analyzeDir, sprintf(mask, datasets{p})));
                field=fieldnames(D);
                maskdata = logical(D.(field{1}));
            end
            
            % Load flat map and mask
            for h=1:2
                if ~isempty(mask)
                    idx = maskdata(:,h);
                else
                    idx = [];
                end
                if fig;figure('units','centimeters','position',[5,5,15,15]);end;
                
                groupDir = fullfile(caretDir,'fsaverage_sym',hemName{h});
                cd(groupDir);
                
                % -- flat map
                Flat = caret_load([hem{h} '.FLAT.coord']);
                xrange      = [floor(min(Flat.data(:,1))),ceil(max(Flat.data(:,1)))];
                yrange      = [floor(min(Flat.data(:,2))),ceil(max(Flat.data(:,2)))];
                
                % -- border files
                Border = caret_load([hem{h} '.display.borderproj']);
                
                % Display mask on flat map
                switch fig
                    case 1
                        % flat map
                        plot(Flat.data(:,1),Flat.data(:,2),'.','color',[0.8 0.8 0.8]);
                        hold on; axis equal;
                        
                        % border dots
                        switch (h)
                            case 1
                                borderidx   = cat(1,Border.Border.vertex);
                                borderX     = Flat.data(borderidx(:),1);
                                borderY     = Flat.data(borderidx(:),2);
                            case 2
                                borderX     = -borderX;
                        end
                        plot(borderX,borderY,'k.');
                        
                        caxis([0.3,0.7]);
                        set(gca,'xlim',xrange,'ylim',yrange);
                        axis off
                        title(hemName{h})
                end
                
                % Generate hexagonal grid
                [nodeidx, edgeidx, sp, XV, YV] = clusterRDM_hexagon_flat('Flat', Flat, 'separation', separation,...
                    'mask', idx);
                
                switch fig
                    case 1
                        plot(XV,YV,'r');
                        for i=1:numel(nodeidx)
                            fillHex(Flat.data(nodeidx(i),1),Flat.data(nodeidx(i),2),1,...
                                'C',[0 0 1]);
                            if shownum
                                text(Flat.data(nodeidx(i),1),Flat.data(nodeidx(i),2),num2str(i),...
                                    'horizontalalignment','center','verticalalignment','middle',...
                                    'fontsize',10);
                            end
                        end
                        set(gca,'xlim',xrange,'ylim',yrange);
                        axis equal;
                end;
                
                Nodes{h} = nodeidx;
                Edges{h} = edgeidx;
                Separations{h} = sp;
            end
            % save result
            save(fullfile(surfaceDir, sprintf('%s_hexNodes.%dmmSeparation.mat',project, round(separation))),...
                'Nodes', 'Edges', 'Separations');
        end
        
        varargout = {Nodes,Separations};
    case 'SL_define_searchlight' % Re-define sparse searchlight based on re-defined node
        pj      = varargin{1}; % project
        overlap = 0.2; % 0~1 (1 was too large)
        separation = 7; % mm
        mask = []; % file name (.img or .nii, saved under 1st-level GLM directory for each subj)
        nVox = 150; % if not explicitly given, automatically calculated using separation and overlap parameters
        vararginoptions(varargin(2:end),{'overlap','separation','nVox'});
        
        % loop over projects
        for p=pj
            project = datasets{p};
            sn = subjects{p};
            subj_name = subj_names{p};
            
            % Get project-wise directory structure
            str=sprintf('Dir=%s;', datasets{p});
            eval(str);
            caretDir = Dir.caretDir;
            glmDir  = Dir.glmDir;
            
            % Load ultra-reduced node indices (created in 'define_nodes')
            load(fullfile(surfaceDir, sprintf('%s_hexNodes.%dmmSeparation.mat',project, round(separation))));
            % Nodes, Edges, Separations
            
            % Calc radius and nVox
            radius  = separation*(0.5*overlap+0.5);
            if isempty(nVox)
                nVox    = 160 * (radius/12)^2; % 2D approximation
                % adjust digit
                nVox = 10*ceil(nVox/10);
            end
            
            % Loop over subjects
            for s=sn
                for h=1:numel(hemtags)
                    % load caret files
                    caret_subjDIR   = fullfile(caretDir,[atlasA{1},subj_name{s},''],hemispheres{h});
                    caret_subjDIRr  = fullfile(caretDir,[atlasA{1},subj_name{s},'_reduced'],hemispheres{h});
                    coord_pial      = caret_load(fullfile(caret_subjDIR, sprintf('%s.PIAL.coord',hemtags{h})));
                    coord_white     = caret_load(fullfile(caret_subjDIR, sprintf('%s.WHITE.coord',hemtags{h})));
                    topo            = caret_load(fullfile(caret_subjDIR, sprintf('%s.CLOSED.topo',hemtags{h})));
                    
                    % get surface roi
                    
                    % resolve mask
                    if isempty(mask) % try default mask created through 1st-level GLM
                        ref = fullfile(glmDir, subj_name{s},'mask.img');
                        if ~exist(ref, 'file'); % check if we have either .img or .nii file
                            ref = fullfile(glmDir, subj_name{s},'mask.nii');
                        end
                    else
                        ref = fullfile(glmDir, subj_name{s},mask);
                    end
                    if ~exist(ref, 'file');
                        warning('%s does not exist.', ref);
                        keyboard();
                    end
                    refV            = spm_vol(ref);
                    epiInfo.mat     = refV.mat;
                    epiInfo.dim     = refV.dim;
                    epiInfo.mask    = spm_read_vols(refV);
                    centernode      = Nodes{h};
                    
                    num_nodes       = coord_white.num_nodes;
                    
                    %- run surfing_voxelselection
                    [LI,voxmin,voxmax,vORr] = ...
                        surfing_voxelselection(coord_white.data',...
                        coord_pial.data',...
                        topo.data',...
                        [radius nVox],...
                        epiInfo,...
                        centernode,...
                        [5,0,1]);
                    % check validity
                    badidx  = false(size(vORr));%isnan(vORr); don't change SL here otherwise you'll pass unequal SL nodes to the searchlight function.
                    LI      = LI(~badidx);
                    voxmin  = voxmin(~badidx,:);
                    voxmax  = voxmax(~badidx,:);
                    vORr    = vORr(~badidx,:);
                    nodeID  = centernode(~badidx)';
                    
                    %- save searchlight def in .mat format
                    if ~exist(caret_subjDIRr);
                        mkdir(caret_subjDIRr);
                    end
                    LI = LI';
                    save(fullfile([caret_subjDIRr], sprintf('clusterRDM_masked_roi_%dmmSeparation_%dvox.mat',separation,nVox)),...
                        'LI','voxmin','voxmax','vORr','nodeID','num_nodes');
                    
                    %- save searchlight def in .metric format
                    M                   = caret_struct('metric','data',vORr);
                    M.data              = NaN(size(coord_white.data,1),1);
                    M.data(nodeID)      = vORr;
                    M.num_rows          = length(M.data);
                    SLname = sprintf('clusterRDM_masked_roi_%dmmSeparation_%dvox.metric',separation,nVox);
                    caret_save(fullfile([caret_subjDIRr], SLname),M);
                end
            end
        end
        
        varargout = {SLname};
        
    case 'distance_import_design'        % Import task design (condition, partiton, etc.) information
        pj = varargin{1};
        vararginoptions(varargin(2:end),{''});
        fileform = 'ROI_pwhBeta.glm%d.mat';
        
        % loop over projects
        for p=pj
            project = datasets{p};
            sn = subjects{p};
            subj_name = subj_names{p};
            distanceform = distanceforms{p};
            glmtype = glmtypes(p);
            % Get project-wise directory structure
            str=sprintf('Dir=%s;', datasets{p});
            eval(str);
            caretDir = Dir.caretDir;
            caretGroupDir = fullfile(caretDir, 'fsaverage_sym');
            glmDir  = Dir.glmDir;
            
            % Load prewhitened beta
            pwhname = sprintf(fileform, glmtype);
            load(fullfile(Dir.regDir,pwhname));
            
            % Get design and give common variable names
            switch (project)
                case 'sh1'
                    design.condition = Design.condition;
                    design.partition = Design.partition;
                    design.subj = Design.subj;
                case 'sh2'
                    design.condition = Design.sequence;
                    design.partition = Design.partition;
                    design.subj = Design.subj;
                    design.activepassive = Design.condition;
                case 'bmw1'
                    design.condition = Design.condition;
                    design.partition = Design.partition;
                    design.subj = Design.subj;
                    design.hand = Design.hand;
            end
            
            % Save
            fname = sprintf('%s_design.mat',project);
            save(fullfile(analyzeDir, fname), '-struct', 'design');
        end
    case 'distance_recalc'      % Re-calculate RDMs only on the sparse searchlight nodes
        pj = varargin{1};
        SLname  = 'clusterRDM_masked_roi_7mmSeparation_150vox.mat'; %'masked_roi_250vox'; % 90
        normmode = 'overall'; % multivariate noise-normalization parameter ('runwise' or 'overall')
        blocksize= 5e07; % default value is 7e05
        vararginoptions(varargin(2:end),{'SLname','blocksize','normmode'});
        
        SLfun = @clusterRDM_calcG_SL;
        
        home = pwd;
        % loop over projects
        for p=pj
            project = datasets{p};
            sn = subjects{p};
            subj_name = subj_names{p};
            % Get project-wise directory structure
            str=sprintf('Dir=%s;', datasets{p});
            eval(str);
            caretDir = Dir.caretDir;
            caretGroupDir = fullfile(caretDir, 'fsaverage_sym');
            glmDir  = Dir.glmDir;
            
            % load design
            Design = load(fullfile(analyzeDir, sprintf('%s_design.mat', project)));
            Ncondition = length(unique(Design.condition));
            Ng = Ncondition*(Ncondition-1)/2+Ncondition; % for second-moment
            Nd = Ncondition*(Ncondition-1)/2; % for squared euclidean distance
            
            for s=sn
                glmDirSubj = fullfile(glmDir, subj_name{s});
                caretDirSubj = fullfile(caretDir, ['x', subj_name{s}, '_reduced']); % sparce surface def
                cd(glmDirSubj);
                
                % Load SPM.mat
                load(fullfile(glmDirSubj,'SPM.mat'));
                % if you are working on different directory structure from
                % which SPM.mat is originally created;
                % SPM = spmj_move_rawdata(SPM,fullfile(baseDir,'imaging_data',subj_name{s},'sess1'));
                
                % Load condition and partition data
                design = getrow(Design, Design.subj==s);
                Cond = design.condition';
                Part = design.partition';
                
                % Cond: column, Part: column
                params = {SPM,Cond,Part,'normmode',normmode};
                
                % Run searchlighg for each hemisphere, separately
                for h=1:2
                    fprintf('%s %s %s ...', project, subj_name{s}, hemName{h})
                    
                    % Define output files
                    metric_out = {};
                    for i=1:Ng
                        metric_out{i} = fullfile(caretDirSubj, hemName{h}, sprintf('%s_clusterRDM_G.metric', subj_name{s}));
                    end
                    for i=1:Nd
                        metric_out{end+1} = fullfile(caretDirSubj, hemName{h}, sprintf('%s_clusterRDM_D.metric', subj_name{s}));
                    end
                    % Load searchlight definition and run
                    SL  = load(fullfile(caretDirSubj, hemName{h}, SLname));
                    clusterRDM_searchlight(SL, SPM.xY.P, metric_out, SLfun,...
                        'params', params, 'isNP',1, 'idealblock', blocksize);
                    
                    fprintf('... done.\n');
                end
            end
        end
        cd(home);
    case 'distance_embed_fullsearchlight' % Embed reduced-searchlight result into full-searchlight (may waste memory)
        pj = varargin{1};
        SLname  = 'clusterRDM_masked_roi_7mmSeparation_150vox.mat'; %'masked_roi_250vox'; % 90
        vararginoptions(varargin(2:end),{'SLname','blocksize','normmode'});
        
        home = pwd;
        % loop over projects
        for p=pj
            project = datasets{p};
            sn = subjects{p};
            subj_name = subj_names{p};
            % Get project-wise directory structure
            str=sprintf('Dir=%s;', datasets{p});
            eval(str);
            caretDir = Dir.caretDir;
            caretGroupDir = fullfile(caretDir, 'fsaverage_sym');
            glmDir  = Dir.glmDir;
            
            for s=sn
                caretDirSubj = fullfile(caretDir, ['x', subj_name{s}]); % full surface def
                caretDirSubj_reduced = fullfile(caretDir, ['x', subj_name{s}, '_reduced']); % sparce surface def
                
                %
                for h=1:2
                    fprintf('%s %s %s ...', project, subj_name{s}, hemName{h})
                    
                    % Load searchlight definition
                    SL  = load(fullfile(caretDirSubj_reduced, hemName{h}, SLname));
                    
                    % Load Flat.coord
                    flat = fullfile(caretGroupDir, hemName{h}, sprintf('%s.Flat.coord', hem{h}));
                    F = caret_load(flat);
                    
                    % Load distance file
                    metricin = fullfile(caretDirSubj_reduced, hemName{h}, sprintf('%s_clusterRDM_D.metric',subj_name{s}));
                    if ~exist(metricin, 'file');
                        warning('\n%s does not exist. Stopping operation.\n', metricin);
                        return;
                    end
                    R = caret_load(metricin);
                    
                    % Define output files
                    metricout = fullfile(caretDirSubj, hemName{h}, sprintf('%s_clusterRDM_D.metric', subj_name{s}));
                    O = R;
                    O.num_rows = F.num_nodes;
                    O.index = F.index;
                    O.data = NaN(size(F.data,1), O.num_cols);
                    O.data(SL.nodeID,:) = R.data;
                    
                    % Embed data
                    caret_savemetric(metricout, O);
                    
                    fprintf('... done.\n');
                end
            end
        end
        cd(home);
    case 'distance_extract'     % Extract RDMs for sparse searchlight nodes fron existing .metric files (if they exist)
        pj = varargin{1}; %
        SLname  = 'clusterRDM_masked_roi_7mmSeparation_150vox'; %'masked_roi_250vox'; % 90
        alpha = 0.7; % proportion of subjects
        vararginoptions(varargin(2:end),{'SLname','alpha'});
        
        % loop over projects
        for p=pj
            project = datasets{p};
            sn = subjects{p};
            subj_name = subj_names{p};
            distanceform = distanceforms{p};
            distlocation = distloci{p};
            % Get project-wise directory structure
            str=sprintf('Dir=%s;', datasets{p});
            eval(str);
            caretDir = Dir.caretDir;
            caretGroupDir = fullfile(caretDir, 'fsaverage_sym');
            glmDir  = Dir.glmDir;
            
            T = []; oldnroi=0;
            for h=1:2
                % load surface info to add x-y coordinate info
                flatname  = [hem{h} '.FLAT.coord'];
                F = caret_load(fullfile(caretGroupDir, hemName{h}, flatname));
                % load shape info to add depth info
                shape = [hem{h} '.surface_shape'];
                S = caret_load(fullfile(caretGroupDir, hemName{h}, shape));
                Nnodes_full = size(F.data,1);
                
                % get node membership to Desikan atlas
                Desikan = caret_load(fullfile(caretGroupDir, hemName{h}, [hem{h} '.desikan.paint']));
                for i=1:numel(Desikan.paintnames);
                    Desikan.paintnames{i} = [hem{h}, '.', Desikan.paintnames{i}];
                end
                
                for s=sn
                    % Searchlight definition
                    caret_subjDIR   = fullfile(caretDir,[atlasA{1},subj_name{s},'_reduced'],hemispheres{h});
                    Searchlight        = fullfile(caret_subjDIR, [SLname,'.mat']);
                    
                    % Distance .metric file for individuals (calculated either on full- or reduced-searchlight)
                    caret_subjDIR   = fullfile(caretDir,[atlasA{1},subj_name{s}, distlocation],hemispheres{h});
                    Metric = fullfile(caret_subjDIR,sprintf(distanceform,subj_name{s}));
                    
                    % Load data
                    SL  = load(Searchlight);
                    M   = caret_load(Metric);
                    
                    % find nodes where result exists
                    nodeID = SL.nodeID;
                    roi = [1:length(nodeID)]+oldnroi;
                    
                    % Check size
                    Nnodes = length(M.index);
                    Nnodes_SL = length(SL.nodeID);
                    if Nnodes<Nnodes_full
                        if Nnodes~=Nnodes_SL
                            error('Invalid size (%d nodes in searchlight definition, but %d nodes in .metric file', Nnodes_SL, Nnodes);
                        end
                        M.index = SL.nodeID;
                        t.data = ssqrt(M.data);
                    elseif Nnodes==Nnodes_full
                        t.data = ssqrt(M.data(nodeID, :)); % dist
                    end
                    
                    % concat data
                    t.nodeID = nodeID;
                    t.coord = [F.data(nodeID,1:2), S.data(nodeID,2)];
                    t.roi = roi';
                    t.roiname = Desikan.paintnames(Desikan.data(nodeID))';
                    t.desikanid = Desikan.data(nodeID);
                    t.hemis = repmat(h, length(t.nodeID),1);
                    t.subj = repmat(s, length(t.nodeID),1);
                    
                    T = addstruct(T,t);
                end
                oldnroi = max(roi);
            end;
            
            % save
            fname = fullfile(analyzeDir, sprintf('%s_distance_%s.mat', project, SLname));
            save(fname, '-struct', 'T');
        end
        
        varargout = {T};
    case 'distance_split_sh2' % Split distance into active and , passive conditions
        SLname  = 'clusterRDM_masked_roi_7mmSeparation_150vox'; %'masked_roi_250vox'; % 90
        vararginoptions(varargin(1:end),{'SLname'});
        
        % load distance
        fname = fullfile(analyzeDir, sprintf('sh2_distance_%s.mat', SLname));
        T = load(fname);
        
        % split and save
        indices = {[1:5],... % active single digit sequences
            [6:11],... % active multi digit sequences
            [12:16],... % passive single digit sequences
            [17:22]}; % passive multi digit sequences
        names = {'activesingle','activemulti','passivesingle','passivemulti'};
        for cond=1:numel(names)
            idx = false(22);
            idx(indices{cond}, indices{cond}) = true;
            M = T; % copy
            M.data = [];
            for i=1:size(T.data,1)
                dist = squareform(T.data(i,:));
                M.data(i,:) = squareform(reshape(dist(idx), numel(indices{cond}), numel(indices{cond})));
            end
            fname = fullfile(analyzeDir, sprintf('sh2_distance_%s_%s.mat', SLname, names{cond}));
            save(fname, '-struct', 'M');
        end
    case 'distance_split_bmw1' % Split distance into unimanual and bimanual parts
        SLname  = 'clusterRDM_masked_roi_7mmSeparation_150vox'; %'masked_roi_250vox'; % 90
        vararginoptions(varargin(1:end),{'SLname'});
        
        % load distance
        fname = fullfile(analyzeDir, sprintf('bmw1_distance_%s.mat', SLname));
        T = load(fname);
        
        % split and save
        indices = {[1:12],... % unimanual movements
            [13:48]}; % bimanual movements
        names = {'unimanual','bimanual'};
        for cond=1:numel(names)
            idx = false(48);
            idx(indices{cond}, indices{cond}) = true;
            M = T; % copy
            M.data = [];
            for i=1:size(T.data,1)
                dist = squareform(T.data(i,:));
                M.data(i,:) = squareform(reshape(dist(idx), numel(indices{cond}), numel(indices{cond})));
            end
            fname = fullfile(analyzeDir, sprintf('bmw1_distance_%s_%s.mat', SLname, names{cond}));
            save(fname, '-struct', 'M');
        end
    case 'distance_visualize' % visualize clustering result on surface
        pj = varargin{1}; % project index
        SLname  = 'clusterRDM_masked_roi_7mmSeparation_150vox'; %'masked_roi_250vox'; % 90
        overlap = 0;
        condition = ''; % '_active' or '_passive' for sh2, '_unimanual' or 'bimanual' for bmw1
        isSave = 0;
        resolution = '-r200';
        threshold=[];
        smooth = 0;
        vararginoptions(varargin(2:end),{'SLname','condition','isSave','resolution','threshold','smooth'});
        
        for p=pj;
            project = datasets{p};
            sn = subjects{p};
            subj_name = subj_names{p};
            Nsubj = numel(sn);
            
            % Get project-wise directory structure
            str=sprintf('Dir=%s;', datasets{p});
            eval(str);
            caretDir = Dir.caretDir;
            glmDir  = Dir.glmDir;
            caretGroupDir = fullfile(caretDir, 'fsaverage_sym');
            
            % Define plot range
            % Set mask (background region on the basic flat map)
            for pp=1:length(datasets)
                % Load distance mask for background
                %mask = load(fullfile(analyzeDir, 'joint_distmask.mat'));
                mask = load(fullfile(analyzeDir, sprintf('%s_distmask.mat', datasets{pp})));
                fn = fieldnames(mask);
                bgmask(:,:,pp) = mask.(fn{1});
            end;
            masks = bgmask;
            bgmask=sum(bgmask,3)>0;
            for h=1:2
                [Flat] = getSurface(fullfile(Dir.caretDir,'fsaverage_sym'), h);
                x(1,h) = min(Flat.data(bgmask(:,h),1));
                x(2,h) = max(Flat.data(bgmask(:,h),1));
                y(1,h) = min(Flat.data(bgmask(:,h),2));
                y(2,h) = max(Flat.data(bgmask(:,h),2));
            end;
            xrange = [floor(min([x(1,1),-x(2,2)])), ceil(max([x(2,1), -x(2,1)]))];
            yrange = [floor(min([y(1,1),y(1,2)])), ceil(max([y(2,1), y(2,2)]))];
            plotrange = {{xrange,yrange},{-xrange([2,1]),yrange}};
            
            % Load distance data
            distfname = fullfile(analyzeDir, sprintf('%s_distance_%s%s.mat', project, SLname,condition));
            T = load(distfname);
            
            % Load node assignment
            overlapstr = {'unoverlapped','overlappled'};
            fname = fullfile(analyzeDir, sprintf('%s_%s_sharednodes_%s.mat',project, SLname,overlapstr{overlap+1}));
            load(fname);
            
            % Threshold
            if isempty(threshold)
                switch (condition)
                    case '_activesingle'
                        i=1;
                    case '_activemulti'
                        i=2;
                    case '_passivesingle'
                        i=1;
                    case '_passivemulti'
                        i=2;
                    case '_unimanual'
                        i=1;
                    case '_bimanual'
                        i=3;
                    otherwise
                        i=1;
                end
                threshold = [maskthresholds{p}(i)]*0.4;
            end
            
            % Map mean distance
            T.data = nanmean(T.data,2);
            [maskout, fighandle, scalehandle, window]=clusterRDM_imana('vis_map', T, SLnodes, caretGroupDir,...
                'plotUnderlay', 1,...
                'makeFig', 1,...
                'MAP', parula(100),...
                'modelname', project,...
                'printScale',1,...
                'threshold',threshold,...
                'plotrange',plotrange,...
                'bgmask',masks(:,:,p),...
                'mask',masks(:,:,p),...
                'smooth',smooth);
            
            % Save map&scale
            figname = fullfile(figDir, project, sprintf('%s%s_DistanceMap_%s',project, condition, SLname));
            mySaveFig(fighandle, figname, isSave, '-dpng',resolution);
            figname = fullfile(figDir, sprintf('%s%s_DistanceScale_%s',project, condition, SLname));
            mySaveFig(scalehandle, figname, isSave, '-dpng',resolution);
            mySaveFig(scalehandle, figname, isSave, '-dpsc2',resolution);
        end
    case 'distance_visualize_all' % create and save map for multiple projects/conditions with one go
        pj = varargin{1};
        conds = {{''},...
            {'_activesingle','_activemulti','_passivesingle'},... % ,'_passivemulti' 
            {'_unimanual','_bimanual'}};
        % pasivemulti is useless as it has almost no discriminability        
        
        for p=pj;
            for con=1:numel(conds{p})
                message=sprintf('Now printing: %s %s...\n', datasets{p},conds{p}{con});
                figure('unit','centimeters','position',[25,25,30,10]);
                title({'',message,''}); axis off; set(gca, 'fontsize', 20); drawnow;
                
                clusterRDM_imana('distance_visualize', p,...
                    'condition', conds{p}{con},...
                    'threshold', 0, ...
                    'smooth', 25,...
                    'isSave', 1,...
                    'resolution','-r200');
            end
            close all;
        end
        
    case 'cluster_estimate' % Run spectral clustering algorithm on sparse searchlight
        pj = varargin{1}; % project index
        % options
        SLname  = 'clusterRDM_masked_roi_7mmSeparation_150vox'; %'masked_roi_250vox'; % 90
        alpha = 0.5; % proportion of subjects
        killnegdist=0; % if we force negative distance to be zero
        estWindivid=1; % if we individually calculate W and then average
        disttype = 'correlation';%'euclidean'; % distance metric to calculate within-subject RDM-similarity matrix
        maxclust = 10; % 10
        usemeandist=0;
        cutoff = 0.5;
        normalization = 'JordanWeiss'; % 'ShiMalik', 'non'
        condition = ''; % active or passive (sh2), unimanual or bimanual (bmw1)
        similarity = 'average'; %average'' 'median' 'tvalue'
        isSave = 1;
        
        % coloring options
        colorspace = 'rgb';
        cweight = 0.5;
        order = [1,2,3];
        
        vararginoptions(varargin(2:end),{'SLname','alpha','killnegdist',...
            'estWindivid','disttype','maxclust','usemeandist','normalization','condition',...
            'colorspace','cweight','order','similarity'});
        
        % loop over projects
        for p=pj
            project = datasets{p};
            sn = subjects{p};
            subj_name = subj_names{p};
            % Get project-wise directory structure
            str=sprintf('Dir=%s;', datasets{p});
            eval(str);
            caretDir = Dir.caretDir;
            caretGroupDir = fullfile(caretDir, 'fsaverage_sym');
            
            % Load distance data
            distfname = fullfile(analyzeDir, sprintf('%s_distance_%s%s.mat', project, SLname,condition));
            T = load(distfname);
            
            % choose good nodes
            G = pivottablerow(T.roi, ~all(isnan(T.data),2), 'nansum');
            goodrois = find(G>alpha*length(sn)); % [1:length(G)];
            fprintf('\n%d/%d nodes has more than %d percent of subject.\n', length(goodrois), length(G), alpha*100);
            Nroi = length(unique(T.roi));
            
            % Evaluate mean-distance and apply mask?
%             switch (condition)
%                 case '_activesingle'
%                     keyboard(); % Not well adjusted?
%                     i=1;
%                 case '_activemulti'
%                     i=2;
%                 case '_passivesingle'
%                     i=1;
%                 case '_passivemulti'
%                     i=2;
%                 case '_unimanual'
%                     i=1;
%                 case '_bimanual'
%                     i=3;
%                 otherwise
%                     i=1;
%             end
%             threshold = [maskthresholds{p}(i)]*0.4;
%             meanD=pivottablerow(T.roi, nanmean(T.data,2), 'nanmean');
            
            
            % Do clustering
            tic; W = [];
            for s=sn % loop over subjects
                % calc adjacency matrix
                D = getrow(T, T.subj==s);
                D = tapply(D, {'nodeID', 'roi', 'roiname', 'desikanid','hemis','subj'},...
                    {'data','nanmean(x,1)','name','data'},...
                    {'coord','nanmean(x,1)','name','coord'});
                gidx = ismember(D.roi, goodrois);
                P = getrow(D, gidx);
                w = clusterRDM_calcAdjacencyMatrix(P.data, 'disttype', disttype,'cutoff', cutoff);
                W = cat(3,W, full(w));
            end
            % average W
            switch (similarity)
                case 'average' % simple average
                    Wavrg = nanmean(W,3);                    
                case 'median' % median
                    Wavrg = nanmedian(W,3);                    
                case 'tvalue' % take inter-subject variability into account
                    Wavrg = nanmean(W,3);
                    Wavrg = Wavrg./(nanstd(W,[],3)/sqrt(size(W,3)));
                    Wavrg = Wavrg/max(Wavrg(isfinite(Wavrg)));                    
            end
            Wavrg(isnan(Wavrg)) = eps; % prevent NaN
            Wavrg(logical(eye(size(Wavrg,1)))) = eps; % kill similarity to oneself
            
            % feed to spectral clustering
            C = zeros(Nroi,1); L = NaN(Nroi); U = NaN(Nroi, maxclust);
            [C_, L_, U_] = clusterRDM_spectralClustering(Wavrg, maxclust, 'normalization', normalization); toc;
            
            % put back to original roi size
            C(gidx,1) = C_;
            L(gidx, gidx) = L_;
            U(gidx, :) = U_;
            
            % -------------------------------------------------------------
            % do hierarchical clustering of result and create color
            Nclusters = numel(unique(C_));
            
            % Evaluate similarity graph (w)
            [Csort, idx] = sort(C_); % sort wrt cluster
            Wsort = Wavrg(idx,:); Wsort=Wsort(:,idx);
            
            % Laplacian eigen vector
            Lsort = L_(idx,:); Lsort=Lsort(:,idx);
            Usort = U_(idx,:);
            
            % Merge similarities according to clusters
            Wsort(logical(eye(size(Wsort)))) = NaN;
            for i=1:Nclusters
                for j=1:Nclusters
                    idxi=Csort==i;
                    idxj=Csort==j;
                    rWsort(i,j) = nanmean(vec(Wsort(idxi,idxj)));
                    rLsort(i,j) = nanmean(vec(Lsort(idxi,idxj)));
                end
                rUsort(i,:) = nanmean(Usort(idxi,:),1); % cluster-mean of Laplacian eigen vectors
            end
            
            % Get dendrogram linkage
            Z = linkage(real(rUsort),'ward','euclidean'); % use cluster-mean of Laplacian eigen vectors, Uc
            % because the Laplacian eigen vectors (U) are responsible for spectral clustering result
            
            % Get colors
            color = colorDendrogram(Z,size(rUsort,1),'colorspace', colorspace,'order',order, 'fig',0,'weight',cweight);
            
            % show diagnostic plot
            hdiag = myFigure([25,25]); colormap('jet');
            [Csort,idx]=sort(C_);
            [V,lambda] = eig(L_);
            subplot(2,2,3); imagesc(Wavrg, [0,1]); axis square
            subplot(2,2,4); imagesc(Wavrg(idx,idx), [0,1]); axis square
            subplot(2,2,[1,2]); bar(diag(lambda)); set(gca,'xlim', [1,2*maxclust],'ylim',[min(diag(lambda(2:end,2:end))), max(lambda(:))]);
            ylabel('Laplacian eigen values'); xlabel('# eigen values');
            drawline(maxclust,'dir','vert'); title(sprintf('%s (%d clusters)', project, maxclust));
            
            % save result
            param.killnegdist=killnegdist;
            param.estWindivid=estWindivid;
            param.disttype=disttype;
            param.maxclust=maxclust;
            param.usemeandist=usemeandist;
            param.noamalization = normalization;
            param.cutoff = cutoff;
            
            cluster.idx = gidx;
            %cluster.C=C_;
            cluster.C=C;
            cluster.U=U;
            cluster.L=L;
            cluster.nodeID = D.nodeID;
            cluster.roi = D.roi;
            cluster.hemis = D.hemis;
            cluster.roiname = D.roiname;
            
            W = W;
            Z = Z;
            color = color;
            
            % Save result
            switch lower(normalization)
                case 'shimalik'
                    postfix = 'ShiMalik';
                case 'non'
                    postfix = 'Unnormalized';
                otherwise
                    postfix = '';
            end
            fname = sprintf('%s%s_cluster_%s%d%s.mat', project,condition,disttype,maxclust,postfix);
            if isSave; save(fullfile(analyzeDir,fname),'param','cluster','W','Z','color');  end; % result
            
            figProjDir = fullfile(figDir, project); dircheck(figProjDir);
            figname = sprintf('%s%s_cluster_%s%d%s', project,condition,disttype,maxclust,postfix);
            mySaveFig(hdiag, fullfile(figProjDir, figname), isSave, '-dpng', '-r150');
        end
        varargout = {cluster};
    case 'cluster_estimate_all' % create and save map for multiple projects/conditions with one go
        pj = [3];
        conds = {{''},{'_activesingle','_activemulti','_passivesingle'},... % ,'_passivemulti'
            {'_unimanual','_bimanual'}};
        clusters = [4:20];
        
        for p=pj;
            for c=clusters
                for con=1:numel(conds{p})
                    message=sprintf('Now estimating: %s %s Ncluster=%d...\n', datasets{p},conds{p}{con},c);
                    figure('unit','centimeters','position',[25,25,30,10]);
                    title({'',message,''}); axis off; set(gca, 'fontsize', 20); drawnow;
                    
                    clusterRDM_imana('cluster_estimate', p,...
                        'condition', conds{p}{con},...
                        'maxclust', c);
                end
                close all;
            end
        end
    case 'cluster_visualize' % visualize clustering result on surface
        pj = varargin{1}; % project index
        SLname  = 'clusterRDM_masked_roi_7mmSeparation_150vox'; %'masked_roi_250vox'; % 90
        disttype = 'correlation';%'euclidean'; % distance metric to calculate within-subject RDM-similarity matrix
        maxclust = 10; % 10
        overlap = 0;
        bgcolor = 'w';
        condition = '';
        isSave = 0;
        resolution = '-r400';
        normalization = 'JordanWeiss';
        vararginoptions(varargin(2:end),{'SLname','alpha',...
            'killnegdist','condition','isSave',...
            'estWindivid','disttype','maxclust',...
            'usemeandist','normalization','overlap',...
            'bgcolor','resolution','normalization'});
        switch lower(normalization)
            case 'shimalik'
                postfix = 'ShiMalik';
            case 'non'
                postfix = 'Unnormalized';
            otherwise
                postfix = '';
        end
        
        % Define plot range
        str=sprintf('Dir=%s;', datasets{1});
        eval(str);        
        % Set mask (background region on the basic flat map)
        for pp=1:length(datasets)
            % Load distance mask for background
            %mask = load(fullfile(analyzeDir, 'joint_distmask.mat'));
            mask = load(fullfile(analyzeDir, sprintf('%s_distmask.mat', datasets{pp})));
            fn = fieldnames(mask);
            bgmask(:,:,pp) = mask.(fn{1});
        end;
        masks=bgmask;
        bgmask=sum(bgmask,3)>0;
        for h=1:2
            [Flat] = getSurface(fullfile(Dir.caretDir,'fsaverage_sym'), h);
            x(1,h) = min(Flat.data(bgmask(:,h),1));
            x(2,h) = max(Flat.data(bgmask(:,h),1));
            y(1,h) = min(Flat.data(bgmask(:,h),2));
            y(2,h) = max(Flat.data(bgmask(:,h),2));
        end;
        xrange = [floor(min([x(1,1),-x(2,2)])), ceil(max([x(2,1), -x(2,1)]))];
        yrange = [floor(min([y(1,1),y(1,2)])), ceil(max([y(2,1), y(2,2)]))];
        plotrange = {{xrange,yrange},{-xrange([2,1]),yrange}};
            
        for p=pj;
            project = datasets{p};
            sn = subjects{p};
            subj_name = subj_names{p};
            Nsubj = numel(sn);
            
            % Get project-wise directory structure
            str=sprintf('Dir=%s;', datasets{p});
            eval(str);
            caretDir = Dir.caretDir;
            glmDir  = Dir.glmDir;
            caretGroupDir = fullfile(caretDir, 'fsaverage_sym');
            
            % Load cluster result
            fname = sprintf('%s%s_cluster_%s%d%s.mat', project,condition,disttype,maxclust,postfix);
            %C = load(fullfile(analyzeDir, fname));
            load(fullfile(analyzeDir, fname)); % param, cluster, W, Z, color
            
            % Load node assignment
            overlapstr = {'unoverlapped','overlappled'};
            fname = fullfile(analyzeDir, sprintf('%s_%s_sharednodes_%s.mat',project, SLname,overlapstr{overlap+1}));
            load(fname);
                        
            % Get color
            %[colors, G, Uc] = sh1_pcm('clusterRDM_defineColor',nanmean(C.W,3),C,'cweight',0.5+0.02*maxclust);
            colors=color;
            
            % Loop over clusters and overlay
            C = cluster;
            C.data = indicatorMatrix('identity_p', C.C);
            fighandle = []; window = [];
            for c=1:maxclust
                T = C;
                T.data = T.data(:,c);
                MAP = colors(c,:);
                modelname = '';
                maptex = sprintf('C%d',c);
                
                [maskout, fighandle, scalehandle, window]=clusterRDM_imana('vis_map', T, SLnodes, caretGroupDir,...
                    'plotUnderlay', c==1,...
                    'bgmask', masks(:,:,p),...
                    'makeFig', c==1,...
                    'plotrange', plotrange,...
                    'MAP', MAP,...
                    'modelname', modelname,...
                    'printScale',0,...
                    'threshold',0.8,...
                    'figurehandle', fighandle,...
                    'figwindow', window,...
                    'bgcolor', bgcolor,...
                    'printScale',0);
            end
            clusterstr = sprintf('%s%dclust%s',disttype,maxclust,postfix);
            
            % Save map&scale
            figProjDir = fullfile(figDir, project); dircheck(figProjDir);
            figname = fullfile(figProjDir, sprintf('%s%s_ClusterMap_%s_%s',project, condition, SLname, clusterstr));
            mySaveFig(fighandle, figname, isSave, '-dpng',resolution);
        end
    case 'cluster_visualize_all' % create and save map for multiple projects/conditions with one go
        pj = varargin{1};
        conds = {{''},...
            {'_activesingle','_activemulti','_passivesingle'},... % ,'_passivemulti' 
            {'_unimanual','_bimanual'}};
        % pasivemulti is useless as it has almost no discriminability
        clusters = [4,5:5:20];
        
        for c=clusters
            for p=pj;
                for con=1:numel(conds{p})
                    message=sprintf('Now printing: %s %s Ncluster=%d...\n', datasets{p},conds{p}{con},c);
                    figure('unit','centimeters','position',[25,25,30,10]);
                    title({'',message,''}); axis off; set(gca, 'fontsize', 20); drawnow;
                    
                    clusterRDM_imana('cluster_visualize', p,...
                        'condition', conds{p}{con},...
                        'maxclust', c, ...
                        'isSave', 1,...
                        'resolution','-r200');
                end
                close all;
            end
        end        
    case 'cluster_showRDM' % show cluster-wise RDMs
        % UNDER CONSTRUCTION %
    case 'cluster_dendrogram' % apply agglomerative clustering to the clustering result
        pj = varargin{1};
        maxclust = varargin{2};
        condition = ''; % additional condition (e.g., 'activesingle' for sh2, 'unimanual', for bmw1)
        disttype = 'correlation';
        isSave = 0;
        vararginoptions(varargin(3:end),{'condition','distatype','isSave'});
        
        % loop over projects
        for p=pj
            project = datasets{p};
            sn = subjects{p};
            subj_name = subj_names{p};
            distanceform = distanceforms{p};
            glmtype = glmtypes(p);
            % Get project-wise directory structure
            str=sprintf('Dir=%s;', datasets{p});
            eval(str);
            caretDir = Dir.caretDir;
            caretGroupDir = fullfile(caretDir, 'fsaverage_sym');
            glmDir  = Dir.glmDir;
            
            % load cluster result (Z)
            fname = sprintf('%s%s_cluster_%s%d.mat', project,condition,disttype,maxclust);
            load(fullfile(analyzeDir, fname));
            
            %fig=myFigure([35,20]*0.5); % 0.7
            fig=myFigure([20,15+5*maxclust]*0.5); % 0.7
            subplot(2,1,1);
            [h,t,per]=dendrogram(Z,'Orientation','right');
            ylabel('Clusters'); xlabel('Dissimilarity');title('Dissimilarity');
            hold on; axis off;
            xlim=get(gca,'xlim');range=diff(xlim);
            for c=1:numel(per)
                plot(max([0, xlim(1)-0.2*range]),c, 'o','markersize',15,...
                    'markeredgecolor','k',...
                    'markerfacecolor',color(per(c),:)); hold on;
                %text(c,ylim(1)-0.3*yrange,c,maptex{per(c)},'horizontalalignment','center','fontsize',15);
            end;
            for i=1:numel(h)
                set(h(i),'color',[0,0,0],'linewidth',1.0);
            end
            set(gca,'xlim',[max([0, xlim(1)-0.2*range]),xlim(2)]);
            set(gca,'ytick',[1:10],'yticklabel',flip([1:10]));
            set(gca,'tickdir','out','ticklength',[0.02,0.02])
            
            subplot(2,1,2);
            for c=1:maxclust
                plot(1, c, 's','markersize',15,...
                    'markeredgecolor','k',...
                    'markerfacecolor',color(c,:)); hold on;
                text(2,c,sprintf('C%d',c),'horizontalalignment','left','fontsize',15);
            end; axis equal;axis off;
            set(gca,'ydir','reverse');
            title({'Original labeling', ''});
            
            figProjDir = fullfile(figDir, project); dircheck(figProjDir);
            figname = sprintf('%s%s_dendrogram_%s%d', project,condition,disttype,maxclust);
            ffigname = fullfile(figProjDir, figname);
            mySaveFig(fig, ffigname, isSave, '-dpng','-r150');
            mySaveFig(fig, ffigname, isSave, '-dpsc2','-r400');
            varargout = {per,color};
        end
    case 'cluster_dendrogram_all'
        pj = varargin{1};
        conds = {{''},...
            {'_activesingle','_activemulti','_passivesingle'},... % ,'_passivemulti'
            {'_unimanual','_bimanual'}};
        clusters = [4:20];
        
        for c=clusters
            for p=pj;
                for con=1:numel(conds{p})
                    message=sprintf('Now printing: %s %s Ncluster=%d...\n', datasets{p},conds{p}{con},c);
                    figure('unit','centimeters','position',[25,25,30,10]);
                    title({'',message,''}); axis off; set(gca, 'fontsize', 20); drawnow;
                    
                    clusterRDM_imana('cluster_dendrogram', p,c,...
                        'condition', conds{p}{con},...
                        'isSave', 1);
                end
                close all;
            end
        end
    case 'cluster_laplacian' % apply agglomerative clustering to the clustering result
        pj = varargin{1};
        maxclust = 5;
        condition = ''; % additional condition (e.g., 'activesingle' for sh2, 'unimanual', for bmw1)
        disttype = 'correlation';
        isSave = 0;
        vararginoptions(varargin(2:end),{'condition','distatype','isSave'});
        
        % loop over projects
        for p=pj
            project = datasets{p};
            % Get project-wise directory structure
            str=sprintf('Dir=%s;', datasets{p});
            eval(str);
            caretDir = Dir.caretDir;
            
            % load cluster result (L)
            fname = sprintf('%s%s_cluster_%s%d.mat', project,condition,disttype,maxclust);
            load(fullfile(analyzeDir, fname));
            
            % calc eigenvalue
            [v, lam] = eigs(cluster.L(cluster.idx,cluster.idx), 25, eps);
            
            fig=myFigure([20,20/sqrt(2)]*0.5);
            
            bar(diag(lam), 0.8, 'facecolor', [0.8,0.8,0.8]); grid on;
            set(gca, 'ylim', [0.8, 1], 'ytick', [0.5, 0.8, 0.9, 0.95, 1]);
            set(gca, 'xlim', [0, 26], 'xtick', [1,5,10,15,20,25]);
            title([project,condition],'interpreter','none'); xlabel('#Clusters'); ylabel('Eigenvalue');
            
            figProjDir = fullfile(figDir, project); dircheck(figProjDir);
            figname = sprintf('%s%s_LaplacianEigenValues_%s', project,condition,disttype);
            ffigname = fullfile(figProjDir, figname);
            mySaveFig(fig, ffigname, isSave, '-dpng','-r150');
            mySaveFig(fig, ffigname, isSave, '-dpsc2','-r350');
        end
    case 'cluster_laplacian_all'
        pj = varargin{1};
        conds = {{''},...
            {'_activesingle','_activemulti','_passivesingle'},... % ,'_passivemulti'
            {'_unimanual','_bimanual'}};
        for p=pj;
            for con=1:numel(conds{p})
                message=sprintf('Now printing: %s %s...\n', datasets{p},conds{p}{con});
                figure('unit','centimeters','position',[25,25,30,10]);
                title({'',message,''}); axis off; set(gca, 'fontsize', 20); drawnow;
                
                clusterRDM_imana('cluster_laplacian', p,...
                    'condition', conds{p}{con},...
                    'isSave', 1);
            end
            close all;
        end
    
    case 'cluster_stability_subj' % Assess stabilities of clustering result using bootstrap/subsampling of participants
        pj = varargin{1}; % project index
        % options
        SLname  = 'clusterRDM_masked_roi_7mmSeparation_150vox'; %'masked_roi_250vox'; % 90
        alpha = 0.5; % proportion of subjects
        killnegdist=0; % if we force negative distance to be zero
        estWindivid=1; % if we individually calculate W and then average
        disttype = 'correlation';%'euclidean'; % distance metric to calculate within-subject RDM-similarity matrix
        maxclust = 10; % 10
        usemeandist=0;
        cutoff = 0.05;
        normalization = 'JordanWeiss'; % 'ShiMalik', 'non'
        condition = ''; % active or passive (sh2), unimanual or bimanual (bmw1)
        isSave = 1;
        Niter = 1000;
        similarity = 'average';
        
        % coloring options
        colorspace = 'rgb';
        cweight = 0.5;
        order = [1,2,3];
        
        vararginoptions(varargin(2:end),{'SLname','alpha','killnegdist',...
            'estWindivid','disttype','maxclust','usemeandist','normalization','condition',...
            'colorspace','cweight','order','Niter','average'});
        
        % loop over projects
        for p=pj
            project = datasets{p};
            sn = subjects{p};
            subj_name = subj_names{p};
            % Get project-wise directory structure
            str=sprintf('Dir=%s;', datasets{p});
            eval(str);
            caretDir = Dir.caretDir;
            caretGroupDir = fullfile(caretDir, 'fsaverage_sym');
            
            % Load distance data
            distfname = fullfile(analyzeDir, sprintf('%s_distance_%s%s.mat', project, SLname,condition));
            T = load(distfname);
            
            % choose good nodes
            G = pivottablerow(T.roi, ~all(isnan(T.data),2), 'nansum');
            goodrois = find(G>alpha*length(sn)); % [1:length(G)];
            fprintf('\n%d/%d nodes has more than %d percent of subject.\n', length(goodrois), length(G), alpha*100);
            Nroi = length(unique(T.roi));
            
            % Evaluate mean-distance and apply mask?
%             switch (condition)
%                 case '_activesingle'
%                     keyboard(); % Not well adjusted?
%                     i=1;
%                 case '_activemulti'
%                     i=2;
%                 case '_passivesingle'
%                     i=1;
%                 case '_passivemulti'
%                     i=2;
%                 case '_unimanual'
%                     i=1;
%                 case '_bimanual'
%                     i=3;
%                 otherwise
%                     i=1;
%             end
%             threshold = [maskthresholds{p}(i)]*0.4;
%             meanD=pivottablerow(T.roi, nanmean(T.data,2), 'nanmean');
            
            
            % Do clustering
            tic; W = [];
            for s=sn % loop over subjects
                % calc adjacency matrix
                D = getrow(T, T.subj==s);
                D = tapply(D, {'nodeID', 'roi', 'roiname', 'desikanid','hemis','subj'},...
                    {'data','nanmean(x,1)','name','data'},...
                    {'coord','nanmean(x,1)','name','coord'});
                gidx = ismember(D.roi, goodrois);
                P = getrow(D, gidx);
                w = clusterRDM_calcAdjacencyMatrix(P.data, 'disttype', disttype,'cutoff', cutoff);
                W = cat(3,W, full(w));
            end
            
            % Run resampling
            CV = []; N=round(numel(sn)*0.8);
            for c=[2:10, 12:2:20];
                csim = zeros(numel(sn),1);
                fprintf('#Cluster=%d, iteration=',c);
                Idx1 = ceil(rand(N,Niter)*numel(sn)); % otherwise random seed is reset in spectralclustering
                Idx2 = ceil(rand(N,Niter)*numel(sn)); 
                for s=1:Niter
                    switch (similarity)
                        case 'average' % simple average
                            W1 = nanmean(W(:,:,Idx1(:,s)),3); 
                            W2 = nanmean(W(:,:,Idx2(:,s)),3);
                            postfix = '';
                        case 'median' % median
                            W1 = nanmedian(W(:,:,Idx1(:,s)),3); 
                            W2 = nanmedian(W(:,:,Idx2(:,s)),3);
                            postfix = 'median';
                        case 'tvalue' % take inter-subject variability into account
                            W1 = nanmean(W(:,:,Idx1(:,s)),3); 
                            W1 = W1./(nanstd(W(:,:,Idx1(:,s)),[],3)/sqrt(size(W(:,:,Idx1(:,s)),3)));
                            W1 = W1/max(W1(isfinite(W1)));
                            W2 = nanmean(W(:,:,Idx2(:,s)),3); 
                            W2 = W2./(nanstd(W(:,:,Idx2(:,s)),[],3)/sqrt(size(W(:,:,Idx2(:,s)),3)));
                            W2 = W2/max(W2(isfinite(W2)));
                            postfix = 'tvalue';
                    end
                    W1(isnan(W1)) = eps;
                    W1(logical(eye(size(W1,1)))) = eps;
                    W2(isnan(W2)) = eps;
                    W2(logical(eye(size(W2,1)))) = eps;
                    
                    % feed to spectral clustering
                    C1 = zeros(Nroi,1); 
                    [C_] = clusterRDM_spectralClustering(W1, c, 'normalization', normalization,...
                        'maxiter', 50, 'replicates',20); 
                    C1(gidx,1) = C_;
                    
                    C2 = zeros(Nroi,1); 
                    [C_] = clusterRDM_spectralClustering(W2, c, 'normalization', normalization,...
                        'maxiter', 50, 'replicates',20); 
                    C2(gidx,1) = C_;
                    
                    % compare similarity between C1 and C2
                    warning off
                    csim(s) = clusterRDM_clusterSimilarity(C1, C2, 'similarity', 'jaccard');
                    
                    fprintf('%d',s);
                end
                cv.similarity = csim;
                cv.ncluster = repmat(c,size(csim));
                
                CV = addstruct(CV, cv);
                
                fprintf('\n');
            end
            
            fname = sprintf('%s%s_clusterstability_%s%d%s.mat', project,condition,disttype,Niter,postfix);
            if isSave; save(fullfile(analyzeDir,fname),'-struct','CV');  end; % result
            
        end
        varargout = {CV};
    case 'cluster_stability_analyze_subj' % Assess stabilities of clustering result using bootstrap/subsampling
        pj = varargin{1}; % project index
        % options
        SLname  = 'clusterRDM_masked_roi_7mmSeparation_150vox'; %'masked_roi_250vox'; % 90
        alpha = 0.5; % proportion of subjects
        killnegdist=0; % if we force negative distance to be zero
        estWindivid=1; % if we individually calculate W and then average
        disttype = 'correlation';%'euclidean'; % distance metric to calculate within-subject RDM-similarity matrix
        maxclust = 10; % 10
        usemeandist=0;
        cutoff = 0.05;
        normalization = 'JordanWeiss'; % 'ShiMalik', 'non'
        condition = ''; % active or passive (sh2), unimanual or bimanual (bmw1)
        isSave = 1;
        Niter = 1000;
        similarity = 'average'; 
        
        % coloring options
        colorspace = 'rgb';
        cweight = 0.5;
        order = [1,2,3];
        
        vararginoptions(varargin(2:end),{'SLname','alpha','killnegdist',...
            'estWindivid','disttype','maxclust','usemeandist','normalization','condition',...
            'colorspace','cweight','order','Niter','similarity'});
        
        % loop over projects
        for p=pj
            project = datasets{p};
            sn = subjects{p};
            subj_name = subj_names{p};
            % Get project-wise directory structure
            str=sprintf('Dir=%s;', datasets{p});
            eval(str);
            caretDir = Dir.caretDir;
            caretGroupDir = fullfile(caretDir, 'fsaverage_sym');
            
            % Load stablity result
            fname = sprintf('%s%s_clusterstability_%s%d%s.mat', project,condition,disttype,Niter,similarity);
            T = load(fullfile(analyzeDir, fname));
            
			% Sort and calc cumulative distribution
			S = []; edge = linspace(0,1,Niter/10); leg={};
			for c=unique(T.ncluster)'
				D = getrow(T, T.ncluster==c);
				s.count = histc(D.similarity,edge);
				s.cumu = cumsum(s.count)/Niter;
				s.edge = edge';
				s.ncluster = repmat(c, size(s.count));
				S = addstruct(S,s);
                leg{end+1} = sprintf('%d',c);
			end
			
			% Plot
			Ncluster = numel(unique(S.ncluster));
            color = cool(Ncluster);
			colors = mat2cell(color, ones(Ncluster,1),3);
			h = myFigure([45,15],'name',sprintf('%s%s', project, condition));
            
            subplot(1,2,1);
			lineplot(S.edge, S.cumu, 'split', S.ncluster, 'leg', leg, 'leglocation', 'northwest',...
				'markertype', 'non','linewidth', 3,'linecolor', colors,'errorcap',0.0001);
			xlabel('Similarity'); ylabel('Cumulative count');
            set(gca,'xtick',[0:0.1:1]);
            
            subplot(1,2,2);
			lineplot(S.ncluster, S.edge, 'subset', S.cumu>=0.89&S.cumu<=0.91,...
				'markertype', 'o','linewidth', 3,'linecolor', colors,'errorcap',0.0001,...
                'markercolor', colors, 'markerfill',[1,1,1],'markersize', 5);
			xlabel('Ncluster'); ylabel('Similarity at c=0.9');
            set(gca,'xtick',[2:max(S.ncluster)]);
            
        end
        varargout = {};
    case 'cluster_stability_node' % Assess stabilities of clustering result using bootstrap/subsampling of participants
        pj = varargin{1}; % project index
        % options
        SLname  = 'clusterRDM_masked_roi_7mmSeparation_150vox'; %'masked_roi_250vox'; % 90
        alpha = 0.5; % minimally necessary proportion of subjects
        killnegdist=0; % if we force negative distance to be zero
        estWindivid=1; % if we individually calculate W and then average
        disttype = 'correlation';%'euclidean'; % distance metric to calculate within-subject RDM-similarity matrix
        maxclust = 10; % 10
        usemeandist=0;
        cutoff = 0.05;
        normalization = 'JordanWeiss'; % 'ShiMalik', 'non'
        condition = ''; % active or passive (sh2), unimanual or bimanual (bmw1)
        isSave = 1;
        Niter = 1000;
		Pdropout = 0.1;
        similarity = 'average';
        
        % coloring options
        colorspace = 'rgb';
        cweight = 0.5;
        order = [1,2,3];
        
        vararginoptions(varargin(2:end),{'SLname','alpha','killnegdist',...
            'estWindivid','disttype','maxclust','usemeandist','normalization','condition',...
            'colorspace','cweight','order','Niter','average'});
        
        % loop over projects
        for p=pj
            project = datasets{p};
            sn = subjects{p};
            subj_name = subj_names{p};
            % Get project-wise directory structure
            str=sprintf('Dir=%s;', datasets{p});
            eval(str);
            caretDir = Dir.caretDir;
            caretGroupDir = fullfile(caretDir, 'fsaverage_sym');
            
            % Load distance data
            distfname = fullfile(analyzeDir, sprintf('%s_distance_%s%s.mat', project, SLname,condition));
            T = load(distfname);
            
            % choose good nodes
            G = pivottablerow(T.roi, ~all(isnan(T.data),2), 'nansum');
            goodrois = find(G>alpha*length(sn)); % [1:length(G)];
            fprintf('\n%d/%d nodes has more than %d percent of subject.\n', length(goodrois), length(G), alpha*100);
            Nroi = length(unique(T.roi));
            
            % Evaluate mean-distance and apply mask?
%             switch (condition)
%                 case '_activesingle'
%                     keyboard(); % Not well adjusted?
%                     i=1;
%                 case '_activemulti'
%                     i=2;
%                 case '_passivesingle'
%                     i=1;
%                 case '_passivemulti'
%                     i=2;
%                 case '_unimanual'
%                     i=1;
%                 case '_bimanual'
%                     i=3;
%                 otherwise
%                     i=1;
%             end
%             threshold = [maskthresholds{p}(i)]*0.4;
%             meanD=pivottablerow(T.roi, nanmean(T.data,2), 'nanmean');
            
            
            % Do clustering
            tic; W = [];
            for s=sn % loop over subjects
                % calc adjacency matrix
                D = getrow(T, T.subj==s);
                D = tapply(D, {'nodeID', 'roi', 'roiname', 'desikanid','hemis','subj'},...
                    {'data','nanmean(x,1)','name','data'},...
                    {'coord','nanmean(x,1)','name','coord'});
                gidx = ismember(D.roi, goodrois);
                P = getrow(D, gidx);
                w = clusterRDM_calcAdjacencyMatrix(P.data, 'disttype', disttype,'cutoff', cutoff);
                W = cat(3,W, full(w));
            end
            
            % Run resampling nodes
			CV = []; 
			Nnodes = size(W,1);
			N=round(Nnodes*Pdropout); % number of nodes to be dropped
            for c=[2:10, 12:2:20];
                csim = zeros(Niter,1);
                fprintf('#Cluster=%d, iteration=',c);
                Idx1 = rand(Nnodes,Niter); 
                Idx2 = rand(Nnodes,Niter);
                for s=1:Niter
					idx1 = Idx1(:,s)>Pdropout;
					idx2 = Idx2(:,s)>Pdropout;
                    switch (similarity)
                        case 'average' % simple average
                            W1 = nanmean(W(idx1,idx1,:),3); 
                            W2 = nanmean(W(idx2,idx2,:),3);
                            postfix = '';
                        case 'median' % median
                            W1 = nanmedian(W(idx1,idx1,:),3); 
                            W2 = nanmedian(W(idx1,idx1,:),3);
                            postfix = 'median';
                        case 'tvalue' % take inter-subject variability into account
                            W1 = nanmean(W(idx1,idx1,:),3); 
                            W1 = W1./(nanstd(W(idx1,id1,:),[],3)/sqrt(size(W(idx1,idx1,:),3)));
                            W1 = W1/max(W1(isfinite(W1)));
                            W2 = nanmean(W(idx2,idx2,:),3); 
                            W2 = W2./(nanstd(W(idx2,idx2,:),[],3)/sqrt(size(W(idx2,idx2,:),3)));
                            W2 = W2/max(W2(isfinite(W2)));
                            postfix = 'tvalue';
                    end
                    W1(isnan(W1)) = eps;
                    W1(logical(eye(size(W1,1)))) = eps;
                    W2(isnan(W2)) = eps;
                    W2(logical(eye(size(W2,1)))) = eps;
                    
                    % feed to spectral clustering
                    C1 = zeros(Nroi,1); C1_=zeros(Nnodes,1);
					try
                    [C_] = clusterRDM_spectralClustering(W1, c, 'normalization', normalization,...
                        'maxiter', 50, 'replicates',20); 
					catch
					end
                    C1_(idx1,1) = C_;
					C1(gidx,1) = C1_;
                    
                    C2 = zeros(Nroi,1); C2_=zeros(Nnodes,1);
					try
                    [C_] = clusterRDM_spectralClustering(W2, c, 'normalization', normalization,...
                        'maxiter', 50, 'replicates',20); 
					catch
					end
					C2_(idx2,:) = C_;
                    C2(gidx,1) = C2_;
                    
                    % compare similarity between C1 and C2
                    warning off
                    csim(s) = clusterRDM_clusterSimilarity(C1, C2, 'similarity', 'jaccard');
                    
					if mod(s,round(Niter/10))==0
						fprintf(' %d ',s);
					end
                end
                cv.similarity = csim;
                cv.ncluster = repmat(c,size(csim));
                
                CV = addstruct(CV, cv);
                
                fprintf('\n');
            end
            
            fname = sprintf('%s%s_clusterstability_node_%s%d%s%d.mat', project,condition,disttype,Niter,postfix,100*Pdropout);
            if isSave; save(fullfile(analyzeDir,fname),'-struct','CV');  end; % result
            
        end
        varargout = {CV};
    case 'cluster_stability_analyze_node' % Assess stabilities of clustering result using bootstrap/subsampling
        pj = varargin{1}; % project index
        % options
        SLname  = 'clusterRDM_masked_roi_7mmSeparation_150vox'; %'masked_roi_250vox'; % 90
        alpha = 0.5; % proportion of subjects
        killnegdist=0; % if we force negative distance to be zero
        estWindivid=1; % if we individually calculate W and then average
        disttype = 'correlation';%'euclidean'; % distance metric to calculate within-subject RDM-similarity matrix
        maxclust = 10; % 10
        usemeandist=0;
        cutoff = 0.05;
        normalization = 'JordanWeiss'; % 'ShiMalik', 'non'
        condition = ''; % active or passive (sh2), unimanual or bimanual (bmw1)
        isSave = 1;
        Niter = 1000;
        similarity = 'average'; 
        postfix='';
        Pdropout = 0.1;
        P = 0.8; 
        
        % coloring options
        colorspace = 'rgb';
        cweight = 0.5;
        order = [1,2,3];
        
        vararginoptions(varargin(2:end),{'SLname','alpha','killnegdist',...
            'estWindivid','disttype','maxclust','usemeandist','normalization','condition',...
            'colorspace','cweight','order','Niter','similarity','postfix','Pdropout','P'});
        
        % loop over projects
        for p=pj
            project = datasets{p};
            sn = subjects{p};
            subj_name = subj_names{p};
            % Get project-wise directory structure
            str=sprintf('Dir=%s;', datasets{p});
            eval(str);
            caretDir = Dir.caretDir;
            caretGroupDir = fullfile(caretDir, 'fsaverage_sym');
            
            % Load stablity result
            fname = sprintf('%s%s_clusterstability_node_%s%d%s%d.mat', project,condition,disttype,Niter,postfix,100*Pdropout);
            T = load(fullfile(analyzeDir, fname));
            
			% Sort and calc cumulative distribution
			S = []; edge = linspace(0,1,Niter/10); leg={};
			for c=unique(T.ncluster)'
				D = getrow(T, T.ncluster==c);
				s.count = histc(D.similarity,edge);
				s.cumu = cumsum(s.count)/Niter;
				s.edge = edge';
				s.ncluster = repmat(c, size(s.count));
				S = addstruct(S,s);
                leg{end+1} = sprintf('%d',c);
			end
			
			% Plot
			Ncluster = numel(unique(S.ncluster));
            color = cool(Ncluster);
			colors = mat2cell(color, ones(Ncluster,1),3);
			h = myFigure([45,15],'name',sprintf('%s%s', project, condition));
            
            subplot(1,3,1);
			lineplot(S.edge, S.cumu, 'split', S.ncluster, 'leg', leg, 'leglocation', 'northwest',...
				'markertype', 'non','linewidth', 3,'linecolor', colors,'errorcap',0.0001);
			xlabel('Similarity'); ylabel('Cumulative count');
            set(gca,'xtick',[0:0.1:1]);
            
            subplot(1,3,2);
            subset = false(size(S.edge));
            [~, idx] = findpeaks(-abs(S.cumu-P));
            subset(idx) = true;
            
			[x,y]=lineplot(S.ncluster, S.edge, 'subset', subset,...
				'markertype', 'o','linewidth', 3,'linecolor', colors,'errorcap',0.0001,...
                'markercolor', colors, 'markerfill',[1,1,1],'markersize', 5);
			xlabel('Ncluster'); ylabel(sprintf('Similarity at p=%1.1f',P));
            set(gca,'xtick',[2:max(S.ncluster)]);
            
            subplot(1,3,3);
            lineplot(x(1:end-1), -diff(y'),...
            	'markertype', 'o','linewidth', 3,'linecolor', colors,'errorcap',0.0001,...
                'markercolor', colors, 'markerfill',[1,1,1],'markersize', 5);
			xlabel('Ncluster'); ylabel(sprintf('\DeltaF(p=%1.1f)',P));
            set(gca,'xtick',[2:max(S.ncluster)]);    
            drawline(0,'dir','horz');
        end
        varargout = {};
    
        
    case 'pcm_extract_prewhitenedbeta' % Extract and compute noise-normalized beta at each SL
        % UNDER CONSTRUCTION %
    case 'pcm_noiseceiling' % Run noise-ceiling estimation on sparse searchlight
        % UNDER CONSTRUCTION %
    case 'pcm_bms' % calculate protected exceedance probability using spm_bms
        % UNDER CONSTRUCTION %
        
    case 'vis_assign_nodes' % Create SL-center-nodes to all SL-member-nodes mapping (necessary for visualization)
        pj=varargin{1};
        SLname  = 'clusterRDM_masked_roi_7mmSeparation_150vox'; %'masked_roi_250vox'; % 90
        overlap = 0;
        mask = '%s_distmask.mat';%'joint_distmask.mat';
        vararginoptions(varargin(2:end), {'SLname','overlap','map','mask'});
        
        for p=pj;
            project = datasets{p};
            sn = subjects{p};
            subj_name = subj_names{p};
            Nsubj = numel(sn);
            
            % Get project-wise directory structure
            str=sprintf('Dir=%s;', datasets{p});
            eval(str);
            caretDir = Dir.caretDir;
            glmDir  = Dir.glmDir;
            caretGroupDir = fullfile(caretDir, 'fsaverage_sym');
            
            % load mask if specified
            nodeidx=[];
            if ~isempty(mask)
                D=load(fullfile(analyzeDir, sprintf(mask,project)));
                field=fieldnames(D);
                nodeidx=D.(field{1});
            end
            
            SLnodes = [];
            % Find common nodes that are shared across all the subjects
            for h=1:2
                if isempty(nodeidx)
                    idx = [];
                else
                    idx = nodeidx(:,h);
                end
                % Load surface info
                flat = fullfile(caretGroupDir, hemName{h}, sprintf('%s.FLAT.coord', hem{h}));
                Flat = caret_load(flat);
                
                c=1; tic(); Allnodes = [];
                for s=sn
                    % Searchlight definition
                    caret_subjDIR   = fullfile(caretDir,[atlasA{1},subj_name{s},'_reduced'],hemispheres{h});
                    Searchlight        = fullfile(caret_subjDIR, [SLname,'.mat']);
                    SL  = load(Searchlight);
                    
                    % assume size of SL at each hemisphere is the same across the subjects
                    Snodes = sh1_assignNode(SL,Flat,'overlap',overlap,'nodeidx', idx);
                    for node=1:numel(SL.nodeID);
                        Allnodes{node}{c} = Snodes{node};
                        %scatter(Flat.data(Snodes{node},1), Flat.data(Snodes{node},2), 5, 'r'); hold on; axis equal
                    end
                    c=c+1;s
                end
                for node=1:numel(SL.nodeID); % stupid way to count shared nodes
                    tmpnodes = cat(1,Allnodes{node}{:});
                    [nodes] = unique(tmpnodes');
                    count = zeros(size(nodes));
                    for i=1:numel(nodes)
                        count(i)=sum(tmpnodes==nodes(i));
                    end
                    SLnodes{h}{node} = nodes(count==Nsubj);
                end
                toc();
            end
            overlapstr = {'unoverlapped','overlappled'};
            fname = fullfile(analyzeDir, sprintf('%s_%s_sharednodes_%s.mat',project, SLname,overlapstr{overlap+1}));
            save(fname, 'SLnodes');
        end
        varargout = {SLnodes};
    case 'vis_map'               % Read SL-center values, assign into surrounding nodes, and directly plot nodes
        T = varargin{1}; % T.data ([Nsubj x Nhemi x Nroi]-by-1 vector) is to be plotted on the surface
        SLnodes = varargin{2};
        caretGroupDir = varargin{3};
        isSave = 0;
        SLname  = 'clusterRDM_masked_roi_7mmSeparation_150vox'; %'masked_roi_250vox'; % 90
        hemis   = [1,2];
        threshold = 0; % threshold for T.data
        overlap = 0;
        mask = [];
        maskout = [];
        bgmask = [];
        
        plotUnderlay = 1;
        printScale = 0;
        valname = 'value';
        modelname = 'modelname';
        
        MAP = [];
        bgcolor = 'w';'k';
        lighting = 1;
        alpha = 1;
        bordersize = 10;
        smooth = 13;
        
        deleteMetric = 1;
        makeFig = 1;
        figurehandle = [];
        figScale = [];
        figwindow = [];
        condition = '';
        ceildata = 0;
        shownum = 0;
        resolution = '-r400';
        
        regionborder = [];
        
        plotrange = {{[-79.8252,104.7103],[-64.1677,81.5173-10]},...
            {[-104.7103,79.8252],[-64.1677,81.5173-10]}};
        
        vararginoptions(varargin(4:end),{'MAP','bgcolor','Nlevel',...
            'threshold','fontsize', 'printScale','valname','modelname',...
            'lighting','overlap','SLname','smooth','isSave','deleteMetric',...
            'alpha','makeFig','mask','postfix','ceildata','shownum',...
            'resolution','regionborder','plotUnderlay',...
            'figurehandle','figScale','figwindow','plotrange','bgmask'});
        home=pwd;
        if smooth==1
            prefix='s';
        else
            prefix='';
        end
        if isempty (MAP)
            MAP = parula(100);
        end
        Nlevel = size(MAP,1);
        
        % Make figure window
        aspectratio = diff(plotrange{1}{2}) / diff(plotrange{1}{1}) /2;
        width = 40;
        heightMap= width*aspectratio*1.05; % 20
        heightScale=20-heightMap;
        if makeFig==1
            figurehandle = myFigure([width heightMap]);
            figurehandle.Color = bgcolor;
            figurehandle.InvertHardcopy = 'off';
            if printScale==1
                figscale = myFigure([width heightScale]);
                figscale.Color = bgcolor;
                figscale.InvertHardcopy = 'off';
            end
        end
        topmargin = 0.01; %0.2
        colorbarsize = [0.2 0.2]; % width height
        centermargin = 0.005;
        axespositions = {[0.0 topmargin 0.5-centermargin 1.0-topmargin],...
            [0.5+centermargin topmargin 0.5-centermargin 1.0-topmargin], ...
            [0.5-colorbarsize(1)/2 0.5-colorbarsize(2)/2 colorbarsize(1), colorbarsize(2)]};
        % third one is for color scale
        
        Tab = []; mapdata=[]; figure(figurehandle);
        for h=hemis
            % Load all surface-related information
            [Flat{h}, Infl, Pial, White, Topo_cut, Topo_closed, Shape{h}, Border{h}, fname{h}] = getSurface(caretGroupDir, h);
            
            % Set axis
            if isempty(figwindow)
                ax(h)=axes('position',axespositions{h});
            else
                ax(h) = figwindow(h);
            end
            
            idx(:,h) = true(length(Flat{h}.index),1);
            
            % Display flatmap (underlaid image)
            xrange      = [floor(min(Flat{h}.data(:,1))),ceil(max(Flat{h}.data(:,1)))];
            yrange      = [floor(min(Flat{h}.data(:,2))),ceil(max(Flat{h}.data(:,2)))];
            depth = Shape{h}.data(:,2);
            if plotUnderlay
                caret_plotflatmap('coord',fname{h}.flat,'topo',fname{h}.cut,...
                    'data',depth,'cscale',[min(depth),max(depth)*lighting],...
                    'xlims',xrange,'ylims',yrange);colormap('gray'); hold on;
                if ~isempty(bgmask);
                    caret_plotflatmap_ay('coord',fname{h}.flat,'topo',fname{h}.cut,...
                        'data',double(idx(:,h)),'idx',bgmask(:,h),...
                        'xlims',xrange,'ylims',yrange,'cscale',[0, 100],...
                        'map', 'hot','Nlevel',100,'alpha',0.8); hold on;
                end
            end
            
            % Define plot region
            xlim{h} = [min(Flat{h}.data(idx(:,h),1)), max(Flat{h}.data(idx(:,h),1))];
            ylim{h} = [min(Flat{h}.data(idx(:,h),2)), max(Flat{h}.data(idx(:,h),2))];
            
            % get all roi from one hemisphere
            H = getrow(T, T.hemis==h);
            
            % Calc summary data
            meanDataAll = pivottablerow(T.nodeID, T.data, 'nanmean(x,1)'); %nanmean(H.data(:,c,:),3); % mean data for plot
            meanData    = pivottablerow(H.roi, H.data, 'nanmean(x,1)'); %nanmean(H.data(:,c,:),3); % mean data for plot
            % SLnodes is based-on roi-based ordering
            
            % Assign values to node-groups
            hNode = SLnodes{h};
            nodedata{h} = NaN(Flat{h}.num_nodes,numel(hNode));
            for i=1:numel(hNode)
                nodedata{h}(hNode{i}, i) = meanData(i);
            end
            nodedata{h} = nanmean(nodedata{h}, 2);
            
            % Export as .metric file and apply smoothing?
            if (smooth>0)
                M = Shape{h}; % meandist
                M.data = nodedata{h};
                M.data(isnan(M.data)) = 0;
                M.num_cols = 1;
                M.column_name = {sprintf('%s-%s',valname,modelname)};
                metricName = fullfile(caretGroupDir, hemName{h}, sprintf('%s-%s-%s.metric',SLname,valname,modelname));
                caret_savemetric(metricName, M);
                S = caret_smooth(metricName, 'coord', fname{h}.flat, 'topo', fname{h}.cut,...
                    'iterations', smooth, 'fwhm', smooth,'algorithm','AN');%25
                sM = caret_load(S{1});
                nodedata{h} = sM.data; % swap node data with smoothed data
                % delete .metric file
                if deleteMetric
                    delete(metricName); delete(S{1});
                end
                if ceildata
                    nodedata{h}=ceil(nodedata{h});
                end
            end
            % update nodes to be mapped
            if ~isempty(threshold)
                mapidx(:,h) = logical(idx(:,h))&(nodedata{h}>threshold(1));
            end
            mapdata = [mapdata; nodedata{h}];
            mapset(:,h) = idx(:,h);
        end
        
        for h=hemis
            axes(ax(h)); % make current
            groupDir = fullfile(caretGroupDir,hemName{h});
            cd(groupDir);
            
            % Show by using caret map function
            cmin = min(mapdata); %min(meanDataAll);
            cmax = max(mapdata); %max(meanDataAll);
            
            % Overlay result taking depth into account
            %depth = Shape.data(:,2);
            %depth = (depth-min(depth))/(max(depth)-min(depth));
            xrange = [min(Flat{h}.data(:,1)), max(Flat{h}.data(:,1))];
            yrange = [min(Flat{h}.data(:,2)), max(Flat{h}.data(:,2))];
            if isempty(threshold)
                scaling = [cmin cmax];
            elseif length(threshold)==1
                scaling = [threshold cmax];
            elseif length(threshold)==2
                scaling = [threshold(1),threshold(2)];
            elseif length(threshold)==3
                scaling = [threshold(1), threshold(3)];
                nodedata{h}(nodedata{h}<threshold(2))=threshold(1); % forece below-threshold ones to floor
            end
            if isempty(mask)
                %mapidx=idx(:,h);
            else
                mapidx(:,h)=mapidx(:,h)&logical(mask(:,h));
            end
            caret_plotflatmap_ay('coord',fname{h}.flat,'topo',fname{h}.cut,...
                'data',nodedata{h},'cscale',scaling,...
                'xlims',xrange,'ylims',yrange, 'idx', mapidx(:,h),'map', MAP,...
                'alpha',alpha,'Nlevel', Nlevel);
            maskout(:,h) = mapidx(:,h);
            
            title(''); % remove title
            %set(gca,'color',bgcolor,'view',[90+180*(1-h), 90]);
            set(gcf,'color',bgcolor);
            %set(gca, 'xlim', [min(Flat.data(:,1)), max(Flat.data(:,1))]);
            %set(gca, 'ylim', [min(Flat.data(:,2)), max(Flat.data(:,2))]);
            set(gca, 'xlim', xlim{h}, 'ylim', ylim{h});
            axis equal;
            
            % ============================= %
            % Draw border
            % ============================= %
            % Draw region border if specified
            if ~isempty(regionborder)
                caret_drawborder_flatmap('coord',fname{h}.flat,'topo',fname{h}.cut,...
                    'xlims',xrange,'ylims',yrange, 'idx', regionborder(:,h)&mapset(:,h),'alpha',alpha,'bordercolor',[1,1,1]);
            end
            
            % Plot border dots
            switch (h)
                case 1
                    borderidx   = cat(1,Border{h}.Border.vertex);
                    borderX{h}     = Flat{h}.data(borderidx(:,1),1);
                    borderY{h}     = Flat{h}.data(borderidx(:,1),2);
                case 2
                    borderidx   = cat(1,Border{h}.Border.vertex);
                    borderX{h}     = Flat{h}.data(borderidx(:,1),1);
                    borderY{h}     = Flat{h}.data(borderidx(:,1),2);
                    %borderX     = -borderX;
            end
            plot(borderX{h},borderY{h},'w.','markersize',bordersize);
            
            % Sow roi number if requested
            if shownum
                roinum=unique(H.roi);
                for roi=roinum'
                    nodeID = H.nodeID(H.roi==roi);
                    nodeID = nodeID(1);
                    x = Flat{h}.data(nodeID,1);
                    y = Flat{h}.data(nodeID,2);
                    text(x,y,sprintf('%d',roi),...
                        'horizontalalignment','center',...
                        'verticalalignment','middle',...
                        'color','w','fontsize',13);
                end
            end
            
            if ~isempty(plotrange{h}{1})
                xrange = plotrange{h}{1};
                yrange = plotrange{h}{2};
            end
            axis normal; axis off;
            axis equal;
            set(gca,'xlim',xrange,'ylim',yrange,'color',bgcolor);
            set(gcf,'color',bgcolor);
        end;
        
        % ============================================================= %
        % draw color scale
        % ============================================================= %
        if printScale
            figure(figscale);
            hlabel = axes('position',axespositions{3});
            
            % handle scaling
            xtickls=[scaling(1), scaling(2)];
            xticks=[1, Nlevel];
            % dealwith flooring threshold
            if length(threshold)==3
                edge=linspace(scaling(1),scaling(2),Nlevel);
                p1=discretize(threshold(2),edge);
                xtickls = sort([xtickls,threshold(2), threshold(2)+(threshold(3)-threshold(2))/2]);
                xticks = sort([xticks,p1,p1+ceil(Nlevel-p1)/2]);
            else
                xtickls=[scaling(1), scaling(1)+diff(scaling)/2, scaling(2)];
                xticks=[1,ceil(Nlevel/2),Nlevel];
            end
            xticklabel = cell(length(xticks),1);
            for i=1:length(xticks)
                xticklabel{i} = sprintf('%1.2f',xtickls(i));
            end
            
            switch bgcolor
                case 'w'
                    axcolor = [0 0 0];%[0.3,0.3,0.3];
                case 'k'
                    axcolor = [0.8,0.8,0.8];
            end
            
            % directly draw color bar
            for i=1:Nlevel
                patch([i-1,i-1,i,i],[0 1 1 0],MAP(i,:),'edgecolor','non');hold on;
            end;
            if exist('p1','var');
                drawline(p1,'dir','vert','color',[1,1,1]);
                for i=1:p1
                    patch([i-1,i-1,i,i],[0 1 1 0],MAP(1,:),'edgecolor','non');hold on;
                end
            end
            ylabel(modelname); %xlabel(valname);
            title(valname,'color',axcolor);
            set(gca,'xlim',[1,Nlevel],'xtick',xticks,...
                'xticklabel',xticklabel,...
                'ylim',[0,1],'ytick',[],'tickdir','out','ticklength',[0.05,0.1]);
            set(gca,'xcolor',axcolor,'ycolor',axcolor);
            
            %axis off;
            set(gcf,'color',bgcolor);
        end
        
        cd(home);
        if ~exist('figurehandle','var');
            figurehandle = [];
        end;
        if ~exist('figscale', 'var');
            figscale = [];
        end
        varargout = {maskout, figurehandle, figscale, ax};
        
    otherwise
        warning('No such case.');
end
end

%=== Local fun ======================== %
function dircheck(name)
if ~exist(name,'dir')
    mkdir(name);
end
end
function D = getDirectory(project, dataDir)
switch (project)
    case 'sh1'
        baseDir = fullfile(dataDir, 'SequenceLearning', 'sh1');
        D.glmDir = fullfile(baseDir, 'GLM_firstlevel_3');
    case 'sh2'
        baseDir = fullfile(dataDir, 'SequenceLearning', 'sh2');
        D.glmDir = fullfile(baseDir, 'GLM_firstlevel_fastnohpf_2exe');
    case 'bmw1'
        baseDir = fullfile(dataDir, 'BimanualWrist_MR', 'bmw1');
        D.glmDir = fullfile(baseDir, 'GLM_firstlevel_6direction_fastnohpf');
    otherwise
        warning('%s is not registered.', project);
end
D.freesurferDir   = fullfile(baseDir, 'surfaceFreesurfer');
D.caretDir        = fullfile(baseDir, 'surfaceCaret');
D.regDir          = fullfile(baseDir, 'RegionOfInterest');
%D.cerebDir        = fullfile(anatomicalDir,'SUIT');
%D.BGDir           = fullfile(anatomicalDir,'basal_ganglia');
end
