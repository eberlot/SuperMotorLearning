function varargout = extract_ts(what,varargin)
baseDir         ='/Users/eberlot/Documents/Data/SuperMotorLearning';
%baseDir         ='/Volumes/MotorControl/data/SuperMotorLearning';
imagingDir      =[baseDir '/imaging_data'];                             
regDir          =[baseDir '/RegionOfInterest/']; 
glmSessDir      ={[baseDir '/glmSess/glmSess1'],[baseDir '/glmSess/glmSess2'],[baseDir '/glmSess/glmSess3'],[baseDir '/glmSess/glmSess4']}; % one glm per session
  
subj_name  = {'s01','s02','s03','s04','s05','s06','s07','s08','s09','s10',...
              's11','s12','s13','s14','s15','s16','s17','s18','s19','s20',...
              's21','s22','s23','s24','s25','s26','s27','s28','s29','s30','s31'}; 
sn = [5:9,11:31]; % good subjects          
switch what
case 'get_ts'
        sessN=1;
        sn = 11:20;
        reg = 1:3;
        vararginoptions(varargin,{'sessN','sn','reg'});
        OO = []; YY = [];
        for s=sn
            % load SPM and all regions
            load(fullfile(glmSessDir{sessN},subj_name{s},'SPM.mat'));
            SPM=spmj_move_rawdata(SPM,fullfile(imagingDir, subj_name{s}));
            load(fullfile(regDir,sprintf('%s_Brodmann_regions.mat',subj_name{s})));
            % extract data and onsets
            data = region_getdata(SPM.xY.VY,R(reg)); % here get the data
            O = spmj_get_ons_struct(SPM);     % Returns onsets in TRs, not secs
            O.sn = ones(size(O.ons))*s;
            % maybe consider filtering here to get y_adj
            for r=1:size(data,2)
                Y.y_raw     = data(r);
                Y.y_filt    = {spm_filter(SPM.xX.K,SPM.xX.W*Y.y_raw{:})};
                Y.B         = {SPM.xX.pKX*Y.y_filt{:}};
                Y.y_res     = {spm_sp('r',SPM.xX.xKXs,Y.y_filt{:})};
                Y.y_hat     = {SPM.xX.xKXs.X*Y.B{:}};
                Y.y_adj     = {Y.y_hat{:} + Y.y_res{:}};
                Y.sn        = s;
                Y.reg       = r;
                Y.xyz       = {R{r}.data};
                Y.flatcoord = {R{r}.flatcoord};
                YY = addstruct(YY,Y);
                fprintf('Done %s - reg %d\n',subj_name{s},r);
            end
            OO = addstruct(OO,O); 
            clear data;
        end
        % save new data structures
        save('task_structure.mat','-struct','OO');
        save('timeseries.mat','-struct','YY')
end