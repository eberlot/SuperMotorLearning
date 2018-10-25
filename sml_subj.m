function ANA=sml_subj(subjname,fig,block,trial)
if nargin<2
    fig=0;
end;
datafilename = ['data/sml1_' subjname '.dat'];
outfilename  = ['analyze/sml1_' subjname '.mat'];

ANA=[];
D=dload(datafilename);
if (nargin<3) 
    s=1;
else
    if (nargin<4)
        s=find(D.BN==block & D.TN==1); 
    else 
        s=find(D.BN==block & D.TN==trial); 
    end;
end;

%define  number of trials
trials=[s:length(D.BN)];
oldblock=-1;
for i=trials % loop over all trials
   if (oldblock~=D.BN(i))
        oldblock=D.BN(i);
        MOV=movload(['data/sml1_' subjname '_' num2str(D.BN(i),'%02d') '.mov']); % all trials of a block
   end;
   fprintf('%d %d\n',D.BN(i),D.TN(i));
   [C]=sml_trial(MOV{D.TN(i)},getrow(D,i),fig);
   C.BN = D.BN(i);
   % update TRA
   ANA=addstruct(ANA,C);
end;
D = ANA;
save(outfilename,'D');