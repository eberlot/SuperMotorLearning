function [D]=sml_trial_old(MOV,D,fig,varargin)


sample		= 5;
forceThres	= 0.25;
minDist 	= 100/sample;

% ------------------------------------------------
% extract data
if (isempty(MOV))
    return;
end;

fingernames = {'Thumb','Index','Middle','Ring','Pinky'};
handnames = {'LH','RH'};

state		= MOV(:,1);
%TR          = MOV(:,2); need to add
TRtime		= MOV(:,2); % 3
t           = MOV(:,3); % 4
sampfreq	= 1000/sample;
% harvest force trace & smooth
F			= (MOV(:,4:end));
F			= smooth_kernel(F,4);
vF			= velocity_discr(F,2);

indx = [6 7 8 9 10; 1 2 3 4 5];  % row: left/right hand; column: 1-5 fing
% need to switch how hands are recorded


% find out the number of finger presses by counting "response" fields
A = fieldnames(D);
AllResp = [];
AllInstr = [];
presscnt = 0;   % press counter - responses
instrcnt = 0;   % press counter - instructions

% compare response, press - if incomplete
for i = 1:length(A)
    if length(A{i}) > 8;
        if strcmp(A{i}(1:8), 'response');
            presscnt = presscnt + 1;
            eval(['AllResp = [AllResp D.',A{i} , '];']);
        end
    elseif length(A{i}) > 5;
        if strcmp(A{i}(1:5), 'press');
            if instrcnt == 9
            else
                instrcnt = instrcnt + 1;
                eval(['AllInstr = [AllInstr D.',A{i} , '];']); 
            end
        end
    end
    
end

press_tot = length(find(AllResp>0));    % number of fingers pressed 
instr_tot = length(find(AllInstr>0));   % number of finger presses instructed 

D.numPress=press_tot;

if (press_tot == 0) || (press_tot ~= instr_tot) || (MOV(end,1) == 4 || MOV(end,1) ==5)   % if forces not recorded until end of trial - state 4 or 5 (instead of 7/8)
   
    D.peakForceTime = NaN(1,9);
    D.peakForceAmount = NaN(1,9);
    D.meanForce=NaN;
else
    
    for press = 0 : (press_tot - 1)
        
        pressName   = ['D.pressTime' , num2str(press) ];
        releaseName = ['D.releaseTime' , num2str(press) ];
        FingerName  = ['D.press' , num2str(press) ];
        
        PressIdx(:,press+1)   = .5* eval(pressName);
        ReleaseIdx(:,press+1) = .5* eval(releaseName);
        
    end
    
    for press = 1 : press_tot
        fing_indx =  AllResp(press);
        hand_indx =  D.hand;
        fing = indx(hand_indx,fing_indx);
        
        % Find the peaks in the force trace and record average force and pacing
        [i j] = max(F(PressIdx(press):ReleaseIdx(press),fing)); 
        % peak force - time
        k = find(F(:,fing)==i);
        
        if isempty(k)
            k=0;
        end
        
        peakHeight(press) = i;
        peakTime(press) = k*2;
        
    end
    
    if press_tot == 3    % for chunks - so it is the same for seq and chunks
        peakHeight(4:9)=0;
        peakTime(4:9)=0;
    end
    
    % extract peak forces time and amount
    D.peakForceTime = peakTime;
    D.peakForceAmount = peakHeight;    
    D.meanForce=mean(peakHeight);

    
    % ------------------------------------------------
    
    % Display trial
    if (fig>0)
        plot(t,F);
        hold on;
        plot(peakTime,peakHeight,'k*');
        %plot(t(peakTime(1:press_tot)),peakHeight(1:press_tot),'k*');
        
        title(['Force for presses of' handnames(D.hand)]);
        hold off;
        keyboard;
    end
end;
