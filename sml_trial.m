function [D]=sml_trial(MOV,D,fig,varargin)

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

indx = [1:5; 6:10];  % row: left/right hand; column: 1-5 fing


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

%D.numPress=press_tot;

    
    for press = 0 : (press_tot - 1)
        
        pressName   = ['D.pressTime' , num2str(press) ];
        FingerName  = ['D.press' , num2str(press) ];        
        PressIdx(:,press+1)   = eval(pressName);
        
    end
    
    % ------------------------------------------------
    
    % Display trial
    if (fig>0)
        figure(1)
        subplot(2,1,1)
        
        if D.hand==1
            plot(t,F(:,1:5),'LineWidth',2);
            title('Force traces for presses of LEFT hand');
        else
            plot(t,F(:,6:10),'LineWidth',2);
            title('Force traces for presses of RIGHT hand');
        end
       
        hold on;
        drawline(PressIdx,'dir','vert');
        
        legend({'thumb','index','middle','ring','little'})
        hold off;
        
        subplot(2,1,2)
        plot(t,state,'LineWidth',2);
        ylim([1 8]);
        keyboard;
    end
end
