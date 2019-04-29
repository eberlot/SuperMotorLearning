function [canvas,colours,opt] = sml_style(styID)
%% Description
%   Default plotting parameters specified by the style styID
%   Define additional styles within the switch statement. Users can create
%   their own styles for their projects and use the style package to point
%   towards the user-defined style
%
%
% Author
%   Naveed Ejaz (ejaz.naveed@gmail.com)

canvas              = 'blackonwhite';
opt                 = [];
opt.save.journal    = 'brain';

% personalized colours
gray=[120 120 120]/255;
lightgray=[170 170 170]/255;
black=[0 0 0]/255;
blue=[49,130,189]/255;
mediumblue=[128,207,231]/255;
lightblue=[158,202,225]/255;
red=[222,45,38]/255;
mediumred=[237,95,76]/255;
%lightred=[252,146,114]/255;
lightred=[251,177,168]/255;
ms=7;

switch(styID)
    case 'default'
        colours                 = {'blue','green','red','orange','aqua','magenta','yellow','black'};
        opt.general.linestyle   = {'-','-','-','-','-','-','-','-',...
                                   '--','--','--','--','--','--','--','--'};
        canvas                  = 'blackonwhite';
        opt.save.journal        = 'brain';
        opt.legend.leglocation  = 'eastoutside';
    case 'gray'
        colours                 = {'black','lightgray','darkgray','black','lightgray','darkgray'};
        canvas                  = 'blackonwhite';
        opt.save.journal        = 'brain';
        opt.general.markertype  = {'o','v','^'};
        opt.general.linestyle   = {'-','-','-','--','--','--'};    
    case 'Seq'
        % set different shades for colours - colour-blind friendly
        colours                 = {red,blue};
        canvas                  = 'blackonwhite';
        opt.save.journal        = 'brain';
        opt.general.markertype  = {'o','v'};
        opt.general.markerfacecolor = colours;
        opt.general.linestyle   = {'-','--'};
        opt.general.markersize  = ms;
    case 'SeqShade'
        colours                 = {red,blue};
        canvas                  = 'blackonwhite';
        opt.save.journal        = 'brain';
        opt.general.markertype  = {'o','v'};
        opt.general.linestyle   = {'-','--'};
        opt.general.errorbars   = 'shade';
    case 'SeqSess'
        colours                 = {red,blue,mediumred,mediumblue,lightred,lightblue};
        canvas                  = 'blackonwhite';
        opt.save.journal        = 'brain';
        opt.general.markertype  = {'o','s'};
        opt.general.linestyle   = {'-','--'};
        opt.general.markersize  = ms;
    case 'SessSeq'
        colours                 = {red,mediumred,blue,mediumblue};
        canvas                  = 'blackonwhite';
        opt.save.journal        = 'brain';
        opt.general.markertype  = {'o','s'};
        opt.general.linestyle   = {'-','--'};
        opt.general.markersize  = ms;
    case 'Sess'
        % set different shades for colours - colour-blind friendly
        colours                 = {black,gray,lightgray,blue};
        canvas                  = 'blackonwhite';
        opt.save.journal        = 'brain';
        opt.general.markertype  = {'o','v','^','s'};
        opt.general.linestyle   = {'-'};
        opt.general.markersize  = ms;
    case 'Trained'
        colours                 = {red,mediumred,lightred};
        canvas                  = 'blackonwhite';
        opt.save.journal        = 'brain';
        opt.general.markertype  = {'o','^','s'};
        opt.general.linestyle   = {'-','--','-.'};
        opt.general.markersize  = ms;
    case 'Untrained'
        colours                 = {blue,mediumblue,lightblue};
        canvas                  = 'blackonwhite';
        opt.save.journal        = 'brain';
        opt.general.markertype  = {'o','^','s'};
        opt.general.linestyle   = {'-','--','-.'};
        opt.general.markersize  = ms;
    case 'SeqShade_small'
        colours                 = {red,blue};
        canvas                  = 'blackonwhite';
        opt.save.journal        = 'brain';
        opt.general.markertype  = {'.','.'};
        opt.general.linestyle   = {'-','-'};
        opt.general.errorbars   = 'shade';
    case 'SeqShade_light'
        colours                 = {lightred,lightblue};
        canvas                  = 'blackonwhite';
        opt.save.journal        = 'brain';
        opt.general.markertype  = {'.','.'};
        opt.general.linestyle   = {'-','-'};
        opt.general.errorbars   = 'shade';
    case 'SessShade_small'
        colours                 = {'darkgray','lightgray'};
        canvas                  = 'blackonwhite';
        opt.save.journal        = 'brain';
        opt.general.markertype  = {'.','.'};
        opt.general.linestyle   = {'-','-'};
        opt.general.errorbars   = 'shade';
    case 'SessTrained_shade'
        colours                 = {red,lightred};
        canvas                  = 'blackonwhite';
        opt.save.journal        = 'brain';
        opt.general.markertype  = {'.','.'};
        opt.general.linestyle   = {'-','-'};
        opt.general.errorbars   = 'shade';
    case 'SessUntrained_shade'
        colours                 = {blue,lightblue};
        canvas                  = 'blackonwhite';
        opt.save.journal        = 'brain';
        opt.general.markertype  = {'.','.'};
        opt.general.linestyle   = {'-','-'};
        opt.general.errorbars   = 'shade';    
    case 'Region'
        colours                 = {[84 13 100]/255,[238 66 102]/255,[14 173 105]/255,[59 206 172]/255,[255 210 63]/255,[78 164 220]/255,[10 10 180]/255,[170 170 170]/255};
        canvas                  = 'blackonwhite';
        opt.save.journal        = 'brain';
        opt.general.linestyle   = {'-'};
        opt.general.markertype  = {'o'};
        opt.legend.leglocation  = 'eastoutside';
       % opt.general.errorbars   = 'shade';
    case 'Region_shade'
        colours                 = {[84 13 100]/255,[238 66 102]/255,[14 173 105]/255,[59 206 172]/255,[255 210 63]/255,[78 164 220]/255,[10 10 180]/255,[170 170 170]/255};
        canvas                  = 'blackonwhite';
        opt.save.journal        = 'brain';
        opt.general.linestyle   = {'-'};
        opt.general.markertype  = {'o'};
        opt.legend.leglocation  = 'eastoutside';
        opt.general.errorbars   = 'shade';
    case 'Region_ceiling'
        colours                 = {[84 13 100]/255,[238 66 102]/255,[14 173 105]/255,[59 206 172]/255,[255 210 63]/255,[78 164 220]/255,[10 10 180]/255,[170 170 170]/255};
        canvas                  = 'blackonwhite';
        opt.save.journal        = 'brain';
        opt.general.linestyle   = {'--'};
        opt.general.markertype  = {'o'};
        opt.legend.leglocation  = 'eastoutside';
        opt.general.errorbars   = 'shade';
        opt.general.facealpha   = 0.1;
end;
