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
lightred=[252,146,114]/255;
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
end;
