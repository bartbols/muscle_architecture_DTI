function varargout = MyHist( data,binsize,varName,units,varargin)
%MYHIST Plot a histogram and adds mean/sd and median to the plot.
%
% Bart Bolsterlee, Neuroscience Research Australia (NeuRA)
% February 2017
%
% ----------------- INPUT -----------------
% data    - nx1 or 1xn vector with data values of which the distribution will be plotted
% binsize - binsize of the histogram
% varName - variable name (to label x-axis)
% units   - variable unit
%
% Optional parameters, provided as 'ParameterName',<value> pairs:
% Type    -  type of plot: 'bar' or 'line'. Default: 'bar'
% Color   -  color of bar/line
% Fraction - logical (0/1) indicating whether the absolute (0 / false) or
%            fraction of the total (1 / true) number of fibres are plotted.
% Text     - 'mean' or 'all'. Displays the mean ('mean') or mean, median and
%            total number of fibres ('all') in a text box
% TextPos - two-element vector with the location of the text box. Default:
%            [0.95 0.95]
% BarAlpha - alpha-value of the bars

p = inputParser;
% addRequired(p,'DTItracts',@isstruct)
addRequired(p,'data')
addRequired(p,'binsize')
addRequired(p,'varName')
addRequired(p,'units')
addParameter(p,'Type','bar',@(x) any(strcmp(x,{'bar','line'})))
addParameter(p,'Color','b')
addParameter(p,'EdgeColor','none')
addParameter(p,'TextColor',[])
addParameter(p,'Fraction',true,@(x) x==0 || x==1 || islogical(x) )
addParameter(p,'Text','all',@(x) any(strcmp(x,{'mean','all'})));
addParameter(p,'TextPos',[0.95 0.95],@(x) isnumeric(x) && numel(x) == 2)
addParameter(p,'BarAlpha',1,@(x) isnumeric(x) && x >= 0 && x <= 1)
addParameter(p,'FontWeight','normal',@(x) strcmp(x,'normal') || strcmp(x,'bold'))
parse(p,data,binsize,varName,units,varargin{:})

color    = p.Results.Color;
plottype = p.Results.Type;
fraction = p.Results.Fraction;
TextPos  = p.Results.TextPos;
Text     = p.Results.Text;
BarAlpha = p.Results.BarAlpha;
fontweight = p.Results.FontWeight;
EdgeColor = p.Results.EdgeColor;
TextColor = p.Results.TextColor;
if isempty(TextColor)
    TextColor = color;
end

% Remove the Inf's
data(isinf(data)) = [];
first  = floor(min(data/binsize)) * binsize;
last   = ceil(max(data/binsize)) * binsize;
[n,c] = hist(data, first : binsize : last);
if numel(c) == 1
    c = c +[-1 0 1]*binsize;
    n = [0 n 0];
end
if nargout < 2
    % If no or one output is requested, plot the distribution as a bar or
    % line plot
    switch plottype
        case 'bar'
            if fraction == true
                handle = bar(c,n / sum(n),'FaceColor',color,...
                    'FaceAlpha',BarAlpha,...
                    'EdgeColor',EdgeColor);
            else
                handle = bar(c,n,'FaceColor',color,...
                    'FaceAlpha',BarAlpha,...
                    'EdgeColor',EdgeColor);
            end
        case 'line'
            if fraction == true
                handle = plot(c,n/sum(n),'Color',color,'LineWidth',2);
            else
                handle = plot(c,n,'Color',color,'LineWidth',2);
            end
    end
    if nargout == 1
        % Return the handle if one output argument is requested.
        varargout{1} = handle;
    end
    % Add title
    title(sprintf('%s, binsize = %.2f',varName,binsize),'FontSize',10)
    xlabel([varName ' (' units ')'])
    if fraction == true
        ylabel('fraction of total')
    else
        ylabel('number of fibres')
    end
    switch Text
        case 'mean'
            txt = sprintf('mean (sd) = %.2f (%.2f), median = %.2f',nanmean(data),nanstd(data),nanmedian(data));
        case 'all'
        txt = {sprintf('mean (sd) = %.2f (%.2f)',nanmean(data),nanstd(data)),...
               sprintf('median = %.2f%.2f)',nanmedian(data)),...
               sprintf('n = %d',sum(~isnan(data)))};
    end
    text(TextPos(1),TextPos(2),txt,...
        'HorizontalAlignment','Right','VerticalAlignment','Top',...
        'Color',TextColor,'Units','Normalized','FontWeight',fontweight,'FontSize',10)
elseif nargout == 2
    % If two outputs are requested, don't plot the data but return the bin
    % centres (c) and counts per bin (n) as outputs.
    varargout{1} = c;
    varargout{2} = n;
    
end

end % of function

