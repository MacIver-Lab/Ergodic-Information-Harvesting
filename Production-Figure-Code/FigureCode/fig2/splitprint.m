function [varargout] = splitprint(varargin)

% SPLITPRINT Divides components into different figures for printing
% SPLITPRINT will divide the current figure into two figures, one 
%   containing lines, text objects, and the containing axes.  It will 
%   be exported using the -depsc2 & -adobecset option, to 'splitprint.eps'.
%   -adobecset prevents certain errors, such as dropped minus signs
%   when the resulting eps file is imported into Adobe Illustrator.
%   The other will contain the patch and surface objects, with their axes' 
%   'Visible' property set to 'off'.  It will be exported using the 
%   -dtiff -r300 option to 'splitprint.tif'.
% SPLITPRINT(FIG) will divide the figure FIG into the default grouping
%   explained above.
% SPLITPRINT(FILENAME) will use the methods described above, but will 
%   print to files with names beginning with the string specified in
%   FILENAME.  FILENAME could also be a cell array of strings.
% SPLITPRINT(FIG,FILENAME) will divide the figure FIG and print to files
%   FILENAME.eps and FILENAME.tif.  If FILENAME is empty, then 
%   the figure will only be split, not printed to a file.
% SPLITPRINT(FIG,FILENAME,TYPES) will divide the objects into the groups
%   specified in TYPES.
% SPLITPRINT(FIG,FILENAME,TYPES,FORMATS) will print to the corresponding 
%   formats specified in the FORMATS cell array.
% SPLITPRINT(FIG,FILENAME,TYPES,FORMATS,ALIGN) wll create a registration
%   mark in the upper right corner to use for realigning the output files
%   in an image editing aplication.
% SPLITPRINT(FIG,FILENAME,TYPES,FORMATS,ALIGN,VISIBLE) will make the axes 
%   in the specified figure visible.
% SPLITPRINT(FIG,FILENAME,TYPES,FORMATS,ALIGN,VISIBLE,OPTIONS) will print 
%   using the additional PRINT options specified
% SPLITPRINT(FIG,FILENAME,TYPES,FORMATS,ALIGN,VISIBLE,OPTIONS,CMAP) will print 
%   using the colormap in CMAP
%
%  Example:
%
%  splitprint(gcf,... %separate the current figure
%     'loco)117',... %filenames will begin with 'disp2'
%     {{'line';'text'},{'surface';'patch';'image'}}, ...% types of objects
%     {'-depsc2','-dtiff'},... %file formats
%     0,... %alignment mark will not be added
%     [1 0],... %axes in first figure will be visible
%     {'-adobecset','-r400'})
%,'-cmyk'}})
   % , ... %second figure is printed in cmyk with -r400 option
    %cmap); % colormap
%
%  Separates the current figure (GCF) into line & text objects and surface & 
%  patch objects.  The new figures will be printed to files whose names begin 
%  with 'disp2', their formats will be printed to encapsulated color PostScript
%  and tiff, respectively.  The alignment mark will be added.  The axes in 
%  the line/text figure will be visible.
%

error(nargchk(0,8,nargin))
error(nargoutchk(0,2,nargout))

% Save temporary copy of figure; will be deleted later
saveas(gcf, [ 'temp.fig'], 'fig');

%default values, modify to set your defaults
fig = gcf;
printtofile = logical(1); %print to file by default, if this is 0, just splits figure
filename = 'splitprint';
types = {{'text';'line'}, {'patch';'surface';'image'}};
fmts = {'-depsc2','-dtiff','-dpng'};
alignmarks = logical(0);
visaxes = logical([1 0]);
options = {'-dtiff',{'-r400','-cmyk'}}; % default eps for first fig; 400 dpi for second


switch nargin
case 1
    if ishandle(varargin{1}) & strcmpi(get(varargin{1},'Type'),'figure')
        fig = varargin{1};
        printtofile = logical(0);
    elseif isstr(varargin{1}) | iscellstr(varargin{1})
        filename = varargin{1};
    else
        error('Single input case, but input is not a figure handle, string, or cell array of strings.')
    end    
case 2
    [fig,filename] = deal(varargin{:});
case 3
    [fig,filename,types] = deal(varargin{:});
    
    %if the number of types is greater than the default number of formats (2)
    %pad the 'fmts' array with '-depsc2'
    if prod(size(types)) > 2,
        l = prod(size(types));
        [fmts{3:l}] = deal('-depsc2');
    end

case 4
    [fig,filename,types,fmts] = deal(varargin{:});
    visaxes = logical([1 zeros(1,prod(size(types))-1)]);
case 5
    [fig,filename,types,fmts,alignmarks] = deal(varargin{:});
    visaxes = logical([1 zeros(1,prod(size(types))-1)]);
case 6
    [fig,filename,types,fmts,alignmarks,visaxes] = deal(varargin{:});
case 7
    [fig,filename,types,fmts,alignmarks,visaxes,options] = deal(varargin{:});
case 8
    [fig,filename,types,fmts,alignmarks,visaxes,options,cmap] = deal(varargin{:});
end

%check figure handle
if ~ishandle(fig)
    error('First input is not a handle.')
elseif ~strcmpi(get(fig,'Type'),'figure')
    error('First input is a handle, but not the handle of a figure.')
elseif prod(size(fig)) > 1,
    error('Too many figure handles in input')
end

%check filename
if ~isempty(filename) & printtofile
    if ~isstr(filename) & ~iscellstr(filename)
        error('Second input is not a string or cell array of strings.')
    end
elseif isempty(filename)
    printtofile = logical(0);
end

%check types
actualtypes = {'image';'light';'line';'patch';'surface';'text'};
for i = 1:length(types)
    for j = 1:prod(size(types{i})),
        if ~any(strcmp(types{i}{j},actualtypes))
            error(['Cell ' num2str(i) ' of the TYPES input does not match any of the accepted types' ])
        end
    end
end

%check fmts
if printtofile & length(fmts) ~= length(types)
    error('The number of print formats does not match the number of groups of types.')
elseif any(~strncmp(fmts,'-',1))
    error('The formats must begin with ''-''.')    
end


%determine number of new figures
nfigs = prod(size(types));

%determine the number og types of objects in ecah figure
ntypes = cellfun('prodofsize',types);

%vector for handles of new figures
newfig = zeros(nfigs,1);

origpos = get(fig,'Position');
origunits = get(fig,'Units');

for i = 1: nfigs,
    close(fig); fig = open([ '/temp.fig']);
    %get the handles for all of the axes
    origaxes = findall(fig,'type','axes');
    
    %create empty figures, same size as original
    newfig(i) = figure('Units',origunits,'Position',origpos);
    set(gcf,'PaperPositionMode',get(fig,'PaperPositionMode'));
    %copy all of the axes into the new figure

    %see solution # 30292 - when copying multiple objects
    %the order of the children is reversed; using flipud is the
    %workaround
    % The following line is needed for some figures (eg three way compare);
    %newaxes{i} = copyobj(flipud(origaxes),newfig(i));

    % Sometimes flipud causes problems; comment the above
    % line and uncomment the following on
    newaxes{i} = copyobj(origaxes,newfig(i));

    %make sure the axes don't rescale when children are deleted 
    set(newaxes{i},...
        'XLimMode','manual',...
        'YLimMode','manual',...
        'ZLimMode','manual');
    %if the corresponding element of 'visaxes' is zero, make the 
    %axes invisible
    if ~visaxes(i),
        set(newaxes{i},'Visible','off')
    end
    
    % determine which children of the axes are the types to keep
    ch = get(newaxes{i},'Children');
    if iscell(ch),
        ch = vertcat(ch{:});
    end
    
    istype = zeros(length(ch),ntypes(i));    
    for j = 1:ntypes(i),
        istype(:,j) = strcmp(types{i}{j},get(ch,'Type')) | strcmp('light',get(ch,'Type'));        
    end
    keep = any(istype,2);
    
    %delete the objects which don't match the possible types
    delete(ch(~keep));
        
end

if nargin==8
    if ~isempty(cmap)
      colormap(cmap);
    end
end


if alignmarks,   
    hidaxes = zeros(nfigs,1);
    for i = 1:nfigs
        %create an invisible axes to position the alignment marks in
        hidaxes(i) = axes('Parent',figure(newfig(i)),...
                            'Units','normalized',...
                            'Position',[0 0 1 1],...
                            'Visible','off');
        allaxes{i} = findall(newfig(i),'Type','axes');
        set(newfig(i),'Children',[newaxes{i};hidaxes(i)]);
        textreg = text('Units','normalized',...
                              'String','+',...
                              'Position',[.975 .975 0]);

    end
end

if printtofile,
    for i = 1:nfigs
        figure(newfig(i));
        pstr = ['print ' fmts{i}];
        if ~isempty(options{i})
            tmp=options{i};
            if iscell(tmp)
              for z = 1:size(tmp,2)
                pstr=[pstr ' ' tmp{z}];
              end     
            else
                pstr=[pstr ' ' tmp];
            end
        end
        if size(filename,1) > 1,
            % Additional quotes needed when there are 
            % spaces in the pathname
            pstr = [pstr ' ' '''' filename{i} ''''];        
        else
            pstr = [pstr ' ' '''' filename ''''];        
        end
        eval(pstr)

    end
end

for i = 1:nfigs
    close(figure(newfig(i)))
end
% open([ 'temp.fig']);
delete([ 'temp.fig']);
