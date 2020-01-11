function varargout = viewover(varargin)
% VIEWOVER displays a fused overlay image and adds a few user controls:
%   Right mouse button win/lev color controls (on overlay)
%   Left mouse button win/lev transparency controls (on overlay)
%   o key makes overlay completely opaque
%   u key makes overlay completely transparent
%   w key resets overlay win/lev for color and transparency to default
%   Up/Down, J/K keys scroll through repetitions (T)
%   Left/Right, H/L scroll through slices (Z)
%   +/- keys zoom in/out
%
% Usage: viewover(im, bkg, 'Name', Value, 'Name', Value ...)
%        viewover(im, bkg, 'Title')
%        viewover(im, 'Name', Value, 'Name', Value, ...)
%
%   where im is a 3D or 4D array to overlay, [X, Y, (Z or T,) T]
%           (if complex magnitude is used)
%         bkg is a 2D or 3D array for underlay, [X, Y, (Z)]
%           (if complex magnitude is used)
%           (if omitted or empty, im is shown without underlay)
%         Name and Value are optional key-value pairs chosen from: 
%           - 'title',    <string>
%           - 'zoom',     <a positive real number greater than 0>
%           - 'colormap', <string specifying a valid colormap function>
%           - 'colorbar', <logical>
%           - 'interp',   <'none', 'bilinear' or 'bicubic'>
%               (interpolates to bkg size)
%           - 'bkgSize',  <a positive real number greater than 0>
%               (only used if bkg array is not provided as second input)
%
% Modified by Jack Miller & Angus Lau 
% 08/2019, Keith Michel


% Some versions of matlab may need to uncomment the below: 
%warning('off','all')
%feature('usehg2',0);
%warning('on','all')
global ImageTag;
if ~nargin, help(mfilename); return; end
if ischar(varargin{1})
    % If first argument passed is a character string, invoke callback routines
    try
        if (nargout)
            [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
        else
            feval(varargin{:}); % FEVAL switchyard
        end
    catch ME
        %disp(ME.identifier);
        disp(ME.message);
        disp(ME.cause);
    end
else
    % Otherwise, open a new figure and post the image
    % For complex data, use the magnitude
    imstack = squeeze(varargin{1});
    args    = {};
    bkg     = [];
    if nargin > 1
        if ischar(varargin{2})
            args = {varargin{2:end}};
        elseif ~ischar(varargin{2})
            bkg = nscale(squeeze(varargin{2}));
        end
    end
    
    ImageTag = 'MyImage';
    if nargin == 2
        if ischar(varargin{2})
            ImageTag = varargin{2};
        end
    elseif nargin == 3
        if ischar(varargin{3}) && ~ischar(varargin{2})
            ImageTag = varargin{3};
        end
    elseif nargin > 3 && isempty(args)
        args = {varargin{3:end}};
    end
    
    if numel(args) ~= 2*floor(numel(args)/2)
        fprintf('Invalid number of inputs to %s \n\n', mfilename)
        help(mfilename)
        return
    end
    
    changeSomething = false;
    if ~isempty(args)
        changeSomething = true;
        p = inputParser;
        if isempty(bkg)
            defaultColormap = 'gray';
        else
            defaultColormap = 'hot';
        end
        expectedColormaps = {'cubehelix', 'linspecer', ...
            'jet','hot','hsv', 'gray','cool','spring', ...
            'summer','autumn','winter','bone','copper','pink', ...
            'Blues', 'BuGn', 'BuPu', 'GnBu', 'Greens', 'Greys', 'OrRd', ...
            'Oranges', 'PuBu', 'PuBuGn', 'PuRd', 'Purples', 'RdPu', ...
            'Reds', 'YlGn', 'YlGnBu', 'YlOrBr', 'YlOrRd'};
        
        defaultTitle    = inputname(1);
        defaultZoom     = 2;
        defaultColorbar = false; 
        defaultInterp   = 'none'; 
        expectedInterps = {'none', 'bilinear', 'bicubic'};
        if isdeployed
            defaultBkgSize = [256, 256];
        else
            defaultBkgSize = size(imstack);
            defaultBkgSize = defaultBkgSize(1:2);
        end
        
        addOptional(p,'colormap',defaultColormap,@(x) any(validatestring(x,expectedColormaps)));
        addOptional(p,'title',defaultTitle,@ischar);
        addOptional(p,'colorbar',defaultColorbar,@islogical)
        addOptional(p,'zoom',defaultZoom)
        addOptional(p,'interp',defaultInterp,@(x) any(validatestring(x,expectedInterps)))
        addOptional(p,'bkgSize',defaultBkgSize)
        
        parse(p, args{:});
        
        usedColormap = p.Results.colormap;
        useColorbar  = p.Results.colorbar;
        usedZoom     = p.Results.zoom; 
        useInterp    = p.Results.interp;
        useBkgSize   = p.Results.bkgSize;
        
        if isempty(usedZoom)
            zoomOnStartup = false;
        else
            zoomOnStartup      = true;
            zoomOnStartupLevel = num2str(usedZoom);
        end
        
        ImageTag = p.Results.title;
        if isempty(ImageTag)
            ImageTag = 'MyImage';
        end
    else
        if isdeployed
            useBkgSize = [256, 256];
        else
            useBkgSize = size(imstack);
            useBkgSize = useBkgSize(1:2);
        end
        useInterp = 'none';
        usedZoom  = 2;
        zoomOnStartupLevel = num2str(usedZoom);
    end
    
    
    % KAM 8/2019 - multidimensional overlays (3D or 4D) and underlays (2D or 3D)
    if ~isreal(imstack), imstack = abs(imstack); end
%     raw_images = imstack;
    imgSize = size(imstack);
    if isempty(bkg)
        switch ndims(imstack)
            case 3
                bkg = zeros([useBkgSize, imgSize(3)]);
            case 4
                bkg = zeros([useBkgSize, imgSize(3)]);
        end
    end
    bkgSize = size(bkg);
    if numel(bkgSize) > 3 || numel(bkgSize) < 2
        error('viewover:bkgSize', 'Background must be 2D or 3D')
    end
    if numel(imgSize) < 3
        error('viewover:imgSize', 'Overlay image must have at least 3 dimenstions')
    elseif numel(imgSize) == 4
        % 4D overlay => [X, Y, Z, T]
        if numel(bkgSize) == 3
            if bkgSize(3) ~= imgSize(3)
                error('viewOver:imgBkgMismatchZ', ...
                '3D background must match 3rd dim of overlay image')
            end
        end
        num_z      = imgSize(3);
        num_t      = imgSize(4);
        num_images = num_t;
        raw_images = imstack;
        imstack    = squeeze(imstack(:,:,1,:));
    elseif numel(imgSize) == 3 && numel(bkgSize) == 2
        % 3D overlay => [X, Y, T], 2D underlay => [X, Y]
        num_z      = 1;
        num_t      = imgSize(3);
        num_images = num_t;
        raw_images = reshape(imstack, imgSize(1), imgSize(2), num_z, num_t);
    elseif numel(imgSize) == 3 && numel(bkgSize) == 3
        % 3D overlay => [X, Y, Z], 3D underlay => [X, Y, Z]
        if bkgSize(3) ~= imgSize(3)
            error('viewOver:imgBkgMismatchZ', ...
                '3D background must match 3rd dim of overlay image')
        end
        num_z      = imgSize(3);
        num_t      = 1;
        num_images = num_z;
        raw_images = imstack;
    end
    
    if ~strcmpi(useInterp, 'none')
        tmp = zeros([bkgSize(1:2), num_images]);
        for ii = 1:num_images
            switch useInterp
                case 'bilinear'
                    tmp(:,:,ii) = imresize(imstack(:,:,ii), bkgSize(1:2), 'bilinear');
                case 'bicubic'
                    tmp(:,:,ii) = imresize(imstack(:,:,ii), bkgSize(1:2), 'bicubic');
            end
        end
        imstack = tmp;
        clear tmp
    end
    
    xresUndr = bkgSize(2);
    yresUndr = bkgSize(1);
    xresOver = size(imstack, 2);
    yresOver = size(imstack, 1);
    xres     = max([xresUndr, xresOver]);
    yres     = max([yresUndr, yresOver]);
    
    % Create the figure graphics object, set the basic characteristics
    
    % JJM
    % Detect primary monitor
    monitorPos = get(0,'MonitorPositions');
    monitorPri = monitorPos(1,:);
    monitorPri = monitorPri + [monitorPri(3)*0.10 monitorPri(4)*0.10 -monitorPri(3)*0.10 -monitorPri(4)*0.10];
    
    % AZL Jul-08-2009 changed toolbar to figure from 'none'
    hfig = figure(  'Tag', 'MyFigure', 'Toolbar', 'figure', ...
        'Resize', 'on', 'NumberTitle', 'off', ...
        'Name', sprintf('Viewimage: %s', ImageTag), ...
        'DoubleBuffer', 'on', ...
        'Position', [monitorPri([1 2]) xres yres] );
    if nargin < 2
        set(hfig, 'Colormap', gray(256))
    else
        set(hfig, 'Colormap', hot(256));
    end

    % Set up the callbacks for right mouse button Win/Lev controls
    set( hfig, 'WindowButtonDownFcn',   'viewover(''Press'', gcf, gca)' );
    set( hfig, 'WindowButtonUpFcn',     'viewover(''Release'', gcf, gca)' );
    
    % JJM - 27-10-14 - Modified to mousewheel through images 
    set( hfig, 'WindowScrollWheelFcn', @Scroll);
    
    % KAM 8/2019 - Add callbacks for keyboard buttons
    set(hfig, 'KeyReleaseFcn', @keyPress);
    
    % Hide the View, Insert, Tools, Window, and Help menus in the figure
    set(0, 'ShowHiddenHandles', 'on');
    set( findobj( hfig, 'Tag',   'figMenuView'  ), 'Visible', 'off');
    set( findobj( hfig, 'Tag',   'figMenuInsert'), 'Visible', 'off');
    set( findobj( hfig, 'Tag',   'figMenuTools' ), 'Visible', 'off');
    set( findobj( hfig, 'Label', '&Window'      ), 'Visible', 'off');
    set( findobj( hfig, 'Label', '&Help'        ), 'Visible', 'off');
    set(0, 'ShowHiddenHandles', 'off')
    
    % Add ZOOM menu with some sensible defaults 
    hmenu1 = uimenu( hfig, 'Label','Zoom');
    uimenu(hmenu1, 'Label', '1X', 'Callback', 'viewover(''zoomit'', gcf, gca, 1)', ...
        'Checked', 'on' );
    uimenu(hmenu1, 'Label', '2X', 'Callback', 'viewover(''zoomit'', gcf, gca, 2)' );
    uimenu(hmenu1, 'Label', '3X', 'Callback', 'viewover(''zoomit'', gcf, gca, 3)' );
    uimenu(hmenu1, 'Label', '4X', 'Callback', 'viewover(''zoomit'', gcf, gca, 4)' );
    uimenu(hmenu1, 'Label', '6X', 'Callback', 'viewover(''zoomit'', gcf, gca, 6)' );
    uimenu(hmenu1, 'Label', '8X', 'Callback', 'viewover(''zoomit'', gcf, gca, 8)' );
    
    % Add WINDOW/LEVEL menu
    hmenu2 = uimenu( hfig, 'Label','W/L');
    uimenu(hmenu2, 'Label', 'AutoW/L', 'Tag', 'AutoWL', 'Checked', 'on', 'Callback', 'viewover(''autowl'', gcf, gca)' );
    uimenu(hmenu2, 'Label', 'To Workspace', 'Callback',   'viewover(''wl_to_ws'', gcf, gca)' );
    uimenu(hmenu2, 'Label', 'From Workspace', 'Callback',   'viewover(''wl_from_ws'', gcf, gca)' );
    
    % Add PLOT menu
    hmenu3 = uimenu( hfig, 'Label','Plot');
    uimenu(hmenu3, 'Label', 'Plot Row',    'Callback',   'viewover(''plot_rc'', gcf, gca, 2)' );
    uimenu(hmenu3, 'Label', 'Plot Column', 'Callback',   'viewover(''plot_rc'', gcf, gca, 1)' );
    if num_images > 1
        uimenu(hmenu3, 'Label', 'Plot Through', 'Callback', 'viewover(''plot_through'', gcf, gca)' );
        uimenu(hmenu3, 'Label', 'Plot ROI Through', 'Callback', 'viewover(''plot_roi_through'', gcf, gca)' );
    end
    
    % Add MATH menu
    roi_exist = 0;
    hmenu4 = uimenu( hfig, 'Label', 'Math' );
    uimenu(hmenu4, 'Label', 'Signal ROI', 'Callback', 'viewover(''get_roi'', 1)' );
    uimenu(hmenu4, 'Label', 'Noise ROI',  'Callback', 'viewover(''get_roi'', 2)' );
    %uimenu(hmenu4, 'Separator', 'on' );
    uimenu(hmenu4, 'Label', 'Calc SNR', 'Callback', 'viewover(''calc_snr'',gca,roi1,roi2)');
    
    % Add COLORMAP menu
    hmenu5 = uimenu( hfig, 'Label', 'Colormap' );
    uimenu(hmenu5, 'Label', 'Cubehelix', 'Callback', 'viewover(''set_colormap'', gcf, gca, ''cubehelix'')' );
    uimenu(hmenu5, 'Label', 'Linspecer', 'Callback', 'viewover(''set_colormap'', gcf, gca, ''linspecer'')' );
    uimenu(hmenu5, 'Label', 'Jet',       'Callback', 'viewover(''set_colormap'', gcf, gca, ''jet'')' );
    uimenu(hmenu5, 'Label', 'Hot',       'Callback', 'viewover(''set_colormap'', gcf, gca, ''hot'')' );
    uimenu(hmenu5, 'Label', 'hsv',       'Callback', 'viewover(''set_colormap'', gcf, gca, ''hsv'')' );
    uimenu(hmenu5, 'Label', 'gray',      'Callback', 'viewover(''set_colormap'', gcf, gca, ''gray'')' );
    uimenu(hmenu5, 'Label', 'cool',      'Callback', 'viewover(''set_colormap'', gcf, gca, ''cool'')' );
    uimenu(hmenu5, 'Label', 'spring',    'Callback', 'viewover(''set_colormap'', gcf, gca, ''spring'')' );
    uimenu(hmenu5, 'Label', 'summer',    'Callback', 'viewover(''set_colormap'', gcf, gca, ''summer'')' );
    uimenu(hmenu5, 'Label', 'autumn',    'Callback', 'viewover(''set_colormap'', gcf, gca, ''autumn'')' );
    uimenu(hmenu5, 'Label', 'winter',    'Callback', 'viewover(''set_colormap'', gcf, gca, ''winter'')' );
    uimenu(hmenu5, 'Label', 'bone',      'Callback', 'viewover(''set_colormap'', gcf, gca, ''bone'')' );
    uimenu(hmenu5, 'Label', 'copper',    'Callback', 'viewover(''set_colormap'', gcf, gca, ''copper'')' );
    uimenu(hmenu5, 'Label', 'pink',      'Callback', 'viewover(''set_colormap'', gcf, gca, ''pink'')' );
    if ~verLessThan('matlab', '8.4')
        uimenu(hmenu5, 'Label', 'Parula',       'Callback', 'viewover(''set_colormap'', gcf, gca, ''parula'')' );
    end
    uimenu(hmenu5, 'Label', 'Blues',    'Callback', 'viewover(''set_colormap'', gcf, gca, ''Blues'')' );
    uimenu(hmenu5, 'Label', 'BuGn',     'Callback', 'viewover(''set_colormap'', gcf, gca, ''BuGn'')' );
    uimenu(hmenu5, 'Label', 'BuPu',     'Callback', 'viewover(''set_colormap'', gcf, gca, ''BuPu'')' );
    uimenu(hmenu5, 'Label', 'GnBu',     'Callback', 'viewover(''set_colormap'', gcf, gca, ''GnBu'')' );
    uimenu(hmenu5, 'Label', 'Greens',   'Callback', 'viewover(''set_colormap'', gcf, gca, ''Greens'')' );
    uimenu(hmenu5, 'Label', 'Greys',    'Callback', 'viewover(''set_colormap'', gcf, gca, ''Greys'')' );
    uimenu(hmenu5, 'Label', 'OrRd',     'Callback', 'viewover(''set_colormap'', gcf, gca, ''OrRd'')' );
    uimenu(hmenu5, 'Label', 'Oranges',  'Callback', 'viewover(''set_colormap'', gcf, gca, ''Oranges'')' );
    uimenu(hmenu5, 'Label', 'PuBu',     'Callback', 'viewover(''set_colormap'', gcf, gca, ''PuBu'')' );
    uimenu(hmenu5, 'Label', 'PuBuGn',   'Callback', 'viewover(''set_colormap'', gcf, gca, ''PuBuGn'')' );
    uimenu(hmenu5, 'Label', 'PuRd',     'Callback', 'viewover(''set_colormap'', gcf, gca, ''PuRd'')' );
    uimenu(hmenu5, 'Label', 'Purples',  'Callback', 'viewover(''set_colormap'', gcf, gca, ''Purples'')' ); 
    uimenu(hmenu5, 'Label', 'RdPu',     'Callback', 'viewover(''set_colormap'', gcf, gca, ''RdPu'')' );
    uimenu(hmenu5, 'Label', 'Reds',     'Callback', 'viewover(''set_colormap'', gcf, gca, ''Reds'')' );
    uimenu(hmenu5, 'Label', 'YlGn',     'Callback', 'viewover(''set_colormap'', gcf, gca, ''YlGn'')' );
    uimenu(hmenu5, 'Label', 'YlGnBu',   'Callback', 'viewover(''set_colormap'', gcf, gca, ''YlGnBu'')' ); 
    uimenu(hmenu5, 'Label', 'YlOrBr',   'Callback', 'viewover(''set_colormap'', gcf, gca, ''YlOrBr'')' );
    uimenu(hmenu5, 'Label', 'YlOrRd',   'Callback', 'viewover(''set_colormap'', gcf, gca, ''YlOrRd'')' );
    
    % Add INTERP menu
    hmenu6 = uimenu( hfig, 'Label', 'Interp' );
    uimenu(hmenu6, 'Label', 'none',     'Callback', 'viewover(''interp_image'', gcf, gca, ''none'', true)', ...
        'Checked', 'on');
    uimenu(hmenu6, 'Label', 'bilinear', 'Callback', 'viewover(''interp_image'', gcf, gca, ''bilinear'', true)' );
    uimenu(hmenu6, 'Label', 'bicubic',  'Callback', 'viewover(''interp_image'', gcf, gca, ''bicubic'', true)' );
    
    % Initialize images
    im  = imstack(:,:,1);
    bkg = uint8(round(255*bkg));
    if any(bkg(:) > 0)
        imalpha  = imstack;
    else
        imalpha  = true(size(imstack));
    end
    imstruct = struct('imstack', imstack, 'imalpha', imalpha);
    alphas   = imalpha(:,:,1);
    
    % Create the background axes
    baxes = axes( 'Tag', 'MyBkgAxes', ...
        'Position', [ 0 0 1 1], ...
        'Visible', 'off', ...
        'CLimMode', 'auto', ...
        'Xlim', [1-.5 xresUndr+.5], 'Ylim', [1-.5 yresUndr+.5], ...
        'YDir', 'reverse' );
    if verLessThan('matlab', '8.4')
        bimage = subimage(bkg(:,:,1), gray(256));
        set(bimage, 'Tag', 'Background')
    else
        bimage = image('Tag', 'Background', 'CData', bkg(:,:,1), ...
            'CDataMapping', 'direct');
        colormap(baxes, gray(256));
    end
    
    % Create the axes graphics object
    haxes = axes(   'Tag', 'MyAxes', ...
        'Position', [ 0 0 1 1], ...
        'Visible', 'off', ...
        'CLimMode', 'auto', ...
        'Xlim', [1-.5 xresOver+.5], 'Ylim', [1-.5 yresOver+.5], ...
        'YDir', 'reverse' );
    setappdata(haxes, 'position_listener', linkprop([baxes, haxes],'Position'))
    
    % Create the image graphics object
    himage = image( 'Tag', ImageTag, 'CData', im, 'UserData', imstruct, ...
        'CDataMapping', 'scaled');
    if verLessThan('matlab','8.4.0') 
    else
        setappdata(hfig.Number,'CData',im); 
        setappdata(hfig.Number,'UserData',imstruct);
    end 
    
    % If there's no underlay, disable transparency
    if any(bkg(:) > 0)
        set(haxes,  'ALim', [min(alphas(:)), max(alphas(:))])
        set(himage, 'AlphaData', alphas, 'AlphaDataMapping', 'scaled')
    else
        set(haxes, 'ALim', [0 1])
        set(himage, 'AlphaData', alphas, 'AlphaDataMapping', 'none')
    end
    
    axes_UD.himage     = himage;
    axes_UD.baxes      = baxes;
    axes_UD.raw_images = raw_images;
    axes_UD.useInterp  = useInterp;
    axes_UD.useZoom    = usedZoom;
    axes_UD.bkg        = bkg;
    axes_UD.num        = 1;          % t dim current index
    axes_UD.num2       = 1;          % z dim current index
    axes_UD.count      = num_images; % t dim total size
    if numel(imgSize) == 4
        axes_UD.count2 = num_z;      % z dun total size
    else
        axes_UD.count2 = 1;
    end
    axes_UD.himage     = himage;
    axes_UD.bimage     = bimage;
    set( haxes, 'UserData', axes_UD );
    
    if changeSomething
        colormap(haxes, usedColormap);
        if useColorbar
            colorbar('peer',haxes);
        end
    end
    
    zoomOnStartup = true;
    if zoomOnStartup
        eval(['viewover(''zoomit'', gcf, gca, ',zoomOnStartupLevel,')']);
    end
    movegui(gcf,'onscreen');
    
    
    %If you uncommented the things at the begining, you probably want to re-enable hg2 later. 
    %feature('usehg2',1);
end

end


%=======================================================================
% CALLBACK  and support functions
%=======================================================================

% Repaint the image
% ----------------------------------------------------
function varargout = update_image( hfig, haxes, himage, varargin)
axes_UD   = get( haxes, 'UserData');
im_num    = axes_UD.num;
im_num2   = axes_UD.num2;
bimage    = axes_UD.bimage;

imstruct = get( himage,'UserData');
imstack  = imstruct.imstack;
imalpha  = imstruct.imalpha;
bkg      = axes_UD.bkg;
imgSize  = size(axes_UD.raw_images);
im       = double(imstack(:,:,im_num));
alphas   = imalpha(:,:,im_num);

set(himage, 'CData', im, 'AlphaData', alphas);
if size(bkg, 3) > 1
    if numel(imgSize) == 3
        bkgNum = im_num;
    else
        bkgNum = im_num2;
    end
    curBkg = bkg(:,:,bkgNum);
    if verLessThan('matlab', '8.4')
        curBkg = ind2rgb(curBkg, gray(256));
    end
    set(bimage, 'CData', curBkg)
end


baxes = axes_UD.baxes;
% set(baxes, 'position', get(haxes, 'position'))
axis([baxes, haxes], 'tight')

end

% Keypress functions
% ----------------------------------------------------------
function varargout = keyPress(~, event)

hfig     = gcf;
haxes    = gca;
axes_UD  = get(haxes, 'UserData');
himage   = axes_UD.himage; 
imstruct = get(himage, 'UserData');
imalpha  = imstruct.imalpha;

switch event.Key
    case {'downarrow', 'j'}
        viewover('DoScroll', hfig, haxes, -1);
    case {'uparrow', 'k'}
        viewover('DoScroll', hfig, haxes, +1);
    case 'pagedown'
        viewover('DoScroll', hfig, haxes, -10);
    case 'pageup'
        viewover('DoScroll', hfig, haxes, +10);
    case {'l', 'rightarrow'}
        viewover('zscroll_image', hfig, haxes, +1);
    case {'h', 'leftarrow'}
        viewover('zscroll_image', hfig, haxes, -1);
    case {'add', 'equal'}
        viewover('zoomit', hfig, haxes, [], +1);
    case {'subtract', 'hyphen'}
        viewover('zoomit', hfig, haxes, [], -1);
    case 'w'
        h_autowl = findobj( hfig, 'Tag', 'AutoWL');
        set( h_autowl, 'Checked', 'off')
        viewover('autowl', hfig, haxes)
    case 'o'
        if strcmpi(get(himage, 'AlphaDataMapping'), 'scaled')
            set(haxes, 'ALim', min(imalpha(:)) - [1 0])
        end
        viewover('update_wl', hfig, haxes)
    case 'u'
        if strcmpi(get(himage, 'AlphaDataMapping'), 'scaled')
            set(haxes, 'ALim', max(imalpha(:)) + [0 1])
        end
        viewover('update_wl', hfig, haxes)
end

end

% Z-scroll function
% ----------------------------------------------------------
function varargout = zscroll_image( hfig, haxes, sc, varargin)
axes_UD = get(haxes, 'UserData');
num2    = axes_UD.num2 + sc;
count2  = axes_UD.count2;
num2    = max([num2, 1]);
num2    = min([num2, count2]);
axes_UD.num2 = num2;
set(haxes, 'UserData', axes_UD)
useInterp = axes_UD.useInterp;
set(hfig,  'Name', sprintf('Viewimage: Im(z)# %d/%d', num2, count2));
viewover('interp_image', hfig, haxes, useInterp, false);

end

% Interp menu function
% ----------------------------------------------------------
function varargout = interp_image( hfig, haxes, useInterp, changeInterp, varargin)
axes_UD  = get(haxes, 'UserData');
himage   = axes_UD.himage;
bimage   = axes_UD.bimage;
xdata    = get( bimage, 'XData');
ydata    = get( bimage, 'YData');
xresUndr = xdata(2)-xdata(1)+1;
yresUndr = ydata(2)-ydata(1)+1;
bkgSize  = [yresUndr, xresUndr];

imstack = squeeze(axes_UD.raw_images);
if ndims(imstack) > 3
    imstack = squeeze(imstack(:,:,axes_UD.num2,:));
end
num_images = size(imstack, 3);
if ~strcmpi(useInterp, 'none')
    tmp = zeros([bkgSize, num_images]);
    for ii = 1:num_images
        switch useInterp
            case 'bilinear'
                tmp(:,:,ii) = imresize(imstack(:,:,ii), bkgSize, 'bilinear');
            case 'bicubic'
                tmp(:,:,ii) = imresize(imstack(:,:,ii), bkgSize, 'bicubic');
        end
    end
    imstack = tmp;
    clear tmp
end
imstruct = get(himage, 'UserData');
imalpha  = imstruct.imalpha;
if all(imalpha(:) == 1)
    imalpha = true(size(imstack));
else
    imalpha = imstack;
end
imstruct = struct('imstack', imstack, 'imalpha', imalpha);
set(himage, 'UserData', imstruct)
axes_UD.useInterp = useInterp;
set(haxes, 'UserData', axes_UD)
viewover('update_image', hfig, haxes, himage);
if changeInterp
    autocheckmark( gcbo );
end

end

% Zoom menu functions
% ----------------------------------------------------------
function varargout = zoomit( hfig, haxes, zoomval, varargin)
global ImageTag;
himage   = findobj( haxes, 'Tag', ImageTag);
xdata    = get( himage, 'XData');
ydata    = get( himage, 'YData');
xresOver = xdata(2)-xdata(1)+1;
yresOver = ydata(2)-ydata(1)+1;

axes_UD  = get(haxes, 'UserData');
if ~isempty(zoomval)
    axes_UD.useZoom = zoomval;
    autocheckmark( gcbo );
else
    zoomval = axes_UD.useZoom + varargin{1};
    axes_UD.useZoom = zoomval;
end
bimage   = axes_UD.bimage;
xdata    = get( bimage, 'XData');
ydata    = get( bimage, 'YData');
xresUndr = xdata(2)-xdata(1)+1;
yresUndr = ydata(2)-ydata(1)+1;

xres = max([xresUndr, xresOver]);
yres = max([yresUndr, yresOver]);

position = get(hfig, 'Position');
upper_edge = position(2)+position(4);
position(3) = xres*zoomval;
position(4) = yres*zoomval;
position(2) = upper_edge - position(4);
set(hfig, 'Position', position);
set(haxes, 'UserData', axes_UD);
% autocheckmark( gcbo );

end

% Menu based W/L functions
% ---------------------------------------------------------------
function varargout = autowl( hfig, haxes, varargin)
h_autowl = findobj( hfig, 'Tag', 'AutoWL');
switch get( h_autowl, 'Checked')
    case 'on'
        set( h_autowl, 'Checked', 'off');
        set( haxes, 'CLimMode', 'manual', 'ALimMode', 'manual');
    case 'off'
        set( h_autowl, 'Checked', 'on');
        set( haxes, 'CLimMode', 'auto', 'ALimMode', 'auto');
        viewover( 'update_wl', hfig, haxes);
end

end

function varargout = wl_to_ws( hfig, haxes, varargin)
lims = get(haxes, 'CLim');
win  = lims(2) - lims(1);
lev  = (lims(1) + lims(2)) / 2;
lims = get(haxes, 'Alim');
awin = lims(2) - lims(1);
alev = (lims(1) + lims(2)) / 2;
assignin( 'base', 'win', win);
assignin( 'base', 'lev', lev);
assignin( 'base', 'awin', awin);
assignin( 'base', 'alev', alev);

end

function varargout = wl_from_ws( hfig, haxes, varargin)
win  = evalin( 'base', 'win');
lev  = evalin( 'base', 'lev');
awin = evalin( 'base', 'awin');
alev = evalin( 'base', 'alev');
set(haxes, 'CLim', [lev-win/2 lev+win/2], ...
           'Alim', [alev-awin/2 alev+awin/2])
viewover( 'update_wl', hfig, haxes);

end

% Update the W/L display
% -----------------------------
function varargout = update_wl( hfig, haxes, varargin)
clims = get( haxes, 'CLim');
alims = get( haxes, 'ALim');
set( hfig, 'Name', sprintf('Viewimage W:%.3g L:%+.3g, Wt:%.3g Lt:%+.3g', ...
    clims(2)-clims(1), (clims(1)+clims(2))/2, ...
    alims(2)-alims(1), (alims(1)+alims(2))/2));

end

% Mouse based W/L and image scroll functions
% -----------------------------------------------------------------------
function varargout = Press(hfig, haxes, varargin)
set( hfig, 'WindowButtonMotionFcn', 'viewover(''Move'', gcf, gca, 0)' );
viewover('Move', gcf, gca, 1);     % Store initial mouse pointer location by calling 'Move' routine

end

% Scrollwheel function 
% -----------------------------------------------------------------------
function varargout = Scroll(object, event) 
viewover('DoScroll',gcf,gca,event.VerticalScrollCount);

end

function varargout = DoScroll(hfig, haxes, sc, varargin)
persistent LASTX LASTY MOVED 
ptr=get(hfig, 'CurrentPoint'); 
        LASTX = ptr(1); 
        LASTY = ptr(2); 

dx=ptr(1)-LASTX; 
dy=sc; 

axes_UD     = get( haxes, 'UserData');
himage      = axes_UD.himage;
imstruct    = get( himage,'UserData');
imstack     = imstruct.imstack;
imalpha     = imstruct.imalpha;
im_num      = axes_UD.num;
nimages     = axes_UD.count;
new_im_num  = im_num - dy;
new_im_num  = max(new_im_num,1);
new_im_num  = min(new_im_num, nimages);
axes_UD.num = new_im_num;
set( hfig,  'Name', sprintf('Viewimage: Im(t)# %d/%d', axes_UD.num, axes_UD.count ));
set( haxes, 'UserData', axes_UD );
viewover( 'update_image', hfig, haxes, himage);
   
end

% Callback to do nothing on mouse release 
% --------------------
function varargout = Release(hfig, haxes, varargin)
set( hfig, 'WindowButtonMotionFcn', '' );    % This makes the callback a do-nothing operation
if verLessThan('matlab','8.4.0')
    axes_UD     = get( haxes, 'UserData');
    ImageTag = get( axes_UD.himage, 'Tag' );
    set( hfig, 'Name', sprintf('Viewimage: %s', ImageTag) );

else
    global ImageTag;
    set( hfig, 'Name', sprintf('Viewimage: %s', ImageTag) );

end

end

function varargout = Move(hfig, haxes, init, varargin)
persistent LASTX LASTY MOVED
ptr = get(hfig,  'CurrentPoint');
if init == 1
    LASTX = ptr(1);
    LASTY = ptr(2);
    
else
    dx = ceil(ptr(1)-LASTX);
    dy = ceil(ptr(2)-LASTY);
    
    switch get(hfig,  'SelectionType')
        
        case 'normal'    % Mouse moved with LEFT button pressed -- change overlay alpha
            alims = get( haxes, 'ALim');
            win = alims(2)-alims(1);
            lev = (alims(1)+alims(2))/2;
            newwin = win * 10^(dx/400);
            newlev = lev - (win*dy)/200;
            set(haxes, 'ALim', [newlev-newwin/2 newlev+newwin/2]);
            set( findobj( hfig, 'Tag', 'AutoWL'), 'Checked', 'off');
            viewover( 'update_wl', hfig, haxes);

        case 'alt'    % Mouse moved with RIGHT button pressed -- change overlay W/L
            clims = get( haxes, 'CLim');
            win = clims(2)-clims(1);
            lev = (clims(1)+clims(2))/2;
            newwin = win * 10^(dx/400);
            newlev = lev - (win*dy)/200;
            set(haxes, 'CLim', [newlev-newwin/2 newlev+newwin/2]);
            set( findobj( hfig, 'Tag', 'AutoWL'), 'Checked', 'off');
            viewover( 'update_wl', hfig, haxes);
            
    end
    LASTX = ptr(1);
    LASTY = ptr(2);
end

end

% Menu based Colormap functions
% ---------------------------------------------------------------
% AZL 11.02.2013
function varargout = set_colormap( hfig, haxes, map, varargin )
switch map
    case {'jet','hot','hsv','gray','cool','spring','summer','autumn', ...
            'winter','bone','copper','pink','parula','cubehelix','linspecer'}
        try
            eval(sprintf('cmap = %s(256);', map))
        catch
            warning('viewover:set_colormap', 'Failed to change colormap')
            cmap = hot(256);
        end
    otherwise
        try
            map  = ['*', map];
            cmap = feval('brewermap', 256, map);
        catch
            warning('viewover:set_colormap', 'brewermap colormap failed')
            cmap = hot(256);
        end
end
colormap(haxes, cmap);

end

% Menu based Math functions
% ---------------------------------------------------------------
function varargout = get_roi( roi_val, varargin)
%roi = roipoly;
roi_free = imfreehand;
roi = roi_free.createMask();
assignin( 'base', strcat(['roi' num2str(roi_val)]), roi);

end

function varargout = calc_snr( haxes, roi1, roi2, varargin )
axes_UD = get( haxes, 'UserData');
im_num = axes_UD.num;
himage = axes_UD.himage;
imstruct = get( himage,'UserData');
imstack  = imstruct.imstack;
im = double(imstack(:,:,im_num));

s = abs(mean(im(roi1)));
n = abs(std(im(roi2)));
snr = s/n;

assignin( 'base', 'snr', snr);

fprintf('------------------\n')
fprintf('signal = %.4f\n',s)
fprintf('noise  = %.4f\n',n)
fprintf('snr    = %.4f\n',snr)
fprintf('------------------\n')

end

% Automatic checkmarks for menu items: checkmark the 'hitem', uncheck the others
% ------------------------------------------------------------------------------
function varargout = autocheckmark( hitem,  varargin)
hmenu     = get( hitem, 'Parent' );
hchildren = get( hmenu, 'Children');
for i = 1:length(hchildren)
    set( hchildren(i), 'Checked', 'off');    % Uncheck all the menu items
end
set( hitem, 'Checked', 'on');    % Check the menu item that made the callback

end

% Plot a row of column of data from the image
% ----------------------------------------------------
function varargout = plot_rc( hfig, haxes, direction, varargin)
[x, y] = ginput(1);
x = round(x);
y = round(y);

axes_UD = get(haxes, 'UserData');
himage  = axes_UD.himage;
im      = get(himage, 'CData');
num     = axes_UD.num;
num2    = axes_UD.num2;

figure;
switch direction
    case 2
        plot( im(round(y),:) );
        title(sprintf('Row %d for z=%d, t=%d', y, num2, num));
        xlabel('X Position')
    case 1
        plot( im(:,round(x)) );
        title(sprintf('Column %d for z=%d, t=%d', x, num2, num));
        xlabel('Y Position')
end

end

% Plot through-plane data from the image
% ----------------------------------------------------
function varargout = plot_through( hfig, haxes, varargin)
global ImageTag;

[voxelXidx, voxelYidx] = ginput(1);
voxelXidx = round(voxelXidx);
voxelYidx = round(voxelYidx);

axes_UD   = get( haxes, 'UserData');
voxelZidx = axes_UD.num2;
himage    = axes_UD.himage;
imstruct  = get( himage,'UserData');
imstack   = imstruct.imstack;

figure;

if strcmpi(ImageTag, 'Pyruvate and Lactate(x5)') && size(imstack, 2) == 2*size(imstack,1)
    % KAM faFuseIms shows pyr and lac(x5) side by side
    dx = ceil(size(imstack,2)/2);
    if voxelXidx > dx
        xx        = voxelXidx;
        voxelXidx = voxelXidx - dx;
    else
        xx = voxelXidx + dx;
    end
    pyrVoxel = squeeze(imstack(voxelYidx,voxelXidx,:));
    plot( pyrVoxel, '.-', 'MarkerSize', 9, 'LineWidth', 1.4, ...
        'Color', [0    0.4470    0.7410]);
    hold on
    lacVoxelx5 = squeeze(imstack(voxelYidx,xx,:));
    plot( lacVoxelx5, '.-', 'MarkerSize', 9, 'LineWidth', 1.4, ...
        'Color', [0.8500    0.3250    0.0980]);
    legend('Pyruvate', 'Lactate(x5)');
    if ~isdeployed
        assignin( 'base', 'pyrVoxel',   pyrVoxel);
        assignin( 'base', 'lacVoxelx5', lacVoxelx5);
        assignin( 'base', 'voxelXidx',  voxelXidx);
        assignin( 'base', 'voxelYidx',  voxelYidx);
        assignin( 'base', 'voxelZidx',  voxelZidx);
    else
        % save(sprintf('PyrLacVoxelDynamics-%s', datestr(now, 30)), ...
        %     'pyrVoxel', 'lacVoxelx5', 'voxelXidx', 'voxelYidx', 'voxelZidx')
        saveVars = {'pyrVoxel', 'lacVoxelx5', 'voxelXidx', 'voxelYidx', ...
            'voxelZidx'};
        uisave(saveVars, sprintf('PyrLacVoxelDynamics-%s', datestr(now, 30)));
    end
else
    dynVoxel = squeeze(imstack(voxelYidx,voxelXidx,:));
    plot( dynVoxel, '.-', 'markersize', 9, 'linewidth', 1.4);
    if ~isdeployed
        assignin( 'base', 'dynVoxel',   dynVoxel);
        assignin( 'base', 'voxelXidx',  voxelXidx);
        assignin( 'base', 'voxelYidx',  voxelYidx);
        assignin( 'base', 'voxelZidx',  voxelZidx);
    else
        % save(sprintf('VoxelDynamics-%s', datestr(now, 30)), ...
        %     'dynVoxel', 'voxelXidx', 'voxelYidx', 'voxelZidx')
        saveVars = {'dynVoxel', 'voxelXidx', 'voxelYidx', 'voxelZidx'};
        uisave(saveVars, sprintf('VoxelDynamics-%s', datestr(now, 30)));
    end
end
ylabel(sprintf('Magnitude at (%d, %d, %d)', voxelXidx, voxelYidx, voxelZidx));
grid on
set(gca, 'fontsize', 13)

end

% Draw an ROI and plot it through the image 
% ----------------------------------------------------
function varargout = plot_roi_through( hfig, haxes, varargin)
global ImageTag;

maskRoi = roipoly();

axes_UD   = get( haxes, 'UserData');
maskZidx  = axes_UD.num2;
himage    = axes_UD.himage;
imstruct  = get(himage, 'UserData');
imstack   = imstruct.imstack;

figure;

if strcmpi(ImageTag, 'Pyruvate and Lactate(x5)') && size(imstack, 2) == 2*size(imstack,1)
    % KAM faFuseIms shows pyr and lac(x5) side by side
    [~,mx] = find(maskRoi, 1, 'last');
    dx = ceil(size(imstack,2)/2);
    if mx > dx
        maskRoi = circshift(maskRoi, dx, 2);
        mask    = repmat(maskRoi, [1,1, size(imstack, 3)]);
    else
        mask = repmat(maskRoi, [1,1, size(imstack, 3)]);
    end
    pyrRoi   = squeeze(sum(sum(mask.*imstack))) ./ nnz(mask(:,:,1));
    mask     = circshift(mask, [0 dx]);
    lacRoix5 = squeeze(sum(sum(mask.*imstack))) ./ nnz(mask(:,:,1));
    plot( pyrRoi, '.-', 'MarkerSize', 9, 'LineWidth', 1.4, ...
        'Color', [0    0.4470    0.7410]);
    hold on
    plot( lacRoix5, '.-', 'MarkerSize', 9, 'LineWidth', 1.4, ...
        'Color', [0.8500    0.3250    0.0980]);
    legend('Pyruvate', 'Lactate(x5)');
    if ~isdeployed
        assignin( 'base', 'pyrRoi',   pyrRoi);
        assignin( 'base', 'lacRoix5', lacRoix5);
        assignin( 'base', 'maskRoi',  maskRoi(:,1:dx));
        assignin( 'base', 'maskZidx', maskZidx);
    else
        % save(sprintf('PyrLacRoiDynamics-%s', datestr(now, 30)), ...
        %     'pyrRoi', 'lacRoix5', 'maskRoi', 'maskZidx')
        saveVars = {'pyrRoi', 'lacRoix5', 'maskRoi', 'maskZidx'};
        uisave(saveVars, sprintf('PyrLacRoiDynamics-%s', datestr(now, 30)));
    end
else
    mask   = repmat(maskRoi, [1,1, size(imstack, 3)]);
    dynRoi = squeeze(sum(sum(mask.*imstack))) ./ nnz(maskRoi);
    plot( dynRoi, '.-', 'markersize', 9, 'linewidth', 1.4);
    if ~isdeployed
        assignin( 'base', 'dynRoi',   dynRoi);
        assignin( 'base', 'maskRoi',  maskRoi);
        assignin( 'base', 'maskZidx', maskZidx);
    else
        % save(sprintf('RoiDynamics-%s', datestr(now, 30)), ...
        %     'dynRoi', 'maskRoi', 'maskZidx')
        saveVars = {'dynRoi', 'maskRoi', 'maskZidx'};
        save(saveVars, sprintf('RoiDynamics-%s', datestr(now, 30)));
    end
end
ylabel('Mean ROI Magnitude');
grid on
set(gca, 'fontsize', 13)

end

%---------------------------------------------------------------
% Again, the below might be required on some versions of matlab 


% function im=getImage(haxes)
% %Required helper function for HG2 graphics 
% try
%     hg1=graphicsversion('handlegraphics'); 
% catch
%     hg1=99;
% end
% 
% if ~hg1
%     
% else %Old graphics
%     himage      = findobj( haxes, 'Tag', 'MyImage');
%     im= get(himage, 'CData');
% end

function [x,n] = nscale(x)
%NSCALE scales data to range of [0, 1]
%
%   Usage: [y, n] = nscale(x)
%
%       where x is array with data to be rescaled
%             y is the normalized magnitude of x
%             n is a two-element vector containing the offset and scale
%               if x is real, y = (x - n(1)) / n(2)
%
%   See also SCAL, AUTOSC, UNSCALE
%
%   06/2019, Keith Michel

%% Parse inputs
if ~nargin,      help(mfilename); return; end
if isempty(x),   return; end
if ~isreal(x),   x = abs(x); end

%% Normalize
n(1) = min(x(:));
x    = x - min(x(:));
n(2) = max(x(:));
x    = x / max(x(:));

end
