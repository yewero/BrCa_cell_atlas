function [kde,xgrid] = kdeSMEZ(data) 
% Modified by Everett Zhou based on Steve Marron's function [kde,xgrid,mker] = kdeSM(data,paramstruct) 
% 3/31/2005

% KDESM, Kernel Density Estimate (1-d, Gaussian Kernel)
%   Steve Marron's matlab function
%     Does 1-d kernel density estimation, using binned (default) or 
%     direct (either matrix, or loops for bigger data sets), 
%     implementations, with the bandwidth either user specified 
%     (can be vector), or data driven (SJPI, Normal Reference, 
%     Silverman's ROT2, Oversmoothed).
%
% Inputs:
%   data        - either n x 1 column vector of 1-d data
%                     or vector of bincts, when imptyp = -1
%
%   paramstruct - a Matlab structure of input parameters
%                    Use: "help struct" and "help datatypes" to
%                         learn about these.
%                    Create one, using commands of the form:
%
%       paramstruct = struct('field1',values1,...
%                            'field2',values2,...
%                            'field3',values3) ;
%
%                          where any of the following can be used,
%                          these are optional, misspecified values
%                          revert to defaults
%
%    fields            values
%
%    vh               vector of bandwidths, or specifies data driven:
%                       0 (or not specified)  -  Sheather Jones Plug In
%                       -1  -  Simple Normal Reference
%                       -2  -  Silverman's Rule Of Thumb 2 
%                                 (20% Smaller than min of sd and IQR)
%                       -3  -  Oversmoothed
%                           Note: <0 only works for imptype = 0
%                       >0  -  Use input number (numbers if vector)
%                           Note: this MUST be >0 for imptyp >= 1
%                                      (the direct implementations)
%
%    vxgrid           vector of parameters for, or values of, grid to evaluate at:
%                       0 (or not specified)  -  use endpts of data and 401 bins
%                       [le; lr]  -  le is left end, re is right, 401 bins
%                              (get error message and no return if le > lr)
%                       [le; lr; nb] - le left, re right, and nb bins
%                       xgrid  -  Use input values 
%                                 Note:  need to have more than 3 entries,
%                                      and only works when imptyp = 1 or 2
%
%    imptyp           flag indicating implementation type:
%                      -1  -  binned version, and "data" is assumed to be
%                                        bincounts of prebinned data
%                       0 (or not specified)  -  linear binned version
%                                        and bin data here
%                       1  -  Direct matrix implementation
%                       2  -  Slow looped implementation (only useful when
%                                 1 creates matrices that are too large)
%
%    eptflag          endpoint truncation flag (only has effect when imptyp = 0):
%                       0 (or not specified)  -  move data outside range to
%                                        nearest endpoint
%                       1  -  truncate data outside range
%
%    ibdryadj         index of boundary adjustment
%                       0 (or not specified)  -  no adjustment
%                       1  -  mirror image adjustment
%                       2  -  circular design
%
%    ndataoverlay     overlay raw data as a jitterplot  (requires  imptyp >= 0)
%                                    (only has effect when plot is made here)
%                       0  -  (default) no data plot
%                       1  -  overlay up to 1000 points (random choice, when more)
%                       2  -  overlay full data set
%                       n > 2   -  overlay n random points
%
%    dolcolor         data overlay color
%                                    (only has effect when plot is made here)
%                           default is 'g'
%
%    ibigdot          0  (default)  use Matlab default for dot sizes
%                     1  force large dot size in prints (useful since some
%                              postscript graphics leave dots too small)
%                              (Caution: shows up as small in Matlab view)
%
%    dolhtseed        seed for random heights used in data overlay jitter plot
%                           default is [] (for using current Matlab seed) 
%                                    (should be an integer with <= 8 digits)
%
%    linewidth        width of lines (only has effect when plot is made here)
%                           default is 2, for length(vh) <= 3,
%                           default is 0.5, for length(vh) > 3,
%
%    linecolor        string with color for lines
%                                    (only has effect when plot is made here)
%                           default is 'b'
%                           use the empty string, '', for standard Matlab colors
%
%    titlestr         string with title (only has effect when plot is made here)
%                           default is empty string, '', for no title
%
%    titlefontsize    font size for title
%                                    (only has effect when plot is made here,
%                                     and when the titlestr is nonempty)
%                           default is empty [], for Matlab default
%
%    xlabelstr        string with x axis label
%                                    (only has effect when plot is made here)
%                           default is empty string, '', for no xlabel
%
%    ylabelstr        string with y axis label
%                                    (only has effect when plot is made here)
%                           default is empty string, '', for no ylabel
%
%    labelfontsize    font size for axis labels
%                                    (only has effect when plot is made here,
%                                     and when a label str is nonempty)
%                           default is empty [], for Matlab default
%
%    plotbottom      bottom of plot window,
%                                    use [] to get 0 - .05 * range (default)
%
%    plottop         top of plot window,
%                                    use [] to get max + .05 * range (default)
%
%    iplot            1  -  plot even when there is numerical output
%                     0  -  (default) only plot when no numerical output
%
% Outputs:
%     (none)  -  Draws a graph of the result (in the current axes)
%     kde     -  col vector of heights of kernel kernel density estimate,
%                    unless vh is a vector (then have matrix, with
%                    corresponding cols as density estimates)
%     xgrid   -  col vector grid of points at which estimate(s) are 
%                    evaluted (useful for plotting), unless grid is input,
%                    can also get this from linspace(le,re,nb)'  
%     mker    -  matrix (vector) of kernel functions, evaluated at xgrid,
%                    which can be plotted to show "effective window 
%                    sizes" (currently scaled to have mass 0.05, may
%                    need some vertical rescaling)
%
% Assumes path can find personal functions:
%    vec2matSM.m
%    lbinrSM.m
%    bwsjpiSM.m
%    bwosSM.m
%    bwrotSM.m
%    bwsnrSM.m

%    Copyright (c) J. S. Marron 1996-2001

%  First set all parameters to defaults
vh = 0 ;      %  use default SJPI
vxgrid = 0 ;
imptyp = 0 ;
eptflag = 0 ;
ibdryadj = 0 ;
ndataoverlay = 0 ;
dolcolor = 'g' ;
ibigdot = 0 ;
dolhtseed = [] ;
linewidth = 2 ;
linecolor = 'b' ;
titlestr = '' ;
titlefontsize = [] ;
xlabelstr = '' ;
ylabelstr = '' ;
labelfontsize = [] ;
plotbottom = [] ;
plottop = [] ;
iplot = 0 ;


%  Now update parameters as specified,
%  by parameter structure (if it is used)
%
% if nargin > 1 ;   %  then paramstruct has been added
% 
%   if isfield(paramstruct,'vh') ;    %  then change to input value
%     vh = getfield(paramstruct,'vh') ; 
%   end ;
% 
%   if isfield(paramstruct,'vxgrid') ;    %  then change to input value
%     vxgrid = getfield(paramstruct,'vxgrid') ; 
%   end ;
% 
%   if isfield(paramstruct,'imptyp') ;    %  then change to input value
%     imptyp = getfield(paramstruct,'imptyp') ; 
%   end ;
% 
%   if isfield(paramstruct,'eptflag') ;    %  then change to input value
%     eptflag = getfield(paramstruct,'eptflag') ; 
%   end ;
% 
%   if isfield(paramstruct,'ibdryadj') ;    %  then change to input value
%     ibdryadj = getfield(paramstruct,'ibdryadj') ; 
%   end ;
% 
%   if isfield(paramstruct,'ndataoverlay') ;    %  then change to input value
%     ndataoverlay = getfield(paramstruct,'ndataoverlay') ; 
%   end ;
% 
%   if isfield(paramstruct,'dolcolor') ;    %  then change to input value
%     dolcolor = getfield(paramstruct,'dolcolor') ; 
%   end ;
% 
%   if isfield(paramstruct,'ibigdot') ;    %  then change to input value
%     ibigdot = getfield(paramstruct,'ibigdot') ; 
%   end ;
% 
%   if isfield(paramstruct,'dolhtseed') ;    %  then change to input value
%     dolhtseed = getfield(paramstruct,'dolhtseed') ; 
%   end ;
% 
%   if isfield(paramstruct,'linewidth') ;    %  then change to input value
%     linewidth = getfield(paramstruct,'linewidth') ; 
%   end ;
% 
%   if isfield(paramstruct,'linecolor') ;    %  then change to input value
%     linecolor = getfield(paramstruct,'linecolor') ; 
%   end ;
% 
%   if isfield(paramstruct,'titlestr') ;    %  then change to input value
%     titlestr = getfield(paramstruct,'titlestr') ; 
%   end ;
% 
%   if isfield(paramstruct,'titlefontsize') ;    %  then change to input value
%     titlefontsize = getfield(paramstruct,'titlefontsize') ; 
%   end ;
% 
%   if isfield(paramstruct,'xlabelstr') ;    %  then change to input value
%     xlabelstr = getfield(paramstruct,'xlabelstr') ; 
%   end ;
% 
%   if isfield(paramstruct,'ylabelstr') ;    %  then change to input value
%     ylabelstr = getfield(paramstruct,'ylabelstr') ; 
%   end ;
% 
%   if isfield(paramstruct,'labelfontsize') ;    %  then change to input value
%     labelfontsize = getfield(paramstruct,'labelfontsize') ; 
%   end ;
% 
%   if isfield(paramstruct,'plotbottom') ;    %  then change to input value
%     plotbottom = getfield(paramstruct,'plotbottom') ; 
%   end ;
% 
%   if isfield(paramstruct,'plottop') ;    %  then change to input value
%     plottop = getfield(paramstruct,'plottop') ; 
%   end ;
% 
%   if isfield(paramstruct,'iplot') ;    %  then change to input value
%     iplot = getfield(paramstruct,'iplot') ; 
%   end ;
% 
% 
% end ;  %  of resetting of input parameters



%  Calculate kde
%
if imptyp > 0 ;    %  Then do direct implementation

  if min(vh) > 0 ;    %  Then have valid bandwidths, so proceed

    n = length(data) ;

    if length(vxgrid) > 3 ;  %  Then use input grid
      xgrid = vxgrid ;
      nbin = length(xgrid) ;
      lend = min(xgrid) ;
      rend = max(xgrid) ;
    else ;                    %  Need to generate a grid
      nbin = 401 ;         %  Default
      lend = min(data) ;   %  Default
      rend = max(data) ;   %  Default
      if length(vxgrid) >= 2 ;      %  use input endpoints
        lend = vxgrid(1) ;
        rend = vxgrid(2) ;
      end ;
      if length(vxgrid) == 3 ;      %  use number of grid points
        nbin = vxgrid(3) ;
      end ;

      if lend > rend ;    %  Then bad range has been input
        disp('!!!   Error in kdeSM: invalid range chosen  !!!') ;
        xgrid = [] ;
      else ;
        xgrid = linspace(lend,rend,nbin)' ;
      end ;
    end ;


    %  do boundary adjustment if needed
    %
    if ibdryadj == 1 ;    %  then do mirror image adjustment
      badata = [(lend - (data - lend)); data; (rend + (rend - data))] ;
    elseif ibdryadj == 2 ;    %  then do circular design adjustment
      badata = [(data - (rend - lend)); data; (rend - lend + data)] ;
    else ;
      badata = data ;
    end ;


    %  Loop through bandwidths
    kde = [] ;
    for ih = 1:length(vh) ;
      h = vh(ih) ;

      if imptyp ~= 2 ;  %  Then do direct matrix implementation
        kdeh = vec2matSM((badata ./ h),nbin) - vec2matSM((xgrid' ./ h),n) ;
          %  efficient way to divide all dif's by h
          %  variable name "kde" is used to avoid creating too many biggies
        kdeh = exp(-(kdeh .^2) / 2) ;
          %  exponential part of Gaussian density
        kdeh = sum(kdeh)' ;
          %  sum part of kde, and make result a column vector
        kdeh = kdeh / (n * h * sqrt(2 * pi)) ;
          %  normalize, and mult by Gaussain density constant
        kde = [kde kdeh] ;
      else ;   %  Do slower looped implementation
        kdeh = [] ;
        for ixg = 1:nbin ;    %  Loop through grid points
          kdehx = (badata - xgrid(ixg)) / h ;
          kdehx = sum(exp(-(kdehx .^2) / 2)) ;
          kdeh = [kdeh; kdehx] ;
        end ;
        kdeh = kdeh / (n * h * sqrt(2 * pi)) ;
        kde = [kde kdeh] ;
      end ;
    end ;

  else ;    %  Have invalid bandwidths

    disp('!!!   Error in kdeSM: A bandwidth is invalid   !!!') ;
    disp('    (Note: cannot use data driven, with direct impl''s)') ;

  end ;

else ;     %  Then do binned implementation

  if imptyp == -1 ;   %  Then data have already been binned

    if (length(vxgrid) == 1) | (length(vxgrid) > 3) ;
                         %  Then can't proceed because don't have bin ends
      disp('!!!   Error: kdeSM needs to know the endpoints   !!!') ;
      disp('!!!            to use this implementation        !!!') ;
      bincts = [] ;
    else ;
      bincts = data ;

      nbin = 401 ;
      lend = vxgrid(1) ;
      rend = vxgrid(2) ;
      if length(vxgrid) == 3 ;          %  then use number of grid points
        nbin = vxgrid(3) ;
      end ;

      if nbin ~= length(bincts) ;    %  Then something is wrong
        disp('!!!   Warning: kdeSM was told the wrong number of bins   !!!') ;
        disp('!!!            will just use the number of counts.       !!!') ;
        nbin = size(bincts,1) ;
      end ;
    end ;

  else ;               %  Then need to bin data

    if length(vxgrid) > 3 ;  %  Then need to warn of change to default
      disp('!!!   Warning: kdeSM was given an xgrid, and also   !!!') ;
      disp('!!!       asked to bin; will bin and ignore xgrid   !!!') ;
    end ;

    %  Specify grid parameters
    nbin = 401 ;         %  Default
    lend = min(data) ;   %  Default
    rend = max(data) ;   %  Default
    if (length(vxgrid) == 2) | (length(vxgrid) == 3) ;
                                     %  then use input endpoints
      lend = vxgrid(1) ;
      rend = vxgrid(2) ;
    end ;
    if length(vxgrid) == 3 ;          %  then use number of grid points
      nbin = vxgrid(3) ;
    end ;

    if lend > rend ;    %  Then bad range has been input
      disp('!!!   Error in kdeSM: invalid range chosen  !!!') ;
      bincts = [] ;
    else ;
      bincts = lbinrSM(data,[lend,rend,nbin],eptflag) ;
    end ;

    %  Can do data-based bandwidth selection here, if specified
    if vh == -1 ;        %  Then use Simple Normal Reference
      vh = bwsnrSM(data) ;
    elseif vh == -2 ;    %  Then use Silverman's Rule Of Thumb 2 
                          %  (~10% Smaller than min of sd and IQR)
      vh = bwrotSM(data) ;
    elseif vh == -3 ;    %  Then use Terrell's Oversmoother
      vh = bwosSM(data) ;
    elseif min(vh) <= 0 ;     %  Then be sure to use default SJPI 
                          %    (in case an unsupported value was input)
      vh = 0 ;
    end ;

  end ;
  n = round(sum(bincts)) ;
          %  put this here in case of truncations during binning

  %  do boundary adjustment if needed
  %
  if ibdryadj == 1 ;    %  then do mirror image adjustment
    babincts = [flipud(bincts); bincts; flipud(bincts)] ;
  elseif ibdryadj == 2 ;    %  then do circular design adjustment
    babincts = [bincts; bincts; bincts] ;
  else ;
    babincts = bincts ;
  end ;


  %  Get bandwidth (if still not yet specified)
  if vh == 0 ;    %  Then use SJPI bandwidth
      vh = bwsjpiSM(bincts,[lend; rend; nbin],0,-1) ;
  end ;


  %  Loop through bandwidths
  kde = [] ;
  for ih = 1:length(vh) ;
    h = vh(ih) ;

    %  Create vector of kernel values, at equally spaced grid
    delta = (rend - lend) / (nbin - 1) ;    %  binwidth
    k = nbin - 1 ;    %  index of last nonzero entry of kernel vector
    arg = linspace(0,k * delta / h,k + 1)' ;
    kvec = exp(-(arg.^2) / 2) / sqrt(2 * pi) ;
    kvec = [flipud(kvec(2:k+1)); kvec] ;

    %  Do actual kernel density estimation
    kdeh = conv(babincts,kvec) ;

    if  ibdryadj == 1  |  ibdryadj == 2 ;    %  then did boundary adjustment
      kdeh = kdeh(nbin+k+1:k+2*nbin) / (n * h) ;
    else ;
      kdeh = kdeh(k+1:k+nbin) / (n * h) ;
    end ;

    if h < 3 * delta ;    %  Then need to normalize
                             %  to make numerical integral roughly 1
      kdeh = kdeh / (sum(kdeh) * delta) ;
    end ;

    kde = [kde kdeh] ;
  end ;

  xgrid = linspace(lend,rend,nbin)' ;

end ;


%  Create matrix of kernels, if this is needed
%
if nargout == 3 ;
  cent = mean([lend; rend]) ;
          %  centerpoint of evaluation grid
  if length(vh) > 1 ;
    mih = vec2matSM(1 ./ vh',nbin) ;
    mker = vec2matSM((xgrid - cent),length(vh)) .* mih;
          %  argument of gaussian kernel
  else ;
    mih = 1 / vh ;
    mker = (xgrid - cent) .* mih;
          %  argument of gaussian kernel
  end ;
  mker = exp(-mker.^2 / 2) .* mih / sqrt(2 * pi) ;
          %  Gaussian kernels with mass 1
  mker = 0.05 * mker ;
          %  Make masses = 0.05
end ;


%  Make plots if no numerical output requested, or if plot requested
%
% if  nargout == 0  | ...
%       iplot == 1  ;  %  Then make a plot
% 
% 
%   if  length(vh) > 3  &  ~isfield(paramstruct,'linewidth')  ;
%                               %  then need to change default value of linewidth
%     linewidth = 0.5 ;
%   end ;
% 
% 
%   if isempty(linecolor) ;
%     plot(xgrid,kde,'LineWidth',linewidth) ;
%   else ;
%     plot(xgrid,kde,'LineWidth',linewidth,'Color',linecolor) ;
%   end ;
% 
%   if  isempty(plottop)  &  isempty(plotbottom)  ;    %  then adjust top and bottom
%     plotbottom = 0 ;
%     plottop = max(max(kde)) ;
%     plotrange = plottop - plotbottom ;
%     plotbottom = plotbottom - 0.05 * plotrange ;
%     plottop = plottop + 0.05 * plotrange ;
%   elseif isempty(plottop) ;                          %  then only adjust top
%     plottop = max(max(kde)) ;
%     plotrange = plottop - plotbottom ;
%     plottop = plottop + 0.05 * plotrange ;
%   elseif isempty(plotbottom) ;                       %  then only adjust bottom
%     plotbottom = 0 ;
%     plotrange = plottop - plotbottom ;
%     plotbottom = plotbottom - 0.05 * plotrange ;
%   end ;
% 
%   vax = [lend,rend,plotbottom,plottop] ;
%   axis(vax) ;
% 
%   if  ndataoverlay > 0  &  imptyp >= 0  ;
%                             %  then overlay jitterplot of data
%   
%     if ~isempty(dolhtseed) ;
%       rand('seed',dolhtseed) ;
%     end ;
% 
%     if ndataoverlay == 1 ;
%       ndo = min(n,1000) ;
%     elseif ndataoverlay == 2 ;
%       ndo = n ;
%     else ;
%       ndo = min(n,ndataoverlay) ;
%     end ;
% 
%     flagleft = (data < lend) ;
%         %  ones where data below left end
%     flagright = (data > rend) ;
%         %  ones where data above right end
%     nleft = sum(flagleft) ;
%     nright = sum(flagright) ;
%     if nleft + nright > 0 ;    %  then need to deal with points outside range
% 
%       if eptflag == 1 ;    %  the truncate data outside range
%         datatrunc = data(~(flagleft | flagright)) ;
%             %  keep data that is not (outside left or outside right)
% 
%       else ;    %  then move outside points to nearest end
%         datatrunc = data ;
%         if nleft > 0 ;    %  then replace those points with lend
%           datatrunc(flagleft) = lend * ones(nleft,1) ;
%         end ;
%         if nright > 0 ;    %  then replace those points with rend
%           datatrunc(flagright) = rend * ones(nright,1) ;
%         end ;
% 
%       end ;
% 
%     else ;
% 
%       datatrunc = data ;
% 
%     end ;
% 
% 
%     if ndo < n ;    %  then need to subsample
%       [temp,randperm] = sort(rand(n,1)) ;
%             %  randperm is a random permutation of 1,2,...,n
%       vkeep = randperm(1:ndo) ;
%             %  indices of elements to keep for display
%       dataol = datatrunc(vkeep,:) ;
%     else ;    %  overlay full data set
%       dataol = datatrunc ;
%     end ;
% 
%     %  overlay selected data
%     %
%     vax = axis ;
%     plotbottom = vax(3) ;
%     ftop = vax(4) ;
%     hold on ;
%       yrand = plotbottom + (0.7 + 0.2 * rand(ndo,1)) ...
%                                             * (plottop - plotbottom) ;
%              %  y coords for jittering
% 
%       if ibigdot == 1 ;   %  plot deliberately large dots
%         plot(dataol,yrand,[dolcolor 'o'],'MarkerSize',1,'LineWidth',2) ;
%       else ;    %  use matlab default dots
%         plot(dataol,yrand,[dolcolor '.']) ;
%       end ;
% 
%     hold off ;
% 
% 
%     %  reorder, so that overlaid data appears under the smooth
%     %
%     vachil = get(gca,'Children') ;
%     nchil = length(vachil) ;
%     vachil = [vachil(2:nchil); vachil(1)] ;
%     set(gca,'Children',vachil) ;
%     
% 
%   end ;
% 
% 
%   if ~isempty(titlestr) ;
%     if isempty(titlefontsize) ;
%       title(titlestr) ;
%     else ;
%       title(titlestr,'FontSize',titlefontsize) ;
%     end ;
%   end ;
% 
% 
%   if ~isempty(xlabelstr) ;
%     if isempty(labelfontsize) ;
%       xlabel(xlabelstr) ;
%     else ;
%       xlabel(xlabelstr,'FontSize',labelfontsize) ;
%     end ;
%   end ;
% 
% 
%   if ~isempty(ylabelstr) ;
%     if isempty(labelfontsize) ;
%       ylabel(ylabelstr) ;
%     else ;
%       ylabel(ylabelstr,'FontSize',labelfontsize) ;
%     end ;
%   end ;
% 
% 
% end ;

