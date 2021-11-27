function BatchAdjustSM() 

% The following comment is not correct anymore for creating executable file
% They are retained here due to the further modification need ---Everett Zhou (EZ)

% function adjdata = BatchAdjustSM(rawdata,batchlabels,paramstruct) 

% BATCHADJUSTSM, of adjustment of "centerpoint" of subpopulations,
%   Steve Marron's matlab function
%     Intended for simple adjustment of subpopulation mean effects
%     in microarray data analysis.
%     Allows output of graphics, including "before" and "after"
%     PCA 2-d Draftsman's plots.
%     Algorithm uses DWD to find an effectrive direction for adjustment,
%     projects the data in this direction, and subtracts the 
%     subpopulation means.
%     This works pairwise, i.e. for a pair of subpopulations
%     For more subpopulations, apply this several times, to groups
%     of subpopulations.  Good choice of groups can be done by first
%     looking at the population structure, e.g. using a class colored
%     version of the PCA 2-d Draftsman's plot given by curvdatSM.
% Inputs:
%   rawdata     - d x n matrix of log gene expression data, 
%                          columns are cases, rows are genes 
%   batchlabels - 1 x n vector of vector of batch labels,
%                          for each case, must be +-1
%   paramstruct - a Matlab structure of input parameters
%                    Use: "help struct" and "help datatypes" to
%                         learn about these.
%                    Create one, using commands of the form:
%
%       paramstruct = struct('field1',values1, ...
%                            'field2',values2, ...
%                            'field3',values3) ;
%
%                          where any of the following can be used,
%                          these are optional, misspecified values
%                          revert to defaults
%
%    fields            values
%
%    viplot           vector of zeros and ones, indicating which plots to make
%                         1st entry:  1 (default) makes "before"
%                                           PCA 2-d Draftsman's plot
%                         2nd entry:  1 (default) makes "projection plot"
%                                           showing DWD performance
%                         3rd entry:  1 makes "afterwards projection plot",
%                                           of DWD applied to adjusted data
%                                           (default is 0, no plot) 
%                         4th entry:  1 (default) makes "after"
%                                           PCA 2-d Draftsman's plot
%                             (use zeros(4,1) for no plots)
%
%    savestr          string controlling saving of output,
%                         either a full path, or a file prefix to
%                         save in matlab's current directory
%                     unspecified:  results only appear on screen
%                     result:  add various plot names (depending 
%                                 on viplot) and add .ps
%
%    titlestr         String for Title of Projection plots 
%                           (will add to this depending on plot)
%                           (leave empty for no title at all)
%                           (default is "Batch Adjustment")
%
%    titlefontsize    Font Size for titles (uses Matlab default)
%                                   (18 is "fairly large")
%
%    legcellstr       Legend Cell String
%                     Use to apply labels to batches
%                     E.g.    legcellstr = {{'Batch 1' 'Batch 2'}} ;
%                         (this "cell within a cell" structure seems
%                          needed to pass in a cell array of strings,
%                          for Matlab 6.0, may be different for other
%                          Matlab versions)
%
%    iscreenwrite     0  (default)  no screen writes
%                     1  write to screen to show progress
%
%    minproj          left end of range in projection plots
%                               (default is output of axisSM)
%
%    maxproj          right end of range in projection plots
%                               (default is output of axisSM)
%
%    npc              Number of Principal Components to use in 
%                     draftsman's plot of 2-d PC projections
%                     (default = 6)
%
%
%1/18/05
%    imeantype        0   move both projected populations to 0
%                                  (sensible for cDNA and other 
%                                   differentially expressed data)
%                     +1  move -1 data set so projected mean is same as 
%                                    +1 data set
%                     -1  move +1 data set so projected mean is same as 
%                                    -1 data set
% Output:
%   adjdata     - d x n matrix of adjusted data, 
%                          columns are cases, rows are genes 
%    

%    Copyright (c) J. S. Marron 2003



%  First set paths to find needed subroutines
%

% Need to comment out the path --EZ

%addpath subroutines -end ;
%addpath subroutines\general -end ;
%addpath subroutines\smoothing -end ;



%  Set all input parameters to defaults
%
% These parameters are for plot
% viplot = [1; 1; 0; 1] ;
% savestr = [] ;
% titlestr = ['Batch Adjustment'] ;
% titlefontsize = [] ;
% legcellstr = {} ;
% iscreenwrite = 0 ;
% minproj = [] ;
% maxproj = [] ;
% npc = 6 ;

%  Now update parameters as specified,
%  by parameter structure (if it is used)
%
% if nargin > 2 ;   %  then paramstruct is an argument

 %  if isfield(paramstruct,'viplot') ;    %  then change to input value
 %    viplot = getfield(paramstruct,'viplot') ; 
 %    viplot                          % Display one of fields ---EZ
 %  end ;

 %  if isfield(paramstruct,'savestr') ;    %  then use input value
 %    savestr = getfield(paramstruct,'savestr') ; 
 %    if ~ischar(savestr) ;    %  then invalid input, so give warning
 %      disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
 %      disp('!!!   Warning from BatchAdjustSM.m:    !!!') ;
 %      disp('!!!   Invalid savestr,                 !!!') ;
 %      disp('!!!   using default of no save         !!!') ;
 %      disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
 %      savestr = [] ;
% end ;
%   end ;

%   if isfield(paramstruct,'titlestr') ;    %  then change to input value
%     titlestr = getfield(paramstruct,'titlestr') ; 
%   end ;

%   if isfield(paramstruct,'titlefontsize') ;    %  then change to input value
%     titlefontsize = getfield(paramstruct,'titlefontsize') ; 
%   end ;

%   if isfield(paramstruct,'legcellstr') ;    %  then change to input value
%     legcellstr = getfield(paramstruct,'legcellstr') ; 
%   end ;

%   if isfield(paramstruct,'iscreenwrite') ;    %  then change to input value
%     iscreenwrite = getfield(paramstruct,'iscreenwrite') ; 
%   end ;

%   if isfield(paramstruct,'minproj') ;    %  then change to input value
%     minproj = getfield(paramstruct,'minproj') ; 
%   end ;

%   if isfield(paramstruct,'maxproj') ;    %  then change to input value
%     maxproj = getfield(paramstruct,'maxproj') ; 
%   end ;

%   if isfield(paramstruct,'npc') ;    %  then change to input value
%     npc = getfield(paramstruct,'npc') ; 
%   end ;

%end ;    %  of resetting of input parameters



%  Set internal parameters
%
% npixshift = 20 ;
    %  number of pixels to shift new figures by

% The above comment is not correct anymore for creating executable file
% The OpenFile function will return two matries used for further processing
% --EZ
[rawdata, batchlabels] = ReadDWDInput;

d = size(rawdata,1);
n = size(rawdata,2);
errflag = logical(0) ;
if  (n ~= size(batchlabels,2))  | ...
    (1 ~= size(batchlabels,1))  ;

  errflag = logical(1) ;

  errstr = ['input "batchlabels" must be a row vector ' ...
                     'of length ' num2str(n)] ;
  %disp('Here-1');
elseif  (sum(batchlabels == 1) + sum(batchlabels == -1))  ~=  n  ;

  errflag = logical(1) ;

  errstr = 'entries in "batchlabels" must all be +- 1' ;
  %disp('Here-2');
end ;

% if ~iscell(legcellstr) ;
%   disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
%   disp('!!!   Warning from BatchAdjustSM:          !!!') ;
%   disp('!!!   legcellstr is not a cell string      !!!') ;
%   disp('!!!   Will turn off the Legend Cell String !!!') ;
%   disp('!!!   i.e. the text with batch labels      !!!') ;
%   disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
%   legcellstr = {} ;

% else ;

%   if length(legcellstr) > 2 ;
%     disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
%     disp('!!!   Warning from BatchAdjustSM:       !!!') ;
%     disp('!!!   legcellstr is too big,            !!!') ;
%     disp('!!!   Will use only first two entries   !!!') ;
%     disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
%     legcellstr = legcellstr([1,2]) ;

%   end ;

% end ;



if errflag ;    %  then had a fatal input error

  disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
  disp('!!!   Error from BatchAdjustSM.m:                   !!!') ;
  disp(['!!!   ' errstr]) ;
  disp('!!!   Terminating execution, with an empty return   !!!') ;
  disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;

  adjdata = [] ;


else ;    %  inputs OK, so do serious work


 %  if length(viplot) ~= 4 ;

 %    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
 %    disp('!!!   Warning from BatchAdjustSM.m:   !!!') ;
 %    disp('!!!   Invalid size of viplot,         !!!') ;
 %    disp('!!!   reverting to default            !!!') ;
 %    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;

 %    viplot = [1; 1; 0; 1] ;
  
 %  end ;

 %  viplot = logical(viplot) ;
      %  make sure this is logical, 
      %  for use in if statements below



  %  Do DWD Batch adjustment
  %
  
  flagp = (batchlabels == 1) ;
  flagn = (batchlabels == -1) ;


  %  Find DWD direction
  
  %disp('Here-3');
  dirvec = DWD1SM(rawdata(:,flagp),rawdata(:,flagn)) ;
 

  %disp('Here-4');
  
    %  Project data
  %
  vprojp = rawdata(:,flagp)' * dirvec ;
  vprojn = rawdata(:,flagn)' * dirvec ;

  meanprojp = mean(vprojp) ;
  meanprojn = mean(vprojn) ;


  %  Subtract respective class means
  %
  %1/18/05
%    imeantype        0   move both projected populations to 0
%                                  (sensible for cDNA and other 
%                                   differentially expressed data)
%                     +1  move -1 data set so projected mean is same as 
%                                    +1 data set
%                     -1  move +1 data set so projected mean is same as 
%                                    -1 data set
%  adjdata(:,flagp) = rawdata(:,flagp) - vec2matSM(meanprojp * dirvec,length(vprojp)) ;
%  adjdata(:,flagn) = rawdata(:,flagn) - vec2matSM(meanprojn * dirvec,length(vprojn)) ;

imeantype = ReadDWDMeanAdjustType;
%imeantype=-1
if imeantype == -1 ;    %  move +1 data, so projected mean is same as -1 data  
    disp ('Mean Adjust Type -1');
     % change length(vprojn) to length(vprojp)
     %adjdata(:,flagp) = rawdata(:,flagp) - vec2matSM(meanprojp * dirvec,length(vprojp)) ...
     %                                   + vec2matSM(meanprojn * dirvec,length(vprojn)) ;
     adjdata(:,flagp) = rawdata(:,flagp) - vec2matSM(meanprojp * dirvec,length(vprojp)) ...
                                        + vec2matSM(meanprojn * dirvec,length(vprojp)) ;
     adjdata(:,flagn) = rawdata(:,flagn) ;
elseif imeantype == 1 ;    %  move -1 data, so projected mean is same as +1 data
    disp ('Mean Adjust Type 1');
     adjdata(:,flagp) = rawdata(:,flagp) ;
     % change length(vprojp) to length(vprojn)
     %adjdata(:,flagn) = (rawdata(:,flagn) - vec2matSM(meanprojn * dirvec,length(vprojn)) )...
     %                                 + vec2matSM(meanprojp * dirvec,length(vprojp)) ;     
     adjdata(:,flagn) = (rawdata(:,flagn) - vec2matSM(meanprojn * dirvec,length(vprojn)) )...
                                      + vec2matSM(meanprojp * dirvec,length(vprojn)) ;
elseif imeantype == 0;    %  move both projected populations to 0
    disp ('Mean Adjust Type 0');
    adjdata(:,flagp) = rawdata(:,flagp) - vec2matSM(meanprojp * dirvec,length(vprojp)) ;
    adjdata(:,flagn) = rawdata(:,flagn) - vec2matSM(meanprojn * dirvec,length(vprojn)) ;
    % [n1 d1]= size (rawdata(:,flagn))
    % [n2 d2]= size (vec2matSM(meanprojn * dirvec,length(vprojn)))
    % [n3 d3]= size (vec2matSM(meanprojp * dirvec,length(vprojp)))

end ;
  
  %adjdata
  
  %dirvec1 = dirvec(1:30)
  %dirvec2 = dirvec(3974:end)  
  
  % Save dirvec back to a file
  WriteDWDVec(dirvec');
  % The result of adjdata will be written into a text file by WriteFile
  WriteDWDOutput (adjdata, batchlabels);
end;





