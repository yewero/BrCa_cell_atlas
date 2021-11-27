function WriteKDEOutput (xgridkde);

% Write data from matrix adjdtat back to the textfile

%dlmwrite ('DWD_cut_out_ez.txt', adjdata, 'delimiter', '\t', 'precision', '%6.3f');
%dlmwrite ('DWD_cut_out_ez.txt', adjdata, 'delimiter', '\t');

% Output tab delimited text file

  fileName = 'C:\DWD\DWDdata\DWD_KDE_Output.txt';
  fid = fopen(fileName,'wt') ;

  if (fid<0)
      error ('Could not open the file DWD_Output.txt for saving data');
  end;
% 'wt' is for "delete contents of this file and open 
% for writing" (with 't' for "text").

  %n = length(vBatch);
  %d = length(batchlabels);
  [n d]= size (xgridkde);
  
  for i = 1:n;
      for j = 1:(d-1) ;
          cntbytes = fprintf(fid,'%6.4f\t',xgridkde(i,j)) ;
      end ;
      cntbytes = fprintf(fid,'%6.4f\n',xgridkde(i,d)) ;
  end ;
  
  fclose(fid) ;