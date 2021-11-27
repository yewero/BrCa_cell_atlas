function WriteDWDVec (dirvec);

% Write data from matrix dirvec back to the textfile

% Output tab delimited text file

  fileName ='C:\DWD\DWDdata\DWD_Vec.txt';
  fid = fopen(fileName,'wt') ;

  if (fid<0)
      error ('Could not open the file DWD_Vec.txt for saving data');
  end;
% 'wt' is for "delete contents of this file and open 
% for writing" (with 't' for "text").
 
  trDirvec = transpose (dirvec);
  [n d]= size (trDirvec);


  for i = 1:n;
      for j = 1:(d-1) ;
          cntbytes = fprintf(fid,'%6.3f\t',trDirvec(i,j)) ;
      end ;
      cntbytes = fprintf(fid,'%6.3f\n',trDirvec(i,d)) ;
  end ;

  fclose(fid) ;