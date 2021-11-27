function [meanAdjustType]= ReadDWDMeanAdjustType();
fid = fopen('C:\DWD\DWDdata\DWD_MeanAdjustType.txt','rt') ;
          %  'rt' is for "read only" and "text"

fstr = '%f' ;

[meanAdjustType,cnt] = fscanf(fid,fstr); 
          %  can add a "size" parameter when don't want to read all.
          %  cnt tells how many reads were done (2 * 3 fields = 6)         
fclose(fid) ;

