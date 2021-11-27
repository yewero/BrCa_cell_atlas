function [mdata]= ReadKDEInput();
% This file is going to open a middle step file saved from Java graph
% It will return one matrix for the further manipulation

%disp('Open a file in this function') ;


  
%fileName = 'C:\DWD\DWDdata\DWD_KDE.txt';
fileName = 'C:\DWD\DWDdata\DWD_KDE_Input.txt';
% To find the columns of each line
%fid = fopen ('Combined_cut.txt','r');
fid = fopen (fileName,'r');
if (fid<0)
    error ('Could not open the file DWD_KDE.txt');
end;

s=fgetl (fid);

col=1;
len = length(s);
for i=1:len
    ch =s(i);
    if ch =='	';
        col=col+1;
    end;
end; 
%col

fclose (fid);

% To find the rows of data file
%fid = fopen ('Combined_cut.txt','r');
fid = fopen (fileName,'r');
if (fid<0)
    error ('Could not open the file DWD_Input.txt');
end;
row=0;
s=fgetl (fid);
while (ischar(s))
    row = row+1;
    s=fgetl (fid);
end;
%row

fclose (fid);

% To read data into an array
%fid = fopen ('Combined_cut.txt','r');
fid = fopen (fileName,'r');
if (fid<0)
    error ('Could not open the file DWD_Input.txt');
end;

text = fscanf (fid, '%f', [col, row]);
fclose (fid);


mdata = text';

%  try 
%      fid = fopen ('Combined_cut.txt','r');
%      text = fscanf (fid, '%f', [18, 8]);
%      fclose (fid);
%      textin = text';
      %[r c] = size(textin)
      %textin = dlmread(fileName,'\t');
      %  catch
%    errmsg = lasterr;
%    disp('** ERROR: The input file was not properly opened.');
%    disp('This was probably because the input file was in the wrong format.');
%    disp('Make sure the input file is formatted as Microsoft Excel 5.0.');
%    disp(errmsg);
%  end
  %textin
  
  %mdata = textin(2:end,1:end);
  %disp(' ') ;
  %disp('  This should give the 1st data value:') ;
  %mdata(1,1)
  %disp(' ') ;
  %disp('  This should give the last data value:') ;
  %mdata(end,end)
  %vGeneName = textin(2:end,1);

  %vBatch = textin(1,1:end);
  %disp(' ') ;
  %disp('  This should give the 1st batch value:') ;
  %vBatch(1)
  %disp(' ') ;
  %disp('  This should give the last batch value:') ;
  %vBatch(end)
  
% vSource has values 1 or 2,
% this maps these to -1 or + 1
% Want to display batchlabels and legcellstr -- EZ
  %batchlabels = 2 * vBatch - 3; 
  %disp ('batchlabels ');
  % fclose (fid);
  