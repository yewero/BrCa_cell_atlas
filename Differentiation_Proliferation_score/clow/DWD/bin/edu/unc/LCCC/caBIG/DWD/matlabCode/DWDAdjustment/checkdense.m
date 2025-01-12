%%********************************************************************
%% checkdense : identify the dense columns of a matrix
%% 
%% SDPT3: version 3.1 
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 15 Sep 2004
%%********************************************************************

  function  idxden = checkdense(A);

   m = size(A,1);
   idxden = []; 
   nzratio = 1;
   if (m > 1000); nzratio = 0.20; end;
   if (m > 2000); nzratio = 0.10; end;
   if (m > 5000); nzratio = 0.05; end;
   if (nzratio < 1)
      nzcolA = sum(spones(A)); 
      idxden = find(nzcolA > nzratio*m);
   end
%%********************************************************************
