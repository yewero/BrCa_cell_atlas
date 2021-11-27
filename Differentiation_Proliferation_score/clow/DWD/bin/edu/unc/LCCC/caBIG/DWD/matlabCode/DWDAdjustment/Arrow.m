%%********************************************************
%% Arrow: 
%%  
%% Fx = Arrow(pblk,f,x,options); 
%%
%% if options == 0; 
%%    Fx = Arr(F)*x
%% if options == 1; 
%%    Fx = Arr(F)^{-1}*x 
%%
%% SDPT3: version 3.1 
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 15 Sep 2004
%%********************************************************

   function   Fx = Arrow(pblk,f,x,options); 

   if nargin == 3; options = 0; end;

   s = 1 + [0, cumsum(pblk{2})]; 
   idx1 = s(1:length(pblk{2})); 
   if options == 0
      inprod = mexqops(pblk,f,x,1);  
      Fx  = mexqops(pblk,f(idx1),x,3) + mexqops(pblk,x(idx1),f,3); 
      Fx(idx1) = inprod; 
   else
      gamf2 = mexqops(pblk,f,f,2);
      gamprod = mexqops(pblk,f,x,2);
      alpha = gamprod./gamf2; 
      Fx = mexqops(pblk,1./f(idx1),x,3) - mexqops(pblk,alpha./f(idx1),f,3); 
      Fx(idx1) = alpha;
   end
%%
%%********************************************************






