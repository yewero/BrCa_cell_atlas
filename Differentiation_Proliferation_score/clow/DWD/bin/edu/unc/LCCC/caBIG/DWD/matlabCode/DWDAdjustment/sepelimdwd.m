
   function [w,beta,residp,residn,alp,totalviolation,dualgap,flag]...
      = sep(Xp,Xn,penalty);

   % This subroutine is designed to calculate a linear discriminator
   % based on the training data Xp and Xn and the violation cost penalty.
   % It uses the DWD method of J.S. Marron and M.J. Todd,
   % "Distance Weighted Discrimination," TR 1339, 
   % School of Operations Research and Industrial Engineering,
   % Cornell University, Ithaca, New York, available at
   % ftp://ftp.orie.cornell.edu/pub/techreps/TR1339.pdf.

   % The user will need to get the SDPT3 optimization package,
   % available at http://www.math.nus.edu.sg/~mattohkc/sdpt3.html,
   % and install it. We recommend that sep.m be run from the SDPT3
   % directory.

   % SDPT3 includes an m-file startup.m, which sets the path,
   % default parameters, and global environment variables.
   % If the user has his/her own startup routine, this may not
   % perform the necessary tasks. Either the instructions
   % in SDPT3's startup.m should be appended, or SDPT3's
   % startup.m be renamed startupSDPT3.m, say, and this called
   % at the appropriate place.

   % The subroutine contains various commented-out statements, which
   % can be reinstated by those users wishing to see the time taken in
   % various steps, the numbers of variables and constraints, etc.

   % Input:  Matrices Xp and Xn, whose columns contain the training
   %            data for positive and negative instances respectively;
   %         scalar penalty, the cost for perturbing each residual.

   % Output: w and beta, such that x with w'x + beta positive (resp.,
   %            negative) is classified in the "+1" (resp. "-1") class.
   %         residp, the values of w'x + beta for x a column of Xp,
   %         residn, the values of w'x + beta for x a column of Xn.
   %         alp, the dual solution (usually, w = [Xp, - Xn]*alp
   %          normalized).
   %         totalviolation, the total amount added to the residuals.
   %         dualgap, the duality gap, a measure of the accuracy
   %            of the solutions found.
   %         flag, an indication of the success of the computation:
   %             0, success;
   %            -1, inaccurate solution;
   %            -2, problem infeasible or unbounded.

   flag = 0;

   % Find the dimensions of the data.

   [dp,np] = size(Xp);
   [dn,nn] = size(Xn);
   if (dn ~= dp), error('The dimensions are incompatible.'), end;
   d = dp;

   % Do the dimension reduction if in HDLSS setting.

   XpnY = [Xp, -Xn];
   XpnY11 = XpnY(1,1);
   n = np + nn; 
   if (d > n),
      [Q,RpnY] = qr(XpnY,0);
      dnew = n;
   else,
      RpnY = XpnY;
      dnew = d; 
   end;
   y = [ones(np,1); -ones(nn,1)];
   ym = y(2:n);

   % nv is the number of variables (eliminating beta), 
   % nc the number of constraints.

   nv = 1 + dnew + 3*n + n;
   nc = 2*n;

   % Set up the block structure, constraint matrix, rhs, and cost vector.

   blk = cell(2,2);
   blk{1,1} = 'q';
   blk{1,2} = [dnew+1, 3*ones(1, n)];
   blk{2,1} = 'l';
   blk{2,2} = n;
   blk;

   Avec = cell(2,1);
   A = zeros(nc,nv-n);
   col1 = RpnY(:,1);
   A(1:n-1,2:dnew+1) = (RpnY(:,2:n) - col1*ym')';
   A(1:n-1,dnew+5:3:dnew+1+3*n) = - speye(n-1,n-1);
   A(1:n-1,dnew+6:3:dnew+2+3*n) = + speye(n-1,n-1);
   A(1:n-1,dnew+2) = ym;
   A(1:n-1,dnew+3) = -ym;
   A(n,1) = 1;
   A(n+1:n+n,dnew+4:3:dnew+3+3*n) = speye(n,n);
   Avec{1,1} = A';
   Avec{2,1} = [[-ym,speye(n-1,n-1)];zeros(1+n,n)]';

   b = [zeros(n-1,1);ones(1+n,1)];

   C = cell(2,1);
   c = zeros(nv-n,1);
   c(dnew+2:3:dnew+1+3*n) = ones(n,1);
   c(dnew+3:3:dnew+2+3*n) = ones(n,1);
   C{1,1} = c;
   C{2,1} = penalty*ones(n,1);

   % Due to modification of parameters, variable OPTIONS is needed -- EZ
   % parameters ;
   [OPTIONS] = sqlparameters ;
   OPTIONS.maxit = 40;
   %disp ('sepelimdwd-sqlp in');
   [X0,lambda0,Z0] = infeaspt(blk,Avec,C,b);
   [obj,X,lambda,Z,info] = sqlp(blk,Avec,C,b,OPTIONS,X0,lambda0,Z0);
   %disp ('sepelimdwd-sqlp out');
   % If infeasible or unbounded, break.

   if (info.termcode > 0), flag = -2; return, end;

   % Compute the normal vector w and constant term beta.

   X1 = X{1}; X2 = X{2};
   barw = X1(2:dnew+1);
   if (d>n),
      w = Q*barw;
   else,
      w = barw;
   end;
   beta = X1(dnew+2) - X1(dnew+3) - X2(1) - col1'*barw;
   normw = norm(w);
   if normw < 1 - 1e-3, normw, end;
   normwm1 = 0;
   if normw > 1 - 1e-3,
      w = w / normw;
      normwm1 = norm(w)-1;
      beta = beta / normw;
   end;

   % Compute the minimum of the supposedly positive 
   % and the maximum of the supposedly negative residuals.
   % Refine the primal solution and print its objective value.

   residp = Xp'*w + beta*ones(np,1);
   residn = Xn'*w + beta*ones(nn,1);
   minresidp = min(residp);
   maxresidn = max(residn);
   res = XpnY'*w + beta*y;
   rsc = 1/sqrt(penalty);
   xi = rsc - res;
   xi = max(xi,0);
   totalviolation = sum(xi);
   minresidpmod = min(residp+xi(1:np)); 
   maxresidnmod = max(residn-xi(np+1:n));
   minxi = min(xi);
   maxxi = max(xi);
   resn = res + xi;
   rresn = 1./resn;
   format long
   primalobj = penalty*sum(xi) + sum(rresn);
   %devprimalobj = primalobj - obj(1);
   %if abs(devprimalobj) > 1e-3,
   %   devprimalobj,
   %   objxioff = penalty*(sum(xi)-sum(X2)),
   %end;

   % Compute the dual solution alp and print its objective value.

   alp = zeros(n,1);
   lambda1 = lambda(1:n-1);
   alp(1) = -ym'*lambda1;
   alp(2:n) = lambda1;
   alp = max(alp,0);
   sump = sum(alp(1:np));
   sumn = sum(alp(np+1:n));
   sum2 = (sump + sumn)/2;
   alp(1:np) = (sum2/sump)*alp(1:np);
   alp(np+1:n) = (sum2/sumn) * alp(np+1:n);
   maxalp = max(alp);
   if (maxalp > penalty | maxxi > 1e-3), 
      alp = (penalty/maxalp) * alp; 
   end;
   minalp = min(alp);
   p = RpnY*alp;
   eta = - norm(p);
   gamma = 2 * sqrt(alp);
   dualobj = eta + sum(gamma);
   %devdualobj = dualobj - obj(2);
   %if abs(devdualobj) > 1e-3,
   %   devdualobj,
   %end;
   format short

   % dualgap is the duality gap, a measure of the accuracy of the solution.

   dualgap = primalobj - dualobj;
   %if normw > 1 - 1e-3,
   %   wfromdual = -(XpnY*alp)/eta;
   %   normdifwprimwdual = norm(w-wfromdual);
   %   wd17 = wfromdual(1:7)';
   %end;
   if (dualgap > 1e-4),
      flag = -1;
   end;
   %arccoswe1 = w(1);
   %if (d>n),
   %   e1 = zeros(d,1); e1(1) = 1;
   %   pe1 = Q*(Q'*e1);
   %   npe1 = sqrt(pe1(1));
   %   arccoswpe1 = w'*pe1 / npe1,
   %   arccose1pe1 = npe1;
   %end;
   return
