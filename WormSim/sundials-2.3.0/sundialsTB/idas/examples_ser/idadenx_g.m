function [g, flag, new_data] = idadenx_g(t,y,yp,data)
%IDADENX_G - Root-finding function for the IDADENX example problem.
%
%   See also: idadenx, IDARootFn 

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date: 2006/07/17 16:49:50 $


g(1) = y(1) - 0.0001;
g(2) = y(3) - 0.01;

flag = 0;
new_data = [];
