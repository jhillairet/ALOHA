function [h]=matrix_pcolor(gp);
% Author: O.Meneghini
A=gp;
A=[A,A(:,end)];
A=[A;A(end,:)];
h=pcolor(abs(A));
view(0,-90)
%set(h,'Edgecolor','none','FaceColor','interp') 