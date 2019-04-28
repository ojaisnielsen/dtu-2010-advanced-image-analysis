function v=MarkShape(go);
go =1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% shape=MarkShape
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Usage:
%
% Make points with the left mouse button in an activ figure.
% Undo last point with  'u'.
% Stop by pressing 'e'.
%
% Note: When plotting using 'DrawShape' the shape will 
% automatically be closed by connecting the first and the 
% last point marked.
%
% See: DrawShape

stopinp=0;
v=[];
disp('left-click : Mark a snake point');
disp('   ''u''     : Undo last point')
disp('   ''e''     : End');
while stopinp==0,
   
[vc,vr,but]=ginput(1);
if but==1
   v=[v; vr vc];
	line( v(:,2), v(:,1), 'Marker', '+','Color','red', 'LineWidth', 2 );
elseif but==117
   line( v(:,2), v(:,1), 'Marker', '+','Color','green', 'LineWidth', 2 );
   v=v(1:end-1,:);
   line( v(:,2), v(:,1), 'Marker', '+','Color','red', 'LineWidth', 2 );
elseif but==101
   stopinp=1;
end
end
