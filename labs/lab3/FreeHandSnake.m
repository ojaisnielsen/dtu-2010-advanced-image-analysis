function v=FreeHandSnake(go);
go =1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% vFH=FreeHandSnake
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Usage:
%
% Make snake points with the left mouse button.
% Stop by pressing the right mouse button.
%
% Note: the snake will automatically be closed by connecting
% the first and the last point marked.

stopinp=0;
v=[];
disp('left-click : Mark a snake point');
disp('   ''u''     : Undo last point')
disp('   ''e''     : End');
while stopinp==0,
   
[vx,vy,but]=ginput(1);
if but==1
   v=[v;ceil(vx) ceil(vy)];
	line( v(:,1), v(:,2), 'Marker', '+','Color','red', 'LineWidth', 2 );
elseif but==117
   v=v(1:end-1,:);
elseif but==101
   stopinp=1;
end
end
