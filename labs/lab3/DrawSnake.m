function DrawSnake(v,color);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DrawSnake(v,colorcode);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Draw snake
%
% arguments:
% v = snake list
% color = (0,1,2,3) -> (black,red,green,blue)
if nargin == 0, error('Not enough input arguments.'); end
if nargin>2, error('Too many input arguments.'); end

if nargin ==1, color=1; end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% line color setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch color
   
case 0, p_color = 'black';
case 1, p_color = 'red';
case 2, p_color = 'green';
case 3, p_color = 'blue';
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% draw initial and optimal snake
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v = [v; v(1,:)];  % close curve by duplicating the first point
line( v(:,1), v(:,2), 'Marker', '+','Color',p_color, 'LineWidth', 2 );
