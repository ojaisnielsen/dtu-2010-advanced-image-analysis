function drawCarShape(v, color)
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

v = [imag(v), real(v)];
drawLine(v(1,:), v(2,:), p_color);
drawLine(v(3,:), v(4,:), p_color);
drawLine(v(4,:), v(1,:), p_color);
drawLine(v(4,:), v(5,:), p_color);
drawLine(v(5,:), v(6,:), p_color);
drawLine(v(6,:), v(7,:), p_color);
drawLine(v(7,:), v(2,:), p_color);
drawLine(v(3,:), v(8,:), p_color);
drawLine(v(8,:), v(9,:), p_color);
drawLine(v(9,:), v(14,:), p_color);
drawLine(v(14,:), v(10,:), p_color);
drawLine(v(10,:), v(11,:), p_color);
drawLine(v(11,:), v(15,:), p_color);
drawLine(v(15,:), v(12,:), p_color);
drawLine(v(12,:), v(13,:), p_color);
drawLine(v(13,:), v(7,:), p_color);
end
