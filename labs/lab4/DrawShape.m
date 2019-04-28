function DrawShape(v,color);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DrawShape(shape,colorcode);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% arguments:
% shape = shape list, size=(np,2), where np is the number of 
%         points in the shape
% color = (0,1,2,3) -> plotting using (black,red,green,blue)
%
% Note: the shape will automatically be closed by connecting
% the first and the last point marked.
%
% See: MarkShape

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

v=[v; v(1,:)];
line( v(:,2), v(:,1), 'Marker', '+','Color',p_color, 'LineWidth', 2 );
