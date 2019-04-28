function vcircle=CircleSnake(c,nbPoints,radius);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% vcircle=CircleSnake(c,nbPoints,radius);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate circle snake (center in the point c)
%
% arguments:
% c = Center-vector for the snake.
% nbPoints = number of points on the snake
% radius = radius of the circle (in pixels)

if nargin == 0, error('Not enough input arguments.'); end
if nargin>3, error('Too many input arguments.'); end
if length(c)~=2, error('Length of center vector != 2.'); end

if nargin == 1,
	nbPoints = 16;
   radius   = 100; 
elseif nargin == 2,
   radius=100;
end

   
	offSet   = (2*pi)/nbPoints;
	
	v   = [ c(1)+radius c(2) ];
	rad = offSet;
	
	for i=1:nbPoints-1
	
	    v = [ v; c(1)+radius*cos(rad) c(2)+radius*sin(rad) ];
	
	    rad = rad + offSet;
	end
	
	vcircle = round( v );

