function [x,y,P,Lines]=Film(nFrames,Draw)

if(Draw)
   close all
end


%[P,Lines]=Hus(Draw);
[P,Lines]=Hus;
nPoints=size(P,2);

x=zeros(nFrames,nPoints);
y=zeros(nFrames,nPoints);

Tfrom=[-10;-30;2];
Tto=[-30;-10;2];
T=zeros(3,nFrames);


for cFrame=1:nFrames
   alpha=(cFrame-1)/(nFrames-1);
   T(:,cFrame)=Tfrom+(Tto-Tfrom)*alpha+[-20;-20;-20]*(alpha*alpha-alpha);
   
   V=[1.5 1.2 2]+[0 0.1 0]*alpha+[1 0 0]*(alpha^2-alpha);  
	R=Rot(V*pi);
   [x(cFrame,:),y(cFrame,:)]=Project(P,T(:,cFrame),R);
end

for cFrame=1:nFrames
  	alpha=(cFrame-1)/(nFrames-1);
     
   V=[1.5 1.2 2]+[0 0.1 0]*alpha+[1 0 0]*(alpha^2-alpha);  
	R=Rot(V*pi);
	CamZ=-R(3,:);
end
   
if(Draw)
	for cFrame=1:nFrames
		figure
		DrawProjFigure(x(cFrame,:),y(cFrame,:),Lines)
   end
end




function [x,y]=Project(P,t,R)
temp=[R,-R*t]*[P;ones(1,size(P,2))];
x=temp(1,:)./temp(3,:);
y=temp(2,:)./temp(3,:);


function [P,Lines]=Hus()

P=[
   0 0 0
   0 1 0
   0 1 1
   0 0 1			%4
   1 0 0
   1 0 1			%6
   0 0.5 1.5
   1 0.5 1.5	%8
   0 1/3 0
   0 1/3 2/3
   0 2/3 0
   0 2/3 2/3	%12
   1/5 0 2/3
   2/5 0 2/3
   1/5 0 1/3
   2/5 0 1/3	%16
   3/5 0 2/3
   4/5 0 2/3
   3/5 0 1/3
   4/5 0 1/3	%20
];

%scale building
P=[1 0 0
   0 1 0
   0 0 1
   ]*P';

Lines=[
   1 2
   2 3 
   3 4
   4 1
   1 5
   5 6
   6 4
   3 7
   4 7
   7 8
   6 8
   9 10
   10 12
   12 11
   13 14
   14 16
   16 15
   15 13
   17 18
   18 20
   20 19
   19 17
];

Lines=Lines';


function R=Rot(V);

c=cos(V);
s=sin(V);

R=zeros(3,3);

R(1,1)=c(2)*c(3);
R(2,1)=-c(2)*s(3);
R(3,1)=s(2);
R(1,2)=c(1)*s(3)+s(1)*s(2)*c(3);
R(2,2)=c(1)*c(3)-s(1)*s(2)*s(3);
R(3,2)=-s(1)*c(2);
R(1,3)=s(1)*s(3)-c(1)*s(2)*c(3);
R(2,3)=s(1)*c(3)+c(1)*s(2)*s(3);
R(3,3)=c(1)*c(2);

function DrawProjFigure(x,y,Lines)

plot(x,y,'*');

hold on
for i=1:size(Lines,2),
	plot([x(Lines(1,i)) x(Lines(2,i))],[y(Lines(1,i)) y(Lines(2,i))],'-');
end

hold off

