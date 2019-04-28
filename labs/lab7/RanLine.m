%Generates a point set wher I points lay on a noisy line, and Out are random points
function x=RanLine(In,Out)

a=(rand(1,In)-0.5)*10;
b=a*0.5+randn(1,In)*0.25;
c=randn(2,Out)*2;
x=[[a;b],c];
perm=randperm(In+Out);
x=x(:,perm);

