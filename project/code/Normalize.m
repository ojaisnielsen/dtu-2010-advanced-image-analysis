function [X,T] = normalize(X)

X = X ./ (ones(3,1) * X(3,:));
Mean=mean(X, 2);
X(1,:)=X(1,:)-Mean(1);
X(2,:)=X(2,:)-Mean(2);
S=mean(sqrt(X(1,:).^2 + X(2,:).^2))/sqrt(2);
X(1:2,:)=X(1:2,:)/S;
T=[eye(2)/S,-Mean(1:2)/S;0 0 1];

end
