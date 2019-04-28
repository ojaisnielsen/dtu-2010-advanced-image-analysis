function [X,T] = normalize(X)

Mean=mean(X, 2);
X(1,:)=X(1,:)-Mean(1);
X(2,:)=X(2,:)-Mean(2);
S=mean(sqrt(diag(X'*X)))/sqrt(2);
X(1:2,:)=X(1:2,:)/S;
T=[eye(2)/S,-Mean(1:2)/S;0 0 1];

end
