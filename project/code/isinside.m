function result = isinside(S, P)

sign = 0;
for i=1:size(S,2),
    result = sign * det([P - S(:,i), S(:,mod(i, size(S,2))+1) - S(:,i)]) >= 0;
    if not(result),
        break
    end
    sign = det([P - S(:,i), S(:,mod(i, size(S,2))+1) - S(:,i)]);
end

end