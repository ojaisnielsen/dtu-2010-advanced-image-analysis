function result = init(cost)

result = zeros(size(cost, 1), 1);
mins = min(cost, [], 2);
for k=1:size(cost, 2)
    result = result + k * (result == 0) .* (cost(:,k) <= mins);
end

end
