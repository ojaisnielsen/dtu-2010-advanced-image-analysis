function [best, F, H] = robustmatches(ransacprob, thres, homogconstraint, points, matches)

for i=1:size(points,2)-1,

    maxScore=0;
    n = 0;
    N = inf; 
    while N > n,
        perm = randperm(size(matches{i}, 2));
        
        if homogconstraint,
            currentMatches = matches{i}(:,perm(1:4));
            currentH = homography(points{i}(:, currentMatches(1,:)), points{i+1}(:, currentMatches(2,:)), true);
            H1 = currentH * points{i}(:, matches{i}(1,:));
            H1 = H1 ./ (ones(3,1)*H1(3,:));
            invH2 = inv(currentH) * points{i+1}(:, matches{i}(2,:));
            invH2 = invH2 ./ (ones(3,1)*invH2(3,:)); 
            
            d = sum((points{i}(:, matches{i}(1,:)) - invH2).^2 +(points{i+1}(:, matches{i}(2,:)) - H1).^2);
        else
            currentMatches = matches{i}(:,perm(1:8));
            currentF = fundamentalmat(points{i}(:, currentMatches(1,:)), points{i+1}(:, currentMatches(2,:)), false);
            points2tFpoints1 = zeros(1,size(matches{i}, 2));
            for n = 1:size(matches{i}, 2),
                points2tFpoints1(n) = points{i+1}(:,matches{i}(2,n))'*currentF*points{i}(:,matches{i}(1,n));
            end

            Fpoints1 = currentF*points{i}(:,matches{i}(1,n));
            Ftpoints2 = currentF'*points{i+1}(:,matches{i}(2,n));
            d =  points2tFpoints1.^2 ./ (Fpoints1(1,:).^2 + Fpoints1(2,:).^2 + Ftpoints2(1,:).^2 + Ftpoints2(2,:).^2);            
        end



        inliers = (abs(d) < thres);

        if (sum(inliers) > maxScore),
            maxScore = sum(inliers);
            best{i} = matches{i}(:, inliers);
            if homogconstraint,
                bestH = currentH;
            else
                bestF = currentF;
            end
        end

        eps = 1-(sum(inliers)/size(matches{i}, 2));
        if homogconstraint,        
            N = log(1-ransacprob)/log(1 - (1 - eps)^4);    
        else
            N = log(1-ransacprob)/log(1 - (1 - eps)^8);              
        end
        n = n+1;
    end
    n
    
    %H(:,:,i) = bestH;
    H(:,:,i) = homography(points{i}(:, best{i}(1,:)), points{i+1}(:, best{i}(2,:)), true);   
    %F(:,:,i) = bestF;
    F(:,:,i) = fundamentalmat(points{i}(:, best{i}(1,:)), points{i+1}(:, best{i}(2,:)), false); 
    
end

% Find tuples
for i=1:size(points,2)-2,
    [common, i1, i2] = intersect(best{1}(end,:), best{i+1}(1,:));
    best{1} = [best{1}(:,i1); best{i+1}(2,i2)];
end
best = best{1};

end
