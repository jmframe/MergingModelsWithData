%%%%%%%%%Jonathan Frame%%%%%%%%%%%%%%%%%%%
%%%%%%%%%CEE 298 Winter 2011%%%%%%%%%%%%%%%%
%%%%%%%%%SCE-UA%%%%%%%%%%%%%%%%%

m = 10; n = 3;
t = [0.6, 2.3, 3.1, 4.4, 5.1, 6, 7, 7.8, 8.4, 9.9];
X = zeros(m,n);  SimplexError = 1001;
PlusOrminus = [0;0;0]; count = 0;
P1bound = 50; P2bound = 150; P3bound = 50;
alpha = 1; gamma = alpha * 3; 
betaplus = alpha * (3/4); betaminus = alpha*(1/4);
RECC = [alpha, gamma, betaplus, betaminus];
dm = [105; 195; 270; 335; 395; 435; 460; 500; 515; 505];

for u=1:m; for v=1:n
        if v == 1; X(u,v)=1; end
        if v == 2; X(u,v)=u; end
        if v == 3; X(u,v)=-(1/2)*u^2; end
end; end

shc = zeros(n+1,n, n+1);
for si = 1:n+1; for sj = 1:n; for sk = 1:n+1
            if sj == 1; shc(si, sj, sk) = rand*P1bound; end
            if sj == 2; shc(si, sj, sk) = rand*P2bound; end
            if sj == 3; shc(si, sj, sk) = rand*P3bound; end
end; end; end

for go = 1:5
    
    for i = 1:n+1
        s1 = shc(i,:,1)'; s2 = shc(i,:,2)'; s3 = shc(i,:,3)'; s4 = shc(i,:,4)';
        badreflection = 0; 
    
        while badreflection < n + 2 
            S = [s1';s2';s3';s4'];
            g = [mean(S(:,1));mean(S(:,2));mean(S(:,3))];
    
            SSEs1 = transpose(dm-X*s1)*(dm-X*s1);
            SSEs2 = transpose(dm-X*s2)*(dm-X*s2);
            SSEs3 = transpose(dm-X*s3)*(dm-X*s3);
            SSEs4 = transpose(dm-X*s4)*(dm-X*s4);
            SSEs = [SSEs1, SSEs2, SSEs3, SSEs4];
            W = max(SSEs);
    
            if W == SSEs1; W = s1; end
            if W == SSEs2; W = s2; end
            if W == SSEs3; W = s3; end
            if W == SSEs4; W = s4; end
    
            if badreflection >0
                if W == s1; W = s2; end
                if W == s2; W = s3; end
                if W == s3; W = s4; end
                if W == s4; W = s1; end
            end
            SSEW = transpose(dm-X*W)*(dm-X*W);

            for j = 1:4
                reflection = g + RECC(1,j)*(g-W);
        
                for k = 1:n;
                    if reflection(k,1) < 0; reflection(k,1) = -reflection(k,1);end;
                end
        
                SSEreflection = transpose(dm-X*reflection)*(dm-X*reflection);
        
                if SSEreflection < SSEW
                    badreflection = 0;
                    if W == s1; s1 = reflection; end
                    if W == s2; s2 = reflection; end
                    if W == s3; s3 = reflection; end
                    if W == s4; s4 = reflection; end
                    break
                end
                badreflection =  badreflection + 1;
            end
            SimplexError = transpose(dm-X*g)*(dm-X*g);
            count = count + 1;
        end
        shc(i,:,1) = s1'; shc(i,:,2) = s2'; shc(i,:,3) = s3'; shc(i,:,4) = s4';
    end

    for si = 1:n+1; for sj = 1:n; for sk = 1:n+1
                shc(:,sj,sk) = shc(randperm(n+1),sj,sk);
    end;end;end
    
end

minMdir = inv(transpose(X)*X)*transpose(X)*dm; 
var = SimplexError/(m - n);
Covariance = var*inv(transpose(X)*X); tscore = 2.8;
for i = 1:n; PlusOrMinus(i,1) = sqrt(Covariance(i,i))*tscore; end
g, PlusOrMinus, SimplexError, count
scatter(t,dm,'r'), line(1:m,X*g),
title('SCE-UA'), legend('data', 'model', 'Location', 'SouthEast')