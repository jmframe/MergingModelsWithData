%PROBLEM 1%%%%%%%%%%%%%%%%%%%


minSSR = 100000000; minM = 1:2;
m = 10; n = 2; 
t = 1:m; T = zeros(m,n); Y = zeros(1,t);
intercept = 21; slope = 45; PlusOrMinus = [0;0];
SSR_ = zeros(1,1000);
for i = t
    Error1 = (0.5-rand) *20; Error2 = (0.5-rand) * 10;  
    Y(i) = (slope + Error1 ) * i + intercept + Error2;
    Y = transpose(Y);
end
for u=1:m
    for v=1:n
        if v == 1
            T(u,v)=1;
        end
        if v == 2;
            T(u,v)=u;
        end
    end
end

M = inv(transpose(T)*T)*transpose(T)*Y;
Yhat = T*M;
E = Y - Yhat;
SSR = transpose(E) * E;
subplot(3,1,1), scatter(t,Y,'r'), line(t,Yhat)
title('linear'); legend('data', 'Model',  'Location', 'SouthEast')
var = SSR/(m - n);
C = var*inv(transpose(T)*T);
tscore = 2.1;
PlusOrMinus(1,1) = sqrt(C(1,1))*tscore;
PlusOrMinus(2,1) = sqrt(C(2,2))*tscore;
fprintf('------------------------problem 1-----------------------------------')
M, PlusOrMinus, SSR


%%%%%%%%%%%%%%%%%%%%%%%%%
%PROBLEM 2%%%%%%%%%%%PROBLEM 2

m = 10; n = 3; count = 0; countmin = 0; ConvergeThresh = 40;
t = [0.6, 2.3, 3.1, 4.4, 5.1, 6, 7, 7.8, 8.4, 9.9];%t = 1:m; 
X = zeros(m,n); SSError = 1000;
delta_m = zeros(1,n);
dm = [105; 195; 270; 335; 395; 435; 460; 500; 515; 505];

for u=1:m
    for v=1:n
        if v == 1
            X(u,v)=1;
        end
        if v == 2;
            X(u,v)=u;
        end
        if v == 3;
            X(u,v)=-(1/2)*u^2;
        end
    end
end

M = [-20; 500; 50]; dp = X * M;
Om = transpose(dm)*dm - 2*transpose(M)*transpose(X)*dm + transpose(X * M)*(X*M);
gm = -2 * transpose(X) * dm + 2 * transpose(X)*X*M;
J = gm*transpose(dm);
Hm = 2*transpose(X)*X; alpha = inv(Hm);

 for i = 0:10
    delta_m = -alpha*gm;
    M = M + delta_m;    
    dp = X * M; 
    count = count +1;
    gm = -2 * transpose(X) * dp + 2 * transpose(X)*X*M; 
    J = gm*transpose(dp);
    SSError = transpose(dm - dp)*(dm - dp);
 end
 

subplot(3,1,2), scatter(t,dm,'r'), line(1:m,dp), 
title('Newtons Method'), legend('data', 'model', 'Location', 'SouthEast')
var = SSError/(m - n);
Covariance = var*inv(transpose(X)*X);
tscore = 2.1;
PlusOrMinus(1,1) = sqrt(Covariance(1,1))*tscore;
PlusOrMinus(2,1) = sqrt(Covariance(2,2))*tscore;
PlusOrMinus(3,1)= sqrt(Covariance(3,3))*tscore;
fprintf('-------------------------problem 2-------------------------')
M, PlusOrMinus, SSError


%%%%%%%%%%%%%%%%%%%%%%%%%
%PROBLEM 3%%%%%%%%%%%%%%%%%%%


m = 10; n = 3; count = 0;
t = [0.6, 2.3, 3.1, 4.4, 5.1, 6, 7, 7.8, 8.4, 9.9];
X = zeros(m,n);  SimplexError = 1001;
PlusOrminus = [0;0;0];
badreflection = 0;
alpha = 1; gamma = alpha * 2; 
betaplus = alpha * (3/4); betaminus = alpha*(1/2);
RECC = [alpha, gamma, betaplus, betaminus];
dm = [105; 195; 270; 335; 395; 435; 460; 500; 515; 505];


for u=1:m
    for v=1:n
        if v == 1
            X(u,v)=1;
        end
        if v == 2;
            X(u,v)=u;
        end
        if v == 3;
            X(u,v)=-(1/2)*u^2;
        end
    end
end

s1 = [rand*30; rand*200; rand*20];
s2 = [rand*30; rand*200; rand*20];
s3 = [rand*30; rand*200; rand*20];
s4 = [rand*30; rand*200; rand*20];
InitialSimplex = [s1,s2,s3, s4];

while SimplexError >= 300 
    
    g = [...
        (s1(1,1)+s2(1,1)+s3(1,1)+s4(1,1))/4;
        (s1(2,1)+s2(2,1)+s3(2,1)+s4(2,1))/4;
        (s1(3,1)+s2(3,1)+s3(3,1)+s4(3,1))/4 
    ];
    
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

    if badreflection > 4
    s1 = [rand*30; rand*200; rand*20];
    s2 = [rand*30; rand*200; rand*20];
    s3 = [rand*30; rand*200; rand*20];
    s4 = [rand*30; rand*200; rand*20];
    InitialSimplex = [s1,s2,s3, s4];
    badreflection = 0;
    end
    
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

count; 
subplot(3,1,3), scatter(t,dm,'r'), line(1:m,X*g),
title('Simplex Method'), legend('data', 'model', 'Location', 'SouthEast')
minMdir = inv(transpose(X)*X)*transpose(X)*dm; 
var = SimplexError/(m - n);
Covariance = var*inv(transpose(X)*X); tscore = 2.1;
PlusOrMinus(1,1) = sqrt(Covariance(1,1))*tscore;
PlusOrMinus(2,1) = sqrt(Covariance(2,2))*tscore;
PlusOrMinus(3,1)= sqrt(Covariance(3,3))*tscore;
fprintf('-----------------problem 3------------------')
g, PlusOrMinus, SimplexError