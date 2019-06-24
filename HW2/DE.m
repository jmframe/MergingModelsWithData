


m = 10; %DATA POINTS
n = 3; %PARAMETERS
P = 100; %POPULATION SIZE
count = 0; DEerror = 100000; ErrorThresh = 1000; 
PlusOrminus = [0;0;0]; %FOR ERROR ESTIMATE
G = zeros(P,n); Gmut = zeros(P,n);
r1 = 0; r2 = 0; r3 = 0;

%DATA
t = [0.6, 2.3, 3.1, 4.4, 5.1, 6, 7, 7.8, 8.4, 9.9];
dm = [105; 195; 270; 335; 395; 435; 460; 500; 515; 505];

% FOR OBJECTIVE FUNCTION
X = zeros(m,n); 
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

%INITIAL POPULATION
for i = 1:P
    for j = 1:n
        G(i,j) = rand*100;
    end
end

%MAIN LOOP
%while DEerror > ErrorThresh

    %Set Mutant Matrix
    MI = transpose(randperm(P)); %Mutant Index
    for i = 1:P
        for  r1i = 1:P; r1 = round(rand(1)*100) ; if r2 ~= i; break; end; end
        for  r2i = 1:P; r2 = round(rand(1)*100) ; if r2 ~= r1; break; end; end
        for  r3i = 1:P; r3 = round(rand(1)*100) ; if r3 ~= r2; break; end; end        
    end
    












%
