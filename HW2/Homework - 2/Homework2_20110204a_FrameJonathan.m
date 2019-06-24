%
% ------------------- HOMEWORK 2: DIFFERENTIAL EVOLUTION  -----------------
%
%                      HERE IS A SET OF TEST FUNCTIONS
% 
% THE TEST SET IS DIMENSION FREE. YOU CAN RUN THESE PROBLEMS IN ANY DIMENSION!
%
% THIS IS A TEST SET OF 25 PROBLEMS. 
% 
% HOMEWORK: PLEASE RUN PROBLEMS 1, 2, 5, 8, 10 and 12
% HOMEWORK: PLEASE RUN EACH OF THESE PROBLEMS IN DIMENSION 2, 5 and 10

% Information:
%
% minn -- minimum value of each parameter (first column is problem no. ; second column is minimum value)
% maxn -- maximum parameter value (first column is problem no.; second column is minimum value)
% accuracy -- if your objective function value is smaller than this you found the solution!
% minn / maxn: if you run in 2, 5 or 10 dimensions then all parameters have
% the same minimum and maximum values 

% Minimum values of each parameter (first column is function number)
minn = [1 -100; 2 -100; 3 -100; 4 -100; 5 -100; 6 -100; 7 0;   8 -32; 9 -5;
    10 -5; 11 -0.5; 12 -pi; 13 -3; 14 -100; 15 -5; 16 -5; 17 -5; 18 -5;
    19 -5; 20 -5; 21 -5; 22 -5; 23 -5; 24 -5; 25 -2];

% Maximun values of each parameter (first column is function number)
maxn = [1  100; 2  100; 3  100; 4  100; 5  100; 6  100; 7 600; 8  32; 9  5;
    10  5; 11  0.5; 12  pi; 13  1; 14  100; 15  5; 16  5; 17  5; 18  5;
    19  5; 20  5; 21  5; 22  5; 23  5; 24  5; 25 5];

% Which dimension do you like to run? (2, 5 or 10) - example is 5
dim = 10;

% Which function (problem) do you want to run? (minimize the function)
optfunction = 1;

       
 % Define initial_flag (dont change this!) -- needed to run program
global initial_flag; initial_flag = 0;

% Lets start with a random value to illustrate how to evaluate the function
% (always make sure x is within bounds!)

x = rand(1,dim);

% You can evaluate the objective function (OF) using:
OF = benchmark_func(x,optfunction); % (x is a horizontal vector with "dim" elements (the parameter values) )

% Define the Parameter Ranges (for initial sample with diferential evolution)
DE.minn = minn(optfunction,2) * ones(1,dim); DE.maxn = maxn(optfunction,2) * ones(1,dim);

% -------------------- Write you DE code here ---------------------

%%%%%%%%%Jonathan Frame%%%%%%%%%%%%%%%%%%%
%%%%%%%%%CEE 298 Winter 2011%%%%%%%%%%%%%%%%
%%%%%%%%%Differential Evolution%%%%%%%%%%%%%%%%%
DEN = 100; Converges = 0; FSR = 0; minCost = 1000000; minSet = ones(1,dim);
Convergance_ = ones(DEN-1,1)*100;
count = 0; ConvergeThresh = 1; accuracy = 0.001; 

% POPULATION, MUTANT POPULATION, TRIAL
x = zeros(DEN,dim); xmut = zeros(DEN,dim); U = zeros(DEN,dim);
CR = 0; r1 = 0; r2 = 0; r3 = 0; F = 0.4; K = 0.6;

fprintf('----------------BEGINNING DIFFERENTIAL EVOLUTION--------------\n'), optfunction

%INITIAL POPULATION
for i = 1:DEN
    for j = 1:dim
        x(i,j) = 2*((DE.maxn(1,j)/2) - rand*DE.maxn(1,j));
    end
end, 
tic,

%MAIN LOOP
while OF > accuracy

    MI = transpose(randperm(DEN)); %Muntant index
    for i = 1:DEN
        for  r1i = 1:DEN; r1 = max(round(rand(1)*DEN),1);
            if r2 ~= i; break; end;
        end
        for  r2i = 1:DEN; r2 = max(round(rand(1)*DEN),1); 
            if r2 ~= r1 && r2 ~= i; break; end;
        end
        for  r3i = 1:DEN; r3 = max(round(rand(1)*DEN),1);
            if r3 ~= r2 && r3 ~= r1 && r3 ~= i; break; end;
        end
        xmut(i,:) = x(i,:) + K * (x(r1,:) - x(i,:)) + F * (x(r2,:) - x(r3,:));        
    end

    for i = 1:DEN   % STAY IN BOUNDS
        for j = 1:dim
            if xmut(i,j) <DE.minn(1,j) 
                xmut(i,j) = DE.minn(1,j)-xmut(i,j); 
            end; 
        end
    end
    for i = 1:DEN 
        for j = 1:dim
            if xmut(i,j) >DE.maxn(1,j)
                xmut(i,j) = DE.maxn(1,j)-xmut(i,j); 
            end 
        end
    end
    
    for i = 1:DEN         %CROSSOVER
        for j = 1:dim
            if round(rand) == CR
                U(i,j) = xmut(i,j);
                else U(i,j) = x(i,j); 
            end
        end
    end

    for i = 1:DEN  %TRIAL VECTOR
        Gcost = benchmark_func(x(i,:),optfunction);
        if Gcost < minCost; minCost = Gcost; minSet = x(i,:); end
        Ucost = benchmark_func(U(i,:),optfunction);
        if Ucost < minCost; minCost = Ucost; minSet = x(i,:); end
        if Ucost < Gcost
            x(i,:) = U(i,:);
        end
    end

    for i = 1:DEN-1
        Convergance_(i,1) = sum(x(i+1,:) - x(i,:))^2;
    end

    OF = min(benchmark_func(mean(x),optfunction),benchmark_func(minSet,optfunction));
    Converge = mean(Convergance_);
    if Converge < ConvergeThresh
        Converges = Converges + 1
        for i = 1:DEN
            for j = 1:dim
                x(i,j) = 2*((DE.maxn(1,j)/2) - rand*DE.maxn(1,j));
            end
        end
    end
    %if Converges > 15;fprintf('\nCONVERGED, BUT DID NOT MEET ACCURACY\n'),break;end
    %accuracy = accuracy * 1.0005; 
    count = count + 1;
fprintf('--------------ITERATIONS---------'), count
fprintf('-----OBJECTIVE FUNCTION-------'), minCost
fprintf('------Optimum Parameters-------'), minSet
end
fprintf('--------------ITERATIONS---------'), count
fprintf('-----OBJECTIVE FUNCTION-------'), minCost
fprintf('------Optimum Parameters-------'), minSet
fprintf('----------------END DIFFERENTIAL EVOLUTION--------------------------\n')
toc,    

% HELP: You can use DE as structure: 
% DEN : number of individuals (population size), 20, 50 or 100 or more?
% DE.minn : minimum values of each parameter (dimension [1,dim])
% DE.maxn : maximum values of each parameter (dimension [1,dim])\