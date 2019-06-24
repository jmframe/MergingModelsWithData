%%%%%%%%%Jonathan Frame%%%%%%%%%%%%%%%%%%%
%%%%%%%%%CEE 298 Winter 2011%%%%%%%%%%%%%%%%
%%%%Multi Objective Opt || Pareto Ranking ||Differential Evolution%%%

dim = 2; DEN = 100; numOF = 3; maxn = 5; minn = 0;
ObjVals = ones(2*DEN, numOF);
CR = 0.5; r1 = 0; r2 = 0; r3 = 0; F = 0.1; K = 0.9;
fprintf('----------------BEGINNING DIFFERENTIAL EVOLUTION--------------\n')

%INITIAL POPULATION
x = rand(2*DEN,dim) * maxn; 
xmut = zeros(DEN,dim); U = zeros(DEN,dim);

%MAIN LOOP
for ii = 1:50

    %Muntant Matrix
    for i = 1:DEN
        for  r1i = 1:DEN; r1 = max(round(rand*DEN),1);
            if r2 ~= i; break; end;
        end
        for  r2i = 1:DEN; r2 = max(round(rand*DEN),1); 
            if r2 ~= r1 && r2 ~= i; break; end;
        end
        for  r3i = 1:DEN; r3 = max(round(rand*DEN),1);
            if r3 ~= r2 && r3 ~= r1 && r3 ~= i; break; end;
        end
        xmut(i,:) = x(i,:) + K * (x(r1,:) - x(i,:)) + F * (x(r2,:) - x(r3,:));        
    end

    for i = 1:DEN   % STAY IN BOUNDS
        for j = 1:dim
            if xmut(i,j) <minn
                xmut(i,j) = minn - xmut(i,j);
            end; 
            if xmut(i,j) >maxn
                xmut(i,j) = maxn + (maxn-xmut(i,j)); 
            end 
            %CROSSOVER
            if rand <= CR
                U(i,j) = xmut(i,j);
                else U(i,j) = x(i,j); 
            end
        end
    end

    for i = 1:DEN  %Objective Values
        ObjVals(i,1) = cs41OF1(x(i,:));
        ObjVals(i,2) = cs41OF2(x(i,:));
        ObjVals(i,2) = cs41OF3(x(i,:));
        % Trial matrix Objective Values
        ObjVals(i+DEN,1) = cs41OF1(U(i,:));
        ObjVals(i+DEN,2) = cs41OF2(U(i,:));
        ObjVals(i+DEN,3) = cs41OF3(U(i,:));
        %Update population
        x(i+DEN,:) = U(i,:);
    end
    
    %Pareto Ranking
    xrank = ParetoRanking(ObjVals);
    
    for i = 1:DEN  %%%%%%%%%TRIAL INDIVIDUAL
        if xrank(i+DEN,:) < xrank(i,:);
            x(i,:) = x(i+DEN,:);
        end
    end

    scatter(x(1:DEN,1),x(1:DEN,2), 'b'); hold on...
        ;scatter(x(DEN:DEN*2,1)...
        ,x(DEN:DEN*2,2), 4 , 'r'),hold off
     legend('population', 'trial pop')
     pause(0.05),
     
     mean_rank = mean(xrank(1:DEN, :))
     if sum(xrank(1:DEN, :)) == DEN; break; end
end
