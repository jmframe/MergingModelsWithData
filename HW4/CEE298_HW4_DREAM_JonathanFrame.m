%%%%%%%%%%%%%%%%%%%%%%%%%
% Jonathan Frame 27288822
% UCI CEE 298 Winter 2011
% Code for plotting results was adapted from... 
%   http://psiexp.ss.uci.edu/research/teachingP205C/205C.pdf

for HWi = 1:3
    
fprintf('---------------------------------------------------')
HW = HWi, min = -10; max = 10;
N = 5; acc = 0; rej = 0; its = 50000;
if HW < 3; n = 1; elseif HW ==3; n = 10; end 
z = zeros(N, n); gamma = 2.4/sqrt(2*n); e = 1;
Chain = zeros(its, N, n); sigma = 1;
M = rand(N, n)*max + rand(N, n)*min;

for i = 1:its

    for ii = 1:N
        idx = randperm(N);
        idx(idx == ii) = [];
        r1 = idx(1); r2 = idx(2);
        e =  rand(1, n) - rand(1, n);
        z(ii, :) = M(ii, :) + gamma * (M(r1, :) - M(r2, :)) + e;
    end
    
    for ii = 1:N
        
        if HW == 1
            
            proposal = normpdf(z(ii), 0, sigma);
            current = normpdf(M(ii), 0, sigma);
        
        elseif HW == 2
        
            proposal = (1/3)*normpdf(z(ii), -5, sigma) +...
                (2/3)*normpdf(z(ii), 5, sigma);
            current = (1/3)*normpdf(M(ii), -5, sigma) +...
                (2/3)*normpdf(M(ii), 5, sigma);
        
        elseif HW == 3
        
            A = 0.5*eye(n) + 0.5*ones(n);
            for iii=1:n
               for jjj=1:n
                   C(iii,jjj) = A(iii,jjj) * sqrt(iii * jjj);
               end
            end
            proposal = mvnpdf(z(ii, :), 0, C); 
            current = mvnpdf(M(ii, :), 0, C);
        
        end
        
        alpha = proposal/current;

        if alpha >= rand; 
            M(ii, :) = z(ii, :); acc = acc+1; else rej = rej +1;
        end
        
        if HW < 3
            Chain(i,ii) = M(ii);
        end
        
        if HW ==3
            Chain(i,ii, :) = M(ii, :);
        end
        
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%_____Plotting Results____________________
% http://psiexp.ss.uci.edu/research/teachingP205C/205C.pdf
%%%%%%%%%%%%%%%%%%%%%%%%^%%

    figure(HWi); clf; 
    if HWi < 3; 
        L = length(Chain);
        ChainPlot = Chain(round(L-L*0.8):L,1);
        subplot( 3,1,1:2 );
        bins = linspace( -10 , 10 , 50 );
        counts = hist( ChainPlot , bins );
        bar( bins , counts/sum(counts) , 'r' );
        xlim( [ -10 10 ] );
        xlabel( 'x' ); ylabel( 'p(x)' );

        if HW ==1
            normal = normpdf(bins, 0, 1);
        elseif HW == 2
            normal = (1/3)*normpdf(bins, -5, sigma) +...
                (2/3)*normpdf(bins, 5, sigma);
        end
        
        hold on;
        plot( bins , normal/sum(normal) , 'b--' , 'LineWidth' , 3 );
        set( gca , 'YTick' , [] );

        subplot( 3,1,3 );
        stairs( Chain )
            
    elseif HWi ==3; 
        for ijk = 1:9
            L = length(Chain);
            ChainPlot = Chain(round(L-L*0.8):L,1, ijk);
            subplot( 3,3,ijk );
            bins = linspace( -5 , 5 , 50 );
            counts = hist( ChainPlot , bins );
            bar( bins , counts/sum(counts) , 'r' );
            xlim( [ -5 5 ] );
            xlabel( ijk ); ylabel( 'p(x)' );
            normal = normpdf(bins, 0, 1);
            hold on;
            plot( bins , normal/sum(normal) , 'b--' , 'LineWidth' , 3 );
            set( gca , 'YTick' , [] );
        end
        
    end
    
    hold off;

    
AcceptanceRatio = acc / (acc + rej)
fprintf('---------------------------------------------------')

end