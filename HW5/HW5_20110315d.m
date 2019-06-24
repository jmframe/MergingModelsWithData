tic,
[A,B,C,D,E,F,G,H] = textread('data.txt',...
    '%f %f %f %f %f %f %f %f');
[Y] = textread('Y.txt', '%f');
N = 8
Len = length(Y);
n = 365*8; % Calibration Time
Np = [3,4,5,8,8,9,9,10];
X = [A,B,C,D,E,F,G,H];
XC = X(1:n,:); YC = Y(1:n,:);
XF = X(n+1:Len,:); YF = Y(n+1:Len,:);

BGA_F = zeros(Len - n,1);
ICA_F = zeros(Len - n,1);
GRA_F =  zeros(Len - n,1);
EW_F =  zeros(Len - n,1);
betaBGA = zeros(N,1);
betaICA = zeros(N,1);
ForVar = zeros(1,N);
InvForVar = zeros(1,N);
LL = zeros(1,N);
I = zeros(1,N);
eToTheIoverTWO = zeros(1,N);

%%%%% SSE %%%%%%%%%%%%%%%%%%
for i = 1:N
    SSE(i) = (XF(:,i)'*XF(:,i))\XF(:,i)'*YF;
end
SSE

%% EQUAL WEIGHTS%%%%%%%%%%%%%%%

betaEW = ones(N,1)/N;

for i = 1:n
    EW_F(i,1) = XF(i,:)*betaEW;
end

SSE_EW = (EW_F'*EW_F)\EW_F'*YF

%% BATES-GRANGER%^%%%%%%%^^%&

for i = 1:N
    ForVar(i) = var(XC(:,i)-YC);
    InvForVar(i) = 1/ForVar(i);
end

for i = 1:N
    betaBGA(i,1) = (1/ForVar(i))/sum(InvForVar);
end

for i = 1:n
    BGA_F(i,1) = XF(i,:)*betaBGA;
end

SSE_BGA = (transpose(BGA_F)*BGA_F)\transpose(BGA_F)*YF

%%%%Information criterion averaging%%%%%%

for i = 1:N
    LL(i) = n * log(ForVar(i)) + n;
    I(i) = LL(i) + 2 * Np(i);
    eToTheIoverTWO(i) = -I(i)/2; %exp(-I(i)/2);
end
for i = 1:N
    betaICA(i,1) = eToTheIoverTWO(i) / sum(eToTheIoverTWO);
end

for i = 1:n
    ICA_F(i,1) = XF(i,:)*betaICA;
end

SSE_BGA = (BGA_F'*ICA_F)\ICA_F'*YF

%%%%%Granger-Ramanathan averaging%%%%%%

betaGRA = (transpose(XC)*XC)\transpose(XC)*YC
for i = 1:n
    GRA_F(i,1) = XF(i,:)*betaGRA;
end

toc,