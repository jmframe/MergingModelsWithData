sshc = shc;

r = randperm(n+1);
for si = 1:n+1
    for sj = 1:n
        for sk = 1:n+1
            sshc(:,sj,sk) = sshc(randperm(n+1),sj,sk)
        end
    end
end
sshc