function v = allVL1eq(n, L1)

global MaxCounter;

n = feval(class(L1),n);
s = n+L1;
sd = double(n)+double(L1);
notoverflowed = double(s)==sd;
if isinf(MaxCounter) && notoverflowed
    v = allVL1nonrecurs(n, L1);
else
    v = allVL1recurs(n, L1);
end

end % allVL1eq