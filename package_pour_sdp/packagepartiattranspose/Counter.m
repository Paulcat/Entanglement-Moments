function res=Counter(newval)
% un simple compteur
    persistent counter;
    if nargin>=1
        counter = newval;
        res = counter;
    else
        res = counter;
        counter = counter+1;
    end
end

