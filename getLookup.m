h1 = load("lo_corrected.mat");
h2 = load("loo_corrected.mat");

s1 = orderfields(h1.files);
s2 = orderfields(h2.files);

for i = 1:length(s2)
    s2(i).const = s2(i).const+2;
end

for i = 1:(length(s1)+length(s2))
    if i <= length(s1)
        S(i) = s1(i);
    else
        S(i) = s2(i-length(s1));
    end
end


S_trimmed = rmfield(S, 'stim');
save("lookup_all.mat", 'S');
