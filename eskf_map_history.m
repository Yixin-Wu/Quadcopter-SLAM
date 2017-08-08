function output = eskf_map_history(data)
A = zeros(10, 864);
for i = 1 : 864
    A(:, i) = eskf_map(data(i));
end
output = A;
end