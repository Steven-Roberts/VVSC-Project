meshSizes = [64, 96, 128];
numMeshes = length(meshSizes);

sols = cell(numMeshes, 1);

for i = 1:numMeshes
    meshSize = meshSizes(i);
    sols{i} = solve_allen_cahn(meshSize, meshSize, 1, 0.25);
end

disErr = zeros(numMeshes - 1);
cgi = zeros(numMeshes - 1);
Fs = 3;
for i = 1:numMeshes-1
    for j = i+1:numMeshes
        disErr(i, j - 1) = (sols{j}.srq - sols{i}.srq) / ((meshSizes(j) / meshSizes(i))^2 - 1);
        cgi(i, j - 1) = Fs * abs(sols{i}.srq - sols{j}.srq) / ((meshSizes(j) / meshSizes(i))^2 - 1);
    end
end


% pHat = log((sols{1}.srq - sols{3}.srq) / (sols{3}.srq - sols{5}.srq)) / log(2);