n = 2;
numGs = 2;
for i = 2:n
    fprintf("i = %d\n", i);
    [A, B, C, S, Peqz] = BMIConstruction(i,numGs);
    save([num2str(numGs),'G ',num2str(i), 'x', num2str(i), 'BMI system test.mat'], 'A', 'B', 'C', 'S', 'Peqz');
end