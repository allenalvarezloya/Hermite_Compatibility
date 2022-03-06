load('L1errors.ext')
load('L2errors.ext')
load('Linferrors.ext')
m = 3;
Ncells = [5 10 20];

conv1 = log2(L1errors(1:end-1)./L1errors(2:end));
conv2 = log2(L2errors(1:end-1)./L2errors(2:end));
convinf = log2(Linferrors(1:end-1)./Linferrors(2:end));

%% Error Table
fprintf("\\begin{table}[ht]\n ")
fprintf(sprintf("\\\\caption{Example 1 mx = %i} \n",m))
fprintf("\\begin{center} \n ")
fprintf("\\begin{tabular}{|c|c|c|c|c|c|c|} \n")
fprintf("\\hline \n")
fprintf(" n & $L_1$ error & Convergence & $L_2$ error & Convergence & $L_{\\infty}$ error & Convergence \\\\ \n ")
fprintf("\\hline \n")
fprintf(sprintf(' %i & %3.2e & - & %3.2e & - & %3.2e & - \\\\\\\\ \n',Ncells(1),L1errors(1),L2errors(1),Linferrors(1)))
for z = 2:length(Ncells)
    fprintf(sprintf('%i & %3.2e & %3.2f & %3.2e & %3.2f & %3.2e & %3.2f \\\\\\\\ \n',Ncells(z),L1errors(z),conv1(z-1),L2errors(z),conv2(z-1),Linferrors(z),convinf(z-1)))
end
fprintf("\\hline \n")
fprintf("\\end{tabular} \n")
fprintf("\\end{center} \n")
fprintf("\\end{table} \n")