# RELATIVE-IMPORTANCE-ANALYSIS-VIA-VECTOR-SPACE-THEORY
Code for the paper [A New Interpretation of Relative Importance on an Analysis of Per and Polyfluorinated Alkyl Substances (PFAS) Exposures on Bone Mineral Density](https://www.mdpi.com/1660-4601/20/5/4539)

#Running .m files 
For all matlab files running on 'women.xslx':
1) Read the data: women = readmatrix('women.xslx')
2) Replace all NaN values with 0: women(isnan(women))=0;
3) [Q,R] = function_name(matrix_name);
