function B = binarize(sFC,threshold)

B = sFC;
B(B>=threshold) = 1;
B(B<threshold) = 0;

end