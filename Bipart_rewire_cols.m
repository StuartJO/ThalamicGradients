function A = Bipart_rewire_cols(A,ITER)

[~,Acols]=size(A);

for col = 1:Acols
    rows2rewire = find(A(:,col)>0);
    for iter = 1:ITER
    r = randperm(length(rows2rewire));
    A(rows2rewire,col) = A(rows2rewire(r),col);
    end
end