function A = Bipart_rewire(A,ITER)

[Arows,Acols]=size(A);

K=nnz(A);
ITER=K*ITER;
maxAttempts = 10;
for i = 1:ITER
    cons_met = 0;
    attempt = 1;
    while ~cons_met
    while attempt <= maxAttempts
    r1 = randi(Arows,1);
    c1 = randi(Acols,1);
    v1 = A(r1,c1);
    r2 = randi(Arows,1);
    c2 = randi(Acols,1); 
    v2 = A(r2,c2);
    if r1 ~= r2 && v1 ~= 0 && v2 ~= 0 && c1 ~= c2
        A(r1,c1) = v2;
        A(r1,c1) = v2;
        cons_met = 1;
    end
    attempt = attempt + 1;
    end
    end
    
end