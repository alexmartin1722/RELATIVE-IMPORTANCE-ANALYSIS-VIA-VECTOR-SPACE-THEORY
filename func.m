% 1. Input matrix X
% 2. Size(X)= (rows, col)=(a,b)
% 3. N=0
% 4. N=N+1
% 5. Generate random number between 0 and a, k_N (I think just use a random number generator)
% 6. Delete rows N and k_N from X,
% 7. Size (X(N))=(a-2,b)
% 8. Here we must normalize X transform by X i =(x ji -E(X i ))/||X i || (Sorry, I forgot this)
% 9. A=first 3*b/4 columns of X(N) (here we want about ¾ of the columns)
% 10. B=all columns of X(N) except the last column
% 11. C=last column of X(N)
% 12. (QA,RA)=householderQR(A)
% 13. (QB,RB)=householderQR(B}
% 14. D(N)=QA*QA’*C
% 15. E(N)=QB*QB’*C
% 16. d(N)=||D(N)|| (norm) (|| K||=sqrt(sum(k_i)^2 ))
% 17. e(N)=||E(N)||
% 18. store (e(N)),d(N)) a 2 column matrix J size(J)=(a-2,2)
% 19. if N<a-2 go to 4, otherwise stop and write J

function func(X)
    [a,b] = size(X);
    N = 0;
    while N < a-2
        [a,b] = size(X);
        N = N+1;
        kn = randi([1,a]);
        minimum = min(kn, N);
        maximum = max(kn, N);
    
        tempMin = X(1:minimum-1, :);
        tempBetween = X(minimum+1:maximum-1,:);
        tempMax = X(maximum+1:end, :);
    
        X = cat(1, tempMin, tempBetween);
        X = cat(1, X, tempMax);
    
        [x,y] = size(X);
%         if x == a-2 && b==y
%              fprintf('properly removed rows\n')
%         end
    
    %     normalize X transform by Xsubi = x_{ij} -E(X_{i}/ ||X_{i}||
        T = X;
        for i = 1:x
            for j = 1:y
                T(i,j) = (X(i,j) - mean(X(:, j)))/norm(X(:,j));
            end
        end
        
        A = T(:, 1:(3*b)/4);
        B = T(:, 1:end-1);
        C = T(:, end);
        
        [QA1, RA1] = householderQR(A);
        [QB1, RB1] = householderQR(B);
        QA1*QA1' == QB1*QB1';
        QA*QA' == QB*QB';
        
        D = QA1*QA1' * C;
        E = QB1*QB1' *C;
        
        d = norm(D);
        e = norm(E);
    
        J = [e,d];
    end
    J
end

    

    


