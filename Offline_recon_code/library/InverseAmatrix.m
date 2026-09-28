function A_tik = InverseAmatrix(A,lambdaTikPercent)
[U, S, V] = svd(A, 'econ');
s = diag(S);
lambdaTik = lambdaTikPercent*max(abs(s));
S_tik = diag(s ./ (s.^2 + lambdaTik));
A_tik = V * S_tik * U';
end

