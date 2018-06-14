function returnValue = gaussianQuadrature(f, a, b)
    w = [1, 1];
    ksi = [-1/sqrt(3), 1/sqrt(3)];
    returnValue = (b - a)/2 .* sum(w .* f((b - a)/2 .* ksi + (b + a)/2 ));
end