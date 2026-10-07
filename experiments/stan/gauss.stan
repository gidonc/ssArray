data { int D; matrix[D, D] L; }            // target: theta ~ multi_normal(0, L L')
parameters { vector[D] theta; }
model { target += -0.5 * dot_self(mdivide_left_tri_low(L, theta)); }
