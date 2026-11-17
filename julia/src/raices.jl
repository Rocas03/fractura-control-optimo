# Regula falsi de Illinois: superlineal y siempre dentro del intervalo.
function raiz(g, a, b; xtol = 4eps(), maxit = 200)
    fa, fb = g(a), g(b)
    fa == 0 && return a
    fb == 0 && return b
    sign(fa) != sign(fb) || throw(ArgumentError("raiz: sin cambio de signo en [$a, $b]"))
    lado = 0
    for _ in 1:maxit
        abs(b - a) ≤ xtol * max(1.0, abs(a), abs(b)) && break
        c = (a * fb - b * fa) / (fb - fa)
        (min(a, b) < c < max(a, b)) || (c = (a + b) / 2)
        fc = g(c)
        fc == 0 && return c
        if sign(fc) == sign(fb)
            b, fb = c, fc
            lado == -1 && (fa /= 2)
            lado = -1
        else
            a, fa = c, fc
            lado == 1 && (fb /= 2)
            lado = 1
        end
    end
    (a + b) / 2
end

# Newton 2×2 con Jacobiano por diferencias centradas y búsqueda lineal por mitades.
function jac_fd(F, z; δ = 1e-6)
    h1, h2 = δ * max(1.0, abs(z[1])), δ * max(1.0, abs(z[2]))
    e1, e2 = SVector(h1, 0.0), SVector(0.0, h2)
    hcat((F(z + e1) - F(z - e1)) / 2h1, (F(z + e2) - F(z - e2)) / 2h2)
end

function newton2(F, z0; tol = 1e-14, maxit = 50)
    z = SVector{2, Float64}(z0)
    r = F(z)
    for it in 0:maxit-1
        maximum(abs, r) ≤ tol && return (; u = collect(z), retcode = "Success", iter = it)
        Δ = jac_fd(F, z) \ r
        t = 1.0
        while true
            zn = z - t * Δ
            rn = F(zn)
            if maximum(abs, rn) < maximum(abs, r)
                z, r = zn, rn
                break
            end
            t /= 2
            t < 1e-6 && return (; u = collect(z), retcode = "Stalled", iter = it)
        end
    end
    (; u = collect(z), retcode = "MaxIters", iter = maxit)
end
