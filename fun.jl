function pressure_poisson(p, u, v, dt, dx, dy, ρ, nit)
    I, J = size(p)    
    pn = zeros(I, J)

    for _ = 1:nit
        for i = 2:I - 1
            for j = 2:J - 1
                pn[i, j] = (((p[i, j + 1] + p[i, j - 1]) * dx^2 + (p[i + 1, j] + p[i - 1, j]) * dy^2) / (2(dx^2 + dy^2)) -
                ((ρ * dx^2 * dy^2) / (2(dx^2 + dy^2))) * ((1 / dt) * ((u[i + 1, j] - u[i - 1, j]) / 2dx + (v[i, j + 1] - v[i, j - 1]) / 2dy) -
                ((u[i + 1, j] - u[i - 1, j]) / 2dx)^2 -
                2 * ((v[i + 1, j] - v[i - 1, j]) / 2dx) * ((u[i, j + 1] - u[i, j - 1]) / 2dy) -
                ((v[i,j + 1] - v[i,j - 1]) / 2dy)^2))
                
                if isnan(pn[i, j])
                    throw(ErrorException("NaN encountered"))
                end
            end
        end
    
        pn[:, end] .= 0
        pn[end, :] .= pn[end - 1, :]
        pn[:, 1] .= pn[:, 2]
        pn[1, :] .= pn[2, :]
        
        p = pn
    end

    return pn
end

function momentum_u(p, u, v, dt, dx, dy, ρ, ν)
    I, J = size(p)
    un = zeros(I, J)

    for i = 2:I - 1
        for j = 2:J - 1
            un[i, j] = (u[i, j] - 
            u[i, j] * (dt / dx) * (u[i, j] - u[i - 1, j]) -
            v[i, j] * (dt / dy) * (u[i, j] - u[i, j - 1]) -
            (dt / (2ρ * dx)) * (p[i + 1, j] - p[i - 1, j]) +
            ν * ((dt / dx^2) * (u[i + 1, j] - 2u[i, j] + u[i - 1, j]) + (dt / dy^2) * (u[i,j + 1] - 2u[i,j] + u[i, j - 1])))
        end
    end

    return un
end

function momentum_v(p, u, v, dt, dx, dy, ρ, ν)
    I, J = size(p)
    vn = zeros(I, J)
    
    for i = 2:I - 1
        for j = 2:J - 1
            vn[i, j] = (v[i, j] - 
            u[i, j] * (dt / dx) * (v[i, j] - v[i - 1, j]) -
            v[i, j] * (dt / dy) * (v[i, j] - v[i, j - 1]) -
            (dt / (2ρ * dy)) * (p[i, j + 1] - p[i, j - 1]) +
            ν * ((dt / dx^2) * (v[i + 1, j] - 2v[i, j] + v[i - 1, j]) + (dt / dy^2) * (v[i,j + 1] - 2v[i,j] + v[i, j - 1])))
        end
    end

    return vn
end

function cavity_(p, u, v, dt, dx, dy, ρ, ν, nit)
    u_, v_ = copy(u), copy(v)   

    p = pressure_poisson(p, u_, v_, dt, dx, dy, ρ, nit)
    u = momentum_u(p, u_, v_, dt, dx, dy, ρ, ν)
    v = momentum_v(p, u_, v_, dt, dx, dy, ρ, ν)
    
    u[1, :] .= u[:, 1] .= u[end, :] .= 0
    v[1, :] .= v[:, 1] .= v[end, :] .= v[:, end] .= 0
    u[:, end] .= 1

    return p, u, v
end

function cavity(p, u, v, dt, dx, dy, ρ, ν, nit, nt)
    for i = 1:nt
        u_, v_ = copy(u), copy(v)   

        p = pressure_poisson(p, u_, v_, dt, dx, dy, ρ, nit)
        u = momentum_u(p, u_, v_, dt, dx, dy, ρ, ν)
        v = momentum_v(p, u_, v_, dt, dx, dy, ρ, ν)
    
        u[1, :] .= u[:, 1] .= u[end, :] .= 0
        v[1, :] .= v[:, 1] .= v[end, :] .= v[:, end] .= 0
        u[:, end] .= 1
    end

    return u, v, p
end
