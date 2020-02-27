using PyPlot, Printf

include("fun.jl")
include("anim.jl")

nx, ny = 101, 101
nt = 30000
nit = 100

dx, dy = 1 ./ ([nx, ny] .- 1) 
x, y = collect(range(0, 1, length = nx)), collect(range(0, 1, length = ny))

ρ = 1
ν = 0.0025
dt = 0.0005

u = v = p = zeros(nx, ny)

figure(figsize = (10, 8))
axis("equal")

mkdir("results")

a = 0

for (i, t) = enumerate(range(0, length=nt, step=dt))
    global p, u, v = cavity_(p, u, v, dt, dx, dy, ρ, ν, nit)
    if (i-1) % 50 == 0
        println("Writing frame for t = $t / $(nt*dt)")
        u_, v_, p_ = transpose(u), transpose(v), transpose(p)

        vel = (u_ .^ 2 + v_ .^ 2) .^ 0.5

        contourf(x, y, vel, levels=collect(0:0.001:1.0), cmap="rainbow")
        cbar=colorbar(label="Velocity Magnitude")
        cbar.set_ticks(collect(0:0.1:1))

        quiver(x[1:5:end], y[1:5:end], u_[1:5:end,1:5:end], v_[1:5:end,1:5:end], color="white", scale=20, headlength=4, headaxislength=4)

        xlim((0, 1))
        ylim((0, 1))

        title("Lid Driven Cavity Flow \$Re=400\$")
        
        xlabel("\$x\$")
        ylabel("\$y\$")
        text(0.02, 0.02, @sprintf("\$t=%.2f\$", t), ha = "left", va = "bottom", color="w")

        savefig("results/$a.png", dpi = 300)

        clf()

        global a += 1
    end
    if (i-1) % 100 == 0
        println("$t/$(dt * nt)")
    end
end

u_, v_, p_ = transpose(u), transpose(v), transpose(p)
vel = (u_ .^ 2 + v_ .^ 2) .^ 0.5
contourf(x, y, vel, levels = collect(0:0.001:1.0), cmap = "rainbow")
cbar = colorbar(label="Velocity Magnitude")
cbar.set_ticks(collect(0:0.1:1))

streamplot(x, y, u_, v_, arrowstyle="->", color="white", linewidth=0.5)

xlim((0, 1))
ylim((0, 1))

title("Lid Driven Cavity Flow \$Re=400\$")

xlabel("\$x\$")
ylabel("\$y\$")
savefig("result.png", dpi=300)

generate_gif("results/%d.png", "result.gif", 60)
rm("results", recursive = true)