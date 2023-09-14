##
using Plots
using Statistics
##

##
δ = 1e-6;
x = [0.0, 3.0, 6.0]
y = zeros(3)
a = Animation()
#ymin = minimum(fvals);
#ymax = maximum(fvals);
ymin = -101;
ymax = 101;
for i=1:1:length(times)
    s = spins[i]
    u = real.(s)
    u[abs.(u) .< δ] .= zero(eltype(u))
    v = imag.(s)
    v[abs.(v) .< δ] .= zero(eltype(v))
    uavg = [mean(u)]
    vavg = [mean(v)]
    running_t = times[1:i]
    running_f = fvals[1:i]
    l = @layout [a b ; c d e]
    p1 = plot(running_t, running_f, xlims=(0, times[end]), ylims=(ymin, ymax),
            label=false, title="Z₁ driving amplitude",
            xlabel="t", ylabel="h₁(t)")
    p2 = quiver([0], [0], quiver=(uavg, vavg), xlims=(-2, 2), ylims=(-2, 2), 
                title="(ψ₁+ψ₂+ψ₃) / 3", c=:orange, xlabel="Re",
                ylabel="Im")
    p3 = scatter([x[1]], [y[1]], xlims=(-2, 2), ylims=(-2, 2), 
                legend=false, title="ψ₁",
                markersize=0, markershape=:utriangle,
                xlabel="Re", ylabel="Im")
    p3 = quiver!([x[1]], [y[1]], quiver=([u[1]],[v[1]]));
    """
    p4 = scatter([x[1]], [y[1]], xlims=(-2, 2), ylims=(-2, 2), 
                legend=false, title="ψ₂",
                markersize=0, markershape=:utriangle,
                xlabel="Re", ylabel="Im")
    p4 = quiver!([x[1]], [y[1]], quiver=([u[2]],[v[2]]));
    p5 = scatter([x[1]], [y[1]], xlims=(-2, 2), ylims=(-2, 2), 
                legend=false, title="ψ₃",
                markersize=0, markershape=:utriangle,
                xlabel="Re", ylabel="Im")
    p5 = quiver!([x[1]], [y[1]], quiver=([u[3]],[v[3]]));
    plt = plot(p1, p2, p3, p4, p5, layout = l)
    """
    plt = plot(p1, p2, p3)
    frame(a, plt)
end

gif(a)
##
 
##
l = @layout [
    a{0.3w} [grid(3,3)
             b{0.2h}  ]
]
plot(
    rand(10, 11),
    layout = l, legend = false, seriestype = [:bar :scatter :path],
    title = ["($i)" for j in 1:1, i in 1:11], titleloc = :right, titlefont = font(8)
)
##

##
i = 25
s = spins[i]
u = real.(s)
u[abs.(u) .< δ] .= zero(eltype(u))
v = imag.(s)
v[abs.(v) .< δ] .= zero(eltype(v))
xavg = [7.5]
yavg = [2]
uavg = [mean(u)]
vavg = [mean(v)]
running_t = times[1:i]
running_f = fvals[1:i]

scatter(x, y, xlims=(-3, 9), ylims=(-2, 4), legend=false)
quiver!(x, y, quiver=(u,v));
quiver!(xavg, yavg, quiver=(uavg, vavg))
plt = plot!(running_t, running_f, inset = (1, bbox(0.05, 0.05, 0.5, 0.25, (0, 1))), subplot=2,
                xlims=(0, times[end]), ylims=(ymin, ymax), 
                legend=false, xlabel="time", ylabel="h₁(t)")
##

##
l = @layout [a b ; c d e]
p1 = plot(running_t, running_f, xlims=(0, times[end]), ylims=(ymin, ymax),
         label=false, title="Z₁ driving amplitude",
         xlabel="t", ylabel="h₁(t)")
p2 = quiver([0], [0], quiver=(uavg, vavg), xlims=(-2, 2), ylims=(-2, 2), 
            title="(ψ₁+ψ₂+ψ₃) / 3", c=:orange, xlabel="Re",
            ylabel="Im")
p3 = scatter([x[1]], [y[1]], xlims=(-2, 2), ylims=(-2, 2), 
            legend=false, title="ψ₁",
            markersize=0, markershape=:utriangle,
            xlabel="Re", ylabel="Im")
p3 = quiver!([x[1]], [y[1]], quiver=([u[1]],[v[1]]));
p4 = scatter([x[1]], [y[1]], xlims=(-2, 2), ylims=(-2, 2), 
            legend=false, title="ψ₂",
            markersize=0, markershape=:utriangle,
            xlabel="Re", ylabel="Im")
p4 = quiver!([x[1]], [y[1]], quiver=([u[2]],[v[2]]));
p5 = scatter([x[1]], [y[1]], xlims=(-2, 2), ylims=(-2, 2), 
            legend=false, title="ψ₃",
            markersize=0, markershape=:utriangle,
            xlabel="Re", ylabel="Im")
p5 = quiver!([x[1]], [y[1]], quiver=([u[3]],[v[3]]));
plot(p1, p2, p3, p4, p5, layout = l)
##


##
##
x = real.(spins)
y = real.(spins)
GR.setarrowsize(1)
plot(x,y, marker =:circle, arrow=(:closed, 2.0))
##

##
x = [0.0, 3.0, 6.0]
y = zeros(3)
u = real.(spins)
u[abs.(u) .< δ] .= zero(eltype(u))
v = imag.(spins)
v[abs.(v) .< δ] .= zero(eltype(v))
scatter(x,y, xlims=(-3, 9), ylims=(-2, 2), legend=false)
quiver!(x,y,quiver=(u,v))
##
##