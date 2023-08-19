# using CairoMakie
using GLMakie
using LinearAlgebra

function r(i, j)
    q = 1.0

    a1 = 2pi .* [1.0, 1 / sqrt(3)] ./ q
    a2 = 2pi .* [0.0, 2 / sqrt(3)] ./ q

    a1 ./= norm(a1)
    a2 ./= norm(a2)

    return i .* a1 .+ j .* a2 .- a1 ./ 2 .- a2 ./ 2
end

function rPolar(r, theta)
    q = 1.0

    a1 = 2pi .* [1.0, 1 / sqrt(3)] ./ q
    a2 = 2pi .* [0.0, 2 / sqrt(3)] ./ q

    a1 ./= norm(a1)
    a2 ./= norm(a2)

    r0 = [
        r * cos(theta),
        r * sin(theta),
    ]

    return r0 # .- a1 ./ 2 .- a2 ./ 2
end

function rPolar(r, θ, N)
    q = 1.0

    a1 = 2pi .* [1.0, 1 / sqrt(3)] ./ q
    a2 = 2pi .* [0.0, 2 / sqrt(3)] ./ q

    # a1 ./= norm(a1)
    # a2 ./= norm(a2)

    r0 = [
        cos(θ),
        sin(θ),
    ]

    # r0 ./= norm(r0)
    # r0 .*= (norm(a1) / 2 * (r / N))
    r0 .*= (norm(a2 ./ 2) * (r / N))

    return r0 # .- a1 ./ 2 .- a2 ./ 2
end

function map_to_sphere(points::Vector)
    # マップされた3次元データ点の座標を保存するための行列を作成
    sphere_points = zeros(eltype(points), size(points))

    # display(hcat(points[:]...)[1, :])

    max = maximum(hcat(points[:]...)[1, :])

    for i in 1:length(points[:])
        x = points[i][1] / max
        y = points[i][2] / max

        # 2次元データ点を3次元球面上にマップする
        r = sqrt(x^2 + y^2)
        # theta = atan(r)
        theta = π * r
        phi = atan(y, x)

        sphere_points[i] = Point3f(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta))
        # sphere_points[i] = sin(theta) * sin(phi)
        # sphere_points[i] = cos(theta)
    end

    return sphere_points
end

function S(i, j, m, phi)
    q = 1.0
    q1 = [q, 0.0]
    q2 = [-q / 2, sqrt(3) * q / 2]
    q3 = [-q / 2, -sqrt(3) * q / 2]

    a1 = 2pi .* [1.0, 1 / sqrt(3)] ./ q
    a2 = 2pi .* [0.0, 2 / sqrt(3)] ./ q

    r0 = r(i, j)

    Q1 = q1 ⋅ r0 + phi[1]
    Q2 = q2 ⋅ r0 + phi[2]
    Q3 = q3 ⋅ r0 + phi[3]

    # psiC = 1/sqrt(3)
    # psiS = 1/sqrt(3)

    s = [
        sqrt(3) / 2 * (sin(Q2) - sin(Q3))
        -sin(Q1) + 1 / 2 * (sin(Q2) + sin(Q3))
        cos(Q1) + cos(Q2) + cos(Q3) + sqrt(3) * m
    ]

    return s ./ norm(s)
end

function S(r, m, phi)
    q = 1.0
    q1 = [q, 0.0]
    q2 = [-q / 2, sqrt(3) * q / 2]
    q3 = [-q / 2, -sqrt(3) * q / 2]

    a1 = 2pi .* [1.0, 1 / sqrt(3)] ./ q
    a2 = 2pi .* [0.0, 2 / sqrt(3)] ./ q

    # r0 = r(i, j)

    Q1 = q1 ⋅ r + phi[1]
    Q2 = q2 ⋅ r + phi[2]
    Q3 = q3 ⋅ r + phi[3]

    # psiC = 1/sqrt(3)
    # psiS = 1/sqrt(3)

    s = [
        sqrt(3) / 2 * (sin(Q2) - sin(Q3))
        -sin(Q1) + 1 / 2 * (sin(Q2) + sin(Q3))
        cos(Q1) + cos(Q2) + cos(Q3) + sqrt(3) * m
    ]

    return s ./ norm(s)
end

function SPolar(r0, m, phi, theta)
    q = 1.0
    q1 = [q, 0.0]
    q2 = [-q / 2, sqrt(3) * q / 2]
    q3 = [-q / 2, -sqrt(3) * q / 2]

    a1 = 2pi .* [1.0, 1 / sqrt(3)] ./ q
    a2 = 2pi .* [0.0, 2 / sqrt(3)] ./ q

    # r0 = r(i, j)
    # r0 = [
    #     r * cos(theta),
    #     r * sin(theta),
    # ]
    # r0 = rPolar(r, theta)

    Q1 = q1 ⋅ r0 + phi[1]
    Q2 = q2 ⋅ r0 + phi[2]
    Q3 = q3 ⋅ r0 + phi[3]

    # psiC = 1/sqrt(3)
    # psiS = 1/sqrt(3)

    s = [
        sqrt(3) / 2 * (sin(Q2) - sin(Q3))
        -sin(Q1) + 1 / 2 * (sin(Q2) + sin(Q3))
        cos(Q1) + cos(Q2) + cos(Q3) + sqrt(3) * m
    ]
    A = [cos(theta) -sin(theta); sin(theta) cos(theta)]

    s[1:2] .= A * s[1:2]

    return s ./ norm(s)
end

function SPolar(r0, m, phi)
    q = 1.0
    q1 = [q, 0.0]
    q2 = [-q / 2, sqrt(3) * q / 2]
    q3 = [-q / 2, -sqrt(3) * q / 2]

    a1 = 2pi .* [1.0, 1 / sqrt(3)] ./ q
    a2 = 2pi .* [0.0, 2 / sqrt(3)] ./ q

    # r0 = r(i, j)
    # r0 = [
    #     r * cos(theta),
    #     r * sin(theta),
    # ]
    # r0 = rPolar(r, theta)

    Q1 = q1 ⋅ r0 + phi[1]
    Q2 = q2 ⋅ r0 + phi[2]
    Q3 = q3 ⋅ r0 + phi[3]

    # psiC = 1/sqrt(3)
    # psiS = 1/sqrt(3)

    s = [
        sqrt(3) / 2 * (sin(Q2) - sin(Q3))
        -sin(Q1) + 1 / 2 * (sin(Q2) + sin(Q3))
        cos(Q1) + cos(Q2) + cos(Q3) + sqrt(3) * m
    ]

    return s ./ norm(s)
end

function writeSpinTexture(p)

    (; Nₛ, m, phi) = p
    # idx(i, j, k) = 3Nₛ * (i - 1) + 3(j - 1) .+ (k)

    nor = sqrt(16pi^2 / 3)

    # r = [i, j]
    ra = 0.5
    cut = 0.5

    ps = [Point3f(r(ra * i, ra * j)[1], r(ra * i, ra * j)[2], 0.0) for i ∈ -Nₛ:Nₛ for j ∈ -Nₛ:Nₛ] .* 2 ./ Nₛ
    # ps = [Point3f(i, j, 0) for i ∈ 1:Nₛ for j ∈ 1:Nₛ] ./ Nₛ
    # append!(ps, [Point3f(0, 0, -10), Point3f(0, 0, -10)])
    append!(ps, [Point3f(10.0, 10.0, -10.0), Point3f(10.0, 10.0, -10.0)])
    # ns = [Vec3f(S(i/10, j/10, m, phi)) for i ∈ 1:Nₛ for j ∈ 1:Nₛ] ./ Nₛ
    ns = [Vec3f(S(r(ra * i, ra * j), m, phi)) for i ∈ -Nₛ:Nₛ for j ∈ -Nₛ:Nₛ] .* 0.6 ./ Nₛ
    append!(ns, [Vec3f(0.0, 0.0, 1 / (2Nₛ)), Vec3f(0.0, 0.0, -1 / (2Nₛ))])

    x = [j for i ∈ -Nₛ/2:Nₛ/2 for j ∈ -Nₛ/2:Nₛ/2] ./ Nₛ
    y = [i for i ∈ -Nₛ/2:Nₛ/2 for j ∈ -Nₛ/2:Nₛ/2] ./ Nₛ
    # x = [r(i, j)[1] for i ∈ 1:Nₛ for j ∈ 1:Nₛ] ./ (nor*Nₛ)
    # y = [r(i, j)[2] for i ∈ 1:Nₛ for j ∈ 1:Nₛ] ./ (nor*Nₛ)
    # z = [df[1, idx(i, j, 3)+1] for i ∈ 1:Nₛ, for j ∈ 1:Nₛ]

    fig = Figure(resolution=(600 * 2, 550 * 2))
    ax = Axis3(fig[1, 1],
        # title="t = 0.0",
        # titlegap=-10,
        protrusions=(0, 0, 0, 0),
        viewmode=:fitzoom,
        limits=(-0.5(1 + 1 / Nₛ), 0.5(1 + 1 / Nₛ), -0.5(1 + 1 / Nₛ), 0.5(1 + 1 / Nₛ), 0, 0.8),
        # limits=(0, 1.0, 0, 1.0, 0, 0.8),
        # elevation=π / 2,
        elevation=π / 6,
        azimuth=0,
        # azimuth=π / 6,
        xypanelcolor=(:gray, 0.2)
    )

    lengths = [S([i / 2, j / 2], m, phi)[3] for i ∈ -Nₛ/2:Nₛ/2 for j ∈ -Nₛ/2:Nₛ/2]# ./ 2
    hm = heatmap!(
        ax, x, y, lengths;
        colormap=(:darkrainbow, 0.4),
        colorrange=(-1, 1),
        interpolate=true
    )

    lengths = [S(r(ra * i, ra * j), m, phi)[3] for i ∈ -Nₛ:Nₛ for j ∈ -Nₛ:Nₛ]
    append!(lengths, [1.0, -1.0])
    ar = arrows!(
        ax, ps[norm.(ps).<cut], ns[norm.(ps).<cut], fxaa=true, # turn on anti-aliasing
        color=lengths[norm.(ps).<cut],
        colormap=:darkrainbow,
        linewidth=0.1 / Nₛ, arrowsize=Vec3f(0.2, 0.2, 0.3) ./ Nₛ,
        align=:center,
        quality=128
    )

    Colorbar(fig[:, end+1], hm)

    hidedecorations!(ax)
    hidespines!(ax)

    fig
end

function writeSpinTexturePolar(p)

    (; Nₛ, m, phi) = p
    # idx(i, j, k) = 3Nₛ * (i - 1) + 3(j - 1) .+ (k)

    nor = sqrt(16pi^2 / 3)

    # r = [i, j]
    # ra = 1.65
    a1 = 2pi .* [1.0, 1 / sqrt(3)]
    ra = 1.0
    cut = 1.0
    Nr = 13
    f = 5

    ps = [Point3f(rPolar(r, theta, Nr)[1], rPolar(r, theta, Nr)[2], 0.0) for r ∈ 1:Nr÷2 for theta ∈ range(0, 2pi, length=f * r)] ./ norm(a1) # .* ra / 2 ./ Nₛ
    append!(ps, [Point3f(rPolar(r, theta, Nr)[1], rPolar(r, theta, Nr)[2], 0.0) for r ∈ Nr÷2+1:Nr-1 for theta ∈ range(0, 2pi, length=f * (Nr - r))] ./ norm(a1)) # .* ra / 2 ./ Nₛ)

    append!(ps, [Point3f(rPolar(0, 0, Nr)[1], rPolar(0, 0, Nr)[2], 0.0), Point3f(rPolar(Nr, 0, Nr)[1], rPolar(Nr, 0, Nr)[2], 0.0) ./ norm(a1)])# * ra / 2 / Nₛ])
    sphere_points = map_to_sphere(ps)

    ns = [Vec3f(SPolar(rPolar(r, theta, Nr), m, phi)) for r ∈ 1:Nr÷2 for theta ∈ range(0, 2pi, length=f * r)] .* 0.6 ./ Nₛ
    append!(ns, [Vec3f(SPolar(rPolar(r, theta, Nr), m, phi)) for r ∈ Nr÷2+1:Nr-1 for theta ∈ range(0, 2pi, length=f * (Nr - r))] .* 0.6 ./ Nₛ)
    append!(ns, [Vec3f(SPolar(rPolar(0, 0, Nr), m, phi)), Vec3f(SPolar(rPolar(Nr, 0, Nr), m, phi))] .* 0.6 ./ Nₛ)

    # ps = [Point3f(rPolar(r, theta)[1], rPolar(r, theta)[2], 0.0) for r ∈ 1:Nr÷2 for theta ∈ range(pi, 2pi, length=3r)] .* ra / 2 ./ Nₛ
    # append!(ps, [Point3f(rPolar(r, theta)[1], rPolar(r, theta)[2], 0.0) for r ∈ Nr÷2+1:Nr-1 for theta ∈ range(pi, 2pi, length=3(Nr - r))] .* ra / 2 ./ Nₛ)

    # append!(ps, [Point3f(rPolar(0, 0)[1], rPolar(0, 0)[2], 0.0), Point3f(rPolar(Nr, 0)[1], rPolar(Nr, 0)[2], 0.0) * ra / 2 / Nₛ])
    # sphere_points = map_to_sphere(ps)

    # ns = [Vec3f(SPolar(rPolar(r / ra, theta), m, phi)) for r ∈ 1:Nr÷2 for theta ∈ range(pi, 2pi, length=3r)] .* 0.6 ./ Nₛ
    # append!(ns, [Vec3f(SPolar(rPolar(r / ra, theta), m, phi)) for r ∈ Nr÷2+1:Nr-1 for theta ∈ range(pi, 2pi, length=3(Nr - r))] .* 0.6 ./ Nₛ)
    # append!(ns, [Vec3f(SPolar(rPolar(0, 0), m, phi)), Vec3f(SPolar(rPolar(Nr, 0), m, phi))] .* 0.6 ./ Nₛ)

    r = 0.85
    θ = range(0.0, pi, length=100)
    φ = range(0.0, 2pi, length=100)
    x = [r * sin(θ) * cos(φ) for θ in θ, φ in φ]
    y = [r * sin(θ) * sin(φ) for θ in θ, φ in φ]
    z = [r * cos(θ) for θ in θ, φ in φ]

    # r = range(-1.0, 1.0, length=100)
    # cube = [(x .^ 2 + y .^ 2 + z .^ 2) for x = r, y = r, z = r]
    # cube = cube .* (cube .> 1.0)
    # # cube = cube[norm.(cube).<1.0]
    # display(cube)

    # x = [r(i, j)[1] for i ∈ 1:Nₛ for j ∈ 1:Nₛ] ./ (nor*Nₛ)
    # y = [r(i, j)[2] for i ∈ 1:Nₛ for j ∈ 1:Nₛ] ./ (nor*Nₛ)
    # z = [df[1, idx(i, j, 3)+1] for i ∈ 1:Nₛ, for j ∈ 1:Nₛ]

    fig = Figure(resolution=(600 * 2, 600 * 2))
    ax = Axis3(fig[1, 1],
        # title="t = 0.0",
        # titlegap=-10,
        aspect=(1, 1, 1),
        protrusions=(0, 0, 0, 0),
        viewmode=:fitzoom,
        # viewmode=:stretch,
        limits=(-1(1 + 1 / Nₛ), 1(1 + 1 / Nₛ), -1(1 + 1 / Nₛ), 1(1 + 1 / Nₛ), -1.2, 1.2),
        # limits=(0, 1.0, 0, 1.0, 0, 0.8),
        # elevation=π / 2,
        # elevation=π / 6,
        elevation=π / 20,
        # elevation=0,
        # azimuth=0,
        # azimuth=π / 6,
        azimuth=π / 2,
        # xypanelcolor=(:gray, 0.2)
    )

    # ratio = 5

    # lengths = [SPolar(rPolar(r / ra, theta), m, phi)[3] for r ∈ 1:Nr÷2 for theta ∈ range(pi, 2pi, length=3r)]# ./ 2
    # append!(lengths, [SPolar(rPolar(r / ra, theta), m, phi)[3] for r ∈ Nr÷2+1:Nr-1 for theta ∈ range(pi, 2pi, length=3(Nr - r))])
    # # append!(lengths, [1.0, -1.0])
    # append!(lengths, [SPolar(rPolar(0, 0), m, phi)[3], SPolar(rPolar(Nr, 0), m, phi)[3]])
    # ar = arrows!(
    #     ax, sphere_points[norm.(ps).<cut], 4ns[norm.(ps).<cut], fxaa=true, # turn on anti-aliasing
    #     color=lengths[norm.(ps).<cut],
    #     colormap=:lightrainbow,
    #     linewidth=0.5 / Nₛ, arrowsize=Vec3f(0.4ratio, 0.4ratio, 0.6ratio) ./ Nₛ,
    #     align=:center,
    #     quality=128
    # )

    lengths = [S([i / 128, j / 128], m, phi)[3] for i ∈ -32Nₛ:32Nₛ for j ∈ -32Nₛ:32Nₛ]# ./ 2
    # lengths[sqrt.(x .^ 2 .+ y .^ 2).>cut] .= -2
    hm = surface!(
        ax, x, y, z;
        # ax, cube;
        # ax, x, y, lengths;
        colormap=:grays,
        colorrange=(-1, 1)#,
        # lowclip=:white,
        # interpolate=true
    )
    # volume!(ax, cube)

    # ps = [Point3f(rPolar(r, theta)[1], rPolar(r, theta)[2], 0.0) for r ∈ 1:Nr÷2 for theta ∈ range(0.0, pi, length=3r)] .* ra / 2 ./ Nₛ
    # append!(ps, [Point3f(rPolar(r, theta)[1], rPolar(r, theta)[2], 0.0) for r ∈ Nr÷2+1:Nr-1 for theta ∈ range(0.0, pi, length=3(Nr - r))] .* ra / 2 ./ Nₛ)

    # append!(ps, [Point3f(rPolar(0, 0)[1], rPolar(0, 0)[2], 0.0), Point3f(rPolar(Nr, 0)[1], rPolar(Nr, 0)[2], 0.0) * ra / 2 / Nₛ])
    # sphere_points = map_to_sphere(ps)

    # ns = [Vec3f(SPolar(rPolar(r / ra, theta), m, phi)) for r ∈ 1:Nr÷2 for theta ∈ range(0.0, pi, length=3r)] .* 0.6 ./ Nₛ
    # append!(ns, [Vec3f(SPolar(rPolar(r / ra, theta), m, phi)) for r ∈ Nr÷2+1:Nr-1 for theta ∈ range(0.0, pi, length=3(Nr - r))] .* 0.6 ./ Nₛ)
    # append!(ns, [Vec3f(SPolar(rPolar(0, 0), m, phi)), Vec3f(SPolar(rPolar(Nr, 0), m, phi))] .* 0.6 ./ Nₛ)

    ratio = 6

    lengths = [SPolar(rPolar(r, theta, Nr), m, phi)[3] for r ∈ 1:Nr÷2 for theta ∈ range(0.0, 2pi, length=f * r)]# ./ 2
    append!(lengths, [SPolar(rPolar(r, theta, Nr), m, phi)[3] for r ∈ Nr÷2+1:Nr-1 for theta ∈ range(0.0, 2pi, length=f * (Nr - r))])
    # append!(lengths, [1.0, -1.0])
    append!(lengths, [SPolar(rPolar(0, 0, Nr), m, phi)[3], SPolar(rPolar(Nr, 0, Nr), m, phi)[3]])
    ar = arrows!(
        ax, sphere_points[norm.(ps).<cut], 6ns[norm.(ps).<cut], fxaa=true, # turn on anti-aliasing
        color=lengths[norm.(ps).<cut],
        colormap=:jet1,
        linewidth=1.0 / Nₛ, arrowsize=Vec3f(0.4ratio, 0.4ratio, 0.6ratio) ./ Nₛ,
        align=:center,
        quality=128
    )

    # Colorbar(fig[:, end+1], hm)

    hidedecorations!(ax)
    hidespines!(ax)

    fig
end


function outputFigToFiles(fig, figName, p)
    # origin = "LLG_"

    # figName = "data/" * origin * parameters * ".pdf"

    # fig = writeFig(df11, df12, df21, df22, df3, p)
    save(figName, fig)
end

# 情報出力・ディレクトリ作成
function makeDirectories()

    # ### stdoutに情報を出力 ###
    # println("Julia source file -> " * file)
    # print("Execution host -> ")
    # flush(stdout)
    # run(`hostname`)
    # println("Parameter -> " * parameters)
    # flush(stdout)
    # ### ----------- ###

    ### ディレクトリの存在確認と作成 ###
    if !isdir("data")
        mkdir("data")
    end
    if !isdir("rawData")
        mkdir("rawData")
    end
    if !isdir("figure")
        mkdir("figure")
    end
    ### ----------- ###
end