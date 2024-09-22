using Gaugefields
using Random

function heatbathtest_4D(NX, NY, NZ, NT, β, NC)
    filefooter = "$(NX)$(NY)$(NZ)$(NT)_SU$(NC)_beta$(β).txt"
    polyfile = "hb_poly_" * filefooter
    plaqfile = "hb_plaq_" * filefooter
    fp_poly = open(polyfile, "w")
    fp_plaq = open(plaqfile, "w")


    Dim = 4
    Nwing = 0

    Random.seed!(123)
    U = Initialize_Gaugefields(NC, Nwing, NX, NY, NZ, NT, condition="hot", randomnumber="Reproducible")
    println(typeof(U))

    gauge_action = GaugeAction(U)
    plaqloop = make_loops_fromname("plaquette", Dim=Dim)
    append!(plaqloop, plaqloop')
    βinp = β / 2
    push!(gauge_action, βinp, plaqloop)

    hnew = Heatbath_update(U, gauge_action)

    show(gauge_action)

    temp1 = similar(U[1])
    temp2 = similar(U[1])
    temp3 = similar(U[1])

    comb = 6
    factor = 1 / (comb * U[1].NV * U[1].NC)
    @time plaq_t = calculate_Plaquette(U, temp1, temp2) * factor
    println("plaq_t = $plaq_t")
    poly = calculate_Polyakov_loop(U, temp1, temp2)
    println("polyakov loop = $(real(poly)) $(imag(poly))")

    numhb = 30000
    count = 0
    sumplaq = 0.0
    for itrj = 1:numhb

        heatbath!(U, hnew)

        plaq_t = calculate_Plaquette(U, temp1, temp2) * factor
        poly = calculate_Polyakov_loop(U, temp1, temp2)

        if itrj % 40 == 0
            println("$itrj plaq_t = $plaq_t")
            println("$itrj polyakov loop = $(real(poly)) $(imag(poly))")

        end
        if itrj > 1000
            count += 1
            sumplaq += plaq_t
            println(fp_plaq, "$itrj $plaq_t $(sumplaq/count)")
            println(fp_poly, "$itrj $(real(poly)) $(imag(poly))")
            flush(fp_poly)
            flush(fp_plaq)
        end

    end

    return plaq_t

end

NX = 6
NY = NX
NZ = NX
NT = NX
NC = 3
β = 5.7
#β = 6
heatbathtest_4D(NX, NY, NZ, NT, β, NC)