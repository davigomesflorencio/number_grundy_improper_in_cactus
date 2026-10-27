# Exemplos de uso e comparação das formulações. Execute com:
#     julia --project=. exemplos.jl

# evita redefinir o módulo (e a ambiguidade de nomes) ao reexecutar no mesmo REPL
isdefined(@__MODULE__, :GrundyImproprioCacto) || include(joinpath(@__DIR__, "src", "GrundyImproprioCacto.jl"))
using .GrundyImproprioCacto
using Random, Printf

function main()
    println("== 1. Uso básico ==")
    G = cacto_aleatorio(12; rng = MersenneTwister(1))
    println(G, ", ", num_ciclos(G), " ciclo(s), blocos: ", blocos(G))
    Γ, cor, info = grundy_improprio_cacto(G)                 # k = 1, representantes
    ordem = ordem_da_coloracao(cor)
    println("Γ_1 = $Γ  (ótimo? $(info.otimo), $(round(info.tempo, digits = 2)) s)")
    println("cores: $cor  →  ordem: $ordem")
    println("a heurística nessa ordem reproduz a coloração? ",
            coloracao_gulosa(G, ordem, 1) == cor)

    println("\n== 2. Γ_1 em famílias de cactos ==")
    familias = [("C_$L", ciclo(L)) for L in 3:8]
    append!(familias, [("F_$m (amizade)", amizade(m)) for m in 1:4])
    append!(familias, [("moinho(3, $L)", moinho(3, L)) for L in 4:5])
    append!(familias, [("cadeia($m, 4)", cadeia_de_ciclos(m, 4)) for m in 2:4])
    append!(familias, [("sol($L)", sol(L)) for L in 3:6])
    @printf("  %-16s %-4s %-3s %-5s %-5s\n", "grafo", "n", "Δ", "Γ_0", "Γ_1")
    for (nome, g) in familias
        Γ0 = grundy_pi_atribuicao(g, 0)[1]
        Γ1 = grundy_pi_atribuicao(g, 1)[1]
        @printf("  %-16s %-4d %-3d %-5d %-5d\n", nome, nv(g), grau_maximo(g), Γ0, Γ1)
    end

    println("\n== 3. Distribuição de Γ_1 nos cactos com n ≤ 8 ==")
    for n in 1:8
        cs = cactos_nao_isomorfos(n)
        Γs = [grundy_pi_atribuicao(g, 1)[1] for g in cs]
        dist = join(["Γ=$t: $(count(==(t), Γs))" for t in sort(unique(Γs))], ", ")
        println("  n = $n ($(length(cs)) cactos):  $dist")
    end

    println("\n== 4. Menor cacto com Γ_1 = 3 ==")
    for n in 1:8
        i = findfirst(g -> grundy_pi_atribuicao(g, 1)[1] == 3, cactos_nao_isomorfos(n))
        i === nothing && continue
        g = cactos_nao_isomorfos(n)[i]
        println("  n = $n: arestas ", arestas(g))
        break
    end

    println("\n== 5. Atribuição × representantes (k = 1, limite de 60 s) ==")
    grundy_pi_atribuicao(ciclo(5), 1); grundy_pi_representantes(ciclo(5), 1)   # compilação
    @printf("  %-4s %-4s %-4s %-12s %-12s\n", "n", "Δ", "Γ_1", "atrib. (s)", "repr. (s)")
    rng = MersenneTwister(10)
    for n in (10, 15, 20, 25)
        g = cacto_aleatorio(n; rng)
        a = grundy_pi_atribuicao(g, 1; tempo_limite = 60)
        r = grundy_pi_representantes(g, 1; tempo_limite = 60)
        marca(x) = x[3].otimo ? "" : "*"
        @printf("  %-4d %-4d %-4s %-12.2f %-12.2f\n", n, grau_maximo(g),
                "$(a[1])$(marca(a))", a[3].tempo, r[3].tempo)
    end
    println("  (* = limite de tempo atingido; o valor é um limitante inferior)")
end

main()
