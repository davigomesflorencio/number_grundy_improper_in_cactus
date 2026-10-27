# Testes de correção. Execute com:
#     julia --project=. testes.jl
# ou, de dentro do REPL, include("testes.jl").

# evita redefinir o módulo (e a ambiguidade de nomes) ao reexecutar no mesmo REPL
isdefined(@__MODULE__, :GrundyImproprioCacto) || include(joinpath(@__DIR__, "src", "GrundyImproprioCacto.jl"))
using .GrundyImproprioCacto
using Random, Test

# A000083: número de cactos não isomorfos com n vértices
const A000083 = [1, 1, 2, 4, 9, 23, 63, 188, 596, 1979]

const FORMULACOES = (grundy_pi_atribuicao, grundy_pi_representantes)

# a solução é ótima, igual a Γ, gulosa, e reproduzida pela heurística?
function confere(g, k, Γ, f)
    Γf, cor, info = f(g, k)
    return info.otimo && Γf == Γ && eh_gulosa(g, cor, k) &&
           coloracao_gulosa(g, ordem_da_coloracao(cor), k) == cor
end

@testset "Grundy k-impróprio em cactos" begin

    @testset "geradores e reconhecimento" begin
        for n in eachindex(A000083)
            cs = cactos_nao_isomorfos(n)
            @test length(cs) == A000083[n]
            @test all(eh_cacto, cs)
            @test length(unique(forma_canonica.(cs))) == length(cs)
        end
        rng = MersenneTwister(1)
        for _ in 1:200
            g = cacto_aleatorio(rand(rng, 1:60); p_ciclo = rand(rng), tam_max = rand(rng, 3:8), rng)
            @test eh_cacto(g)
            # a forma canônica não depende dos rótulos
            @test forma_canonica(renomeia(g, randperm(rng, nv(g)))) == forma_canonica(g)
        end
        @test eh_cacto(caminho(5)) && eh_cacto(estrela(6)) && eh_cacto(ciclo(7))
        @test eh_cacto(amizade(4)) && nv(amizade(4)) == 9 && num_ciclos(amizade(4)) == 4
        @test eh_cacto(moinho(3, 5)) && nv(moinho(3, 5)) == 13
        @test eh_cacto(cadeia_de_ciclos(4, 6)) && num_ciclos(cadeia_de_ciclos(4, 6)) == 4
        @test eh_cacto(sol(5)) && nv(sol(5)) == 10
        @test !eh_cacto(completo(4))                                   # bloco K4
        @test !eh_cacto(grafo_de_arestas(4, [(1, 2), (2, 3), (3, 4), (4, 1), (1, 3)]))  # dois ciclos com aresta comum
        @test !eh_cacto(Grafo(2))                                      # desconexo
    end

    @testset "valores conhecidos" begin
        @test grundy_forca_bruta(completo(5), 1)[1] == 3
        @test grundy_forca_bruta(ciclo(6), 1)[1] == 2
        for f in FORMULACOES
            @test f(completo(5), 1)[1] == 3
            @test f(ciclo(6), 1)[1] == 2
            @test [f(caminho(n), 1)[1] for n in 1:6] == [1, 1, 2, 2, 2, 2]
        end
    end

    @testset "formulações = força bruta, todos os cactos com n ≤ 7, k = 0..2" begin
        for n in 1:7, g in cactos_nao_isomorfos(n), k in 0:2
            Γ = grundy_forca_bruta(g, k)[1]
            for f in FORMULACOES
                @test confere(g, k, Γ, f)
            end
        end
    end

    @testset "formulações = força bruta, cactos com n = 8, k = 1" begin
        for g in cactos_nao_isomorfos(8)
            Γ = grundy_forca_bruta(g, 1)[1]
            for f in FORMULACOES
                @test confere(g, 1, Γ, f)
            end
        end
    end

    @testset "formulações = força bruta em grafos gerais" begin
        rng = MersenneTwister(2)
        for _ in 1:100
            g = grafo_aleatorio(rand(rng, 1:7), rand(rng); rng)
            k = rand(rng, 0:3)
            Γ = grundy_forca_bruta(g, k)[1]
            for f in FORMULACOES
                @test confere(g, k, Γ, f)
            end
        end
    end

    @testset "busca de certificados = força bruta (cactos n ≤ 8, grafos gerais)" begin
        valido(g, k, Γ, cor) = maximum(coloracao_gulosa(g, ordem_do_certificado(cor), k)) >= Γ
        for n in 1:8, g in cactos_nao_isomorfos(n), k in 0:2
            Γ, cor = grundy_certificado(g, k)
            @test Γ == grundy_forca_bruta(g, k)[1] && valido(g, k, Γ, cor)
        end
        rng = MersenneTwister(4)
        for _ in 1:200
            g = grafo_aleatorio(rand(rng, 1:8), rand(rng); rng)
            k = rand(rng, 0:3)
            Γ, cor = grundy_certificado(g, k)
            @test Γ == grundy_forca_bruta(g, k)[1] && valido(g, k, Γ, cor)
        end
    end

    @testset "busca de certificados = PI em cactos aleatórios" begin
        rng = MersenneTwister(5)
        for _ in 1:20
            g = cacto_aleatorio(rand(rng, 11:18); p_ciclo = rand(rng), rng)
            @test grundy_certificado(g, 1)[1] == grundy_pi_atribuicao(g, 1)[1]
        end
    end

    @testset "atribuição = representantes em cactos aleatórios" begin
        rng = MersenneTwister(3)
        for _ in 1:20
            g = cacto_aleatorio(rand(rng, 9:14); rng)
            a = grundy_improprio_cacto(g; formulacao = :atribuicao)
            r = grundy_improprio_cacto(g; formulacao = :representantes)
            @test a[3].otimo && r[3].otimo && a[1] == r[1]
            @test a[1] <= grau_maximo(g) + 1
        end
        @test_throws ArgumentError grundy_improprio_cacto(completo(4))
    end
end
