# Menor cacto com Γ_k ≥ t. Execute com:
#     julia --project=. menor_cacto.jl [t] [k] [nmax]
# (padrão: t = 4, k = 1, nmax = 16)
#
# Como Γ_k é monótono em subgrafos induzidos e todo subgrafo induzido conexo
# de um cacto é um cacto, basta percorrer os cactos não isomorfos por número
# de vértices e parar no primeiro n em que algum tem Γ_k ≥ t. Para cada n,
# o script decide Γ_k ≥ t em todos os cactos com n vértices pela busca de
# certificados. Os cactos ficam guardados como listas de arestas compactas,
# e o último nível não é guardado (só as formas canônicas, para descartar
# isomorfos).

isdefined(@__MODULE__, :GrundyImproprioCacto) || include(joinpath(@__DIR__, "src", "GrundyImproprioCacto.jl"))
using .GrundyImproprioCacto
using Printf

const Arestas = Vector{Tuple{UInt8,UInt8}}

compacta(g::Grafo) = Arestas([(UInt8(u), UInt8(v)) for (u, v) in arestas(g)])
expande(n::Int, E::Arestas) = grafo_de_arestas(n, [(Int(u), Int(v)) for (u, v) in E])

# chama f(h) para cada cacto com m vértices, sem repetir classes de isomorfismo
function cada_cacto(f, nivel::Vector{Vector{Arestas}}, m::Int)
    vistos = Set{String}()
    considera(h) = (c = forma_canonica(h); c in vistos || (push!(vistos, c); f(h)))
    for E in nivel[m-1]
        g = expande(m - 1, E)
        for v in 1:m-1
            h = Grafo([copy(a) for a in g.adj]); push!(h.adj, Int[])
            adiciona_aresta!(h, v, m)
            considera(h)
        end
    end
    for L in 3:m, E in nivel[m-L+1]
        g = expande(m - L + 1, E)
        for v in 1:m-L+1
            h = Grafo([copy(a) for a in g.adj])
            prox = nv(h) + 1
            append!(h.adj, [Int[] for _ in 1:L-1])
            ant = v
            for w in prox:prox+L-2
                adiciona_aresta!(h, ant, w); ant = w
            end
            adiciona_aresta!(h, ant, v)
            considera(h)
        end
    end
    return length(vistos)
end

function main(t::Int, k::Int, nmax::Int)
    println("Menor cacto com Γ_$k ≥ $t (até n = $nmax)")
    @printf("  %-4s %-10s %-10s %s\n", "n", "cactos", "Γ ≥ $t", "tempo (s)")
    nivel = Vector{Vector{Arestas}}(undef, nmax)
    nivel[1] = [Arestas()]
    for m in 2:nmax
        guarda = m < nmax
        nivel[m] = Arestas[]
        achados = Grafo[]
        tempo = @elapsed total = cada_cacto(nivel, m) do h
            guarda && push!(nivel[m], compacta(h))
            certificado_grundy(h, k, t) === nothing || push!(achados, h)
        end
        @printf("  %-4d %-10d %-10d %.1f\n", m, total, length(achados), tempo)
        flush(stdout)
        if !isempty(achados)
            println("\nMenor n: $m. $(length(achados)) cacto(s) não isomorfo(s) com Γ_$k ≥ $t:")
            for h in achados
                Γ, cor, info = grundy_pi_atribuicao(h, k)
                println("  arestas: ", arestas(h))
                println("    Γ_$k = $Γ (PI ótimo: $(info.otimo)), Δ = $(grau_maximo(h)), ",
                        "$(num_ciclos(h)) ciclo(s), coloração gulosa: $cor")
            end
            return
        end
    end
    println("Nenhum cacto com até $nmax vértices tem Γ_$k ≥ $t.")
end

args = parse.(Int, ARGS)
main(get(args, 1, 4), get(args, 2, 1), get(args, 3, 16))
