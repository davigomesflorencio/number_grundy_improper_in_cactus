# Geradores de grafos e cactos.

using Random

caminho(n::Integer) = grafo_de_arestas(n, [(i, i + 1) for i in 1:n-1])
ciclo(n::Integer) = (n >= 3 || throw(ArgumentError("ciclo precisa de n ≥ 3"));
                     grafo_de_arestas(n, [[(i, i + 1) for i in 1:n-1]; (n, 1)]))
completo(n::Integer) = grafo_de_arestas(n, [(i, j) for i in 1:n for j in i+1:n])

"Estrela com `n` vértices; o centro é o vértice 1."
estrela(n::Integer) = grafo_de_arestas(n, [(1, i) for i in 2:n])

# liga v, w1, ..., w_{L-1} num ciclo, com wi = prox, ..., prox+L-2
function _ciclo_em!(g::Grafo, v::Int, prox::Int, L::Int)
    anterior = v
    for w in prox:prox+L-2
        _liga!(g, anterior, w)
        anterior = w
    end
    _liga!(g, anterior, v)
    return g
end

"""
    moinho(m, L) -> Grafo

`m` ciclos de comprimento `L` com um único vértice em comum (o vértice 1).
Tem 1 + m(L-1) vértices.
"""
function moinho(m::Integer, L::Integer)
    L >= 3 || throw(ArgumentError("os ciclos precisam de comprimento ≥ 3"))
    g = Grafo(1 + m * (L - 1))
    for i in 0:m-1
        _ciclo_em!(g, 1, 2 + i * (L - 1), Int(L))
    end
    return g
end

"Grafo da amizade F_m: `m` triângulos com um vértice em comum (o vértice 1)."
amizade(m::Integer) = moinho(m, 3)

"""
    cadeia_de_ciclos(m, L) -> Grafo

`m` ciclos de comprimento `L` em fila: cada ciclo compartilha um vértice com
o seguinte, e os dois vértices de corte de um ciclo ficam a distância L÷2.
Tem 1 + m(L-1) vértices.
"""
function cadeia_de_ciclos(m::Integer, L::Integer)
    L >= 3 || throw(ArgumentError("os ciclos precisam de comprimento ≥ 3"))
    g = Grafo(1 + m * (L - 1))
    corte = 1
    for i in 0:m-1
        prox = 2 + i * (L - 1)
        _ciclo_em!(g, corte, prox, Int(L))
        corte = prox + L ÷ 2 - 1
    end
    return g
end

"""
    sol(L) -> Grafo

Ciclo C_L (vértices 1..L) com uma folha pendurada em cada vértice.
"""
function sol(L::Integer)
    g = ciclo(L)
    append!(g.adj, [Int[] for _ in 1:L])
    for i in 1:L
        _liga!(g, i, L + i)
    end
    return g
end

"""
    cacto_aleatorio(n; p_ciclo = 0.5, tam_max = 6, rng) -> Grafo

Cacto conexo com `n` vértices, construído bloco a bloco: a cada passo um
vértice já existente é sorteado e recebe, com probabilidade `p_ciclo`, um
ciclo novo de comprimento uniforme em 3..tam_max (limitado pelos vértices
que faltam) e, caso contrário, uma folha. Todo cacto pode ser gerado assim,
mas a distribuição não é uniforme. Com `p_ciclo = 0` gera árvores
recursivas; com `p_ciclo = 1` e `tam_max = 3`, cactos de triângulos (para n
ímpar).
"""
function cacto_aleatorio(n::Integer; p_ciclo::Real = 0.5, tam_max::Integer = 6,
                         rng::AbstractRNG = Random.default_rng())
    n >= 1 || throw(ArgumentError("n deve ser ≥ 1"))
    g = Grafo(n)
    m = 1
    while m < n
        v = rand(rng, 1:m)
        Lmax = min(tam_max, n - m + 1)
        if Lmax >= 3 && rand(rng) < p_ciclo
            L = rand(rng, 3:Lmax)
            _ciclo_em!(g, v, m + 1, L)
            m += L - 1
        else
            _liga!(g, v, m + 1)
            m += 1
        end
    end
    return g
end

"Grafo aleatório G(n, p) de Erdős–Rényi."
function grafo_aleatorio(n::Integer, p::Real; rng::AbstractRNG = Random.default_rng())
    g = Grafo(n)
    for i in 1:n, j in i+1:n
        rand(rng) < p && _liga!(g, i, j)
    end
    return g
end

"Cópia de `g` com os vértices renomeados: o vértice v passa a ser `perm[v]`."
function renomeia(g::Grafo, perm::AbstractVector{<:Integer})
    h = Grafo(nv(g))
    for (u, v) in arestas(g)
        _liga!(h, Int(perm[u]), Int(perm[v]))
    end
    return h
end
