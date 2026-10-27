# Estrutura de grafo simples, não direcionado, com vértices 1..n
# armazenado como listas de adjacência.

struct Grafo
    adj::Vector{Vector{Int}}
end

"Grafo sem arestas com `n` vértices."
Grafo(n::Integer) = Grafo([Int[] for _ in 1:n])

nv(g::Grafo) = length(g.adj)
ne(g::Grafo) = sum(length, g.adj; init = 0) ÷ 2
vizinhos(g::Grafo, v::Integer) = g.adj[v]
grau(g::Grafo, v::Integer) = length(g.adj[v])
grau_maximo(g::Grafo) = nv(g) == 0 ? 0 : maximum(length, g.adj)
tem_aresta(g::Grafo, u::Integer, v::Integer) = v in g.adj[u]

"""
    adiciona_aresta!(g, u, v) -> Bool

Adiciona a aresta {u, v}. Devolve `false` se ela já existia.
Laços não são permitidos.
"""
function adiciona_aresta!(g::Grafo, u::Integer, v::Integer)
    u == v && throw(ArgumentError("laço em $u não é permitido"))
    (1 <= u <= nv(g) && 1 <= v <= nv(g)) || throw(BoundsError(g.adj, (u, v)))
    tem_aresta(g, u, v) && return false
    push!(g.adj[u], v)
    push!(g.adj[v], u)
    return true
end

# versão sem verificação de duplicidade, usada pelos geradores de árvores
function _liga!(g::Grafo, u::Int, v::Int)
    push!(g.adj[u], v)
    push!(g.adj[v], u)
    return g
end

"Lista de arestas (u, v) com u < v."
arestas(g::Grafo) = [(u, v) for u in 1:nv(g) for v in g.adj[u] if u < v]

"Constrói um grafo com `n` vértices a partir de uma coleção de pares (u, v)."
function grafo_de_arestas(n::Integer, E)
    g = Grafo(n)
    for (u, v) in E
        adiciona_aresta!(g, u, v)
    end
    return g
end

"""
    bfs(g, raiz) -> (pai, ordem)

Busca em largura a partir de `raiz`. `pai[raiz] == 0` e `pai[v] == -1`
para vértices não alcançados. `ordem` lista os vértices alcançados.
"""
function bfs(g::Grafo, raiz::Integer)
    pai = fill(-1, nv(g))
    pai[raiz] = 0
    ordem = [Int(raiz)]
    i = 1
    while i <= length(ordem)
        v = ordem[i]
        for u in g.adj[v]
            if pai[u] == -1
                pai[u] = v
                push!(ordem, u)
            end
        end
        i += 1
    end
    return pai, ordem
end

eh_conexo(g::Grafo) = nv(g) == 0 || length(bfs(g, 1)[2]) == nv(g)
eh_arvore(g::Grafo) = nv(g) >= 1 && ne(g) == nv(g) - 1 && eh_conexo(g)

"""
    le_lista_arestas(caminho) -> Grafo

Lê um arquivo no formato

    n
    u v
    u v
    ...

com vértices numerados de 1 a n. Linhas vazias e linhas iniciadas
por `#` são ignoradas.
"""
function le_lista_arestas(caminho::AbstractString)
    linhas = [strip(l) for l in eachline(caminho)]
    linhas = [l for l in linhas if !isempty(l) && !startswith(l, "#")]
    isempty(linhas) && throw(ArgumentError("arquivo vazio: $caminho"))
    n = parse(Int, linhas[1])
    g = Grafo(n)
    for l in linhas[2:end]
        u, v = parse.(Int, split(l))
        adiciona_aresta!(g, u, v)
    end
    return g
end

"Grava `g` no mesmo formato lido por `le_lista_arestas`."
function escreve_lista_arestas(caminho::AbstractString, g::Grafo)
    open(caminho, "w") do io
        println(io, nv(g))
        for (u, v) in arestas(g)
            println(io, u, " ", v)
        end
    end
    return caminho
end

Base.show(io::IO, g::Grafo) = print(io, "Grafo(", nv(g), " vértices, ", ne(g), " arestas)")
