# Cactos: decomposição em blocos, reconhecimento, forma canônica e
# enumeração dos cactos não isomorfos.
#
# Um cacto é um grafo conexo em que toda aresta está em no máximo um ciclo,
# ou seja, todo bloco (componente biconexa) é uma aresta ou um ciclo.

"""
    blocos_arestas(g) -> Vector{Vector{Tuple{Int,Int}}}

Blocos (componentes biconexas) de `g`, cada um como a lista de suas arestas.
Algoritmo de Tarjan com pilha de arestas, em versão iterativa.
"""
function blocos_arestas(g::Grafo)
    n = nv(g)
    disc, low = zeros(Int, n), zeros(Int, n)
    tempo = 0
    pilha_a = Tuple{Int,Int}[]
    res = Vector{Vector{Tuple{Int,Int}}}()
    for s in 1:n
        disc[s] == 0 || continue
        tempo += 1
        disc[s] = low[s] = tempo
        pilha = [(s, 0, 1)]                # (vértice, pai, próximo vizinho)
        while !isempty(pilha)
            v, p, i = pilha[end]
            if i <= length(g.adj[v])
                pilha[end] = (v, p, i + 1)
                u = g.adj[v][i]
                if disc[u] == 0
                    push!(pilha_a, (v, u))
                    tempo += 1
                    disc[u] = low[u] = tempo
                    push!(pilha, (u, v, 1))
                elseif u != p && disc[u] < disc[v]
                    push!(pilha_a, (v, u))
                    low[v] = min(low[v], disc[u])
                end
            else
                pop!(pilha)
                p == 0 && continue
                low[p] = min(low[p], low[v])
                if low[v] >= disc[p]
                    B = Tuple{Int,Int}[]
                    while true
                        e = pop!(pilha_a)
                        push!(B, e)
                        e == (p, v) && break
                    end
                    push!(res, B)
                end
            end
        end
    end
    return res
end

# vértices de um bloco-ciclo na ordem cíclica (para uma ponte, os dois extremos)
function _ordem_ciclica(B::Vector{Tuple{Int,Int}})
    length(B) == 1 && return [B[1][1], B[1][2]]
    viz = Dict{Int,Vector{Int}}()
    for (u, v) in B
        push!(get!(viz, u, Int[]), v)
        push!(get!(viz, v, Int[]), u)
    end
    all(l -> length(l) == 2, values(viz)) || return nothing   # não é ciclo
    inicio = B[1][1]
    seq = [inicio]
    anterior, atual = inicio, viz[inicio][1]
    while atual != inicio
        push!(seq, atual)
        a, b = viz[atual]
        anterior, atual = atual, (a == anterior ? b : a)
    end
    return length(seq) == length(viz) ? seq : nothing
end

"""
    blocos(g) -> Vector{Vector{Int}}

Blocos de um cacto `g`: cada ponte vem como `[u, v]` e cada ciclo como a
lista de seus vértices na ordem cíclica. Lança `ArgumentError` se `g` tiver
um bloco que não é aresta nem ciclo.
"""
function blocos(g::Grafo)
    res = Vector{Int}[]
    for B in blocos_arestas(g)
        seq = _ordem_ciclica(B)
        seq === nothing && throw(ArgumentError("o grafo tem um bloco que não é aresta nem ciclo"))
        push!(res, seq)
    end
    return res
end

"`true` se `g` é conexo e cada bloco é uma aresta ou um ciclo."
function eh_cacto(g::Grafo)
    nv(g) >= 1 && eh_conexo(g) || return false
    return all(B -> length(B) == 1 || _ordem_ciclica(B) !== nothing, blocos_arestas(g))
end

"Número de ciclos de um cacto (pelo número ciclomático m - n + 1)."
num_ciclos(g::Grafo) = ne(g) - nv(g) + 1

# ----------------------------------------------------------------------
# Forma canônica. Enraizado num vértice r, o cacto é codificado como uma
# árvore de blocos (no estilo AHU): o código de um vértice é a lista
# ordenada dos códigos dos blocos-filhos; o de uma ponte é "e" seguido do
# código do outro extremo; o de um ciclo de comprimento L é "cL[...]" com
# os códigos dos outros L-1 vértices na ordem cíclica, no sentido que der a
# menor string. O mínimo sobre todas as raízes é um invariante completo.
# ----------------------------------------------------------------------

function _codigo(v::Int, pb::Int, B::Vector{Vector{Int}}, bl::Vector{Vector{Int}})
    filhos = String[]
    for b in bl[v]
        b == pb && continue
        cyc = B[b]
        L = length(cyc)
        i = findfirst(==(v), cyc)
        seq = [_codigo(cyc[mod1(i + d, L)], b, B, bl) for d in 1:L-1]
        if L == 2
            push!(filhos, "e" * seq[1])
        else
            push!(filhos, "c$L[" * min(join(seq, ","), join(reverse(seq), ",")) * "]")
        end
    end
    sort!(filhos)
    return "(" * join(filhos) * ")"
end

"""
    forma_canonica(g) -> String

String que identifica o cacto `g` a menos de isomorfismo: dois cactos são
isomorfos se e só se têm a mesma forma canônica. Custo O(n² log n).
"""
function forma_canonica(g::Grafo)
    B = blocos(g)
    bl = [Int[] for _ in 1:nv(g)]
    for (b, vs) in enumerate(B), v in vs
        push!(bl[v], b)
    end
    return minimum(_codigo(r, 0, B, bl) for r in 1:nv(g))
end

# cópia de g com um vértice novo pendurado em v
function _anexa_folha(g::Grafo, v::Int)
    h = Grafo([copy(a) for a in g.adj])
    push!(h.adj, Int[])
    return _liga!(h, v, nv(h))
end

# cópia de g com um ciclo novo de comprimento L passando por v
function _anexa_ciclo(g::Grafo, v::Int, L::Int)
    h = Grafo([copy(a) for a in g.adj])
    prox = nv(h) + 1
    append!(h.adj, [Int[] for _ in 1:L-1])
    return _ciclo_em!(h, v, prox, L)
end

"""
    cactos_nao_isomorfos(n) -> Vector{Grafo}

Todos os cactos com `n` vértices, um de cada classe de isomorfismo
(sequência OEIS A000083): n = 1..15 → 1, 1, 2, 4, 9, 23, 63, 188, 596,
1979, 6804, 24118, 87379, 322652, 1209808.

Todo cacto com n ≥ 2 vértices tem um bloco-folha (uma folha ou um ciclo com
um só vértice de corte); removê-lo deixa um cacto menor. Por isso os cactos
com n vértices são gerados anexando uma folha ou um ciclo a cada vértice
dos cactos menores, e as cópias isomorfas são descartadas pela forma
canônica.
"""
function cactos_nao_isomorfos(n::Integer)
    n >= 1 || return Grafo[]
    nivel = Vector{Vector{Grafo}}(undef, n)
    nivel[1] = [Grafo(1)]
    for m in 2:n
        res, vistos = Grafo[], Set{String}()
        function considera(h)
            c = forma_canonica(h)
            c in vistos && return
            push!(vistos, c)
            push!(res, h)
        end
        for g in nivel[m-1], v in 1:m-1
            considera(_anexa_folha(g, v))
        end
        for L in 3:m, g in nivel[m-L+1], v in 1:m-L+1
            considera(_anexa_ciclo(g, v, L))
        end
        nivel[m] = res
    end
    return nivel[n]
end
