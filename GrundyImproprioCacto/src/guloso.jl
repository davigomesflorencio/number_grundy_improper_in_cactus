# Heurística gulosa k-imprópria, verificações e força bruta (grafos gerais).

"""
    coloracao_gulosa(g, ordem, k) -> Vector{Int}

Executa a heurística gulosa k-imprópria na `ordem` dada (uma permutação
dos vértices). O vértice v recebe a menor cor c tal que
  (i)  v tem no máximo k vizinhos já coloridos com c, e
  (ii) nenhum vizinho de v com cor c já tem k vizinhos com cor c.
"""
function coloracao_gulosa(g::Grafo, ordem::AbstractVector{<:Integer}, k::Integer)
    n = nv(g)
    length(ordem) == n || throw(ArgumentError("a ordem deve conter os $n vértices"))
    cor = zeros(Int, n)
    # cnt[v][c] = número de vizinhos já coloridos de v com a cor c. Como v
    # só fica bloqueado de uma cor se tiver um vizinho com ela, cor[v] ≤
    # grau(v) + 1, e cores maiores nunca são consultadas: memória O(n + m)
    cnt = [zeros(Int, length(g.adj[v]) + 1) for v in 1:n]
    for v in ordem
        cor[v] == 0 || throw(ArgumentError("vértice $v repetido na ordem"))
        c = 1
        while true
            if cnt[v][c] <= k
                ok = true
                for u in g.adj[v]
                    if cor[u] == c && cnt[u][c] >= k
                        ok = false
                        break
                    end
                end
                ok && break
            end
            c += 1
        end
        cor[v] = c
        for u in g.adj[v]
            c <= length(cnt[u]) && (cnt[u][c] += 1)
        end
    end
    return cor
end

"`true` se cada vértice tem no máximo k vizinhos da própria cor."
function eh_k_impropria(g::Grafo, cor::AbstractVector{<:Integer}, k::Integer)
    all(v -> count(u -> cor[u] == cor[v], g.adj[v]) <= k, 1:nv(g))
end

"""
    eh_gulosa(g, cor, k) -> Bool

Verifica, pelo Lema da Caracterização, se `cor` é uma coloração gulosa
k-imprópria: ela é k-imprópria e todo vértice v está bloqueado de toda
cor j < cor[v], isto é, tem ≥ k+1 vizinhos de cor j ou um vizinho de
cor j saturado (com exatamente k vizinhos da própria cor).
"""
function eh_gulosa(g::Grafo, cor::AbstractVector{<:Integer}, k::Integer)
    eh_k_impropria(g, cor, k) || return false
    mesma = [count(u -> cor[u] == cor[v], g.adj[v]) for v in 1:nv(g)]
    for v in 1:nv(g), j in 1:cor[v]-1
        viz_j = [u for u in g.adj[v] if cor[u] == j]
        (length(viz_j) >= k + 1 || any(u -> mesma[u] == k, viz_j)) || return false
    end
    return true
end

# próxima permutação em ordem lexicográfica (in place); false ao terminar
function _proxima_permutacao!(p::Vector{Int})
    i = length(p) - 1
    while i >= 1 && p[i] >= p[i+1]
        i -= 1
    end
    i == 0 && return false
    j = length(p)
    while p[j] <= p[i]
        j -= 1
    end
    p[i], p[j] = p[j], p[i]
    reverse!(p, i + 1, length(p))
    return true
end

"""
    grundy_forca_bruta(g, k) -> (Γ_k, ordem)

Γ_k(G) exato testando todas as n! ordens (use só com n ≤ 10).
Para assim que atinge o limite Δ(G) + 1.
"""
function grundy_forca_bruta(g::Grafo, k::Integer)
    n = nv(g)
    n == 0 && return 0, Int[]
    limite = grau_maximo(g) + 1
    p = collect(1:n)
    melhor, melhor_ordem = 0, copy(p)
    while true
        m = maximum(coloracao_gulosa(g, p, k))
        if m > melhor
            melhor, melhor_ordem = m, copy(p)
            melhor == limite && break
        end
        _proxima_permutacao!(p) || break
    end
    return melhor, melhor_ordem
end

"""
    ordem_da_coloracao(cor) -> Vector{Int}

Vértices por cor crescente. Se `cor` é gulosa k-imprópria, a heurística
gulosa nessa ordem reproduz exatamente `cor`.
"""
ordem_da_coloracao(cor::AbstractVector{<:Integer}) = sortperm(cor)
