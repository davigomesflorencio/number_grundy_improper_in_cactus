# Decisão rápida de "Γ_k(G) ≥ t?" por busca de certificados.
#
# Γ_k é monótono em subgrafos induzidos: uma ordem ótima de G[S] seguida
# dos demais vértices colore S do mesmo jeito em G. Assim Γ_k(G) ≥ t se e
# só se existe um conjunto S e uma coloração de G[S] com cores em 1..t que
# é gulosa k-imprópria em G[S] (Lema da Caracterização) e dá a cor t a
# algum vértice r. Esse par (S, cor) é um *certificado*.
#
# A busca fixa cor[r] = t e só decide vértices que podem ajudar a bloquear
# um vértice de S ainda não bloqueado: seus vizinhos e, como um vizinho
# saturado também bloqueia, os vizinhos de seus vizinhos. Cada vértice é
# indefinido (-1), fora de S (0) ou colorido (1..t); indefinido conta como
# fora de S. A poda usa um limite otimista: cada vizinho indefinido
# bloqueia no máximo uma das cores que ainda faltam.

mutable struct _Busca
    g::Grafo
    k::Int
    t::Int
    cor::Vector{Int}       # -1 indefinido, 0 fora de S, c > 0 cor
    nind::Vector{Int}      # número de vizinhos indefinidos
    proibe_t::Vector{Bool} # vértices que já se sabe que não recebem a cor t
end

# vizinhos coloridos de v com a mesma cor de v
_mesma(b::_Busca, v::Int) = count(u -> b.cor[u] == b.cor[v], b.g.adj[v])

# v está bloqueado da cor j pelos vizinhos já coloridos?
function _bloqueado(b::_Busca, v::Int, j::Int)
    cnt = 0
    for x in b.g.adj[v]
        b.cor[x] == j || continue
        cnt += 1
        (cnt >= b.k + 1 || _mesma(b, x) == b.k) && return true
    end
    return false
end

# v ainda pode ficar bloqueado de j sem novos vizinhos coloridos?
# (conta um vizinho x de cor j que ainda pode saturar)
function _talvez_bloqueado(b::_Busca, v::Int, j::Int)
    cnt = 0
    for x in b.g.adj[v]
        b.cor[x] == j || continue
        cnt += 1
        cnt >= b.k + 1 && return true
        b.k == 0 && return true
        _mesma(b, x) + b.nind[x] >= b.k && return true
    end
    return false
end

_satisfeito(b::_Busca, v::Int) = all(j -> _bloqueado(b, v, j), 1:b.cor[v]-1)

function _viavel(b::_Busca, v::Int)
    _mesma(b, v) <= b.k || return false
    faltam = count(j -> !_talvez_bloqueado(b, v, j), 1:b.cor[v]-1)
    return faltam <= b.nind[v]
end

function _define!(b::_Busca, u::Int, c::Int)
    antes = b.cor[u]
    if antes == -1 && c != -1
        for y in b.g.adj[u]
            b.nind[y] -= 1
        end
    elseif antes != -1 && c == -1
        for y in b.g.adj[u]
            b.nind[y] += 1
        end
    end
    b.cor[u] = c
    return b
end

# a viabilidade só muda a distância ≤ 2 de u
function _consistente(b::_Busca, u::Int)
    b.cor[u] > 0 && !_viavel(b, u) && return false
    for w in b.g.adj[u]
        b.cor[w] > 0 && !_viavel(b, w) && return false
        for y in b.g.adj[w]
            y != u && b.cor[y] > 0 && !_viavel(b, y) && return false
        end
    end
    return true
end

# vértice indefinido que pode ajudar a bloquear v: um vizinho de v ou um
# vizinho de um vizinho x de v, de cor menor, que ainda não saturou
function _a_ramificar(b::_Busca, v::Int)
    for y in b.g.adj[v]
        b.cor[y] == -1 && return y
    end
    for x in b.g.adj[v]
        (0 < b.cor[x] < b.cor[v] && _mesma(b, x) < b.k) || continue
        for y in b.g.adj[x]
            b.cor[y] == -1 && return y
        end
    end
    return 0
end

# cores a tentar em u: primeiro as que faltam a v, por último "fora de S"
function _candidatas(b::_Busca, u::Int, v::Int)
    cmax = min(b.t, grau(b.g, u) + 1)
    b.proibe_t[u] && cmax == b.t && (cmax -= 1)
    faltam = [j for j in b.cor[v]-1:-1:1 if j <= cmax && !_bloqueado(b, v, j)]
    outras = [c for c in cmax:-1:1 if !(c in faltam)]
    return [faltam; outras; 0]
end

function _busca!(b::_Busca)
    # vértice de S ainda não bloqueado, com menos vizinhos indefinidos
    v, menor = 0, typemax(Int)
    for x in eachindex(b.cor)
        if b.cor[x] > 0 && b.nind[x] < menor && !_satisfeito(b, x)
            v, menor = x, b.nind[x]
        end
    end
    v == 0 && return true
    u = _a_ramificar(b, v)
    u == 0 && return false
    for c in _candidatas(b, u, v)
        _define!(b, u, c)
        _consistente(b, u) && _busca!(b) && return true
    end
    _define!(b, u, -1)
    return false
end

"""
    certificado_grundy(g, k, t) -> Union{Vector{Int}, Nothing}

Decide se Γ_k(g) ≥ t. Devolve um certificado `cor` (com `cor[v] = 0` para
os vértices fora de S): uma coloração gulosa k-imprópria de G[S] que usa a
cor `t`. Devolve `nothing` se não existir, o que prova Γ_k(g) < t.
`ordem_do_certificado(cor)` faz a heurística gulosa usar pelo menos `t`
cores em `g`.
"""
function certificado_grundy(g::Grafo, k::Integer, t::Integer)
    k >= 0 || throw(ArgumentError("k deve ser ≥ 0"))
    n = nv(g)
    t <= 0 && return zeros(Int, n)
    b = _Busca(g, Int(k), Int(t), fill(-1, n), [grau(g, v) for v in 1:n], falses(n))
    # a cor t exige grau ≥ t - 1; raízes de grau alto primeiro
    for r in sort(collect(1:n); by = v -> -grau(g, v))
        grau(g, r) >= t - 1 || break
        _define!(b, r, Int(t))
        _busca!(b) && return [max(c, 0) for c in b.cor]
        _define!(b, r, -1)
        b.proibe_t[r] = true        # nenhum certificado dá a cor t a r
    end
    return nothing
end

"""
    grundy_certificado(g, k) -> (Γ_k, cor)

Γ_k(g) exato pela busca de certificados: testa t = 2, 3, ... até que não
haja certificado. Devolve o certificado da maior cor (0 fora de S).
"""
function grundy_certificado(g::Grafo, k::Integer)
    nv(g) == 0 && return 0, Int[]
    melhor = certificado_grundy(g, k, 1)
    Γ = 1
    for t in 2:grau_maximo(g)+1
        cor = certificado_grundy(g, k, t)
        cor === nothing && break
        Γ, melhor = t, cor
    end
    return Γ, melhor
end

"""
    ordem_do_certificado(cor) -> Vector{Int}

Os vértices de S (`cor[v] > 0`) por cor crescente e, depois, os demais.
A heurística gulosa nessa ordem colore S exatamente como `cor`.
"""
function ordem_do_certificado(cor::AbstractVector{<:Integer})
    S = sort([v for v in eachindex(cor) if cor[v] > 0]; by = v -> cor[v])
    return vcat(S, [v for v in eachindex(cor) if cor[v] == 0])
end
