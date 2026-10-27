# Formulações de programação inteira para Γ_k (grafos gerais).
#
# 1. Atribuição: a formulação (5.8)-(5.10), (5.18)-(5.29) da dissertação de
#    E. N. H. D. Rodrigues, "Coloração k-imprópria gulosa" (UFC, 2020),
#    seção 5.2.2. x[v,c] = 1 se v recebe a cor c.
#
# 2. Representantes: adaptação, para a coloração gulosa k-imprópria, da
#    formulação por representantes da coloração (Campêlo, Campos e Corrêa).
#    Cada classe de cor é representada pelo seu vértice de menor índice e as
#    classes são ordenadas por variáveis de precedência. Pelo Lema da
#    Caracterização, basta exigir que cada classe seja k-imprópria e que todo
#    vértice esteja bloqueado de todas as classes anteriores à sua.

using JuMP
import HiGHS

# cria o modelo com o otimizador, o limite de tempo e o modo silencioso
function _modelo(otimizador, tempo_limite, silencioso)
    m = Model(otimizador)
    silencioso && set_silent(m)
    tempo_limite === nothing || set_time_limit_sec(m, tempo_limite)
    return m
end

# resultado comum às duas formulações
function _resultado(m::Model, cor::Vector{Int})
    st = termination_status(m)
    info = (status = st, otimo = st == MOI.OPTIMAL, tempo = solve_time(m),
            limitante = has_values(m) ? objective_bound(m) : NaN,
            variaveis = num_variables(m),
            restricoes = num_constraints(m; count_variable_in_set_constraints = false))
    return maximum(cor; init = 0), cor, info
end

"""
    modelo_atribuicao(g, k; cmax, otimizador, tempo_limite, silencioso) -> (m, x)

Monta a formulação de atribuição da dissertação de Rodrigues (2020, seção
5.2.2). As cores vão de 1 a `cmax` (padrão Δ + 1, que limita Γ_k).
Variáveis binárias:

- `x[v,c]`: v recebe a cor c;
- `y[v,c]`: v tem mais que k vizinhos com a cor c;
- `w[v,c]`: v tem a cor c e está saturado (k vizinhos com a cor c);
- `p[v,c]`: a cor c é proibida para v;
- `z[c]`: a cor c é usada.
"""
function modelo_atribuicao(g::Grafo, k::Integer; cmax::Integer = grau_maximo(g) + 1,
                           otimizador = HiGHS.Optimizer, tempo_limite = nothing,
                           silencioso::Bool = true)
    n, V, C = nv(g), 1:nv(g), 1:cmax
    N(v) = g.adj[v]
    d(v) = length(g.adj[v])
    m = _modelo(otimizador, tempo_limite, silencioso)
    @variable(m, x[V, C], Bin)
    @variable(m, y[V, C], Bin)
    @variable(m, w[V, C], Bin)
    @variable(m, p[V, C], Bin)
    @variable(m, z[C], Bin)
    @objective(m, Max, sum(z))
    @constraint(m, [v in V], sum(x[v, c] for c in C) == 1)                          # (5.9)
    @constraint(m, [c in C], z[c] <= sum(x[v, c] for v in V))                       # (5.10)
    for v in V, c in C
        sx = sum(x[u, c] for u in N(v); init = AffExpr(0))
        @constraint(m, (d(v) - k) * y[v, c] >= sx - k)                               # (5.18)
        @constraint(m, (k + 1) * y[v, c] <= sx)                                      # (5.19)
        @constraint(m, w[v, c] <= x[v, c])                                           # (5.20)
        @constraint(m, k * w[v, c] <= sx)                                            # (5.21)
        @constraint(m, w[v, c] + (1 - x[v, c]) * (d(v) - k + 1) >= sx - k + 1)       # (5.22)
        @constraint(m, p[v, c] <= y[v, c] + sum(w[u, c] for u in N(v); init = AffExpr(0)))  # (5.23)
        @constraint(m, p[v, c] >= y[v, c])                                           # (5.24)
        for u in N(v)
            @constraint(m, p[v, c] >= w[u, c] - x[v, c])                             # (5.25)
        end
        @constraint(m, x[v, c] <= 1 - p[v, c])                                       # (5.26)
        for dd in 1:c-1
            @constraint(m, x[v, c] <= p[v, dd])                                      # (5.27)
        end
    end
    return m, x
end

"""
    grundy_pi_atribuicao(g, k; kwargs...) -> (Γ_k, cor, info)

Resolve a formulação de atribuição e devolve o número de cores, a coloração
gulosa encontrada e `info = (status, otimo, tempo, limitante, variaveis,
restricoes)`. Γ_k só é exato quando `info.otimo` é `true`; com limite de
tempo, `Γ_k` é um limitante inferior e `info.limitante` um superior.
Aceita os argumentos de `modelo_atribuicao`.
"""
function grundy_pi_atribuicao(g::Grafo, k::Integer; kwargs...)
    k >= 0 || throw(ArgumentError("k deve ser ≥ 0"))
    nv(g) == 0 && return 0, Int[], nothing
    m, x = modelo_atribuicao(g, k; kwargs...)
    optimize!(m)
    has_values(m) || return _resultado(m, zeros(Int, nv(g)))
    C = axes(x, 2)
    cor = [findfirst(c -> value(x[v, c]) > 0.5, C) for v in 1:nv(g)]
    return _resultado(m, cor)
end

"""
    modelo_representantes(g, k; otimizador, tempo_limite, silencioso) -> (m, x, o)

Formulação por representantes (adaptada). Variáveis binárias:

- `x[u,v]`, u ≤ v: v está na classe representada por u (o menor vértice da
  classe); `x[u,u] = 1` se u é representante;
- `o[u,v]`: u é representante e a classe de u vem antes da classe de v;
- `y[v,u]`: v tem mais que k vizinhos na classe de u;
- `w[a,u]`: a está na classe de u e tem k vizinhos nela (saturado).

Restrições: cada vértice tem um representante; as classes são k-impróprias;
as precedências formam uma ordem total dos representantes, compatível com
as classes; e se a classe de u vem antes da de v, então v está bloqueado da
classe de u (y[v,u] = 1 ou um vizinho de v em u saturado). O objetivo é
maximizar o número de representantes.
"""
function modelo_representantes(g::Grafo, k::Integer; otimizador = HiGHS.Optimizer,
                               tempo_limite = nothing, silencioso::Bool = true)
    n, V = nv(g), 1:nv(g)
    N(v) = g.adj[v]
    d(v) = length(g.adj[v])
    m = _modelo(otimizador, tempo_limite, silencioso)
    P = [(u, v) for u in V for v in u:n]
    @variable(m, x[P], Bin)
    @variable(m, o[V, V], Bin)
    @variable(m, y[V, V], Bin)
    @variable(m, w[V, V], Bin)
    X(u, v) = u <= v ? x[(u, v)] : AffExpr(0)      # x[u,v] só existe para u ≤ v
    # quantos vizinhos de v estão na classe de u
    viz_em(v, u) = sum(X(u, a) for a in N(v); init = AffExpr(0))
    @objective(m, Max, sum(x[(u, u)] for u in V))
    for v in V
        @constraint(m, sum(x[(u, v)] for u in 1:v) == 1)              # um representante
        @constraint(m, o[v, v] == 0)
    end
    for (u, v) in P
        u < v && @constraint(m, x[(u, v)] <= x[(u, u)])               # u representa a si mesmo
    end
    for u in V, v in V
        # classe k-imprópria: se v está na classe de u, no máximo k vizinhos nela
        if u <= v && d(v) > k
            @constraint(m, viz_em(v, u) <= k + (d(v) - k) * (1 - x[(u, v)]))
        end
        # precedência só entre representantes, e não com a própria classe
        @constraint(m, o[u, v] <= x[(u, u)])
        u <= v && @constraint(m, o[u, v] <= 1 - x[(u, v)])
        # bloqueio: classe de u antes da de v ⇒ v bloqueado da classe de u
        @constraint(m, (k + 1) * y[v, u] <= viz_em(v, u))
        @constraint(m, w[v, u] <= X(u, v))
        @constraint(m, k * w[v, u] <= viz_em(v, u))
        @constraint(m, o[u, v] <= y[v, u] + sum(w[a, u] for a in N(v); init = AffExpr(0)))
        u == v && continue
        # ordem total entre representantes
        u < v && @constraint(m, o[u, v] + o[v, u] >= x[(u, u)] + x[(v, v)] - 1)
        u < v && @constraint(m, o[u, v] + o[v, u] <= 1)
        # v na classe de r ⇒ o[u,v] = o[u,r]
        for r in 1:v
            r == u && continue
            @constraint(m, o[u, v] >= o[u, r] + x[(r, v)] - 1)
            @constraint(m, o[u, v] <= o[u, r] + 1 - x[(r, v)])
        end
    end
    # transitividade: u antes de r e r antes de v ⇒ u antes de v
    for u in V, r in V, v in V
        length(unique((u, r, v))) == 3 || continue
        @constraint(m, o[u, r] + o[r, v] <= 1 + o[u, v])
    end
    return m, x, o
end

"""
    grundy_pi_representantes(g, k; kwargs...) -> (Γ_k, cor, info)

Resolve a formulação por representantes. A cor de cada vértice é 1 mais o
número de classes que vêm antes da sua. Mesma saída de
`grundy_pi_atribuicao` e mesmos argumentos de `modelo_representantes`.
"""
function grundy_pi_representantes(g::Grafo, k::Integer; kwargs...)
    k >= 0 || throw(ArgumentError("k deve ser ≥ 0"))
    nv(g) == 0 && return 0, Int[], nothing
    m, x, o = modelo_representantes(g, k; kwargs...)
    optimize!(m)
    n = nv(g)
    has_values(m) || return _resultado(m, zeros(Int, n))
    rep = [findfirst(u -> value(x[(u, v)]) > 0.5, 1:v) for v in 1:n]
    cor = [1 + count(u -> value(o[u, rep[v]]) > 0.5, 1:n) for v in 1:n]
    return _resultado(m, cor)
end

"""
    grundy_improprio_cacto(G, k = 1; formulacao = :representantes, kwargs...) -> (Γ_k, cor, info)

Γ_k de um cacto (por padrão, o número de Grundy 1-impróprio) por uma das
formulações (`:atribuicao` ou `:representantes`). Lança `ArgumentError` se
`G` não for cacto.
"""
function grundy_improprio_cacto(G::Grafo, k::Integer = 1; formulacao::Symbol = :representantes,
                                kwargs...)
    eh_cacto(G) || throw(ArgumentError("o grafo não é um cacto"))
    formulacao == :atribuicao && return grundy_pi_atribuicao(G, k; kwargs...)
    formulacao == :representantes && return grundy_pi_representantes(G, k; kwargs...)
    throw(ArgumentError("formulação desconhecida: $formulacao"))
end
