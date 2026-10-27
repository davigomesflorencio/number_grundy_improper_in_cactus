# GrundyImproprioCacto.jl

Implementação em Julia do **número de Grundy k-impróprio** Γ_k(G) em
**cactos**, com foco em k = 1: o maior número de cores que a heurística
gulosa 1-imprópria pode usar num cacto G, considerando todas as ordens
possíveis dos vértices.

Segue a estrutura de
[number_grundy_improper_in_three](https://github.com/davigomesflorencio/number_grundy_improper_in_three)
(árvores). O módulo traz:

- estruturas de grafo (listas de adjacência) e leitura/escrita de lista de arestas;
- geradores de cactos (ciclos, moinhos, grafos da amizade, cadeias de ciclos,
  sóis, cactos aleatórios) e a enumeração de todos os cactos não isomorfos;
- decomposição em blocos (Tarjan), reconhecimento de cactos e forma canônica;
- a heurística gulosa k-imprópria, a verificação pelo Lema da Caracterização
  e a força bruta (usada só para validar);
- duas **formulações de programação inteira** para Γ_k, resolvidas com
  JuMP + HiGHS:
  - **atribuição**: a formulação da seção 5.2.2 da dissertação de
    E. N. H. D. Rodrigues, *Coloração k-imprópria gulosa* (UFC, 2020);
  - **representantes**: adaptação, para a coloração gulosa k-imprópria, da
    formulação por representantes da coloração (Campêlo, Campos e Corrêa).

## Sumário

- [Definições](#definições)
- [Formulações](#formulações)
- [Requisitos](#requisitos)
- [Como executar](#como-executar)
- [Uso no REPL](#uso-no-repl)
- [Referência da API](#referência-da-api)
- [Formato de arquivo](#formato-de-arquivo)
- [Estrutura do projeto](#estrutura-do-projeto)
- [Resultados em subclasses de cactos](#resultados-em-subclasses-de-cactos)
- [Desempenho](#desempenho)
- [Changelog](#changelog)

## Definições

Uma coloração é **k-imprópria** quando cada vértice tem no máximo `k`
vizinhos da mesma cor. Um vértice com exatamente `k` vizinhos da sua cor
está **saturado**.

A **heurística gulosa k-imprópria** percorre os vértices numa ordem dada e
atribui a cada vértice `v` a menor cor `c` tal que:

1. `v` tem no máximo `k` vizinhos já coloridos com `c`; e
2. nenhum vizinho de `v` com cor `c` já tem `k` vizinhos com cor `c`.

**Γ_k(G)** é o maior número de cores que essa heurística usa, considerando
todas as ordens. Vale Γ_k(G) ≤ Δ(G) + 1.

**Lema da Caracterização.** Uma coloração é gulosa k-imprópria se e só se é
k-imprópria e todo vértice `v` está bloqueado de toda cor `j < cor(v)`, isto
é, tem pelo menos `k + 1` vizinhos de cor `j` ou um vizinho de cor `j`
saturado. Ordenar os vértices por cor reproduz a coloração.

Um **cacto** é um grafo conexo em que toda aresta está em no máximo um
ciclo, ou seja, todo bloco é uma aresta ou um ciclo.

## Formulações

### Atribuição (Rodrigues, 2020, seção 5.2.2)

Cores `C = {1, …, Δ+1}`. Variáveis binárias `x[v,c]` (v tem a cor c),
`y[v,c]` (v tem mais que k vizinhos com a cor c), `w[v,c]` (v tem a cor c e
está saturado), `p[v,c]` (a cor c é proibida para v) e `z[c]` (a cor c é
usada). Maximiza `Σ z[c]` sujeito a (5.9), (5.10) e (5.18)–(5.27): cada
vértice recebe uma cor, `y` e `w` ficam amarrados às contagens de vizinhos,
`p[v,c] = 1` quando `y[v,c] = 1` ou quando um vizinho saturado em `c` não
depende de `v`, e `x[v,c] ≤ p[v,d]` para todo `d < c` (v recebe a menor cor
permitida).

A dissertação limita as cores por n/(k+1). Aqui usamos Δ + 1, que também é
um limitante válido e, nos cactos, costuma ser bem menor.

### Representantes (adaptação)

Cada classe de cor é representada pelo seu vértice de menor índice.
Variáveis binárias:

- `x[u,v]`, `u ≤ v`: v está na classe representada por u (`x[u,u] = 1` se u é representante);
- `o[u,v]`: u é representante e a classe de u vem antes da classe de v;
- `y[v,u]`: v tem mais que k vizinhos na classe de u;
- `w[a,u]`: a está na classe de u e tem k vizinhos nela (saturado).

Maximiza o número de representantes `Σ x[u,u]` sujeito a:

- cada vértice tem exatamente um representante, e `x[u,v] ≤ x[u,u]`;
- cada classe é k-imprópria: `Σ_{a∈N(v)} x[u,a] ≤ k + (d(v) − k)(1 − x[u,v])`;
- `o` define uma ordem total e transitiva dos representantes, compatível com
  as classes (`o[u,v] = o[u,r]` quando `x[r,v] = 1`);
- bloqueio: `o[u,v] ≤ y[v,u] + Σ_{a∈N(v)} w[a,u]`, com
  `(k+1)·y[v,u] ≤ Σ_{a∈N(v)} x[u,a]`, `w[a,u] ≤ x[u,a]` e
  `k·w[a,u] ≤ Σ_{b∈N(a)} x[u,b]`.

Pelo Lema da Caracterização, as soluções viáveis são exatamente as
colorações gulosas k-impróprias, com cor(v) = 1 + o número de classes antes
da classe de v. Ao contrário da formulação de atribuição, essa não tem a
simetria entre cores iguais com rótulos trocados, mas tem O(n³) restrições
de ordem.

Os testes conferem as duas formulações com a força bruta em todos os
cactos com até 8 vértices e em grafos gerais aleatórios.

## Requisitos

- Julia 1.x (testado na 1.13.2).
- [JuMP](https://jump.dev) e [HiGHS](https://highs.dev), instalados no
  ambiente do projeto (`Project.toml` e `Manifest.toml` em
  `GrundyImproprioCacto/`). Na primeira vez:

```bash
cd GrundyImproprioCacto
julia --project=. -e "using Pkg; Pkg.instantiate()"
```

Para usar outro solver (por exemplo, Gurobi), passe
`otimizador = Gurobi.Optimizer` a qualquer função das formulações.

## Como executar

Rode os scripts de dentro da pasta `GrundyImproprioCacto/`:

```bash
julia --project=. testes.jl     # testes de correção (cerca de 80 s)
julia --project=. exemplos.jl   # exemplos, tabelas e comparação das formulações (cerca de 55 s)
julia --project=. menor_cacto.jl 4 1 16   # menor cacto com Γ_1 ≥ 4 (cerca de 50 min)
```

`menor_cacto.jl [t] [k] [nmax]` procura o menor cacto com Γ_k ≥ t; com
`3 1 10` termina em poucos segundos.

## Uso no REPL

```julia
julia> cd("GrundyImproprioCacto"); using Pkg; Pkg.activate(".")

julia> include("src/GrundyImproprioCacto.jl"); using .GrundyImproprioCacto

julia> G = cacto_aleatorio(15);

julia> Γ, cor, info = grundy_improprio_cacto(G)        # Γ_1, formulação de representantes

julia> grundy_improprio_cacto(G; formulacao = :atribuicao)

julia> coloracao_gulosa(G, ordem_da_coloracao(cor), 1) == cor   # confere por simulação
true
```

Mais exemplos:

```julia
# k qualquer, em qualquer grafo, com limite de tempo
Γ, cor, info = grundy_pi_atribuicao(moinho(3, 4), 2; tempo_limite = 60)
info.otimo       # false se o limite foi atingido; então Γ é um limitante inferior

# todos os cactos não isomorfos com 8 vértices
cs = cactos_nao_isomorfos(8)          # 188 cactos
count(g -> grundy_pi_atribuicao(g, 1)[1] == 3, cs)

# o modelo JuMP, para inspecionar ou alterar
m, x = modelo_atribuicao(ciclo(6), 1)
```

## Referência da API

### Grafos (`src/grafos.jl`)

| Função | Descrição |
|---|---|
| `Grafo(n)` | grafo sem arestas com vértices `1..n` |
| `nv(g)`, `ne(g)` | número de vértices e de arestas |
| `vizinhos(g, v)`, `grau(g, v)`, `grau_maximo(g)` | vizinhança e graus |
| `tem_aresta(g, u, v)`, `adiciona_aresta!(g, u, v)` | consulta e inserção (laços não são permitidos) |
| `arestas(g)`, `grafo_de_arestas(n, E)` | lista de arestas e construção a partir dela |
| `bfs(g, raiz)`, `eh_conexo(g)`, `eh_arvore(g)` | busca em largura e testes estruturais |
| `le_lista_arestas(arq)`, `escreve_lista_arestas(arq, g)` | E/S de arquivo ([formato](#formato-de-arquivo)) |

### Geradores (`src/geradores.jl`)

| Função | Descrição |
|---|---|
| `caminho(n)`, `ciclo(n)`, `completo(n)`, `estrela(n)` | famílias clássicas |
| `moinho(m, L)` | `m` ciclos de comprimento `L` com um vértice comum |
| `amizade(m)` | grafo da amizade F_m = `moinho(m, 3)` |
| `cadeia_de_ciclos(m, L)` | `m` ciclos de comprimento `L` em fila, consecutivos com um vértice comum |
| `sol(L)` | ciclo C_L com uma folha em cada vértice |
| `cacto_aleatorio(n; p_ciclo, tam_max, rng)` | cacto construído anexando folhas e ciclos a vértices sorteados |
| `grafo_aleatorio(n, p; rng)` | G(n, p) de Erdős–Rényi |
| `renomeia(g, perm)` | cópia de `g` com os vértices renomeados |

### Cactos (`src/cacto.jl`)

| Função | Descrição |
|---|---|
| `blocos_arestas(g)` | blocos (componentes biconexas) como listas de arestas |
| `blocos(g)` | blocos de um cacto: pontes `[u, v]` e ciclos na ordem cíclica |
| `eh_cacto(g)`, `num_ciclos(g)` | reconhecimento e número de ciclos |
| `forma_canonica(g)` | string que identifica o cacto a menos de isomorfismo |
| `cactos_nao_isomorfos(n)` | um cacto de cada classe de isomorfismo (OEIS A000083) |

### Heurística gulosa (`src/guloso.jl`)

| Função | Descrição |
|---|---|
| `coloracao_gulosa(g, ordem, k)` | executa a heurística na ordem dada |
| `eh_k_impropria(g, cor, k)`, `eh_gulosa(g, cor, k)` | verificações (Lema da Caracterização) |
| `ordem_da_coloracao(cor)` | vértices por cor crescente, que reproduz uma coloração gulosa |
| `grundy_forca_bruta(g, k)` | Γ_k testando as n! ordens; só para validar, com n ≤ 9 |

### Busca de certificados (`src/certificado.jl`)

| Função | Descrição |
|---|---|
| `certificado_grundy(g, k, t)` | decide Γ_k ≥ t: devolve uma coloração gulosa de um subgrafo induzido que usa a cor `t` (0 fora dele), ou `nothing` |
| `grundy_certificado(g, k)` | Γ_k exato testando t = 2, 3, … e o certificado da maior cor |
| `ordem_do_certificado(cor)` | ordem que faz a heurística gulosa reproduzir o certificado |

É bem mais rápida que o PI para responder "Γ_k ≥ t?" em cactos pequenos,
mas o tempo é exponencial no pior caso.

### Programação inteira (`src/formulacoes.jl`)

| Função | Descrição |
|---|---|
| `modelo_atribuicao(g, k; cmax, otimizador, tempo_limite, silencioso)` | modelo JuMP da formulação de atribuição e as variáveis `x` |
| `grundy_pi_atribuicao(g, k; ...)` | resolve e devolve `(Γ_k, cor, info)` |
| `modelo_representantes(g, k; otimizador, tempo_limite, silencioso)` | modelo JuMP da formulação por representantes e as variáveis `x`, `o` |
| `grundy_pi_representantes(g, k; ...)` | resolve e devolve `(Γ_k, cor, info)` |
| `grundy_improprio_cacto(G, k = 1; formulacao = :representantes, ...)` | Γ_k de um cacto; lança `ArgumentError` se `G` não for cacto |

`info` tem os campos `status`, `otimo`, `tempo`, `limitante`, `variaveis` e
`restricoes`. Γ_k só é exato quando `info.otimo` é `true`.

## Formato de arquivo

O mesmo do projeto de árvores: número de vértices na primeira linha e uma
aresta por linha. Linhas vazias e linhas começando com `#` são ignoradas.

```text
# ciclo com 4 vértices
4
1 2
2 3
3 4
4 1
```

## Estrutura do projeto

```text
.
├── README.md
├── CHANGELOG.md
└── GrundyImproprioCacto/
    ├── Project.toml, Manifest.toml  # dependências (JuMP, HiGHS)
    ├── exemplos.jl                  # exemplos, tabelas e comparação das formulações
    ├── testes.jl                    # testes de correção
    ├── menor_cacto.jl               # menor cacto com Γ_k ≥ t
    ├── resultados/                  # saída de menor_cacto.jl para t = 3 e t = 4
    └── src/
        ├── GrundyImproprioCacto.jl  # módulo e exports
        ├── grafos.jl                # estrutura Grafo, BFS, E/S
        ├── geradores.jl             # geradores de cactos e grafos
        ├── cacto.jl                 # blocos, reconhecimento, forma canônica, enumeração
        ├── guloso.jl                # heurística gulosa e força bruta
        ├── certificado.jl           # decisão rápida de Γ_k ≥ t
        └── formulacoes.jl           # formulações de atribuição e de representantes
```

## Resultados em subclasses de cactos

Valores de Γ_1 calculados pela formulação de atribuição, todos com
otimalidade provada pelo HiGHS. São resultados computacionais para os
tamanhos indicados, não teoremas para todo n.

| Subclasse | Gerador | Tamanhos testados | Γ_1 |
|---|---|---|---|
| Caminho P_n | `caminho(n)` | n = 1..10 | 1 para n ≤ 2; 2 para n ≥ 3 |
| Estrela K_{1,n-1} | `estrela(n)` | n = 2..10 | 1 para n = 2; 2 para n ≥ 3 |
| Ciclo C_L | `ciclo(L)` | L = 3..12 | 2 |
| Amizade F_m (m triângulos com um vértice comum) | `amizade(m)` | m = 1..6 | 2 |
| Sol (C_L com uma folha em cada vértice) | `sol(L)` | L = 3..8 | 2 |
| Moinho de ciclos C_L, L ≥ 4 | `moinho(m, L)` | m = 1..4, L = 4..6 | 2 para m = 1 (é o ciclo); 3 para m ≥ 2 |
| Cadeia de triângulos | `cadeia_de_ciclos(m, 3)` | m = 1..5 | 2 para m ≤ 3; 3 para m ≥ 4 |
| Cadeia de ciclos C_L, L ≥ 4 | `cadeia_de_ciclos(m, L)` | m = 1..5, L = 4..6 | 2 para m = 1; 3 para m ≥ 2 |

Observações:

- **O menor cacto com Γ_1 = 3 tem 7 vértices**: um triângulo 1-2-3 com uma
  folha em 2 e um vértice 4 ligado a 1, com duas folhas. Em árvores, Γ_1 = t
  exige pelo menos 3^(t-1) vértices, ou seja, 9 para t = 3 (o mínimo é
  atingido pela árvore binomial B^1_3).
- Com n = 7, 15 dos 63 cactos têm Γ_1 = 3; com n = 8, 83 dos 188 (seção 3
  de `exemplos.jl`).
- Ciclos, grafos da amizade e sóis ficam em Γ_1 = 2, mesmo com grau máximo
  grande, como no grafo da amizade. Já dois ciclos de comprimento ≥ 4 com um
  vértice comum bastam para Γ_1 = 3. Com triângulos, são necessários quatro
  em cadeia.

### Menor cacto com Γ_1 = t

Γ_1 é monótono em subgrafos induzidos, e todo subgrafo induzido conexo de um
cacto é um cacto. Então o menor cacto com Γ_1 = t é o primeiro que aparece
quando se percorrem todos os cactos não isomorfos por número de vértices.
É o que faz `menor_cacto.jl`, decidindo Γ_1 ≥ t em cada cacto pela busca de
certificados ([`src/certificado.jl`](GrundyImproprioCacto/src/certificado.jl)).
Essa busca foi conferida contra a força bruta e contra o PI em cerca de 4 mil
casos, incluindo todos os cactos com 9 e 10 vértices, sem nenhuma diferença.

| t | menor n (cactos) | cactos não isomorfos com esse n e Γ_1 = t | menor n (árvores) |
|---|---|---|---|
| 1 | 1 | 1 | 1 |
| 2 | 3 | 2 (P_3 e C_3) | 3 |
| 3 | 7 | 15 | 9 |
| 4 | **15** | **72** | 27 |

Para t = 3 (`menor_cacto.jl 3 1 10`, menos de 1 s), nenhum cacto com até 6
vértices tem Γ_1 ≥ 3. Com 7 vértices há 15, todos confirmados pelo PI. Nenhum
deles é árvore, o que é coerente com o mínimo de 9 vértices nas árvores.
Oito têm Δ = 3 e um ciclo; sete têm Δ = 4 e dois ciclos. Os ciclos que
aparecem:

| Ciclos | Cactos |
|---|---|
| um C_3 | 2 |
| um C_4 | 4 |
| um C_5 | 2 |
| dois C_3 | 2 |
| C_3 e C_4 | 3 |
| C_3 e C_5 | 1 |
| dois C_4 (o moinho com dois quadrados) | 1 |

A lista completa está em
[`resultados/menor_cacto_gamma1_3.txt`](GrundyImproprioCacto/resultados/menor_cacto_gamma1_3.txt).

Para t = 4, o script percorreu todos os 1 653 627 cactos com até 15 vértices
(cerca de 50 minutos). Nenhum cacto com até 14 vértices tem Γ_1 ≥ 4, e os
72 de 15 vértices foram confirmados pelo PI (Γ_1 = 4, otimalidade provada).
Todos têm Δ = 4; 48 têm 4 ciclos e 24 têm 5. A lista completa está em
[`resultados/menor_cacto_gamma1_4.txt`](GrundyImproprioCacto/resultados/menor_cacto_gamma1_4.txt).
Um deles:

```text
arestas: (1,2) (1,3) (1,4) (1,5) (2,3) (2,6) (2,7) (3,11) (3,12) (4,5)
         (4,8) (5,13) (6,7) (6,14) (8,9) (8,10) (9,10) (9,15)
cores:   3 2 4 3 1 1 2 2 1 2 1 1 1 1 1   (vértices 1..15)
```

Os mínimos 1, 3, 7, 15 seguem 2^t − 1. Em árvores, eles são 3^(t−1).
Se o padrão continuar, os cactos satisfazem Γ_1 ≤ log₂(n + 1), um limite
bem mais justo que 1 + 2·log₂ n. Isso é só uma **conjectura**, apoiada em
quatro valores. O próximo caso (t = 5, n = 31) está fora do alcance da
enumeração completa.

## Desempenho

Tempos do HiGHS para Γ_1 em cactos aleatórios (seção 5 de `exemplos.jl`):

| n | Δ | Γ_1 | atribuição (s) | representantes (s) |
|---|---|---|---|---|
| 10 | 4 | 3 | 0.06 | 0.24 |
| 15 | 3 | 3 | 0.48 | 1.05 |
| 20 | 6 | 3 | 1.82 | 3.69 |
| 25 | 6 | 3 | 5.78 | 4.37 |

Nos cactos com até 20 vértices, a formulação de atribuição costuma ser mais
rápida. A de representantes cresce como O(n³) em restrições, e fica mais
competitiva quando há muitas cores possíveis (Δ grande). Os dois tempos
crescem rápido com n: para cactos maiores, use `tempo_limite` e confira
`info.otimo`.

## Changelog

O histórico de mudanças está em [CHANGELOG.md](CHANGELOG.md).
