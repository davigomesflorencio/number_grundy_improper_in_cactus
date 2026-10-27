"""
Número de Grundy k-impróprio em cactos: estruturas de grafo, geradores e
enumeração de cactos, heurística gulosa k-imprópria, força bruta e as
formulações de programação inteira de atribuição e de representantes.

Uso (de dentro de GrundyImproprioCacto/, com `julia --project=.`):
    include("src/GrundyImproprioCacto.jl")
    using .GrundyImproprioCacto
"""
module GrundyImproprioCacto

using Random

include("grafos.jl")
include("geradores.jl")
include("cacto.jl")
include("guloso.jl")
include("formulacoes.jl")

# estruturas de grafo
export Grafo, nv, ne, vizinhos, grau, grau_maximo, tem_aresta, adiciona_aresta!,
       arestas, grafo_de_arestas, bfs, eh_conexo, eh_arvore,
       le_lista_arestas, escreve_lista_arestas
# geradores
export caminho, ciclo, completo, estrela, moinho, amizade, cadeia_de_ciclos, sol,
       cacto_aleatorio, grafo_aleatorio, renomeia
# cactos
export blocos_arestas, blocos, eh_cacto, num_ciclos, forma_canonica,
       cactos_nao_isomorfos
# heurística gulosa
export coloracao_gulosa, eh_k_impropria, eh_gulosa, grundy_forca_bruta
# programação inteira
export modelo_atribuicao, grundy_pi_atribuicao, modelo_representantes,
       grundy_pi_representantes, grundy_improprio_cacto, ordem_da_coloracao

end
