# Changelog

Todas as mudanças relevantes deste projeto são registradas aqui.

O formato segue o [Keep a Changelog](https://keepachangelog.com/pt-BR/1.1.0/)
e o projeto adota o [Versionamento Semântico](https://semver.org/lang/pt-BR/).

## [Não lançado]

### Adicionado

- Estrutura `Grafo`, heurística gulosa k-imprópria, verificação pelo Lema da
  Caracterização e força bruta, vindas do projeto de árvores
  (number_grundy_improper_in_three).
- Geradores de cactos: `moinho`, `amizade`, `cadeia_de_ciclos`, `sol` e
  `cacto_aleatorio`, além de `renomeia`.
- Cactos: decomposição em blocos (Tarjan iterativo), `eh_cacto`,
  `forma_canonica` e `cactos_nao_isomorfos` (confere com a OEIS A000083 até
  n = 11).
- Formulação de programação inteira de atribuição da dissertação de
  Rodrigues (2020, seção 5.2.2) e uma formulação por representantes adaptada
  à coloração gulosa k-imprópria, com JuMP + HiGHS.
- `grundy_improprio_cacto`, `ordem_da_coloracao`, `testes.jl` e `exemplos.jl`.

[Não lançado]: #nao-lancado
