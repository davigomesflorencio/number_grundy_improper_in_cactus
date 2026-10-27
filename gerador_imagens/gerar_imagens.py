"""Gera imagens PNG dos cactos listados nos resultados de menor_cacto.jl.

Para cada arquivo de resultados, cria uma pasta com uma imagem por cacto
(vértices pintados com a coloração gulosa encontrada pelo PI) e uma visão
geral com todos os cactos lado a lado.

Uso (da raiz do repositório):
    python gerador_imagens/gerar_imagens.py
    python gerador_imagens/gerar_imagens.py GrundyImproprioCacto/resultados/menor_cacto_gamma1_3.txt
    python gerador_imagens/gerar_imagens.py --layout planar --saida outra_pasta
"""

import argparse
import ast
import math
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

RAIZ = Path(__file__).resolve().parent.parent
RESULTADOS = RAIZ / "GrundyImproprioCacto" / "resultados"
SAIDA = Path(__file__).resolve().parent / "imagens"

# paleta Okabe-Ito (segura para daltonismo), uma cor por cor da coloração
PALETA = ["#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#F0E442", "#0072B2", "#D55E00"]

RE_CABECALHO = re.compile(r"Menor cacto com Γ_(\d+) ≥ (\d+)")
RE_DADOS = re.compile(
    r"Γ_(\d+) = (\d+).*Δ = (\d+), (\d+) ciclo\(s\), coloração gulosa: (\[.*\])"
)


def le_resultados(arquivo):
    """Devolve (k, t, cactos); cada cacto é um dict com arestas e dados."""
    linhas = Path(arquivo).read_text(encoding="utf-8").splitlines()
    k = t = None
    cactos = []
    for i, linha in enumerate(linhas):
        m = RE_CABECALHO.search(linha)
        if m:
            k, t = int(m.group(1)), int(m.group(2))
        if linha.strip().startswith("arestas:"):
            arestas = ast.literal_eval(linha.split(":", 1)[1].strip())
            d = RE_DADOS.search(linhas[i + 1])
            if d is None:
                raise ValueError(f"{arquivo}:{i + 2}: linha de dados inesperada")
            cactos.append({
                "arestas": arestas,
                "k": int(d.group(1)),
                "grundy": int(d.group(2)),
                "delta": int(d.group(3)),
                "ciclos": int(d.group(4)),
                "cor": ast.literal_eval(d.group(5)),
            })
    return k, t, cactos


def grafo(cacto):
    g = nx.Graph()
    g.add_nodes_from(range(1, len(cacto["cor"]) + 1))
    g.add_edges_from(cacto["arestas"])
    return g


def posicoes(g, layout):
    if layout == "planar":
        return nx.planar_layout(g)
    if layout == "spring":
        return nx.spring_layout(g, seed=1)
    return nx.kamada_kawai_layout(g)


def desenha(ax, cacto, layout, rotulos=True, tam_no=420):
    """Desenha um cacto: arestas de ciclo sólidas, pontes tracejadas, e o
    vértice de maior cor com borda grossa."""
    g = grafo(cacto)
    pos = posicoes(g, layout)
    em_ciclo = {frozenset(e) for c in nx.cycle_basis(g) for e in zip(c, c[1:] + c[:1])}
    ciclo = [e for e in g.edges if frozenset(e) in em_ciclo]
    ponte = [e for e in g.edges if frozenset(e) not in em_ciclo]
    nx.draw_networkx_edges(g, pos, edgelist=ciclo, ax=ax, width=2.2, edge_color="#333333")
    nx.draw_networkx_edges(g, pos, edgelist=ponte, ax=ax, width=1.6,
                           edge_color="#777777", style="dashed")
    cor = cacto["cor"]
    nos = list(g.nodes)
    topo = max(cor)
    nx.draw_networkx_nodes(
        g, pos, nodelist=nos, ax=ax, node_size=tam_no,
        node_color=[PALETA[(cor[v - 1] - 1) % len(PALETA)] for v in nos],
        edgecolors=["#000000" if cor[v - 1] == topo else "#444444" for v in nos],
        linewidths=[3.0 if cor[v - 1] == topo else 1.0 for v in nos],
    )
    if rotulos:
        nx.draw_networkx_labels(g, pos, ax=ax, font_size=9,
                                labels={v: f"{v}\n{cor[v - 1]}" for v in nos})
    ax.set_axis_off()
    ax.margins(0.12)


def legenda(fig, t):
    itens = [Patch(facecolor=PALETA[c - 1], edgecolor="#444444", label=f"cor {c}")
             for c in range(1, t + 1)]
    itens += [Line2D([0], [0], color="#333333", lw=2.2, label="aresta de ciclo"),
              Line2D([0], [0], color="#777777", lw=1.6, ls="--", label="ponte")]
    fig.legend(handles=itens, loc="lower center", ncol=len(itens), frameon=False, fontsize=9)


def gera(arquivo, saida, layout, dpi):
    k, t, cactos = le_resultados(arquivo)
    if not cactos:
        print(f"{arquivo}: nenhum cacto encontrado")
        return
    pasta = saida / Path(arquivo).stem
    pasta.mkdir(parents=True, exist_ok=True)
    n = len(cactos[0]["cor"])
    total = len(cactos)
    t = t or max(c["grundy"] for c in cactos)

    for i, cacto in enumerate(cactos, start=1):
        fig, ax = plt.subplots(figsize=(6, 6))
        desenha(ax, cacto, layout)
        ax.set_title(f"Cacto {i}/{total}: n = {n}, Γ_{cacto['k']} = {cacto['grundy']}, "
                     f"Δ = {cacto['delta']}, {cacto['ciclos']} ciclo(s)\n"
                     "rótulo: vértice / cor gulosa", fontsize=10)
        legenda(fig, t)
        fig.savefig(pasta / f"cacto_{i:02d}.png", dpi=dpi, bbox_inches="tight")
        plt.close(fig)

    colunas = math.ceil(math.sqrt(total))
    linhas = math.ceil(total / colunas)
    fig, eixos = plt.subplots(linhas, colunas, figsize=(2.6 * colunas, 2.6 * linhas + 0.8))
    eixos = [eixos] if total == 1 else list(getattr(eixos, "flat", eixos))
    for i, ax in enumerate(eixos):
        if i < total:
            desenha(ax, cactos[i], layout, rotulos=False, tam_no=90)
            ax.set_title(f"{i + 1}", fontsize=9)
        else:
            ax.set_axis_off()
    fig.suptitle(f"Os {total} cactos não isomorfos com {n} vértices e Γ_{k} = {t}", fontsize=14)
    legenda(fig, t)
    fig.savefig(pasta / "visao_geral.png", dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"{arquivo}: {total} imagem(ns) + visao_geral.png em {pasta}")


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("arquivos", nargs="*", type=Path,
                   help="arquivos de resultados (padrão: todos os .txt em GrundyImproprioCacto/resultados)")
    p.add_argument("--saida", type=Path, default=SAIDA, help="pasta das imagens")
    p.add_argument("--layout", choices=["kamada_kawai", "planar", "spring"], default="kamada_kawai")
    p.add_argument("--dpi", type=int, default=120)
    args = p.parse_args()
    arquivos = args.arquivos or sorted(RESULTADOS.glob("*.txt"))
    for arquivo in arquivos:
        gera(arquivo, args.saida, args.layout, args.dpi)


if __name__ == "__main__":
    main()
