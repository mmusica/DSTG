import networkx as nx
from ipy_table import *
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pylab
import itertools as it
import random
import copy
import os

plt.rcParams['figure.facecolor'] = 'white'

from IPython.core.display import HTML
from IPython.core.display import Image

PLATON = {'DODEKAEDAR': {0: [0, 2.466], 1: [-2.345, 0.762], 2: [-1.449, -1.995], 3: [1.449, -1.995], 19: [2.345, 0.762],
                         10: [0, 1.628], 8: [-1.548, 0.503], 6: [-0.957, -1.317], 4: [0.957, -1.317],
                         18: [1.548, 0.503],
                         9: [-0.7, 0.965], 11: [0.7, 0.965], 17: [1.134, -0.368], 5: [0, -1.192], 7: [-1.134, -0.368],
                         13: [-0.302, 0.416], 12: [0.302, 0.416], 16: [0.489, -0.159], 15: [0, -0.514],
                         14: [-0.489, -0.159]},
          'IKOZAEDAR': {1: [0, -0.314], 10: [0, 3.602], 8: [0, -1.165], 4: [0, 0.78], 5: [0.272, 0.157],
                        6: [-0.272, 0.157],
                        7: [3.12, -1.801], 9: [-3.12, -1.801], 0: [0.675, -0.39], 2: [-0.675, -0.39],
                        11: [1.009, 0.583], 3: [-1.009, 0.583]},
          'OKTAEDAR': {0: [-0.316, 0.182], 1: [0.316, 0.182], 2: [0, -0.365], 3: [0, 1.385], 4: [-1.2, -0.693],
                       5: [1.2, -0.693]},
          'KOCKA': {0: [-0.333, -0.333], 4: [-1, -1], 1: [-0.333, 0.333], 7: [-1, 1], 3: [0.333, -0.333], 5: [1, -1],
                    2: [0.333, 0.333], 6: [1, 1]},
          'TETRAEDAR': {0: [-0.866, -0.5], 1: [0.866, -0.5], 2: [0, 1], 3: [0, 0]}}


def ispis(objekt, br):
    objekt = str(objekt)
    d = len(objekt)
    red = d // br
    ostatak = d % br
    za_ispis = ""
    for k in range(red):
        za_ispis = za_ispis + objekt[k * br:(k + 1) * br] + "\n"
    if ostatak == 0:
        za_ispis = za_ispis[:-1]
    else:
        za_ispis = za_ispis + objekt[red * br:]
    print(za_ispis)


def kvocijent(a, b):
    if b == 0:
        return "Error: dijeljenje s nulom"
    elif b > 0:
        return int(np.floor(a / b))
    else:
        return int(np.ceil(a / b))


def ostatak(a, b):
    if b == 0:
        return "Error: dijeljenje s nulom"
    elif b > 0:
        return int(a - np.floor(a / b) * b)
    else:
        return int(a - np.ceil(a / b) * b)


def euklid_rucno(a, b, izlaz='html'):
    qi = []
    xi = [1, 0]
    yi = [0, 1]
    ri = [a, b]
    while True:
        q = kvocijent(ri[-2], ri[-1])
        r = ri[-2] - q * ri[-1]
        if r != 0:
            qi.append(q)
            xi.append(xi[-2] - q * xi[-1])
            yi.append(yi[-2] - q * yi[-1])
            ri.append(ri[-2] - q * ri[-1])
        else:
            break
    rjecnik = {}
    rjecnik['ri'] = ri
    rjecnik['qi'] = qi
    rjecnik['xi'] = xi
    rjecnik['yi'] = yi
    rjecnik['NZM'] = ri[-1]
    if izlaz == 'rjecnik':
        return rjecnik
    else:
        indeksi = list(range(-1, len(qi) + 1))
        red1 = [' '] + indeksi
        red2 = ['r<sub>i</sub>'] + ri
        red3 = ['q<sub>i</sub>', ' ', ' '] + qi
        red4 = ['x<sub>i</sub>'] + xi
        red5 = ['y<sub>i</sub>'] + yi
        tablica = [red1, red2, red3, red4, red5]
        tablica = make_table(tablica)
        apply_theme('basic_both')
        set_global_style(align='center')
        return tablica


def lin_kong_rucno(a, b, n, izlaz='html'):
    qi = []
    yi = [0, 1]
    n1, q, r1, d = a, n // a, n % a, a
    while r1 != 0:
        qi.append(q)
        yi.append(yi[-2] - q * yi[-1])
        if r1 != 0: d = r1
        q = n1 // r1
        n1, r1 = r1, n1 % r1
    if b % d != 0: return False
    a1, b1, n1 = a // d, b // d, n // d
    x1 = (yi[-1] * b1) % n1
    rjesenja = [x1 + k * n1 for k in range(d)]
    rjecnik = {}
    rjecnik['qi'] = qi
    rjecnik['yi'] = yi
    rjecnik['a1b1n1'] = [a1, b1, n1]
    rjecnik['rjesenja'] = rjesenja
    if izlaz == 'rjecnik':
        return rjecnik
    else:
        indeksi = list(range(-1, len(qi) + 1))
        red1 = [' '] + indeksi + ['rjesenja']
        red2 = ['q<sub>i</sub>', ' ', ' '] + qi + [rjesenja]
        red3 = ['y<sub>i</sub>'] + yi + ['(a<sub>1</sub>,b<sub>1</sub>,n<sub>1</sub>)=({},{},{})'.format(a1, b1, n1)]
        tablica = [red1, red2, red3]
        tablica = make_table(tablica)
        apply_theme('basic_both')
        set_global_style(align='center')
        return tablica


def KTO_rucno(A, N, izlaz='html'):
    for t in it.combinations(N, 2):
        if np.gcd(t[0], t[1]) != 1: return "Error: moduli nisu u parovima relativno prosti."
    n = np.prod(N)
    ki = list(map(lambda x: n // x, N))
    xi = []
    for i in range(len(N)):
        xi.append(lin_kong_rucno(ki[i], A[i], N[i], izlaz='rjecnik')['rjesenja'][0])
    x0 = sum(list(map(lambda x: x[0] * x[1], zip(ki, xi))))
    x0mod = x0 % n
    if izlaz == 'rjecnik':
        return dict(zip(['n', 'k_i', 'x_i', 'x0', 'x0_mod'], [n, ki, xi, x0, x0mod]))
    else:
        indeksi = [' '] + list(range(1, len(A) + 1))
        ki = ['k<sub>i</sub>'] + ki
        xi = ['x<sub>i</sub>'] + xi
        x0 = ['x<sub>0</sub>', x0] + [' '] * (len(A) - 1)
        x0mod = ['x<sub>0</sub>-mod', x0mod] + [' '] * (len(A) - 1)
        n = ['n', n] + [' '] * (len(A) - 1)
        tablica = [indeksi, ki, xi, x0, x0mod, n]
        tablica = make_table(tablica)
        apply_theme('basic_both')
        set_global_style(align='center')
        set_cell_style(3, 1, column_span=len(A))
        set_cell_style(4, 1, column_span=len(A))
        set_cell_style(5, 1, column_span=len(A))
        return tablica


def struk(G):
    if G.number_of_selfloops() > 0: return 1
    if G.is_multigraph(): return 2
    strukG = G.number_of_nodes()
    parent = {}
    D = {}
    for v in G.nodes():
        S = []
        R = [v]
        parent[v] = 'null'
        D[v] = 0
        while len(R) != 0:
            x = R[0]
            S.append(x)
            del R[0]
            for y in G.neighbors(x):
                if y == parent[x]:
                    continue
                if not (y in S):
                    parent[y] = x
                    D[y] = D[x] + 1
                    R.append(y)
                else:
                    strukG = min(strukG, D[x] + D[y] + 1)
    return strukG


def flatten(a):
    for elem in a:
        if type(elem) in (tuple, list):
            for i in flatten(elem):
                yield i
        else:
            yield elem


def pretvori(G):
    if len(set(G.edges())) == len(G.edges()): return G
    mul = {}
    bridovi = list(G.edges())
    for (x, y) in bridovi:
        mul[(x, y)] = bridovi.count((x, y))
    H = nx.Graph(G)
    for x in H.edges():
        H[x[0]][x[1]]['weight'] = mul[x]
    return H


def minimum_cut_phase(G1, W, a):
    G = G1.copy()
    V = G.nodes()
    A = [a]
    del W[a]
    for k in W:
        if G.get_edge_data(a, k):
            W[k] += G.get_edge_data(a, k)['weight']

    while len(A) != len(V):
        m = max(W, key=W.get)
        A.append(m)
        del W[m]
        for k in W:
            if G.get_edge_data(m, k):
                W[k] += G.get_edge_data(m, k)['weight']

    rez1 = list(flatten([A[-1]]))
    rez = [rez1, list(set(list(flatten(V))).difference(set(rez1)))]
    tezina_reza = sum(list(map(lambda x: G[A[-1]][x]['weight'], list(G.neighbors(A[-1])))))

    st = A[-2:]
    if st[0] in G.neighbors(st[1]):
        G.remove_edge(st[0], st[1])
    susjedi = set(G.neighbors(st[0])).union(set(G.neighbors(st[1])))
    presjek = set(G.neighbors(st[0])).intersection(set(G.neighbors(st[1])))
    susjedi0 = susjedi.difference(set(G.neighbors(st[1])))
    susjedi1 = susjedi.difference(set(G.neighbors(st[0])))
    for x in susjedi:
        if x in presjek:
            tezina = G[x][st[0]]['weight'] + G[x][st[1]]['weight']
            G.add_weighted_edges_from([(x, tuple(st), tezina)])
        if x in susjedi0:
            tezina = G[x][st[0]]['weight']
            G.add_weighted_edges_from([(x, tuple(st), tezina)])
        if x in susjedi1:
            tezina = G[x][st[1]]['weight']
            G.add_weighted_edges_from([(x, tuple(st), tezina)])
    G.remove_node(st[0])
    G.remove_node(st[1])

    return (tezina_reza, rez, G)


def minimum_edge_cut(G1):
    G = G1.copy()
    G = pretvori(G)
    V = list(G.nodes())
    E = list(G.edges())
    nema_tezine = G[E[0][0]][E[0][1]]
    if not (nema_tezine):
        for x in E:
            G[x[0]][x[1]]['weight'] = 1
    minimalna_tezina = 0
    for b in E:
        minimalna_tezina += G[b[0]][b[1]]['weight']
    a = random.choice(V)
    W = {}
    for x in V:
        W[x] = 0
    broj = len(V)

    while broj > 1:
        rezultat = minimum_cut_phase(G, W, a)
        if rezultat[0] < minimalna_tezina:
            minimalni_rez = rezultat[1]
            minimalna_tezina = rezultat[0]
        G = rezultat[2]
        W = {}
        for x in G.nodes():
            W[x] = 0
        broj = len(G.nodes())

    bridovi_minimalnog_reza = []
    for u in minimalni_rez[0]:
        for v in minimalni_rez[1]:
            if (u, v) in G1.edges():
                bridovi_minimalnog_reza.append((u, v))

    return (minimalna_tezina, minimalni_rez, bridovi_minimalnog_reza)


def CrtajTezinskiGraf(G, pos, omjeri=None, pomaci=None, slika=None, velicinaVrha=300, bojaVrha='y', rubVrha=None,
                      oblikVrha='o', debljinaBrida=1, bojaBrida='k', fontTezine=12, fontVrh=12):
    pozicijaOznake = {}
    if omjeri == None:
        for u in G.edges().__iter__():
            pozicijaOznake[u] = 0.5
    else:
        pozicijaOznake = omjeri
        for u in G.edges().__iter__():
            if not (u in pozicijaOznake.keys()):
                pozicijaOznake[u] = 0.5
    dxdy = {}
    if pomaci == None:
        for u in G.edges().__iter__():
            dxdy[u] = (0, 0)
    else:
        dxdy = pomaci
        for u in G.edges().__iter__():
            if not (u in dxdy.keys()):
                dxdy[u] = (0, 0)
    rjecnik = {}
    for u in G.edges(data=True).__iter__():
        rjecnik[(u[0], u[1])] = str(u[2]['weight'])
    if rubVrha == None:
        nx.draw_networkx_nodes(G, pos, node_size=velicinaVrha, node_color=bojaVrha, node_shape=oblikVrha)
    else:
        nx.draw_networkx_nodes(G, pos, node_size=velicinaVrha, node_color=bojaVrha, node_shape=oblikVrha,
                               edgecolors=rubVrha)
    nx.draw_networkx_edges(G, pos, width=debljinaBrida, edge_color=bojaBrida)
    # nx.draw_networkx_edge_labels(G,pos,edge_labels=rjecnik)
    for (u, v) in G.edges().__iter__():
        plt.text((1 - pozicijaOznake[(u, v)]) * pos[u][0] + pozicijaOznake[(u, v)] * pos[v][0] + dxdy[(u, v)][0],
                 (1 - pozicijaOznake[(u, v)]) * pos[u][1] + pozicijaOznake[(u, v)] * pos[v][1] + dxdy[(u, v)][1],
                 rjecnik[(u, v)], bbox=dict(facecolor='white', edgecolor='white', alpha=1), size=fontTezine)
    nx.draw_networkx_labels(G, pos, font_size=fontVrh)
    # pylab.gcf().gca().set_axis_off()
    plt.axis('off')
    pylab.gcf().set_facecolor('w')
    if slika == None:
        plt.axis('auto')
    else:
        plt.xlim(*slika[0])
        plt.ylim(*slika[1])
    plt.show()
    return


def CrtajNajkraciPut(G, pos, prvi, drugi, tezine=True, omjeri=None, pomaci=None, slika=None, velicinaVrha=300,
                     bojaVrha='y', oblikVrha='o', rubVrha=None, debljinaBrida=1, bojaBrida='k', fontTezine=12,
                     fontVrh=12, alfa=0.5, bojaPuta='r', debljinaPuta=5, bojaVrha2='cyan'):
    # crtanje tezinskog grafa
    vrhovi = list(G.nodes())
    vrhovi.remove(prvi)
    vrhovi.remove(drugi)
    if tezine:
        pozicijaOznake = {}
        if omjeri == None:
            for u in G.edges().__iter__():
                pozicijaOznake[u] = 0.5
        else:
            pozicijaOznake = omjeri
            for u in G.edges().__iter__():
                if not (u in pozicijaOznake.keys()):
                    pozicijaOznake[u] = 0.5
        dxdy = {}
        if pomaci == None:
            for u in G.edges().__iter__():
                dxdy[u] = (0, 0)
        else:
            dxdy = pomaci
            for u in G.edges().__iter__():
                if not (u in dxdy.keys()):
                    dxdy[u] = (0, 0)
        rjecnik = {}
        for u in G.edges(data=True).__iter__():
            rjecnik[(u[0], u[1])] = str(u[2]['weight'])
        for (u, v) in G.edges().__iter__():
            plt.text((1 - pozicijaOznake[(u, v)]) * pos[u][0] + pozicijaOznake[(u, v)] * pos[v][0] + dxdy[(u, v)][0],
                     (1 - pozicijaOznake[(u, v)]) * pos[u][1] + pozicijaOznake[(u, v)] * pos[v][1] + dxdy[(u, v)][1],
                     rjecnik[(u, v)], bbox=dict(facecolor='white', edgecolor='white', alpha=1), size=fontTezine)
    if rubVrha == None:
        nx.draw_networkx_nodes(G, pos, nodelist=vrhovi, node_size=velicinaVrha, node_color=bojaVrha,
                               node_shape=oblikVrha)
    else:
        nx.draw_networkx_nodes(G, pos, nodelist=vrhovi, node_size=velicinaVrha, node_color=bojaVrha,
                               node_shape=oblikVrha, edgecolors=rubVrha)
    nx.draw_networkx_edges(G, pos, width=debljinaBrida, edge_color=bojaBrida)
    nx.draw_networkx_labels(G, pos, font_size=fontVrh)
    # isticanje najkraceg puta
    istaknutiBridovi = []
    if tezine:
        for i in range(len(nx.dijkstra_path(G, prvi, drugi)) - 1):
            istaknutiBridovi.append((nx.dijkstra_path(G, prvi, drugi)[i], nx.dijkstra_path(G, prvi, drugi)[i + 1]))
    else:
        for i in range(len(nx.shortest_path(G, prvi, drugi)) - 1):
            istaknutiBridovi.append((nx.shortest_path(G, prvi, drugi)[i], nx.shortest_path(G, prvi, drugi)[i + 1]))
    nx.draw_networkx_edges(G, pos, edgelist=istaknutiBridovi, width=debljinaPuta, edge_color=bojaPuta, alpha=alfa)
    if rubVrha == None:
        nx.draw_networkx_nodes(G, pos, nodelist=[prvi, drugi], node_size=velicinaVrha, node_color=bojaVrha2,
                               node_shape=oblikVrha)
    else:
        nx.draw_networkx_nodes(G, pos, nodelist=[prvi, drugi], node_size=velicinaVrha, node_color=bojaVrha2,
                               node_shape=oblikVrha, edgecolors=rubVrha)
    plt.axis('off')
    pylab.gcf().set_facecolor('w')
    if slika == None:
        plt.axis('auto')
    else:
        plt.xlim(*slika[0])
        plt.ylim(*slika[1])
    plt.show()
    return


def StabloNajkracihPutova(G, pos, pocetak, tezine=True, omjeri=None, pomaci=None, slika=None, velicinaVrha=300,
                          bojaVrha='y', oblikVrha='o', rubVrha=None, debljinaBrida=1, bojaBrida='k', fontTezine=12,
                          fontVrh=12, alfa=0.5, bojaPuta='r', debljinaPuta=5, bojaVrha2='cyan'):
    # crtanje tezinskog grafa
    vrhovi = list(G.nodes())
    vrhovi.remove(pocetak)
    if tezine:
        pozicijaOznake = {}
        if omjeri == None:
            for u in G.edges().__iter__():
                pozicijaOznake[u] = 0.5
        else:
            pozicijaOznake = omjeri
            for u in G.edges().__iter__():
                if not (u in pozicijaOznake.keys()):
                    pozicijaOznake[u] = 0.5
        dxdy = {}
        if pomaci == None:
            for u in G.edges().__iter__():
                dxdy[u] = (0, 0)
        else:
            dxdy = pomaci
            for u in G.edges().__iter__():
                if not (u in dxdy.keys()):
                    dxdy[u] = (0, 0)
        rjecnik = {}
        for u in G.edges(data=True).__iter__():
            rjecnik[(u[0], u[1])] = str(u[2]['weight'])
        for (u, v) in G.edges().__iter__():
            plt.text((1 - pozicijaOznake[(u, v)]) * pos[u][0] + pozicijaOznake[(u, v)] * pos[v][0] + dxdy[(u, v)][0],
                     (1 - pozicijaOznake[(u, v)]) * pos[u][1] + pozicijaOznake[(u, v)] * pos[v][1] + dxdy[(u, v)][1],
                     rjecnik[(u, v)], bbox=dict(facecolor='white', edgecolor='white', alpha=1), size=fontTezine)
    if rubVrha == None:
        nx.draw_networkx_nodes(G, pos, nodelist=vrhovi, node_size=velicinaVrha, node_color=bojaVrha,
                               node_shape=oblikVrha)
    else:
        nx.draw_networkx_nodes(G, pos, nodelist=vrhovi, node_size=velicinaVrha, node_color=bojaVrha,
                               node_shape=oblikVrha, edgecolors=rubVrha)
    nx.draw_networkx_edges(G, pos, width=debljinaBrida, edge_color=bojaBrida)
    nx.draw_networkx_labels(G, pos, font_size=fontVrh)
    # isticanje stabla najkracih putova
    istaknutiBridovi = []
    if tezine:
        stablo = nx.dijkstra_predecessor_and_distance(G, pocetak)[0]
        del stablo[pocetak]
        for u in stablo:
            istaknutiBridovi.append((u, stablo[u][0]))
    else:
        stablo = nx.predecessor(G, pocetak)
        del stablo[pocetak]
        for u in stablo:
            istaknutiBridovi.append((u, stablo[u][0]))
    nx.draw_networkx_edges(G, pos, edgelist=istaknutiBridovi, width=debljinaPuta, edge_color=bojaPuta, alpha=alfa)
    if rubVrha == None:
        nx.draw_networkx_nodes(G, pos, nodelist=[pocetak], node_size=velicinaVrha, node_color=bojaVrha2,
                               node_shape=oblikVrha)
    else:
        nx.draw_networkx_nodes(G, pos, nodelist=[pocetak], node_size=velicinaVrha, node_color=bojaVrha2,
                               node_shape=oblikVrha, edgecolors=rubVrha)
    plt.axis('off')
    pylab.gcf().set_facecolor('w')
    if slika == None:
        plt.axis('auto')
    else:
        plt.xlim(*slika[0])
        plt.ylim(*slika[1])
    plt.show()
    return


def RootedEmbedding(G, korijen):
    root = korijen
    n = G.number_of_nodes()
    e = dict(G.adjacency())
    pos = dict(zip(G.nodes(), [(-np.ceil(np.sqrt(n)), np.ceil(np.sqrt(n))) for i in range(n)]))
    v = dict(zip(G.nodes(), [(0, 0) for i in range(n)]))
    next = [root]
    done = [root]
    while next != []:
        root = next[0]
        new = []
        for i in e[root]:
            if not (i in done): new.append(i)
        for i in range(len(new)):
            dx = (pos[root][1] - pos[root][0]) / len(new)
            pos[new[i]] = (pos[root][0] + i * dx, pos[root][0] + (i + 1) * dx)
            x = sum(pos[new[i]]) / 2.0
            v[new[i]] = (x, v[root][1] - 1)
        next = next[1:] + new
        done = done + new
    return v


def KorijenskoStabloNajkracihPutova(G, korijen, tezine=True, omjeri=None, pomaci=None, slika=None, velicinaVrha=300,
                                    bojaVrha='y', oblikVrha='o', rubVrha=None, debljinaBrida=1, bojaBrida='k',
                                    fontTezine=12, fontVrh=12):
    istaknutiBridovi = []
    if tezine:
        stablo = nx.dijkstra_predecessor_and_distance(G, korijen)[0]
    else:
        stablo = nx.predecessor(G, korijen)
    del stablo[korijen]
    for u in stablo:
        if (u, stablo[u][0]) in G.edges():
            istaknutiBridovi.append((u, stablo[u][0]))
        else:
            istaknutiBridovi.append((stablo[u][0], u))
    G1 = nx.Graph(istaknutiBridovi)
    pos = RootedEmbedding(G1, korijen)
    pozicijaOznake = {}
    if tezine:
        if omjeri == None:
            for u in G.edges().__iter__():
                pozicijaOznake[u] = 0.5
        else:
            pozicijaOznake = omjeri
            for u in G.edges().__iter__():
                if not (u in pozicijaOznake.keys()):
                    pozicijaOznake[u] = 0.5
        dxdy = {}
        if pomaci == None:
            for u in G.edges().__iter__():
                dxdy[u] = (0, 0)
        else:
            dxdy = pomaci
            for u in G.edges().__iter__():
                if not (u in dxdy.keys()):
                    dxdy[u] = (0, 0)
        rjecnik = {}
        for u in G.edges(data=True).__iter__():
            rjecnik[(u[0], u[1])] = str(u[2]['weight'])
        for (u, v) in G1.edges().__iter__():
            try:
                plt.text(
                    (1 - pozicijaOznake[(u, v)]) * pos[u][0] + pozicijaOznake[(u, v)] * pos[v][0] + dxdy[(u, v)][0],
                    (1 - pozicijaOznake[(u, v)]) * pos[u][1] + pozicijaOznake[(u, v)] * pos[v][1] + dxdy[(u, v)][1],
                    rjecnik[(u, v)], bbox=dict(facecolor='white', edgecolor='white', alpha=1), size=fontTezine)
            except:
                plt.text(
                    (1 - pozicijaOznake[(v, u)]) * pos[v][0] + pozicijaOznake[(v, u)] * pos[u][0] + dxdy[(v, u)][0],
                    (1 - pozicijaOznake[(v, u)]) * pos[v][1] + pozicijaOznake[(v, u)] * pos[u][1] + dxdy[(v, u)][1],
                    rjecnik[(v, u)], bbox=dict(facecolor='white', edgecolor='white', alpha=1), size=fontTezine)
    if rubVrha == None:
        nx.draw_networkx_nodes(G1, pos, node_size=velicinaVrha, node_color=bojaVrha, node_shape=oblikVrha)
    else:
        nx.draw_networkx_nodes(G1, pos, node_size=velicinaVrha, node_color=bojaVrha, node_shape=oblikVrha,
                               edgecolors=rubVrha)
    nx.draw_networkx_edges(G1, pos, width=debljinaBrida, edge_color=bojaBrida)
    nx.draw_networkx_labels(G1, pos, font_size=fontVrh)
    # pylab.gcf().gca().set_axis_off()
    plt.axis('off')
    pylab.gcf().set_facecolor('w')
    if slika == None:
        plt.axis('auto')
    else:
        plt.xlim(*slika[0])
        plt.ylim(*slika[1])
    plt.show()
    return


def Dijkstra_graf(G, v0, v1=None):
    """Poboljsana verzija Dijkstrinog algoritma za najkrace putove u tezinskom grafu G od vrha v0.
    Ako graf ima i negativnih tezina, tada algoritam vraca error.
    Ako je v1=None, na izlazu se daje uredjeni par od dva rjecnika. U prvom rjecniku su dane
    udaljenosti pojedinih vrhova od vrha v0, a u drugom neposredni prethodnici pojedinih vrhova
    na najkracem putu od vrha v0.
    U protivnom se daje najkraci put izmedju vrhova v0 i v1 ukoliko postoji."""
    MM = min(map(lambda x: x[2]['weight'], G.edges(data=True)))
    tezina_brida = {}
    for e in G.edges(data=True):
        if e[2] != {}:
            tezina_brida[(e[0], e[1])] = e[2]['weight']
            tezina_brida[(e[1], e[0])] = e[2]['weight']
        else:
            tezina_brida[(e[0], e[1])] = 1
            tezina_brida[(e[1], e[0])] = 1
    P = {}
    udaljenost = {}
    Q = set(G.nodes())
    for v in Q:
        if v == v0:
            udaljenost[v] = 0
            P[v] = None
        else:
            udaljenost[v] = np.inf
            P[v] = None
    while len(Q) > 0:
        m = min([udaljenost[u] for u in Q])
        if m == np.inf: break
        v = random.choice(list(filter(lambda x: udaljenost[x] == m, Q)))
        Q.remove(v)
        for u in set(G.neighbors(v)).intersection(Q):
            if udaljenost[u] > udaljenost[v] + tezina_brida[(v, u)]:
                udaljenost[u] = udaljenost[v] + tezina_brida[(v, u)]
                P[u] = v
    if MM < 0: print("error: Dijkstra ne radi ispravno na negativnim tezinama")
    if v1 == None: return (udaljenost, P)
    v = v1
    min_put = [v1]
    while v != v0 and v != None:
        v = P[v]
        min_put = [v] + min_put
    if v == None:
        return []
    else:
        return min_put


def to_latex_dijkstra(tablica):
    """Pomocna funkcija za latex prikaz tablice kod rucnih prikaza Dijkstrinih algoritama."""
    for i in range(len(tablica)):
        for j in range(len(tablica[0])):
            if type(tablica[i][j]) in (int, str):
                pass
            elif tablica[i][j][0] == '-':
                if tablica[i][j][1] == np.inf:
                    tablica[i][j] = '$(-,\\infty)$'
                else:
                    tablica[i][j] = '$(-,0)$'
            else:
                tablica[i][j] = '$(\\text{' + str(tablica[i][j][0]) + '},' + str(tablica[i][j][1]) + ')$'
    return tablica


def rucni_Dijkstra_graf(G, v0, poredak=None):
    """Na izlazu daje tablicu kakva se dobiva kod rucnog provodjenja poboljsane verzije
    Dijkstrinog algoritma  na tezinskom grafu s pocetnim vrhom v0.
    Opcija poredak=None daje poredak vrhova u tablici kako se oni javljaju u listi
    G.nodes(), a ako zelimo neki drugi poredak, onda mozemo sami navesti listu sa
    zeljenim poretkom."""
    if poredak == None:
        V = list(G.nodes())
    else:
        V = poredak
    koraci = {}
    tezina_brida = {}
    for e in G.edges(data=True):
        tezina_brida[(e[0], e[1])] = e[2]['weight']
        tezina_brida[(e[1], e[0])] = e[2]['weight']
    koraci[0] = {}
    for v in V:
        if v == v0:
            koraci[0][v] = ('-', 0)
        else:
            koraci[0][v] = ('-', np.inf)
    Q = set(G.nodes())
    k = 1
    while len(Q) > 1:
        m = min([koraci[k - 1][x][1] for x in Q])
        if m == np.inf:
            break
        else:
            koraci[k] = koraci[k - 1].copy()
            w = random.choice(list(filter(lambda x: koraci[k - 1][x][1] == m, Q)))
            Q.remove(w)
            koraci[k][w] = '*****'
            for v in set(G.neighbors(w)).intersection(Q):
                if koraci[k - 1][v][1] > koraci[k - 1][w][1] + tezina_brida[(w, v)]:
                    koraci[k][v] = (w, koraci[k - 1][w][1] + tezina_brida[(w, v)])
        k += 1
    tablica = []
    for i in range(len(koraci)):
        redak = [str(i)]
        for v in V:
            redak.append(koraci[i][v])
        tablica.append(redak)
    tablica = [['vrh | korak'] + V] + tablica
    broj_stupaca = len(koraci) + 1
    tablica = list(map(list, zip(*tablica)))
    tablica = to_latex_dijkstra(tablica)
    tablica = make_table(tablica)
    apply_theme('basic_both')
    set_global_style(align='center')
    return tablica


def FW_html(lista, redak, stupac, korak=' '):
    stupac = [korak] + stupac
    for k in range(len(lista)):
        lista[k] = [redak[k]] + lista[k]
    lista = [stupac] + lista
    for i in range(len(lista)):
        for j in range(len(lista[0])):
            if lista[i][j] != np.inf:
                lista[i][j] = '$' + str(lista[i][j]) + '$'
            else:
                lista[i][j] = '$\\infty$'
    lista = make_table(lista)
    apply_theme('basic_both')
    set_global_style(align='center')
    return lista


def FW(graf, step, ncol, redoslijed_vrhova=None):
    if redoslijed_vrhova == None:
        redoslijed_vrhova = list(graf.nodes())
    bridovi = graf.edges()
    matrica = [[0 for j in range(len(redoslijed_vrhova))] for i in range(len(redoslijed_vrhova))]
    for i in range(len(redoslijed_vrhova)):
        for j in range(i + 1, len(redoslijed_vrhova)):
            if (redoslijed_vrhova[i], redoslijed_vrhova[j]) in bridovi:
                matrica[i][j] = graf.get_edge_data(redoslijed_vrhova[i], redoslijed_vrhova[j])['weight']
                matrica[j][i] = graf.get_edge_data(redoslijed_vrhova[i], redoslijed_vrhova[j])['weight']
            elif (redoslijed_vrhova[j], redoslijed_vrhova[i]) in bridovi:
                matrica[i][j] = graf.get_edge_data(redoslijed_vrhova[j], redoslijed_vrhova[i])['weight']
                matrica[j][i] = graf.get_edge_data(redoslijed_vrhova[j], redoslijed_vrhova[i])['weight']
            else:
                matrica[i][j] = np.inf
                matrica[j][i] = np.inf
    koraci = {0: copy.deepcopy(matrica)}
    for k in range(len(redoslijed_vrhova)):
        for i in range(len(redoslijed_vrhova)):
            for j in range(len(redoslijed_vrhova)):
                matrica[i][j] = min(matrica[i][j], matrica[i][k] + matrica[k][j])
        koraci[k + 1] = copy.deepcopy(matrica)
    lista_tablica = []
    for t in step:
        lista_tablica.append(FW_html(koraci[t], redoslijed_vrhova, redoslijed_vrhova, 'k=' + str(t)))
    prikaz = '<table><tr style="background-color:white;">'
    for i in range(len(step)):
        prikaz = prikaz + '<td>' + lista_tablica[i]._repr_html_() + '</td>'
        if (i + 1) % ncol == 0:
            prikaz = prikaz + '</tr>'
            if i + 1 < len(step):
                prikaz = prikaz + '<tr style="background-color:white;">'
    if len(step) % ncol != 0:
        prikaz = prikaz + '</tr></table>'
    else:
        prikaz = prikaz + '</table>'
    return HTML(prikaz)


def BrojRazapinjucihStabala(H):
    if nx.number_connected_components(H) > 1: return 0
    G = H.copy()
    for (u, v) in list(G.edges()):
        if u == v:
            G.remove_edge(u, v)
    bridovi = list(G.edges())
    vrhovi = list(G.nodes())
    nu = len(vrhovi)
    mat = np.matrix([[0 for i in range(nu)] for j in range(nu)])
    for i in range(nu):
        for j in range(nu):
            if ((vrhovi[i], vrhovi[j]) in bridovi):
                m = bridovi.count((vrhovi[i], vrhovi[j]))
                mat[i, j] += m
                mat[j, i] += m
    for i in range(nu):
        for j in range(nu):
            if i == j:
                mat[i, j] = G.degree(vrhovi[i])
            else:
                mat[i, j] = -mat[i, j]
    mat2 = mat[1:, 1:]
    return int(np.round(np.linalg.det(mat2)))


def DFSstablo(G, pos, pocetak, slika=None, velicinaVrha=300, bojaVrha='y', oblikVrha='o', rubVrha=None, debljinaBrida=1,
              bojaBrida='k', fontVrh=12, alfa=0.5, bojaStabla='r', debljinaStabla=5, bojaVrha2='cyan'):
    vrhovi = list(G.nodes())
    vrhovi.remove(pocetak)
    dfs_bridovi = list(nx.dfs_edges(G, pocetak))
    if rubVrha == None:
        nx.draw_networkx_nodes(G, pos, nodelist=vrhovi, node_size=velicinaVrha, node_color=bojaVrha,
                               node_shape=oblikVrha)
        nx.draw_networkx_nodes(G, pos, nodelist=[pocetak], node_size=velicinaVrha, node_color=bojaVrha2,
                               node_shape=oblikVrha)
    else:
        nx.draw_networkx_nodes(G, pos, nodelist=vrhovi, node_size=velicinaVrha, node_color=bojaVrha,
                               node_shape=oblikVrha, edgecolors=rubVrha)
        nx.draw_networkx_nodes(G, pos, nodelist=[pocetak], node_size=velicinaVrha, node_color=bojaVrha2,
                               node_shape=oblikVrha, edgecolors=rubVrha)
    nx.draw_networkx_labels(G, pos, font_size=fontVrh)
    nx.draw_networkx_edges(G, pos, width=debljinaBrida, edge_color=bojaBrida)
    nx.draw_networkx_edges(G, pos, edgelist=dfs_bridovi, width=debljinaStabla, edge_color=bojaStabla, alpha=alfa)
    plt.axis('off')
    pylab.gcf().set_facecolor('w')
    if slika == None:
        plt.axis('auto')
    else:
        plt.xlim(*slika[0])
        plt.ylim(*slika[1])
    plt.show()
    return


def DFSkorijenskoStablo(G, pocetak, slika=None, velicinaVrha=300, bojaVrha='y', oblikVrha='o', rubVrha=None,
                        debljinaBrida=1, bojaBrida='k', fontVrh=12):
    dfs_stablo = nx.Graph(list(nx.dfs_edges(G, pocetak)))
    pozicije = RootedEmbedding(dfs_stablo, pocetak)
    if rubVrha == None:
        nx.draw_networkx_nodes(dfs_stablo, pozicije, node_size=velicinaVrha, node_color=bojaVrha, node_shape=oblikVrha)
    else:
        nx.draw_networkx_nodes(dfs_stablo, pozicije, node_size=velicinaVrha, node_color=bojaVrha, node_shape=oblikVrha,
                               edgecolors=rubVrha)
    nx.draw_networkx_labels(dfs_stablo, pozicije, font_size=fontVrh)
    nx.draw_networkx_edges(dfs_stablo, pozicije, width=debljinaBrida, edge_color=bojaBrida)
    plt.axis('off')
    pylab.gcf().set_facecolor('w')
    if slika == None:
        plt.axis('auto')
    else:
        plt.xlim(*slika[0])
        plt.ylim(*slika[1])
    plt.show()
    return


def BFSstablo(G, pos, pocetak, slika=None, velicinaVrha=300, bojaVrha='y', oblikVrha='o', rubVrha=None, debljinaBrida=1,
              bojaBrida='k', fontVrh=12, alfa=0.5, bojaStabla='r', debljinaStabla=5, bojaVrha2='cyan'):
    vrhovi = list(G.nodes())
    vrhovi.remove(pocetak)
    bfs_bridovi = list(nx.bfs_edges(G, pocetak))
    if rubVrha == None:
        nx.draw_networkx_nodes(G, pos, nodelist=vrhovi, node_size=velicinaVrha, node_color=bojaVrha,
                               node_shape=oblikVrha)
        nx.draw_networkx_nodes(G, pos, nodelist=[pocetak], node_size=velicinaVrha, node_color=bojaVrha2,
                               node_shape=oblikVrha)
    else:
        nx.draw_networkx_nodes(G, pos, nodelist=vrhovi, node_size=velicinaVrha, node_color=bojaVrha,
                               node_shape=oblikVrha, edgecolors=rubVrha)
        nx.draw_networkx_nodes(G, pos, nodelist=[pocetak], node_size=velicinaVrha, node_color=bojaVrha2,
                               node_shape=oblikVrha, edgecolors=rubVrha)
    nx.draw_networkx_labels(G, pos, font_size=fontVrh)
    nx.draw_networkx_edges(G, pos, width=debljinaBrida, edge_color=bojaBrida)
    nx.draw_networkx_edges(G, pos, edgelist=bfs_bridovi, width=debljinaStabla, edge_color=bojaStabla, alpha=alfa)
    plt.axis('off')
    pylab.gcf().set_facecolor('w')
    if slika == None:
        plt.axis('auto')
    else:
        plt.xlim(*slika[0])
        plt.ylim(*slika[1])
    plt.show()
    return


def BFSkorijenskoStablo(G, pocetak, slika=None, velicinaVrha=300, bojaVrha='y', oblikVrha='o', rubVrha=None,
                        debljinaBrida=1, bojaBrida='k', fontVrh=12):
    bfs_stablo = nx.Graph(list(nx.bfs_edges(G, pocetak)))
    pozicije = RootedEmbedding(bfs_stablo, pocetak)
    if rubVrha == None:
        nx.draw_networkx_nodes(bfs_stablo, pozicije, node_size=velicinaVrha, node_color=bojaVrha, node_shape=oblikVrha)
    else:
        nx.draw_networkx_nodes(bfs_stablo, pozicije, node_size=velicinaVrha, node_color=bojaVrha, node_shape=oblikVrha,
                               edgecolors=rubVrha)
    nx.draw_networkx_labels(bfs_stablo, pozicije, font_size=fontVrh)
    nx.draw_networkx_edges(bfs_stablo, pozicije, width=debljinaBrida, edge_color=bojaBrida)
    plt.axis('off')
    pylab.gcf().set_facecolor('w')
    if slika == None:
        plt.axis('auto')
    else:
        plt.xlim(*slika[0])
        plt.ylim(*slika[1])
    plt.show()
    return


def is_acyclic(G):
    if G.number_of_selfloops() > 0: return False
    if G.is_multigraph(): return False
    parent = {}
    for v in G.nodes():
        S = []
        R = [v]
        parent[v] = 'null'
        while len(R) != 0:
            x = R[0]
            S.append(x)
            del R[0]
            for y in G.neighbors(x):
                if y == parent[x]:
                    continue
                if not (y in S):
                    parent[y] = x
                    R.append(y)
                else:
                    return False
    return True


def rezni_bridovi(G):
    komponente = nx.connected_component_subgraphs(G)
    rezni = []
    for H in komponente:
        rezni_kandidati = nx.dfs_edges(H, list(H.nodes())[0])
        for brid in rezni_kandidati:
            graf = H.copy()
            graf.remove_edge(*brid)
            if nx.number_connected_components(graf) > nx.number_connected_components(H):
                rezni.append(brid)
    return rezni


def dfs_graf(G, v0):
    V = list(G.nodes())
    V.remove(v0)
    vrhovi = [v0] + V
    d = {}
    f = {}
    P = {}
    color = {}
    tree_edges = []
    back_edges = []
    for u in vrhovi:
        color[u] = "WHITE"
        P[u] = None
    time = 0
    for u in vrhovi:
        if color[u] == "WHITE":
            S = [u]
            while S != []:
                x = S[-1]
                del S[-1]
                if color[x] == "WHITE":
                    time += 1
                    d[x] = time
                    color[x] = "GRAY"
                    S.append(x)
                    for v in G.neighbors(x):
                        if color[v] == "WHITE":
                            P[v] = x
                            S.append(v)
                elif color[x] == "GRAY":
                    time += 1
                    f[x] = time
                    color[x] = "BLACK"
    redoslijed_vrhova = sorted(G.nodes(), key=lambda x: d[x])
    tree_edges = []
    for u in redoslijed_vrhova:
        if P[u] != None:
            tree_edges.append((u, P[u]))
    back_edges = []
    for (u, v) in G.edges():
        if ((u, v) in tree_edges) or ((v, u) in tree_edges):
            pass
        else:
            back_edges.append((u, v))
    return (redoslijed_vrhova, d, f, P, tree_edges, back_edges)


def DFS_orijentacija(G, v0):
    """Daje orijentaciju grafa G na temelju DFS algoritma koji je poceo
    s vrhom v0. Ukoliko graf nema reznih bridova, tada je dobivena
    ujedno i jaka orijentacija na grafu G. Na izlazu se daje digraf."""
    D = nx.DiGraph()
    lukovi = []
    vrijeme = dfs_graf(G, v0)[1]
    bridoviDFS = dfs_graf(G, v0)[4]
    for e in G.edges():
        if ((e[0], e[1]) in bridoviDFS) or ((e[1], e[0]) in bridoviDFS):
            if vrijeme[e[0]] < vrijeme[e[1]]:
                lukovi.append((e[0], e[1]))
            else:
                lukovi.append((e[1], e[0]))
        else:
            if vrijeme[e[0]] < vrijeme[e[1]]:
                lukovi.append((e[1], e[0]))
            else:
                lukovi.append((e[0], e[1]))
    D.add_edges_from(lukovi)
    return D


def random_cycle(G):
    v0 = random.choice(list(G.nodes()))
    alg = dfs_graf(G, v0)
    if alg[5] == []: return False
    dfs_forest = nx.Graph(alg[4])
    brid_ciklus = random.choice(alg[5])
    u1 = brid_ciklus[0]
    u2 = brid_ciklus[1]
    ciklus = nx.shortest_path(dfs_forest, u1, u2) + [u1]
    return ciklus


def random_cycle_graph(G, pos, slika=None, velicinaVrha=300, bojaVrha='y', oblikVrha='o', rubVrha=None,
                       debljinaBrida=1, bojaBrida='k', fontVrh=12, alfa=0.5, bojaCiklusa='r', debljinaCiklusa=5):
    vrhovi = G.nodes()
    ciklus = random_cycle(G)
    bridovi_ciklus = []
    for k in range(len(ciklus) - 1):
        bridovi_ciklus.append((ciklus[k], ciklus[k + 1]))
    bridovi_ciklus.append((ciklus[-1], ciklus[0]))
    if rubVrha == None:
        nx.draw_networkx_nodes(G, pos, nodelist=vrhovi, node_size=velicinaVrha, node_color=bojaVrha,
                               node_shape=oblikVrha)
    else:
        nx.draw_networkx_nodes(G, pos, nodelist=vrhovi, node_size=velicinaVrha, node_color=bojaVrha,
                               node_shape=oblikVrha, edgecolors=rubVrha)
    nx.draw_networkx_labels(G, pos, font_size=fontVrh)
    nx.draw_networkx_edges(G, pos, width=debljinaBrida, edge_color=bojaBrida)
    nx.draw_networkx_edges(G, pos, edgelist=bridovi_ciklus, width=debljinaCiklusa, edge_color=bojaCiklusa, alpha=alfa)
    plt.axis('off')
    pylab.gcf().set_facecolor('w')
    if slika == None:
        plt.axis('auto')
    else:
        plt.xlim(*slika[0])
        plt.ylim(*slika[1])
    plt.show()
    return


def KruskalStablo(G, pos, opcija='min', omjeri=None, pomaci=None, slika=None, velicinaVrha=300, bojaVrha='y',
                  oblikVrha='o', debljinaBrida=1, rubVrha=None,
                  bojaBrida='k', fontTezine=12, fontVrh=12, alfa=0.5, bojaStabla='r', debljinaStabla=5):
    # crtanje tezinskog grafa
    pozicijaOznake = {}
    if omjeri == None:
        for u in G.edges().__iter__():
            pozicijaOznake[u] = 0.5
    else:
        pozicijaOznake = omjeri
        for u in G.edges().__iter__():
            if not (u in pozicijaOznake.keys()):
                pozicijaOznake[u] = 0.5
    dxdy = {}
    if pomaci == None:
        for u in G.edges().__iter__():
            dxdy[u] = (0, 0)
    else:
        dxdy = pomaci
        for u in G.edges().__iter__():
            if not (u in dxdy.keys()):
                dxdy[u] = (0, 0)
    rjecnik = {}
    for u in G.edges(data=True).__iter__():
        rjecnik[(u[0], u[1])] = str(u[2]['weight'])
    for (u, v) in G.edges().__iter__():
        plt.text((1 - pozicijaOznake[(u, v)]) * pos[u][0] + pozicijaOznake[(u, v)] * pos[v][0] + dxdy[(u, v)][0],
                 (1 - pozicijaOznake[(u, v)]) * pos[u][1] + pozicijaOznake[(u, v)] * pos[v][1] + dxdy[(u, v)][1],
                 rjecnik[(u, v)], bbox=dict(facecolor='white', edgecolor='white', alpha=1), size=fontTezine)
    if rubVrha == None:
        nx.draw_networkx_nodes(G, pos, node_size=velicinaVrha, node_color=bojaVrha, node_shape=oblikVrha)
    else:
        nx.draw_networkx_nodes(G, pos, node_size=velicinaVrha, node_color=bojaVrha, node_shape=oblikVrha,
                               edgecolors=rubVrha)
    nx.draw_networkx_edges(G, pos, width=debljinaBrida, edge_color=bojaBrida)
    nx.draw_networkx_labels(G, pos, font_size=fontVrh)
    # isticanje minimalnog stabla
    istaknutiBridovi = []
    if opcija == 'min':
        lista = list(nx.minimum_spanning_edges(G, data=False))
    if opcija == 'max':
        lista = list(nx.maximum_spanning_edges(G, data=False))
    # for u in lista:
    #    istaknutiBridovi.append((u[0],u[1]))
    nx.draw_networkx_edges(G, pos, edgelist=lista, width=debljinaStabla, edge_color=bojaStabla, alpha=alfa)
    plt.axis('off')
    pylab.gcf().set_facecolor('w')
    if slika == None:
        plt.axis('auto')
    else:
        plt.xlim(*slika[0])
        plt.ylim(*slika[1])
    plt.show()
    return


def PrimStablo(G, v0, pos, opcija='min', omjeri=None, pomaci=None, slika=None, velicinaVrha=300, bojaVrha='y',
               oblikVrha='o', rubVrha=None,
               debljinaBrida=1, bojaBrida='k', fontTezine=12, fontVrh=12, alfa=0.5, bojaStabla='r', debljinaStabla=5):
    # crtanje tezinskog grafa
    pozicijaOznake = {}
    if omjeri == None:
        for u in G.edges().__iter__():
            pozicijaOznake[u] = 0.5
    else:
        pozicijaOznake = omjeri
        for u in G.edges().__iter__():
            if not (u in pozicijaOznake.keys()):
                pozicijaOznake[u] = 0.5
    dxdy = {}
    if pomaci == None:
        for u in G.edges().__iter__():
            dxdy[u] = (0, 0)
    else:
        dxdy = pomaci
        for u in G.edges().__iter__():
            if not (u in dxdy.keys()):
                dxdy[u] = (0, 0)
    rjecnik = {}
    for u in G.edges(data=True).__iter__():
        rjecnik[(u[0], u[1])] = str(u[2]['weight'])
    for (u, v) in G.edges().__iter__():
        plt.text((1 - pozicijaOznake[(u, v)]) * pos[u][0] + pozicijaOznake[(u, v)] * pos[v][0] + dxdy[(u, v)][0],
                 (1 - pozicijaOznake[(u, v)]) * pos[u][1] + pozicijaOznake[(u, v)] * pos[v][1] + dxdy[(u, v)][1],
                 rjecnik[(u, v)], bbox=dict(facecolor='white', edgecolor='white', alpha=1), size=fontTezine)
    if rubVrha == None:
        nx.draw_networkx_nodes(G, pos, node_size=velicinaVrha, node_color=bojaVrha, node_shape=oblikVrha)
    else:
        nx.draw_networkx_nodes(G, pos, node_size=velicinaVrha, node_color=bojaVrha, node_shape=oblikVrha,
                               edgecolors=rubVrha)
    nx.draw_networkx_edges(G, pos, width=debljinaBrida, edge_color=bojaBrida)
    nx.draw_networkx_labels(G, pos, font_size=fontVrh)
    # isticanje minimalnog stabla
    istaknutiBridovi = []
    if opcija == 'min':
        lista = list(map(lambda x: (x[0], x[1]), prim_alg_min(G, v0)))
    if opcija == 'max':
        lista = list(map(lambda x: (x[0], x[1]), prim_alg_max(G, v0)))
    # for u in lista:
    #    istaknutiBridovi.append((u[0],u[1]))
    nx.draw_networkx_edges(G, pos, edgelist=lista, width=debljinaStabla, edge_color=bojaStabla, alpha=alfa)
    plt.axis('off')
    pylab.gcf().set_facecolor('w')
    if slika == None:
        plt.axis('auto')
    else:
        plt.xlim(*slika[0])
        plt.ylim(*slika[1])
    plt.show()
    return


def KruskalKorijenskoStablo(G, korijen, opcija='min', omjeri=None, pomaci=None, slika=None, velicinaVrha=300,
                            bojaVrha='y', oblikVrha='o', rubVrha=None,
                            debljinaBrida=1, bojaBrida='k', fontTezine=12, fontVrh=12):
    istaknutiBridovi = []
    if opcija == 'min':
        stablo = list(nx.minimum_spanning_edges(G, data=False))
    if opcija == 'max':
        stablo = list(nx.maximum_spanning_edges(G, data=False))
    # for u in stablo:
    #    istaknutiBridovi.append((u[0],u[1]))
    G1 = nx.Graph(stablo)
    pos = RootedEmbedding(G1, korijen)
    pozicijaOznake = {}
    if omjeri == None:
        for u in G.edges().__iter__():
            pozicijaOznake[u] = 0.5
    else:
        pozicijaOznake = omjeri
        for u in G.edges().__iter__():
            if not (u in pozicijaOznake.keys()):
                pozicijaOznake[u] = 0.5
    dxdy = {}
    if pomaci == None:
        for u in G.edges().__iter__():
            dxdy[u] = (0, 0)
    else:
        dxdy = pomaci
        for u in G.edges().__iter__():
            if not (u in dxdy.keys()):
                dxdy[u] = (0, 0)
    rjecnik = {}
    for u in G.edges(data=True).__iter__():
        rjecnik[(u[0], u[1])] = str(u[2]['weight'])
    for (u, v) in G1.edges().__iter__():
        try:
            plt.text((1 - pozicijaOznake[(u, v)]) * pos[u][0] + pozicijaOznake[(u, v)] * pos[v][0] + dxdy[(u, v)][0],
                     (1 - pozicijaOznake[(u, v)]) * pos[u][1] + pozicijaOznake[(u, v)] * pos[v][1] + dxdy[(u, v)][1],
                     rjecnik[(u, v)], bbox=dict(facecolor='white', edgecolor='white', alpha=1), size=fontTezine)
        except:
            plt.text((1 - pozicijaOznake[(v, u)]) * pos[v][0] + pozicijaOznake[(v, u)] * pos[u][0] + dxdy[(v, u)][0],
                     (1 - pozicijaOznake[(v, u)]) * pos[v][1] + pozicijaOznake[(v, u)] * pos[u][1] + dxdy[(v, u)][1],
                     rjecnik[(v, u)], bbox=dict(facecolor='white', edgecolor='white', alpha=1), size=fontTezine)
    if rubVrha == None:
        nx.draw_networkx_nodes(G1, pos, node_size=velicinaVrha, node_color=bojaVrha, node_shape=oblikVrha)
    else:
        nx.draw_networkx_nodes(G1, pos, node_size=velicinaVrha, node_color=bojaVrha, node_shape=oblikVrha,
                               edgecolors=rubVrha)
    nx.draw_networkx_edges(G1, pos, width=debljinaBrida, edge_color=bojaBrida)
    nx.draw_networkx_labels(G1, pos, font_size=fontVrh)
    # pylab.gcf().gca().set_axis_off()
    plt.axis('off')
    pylab.gcf().set_facecolor('w')
    if slika == None:
        plt.axis('auto')
    else:
        plt.xlim(*slika[0])
        plt.ylim(*slika[1])
    plt.show()
    return


def PrimKorijenskoStablo(G, korijen, opcija='min', omjeri=None, pomaci=None, slika=None, velicinaVrha=300, bojaVrha='y',
                         oblikVrha='o', rubVrha=None,
                         debljinaBrida=1, bojaBrida='k', fontTezine=12, fontVrh=12):
    istaknutiBridovi = []
    if opcija == 'min':
        stablo = map(lambda x: (x[0], x[1]), prim_alg_min(G, korijen))
    if opcija == 'max':
        stablo = map(lambda x: (x[0], x[1]), prim_alg_max(G, korijen))
    # for u in stablo:
    #    istaknutiBridovi.append((u[0],u[1]))
    G1 = nx.Graph(stablo)
    pos = RootedEmbedding(G1, korijen)
    pozicijaOznake = {}
    if omjeri == None:
        for u in G.edges().__iter__():
            pozicijaOznake[u] = 0.5
    else:
        pozicijaOznake = omjeri
        for u in G.edges().__iter__():
            if not (u in pozicijaOznake.keys()):
                pozicijaOznake[u] = 0.5
    dxdy = {}
    if pomaci == None:
        for u in G.edges().__iter__():
            dxdy[u] = (0, 0)
    else:
        dxdy = pomaci
        for u in G.edges().__iter__():
            if not (u in dxdy.keys()):
                dxdy[u] = (0, 0)
    rjecnik = {}
    for u in G.edges(data=True).__iter__():
        rjecnik[(u[0], u[1])] = str(u[2]['weight'])
    for (u, v) in G1.edges().__iter__():
        try:
            plt.text((1 - pozicijaOznake[(u, v)]) * pos[u][0] + pozicijaOznake[(u, v)] * pos[v][0] + dxdy[(u, v)][0],
                     (1 - pozicijaOznake[(u, v)]) * pos[u][1] + pozicijaOznake[(u, v)] * pos[v][1] + dxdy[(u, v)][1],
                     rjecnik[(u, v)], bbox=dict(facecolor='white', edgecolor='white', alpha=1), size=fontTezine)
        except:
            plt.text((1 - pozicijaOznake[(v, u)]) * pos[v][0] + pozicijaOznake[(v, u)] * pos[u][0] + dxdy[(v, u)][0],
                     (1 - pozicijaOznake[(v, u)]) * pos[v][1] + pozicijaOznake[(v, u)] * pos[u][1] + dxdy[(v, u)][1],
                     rjecnik[(v, u)], bbox=dict(facecolor='white', edgecolor='white', alpha=1), size=fontTezine)
    if rubVrha == None:
        nx.draw_networkx_nodes(G1, pos, node_size=velicinaVrha, node_color=bojaVrha, node_shape=oblikVrha)
    else:
        nx.draw_networkx_nodes(G1, pos, node_size=velicinaVrha, node_color=bojaVrha, node_shape=oblikVrha,
                               edgecolors=rubVrha)
    nx.draw_networkx_edges(G1, pos, width=debljinaBrida, edge_color=bojaBrida)
    nx.draw_networkx_labels(G1, pos, font_size=fontVrh)
    # pylab.gcf().gca().set_axis_off()
    plt.axis('off')
    pylab.gcf().set_facecolor('w')
    if slika == None:
        plt.axis('auto')
    else:
        plt.xlim(*slika[0])
        plt.ylim(*slika[1])
    plt.show()
    return


def prim_alg_min(G, v1):
    vrhovi = list(G.nodes())
    vrhovi.remove(v1)
    oznake_vrhova = {}
    bridovi_stablo = []
    odabrani_vrhovi = [v1]
    for v in vrhovi:
        if G.get_edge_data(v, v1) != None:
            oznake_vrhova[v] = (v1, G.get_edge_data(v, v1)['weight'])
        else:
            oznake_vrhova[v] = (None, np.inf)
    while vrhovi != []:
        # biranje vrha s najmanjom oznakom
        m = np.inf
        for x in vrhovi:
            if oznake_vrhova[x][1] < m:
                m = oznake_vrhova[x][1]
                u = x
        # dodavanje brida
        bridovi_stablo.append((oznake_vrhova[u][0], u, G.get_edge_data(oznake_vrhova[u][0], u)['weight']))
        # promjena oznaka susjednim neodabranim vrhovima
        for x in G.neighbors(u):
            if (x in oznake_vrhova) and (G.get_edge_data(u, x)['weight'] < oznake_vrhova[x][1]):
                oznake_vrhova[x] = (u, G.get_edge_data(u, x)['weight'])
        # brisanje  odabranog vrha
        vrhovi.remove(u)
    return bridovi_stablo


def prim_alg_max(H, v1):
    max_tezine = map(lambda x: (x[0], x[1], -x[2]['weight']), H.edges(data=True))
    G = nx.Graph()
    G.add_weighted_edges_from(max_tezine)
    vrhovi = list(G.nodes())
    vrhovi.remove(v1)
    oznake_vrhova = {}
    bridovi_stablo = []
    odabrani_vrhovi = [v1]
    for v in vrhovi:
        if G.get_edge_data(v, v1) != None:
            oznake_vrhova[v] = (v1, G.get_edge_data(v, v1)['weight'])
        else:
            oznake_vrhova[v] = (None, np.inf)
    while vrhovi != []:
        # biranje vrha s najmanjom oznakom
        m = np.inf
        for x in vrhovi:
            if oznake_vrhova[x][1] < m:
                m = oznake_vrhova[x][1]
                u = x
        # dodavanje brida
        bridovi_stablo.append((oznake_vrhova[u][0], u, H.get_edge_data(oznake_vrhova[u][0], u)['weight']))
        # promjena oznaka susjednim neodabranim vrhovima
        for x in G.neighbors(u):
            if (x in oznake_vrhova) and (G.get_edge_data(u, x)['weight'] < oznake_vrhova[x][1]):
                oznake_vrhova[x] = (u, G.get_edge_data(u, x)['weight'])
        # brisanje  odabranog vrha
        vrhovi.remove(u)
    return bridovi_stablo


def rucni_Kruskal(G, opcija='min'):
    if opcija == 'min':
        bridovi_stablo = list(nx.minimum_spanning_edges(G))
    else:
        bridovi_stablo = list(nx.maximum_spanning_edges(G))
    prvi_redak = ['korak'] + list(range(1, G.number_of_nodes()))
    drugi_redak = ['brid'] + list(map(lambda x: (x[0], x[1]), bridovi_stablo))
    treci_redak = ['tezina'] + list(map(lambda x: x[2]['weight'], bridovi_stablo))
    tezina = sum(map(lambda x: x[2]['weight'], bridovi_stablo))
    print("tezina stabla: ", tezina)
    tablica = make_table([prvi_redak, drugi_redak, treci_redak])
    apply_theme('basic_both')
    set_global_style(align='center')
    return tablica


def rucni_Prim(G, v0, opcija='min'):
    if opcija == 'min':
        bridovi_stablo = prim_alg_min(G, v0)
    else:
        bridovi_stablo = prim_alg_max(G, v0)
    prvi_redak = [str(v0) + ' | korak'] + list(range(1, G.number_of_nodes()))
    drugi_redak = ['brid'] + list(map(lambda x: (x[0], x[1]), bridovi_stablo))
    treci_redak = ['tezina'] + list(map(lambda x: x[2], bridovi_stablo))
    tezina = sum(map(lambda x: x[2], bridovi_stablo))
    print("tezina stabla: ", tezina)
    tablica = make_table([prvi_redak, drugi_redak, treci_redak])
    apply_theme('basic_both')
    set_global_style(align='center')
    return tablica


def BellmanFord(D, v0, v1=None):
    """Bellman-Fordov algoritam za najkrace putove u tezinskom digrafu D od vrha v0.
    Ako digraf ima usmjereni ciklus negativne tezine koji se moze doseci iz vrha v0,
    tada algoritam vraca error.
    Ako je v1=None, na izlazu se daje uredjeni par od dva rjecnika. U prvom rjecniku su
    dane udaljenosti pojedinih vrhova od vrha v0, a u drugom neposredni prethodnici
    pojedinih vrhova na najkracem putu od vrha v0.
    U protivnom se daje najkraci put izmedju vrhova v0 i v1 ukoliko postoji."""
    tezina_brida = {}
    for e in D.edges(data=True):
        if e[2] != {}:
            tezina_brida[(e[0], e[1])] = e[2]['weight']
        else:
            tezina_brida[(e[0], e[1])] = 1
    P = {}
    udaljenost = {}
    V = D.nodes()
    E = D.edges()
    udaljenost[0] = {}
    for v in V:
        if v == v0:
            udaljenost[0][v] = 0
            P[v] = None
        else:
            udaljenost[0][v] = np.inf
            P[v] = None
    for i in range(1, len(V) + 1):
        udaljenost[i] = udaljenost[i - 1].copy()
        promjena = False
        for v in V:
            for u in D.predecessors(v):
                if udaljenost[i][v] > udaljenost[i - 1][u] + tezina_brida[(u, v)]:
                    udaljenost[i][v] = udaljenost[i - 1][u] + tezina_brida[(u, v)]
                    P[v] = u
                    promjena = True
        if promjena == False:
            break
    if udaljenost[len(udaljenost) - 1] != udaljenost[len(udaljenost) - 2]: print("error: digraf ima negativne cikluse")
    if v1 == None: return (udaljenost[len(udaljenost) - 1], P)
    v = v1
    min_put = [v1]
    while v != v0 and v != None:
        v = P[v]
        min_put = [v] + min_put
    if v == None:
        return []
    else:
        return min_put


def Dijkstra_digraf(D, v0, v1=None):
    """Poboljsana verzija Dijkstrinog algoritma za najkrace putove u tezinskom digrafu D od vrha v0.
    Ako digraf ima i negativnih tezina, tada algoritam vraca error.
    Ako je v1=None, na izlazu se daje uredjeni par od dva rjecnika. U prvom rjecniku su dane
    udaljenosti pojedinih vrhova od vrha v0, a u drugom neposredni prethodnici pojedinih vrhova
    na najkracem putu od vrha v0.
    U protivnom se daje najkraci put izmedju vrhova v0 i v1 ukoliko postoji."""
    MM = min(map(lambda x: x[2]['weight'], D.edges(data=True)))
    tezina_brida = {}
    for e in D.edges(data=True):
        if e[2] != {}:
            tezina_brida[(e[0], e[1])] = e[2]['weight']
        else:
            tezina_brida[(e[0], e[1])] = 1
    P = {}
    udaljenost = {}
    Q = set(D.nodes())
    for v in Q:
        if v == v0:
            udaljenost[v] = 0
            P[v] = None
        else:
            udaljenost[v] = np.inf
            P[v] = None
    while len(Q) > 0:
        m = min([udaljenost[u] for u in Q])
        if m == np.inf: break
        v = random.choice(list(filter(lambda x: udaljenost[x] == m, Q)))
        Q.remove(v)
        for u in set(D.neighbors(v)).intersection(Q):
            if udaljenost[u] > udaljenost[v] + tezina_brida[(v, u)]:
                udaljenost[u] = udaljenost[v] + tezina_brida[(v, u)]
                P[u] = v
    if MM < 0: print("error: Dijkstra ne radi ispravno na negativnim tezinama")
    if v1 == None: return (udaljenost, P)
    v = v1
    min_put = [v1]
    while v != v0 and v != None:
        v = P[v]
        min_put = [v] + min_put
    if v == None:
        return []
    else:
        return min_put


def rucni_BellmanFord(D, v0, poredak=None):
    """Na izlazu daje tablicu kakva se dobiva kod rucnog provodjenja
    Bellman-Fordovog algoritma  na tezinskom digrafu s pocetnim vrhom v0.
    Opcija poredak=None daje poredak vrhova u tablici kako se oni javljaju u listi
    D.nodes(), a ako zelimo neki drugi poredak, onda mozemo sami navesti listu sa
    zeljenim poretkom."""
    koraci = {}
    if poredak == None:
        V = list(D.nodes())
    else:
        V = poredak
    tezina_brida = {}
    for e in D.edges(data=True):
        tezina_brida[(e[0], e[1])] = e[2]['weight']
    koraci[0] = {}
    for v in V:
        if v == v0:
            koraci[0][v] = ('-', 0)
        else:
            koraci[0][v] = ('-', np.inf)
    for i in range(1, len(V) + 1):
        koraci[i] = koraci[i - 1].copy()
        promjena = False
        for v in V:
            for u in D.predecessors(v):
                if koraci[i][v][1] > koraci[i - 1][u][1] + tezina_brida[(u, v)]:
                    koraci[i][v] = (u, koraci[i - 1][u][1] + tezina_brida[(u, v)])
                    promjena = True
        if promjena == False:
            break
    tablica = []
    for i in range(len(koraci)):
        redak = [str(i)]
        for v in V:
            redak.append(koraci[i][v])
        tablica.append(redak)
    tablica = [['vrh | korak'] + V] + tablica
    tablica = to_latex_dijkstra(list(map(list, zip(*tablica))))
    tablica = make_table(tablica)
    apply_theme('basic_both')
    set_global_style(align='center')
    if koraci[len(koraci) - 1] != koraci[len(koraci) - 2]: print("error: digraf ima negativne cikluse")
    return tablica


def rucni_Dijkstra_digraf(D, v0, poredak=None):
    """Na izlazu daje tablicu kakva se dobiva kod rucnog provodjenja poboljsane verzije
    Dijkstrinog algoritma  na tezinskom digrafu s pocetnim vrhom v0.
    Opcija poredak=None daje poredak vrhova u tablici kako se oni javljaju u listi
    D.nodes(), a ako zelimo neki drugi poredak, onda mozemo sami navesti listu sa
    zeljenim poretkom."""
    if poredak == None:
        V = list(D.nodes())
    else:
        V = poredak
    koraci = {}
    tezina_brida = {}
    for e in D.edges(data=True):
        tezina_brida[(e[0], e[1])] = e[2]['weight']
    MM = min(map(lambda x: x[2]['weight'], D.edges(data=True)))
    koraci[0] = {}
    for v in V:
        if v == v0:
            koraci[0][v] = ('-', 0)
        else:
            koraci[0][v] = ('-', np.inf)
    Q = set(D.nodes())
    k = 1
    while len(Q) > 1:
        m = min([koraci[k - 1][x][1] for x in Q])
        if m == np.inf:
            break
        else:
            koraci[k] = koraci[k - 1].copy()
            w = random.choice(list(filter(lambda x: koraci[k - 1][x][1] == m, Q)))
            Q.remove(w)
            koraci[k][w] = '*****'
            for v in set(D.neighbors(w)).intersection(Q):
                if koraci[k - 1][v][1] > koraci[k - 1][w][1] + tezina_brida[(w, v)]:
                    koraci[k][v] = (w, koraci[k - 1][w][1] + tezina_brida[(w, v)])
        k += 1
    tablica = []
    for i in range(len(koraci)):
        redak = [str(i)]
        for v in V:
            redak.append(koraci[i][v])
        tablica.append(redak)
    tablica = [['vrh | korak'] + V] + tablica
    tablica = to_latex_dijkstra(list(map(list, zip(*tablica))))
    tablica = make_table(tablica)
    apply_theme('basic_both')
    set_global_style(align='center')
    if MM < 0: print("error: Dijkstra ne radi ispravno na negativnim tezinama")
    return tablica


def graf_string(G, slika=None, fontV=12, fontE=12, xy=None, bojaVrha="white", bojaBrida="black", debljinaV=1,
                debljinaE=1, bojaTezine="black",
                bojeVrhova=None, bojeBridova=None, tezine=True, dekor=False, d={}, kut={}):
    """Pretvara networkx graf G u graphviz format.
    slika -> dimenzije nacrtane slike na izlazu
    fontV -> velicina fonta oznake vrhova
    fontE -> velicina fonta tezina bridova
    xy -> rjecnik koordinata vrhova grafa G
    bojaVrha -> boje svih vrhova na slici
    bojaBrida -> boje svih bridova na slici
    debljinaV -> debljina linije kod vrhova
    debljinaE -> debljina bridova
    bojaTezine -> boja u kojoj ce biti ispisane tezine bridova
    bojeVrhova -> rjecnik boja za svaki pojedini vrh ako zelimo da vrhovi budu razlicito obojani
    bojeBridova -> rjecnik boja za svaki pojedini brid ako zelimo bridove obojati u razlicitim bojama
    tezine -> da li ce tezine bridova biti ispisane ili nece
    dekor -> dekorirani ili nedekorirani ispis tezina
    d -> udaljenost tezine brida od njegovog kraja
    kut -> kut za koji je rotirana tezina brida oko njegovog kraja"""
    if G.is_directed():
        rijec = 'digraph {\n '
        strelica = '>'
    else:
        rijec = 'graph {\n '
        strelica = '-'
    if slika != None:
        rijec += 'size="{},{}!"\n '.format(*slika)
    rijec += 'node [shape=circle width=0 height=0 margin=0 style=filled fontsize={} penwidth={}];\n '.format(fontV,
                                                                                                             debljinaV)
    if dekor:
        rijec += 'edge[color={},fontcolor={} overlap=false splines=true decorate=true fontsize={} penwidth={}];\n '.format(
            bojaBrida, bojaTezine, fontE, debljinaE)
    else:
        rijec += 'edge[color={},fontcolor={} overlap=false splines=true,fontsize={} penwidth={}];\n '.format(bojaBrida,
                                                                                                             bojaTezine,
                                                                                                             fontE,
                                                                                                             debljinaE)
    if bojeVrhova == None:
        bojeVrhova = {}
        for v in G.nodes():
            bojeVrhova[v] = bojaVrha
    if xy != None:
        for v in G.nodes():
            rijec += '"{}" [label="{}" pos="{},{}!" fillcolor={}];\n '.format(v, v, xy[v][0], xy[v][1], bojeVrhova[v])
    else:
        for v in G.nodes():
            rijec += '"{}" [label="{}" style=filled fillcolor={}];\n '.format(v, v, bojeVrhova[v])
    if bojeBridova == None:
        bojeBridova = {}
        for e in G.edges():
            bojeBridova[(e[0], e[1])] = bojaBrida
    if tezine == False:
        for e in G.edges():
            rijec += '"{}" -{} "{}" [color={}];\n '.format(e[0], strelica, e[1], bojeBridova[(e[0], e[1])])
    if tezine == True:
        if d == {}:
            for e in G.edges(data=True):
                rijec += '"{}" -{} "{}" [color={} label="{}"];\n '.format(e[0], strelica, e[1],
                                                                          bojeBridova[(e[0], e[1])], e[2]['weight'])
        else:
            for e in G.edges(data=True):
                if (e[0], e[1]) in d.keys():
                    rijec += '"{}" -{} "{}" [color={} headlabel="{}" labeldistance={} labelangle={}];\n '.format(e[0],
                                                                                                                 strelica,
                                                                                                                 e[1],
                                                                                                                 bojeBridova[
                                                                                                                     (e[
                                                                                                                          0],
                                                                                                                      e[
                                                                                                                          1])],
                                                                                                                 e[2][
                                                                                                                     'weight'],
                                                                                                                 d[(
                                                                                                                 e[0],
                                                                                                                 e[1])],
                                                                                                                 kut[(
                                                                                                                 e[0],
                                                                                                                 e[1])])
                else:
                    rijec += '"{}" -{} "{}" [color={} label="{}"];\n '.format(e[0], strelica, e[1],
                                                                              bojeBridova[(e[0], e[1])], e[2]['weight'])
    rijec += '}'
    return rijec


def crtaj_graphviz(G, ime, slika=None, fontV=12, fontE=12, xy=None, bojaVrha="white", bojaBrida="black", debljinaV=1,
                   debljinaE=1, bojaTezine="black",
                   bojeVrhova=None, bojeBridova=None, tezine=True, dekor=False, d={}, kut={}):
    """Crta networkx graf G pomocu graphviza.
    ime -> ime datoteke u koju ce biti spremljena slika
    slika -> dimenzije nacrtane slike na izlazu
    fontV -> velicina fonta oznake vrhova
    fontE -> velicina fonta tezina bridova
    xy -> rjecnik koordinata vrhova grafa G
    bojaVrha -> boje svih vrhova na slici
    bojaBrida -> boje svih bridova na slici
    debljinaV -> debljina linije kod vrhova
    debljinaE -> debljina bridova
    bojaTezine -> boja u kojoj ce biti ispisane tezine bridova
    bojeVrhova -> rjecnik boja za svaki pojedini vrh ako zelimo da vrhovi budu razlicito obojani
    bojeBridova -> rjecnik boja za svaki pojedini brid ako zelimo bridove obojati u razlicitim bojama
    tezine -> da li ce tezine bridova biti ispisane ili nece
    dekor -> dekorirani ili nedekorirani ispis tezina
    d -> udaljenost tezine brida od njegovog kraja
    kut -> kut za koji je rotirana tezina brida oko njegovog kraja"""
    rijec = graf_string(G, slika, fontV, fontE, xy, bojaVrha, bojaBrida, debljinaV, debljinaE, bojaTezine, bojeVrhova,
                        bojeBridova, tezine, dekor, d, kut)
    os.system("echo '{}' | dot -Kfdp -Tpng > {}".format(rijec, ime))
    return Image(filename=ime)


def graphviz_najkraci_put(D, v1, v2, ime, algoritam=BellmanFord, slika=None, fontV=12, fontE=12, debljinaV=1,
                          debljinaE=1, xy=None, vrh0="yellow",
                          vrh1="pink", brid0="red", brid1="black", bojaTezine="black", tezine=True, dekor=False, d={},
                          kut={}):
    """Istice najkraci put od vrha v1 do vrha v2 u tezinskom digrafu (ili grafu) D na
    kojeg je primijenjen neki od nasih implementiranih algoritama.

    Opcija 'algoritam' moze imati neku od tri vrijednosti: BellmanFord, Dijkstra_digraf, Dijkstra_graf.

    Slika se crta u graphviz formatu.
    ime -> ime datoteke u koju ce biti spremljena slika
    slika -> dimenzije nacrtane slike na izlazu
    fontV -> velicina fonta oznake vrhova
    fontE -> velicina fonta tezina bridova
    debljinaV -> debljina linije kod vrhova
    debljinaE -> debljina bridova
    xy -> rjecnik koordinata vrhova grafa G
    vrh0 -> boja vrhova v1 i v2
    vrh1 -> boja preostalih vrhova
    brid0 -> boja bridova koji pripadaju najkracem putu
    brid1 -> boja bridova koji ne pripadaju najkracem putu
    bojaTezine -> boja u kojoj ce biti ispisane tezine bridova
    tezine -> da li ce tezine bridova biti ispisane ili nece
    dekor -> dekorirani ili nedekorirani ispis tezina
    d -> udaljenost tezine brida od njegovog kraja
    kut -> kut za koji je rotirana tezina brida oko njegovog kraja"""
    put = algoritam(D, v1, v2)
    ostali_vrhovi = list(set(D.nodes()).difference(set([v1, v2])))
    bridovi_put = []
    for i in range(len(put) - 1):
        bridovi_put.append((put[i], put[i + 1]))
    bojeVrhova = {v1: vrh0, v2: vrh0}
    for v in ostali_vrhovi:
        bojeVrhova[v] = vrh1
    bojeBridova = {}
    for e in D.edges():
        if ((e[0], e[1]) in bridovi_put) or ((e[1], e[0]) in bridovi_put):
            bojeBridova[e] = brid0
        else:
            bojeBridova[e] = brid1
    rijec = graf_string(D, slika, fontV, fontE, xy, vrh1, brid1, debljinaV, debljinaE, bojaTezine, bojeVrhova,
                        bojeBridova, tezine, dekor, d, kut)
    os.system("echo '{}' | dot -Kfdp -Tpng > {}".format(rijec, ime))
    return Image(filename=ime)


def graphviz_stablo_min_putova(D, v0, ime, algoritam=BellmanFord, slika=None, fontV=12, fontE=12, xy=None, debljinaV=1,
                               debljinaE=1, vrh0="yellow",
                               vrh1="pink", brid0="red", brid1="black", bojaTezine="black", tezine=True, dekor=False,
                               d={}, kut={}):
    """Istice stablo najkracih putova od vrha v0 prema svim preostalim vrhovima na tezinskom digrafu (ili grafu) na
    kojeg je primijenjen neki od nasih implementiranih algoritama.

    Opcija 'algoritam' moze imati neku od tri vrijednosti: BellmanFord, Dijkstra_digraf, Dijkstra_graf.

    Slika se crta u graphviz formatu.
    ime -> ime datoteke u koju ce biti spremljena slika
    slika -> dimenzije nacrtane slike na izlazu
    fontV -> velicina fonta oznake vrhova
    fontE -> velicina fonta tezina bridova
    xy -> rjecnik koordinata vrhova grafa G
    debljinaV -> debljina linije kod vrhova
    debljinaE -> debljina bridova
    vrh0 -> boje vrha v0
    vrh1 -> boja preostalih vrhova
    brid0 -> boja bridova koji pripadaju stablu
    brid1 -> boja bridova koji ne pripadaju stablu
    bojaTezine -> boja u kojoj ce biti ispisane tezine bridova
    tezine -> da li ce tezine bridova biti ispisane ili nece
    dekor -> dekorirani ili nedekorirani ispis tezina
    d -> udaljenost tezine brida od njegovog kraja
    kut -> kut za koji je rotirana tezina brida oko njegovog kraja"""
    P = algoritam(D, v0)[1]
    bridovi = []
    for x in P:
        if P[x] != None: bridovi.append((P[x], x))
    ostali_vrhovi = list(D.nodes())
    ostali_vrhovi.remove(v0)
    bojeVrhova = {v0: vrh0}
    for v in ostali_vrhovi:
        bojeVrhova[v] = vrh1
    bojeBridova = {}
    for e in D.edges():
        if ((e[0], e[1]) in bridovi) or ((e[1], e[0]) in bridovi):
            bojeBridova[e] = brid0
        else:
            bojeBridova[e] = brid1
    rijec = graf_string(D, slika, fontV, fontE, xy, vrh1, brid1, debljinaV, debljinaE, bojaTezine, bojeVrhova,
                        bojeBridova, tezine, dekor, d, kut)
    os.system("echo '%s' | dot -Kfdp -Tpng > %s" % (rijec, ime))
    return Image(filename=ime)


def usmjereni_ciklus(D):
    """Daje neki usmjereni ciklus u digrafu ukoliko digraf nije aciklicki.
    Zapravo je implementiran Brentov algoritam na nacin da ako vise
    puta primijenimo algoritam na istom digrafu, opcenito ce se dobiti
    na izlazu razliciti usmjereni ciklusi."""
    if nx.is_directed_acyclic_graph(D): return "nema usmjerenih ciklusa"
    x0 = random.choice(list(D.nodes()))
    f = {}
    for v in D.nodes():
        f[v] = random.choice(list(D.neighbors(v)))
    # racunanje lambde
    power = 1
    lam = 1
    tortoise = x0
    hare = f[x0]
    while tortoise != hare:
        if power == lam:
            tortoise = hare
            power *= 2
            lam = 0
        hare = f[hare]
        lam += 1
    # mu = pozicija prvog ponavljanja duljine lambda
    mu = 0
    tortoise = hare = x0
    for i in range(lam):
        hare = f[hare]
    while tortoise != hare:
        tortoise = f[tortoise]
        hare = f[hare]
        mu += 1
    ciklus = [tortoise]
    for i in range(lam):
        tortoise = f[tortoise]
        ciklus.append(tortoise)
    return ciklus


def graphviz_usmjereni_ciklus(D, ime, slika=None, fontV=12, fontE=12, debljinaV=1, debljinaE=1, xy=None, vrh0="yellow",
                              vrh1="pink",
                              brid0="red", brid1="black", bojaTezine="black", tezine=False, dekor=False, d={}, kut={}):
    """Istice neki usmjereni ciklus u tezinskom digrafu (ili grafu) D na
    kojeg je primijenjen Brentov algoritam.

    Brentov algoritam je implementiran na nacin da ako vise
    puta primijenimo algoritam na istom digrafu, opcenito ce se dobiti
    na izlazu razliciti usmjereni ciklusi.

    Slika se crta u graphviz formatu.
    ime -> ime datoteke u koju ce biti spremljena slika
    slika -> dimenzije nacrtane slike na izlazu
    fontV -> velicina fonta oznake vrhova
    fontE -> velicina fonta tezina bridova
    debljinaV -> debljina linije kod vrhova
    debljinaE -> debljina bridova
    xy -> rjecnik koordinata vrhova grafa G
    vrh0 -> boja vrhova u usmjerenom ciklusu
    vrh1 -> boja preostalih vrhova
    brid0 -> boja bridova koji pripadaju usmjerenom ciklusu
    brid1 -> boja bridova koji ne pripadaju usmjerenom ciklusu
    bojaTezine -> boja u kojoj ce biti ispisane tezine bridova
    tezine -> da li ce tezine bridova biti ispisane ili nece
    dekor -> dekorirani ili nedekorirani ispis tezina
    d -> udaljenost tezine brida od njegovog kraja
    kut -> kut za koji je rotirana tezina brida oko njegovog kraja"""
    ciklus = usmjereni_ciklus(D)
    bridovi_ciklus = []
    for i in range(len(ciklus) - 1):
        bridovi_ciklus.append((ciklus[i], ciklus[i + 1]))
    bojeVrhova = {}
    for v in D.nodes():
        if v in ciklus:
            bojeVrhova[v] = vrh0
        else:
            bojeVrhova[v] = vrh1
    bojeBridova = {}
    for e in D.edges():
        if e in bridovi_ciklus:
            bojeBridova[e] = brid0
        else:
            bojeBridova[e] = brid1
    rijec = graf_string(D, slika, fontV, fontE, xy, vrh1, brid1, debljinaV, debljinaE, bojaTezine, bojeVrhova,
                        bojeBridova, tezine, dekor, d, kut)
    os.system("echo '%s' | dot -Kfdp -Tpng > %s" % (rijec, ime))
    return Image(filename=ime)


def snaga_vrha(T, v0):
    indeks = list(T.nodes()).index(v0)
    matrica = nx.adj_matrix(T) + nx.adj_matrix(T) ** 2
    return int(matrica[indeks].sum())


def rang_lista(turnir, izlaz="html"):
    rezultat1 = map(lambda v: [v, turnir.out_degree(v), snaga_vrha(turnir, v)], turnir.nodes())
    rezultat1 = sorted(rezultat1, key=lambda x: x[2], reverse=True)
    rezultat = list(map(lambda v: [str(v[0]), str(v[1]), str(v[2])], rezultat1))
    if izlaz == "lista":
        return rezultat1
    elif izlaz == "html":
        rezultat = [['igraci', 'broj pobjeda', 'snaga pobjeda']] + rezultat
        tablica = make_table(rezultat)
        apply_theme('basic_both')
        set_global_style(align='center')
        return tablica


def to_latex_protok(tablica):
    """Pomocna funkcija za latex prikaz tablice kod rucnih prikaza protoka."""
    for i in range(len(tablica)):
        for j in range(len(tablica[0])):
            if type(tablica[i][j]) in (int, str):
                pass
            else:
                tablica[i][j] = '$(\\text{' + str(tablica[i][j][0]) + '},\\text{' + str(tablica[i][j][1]) + '})$'
    return tablica


def maksimalni_protok(mreza, v0, v1, kapacitet='weight'):
    protok = nx.maximum_flow(mreza, v0, v1, capacity=kapacitet)
    kapacitet_luk = []
    protok_luk = []
    for e in mreza.edges(data=True):
        kapacitet_luk.append(str(e[2][kapacitet]))
        protok_luk.append(str(protok[1][e[0]][e[1]]))
    lukovi = ['luk'] + list(mreza.edges())
    kapacitet_luk = ['kapacitet luka'] + kapacitet_luk
    protok_luk = ['protok kroz luk'] + protok_luk
    tablica = [lukovi, kapacitet_luk, protok_luk]
    tablica = to_latex_protok(list(map(list, zip(*tablica))))
    tablica = make_table(tablica)
    apply_theme('basic_both')
    set_global_style(align='center')
    print("Vrijednost maksimalnog protoka: ", protok[0])
    return tablica


def minimalni_rez(mreza, i, p, kapacitet='weight'):
    rez = nx.minimum_cut(mreza, i, p, capacity=kapacitet)[1]
    bridovi_rez = []
    for u in rez[0]:
        for v in rez[1]:
            if (u, v) in mreza.edges():
                bridovi_rez.append((u, v))
    return {'particija': (rez[0], rez[1]), 'bridovi': bridovi_rez}


def graphviz_minimalni_rez(D, v1, v2, ime, kapacitet='weight', slika=None, fontV=12, fontE=12, debljinaV=1, debljinaE=1,
                           xy=None,
                           vrh0="yellow", vrh1="pink", brid0="red", brid1="black", bojaTezine="black", tezine=True,
                           dekor=False, d={}, kut={}):
    """Istice minimalni (v1,v2)-rez u transportnoj mrezi D.

    Slika se crta u graphviz formatu.
    ime -> ime datoteke u koju ce biti spremljena slika
    slika -> dimenzije nacrtane slike na izlazu
    fontV -> velicina fonta oznake vrhova
    fontE -> velicina fonta tezina bridova
    debljinaV -> debljina linije kod vrhova
    debljinaE -> debljina bridova
    xy -> rjecnik koordinata vrhova grafa G
    vrh0 -> boja vrhova u clanu particije u kojemu je vrh v1
    vrh1 -> boja vrhova u clanu particije u kojemu je vrh v2
    brid0 -> boja bridova koji pripadaju minimalnom rezu
    brid1 -> boja bridova koji ne pripadaju minimalnom rezu
    bojaTezine -> boja u kojoj ce biti ispisane tezine bridova
    tezine -> da li ce tezine bridova biti ispisane ili nece
    dekor -> dekorirani ili nedekorirani ispis tezina
    d -> udaljenost tezine brida od njegovog kraja
    kut -> kut za koji je rotirana tezina brida oko njegovog kraja"""
    particija = minimalni_rez(D, v1, v2, kapacitet)
    bojeVrhova = {}
    for v in particija['particija'][0]:
        bojeVrhova[v] = vrh0
    for v in particija['particija'][1]:
        bojeVrhova[v] = vrh1
    bojeBridova = {}
    for e in D.edges():
        if e in particija['bridovi']:
            bojeBridova[e] = brid0
        else:
            bojeBridova[e] = brid1
    rijec = graf_string(D, slika, fontV, fontE, xy, vrh1, brid1, debljinaV, debljinaE, bojaTezine, bojeVrhova,
                        bojeBridova, tezine, dekor, d, kut)
    os.system("echo '%s' | dot -Kfdp -Tpng > %s" % (rijec, ime))
    return Image(filename=ime)


def kromatski_broj(G):
    H = G.copy()
    while True:
        nuH = H.number_of_nodes()
        vrhoviH = list(H.nodes())
        nesusjedni_vrhovi = []
        for i in range(nuH):
            for j in range(i + 1, nuH):
                if not ((vrhoviH[i], vrhoviH[j]) in H.edges()):
                    nesusjedni_vrhovi.append((vrhoviH[i], vrhoviH[j]))
        if len(nesusjedni_vrhovi) == 0:
            return nuH
        else:
            broj_zajednickih_susjeda = -1
            for par in nesusjedni_vrhovi:
                br = len(set(H.neighbors(par[0])).intersection(set(H.neighbors(par[1]))))
                if br > broj_zajednickih_susjeda:
                    broj_zajednickih_susjeda = br
                    par_vrhova = par
        H = nx.contracted_nodes(H, *par_vrhova)


def bojanje_vrhova(G, num=2000, odredi_KB=True, check=True):
    OPT = G.number_of_nodes()
    if (check) and (OPT > 120):
        return "Broj vrhova > 120. Stavi check=False ako ipak zelis probati."
    if odredi_KB:
        KR_BR = kromatski_broj(G)
        bojanje = nx.greedy_color(G, 'DSATUR')
        BROJ_BOJA = max(bojanje.values()) + 1
        if BROJ_BOJA == KR_BR:
            return bojanje
        else:
            OPT = BROJ_BOJA
            OPTIMALNO_BOJANJE = bojanje
        bojanje = nx.greedy_color(G, 'largest_first')
        BROJ_BOJA = max(bojanje.values()) + 1
        if BROJ_BOJA == KR_BR:
            return bojanje
        else:
            if BROJ_BOJA < OPT:
                OPT = BROJ_BOJA
                OPTIMALNO_BOJANJE = bojanje
        for i in range(num):
            bojanje = nx.greedy_color(G, 'random_sequential')
            BROJ_BOJA = max(bojanje.values()) + 1
            if BROJ_BOJA == KR_BR:
                return bojanje
            else:
                if BROJ_BOJA < OPT:
                    OPT = BROJ_BOJA
                    OPTIMALNO_BOJANJE = bojanje
        print('Bojanje koristi {} boja, a kromatski broj grafa je {}'.format(OPT, KR_BR))
        return OPTIMALNO_BOJANJE
    else:
        bojanje = nx.greedy_color(G, 'DSATUR')
        BROJ_BOJA = max(bojanje.values()) + 1
        OPT = BROJ_BOJA
        OPTIMALNO_BOJANJE = bojanje
        bojanje = nx.greedy_color(G, 'largest_first')
        BROJ_BOJA = max(bojanje.values()) + 1
        if BROJ_BOJA < OPT:
            OPT = BROJ_BOJA
            OPTIMALNO_BOJANJE = bojanje
        for i in range(num):
            bojanje = nx.greedy_color(G, 'random_sequential')
            BROJ_BOJA = max(bojanje.values()) + 1
            if BROJ_BOJA < OPT:
                OPT = BROJ_BOJA
                OPTIMALNO_BOJANJE = bojanje
        print('Dobiveno bojanje koristi {} boja i mozda nije optimalno.'.format(OPT))
        return OPTIMALNO_BOJANJE


def obojani_vrhovi(G, raspored_vrhova, velicina_vrha=300):
    bojanje = bojanje_vrhova(G)
    br = max(bojanje.values())
    cmap = matplotlib.cm.get_cmap('rainbow')
    boje = [cmap(i / br) for i in range(br + 1)]
    boje_vrhova = []
    for v in G.nodes():
        boje_vrhova.append(boje[bojanje[v]])
    nx.draw(G, pos=raspored_vrhova, node_size=velicina_vrha, edgecolors='k', with_labels=True, node_color=boje_vrhova)


def bojanje_vrha(graf, w, boje):
    obojaniSusjedi = list(set(nx.neighbors(graf, w)).intersection(set(boje.keys())))
    trazenjeBoje = []
    for u in obojaniSusjedi:
        trazenjeBoje.append(boje[u])
    if set(boje.values()).difference(set(trazenjeBoje)) != set([]):
        boje[w] = min(set(boje.values()).difference(set(trazenjeBoje)))
    else:
        boje[w] = max(boje.values()) + 1
    return None


def stupanj_zasicenosti(graf, w, boje):
    obojaniSusjedi = list(set(nx.neighbors(graf, w)).intersection(set(boje.keys())))
    brojBoja = len(set(map(lambda x: boje[x], obojaniSusjedi)))
    return brojBoja


def Welsh_Powell(graf):
    vrhovi = list(graf.nodes())
    vrhovi_sorted = []
    while vrhovi != []:
        D = max(dict(graf.degree(vrhovi)).values())
        max_vrh = []
        for v in vrhovi:
            if graf.degree(v) == D: max_vrh = max_vrh + [v]
        random.shuffle(max_vrh)
        vrhovi_sorted.extend(max_vrh)
        vrhovi = list(set(vrhovi).difference(set(max_vrh)))
    obojaniVrhovi = [vrhovi_sorted[0]]
    boje = {vrhovi_sorted[0]: 1}
    del vrhovi_sorted[0]
    while vrhovi_sorted != []:
        w = vrhovi_sorted[0]
        bojanje_vrha(graf, w, boje)
        obojaniVrhovi.append(w)
        del vrhovi_sorted[0]
    return (obojaniVrhovi, boje, len(set(boje.values())))


def Brelaz(graf):
    vrhovi = list(graf.nodes())
    D = max(dict(graf.degree()).values())
    max_vrh = ()
    for v in vrhovi:
        if graf.degree(v) == D: max_vrh = max_vrh + (v,)
    v1 = random.choice(max_vrh)
    obojaniVrhovi = [v1]
    boje = {v1: 1}
    del vrhovi[vrhovi.index(v1)]
    while vrhovi != []:
        # racunanje stupnjeva zasicenosti
        ds = {}
        for v in vrhovi:
            ds[v] = stupanj_zasicenosti(graf, v, boje)
        # izdvajanje vrhova s maksimalnim stupnjem zasicenosti
        ds_max = max(ds.values())
        ds_maxVrh = ()
        for v in ds:
            if ds[v] == ds_max: ds_maxVrh = ds_maxVrh + (v,)
        w = random.choice(ds_maxVrh)
        # bojanje vrha s maksimlnim stupnjem zasicenosti
        bojanje_vrha(graf, w, boje)
        obojaniVrhovi.append(w)
        del vrhovi[vrhovi.index(w)]
    return (obojaniVrhovi, boje, len(set(boje.values())))


def bridno_kromatski_broj(G):
    return kromatski_broj(nx.line_graph(G))


def bojanje_bridova(G, num=2000, odredi_BKB=True, check=True):
    H = nx.line_graph(G)
    OPT = H.number_of_nodes()
    if (check) and (OPT > 120):
        return "Broj bridova > 120. Stavi check=False ako ipak zelis probati."
    if odredi_BKB:
        BKR_BR = kromatski_broj(H)
        bojanje = nx.greedy_color(H, 'DSATUR')
        BROJ_BOJA = max(bojanje.values()) + 1
        if BROJ_BOJA == BKR_BR:
            return bojanje
        else:
            OPT = BROJ_BOJA
            OPTIMALNO_BOJANJE = bojanje
        bojanje = nx.greedy_color(H, 'largest_first')
        BROJ_BOJA = max(bojanje.values()) + 1
        if BROJ_BOJA == BKR_BR:
            return bojanje
        else:
            if BROJ_BOJA < OPT:
                OPT = BROJ_BOJA
                OPTIMALNO_BOJANJE = bojanje
        for i in range(num):
            bojanje = nx.greedy_color(H, 'random_sequential')
            BROJ_BOJA = max(bojanje.values()) + 1
            if BROJ_BOJA == BKR_BR:
                return bojanje
            else:
                if BROJ_BOJA < OPT:
                    OPT = BROJ_BOJA
                    OPTIMALNO_BOJANJE = bojanje
        print('Bojanje koristi {} boja, a bridno kromatski broj grafa je {}'.format(OPT, BKR_BR))
        return OPTIMALNO_BOJANJE
    else:
        bojanje = nx.greedy_color(H, 'DSATUR')
        BROJ_BOJA = max(bojanje.values()) + 1
        OPT = BROJ_BOJA
        OPTIMALNO_BOJANJE = bojanje
        bojanje = nx.greedy_color(H, 'largest_first')
        BROJ_BOJA = max(bojanje.values()) + 1
        if BROJ_BOJA < OPT:
            OPT = BROJ_BOJA
            OPTIMALNO_BOJANJE = bojanje
        for i in range(num):
            bojanje = nx.greedy_color(H, 'random_sequential')
            BROJ_BOJA = max(bojanje.values()) + 1
            if BROJ_BOJA < OPT:
                OPT = BROJ_BOJA
                OPTIMALNO_BOJANJE = bojanje
        print('Dobiveno bojanje koristi {} boja i mozda nije optimalno.'.format(OPT))
        return OPTIMALNO_BOJANJE


def obojani_bridovi(G, raspored_vrhova, velicina_vrha=300, boja_vrha='y', debljinaBrida=2):
    bojanje = bojanje_bridova(G)
    br = max(bojanje.values())
    cmap = matplotlib.cm.get_cmap('rainbow')
    boje = [cmap(i / br) for i in range(br + 1)]
    boje_bridova = []
    for e in G.edges():
        try:
            boje_bridova.append(boje[bojanje[(e[0], e[1])]])
        except:
            boje_bridova.append(boje[bojanje[(e[1], e[0])]])
    nx.draw(G, pos=raspored_vrhova, node_size=velicina_vrha, edgecolors='k', with_labels=True, node_color=boja_vrha,
            edge_color=boje_bridova, width=debljinaBrida)


def raspored_termina(G, stupac, particija=1, permutacija=None):
    BR = bridno_kromatski_broj(G)
    if len(stupac) != BR:
        return "Error: Broj termina nije jednak bridno kromatskom broju grafa"
    redak = list(nx.bipartite.sets(G)[particija - 1])
    bojanje = bojanje_bridova(G)
    if permutacija == None:
        permutacija = list(range(BR))
    else:
        permutacija = list(map(lambda x: x - 1, permutacija))
    tablica = [[" "] * len(redak) for _ in range(BR)]
    for brid in bojanje:
        if brid[0] in redak:
            tablica[permutacija[bojanje[brid]]][redak.index(brid[0])] = brid[1]
        else:
            tablica[permutacija[bojanje[brid]]][redak.index(brid[1])] = brid[0]
    for i in range(BR):
        tablica[i] = [stupac[i]] + tablica[i]
    tablica = [[" "] + redak] + tablica
    tablica = make_table(tablica)
    apply_theme('basic_both')
    set_global_style(align='center')
    return tablica


def random_maximal_matching(G, M1=[]):
    M = M1[:]
    bridovi = list(G.edges())
    provjera = []
    if M != []:
        for (u, v) in M:
            if not ((u, v) in bridovi):
                return "Error: U sparivanju M su nepostojeci bridovi."
        for (u, v) in M:
            if ((u in provjera) or (v in provjera)):
                return "Error: M nije sparivanje u grafu G."
            else:
                provjera.extend([u, v])
    bridovi = list(set(bridovi).difference(set(M)))
    while bridovi != []:
        b = random.choice(bridovi)
        if not ((b[0] in provjera) or (b[1] in provjera)):
            M.append(b)
            provjera.extend([b[0], b[1]])
        bridovi.remove(b)
    return M


def najvece_sparivanje_graf(G, raspored_vrhova, velicina_vrha=300, boja_vrha='y', debljinaBrida=7, bojaBrida='b',
                            slika=[0.2, 0.2]):
    sparivanje = nx.max_weight_matching(G, maxcardinality=True)
    plt.margins(x=slika[0], y=slika[1])
    plt.axis('off')
    nx.draw_networkx_nodes(G, pos=raspored_vrhova, node_size=velicina_vrha, node_color=boja_vrha, edgecolors='k')
    nx.draw_networkx_edges(G, pos=raspored_vrhova, width=debljinaBrida, edgelist=sparivanje, edge_color=bojaBrida,
                           alpha=0.5)
    nx.draw_networkx_edges(G, pos=raspored_vrhova)
    nx.draw_networkx_labels(G, pos=raspored_vrhova)


def maksimalno_sparivanje_graf(G, raspored_vrhova, M=[], velicina_vrha=300, boja_vrha='y', debljinaBrida=7,
                               bojaBridaM='r', bojaBrida='b', slika=[0.2, 0.2]):
    sparivanje = random_maximal_matching(G, M)
    if M != []:
        sparivanje = list(set(sparivanje).difference(set(M)))
    plt.margins(x=slika[0], y=slika[1])
    plt.axis('off')
    nx.draw_networkx_nodes(G, pos=raspored_vrhova, node_size=velicina_vrha, node_color=boja_vrha, edgecolors='k')
    nx.draw_networkx_edges(G, pos=raspored_vrhova, width=debljinaBrida, edgelist=sparivanje, edge_color=bojaBrida,
                           alpha=0.5)
    nx.draw_networkx_edges(G, pos=raspored_vrhova, width=debljinaBrida, edgelist=M, edge_color=bojaBridaM, alpha=0.5)
    nx.draw_networkx_edges(G, pos=raspored_vrhova)
    nx.draw_networkx_labels(G, pos=raspored_vrhova)