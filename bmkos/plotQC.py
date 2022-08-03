import numpy as np
import pandas as pd
from plotly import graph_objects as go
from plotly import express as px


def plotLen(qlen, lastLen=2000, wid=100):
    if qlen.max() < lastLen + 1:
        h = np.histogram(qlen, bins=range(0, qlen.max() + wid, wid))
        cus = [(i, i + wid) for i in h[1]]
        fig = go.Figure(go.Bar(x=h[1], y=h[0], offset=0, width=wid,
                               hovertemplate="ReadLength: %{customdata[0]}~%{customdata[1]}bp<br>Count: %{y}",
                               name='', customdata=cus))
        fig = fig.update_layout(xaxis={'tickvals': h[1][:-1]})
    else:
        h = np.histogram(qlen, bins=range(0, lastLen + 1, wid))
        cus = [(i, i + wid) for i in h[1]]
        cus[-1] = ('>' + str(int(lastLen)), '')
        h0 = np.append(h[0], (qlen > lastLen).sum())
        h1 = np.append(h[1], h[1][-1] + wid)
        fig = go.Figure(go.Bar(x=h1, y=h0, offset=0, width=wid,
                               hovertemplate="ReadLength: %{customdata[0]}~%{customdata[1]}bp<br>Count: %{y}",
                               name='', customdata=cus))
        fig = fig.update_layout(xaxis={'tickvals': h1[:-1]})
    return fig



def plotQscore(qscore, firstScore=6, lastScore=30):
    h = np.histogram(qscore, bins=range(firstScore, lastScore + 1, 1))
    y = np.concatenate([[(qscore < firstScore).sum(), ], h[0], [(qscore > lastScore).sum(), ]])
    x = np.concatenate([[firstScore - 1, ], h[1], [lastScore + 1, ]])
    cus = [(firstScore - 1, firstScore), ] + [(i, i + 1) for i in h[1]]
    cus[-1] = ('>' + str(int(lastScore)), '')
    fig = go.Figure(go.Bar(x=x, y=y, offset=0, width=1,
                           hovertemplate="Qscore: %{customdata[0]}~%{customdata[1]}<br>Count: %{y}",
                           name='', customdata=cus))
    fig = fig.update_layout(xaxis={'tickvals': x[:-1]})
    return fig











