"""Gráficas Plotly: independientes de Streamlit para facilitar la edición."""
from __future__ import annotations
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from modelo import perfil_en_2d

COLORS = ['#2E7DAB', '#0D9488', '#EAAD5C', '#D35B6B']


def formatear(fig, alto=360):
    fig.update_layout(
        font=dict(family='Arial, sans-serif', size=13, color='#223248'),
        margin=dict(l=15, r=15, t=42, b=22),
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)',
        height=alto,
        legend=dict(orientation='h', y=1.10, x=0),
        hoverlabel=dict(bgcolor='white', font_color='#17324B'),
    )
    fig.update_xaxes(showgrid=True, gridcolor='#E6EDF4', zeroline=False)
    fig.update_yaxes(showgrid=True, gridcolor='#E6EDF4', zeroline=False)
    return fig


def mapa_concentracion(res, idx, pix=111):
    x, z = perfil_en_2d(res.radios, res.valores[:, idx], res.radio_total, n_pix=pix)
    fig = go.Figure()
    fig.add_trace(go.Heatmap(
        x=x, y=x, z=z, zmin=0, zmax=1, colorscale='Viridis',
        colorbar=dict(title='C/Cₛ', len=0.8, thickness=12),
        hovertemplate='x=%{x:.0f} µm<br>y=%{y:.0f} µm<br>C/Cₛ=%{z:.3f}<extra></extra>',
    ))
    theta = np.linspace(0, 2*np.pi, 241)
    fig.add_trace(go.Scatter(
        x=res.radio_total*np.cos(theta), y=res.radio_total*np.sin(theta),
        mode='lines', line=dict(color='#8799A7', width=1.5, dash='dot'),
        showlegend=False, hoverinfo='skip'))
    fig.update_xaxes(title='x (µm)', scaleanchor='y', scaleratio=1, constrain='domain')
    fig.update_yaxes(title='y (µm)', constrain='domain')
    fig.update_layout(title=f'Sección 2D · t = {res.tiempos[idx]/60:.1f} min')
    return formatear(fig, 435)


def animacion_concentracion(res, fotogramas=26, pix=91):
    """Animación real: las imágenes proceden del modelo, no de ilustraciones."""
    indices = np.linspace(0, res.tiempos.size-1, fotogramas, dtype=int)
    x, z0 = perfil_en_2d(res.radios, res.valores[:, indices[0]], res.radio_total, n_pix=pix)
    fig = go.Figure(data=[go.Heatmap(
        z=z0, x=x, y=x, colorscale='Viridis', zmin=0, zmax=1,
        colorbar=dict(title='C/Cₛ', thickness=12),
        hovertemplate='x=%{x:.0f} µm<br>y=%{y:.0f} µm<br>C/Cₛ=%{z:.3f}<extra></extra>')])
    frames = []
    for i in indices:
        _, z = perfil_en_2d(res.radios, res.valores[:, i], res.radio_total, n_pix=pix)
        frames.append(go.Frame(data=[go.Heatmap(z=z, x=x, y=x, zmin=0, zmax=1, colorscale='Viridis')],
                               name=str(i)))
    fig.frames = frames
    fig.update_layout(
        title='Animación temporal de la penetración',
        xaxis=dict(title='x (µm)', scaleanchor='y', scaleratio=1),
        yaxis=dict(title='y (µm)'),
        updatemenus=[dict(type='buttons', showactive=False, y=-0.10, x=0,
          buttons=[dict(label='▶ Reproducir', method='animate',
            args=[None, {'frame': {'duration': 160, 'redraw': True},
                         'transition': {'duration': 0}, 'fromcurrent': True}]),
                   dict(label='⏸ Pausar', method='animate',
                    args=[[None], {'frame': {'duration': 0, 'redraw': False},
                                  'mode': 'immediate'}])])],
        sliders=[dict(steps=[dict(label=f'{res.tiempos[i]/60:.0f} min', method='animate',
               args=[[str(i)], {'mode': 'immediate',
                               'frame': {'duration': 0, 'redraw': True}}])
                              for i in indices],
            active=0, y=-0.23, x=0.05, len=0.9, pad=dict(t=15))]
    )
    return formatear(fig, 540)


def perfiles_tiempo(res, idx):
    fig = go.Figure()
    ti = sorted(set([0, max(0,idx//3), max(0,2*idx//3), idx]))
    for j, i in enumerate(ti):
        fig.add_trace(go.Scatter(
            x=np.r_[0, res.radios, res.radio_total],
            y=np.r_[res.valores[0,i], res.valores[:,i], 1.],
            mode='lines', line=dict(width=2.5, color=COLORS[j]),
            name=f'{res.tiempos[i]/60:.0f} min',
            hovertemplate='r=%{x:.1f} µm<br>C/Cₛ=%{y:.3f}<extra></extra>'))
    fig.update_layout(title='Perfil radial de concentración',
                      xaxis_title='Distancia al centro r (µm)',
                      yaxis_title='Concentración libre C/Cₛ')
    fig.update_yaxes(range=[0, 1.04])
    return formatear(fig)


def curvas_temporales(res):
    fig = go.Figure()
    for val, nombre, color in [
        (res.centro, 'Centro (aprox.)', COLORS[0]),
        (res.media, 'Media volumétrica', COLORS[1])]:
        fig.add_trace(go.Scatter(x=res.tiempos/60, y=val, name=nombre,
            line=dict(width=3, color=color),
            hovertemplate='t=%{x:.1f} min<br>C/Cₛ=%{y:.4f}<extra></extra>'))
    fig.update_layout(title='Evolución de la captación libre',
                      xaxis_title='Tiempo (min)', yaxis_title='C/Cₛ')
    fig.update_yaxes(range=[0, 1.04])
    return formatear(fig)


def comparar_perfiles(casos, tiempo_min):
    fig = go.Figure()
    for j, (nombre, res) in enumerate(casos):
        vals = res.valores[:, -1]
        fig.add_trace(go.Scatter(
          x=np.r_[0,res.radios/res.radio_total, 1],
          y=np.r_[vals[0], vals, 1], name=nombre, mode='lines',
          line=dict(width=3, color=COLORS[j % len(COLORS)])))
    fig.update_layout(title=f'Perfiles de penetración · t = {tiempo_min:g} min',
                      xaxis_title='Posición radial normalizada r/R',
                      yaxis_title='Concentración libre C/Cₛ')
    fig.update_yaxes(range=[0,1.05])
    return formatear(fig, 390)


def figura_imagenes(coords, real, medida, radio_roi):
    mask_real = real.copy()
    mask_med = medida.copy()
    fig=make_subplots(rows=1,cols=2, subplot_titles=['Distribución verdadera (2D)', 'Imagen sintética observada'],
                      horizontal_spacing=0.13)
    for col,z in [(1,mask_real),(2,mask_med)]:
        fig.add_trace(go.Heatmap(x=coords,y=coords,z=z,zmin=0,zmax=1,
            colorscale='Viridis', showscale=(col==2),
            colorbar=dict(title='Relativa',len=0.7,thickness=10),
            hovertemplate='x=%{x:.0f}<br>y=%{y:.0f}<br>señal=%{z:.3f}<extra></extra>'), row=1,col=col)
        a = np.linspace(0,2*np.pi,121)
        fig.add_trace(go.Scatter(x=radio_roi*np.cos(a),y=radio_roi*np.sin(a),mode='lines',
            line=dict(color='white', width=1.5, dash='dot'),showlegend=False,hoverinfo='skip'),row=1,col=col)
        fig.update_xaxes(scaleanchor=('y' if col==1 else 'y2'),scaleratio=1, title_text='x (µm)',row=1,col=col)
        fig.update_yaxes(title_text='y (µm)' if col==1 else '',row=1,col=col)
    fig.update_layout(title='Efectos de resolución y estadística de recuentos')
    return formatear(fig,435)


def histograma_montecarlo(muestras, referencia):
    fig=go.Figure(go.Histogram(x=muestras, nbinsx=24, marker_color=COLORS[0], opacity=0.85))
    fig.add_vline(x=referencia, line_width=2.5, line_color=COLORS[3],
                  annotation_text='Referencia real', annotation_position='top right')
    fig.update_layout(title='Distribución de estimaciones en 200 adquisiciones simuladas',
                      xaxis_title='Señal media en la ROI (unidades relativas)',
                      yaxis_title='Número de simulaciones')
    return formatear(fig,330)
