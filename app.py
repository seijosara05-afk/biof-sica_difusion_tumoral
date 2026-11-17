"""Aplicación web editable de Biofísica · Universidad de Santiago de Compostela.

Estructura:
  01 problema -> 02 modelo -> 03 laboratorio -> 04 experimentos
  -> 05 analogía de imagen molecular -> 06 discusión, referencias.

Modelos físicos: modelo.py. Gráficas: graficos.py.
"""
from __future__ import annotations
import io
import numpy as np
import pandas as pd
import streamlit as st
from scipy.ndimage import gaussian_filter
from modelo import (simular, perfil_estacionario_analitico,
                    porcentaje_volumen_superior_umbral, perfil_en_2d)
from graficos import (mapa_concentracion, animacion_concentracion,
                      perfiles_tiempo, curvas_temporales, comparar_perfiles,
                      figura_imagenes, histograma_montecarlo)

st.set_page_config(page_title='Biofísica · Difusión tumoral', page_icon='🧫',
                   layout='wide', initial_sidebar_state='expanded')

# ---------- Estilo general (editable) ----------
st.markdown('''
<style>
  .block-container {padding-top:1.35rem; padding-bottom:3rem; max-width:1380px}
  h1,h2,h3 {letter-spacing:-.025em}
  h1 {font-size:2.25rem !important; line-height:1.17 !important}
  [data-testid="stSidebar"] {background-color:#F2F6FA}
  .hero {background:linear-gradient(105deg,#0E2940,#126878);
         border-radius:18px; padding:30px 34px; color:white; margin-bottom:20px}
  .hero h1 {color:white !important;margin:0.25rem 0 .8rem 0 !important}
  .hero p {color:#CEE8EE;font-size:1.04rem;margin-bottom:0}
  .eyebrow {letter-spacing:0.12em;text-transform:uppercase;color:#8BDFD6;
            font-weight:800;font-size:.79rem}
  .mini {color:#577186;font-size:.9rem}
  div[data-testid="stMetric"] {background:#F4F8FA; border:1px solid #DCE6EC;
      padding:16px; border-radius:14px}
  div[data-testid="stMetricLabel"]{color:#486275}
  .stTabs [data-baseweb="tab"] {font-weight:700; padding-left:10px; padding-right:10px}
  a {color:#087A87}
</style>
''', unsafe_allow_html=True)

# ---------- Estado y parámetros ----------
if 'D' not in st.session_state: st.session_state.D = 8.5
if 'k' not in st.session_state: st.session_state.k = 0.036
if 'R' not in st.session_state: st.session_state.R = 400
if 'horizonte' not in st.session_state: st.session_state.horizonte = 120


def restaurar():
    st.session_state.D = 8.5
    st.session_state.k = 0.036
    st.session_state.R = 400
    st.session_state.horizonte = 120

with st.sidebar:
    st.markdown('### ⚙️ Parámetros del tejido')
    st.caption('Cambia los valores y observa cómo responde el mismo modelo biofísico.')
    D = st.slider('Difusión D (µm²/s)', 0.5, 30.0, step=0.5, key='D',
                  help='Movilidad efectiva de la sustancia en el tejido.')
    
    k = st.slider(
        'Captación k (s⁻¹)',
        min_value=0.0,
        max_value=0.010,
        value=0.001,
        step=0.0001,
        key='captacion_k',
        format='%.4f',
        help=(
            "Constante que representa la velocidad de retirada del fármaco libre "
            "debido a procesos como la captación celular o la degradación. "
            "Cuanto mayor es k, más rápidamente disminuye la concentración libre "
            "y menor es su penetración característica en el tejido. "
            "Si k = 0, únicamente actúa la difusión. "
            "Unidad: s⁻¹."))


    R = st.slider('Radio R (µm)', 100, 600, step=25, key='R',
                  help='Radio del esferoide esférico idealizado.')
    horizonte = st.slider('Duración total (min)', 15, 480, step=15,
                          key='horizonte')
    st.button('↺ Restaurar valores iniciales', on_click=restaurar, use_container_width=True)
    st.divider()
    ell = np.inf if k==0 else np.sqrt(D/k)
    phi = 0 if k==0 else R/ell
    st.markdown('**Escalas físicas**')
    st.write(f'Longitud característica: **{"∞ - no existe una escala finita de penetración asociada a la captación (k = 0)" if k==0 else f"{ell:.1f} µm"}**')
    st.write(f'Parámetro adimensional: **φ = {phi:.2f}**')
    
@st.cache_data(show_spinner=False, max_entries=28)
def calcular(D, k, R, duracion_min, n=160):
    return simular(D, k, R, np.linspace(0, duracion_min*60, 101), n=n)

res = calcular(D, k, R, horizonte)

st.markdown('''<section class="hero">
<div class="eyebrow">BIOFÍSICA · PROYECTO FINAL · LABORATORIO VIRTUAL</div>
<h1>¿Por qué los fármacos no siempre alcanzan el interior de los tumores?</h1>
<p>Simulación interactiva del transporte molecular en esferoides tumorales ·<br>
Difusión, captación, geometría y cuantificación de distribuciones.</p>
</section>''', unsafe_allow_html=True)

tab1,tab2,tab3,tab4,tab5,tab6 = st.tabs([
    '01 · El problema', '02 · Modelo biofísico', '03 · Laboratorio',
    '04 · Experimentos', '05 · Imagen molecular', '06 · Conclusiones'
])

# ===== 01 =====
with tab1:
    st.header('01 — El problema')
    st.subheader('¿Llegar al tumor es lo mismo que llegar a todas sus células?')
    st.markdown('''
    Para que un fármaco actúe sobre una célula, primero tiene que conseguir llegar hasta ella.
    En los tumores sólidos, la presencia de una sustancia en el entorno del tumor
    no garantiza que alcance por igual todas sus regiones. Su distribución depende
    de las distancias que debe recorrer, de la estructura del tejido y de las
    interacciones que experimenta durante el trayecto.

    En un tumor real intervienen numerosos factores en el proceso como la vascularización,
    la permeabilidad de los vasos sanguíneos, la estructura de la matriz extracelular,...
    Nosotros **proponemos un modelo simplificado que considera únicamente dos mecanismos
    fundamentales: la difusión molecular y la retirada de sustancia libre**. Esta 
    aproximación permite analizar de forma aislada cómo la competencia entre ambos procesos
    condiciona la distribución espacial del fármaco y su capacidad para alcanzar las
    regiones más profundas del tejido, sin introducir toda la complejidad de un sistema
    tumoral real.


    Los **esferoides tumorales multicelulares** son conjuntos tridimensionales de
    células utilizados como modelos experimentales simplificados de tejido tumoral.
    En este trabajo representamos uno de ellos mediante una **esfera de radio $R$**:
    la sustancia se suministra desde el medio exterior, atraviesa la superficie
    y comienza a difundirse hacia el centro.
    
    Vamos a suponer que el exterior del esferoide se mantiene a una concentración fija, 
    mientras que en su interior se difunde una sustancia y parte de
    ella deja de estar libre. La cuestión no es solo cuánto entra, sino hasta dónde llega.
    ''')
    c1,c2,c3=st.columns(3)
    c1.info('**Difusión**\n\nLa agitación molecular permite el transporte neto desde zonas de mayor a menor concentración. Un $D$ más alto favorece la redistribución espacial.')
    c2.info('**Captación / eliminación**\n\nUna parte de las moléculas deja de estar libre para seguir difundiendo. El parámetro $k$ regula esa retirada.')
    c3.info('**Tamaño**\n\nUn esferoide mayor obliga a recorrer distancias mayores. Por eso el radio $R$ también afecta a la penetración.')
    st.markdown(r"""
    ### Los parámetros del simulador

    Los controles de la izquierda permiten cambiar **un mecanismo físico cada vez**
    y comprobar cómo cambia la distribución. Sus valores son **exploratorios**:
    los hemos elegido para visualizar distintos regímenes de nuestro modelo.
    """)
    st.table(pd.DataFrame([
        {'Parámetro':'Difusión D',
         'Unidad':'µm²/s', 'Intervalo':'0,5–30'},
        {'Parámetro':'Captación k',
         'Unidad':'s⁻¹', 'Intervalo':'0–0,060'},
        {'Parámetro':'Radio R',
         'Unidad':'µm', 'Intervalo':'100–600'},
        {'Parámetro':'Duración',
         'Unidad':'min', 'Intervalo':'15–480'},
    ]))
    st.markdown(r"""
    **¿Por qué esos intervalos?** Son rangos de exploración que permiten comparar
    procesos lentos y rápidos dentro del *mismo modelo*.""")

    st.markdown(r"""
    La difusión y la captación compiten en **escalas temporales distintas**.
    Para compararlas, utilizamos una escala difusiva $t_D=R^2/D$, el tiempo
    característico de retirada $t_k=1/k$ (cuando $k>0$) y la longitud
    $\ell=\sqrt{D/k}$. Esta última nos indica cuánto puede penetrar
    característicamente la concentración libre antes de disminuir de forma
    apreciable por captación.
    """)
    t_d = R**2 / D
    t_k = float('inf') if k == 0 else 1 / k
    escala1,escala2,escala3=st.columns(3)
    escala1.metric('Escala difusiva R²/D',f'{t_d/60:.1f} min',
                   help='Escala dimensional; no equivale al tiempo exacto de estabilización.')
    escala2.metric('Escala de captación 1/k',
                   'Sin captación' if k == 0 else f'{t_k:.1f} s',
                   help='Tiempo para caer a 1/e sin aporte ni difusión.')
    escala3.metric('Longitud √(D/k)',
                   'No finita' if k == 0 else f'{ell:.1f} µm',
                   help='Escala de penetración asociada a difusión y retirada.')
    if k == 0:
        st.success('**Con estos valores:** no hay retirada de sustancia libre. Si mantenemos el suministro exterior el tiempo suficiente, la concentración tenderá a uniformizarse en todo el esferoide.')
    elif phi < 1:
        st.info(f'**Con estos valores, φ = {phi:.2f}:** la longitud característica de penetración ({ell:.1f} µm) es mayor que el radio ({R} µm). En el régimen estacionario se espera una concentración relativamente homogénea.')
    elif phi < 3:
        st.info(f'**Con estos valores, φ = {phi:.2f}:** difusión y retirada compiten a escalas comparables. Se espera un gradiente de concentración apreciable entre la superficie y el centro.')
    else:
        st.warning(f'**Con estos valores, φ = {phi:.2f}:** la longitud característica ({ell:.1f} µm) es mucho menor que el radio ({R} µm). Se prevé una concentración libre central muy reducida respecto a la superficie. Prueba a disminuir k o R para comparar.')


    st.markdown(r"""
    ### La pregunta que investigaremos

    **¿Puede una captación celular más intensa producir una concentración
    *libre* menor en el centro, aunque el suministro exterior sea el mismo?**
    El laboratorio permitirá poner a prueba esta hipótesis variando $D$, $k$
    y $R$ por separado, y observando perfiles radiales, mapas de concentración
    y la fracción de volumen que supera un umbral relativo.

    Conviene distinguir **llegada de sustancia libre** y **efecto terapéutico**:
    el modelo no calcula cuánta sustancia se acumula dentro de las células,
    su toxicidad ni la muerte celular. Tampoco reproduce todas las barreras
    presentes en un tumor real, como mencionamos anteriormente. Su finalidad es explicar, mediante un sistema
    controlado, cómo la competición entre difusión y captación puede generar
    distribuciones espaciales heterogéneas.
    """)
   
# ===== 02 =====
with tab2:
    st.header('02 — El modelo biofísico')
    st.markdown(r'''
    La **segunda ley de Fick** describe la evolución de la concentración libre
    $C$ cuando el transporte se debe a difusión:
    ''')
    st.latex(r'\frac{\partial C}{\partial t}=D\nabla^2 C')
    st.markdown(r'''
    Para incluir una retirada de primer orden, introducimos el término $-kC$.
    Suponiendo **simetría esférica** (la concentración depende solo del radio $r$),
    el problema queda:
    ''')
    st.latex(r'\boxed{\frac{\partial C}{\partial t}=\frac{D}{r^2}\frac{\partial}{\partial r}\left(r^2\frac{\partial C}{\partial r}\right)-kC}')
    st.markdown(r'''
    La variable $C$ representa únicamente **sustancia libre**; el modelo no distingue
    moléculas ligadas, internalizadas y degradadas, que quedan agrupadas en un único
    proceso efectivo de retirada.
    ''')
    st.markdown('**Hipótesis y condiciones de contorno**')
    st.markdown(r'''
    - Esferoide homogéneo e isótropo de radio $R$; $D$ y $k$ son constantes.
    - No existe flujo convectivo ni un sistema vascular interno.
    - El interior está inicialmente vacío: $C(r,0)=0$, para $0\le r<R$.
    - La superficie está en contacto con un reservorio mantenido a $C(R,t)=C_s$.
    - En el centro se cumple la condición de simetría $\partial_r C(0,t)=0$.
    - Representamos la concentración **normalizada** $C/C_s$; no trabajamos con una dosis real.
    ''')
    st.markdown('**Competencia entre escalas**')
    st.latex(r'\ell = \sqrt{D/k},\qquad \phi=\frac{R}{\ell}=R\sqrt{k/D}\quad(k>0)')
    st.markdown(r'''
    La longitud $\ell$ es la escala de penetración asociada a la competencia entre
    difusión y retirada. El parámetro adimensional $\phi$ compara el tamaño del
    esferoide con esa longitud. Si $\phi$ es elevado, la concentración puede caer
    de forma muy marcada hacia el centro. Si $k=0$, no existe una longitud de
    eliminación finita y, a tiempo suficientemente largo, la concentración se uniformiza.
    ''')
    with st.expander('Solución analítica en régimen estacionario (para contrastar el cálculo numérico)'):
        st.markdown(r'''
        Al imponer $\partial C/\partial t=0$ y mantener $C(R)=C_s$, se obtiene:
        ''')
        st.latex(r'\frac{C_\infty(r)}{C_s}=\frac{R}{r}\,\frac{\sinh(r/\ell)}{\sinh(R/\ell)}\qquad (r>0)')
        st.markdown(r'''
        con límite central $C_\infty(0)/C_s=\phi/\sinh\phi$.
        Esta solución corresponde al modelo esférico continuo, no a datos clínicos.
        ''')
    with st.expander('Método numérico y verificación física'):
        st.markdown(r'''
        Se discretiza la esfera en **160 capas de volumen finito**, respetando el
        volumen de cada capa y sus superficies de intercambio. El flujo entre capas
        se calcula mediante la primera ley de Fick. Esto reproduce el balance:
        ''')
        st.latex(r'\frac{d}{dt}\int_V C\,dV=J_{\mathrm{entrada}}-k\int_V C\,dV')
        st.markdown(r'''
        El sistema lineal resultante se resuelve con una descomposición espectral,
        exacta en tiempo para la discretización elegida. Esto evita la inestabilidad
        temporal de un esquema explícito, aunque persiste el error espacial de la malla.
        La solución se contrasta con el estacionario analítico y se incluyen pruebas
        automáticas en `tests/test_modelo.py`.
        ''')
        anal = perfil_estacionario_analitico(res.radios, D, k, R)
        num = simular(D,k,R,np.array([1e8]),n=160).valores[:,-1]
        err = np.max(np.abs(num-anal))
        st.metric('Error máximo estacionario (malla vs. teoría)',f'{err:.3g} (C/Cₛ)')
        st.caption('Para esta validación se calcula un tiempo suficientemente largo (independiente del tiempo mostrado en el laboratorio), aislando el error de discretización espacial.')
    st.caption('Conexión con los apuntes: ecuaciones de reacción–difusión de la asignatura de Biofísica, presentación «Modelos matemáticos», pp. 34–36.')

# ===== 03 =====
with tab3:
    st.header('03 — Laboratorio interactivo')
    st.markdown('Modifica **D, k, R y duración** en la barra lateral. Explora cómo cambia el perfil de concentración libre en el interior del esferoide.')
    tiempo_min = st.slider('Instante que queremos observar (min)', 0, horizonte,
                           min(60, horizonte), 1, key='tiempo_mostrado')
    idx = int(np.argmin(np.abs(res.tiempos/60-tiempo_min)))
    mapa, perfil=st.columns([1.06,0.94], gap='large')
    with mapa:
        st.plotly_chart(mapa_concentracion(res,idx), use_container_width=True)
    with perfil:
        st.plotly_chart(perfiles_tiempo(res,idx), use_container_width=True)

    cobertura=porcentaje_volumen_superior_umbral(res.radios,res.valores[:,idx],R,0.1)
    c1,c2,c3,c4=st.columns(4)
    c1.metric('Centro (C/Cₛ)',f'{100*res.centro[idx]:.3g} %',help='Primera capa radial, aprox. r=0')
    c2.metric('Media en volumen',f'{100*res.media[idx]:.2f} %')
    c3.metric('Volumen con C/Cₛ ≥ 0,10',f'{cobertura:.1f} %')
    c4.metric('φ = R/ℓ', f'{res.phi:.2f}')
    st.caption('Los porcentajes se refieren a la concentración superficial normalizada, no a células afectadas ni a eficacia terapéutica.')
    if k > 0 and (R/160) > 0.5*np.sqrt(D/k):
        st.warning('La longitud de penetración es pequeña respecto a la malla de 160 capas: en estos parámetros extremos, conviene refinar la discretización para obtener resultados cuantitativos precisos.')
    st.plotly_chart(curvas_temporales(res), use_container_width=True)
    with st.expander('▶ Ver animación de la penetración (controles de reproducción)'):
        st.plotly_chart(animacion_concentracion(res),use_container_width=True)
        st.caption('Los fotogramas se calculan con el mismo modelo que las demás gráficas. La reproducción se controla dentro de la figura.')

    csv_radial = pd.DataFrame({'radio_um':res.radios,'concentracion_relativa':res.valores[:,idx]}).to_csv(index=False).encode('utf-8-sig')
    csv_temporal = pd.DataFrame({'tiempo_min':res.tiempos/60,
        'concentracion_centro_relativa':res.centro,
        'concentracion_media_volumen_relativa':res.media}).to_csv(index=False).encode('utf-8-sig')
    d1,d2=st.columns(2)
    d1.download_button('⬇ Descargar perfil radial (CSV)',csv_radial,'perfil_radial.csv','text/csv',use_container_width=True)
    d2.download_button('⬇ Descargar evolución temporal (CSV)',csv_temporal,'evolucion_temporal.csv','text/csv',use_container_width=True)

# ===== 04 =====
with tab4:
    st.header('04 — ¿Qué hemos descubierto?')
    st.markdown('Tres experimentos numéricos: **en cada uno variamos solo un parámetro**, manteniendo los demás y el instante final constantes. El botón del laboratorio permite explorar más casos.')
    tf = st.slider('Tiempo de comparación (min)', 15, 480, 120, 15, key='tiempo_comparacion')
    def mostrar_casos(casos_param):
        lista=[]; filas=[]
        for nombre,dd,kk,rr in casos_param:
            salida = simular(dd,kk,rr,np.array([0.,float(tf*60)]),n=160)
            lista.append((nombre,salida))
            filas.append({
                'Escenario':nombre, 'D (µm²/s)':round(dd,3), 'k (s⁻¹)':round(kk,4),
                'R (µm)':int(rr), 'φ':round(salida.phi,2),
                'C central (% Cₛ)':round(100*salida.centro[-1],4),
                'C media (% Cₛ)':round(100*salida.media[-1],2),
                'Volumen ≥10% (%)':round(porcentaje_volumen_superior_umbral(
                    salida.radios,salida.valores[:,-1],rr),1)
            })
        st.plotly_chart(comparar_perfiles(lista,tf),use_container_width=True)
        st.dataframe(pd.DataFrame(filas),hide_index=True,use_container_width=True)
        return filas

    a,b,c=st.tabs(['A · Cambiar la difusión','B · Cambiar la captación','C · Cambiar el radio'])
    with a:
        st.markdown('**Hipótesis.** A igualdad de geometría y captación, aumentar $D$ favorece la penetración de moléculas libres.')
        aa=mostrar_casos([('D/4',D/4,k,R),('D (referencia)',D,k,R),('4D',4*D,k,R)])
        st.markdown(f'**Resultado numérico:** al comparar D/4 y 4D, la fracción de volumen por encima del 10 % superficial cambia de **{aa[0]["Volumen ≥10% (%)"]:.1f} %** a **{aa[2]["Volumen ≥10% (%)"]:.1f} %**.')
        st.caption('Comparación a tiempo fijo. La magnitud del efecto depende de si el régimen es transitorio o estacionario.')
    with b:
        st.markdown('**Hipótesis.** Una retirada rápida de moléculas libres reduce su capacidad de alcanzar las regiones interiores.')
        bb=mostrar_casos([('Sin retirada (k=0)',D,0.,R),('k (referencia)',D,k,R),('4k',D,4*k,R)])
        st.markdown(f'**Resultado numérico:** la concentración media libre es **{bb[0]["C media (% Cₛ)"]:.2f} %** sin captación, frente a **{bb[2]["C media (% Cₛ)"]:.2f} %** para 4k.')
        st.caption('El término -kC representa una desaparición efectiva de sustancia libre; no permite deducir qué fracción queda unida a las células.')
    with c:
        st.markdown('**Hipótesis.** Aumentar la distancia hasta el centro dificulta su abastecimiento en el mismo tiempo de exposición.')
        cc=mostrar_casos([('R/2',D,k,R/2),('R (referencia)',D,k,R),('1,5R',D,k,R*1.5)])
        st.markdown(f'**Resultado numérico:** al pasar de R/2 a 1,5R, la concentración media libre cambia de **{cc[0]["C media (% Cₛ)"]:.2f} %** a **{cc[2]["C media (% Cₛ)"]:.2f} %**.')
        st.caption('Se utiliza r/R para superponer perfiles de esferoides con distintos radios.')

    st.markdown('### Lectura física de los resultados')
    st.markdown(r'''
    Las tres comparaciones muestran que el reparto espacial depende tanto de la **velocidad
    difusiva** como de la **cinética de retirada** y de la **longitud que debe recorrerse**.
    El grupo $\phi=R\sqrt{k/D}$ reúne esas tres dependencias en una sola magnitud.
    Conviene evitar la interpretación simplista «más captación = más fármaco en todo el tumor»:
    el modelo representa la **fracción libre**, que puede disminuir en el interior al
    elevar $k$. Tampoco podemos concluir que esa captación implique una mayor eficacia.
    ''')

# ===== 05 =====
with tab5:
    st.header('05 — Aplicación a imagen molecular')
    st.markdown('''
    **Extensión conceptual inspirada en PET preclínica.** Partimos de una sección 2D
    del campo calculado en el laboratorio y suponemos, *únicamente para este experimento*,
    que su valor es proporcional a una actividad radiactiva hipotética. Simulamos una
    imagen de baja resolución mediante desenfoque gaussiano y ruido de Poisson sobre
    los píxeles. **Esto no es un simulador de adquisición ni de reconstrucción PET.**
    ''')
    x1,x2,x3,x4=st.columns(4)
    with x1:
        fwhm = st.slider('Desenfoque FWHM (µm)', 10, 250, 85, 5)
    with x2:
        ganancia = st.slider('Recuentos por unidad de señal', 10, 1000, 180, 10)
    with x3:
        frac_roi=st.slider('Radio de la ROI / R', 0.40, 1.0, 0.90, 0.05)
    with x4:
        semilla=st.number_input('Semilla aleatoria',min_value=0,max_value=100000,value=7)
    npix=111
    coords,verdad=perfil_en_2d(res.radios,res.valores[:, -1],R,n_pix=npix,normalizado=True)
    dx=coords[1]-coords[0]
    sigma_pix=fwhm/(2*np.sqrt(2*np.log(2))*dx)
    desenf=gaussian_filter(verdad,sigma=sigma_pix,mode='constant',cval=0)
    rng=np.random.default_rng(int(semilla))
    lam=ganancia*desenf
    observacion=rng.poisson(lam)/ganancia
    xx,yy=np.meshgrid(coords,coords)
    roi=(xx**2+yy**2)<=(frac_roi*R)**2
    referencia=float(np.mean(verdad[roi]))
    estimacion=float(np.mean(observacion[roi]))
    m1,m2=st.columns([1.7,1])
    with m1:
        st.plotly_chart(figura_imagenes(coords,verdad,observacion,frac_roi*R),use_container_width=True)
    with m2:
        st.markdown('#### Cuantificación en la ROI')
        st.metric('Señal media verdadera',f'{referencia:.4f}')
        st.metric('Señal media estimada',f'{estimacion:.4f}')
        st.metric('Diferencia absoluta',f'{estimacion-referencia:+.4f}')
        st.caption('Unidades de actividad relativas; no se ha calculado SUV ni un biomarcador clínico.')
    # Repeticiones independientes de recuentos de Poisson para medir precisión.
    muestra = rng.poisson(lam, size=(200,npix,npix))/ganancia
    estimaciones=muestra[:,roi].mean(axis=1)
    sesgo=float(estimaciones.mean()-referencia)
    dispersion=float(estimaciones.std(ddof=1))
    st.plotly_chart(histograma_montecarlo(estimaciones,referencia),use_container_width=True)
    p1,p2,p3=st.columns(3)
    p1.metric('Sesgo estimado (media de 200 − verdad)',f'{sesgo:+.5f}')
    p2.metric('Desviación típica entre réplicas',f'{dispersion:.5f}')
    p3.metric('Número de realizaciones', '200')
    st.markdown(r'''
    **Interpretación.** Cambiar la resolución espacial puede alterar la señal media
    dentro de una región de interés. Cambiar los recuentos modifica la dispersión de
    las estimaciones. El sesgo compara la media de las medidas simuladas con una
    referencia conocida; la desviación típica cuantifica su variabilidad.
    ''')
    st.info('En PET real, la incertidumbre depende de los recuentos detectados, correcciones, reconstrucción, volumen parcial, delimitación de regiones y otros factores. Aquí solo aislamos dos mecanismos mediante un modelo didáctico 2D.')

# ===== 06 =====
with tab6:
    st.header('06 — Conclusiones y referencias')
    st.markdown('### Conclusiones de la simulación')
    st.markdown(f'''
    Con los parámetros seleccionados (**D = {D:g} µm²/s, k = {k:.3f} s⁻¹,
    R = {R:g} µm**), la razón adimensional es **φ = {res.phi:.2f}**.
    En el instante final de **{horizonte} min**, la concentración libre media del esferoide
    equivale al **{100*res.media[-1]:.2f} %** de la concentración superficial, mientras que
    la primera capa central alcanza aproximadamente el **{100*res.centro[-1]:.3g} %**.

    Estos resultados son predicciones de un modelo idealizado. La difusión y la
    retirada de sustancia libre compiten; la geometría esférica introduce distancias
    de transporte variables y hace relevante la relación entre tamaño y longitud
    característica. Comparar escenarios con una sola variable modificada ayuda a
    identificar mecanismos, pero no proporciona parámetros farmacológicos reales.
    ''')
    st.markdown('### Limitaciones')
    st.markdown('''
    - No hay vasos, convección, necrosis, heterogeneidad, receptores saturables ni unión reversible explícita.
    - El baño externo tiene una concentración constante: no se modela una dosis intravenosa real.
    - La sustancia libre retirada no se sigue en un compartimento celular; no se modela la eficacia ni la toxicidad.
    - Las figuras 2D son cortes de un modelo radial 3D; la pestaña de imagen molecular no reconstruye datos PET.
    - Los parámetros son exploratorios; sería necesario contrastarlos con experimentos para estudiar un compuesto concreto.
    ''')
    st.markdown('### Bibliografía y relación con la asignatura')
    st.markdown('''
    1. **Ruso, J. M.** *Modelos matemáticos*, material de Biofísica de la USC,
       especialmente las diapositivas de ecuaciones de reacción–difusión (pp. 34–36).
    2. **Minchinton, A. I. y Tannock, I. F. (2006).** Drug penetration in solid tumours.
       *Nature Reviews Cancer*, 6, 583–592.
       [DOI: 10.1038/nrc1893](https://doi.org/10.1038/nrc1893).
    3. **Kim, M., Gillies, R. J. y Rejniak, K. A. (2013).** Current Advances in Mathematical
       Modeling of Anti-Cancer Drug Penetration into Tumor Tissues.
       *Frontiers in Oncology*, 3, 278.
       [DOI: 10.3389/fonc.2013.00278](https://doi.org/10.3389/fonc.2013.00278).
    4. **Ruso, J. M.** *Dimensiones y constituyentes de los sistemas biológicos*,
       material de Biofísica de la USC (importancia de la escala térmica $k_B T$).
    ''')
    st.markdown('### Reproducibilidad y edición')
    st.markdown('''
    El repositorio incluye `modelo.py` (ecuaciones y cálculo), `graficos.py`
    (visualizaciones), `app.py` (interfaz y contenido), `tests/test_modelo.py`
    (pruebas automáticas) y `README.md` (instrucciones).
    Todos los gráficos y resultados se obtienen de las ecuaciones implementadas,
    sin introducir datos experimentales inventados.
    ''')

st.divider()
st.caption('Biofísica · Simulación educativa · Concentraciones normalizadas · Python / NumPy / SciPy / Plotly / Streamlit')
