# Laboratorio virtual de Biofísica · Difusión tumoral

**Proyecto final de Biofísica (4.º de Física, USC)**

**Pregunta científica:** ¿Por qué los fármacos no siempre alcanzan el interior de los tumores?

Aplicación científica interactiva desarrollada con Streamlit que resuelve un modelo de **reacción–difusión esférica** en un esferoide tumoral idealizado. Incluye un laboratorio de parámetros, tres experimentos comparativos, análisis cuantitativo, una **analogía didáctica** de imagen molecular y bibliografía.

> **Aviso científico:** No se trata de una simulación terapéutica ni de adquisición PET real. El término de captación describe la retirada de *fármaco libre*, no su acumulación intracelular. Los parámetros son exploratorios.

## Qué contiene

| Archivo | Función | Qué puedes editar |
|---|---|---|
| `app.py` | Todas las pantallas y controles del proyecto | Títulos, textos, parámetros iniciales, experimentos, diseño |
| `modelo.py` | Ecuaciones físicas, solución numérica y métricas | Modelo, condiciones de contorno, método |
| `graficos.py` | Gráficas y animación Plotly | Colores, títulos, leyendas, estilos |
| `requirements.txt` | Bibliotecas Python para Streamlit Cloud | Dependencias |
| `.streamlit/config.toml` | Aspecto visual y tema | Colores generales |
| `tests/test_modelo.py` | Comprobaciones de coherencia matemática | Pruebas físicas |
| `README.md` | Documentación y despliegue | Presentación del repositorio |

## Subir a GitHub (desde el navegador)

1. Descarga y descomprime el ZIP entregado.
2. En el repositorio que acabas de crear, haz clic en **Add file → Upload files**.
3. Arrastra **los archivos y carpetas que están dentro de la carpeta del proyecto**, no el ZIP ni la carpeta contenedora.
4. Comprueba que `app.py` y `requirements.txt` estén en la **raíz** del repositorio y pulsa **Commit changes**.
5. En [Streamlit Community Cloud](https://share.streamlit.io/), selecciona **Create app**, tu repositorio, la rama `main` y **Main file path: `app.py`**.
6. La aplicación quedará asociada al repositorio. Cada cambio que guardes con **Commit changes** en GitHub podrá desplegarse automáticamente.

Para editar código sin instalar nada, abre un archivo en GitHub y pulsa el lápiz, o utiliza GitHub Codespaces.

## Ejecución en un ordenador propio (opcional)

```bash
python -m pip install -r requirements.txt
python -m streamlit run app.py
```

### Pruebas automáticas

```bash
python -m pip install pytest
python -m pytest -q
```

## Matemática del modelo

Se supone simetría radial esférica, con radio `R`, coeficiente de difusión `D` y retirada de primer orden `k`:

\[
\partial_t C= D\,r^{-2}\,\partial_r(r^2\partial_r C)-k C.
\]

Condiciones: interior inicialmente vacío, concentración superficial fija `C(R,t)=C_s`, y simetría en el centro `∂C/∂r (0,t)=0`.

La concentración presentada está **normalizada**: `C/C_s`, sin unidades. El radio se expresa en `µm`, el tiempo en `s` y los parámetros en `µm²/s` y `s⁻¹`.

La discretización utiliza volúmenes finitos esféricos y un propagador espectral de un sistema tridiagonal lineal. La forma conservativa respeta el balance de cantidad libre en la esfera:

\[
\frac{dM}{dt}=J_{\mathrm{entrada}}-kM,\qquad M=\int_V C\,dV.
\]

Magnitudes: longitud característica `ell = sqrt(D/k)` (si k>0), módulo `phi = R/ell`, concentración central aproximada, media en volumen y porcentaje del volumen por encima del umbral `0.1 C_s`.

Para `k>0`, el régimen estacionario del modelo **continuo** permite verificar los resultados numéricos mediante:

\[
\frac{C_\infty(r)}{C_s}=\frac{R}{r}\frac{\sinh(r/\ell)}{\sinh(R/\ell)}.
\]

### Importante sobre el PET

La pestaña 05 **no** calcula SUV, ni modela fotones, detectores, correcciones o reconstrucción. Trata una **sección 2D sintética** como una distribución de actividad *hipotética* proporcional a `C/C_s`; simula desenfoque con una PSF gaussiana y recuentos de Poisson simplificados **a nivel píxel**. Esto permite practicar los conceptos de sesgo, dispersión e incertidumbre de medición, sin atribuir a la simulación validez clínica.

## Bibliografía

- Ruso, J. M. *Modelos matemáticos*, material de Biofísica (USC), pp. 34–36 (reacción–difusión).
- Ruso, J. M. *Dimensiones y constituyentes de los sistemas biológicos*, material de Biofísica (USC).
- Minchinton, A. I. y Tannock, I. F. (2006). *Drug penetration in solid tumours*. Nature Reviews Cancer 6, 583–592. https://doi.org/10.1038/nrc1893
- Kim, M., Gillies, R. J. y Rejniak, K. A. (2013). *Current Advances in Mathematical Modeling of Anti-Cancer Drug Penetration into Tumor Tissues*. Frontiers in Oncology 3, 278. https://doi.org/10.3389/fonc.2013.00278

## Aspectos pendientes para la entrega académica

Antes de entregar, incorpora en la app la autoría y las fechas correctas, confirma que las suposiciones físicas están claras, comprueba los tres experimentos con varios conjuntos de parámetros y, si el profesor lo exige, adjunta una breve discusión personal sobre los resultados. La aplicación está diseñada para modificarla y convertirla en vuestro trabajo.
