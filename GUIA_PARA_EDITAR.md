# Guía rápida para personalizar el proyecto en GitHub

No necesitas volver a Google Colab. El proyecto se modifica directamente en los archivos de GitHub y Streamlit despliega los cambios tras el commit.

## ¿Dónde cambio cada cosa?

- **Título y descripción**: abre `app.py` y busca `<section class="hero">`.
- **Textos académicos**: abre `app.py` y busca `# ===== 01 =====`, `# ===== 02 =====` ... hasta `# ===== 06 =====`. Cada bloque pertenece a su pestaña.
- **Parámetros iniciales**: en `app.py`, busca `st.session_state.D = 8.5`, `st.session_state.k = 0.036` y `st.session_state.R = 400`. Ajusta también la función `restaurar()`.
- **Rango de parámetros**: en `app.py`, busca los `st.slider` de la barra lateral.
- **Colores y aspecto de la web**: edita las reglas de CSS de la parte superior de `app.py` y `.streamlit/config.toml`.
- **Gráficas**: en `graficos.py`, edita las funciones `mapa_concentracion`, `perfiles_tiempo`, `comparar_perfiles`, etc.
- **Ecuación física / solución numérica**: en `modelo.py`. Ten en cuenta que modificar los términos o las fronteras requiere volver a comprobar las pruebas.

## Cómo editar desde GitHub

1. Abre el archivo que quieras modificar.
2. Pulsa el icono del **lápiz** (Edit this file).
3. Modifica el texto o código.
4. Pulsa **Commit changes** y confirma.
5. Regresa a la página de Streamlit para ver los cambios publicados. Si aparecen errores, abre el apartado de logs de Streamlit Cloud.

## Recomendaciones para la entrega

- Completa los datos de los autores y el curso si decides añadirlos.
- Prueba los tres experimentos con los parámetros iniciales y con valores alternativos.
- Discute por qué la concentración **libre** puede ser baja en el centro aun cuando la superficie esté expuesta continuamente.
- No afirmes que la concentración libre sea igual a la concentración de medicamento ligada o internalizada.
- No presentes el experimento de imágenes de la pestaña 05 como un PET reconstruido: es una ilustración deliberadamente simplificada.
- Añade capturas y, si lo exige el profesor, exporta los resultados CSV de la pestaña 03 para documentar los experimentos.

## Qué pasa si el simulador deja de funcionar

- Comprueba que `app.py`, `modelo.py`, `graficos.py` y `requirements.txt` estén en la **misma carpeta raíz** del repositorio.
- No subas solo el archivo ZIP: Streamlit necesita los `.py` descomprimidos.
- En Streamlit Community Cloud, selecciona `app.py` como archivo principal.
- Si se rompe una fórmula tras un cambio, vuelve a la versión previa desde el historial de GitHub.
