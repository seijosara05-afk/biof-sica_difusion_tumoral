"""Modelo biofísico radial de reacción–difusión para un esferoide.

Unidades: radio y coordenadas en µm, tiempo en s, D en µm²/s, k en s⁻¹.
C/C_s es adimensional. C describe exclusivamente fármaco *libre*.

Ecuación: ∂C/∂t = D/r² ∂r(r² ∂r C) - k C,
con C(R,t)=C_s, ∂r C(0,t)=0, C(r,0)=0 dentro de la esfera.

Se usa un método conservativo de volúmenes finitos (capas esféricas)
y una solución espectral exacta en tiempo del sistema semidiscreto.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy.linalg import eigh_tridiagonal, solve_banded


@dataclass
class Resultado:
    tiempos: np.ndarray  # s
    radios: np.ndarray  # µm, centros de capas
    valores: np.ndarray  # (N, n_t): C/C_s
    volumenes: np.ndarray  # µm³
    radio_total: float  # µm
    D: float  # µm²/s
    k: float  # s⁻¹
    cs: float  # valor normalizado de superficie = 1

    @property
    def media(self):
        return self.volumenes @ self.valores / self.volumenes.sum()

    @property
    def centro(self):
        # Aproximación mediante la primera capa de radio dr/2.
        return self.valores[0, :]

    @property
    def phi(self):
        return float(self.radio_total * np.sqrt(self.k / self.D))

    @property
    def longitud(self):
        return np.inf if self.k == 0 else float(np.sqrt(self.D / self.k))


def sistema_radial(D: float, k: float, R: float, n: int = 160):
    """Devuelve A tridiagonal, b y volúmenes para dc/dt=A c+b."""
    if not (D > 0 and k >= 0 and R > 0 and n >= 3):
        raise ValueError('Se requiere D>0, k>=0, R>0 y n>=3.')
    caras = np.linspace(0., R, n + 1)
    dr = R / n
    centros = (caras[:-1] + caras[1:]) / 2
    V = 4 * np.pi / 3 * (caras[1:] ** 3 - caras[:-1] ** 3)
    # Conductancia difusiva entre capas: superficie * D / separación.
    G = 4 * np.pi * D * caras[1:-1] ** 2 / dr
    # Última capa y baño exterior: distancia centro-superficie = dr/2.
    Gb = 4 * np.pi * D * R ** 2 / (dr / 2)
    diagonal = np.full(n, -k, dtype=float)
    diagonal[:-1] -= G / V[:-1]
    diagonal[1:] -= G / V[1:]
    diagonal[-1] -= Gb / V[-1]
    superior = G / V[:-1]
    inferior = G / V[1:]
    fuente = np.zeros(n)
    fuente[-1] = Gb / V[-1]  # C_s = 1
    return centros, V, diagonal, superior, inferior, fuente, G, Gb


def simular(D: float, k: float, R: float, tiempos: np.ndarray, n: int = 160) -> Resultado:
    """Resuelve la dinámica (sin restricción CFL gracias al propagador espectral)."""
    t = np.asarray(tiempos, dtype=float)
    if t.ndim != 1 or t.size == 0 or np.any(t < 0):
        raise ValueError('tiempos debe ser un vector no vacío de tiempos ≥0')
    radios, V, diag, sup, inf, b, _, _ = sistema_radial(D, k, R, n)

    # Sistema tridiagonal para el estacionario: A c_est = -b.
    ab = np.zeros((3, n))
    ab[0, 1:] = sup
    ab[1, :] = diag
    ab[2, :-1] = inf
    estacionario = solve_banded((1, 1), ab, -b)

    # Cambio z_i=sqrt(V_i)c_i. La matriz simétrica tiene diagonal 'diag'
    # y subdiagonal G_i/sqrt(V_i V_{i+1}).
    G = sup * V[:-1]
    off = G / np.sqrt(V[:-1] * V[1:])
    autovalores, autovectores = eigh_tridiagonal(diag, off)
    z_inicial = -np.sqrt(V) * estacionario   # c(0)=0.
    coeficientes = autovectores.T @ z_inicial
    modos = coeficientes[:, None] * np.exp(autovalores[:, None] * t[None, :])
    valores = estacionario[:, None] + (autovectores @ modos) / np.sqrt(V)[:, None]
    # Solo suprime error de redondeo (~1e-15) en situaciones de gran captación.
    valores = np.clip(valores, 0., 1.)
    valores[:, t == 0] = 0.
    return Resultado(t, radios, valores, V, R, D, k, 1.)


def perfil_estacionario_analitico(r: np.ndarray, D: float, k: float, R: float):
    """Solución continua estacionaria para concentración exterior unitaria."""
    r = np.asarray(r, dtype=float)
    if k == 0:
        return np.ones_like(r)
    L = np.sqrt(D / k)
    phi = R / L
    # C(r)/C_s = [sinh(r/L)/(r/L)]/[sinh(R/L)/(R/L)].
    x = r / L
    ratio = np.ones_like(x)
    valid = np.abs(x) > 1e-12
    ratio[valid] = np.sinh(x[valid]) / x[valid]
    return ratio * phi / np.sinh(phi)


def porcentaje_volumen_superior_umbral(radios, perfil, R, umbral=0.1):
    """Porcentaje del volumen de la esfera con C/C_s ≥ umbral.

    Se interpola el perfil radial y se supone que es monótono en r.
    """
    if not 0 < umbral < 1:
        raise ValueError('El umbral debe pertenecer a (0,1).')
    rs = np.r_[0.0, radios, R]
    cs = np.r_[perfil[0], perfil, 1.0]
    cs = np.maximum.accumulate(cs)
    if cs[0] >= umbral:
        return 100.0
    r_umbral = np.interp(umbral, cs, rs)
    return float(100. * (1. - (r_umbral / R) ** 3))


def perfil_en_2d(radios, perfil, R, n_pix=111, normalizado=False):
    """Convierte la solución esférica a una sección 2D para representar.

    Fuera de la esfera se fija 0 si normalizado=True (imagen sintética),
    o NaN para ocultar la región exterior en la figura del tejido.
    """
    coords = np.linspace(-R, R, n_pix)
    xx, yy = np.meshgrid(coords, coords)
    rr = np.sqrt(xx**2 + yy**2)
    z = np.interp(rr, np.r_[0, radios, R], np.r_[perfil[0], perfil, 1.0])
    z[rr > R] = 0.0 if normalizado else np.nan
    return coords, z


def diagnostico_balance(D: float, k: float, R: float, c: np.ndarray):
    """Prueba la identidad dM/dt=flujo_entrante-kM en el sistema discreto."""
    n = len(c)
    _, V, diag, sup, inf, b, _, Gb = sistema_radial(D, k, R, n)
    deriv = diag * c + b
    deriv[:-1] += sup * c[1:]
    deriv[1:] += inf * c[:-1]
    lhs = float(V @ deriv)
    rhs = float(Gb * (1. - c[-1]) - k * (V @ c))
    return lhs, rhs
