"""Ejecutar: python -m pytest -q (instalar pytest si es necesario)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
import pytest
from modelo import (
    simular, sistema_radial, perfil_estacionario_analitico,
    porcentaje_volumen_superior_umbral, diagnostico_balance, perfil_en_2d)


def test_inicio_y_valores_fisicos():
    r = simular(8.5, 0.036, 400, np.array([0, 10., 100., 7200.]), n=160)
    assert np.array_equal(r.valores[:, 0], np.zeros(160))
    assert np.all(np.isfinite(r.valores))
    assert np.all((r.valores >= 0) & (r.valores <= 1))
    assert np.all(np.diff(r.media) >= -1e-9)
    assert r.phi == pytest.approx(400 * np.sqrt(0.036/8.5))


def test_difusion_sin_captacion_tiende_a_uniforme():
    r = simular(8.5, 0, 400, np.array([0., 1e5]), n=140)
    assert r.media[-1] == pytest.approx(1., abs=1e-7)
    assert r.centro[-1] == pytest.approx(1., abs=1e-7)


def test_estacionario_con_analitica_esferica():
    for dd, kk, rr in [(8.5, .036, 400), (20., .001, 250), (4., .003, 300)]:
        r = simular(dd, kk, rr, np.array([360000.]), n=220)
        teoria = perfil_estacionario_analitico(r.radios, dd, kk, rr)
        assert np.max(np.abs(r.valores[:, -1] - teoria)) < .003


def test_balance_de_masa_discreto():
    r = simular(8.5, .036, 400, np.array([15.,300.,1200.]), n=140)
    for i in range(r.valores.shape[1]):
        lhs, rhs = diagnostico_balance(8.5, .036, 400, r.valores[:,i])
        # Comparación relativa a un flujo de escala micrométrica y tolerancia
        # absoluta al redondeo de grandes y pequeños términos compensados.
        assert np.isclose(lhs, rhs, atol=1e-7, rtol=1e-10)


def test_captacion_reduce_centro_y_media():
    t = np.array([3600.])
    sin = simular(8.5,0,400,t,n=160)
    con = simular(8.5,.036,400,t,n=160)
    assert sin.media[-1] > con.media[-1]
    assert sin.centro[-1] > con.centro[-1]


def test_difusion_y_radio_tienen_efecto_en_captacion():
    t=np.array([3600.])
    bajo=simular(2,.005,300,t,n=160)
    alto=simular(15,.005,300,t,n=160)
    assert alto.media[-1]>bajo.media[-1]
    pequeno=simular(8,.005,150,t,n=160)
    grande=simular(8,.005,400,t,n=160)
    assert pequeno.media[-1]>grande.media[-1]


def test_volumen_y_mapa():
    r=simular(8,.004,250,np.array([3600.]),n=160)
    assert np.sum(r.volumenes)==pytest.approx(4*np.pi*250**3/3)
    porc=porcentaje_volumen_superior_umbral(r.radios,r.valores[:,-1],250)
    assert 0<=porc<=100
    coords,z=perfil_en_2d(r.radios,r.valores[:,-1],250)
    assert coords.shape==(111,)
    assert z.shape==(111,111)
    assert np.isnan(z[0,0])


def test_argumentos_incorrectos():
    with pytest.raises(ValueError): simular(0,0,100,np.array([0.]))
    with pytest.raises(ValueError): simular(3,-1,100,np.array([0.]))
    with pytest.raises(ValueError): simular(3,0,100,np.array([-1.]))
