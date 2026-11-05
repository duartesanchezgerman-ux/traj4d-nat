"""Fase 3: validación del modelo Poll–Schumann sin BADA, frente a la calculadora EMEP/EEA (EUROCONTROL).

5-11-2026: se quitó la comparación con la ICAO CORSIA CERT. La OACI no respondió a la consulta sobre su licencia
(solo investigación, sin uso comercial, propiedad intelectual de la OACI), así que el artículo usa solo la
referencia de EUROCONTROL. Este script ya no lee la CERT.
Ejecutar:  PYTHONPATH=src python scripts/f3_validacion_modelo.py | tee resultados/resultados_f3_validacion_modelo.txt

Referencias (en datos/referencias/, ver docs/DATOS.md):
  * EEA/EMEP Master emissions calculator 2023 (EUROCONTROL, método AEM basado en BADA): combustible de
    ascenso-crucero-descenso (CCD) más LTO (despegue y aterrizaje) por tipo y distancia (en NM).
Simulación con Poll–Schumann (pycontrails PSFlight): perfil sintético en atmósfera ISA y sin viento.
Ascenso a 1800 ft/min, crucero FL380 a M0,85, descenso a 2000 ft/min, entre 1500 ft y 1500 ft
(sin rodaje ni LTO completo). PSFlight estima la masa iterando el combustible del viaje y las reservas,
con factor de carga de pasajeros por defecto (~0,81).
Diferencia esperable por construcción: la simulación no incluye rodaje ni el LTO completo,
así que infravalora los vuelos cortos.
"""
import warnings; warnings.filterwarnings('ignore')
import numpy as np, pandas as pd, openpyxl
from pycontrails import Flight
from pycontrails.models.ps_model import PSFlight

from traj4d import atmosfera as atm
NM=1852.0
def perfil(dist_km, fl_cru, mach, roc=1800, rod=2000):
    dt=60.0; alt=[1500.0]; s=[0.0]; tas=[]
    h_c=fl_cru*100; D=dist_km*1e3; fase='asc'
    def v(h):
        if h<10000: return 250*atm.KT*(1+h/100000)
        cas_tas = 300*atm.KT*np.sqrt(atm.presion(0)/atm.presion(h*atm.FT))**1.0*np.sqrt(atm.temperatura(h*atm.FT)/288.15)
        return float(min(atm.mach_a_tas(mach,h*atm.FT), cas_tas))
    while True:
        h=alt[-1]; vt=v(h)
        dist_desc=(h-1500)/rod*60*v(max(h*0.5,1500))
        if fase!='des' and D-s[-1]<=dist_desc: fase='des'
        if fase=='asc':
            h=min(h+roc*dt/60,h_c); fase='cru' if h>=h_c else 'asc'
        elif fase=='des': h=max(h-rod*dt/60,1500)
        s.append(s[-1]+vt*dt); alt.append(h); tas.append(vt)
        if fase=='des' and h<=1500: break
    tas.append(tas[-1])
    return np.array(s), np.array(alt), np.array(tas)
def vuelo(tipo, dist_km, fl_cru=380, mach=0.85, lf=None):
    s,alt,tas=perfil(dist_km,fl_cru,mach); n=len(s)
    fl=Flight(longitude=np.degrees(s/6371008.8), latitude=np.zeros(n), altitude_ft=alt,
              time=pd.date_range('2024-01-01',periods=n,freq='60s'), aircraft_type=tipo)
    fl['true_airspeed']=tas; fl['air_temperature']=atm.temperatura(alt*atm.FT)
    if lf is not None: fl.attrs.update(passenger_load_factor=lf[0], cargo_load_factor=lf[1])
    out=PSFlight(fill_low_altitude_with_isa_temperature=True, fill_low_altitude_with_zero_wind=True).eval(fl)
    ff=np.asarray(out['fuel_flow']); dt=np.diff(out['time']).astype('timedelta64[s]').astype(float)
    return float(np.nansum(ff[:-1]*dt)), dict(out.attrs), float(np.asarray(out['aircraft_mass'])[0])
wb2=openpyxl.load_workbook('datos/referencias/1.A.3.a Aviation  -Annex 1 - Master emissions calculator - 2023 - Protected - v1.5_18_09_2024.xlsx',read_only=True,data_only=True)
eea={}
for r in wb2['database'].iter_rows(values_only=True):
    if r and r[2] in ('B789','A359') and r[9] in ('CCD','LTO'):
        eea.setdefault(r[2],{}).setdefault(r[9],{})[r[11]]=r[13]
print('LTO EEA (kg):', {k: eea[k].get('LTO',{}).get(125) for k in eea}, '| claves CCD B789:', sorted(eea['B789']['CCD'])[:3], '...', sorted(eea['B789']['CCD'])[-3:])
filas=[]
for tipo in ('B789','A359'):
    lto=eea[tipo]['LTO'].get(125,np.nan) if 'LTO' in eea[tipo] else np.nan
    for dnm in (1500,2500,3000,3500,4500,5500):
        dkm=dnm*NM/1e3
        f_ps,at,m0=vuelo(tipo,dkm)
        ccd=eea[tipo]['CCD'].get(dnm,np.nan)
        f_eea=ccd+(lto if np.isfinite(lto) else 0)
        filas.append(dict(tipo=tipo,dist_nm=dnm,dist_km=round(dkm),PS=round(f_ps),EEA=round(f_eea),
            PS_vs_EEA=round(100*(f_ps/f_eea-1),1),
            masa_despegue_t=round(m0/1e3,1), lf_pax=at.get('passenger_load_factor')))
t=pd.DataFrame(filas); print(t.to_string()); t.to_csv("datos/procesado/validacion_modelo_ps.csv", index=False)
