#!/usr/bin/env python3
"""
fermentation_benchmark_full.py

One-command CSV generator for a toy ODE-network benchmark of
surrounding-environment--yeast interaction across four domains:

1) bread dough proxy
2) rice wine / simultaneous saccharification-fermentation proxy
3) beer wort fermentation proxy
4) malt / mash-to-wort-to-fermentation proxy

Run once:
    python fermentation_benchmark_full.py

Outputs CSV files only:
    fermentation_benchmark_results/experiment_conditions.csv
    fermentation_benchmark_results/timeseries.csv
    fermentation_benchmark_results/summary.csv
    fermentation_benchmark_results/environment_rankings.csv
    fermentation_benchmark_results/column_dictionary.csv
    fermentation_benchmark_results/manifest.csv

This is NOT a validated industrial fermentation model. It is a lightweight
research scaffold: a Ban-style synchronous ODE network where chemical pools,
enzyme reactions, yeast physiology, and environment variables interact.

No external packages are required.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import sys
from dataclasses import asdict, dataclass, fields
from typing import Dict, Iterable, List, Tuple

EPS = 1e-12

STATE_NAMES = [
    # Non-fermentable / slowly released substrate pools
    "starch_g_L",
    "dextrin_g_L",
    "protein_g_L",
    # Fermentable carbohydrate pools
    "maltose_g_L",
    "maltotriose_g_L",
    "glucose_g_L",
    "fructose_g_L",
    # Yeast nutrients and physical environment
    "nitrogen_g_L",
    "oxygen_mg_L",
    # Organism state
    "biomass_g_L",
    # Products / inhibitors / proxies
    "ethanol_g_L",
    "co2_g_L",
    "retained_co2_g_L",
    "glycerol_g_L",
    "acetate_g_L",
    "acid_mM",
]
IDX = {name: i for i, name in enumerate(STATE_NAMES)}
FERMENTABLE_SUGARS = ["maltose_g_L", "maltotriose_g_L", "glucose_g_L", "fructose_g_L"]
CARBOHYDRATES = ["starch_g_L", "dextrin_g_L"] + FERMENTABLE_SUGARS


@dataclass(frozen=True)
class StrainPhysiology:
    """Yeast-side coefficients; deliberately coarse and interpretable."""

    name: str
    q_sugar_max_g_g_h: float = 2.0

    glucose_affinity_g_L: float = 0.80
    fructose_affinity_g_L: float = 1.40
    maltose_affinity_g_L: float = 2.50
    maltotriose_affinity_g_L: float = 5.00
    maltose_capacity: float = 0.70
    maltotriose_capacity: float = 0.35
    glucose_repression_g_L: float = 3.0

    nitrogen_affinity_g_L: float = 0.05
    sugar_inhibition_g_L: float = 300.0
    ethanol_tolerance_g_L: float = 90.0
    ethanol_inhibition_power: float = 2.0

    temp_opt_C: float = 28.0
    temp_width_C: float = 8.0
    pH_opt: float = 4.8
    pH_width: float = 1.0

    crabtree_strength: float = 0.65
    crabtree_K_g_L: float = 8.0
    acid_pump_strength: float = 1.0
    glycerol_stress_gain: float = 1.0

    death_base_h: float = 0.003
    death_stress_h: float = 0.040


@dataclass(frozen=True)
class YieldParameters:
    biomass_yield_anaerobic_g_g: float = 0.060
    biomass_yield_aerobic_g_g: float = 0.230
    ethanol_yield_g_g: float = 0.470
    co2_yield_g_g: float = 0.490
    glycerol_yield_g_g: float = 0.035
    acetate_yield_g_g: float = 0.012
    nitrogen_yield_biomass_g_g: float = 8.0
    oxygen_mg_per_g_sugar_aerobic: float = 450.0
    acid_mM_per_g_sugar: float = 0.30


@dataclass
class Condition:
    condition_id: str
    environment_type: str
    strain: str

    # Simulation protocol
    hours: float = 120.0
    dt_h: float = 0.05
    record_every_h: float = 1.0

    # Temperature protocol. If pre_saccharification_h > 0, the simulation uses
    # pre_saccharification_temp_C until that time, then fermentation_temp_C.
    fermentation_temp_C: float = 25.0
    pre_saccharification_h: float = 0.0
    pre_saccharification_temp_C: float = 25.0
    inoculation_time_h: float = 0.0
    inoculum_biomass_g_L: float = 0.0

    # Environment
    initial_pH: float = 4.8
    buffer_capacity_mM_per_pH: float = 80.0
    kla_h: float = 0.0
    initial_oxygen_mg_L: float = -1.0  # -1 means infer from kla and temp
    oxygen_sat_mg_L_at_30C: float = 7.6
    evaporation_ethanol_h: float = 0.0

    # Initial state
    starch_g_L: float = 0.0
    dextrin_g_L: float = 0.0
    protein_g_L: float = 0.0
    maltose_g_L: float = 0.0
    maltotriose_g_L: float = 0.0
    glucose_g_L: float = 0.0
    fructose_g_L: float = 0.0
    nitrogen_g_L: float = 0.30
    biomass_g_L: float = 0.30
    ethanol_g_L: float = 0.0
    co2_g_L: float = 0.0
    retained_co2_g_L: float = 0.0
    glycerol_g_L: float = 0.0
    acetate_g_L: float = 0.0
    acid_mM: float = 0.0

    # Enzyme/saccharification module
    amylase_activity: float = 0.0
    glucoamylase_activity: float = 0.0
    protease_activity: float = 0.0
    k_starch_to_dextrin_h: float = 0.02
    k_dextrin_to_maltose_h: float = 0.03
    k_maltose_to_glucose_h: float = 0.00
    k_maltotriose_to_glucose_h: float = 0.00
    k_protein_to_nitrogen_h: float = 0.01
    enzyme_temp_opt_C: float = 35.0
    enzyme_temp_width_C: float = 14.0
    enzyme_pH_opt: float = 5.2
    enzyme_pH_width: float = 1.1

    # Gas handling; useful for bread dough proxy
    co2_retention_fraction: float = 0.10
    co2_escape_h: float = 0.50

    # Extra notes carried to output
    design_axis: str = "default"


STRAINS: Dict[str, StrainPhysiology] = {
    "reference": StrainPhysiology(name="reference"),
    "baker_like": StrainPhysiology(
        name="baker_like",
        q_sugar_max_g_g_h=2.75,
        glucose_affinity_g_L=0.55,
        fructose_affinity_g_L=0.90,
        maltose_capacity=0.55,
        maltotriose_capacity=0.20,
        ethanol_tolerance_g_L=65.0,
        temp_opt_C=31.0,
        temp_width_C=8.5,
        pH_opt=5.4,
        pH_width=1.1,
        crabtree_strength=0.85,
        acid_pump_strength=0.85,
    ),
    "rice_wine_like": StrainPhysiology(
        name="rice_wine_like",
        q_sugar_max_g_g_h=2.05,
        maltose_capacity=0.80,
        maltotriose_capacity=0.35,
        ethanol_tolerance_g_L=115.0,
        temp_opt_C=25.0,
        temp_width_C=9.0,
        pH_opt=4.5,
        pH_width=1.2,
        crabtree_strength=0.70,
        glycerol_stress_gain=1.15,
    ),
    "ale_like": StrainPhysiology(
        name="ale_like",
        q_sugar_max_g_g_h=1.95,
        maltose_affinity_g_L=1.80,
        maltotriose_affinity_g_L=4.00,
        maltose_capacity=1.00,
        maltotriose_capacity=0.55,
        ethanol_tolerance_g_L=95.0,
        temp_opt_C=20.0,
        temp_width_C=6.5,
        pH_opt=4.7,
        pH_width=0.95,
        crabtree_strength=0.62,
    ),
    "lager_like": StrainPhysiology(
        name="lager_like",
        q_sugar_max_g_g_h=1.35,
        maltose_affinity_g_L=1.80,
        maltotriose_affinity_g_L=3.30,
        maltose_capacity=1.05,
        maltotriose_capacity=0.70,
        ethanol_tolerance_g_L=95.0,
        temp_opt_C=12.0,
        temp_width_C=5.5,
        pH_opt=4.7,
        pH_width=0.95,
        crabtree_strength=0.55,
        death_stress_h=0.030,
    ),
    "stress_tolerant": StrainPhysiology(
        name="stress_tolerant",
        q_sugar_max_g_g_h=1.75,
        maltose_capacity=0.75,
        maltotriose_capacity=0.45,
        ethanol_tolerance_g_L=130.0,
        temp_opt_C=27.0,
        temp_width_C=11.0,
        pH_opt=4.6,
        pH_width=1.45,
        glycerol_stress_gain=1.45,
        death_stress_h=0.022,
    ),
    "fast_fermenter": StrainPhysiology(
        name="fast_fermenter",
        q_sugar_max_g_g_h=3.00,
        glucose_affinity_g_L=0.55,
        fructose_affinity_g_L=0.95,
        maltose_capacity=0.72,
        maltotriose_capacity=0.25,
        ethanol_tolerance_g_L=85.0,
        temp_opt_C=30.0,
        temp_width_C=8.0,
        crabtree_strength=0.82,
        acid_pump_strength=1.10,
    ),
}


def clamp(x: float, lo: float, hi: float) -> float:
    return lo if x < lo else hi if x > hi else x


def monod(s: float, k: float) -> float:
    s = max(0.0, s)
    return s / (k + s + EPS)


def gaussian_factor(x: float, opt: float, width: float) -> float:
    return math.exp(-((x - opt) / max(width, EPS)) ** 2)


def current_temperature_C(t_h: float, cond: Condition) -> float:
    if cond.pre_saccharification_h > EPS and t_h < cond.pre_saccharification_h:
        return cond.pre_saccharification_temp_C
    return cond.fermentation_temp_C


def current_phase(t_h: float, cond: Condition) -> str:
    if cond.pre_saccharification_h > EPS and t_h < cond.pre_saccharification_h:
        return "pre_saccharification"
    if t_h < cond.inoculation_time_h - EPS:
        return "pre_inoculation"
    return "fermentation"


def oxygen_saturation_mg_L(temp_C: float, cond: Condition) -> float:
    # Lightweight temperature correction around 30 C; enough for a toy benchmark.
    return max(1.0, cond.oxygen_sat_mg_L_at_30C - 0.20 * (temp_C - 30.0))


def initial_state(cond: Condition) -> List[float]:
    oxygen = cond.initial_oxygen_mg_L
    if oxygen < 0.0:
        if cond.kla_h <= EPS:
            oxygen = 0.05
        else:
            oxygen = 0.50 * oxygen_saturation_mg_L(cond.fermentation_temp_C, cond)
    state = [0.0 for _ in STATE_NAMES]
    for name in STATE_NAMES:
        if name == "oxygen_mg_L":
            state[IDX[name]] = oxygen
        else:
            state[IDX[name]] = float(getattr(cond, name))
    return state


def pH_from_state(y: List[float], cond: Condition) -> float:
    acid = max(y[IDX["acid_mM"]], 0.0)
    return clamp(cond.initial_pH - acid / max(cond.buffer_capacity_mM_per_pH, EPS), 2.2, 7.5)


def fermentable_sugar(y: List[float]) -> float:
    return sum(max(y[IDX[name]], 0.0) for name in FERMENTABLE_SUGARS)


def total_carbohydrate(y: List[float]) -> float:
    return sum(max(y[IDX[name]], 0.0) for name in CARBOHYDRATES)


def enzyme_activity_factor(t_h: float, y: List[float], cond: Condition) -> float:
    temp = current_temperature_C(t_h, cond)
    ph = pH_from_state(y, cond)
    temp_factor = gaussian_factor(temp, cond.enzyme_temp_opt_C, cond.enzyme_temp_width_C)
    ph_factor = gaussian_factor(ph, cond.enzyme_pH_opt, cond.enzyme_pH_width)
    return clamp(temp_factor * ph_factor, 0.0, 1.0)


def rate_factors(t_h: float, y: List[float], phys: StrainPhysiology, cond: Condition) -> Dict[str, float]:
    temp = current_temperature_C(t_h, cond)
    ph = pH_from_state(y, cond)
    sugar = fermentable_sugar(y)
    nitrogen = max(y[IDX["nitrogen_g_L"]], 0.0)
    oxygen = max(y[IDX["oxygen_mg_L"]], 0.0)
    ethanol = max(y[IDX["ethanol_g_L"]], 0.0)

    sugar_lim = monod(sugar, phys.glucose_affinity_g_L)
    n_lim = 0.15 + 0.85 * monod(nitrogen, phys.nitrogen_affinity_g_L)
    substrate_inhib = 1.0 / (1.0 + sugar / max(phys.sugar_inhibition_g_L, EPS))
    ethanol_inhib = clamp(1.0 - (ethanol / max(phys.ethanol_tolerance_g_L, EPS)) ** phys.ethanol_inhibition_power, 0.0, 1.0)
    temp_factor = gaussian_factor(temp, phys.temp_opt_C, phys.temp_width_C)
    ph_factor = gaussian_factor(ph, phys.pH_opt, phys.pH_width)
    oxygen_avail = monod(oxygen, 0.35)

    crabtree = phys.crabtree_strength * monod(sugar, phys.crabtree_K_g_L)
    fermentative_fraction = clamp(max(1.0 - oxygen_avail, crabtree), 0.0, 1.0)
    aerobic_fraction = 1.0 - fermentative_fraction

    activity = sugar_lim * n_lim * substrate_inhib * ethanol_inhib * temp_factor * ph_factor
    stress = stress_index(t_h, y, phys, cond)
    return {
        "temperature_C": temp,
        "pH": ph,
        "fermentable_sugar_g_L": sugar,
        "total_carbohydrate_g_L": total_carbohydrate(y),
        "activity": activity,
        "stress": stress,
        "aerobic_fraction": aerobic_fraction,
        "fermentative_fraction": fermentative_fraction,
        "enzyme_activity_factor": enzyme_activity_factor(t_h, y, cond),
    }


def stress_index(t_h: float, y: List[float], phys: StrainPhysiology, cond: Condition) -> float:
    temp = current_temperature_C(t_h, cond)
    ph = pH_from_state(y, cond)
    sugar = fermentable_sugar(y)
    ethanol = max(y[IDX["ethanol_g_L"]], 0.0)
    temp_bad = 1.0 - gaussian_factor(temp, phys.temp_opt_C, phys.temp_width_C)
    ph_bad = 1.0 - gaussian_factor(ph, phys.pH_opt, phys.pH_width)
    ethanol_bad = clamp((ethanol / max(phys.ethanol_tolerance_g_L, EPS)) ** phys.ethanol_inhibition_power, 0.0, 2.0)
    osmotic_bad = clamp(sugar / max(phys.sugar_inhibition_g_L, EPS), 0.0, 2.0)
    return clamp(0.32 * temp_bad + 0.26 * ph_bad + 0.30 * ethanol_bad + 0.12 * osmotic_bad, 0.0, 2.0)


def sugar_uptake_shares(y: List[float], phys: StrainPhysiology) -> Dict[str, float]:
    glucose = max(y[IDX["glucose_g_L"]], 0.0)
    fructose = max(y[IDX["fructose_g_L"]], 0.0)
    maltose = max(y[IDX["maltose_g_L"]], 0.0)
    maltotriose = max(y[IDX["maltotriose_g_L"]], 0.0)

    repression = 1.0 / (1.0 + glucose / max(phys.glucose_repression_g_L, EPS))
    drives = {
        "glucose_g_L": monod(glucose, phys.glucose_affinity_g_L),
        "fructose_g_L": monod(fructose, phys.fructose_affinity_g_L),
        "maltose_g_L": phys.maltose_capacity * monod(maltose, phys.maltose_affinity_g_L) * repression,
        "maltotriose_g_L": phys.maltotriose_capacity * monod(maltotriose, phys.maltotriose_affinity_g_L) * repression,
    }
    denom = sum(drives.values())
    if denom <= EPS:
        return {name: 0.0 for name in FERMENTABLE_SUGARS}
    return {name: val / denom for name, val in drives.items()}


def rhs(t_h: float, y: List[float], phys: StrainPhysiology, cond: Condition, yields: YieldParameters) -> List[float]:
    y = [max(v, 0.0) for v in y]
    dy = [0.0 for _ in STATE_NAMES]

    # Enzymatic release: starch -> dextrin -> maltose -> glucose, protein -> assimilable N.
    ef = enzyme_activity_factor(t_h, y, cond)
    starch = y[IDX["starch_g_L"]]
    dextrin = y[IDX["dextrin_g_L"]]
    maltose = y[IDX["maltose_g_L"]]
    maltotriose = y[IDX["maltotriose_g_L"]]
    protein = y[IDX["protein_g_L"]]

    starch_to_dextrin = cond.k_starch_to_dextrin_h * cond.amylase_activity * ef * starch
    dextrin_to_maltose = cond.k_dextrin_to_maltose_h * cond.amylase_activity * ef * dextrin
    maltose_to_glucose = cond.k_maltose_to_glucose_h * cond.glucoamylase_activity * ef * maltose
    maltotriose_to_glucose = cond.k_maltotriose_to_glucose_h * cond.glucoamylase_activity * ef * maltotriose
    protein_to_nitrogen = cond.k_protein_to_nitrogen_h * cond.protease_activity * ef * protein

    dy[IDX["starch_g_L"]] -= starch_to_dextrin
    dy[IDX["dextrin_g_L"]] += 0.98 * starch_to_dextrin - dextrin_to_maltose
    dy[IDX["maltose_g_L"]] += 0.92 * dextrin_to_maltose - maltose_to_glucose
    dy[IDX["maltotriose_g_L"]] -= maltotriose_to_glucose
    dy[IDX["glucose_g_L"]] += 1.05 * maltose_to_glucose + 1.05 * maltotriose_to_glucose
    dy[IDX["protein_g_L"]] -= protein_to_nitrogen
    dy[IDX["nitrogen_g_L"]] += 0.16 * protein_to_nitrogen

    # Yeast-mediated metabolism.
    biomass = y[IDX["biomass_g_L"]]
    factors = rate_factors(t_h, y, phys, cond)
    sugar_total = factors["fermentable_sugar_g_L"]
    if biomass <= EPS or sugar_total <= EPS or t_h < cond.inoculation_time_h - EPS:
        sugar_uptake = 0.0
    else:
        q_s = phys.q_sugar_max_g_g_h * factors["activity"]
        sugar_uptake = min(q_s * biomass, sugar_total / 0.05)  # loose stability cap, g/L/h
        sugar_uptake = max(0.0, sugar_uptake)

    shares = sugar_uptake_shares(y, phys)
    for sugar_name, share in shares.items():
        dy[IDX[sugar_name]] -= sugar_uptake * share

    aerobic = factors["aerobic_fraction"]
    fermentative = factors["fermentative_fraction"]
    yx_s = yields.biomass_yield_anaerobic_g_g + (yields.biomass_yield_aerobic_g_g - yields.biomass_yield_anaerobic_g_g) * aerobic
    growth = yx_s * sugar_uptake
    death = (phys.death_base_h + phys.death_stress_h * factors["stress"]) * biomass
    ethanol_prod = yields.ethanol_yield_g_g * sugar_uptake * fermentative
    co2_prod = yields.co2_yield_g_g * sugar_uptake * (0.55 + 0.45 * fermentative)
    glycerol_prod = yields.glycerol_yield_g_g * sugar_uptake * (1.0 + phys.glycerol_stress_gain * factors["stress"])
    acetate_prod = yields.acetate_yield_g_g * sugar_uptake * (0.35 + 0.65 * aerobic)
    nitrogen_use = growth / max(yields.nitrogen_yield_biomass_g_g, EPS)
    oxygen_use = yields.oxygen_mg_per_g_sugar_aerobic * sugar_uptake * aerobic

    temp = factors["temperature_C"]
    oxygen_transfer = cond.kla_h * (oxygen_saturation_mg_L(temp, cond) - y[IDX["oxygen_mg_L"]])
    acid_prod = yields.acid_mM_per_g_sugar * phys.acid_pump_strength * sugar_uptake + 3.0 * acetate_prod

    dy[IDX["nitrogen_g_L"]] -= min(nitrogen_use, y[IDX["nitrogen_g_L"]] / 0.05)
    dy[IDX["oxygen_mg_L"]] += oxygen_transfer - oxygen_use
    dy[IDX["biomass_g_L"]] += growth - death
    dy[IDX["ethanol_g_L"]] += ethanol_prod - cond.evaporation_ethanol_h * y[IDX["ethanol_g_L"]]
    dy[IDX["co2_g_L"]] += co2_prod
    dy[IDX["retained_co2_g_L"]] += cond.co2_retention_fraction * co2_prod - cond.co2_escape_h * y[IDX["retained_co2_g_L"]]
    dy[IDX["glycerol_g_L"]] += glycerol_prod
    dy[IDX["acetate_g_L"]] += acetate_prod
    dy[IDX["acid_mM"]] += acid_prod
    return dy


def euler_step(t_h: float, y: List[float], dt_h: float, phys: StrainPhysiology, cond: Condition, yields: YieldParameters) -> List[float]:
    dy = rhs(t_h, y, phys, cond, yields)
    out = [max(0.0, y[i] + dt_h * dy[i]) for i in range(len(y))]
    return out


def snapshot(t_h: float, y: List[float], phys: StrainPhysiology, cond: Condition) -> Dict[str, float | str]:
    factors = rate_factors(t_h, y, phys, cond)
    row: Dict[str, float | str] = {
        "condition_id": cond.condition_id,
        "environment_type": cond.environment_type,
        "strain": cond.strain,
        "design_axis": cond.design_axis,
        "time_h": round(t_h, 8),
        "phase": current_phase(t_h, cond),
        "temperature_C": factors["temperature_C"],
        "initial_pH": cond.initial_pH,
        "pH": factors["pH"],
        "activity": factors["activity"],
        "stress": factors["stress"],
        "aerobic_fraction": factors["aerobic_fraction"],
        "fermentative_fraction": factors["fermentative_fraction"],
        "enzyme_activity_factor": factors["enzyme_activity_factor"],
        "fermentable_sugar_g_L": factors["fermentable_sugar_g_L"],
        "total_carbohydrate_g_L": factors["total_carbohydrate_g_L"],
    }
    for name in STATE_NAMES:
        row[name] = y[IDX[name]]
    return row


def run_simulation(cond: Condition, yields: YieldParameters | None = None) -> Tuple[List[Dict[str, float | str]], Dict[str, float | str]]:
    if cond.strain not in STRAINS:
        raise ValueError(f"Unknown strain {cond.strain!r}. Available: {', '.join(sorted(STRAINS))}")
    if yields is None:
        yields = YieldParameters()
    phys = STRAINS[cond.strain]
    y = initial_state(cond)
    n_steps = int(math.ceil(cond.hours / cond.dt_h))
    record_every_steps = max(1, int(round(cond.record_every_h / cond.dt_h)))
    records: List[Dict[str, float | str]] = []

    inoculated = cond.inoculation_time_h <= EPS or cond.inoculum_biomass_g_L <= EPS
    if cond.inoculation_time_h <= EPS and cond.inoculum_biomass_g_L > EPS:
        y[IDX["biomass_g_L"]] += cond.inoculum_biomass_g_L
        inoculated = True

    records.append(snapshot(0.0, y, phys, cond))
    for step in range(1, n_steps + 1):
        t_prev = (step - 1) * cond.dt_h
        t_now = step * cond.dt_h
        if (not inoculated) and t_prev < cond.inoculation_time_h <= t_now + EPS:
            y[IDX["biomass_g_L"]] += cond.inoculum_biomass_g_L
            inoculated = True
        y = euler_step(t_prev, y, cond.dt_h, phys, cond, yields)
        if step % record_every_steps == 0 or step == n_steps:
            records.append(snapshot(t_now, y, phys, cond))

    summary = summarize(records, cond, phys)
    return records, summary


def value_at_or_after(records: List[Dict[str, float | str]], metric: str, time_h: float) -> float:
    best = None
    for r in records:
        if float(r["time_h"]) >= time_h:
            best = r
            break
    if best is None:
        best = records[-1]
    try:
        return float(best[metric])
    except Exception:
        return float("nan")


def max_rate(records: List[Dict[str, float | str]], metric: str) -> float:
    best = 0.0
    for prev, cur in zip(records, records[1:]):
        dt = float(cur["time_h"]) - float(prev["time_h"])
        if dt <= EPS:
            continue
        rate = (float(cur[metric]) - float(prev[metric])) / dt
        if rate > best:
            best = rate
    return best


def first_time_below(records: List[Dict[str, float | str]], metric: str, threshold: float) -> float:
    for r in records:
        if float(r[metric]) <= threshold:
            return float(r["time_h"])
    return float("nan")


def mean(values: Iterable[float]) -> float:
    vals = list(values)
    return sum(vals) / len(vals) if vals else float("nan")


def infer_failure_mode(records: List[Dict[str, float | str]], cond: Condition, phys: StrainPhysiology) -> str:
    final = records[-1]
    final_sugar = float(final["fermentable_sugar_g_L"])
    final_starch_dextrin = float(final["starch_g_L"]) + float(final["dextrin_g_L"])
    final_n = float(final["nitrogen_g_L"])
    final_ethanol = float(final["ethanol_g_L"])
    max_stress_val = max(float(r["stress"]) for r in records)
    min_pH_val = min(float(r["pH"]) for r in records)

    # In bread, residual starch is expected; CO2 generation is the target, not
    # complete saccharification. Do not call that a failure.
    if cond.environment_type == "bread_dough":
        if max_stress_val > 0.75:
            return "yeast_stress"
        if final_sugar < 10:
            return "proof_sugar_depleted"
        return "proof_residual_sugar_expected"

    if cond.environment_type == "malt_mash_wort":
        # Malt/wort systems can intentionally retain dextrins. Treat residual
        # starch as a mash-conversion issue, but not residual dextrin alone.
        if float(final["starch_g_L"]) > 30 and cond.amylase_activity > 0.0:
            return "mash_starch_conversion_limited"
    elif final_starch_dextrin > 40 and cond.amylase_activity > 0.0:
        return "saccharification_limited"
    if final_sugar > 20 and final_n < 0.03:
        return "nitrogen_limited"
    if final_sugar > 20 and final_ethanol > 0.70 * phys.ethanol_tolerance_g_L:
        return "ethanol_inhibited"
    if final_sugar > 20 and min_pH_val < 3.3:
        return "acid_stress"
    if final_sugar > 20 and max_stress_val > 0.75:
        return "yeast_stress"
    if final_sugar < 10:
        return "completed"
    return "balanced_or_slow"


def environment_objective(records: List[Dict[str, float | str]], cond: Condition) -> float:
    final = records[-1]
    final_ethanol = float(final["ethanol_g_L"])
    final_sugar = float(final["fermentable_sugar_g_L"])
    final_carb = float(final["total_carbohydrate_g_L"])
    final_acetate = float(final["acetate_g_L"])
    final_glycerol = float(final["glycerol_g_L"])
    max_stress_val = max(float(r["stress"]) for r in records)
    min_pH_val = min(float(r["pH"]) for r in records)

    if cond.environment_type == "bread_dough":
        co2_4h = value_at_or_after(records, "co2_g_L", 4.0)
        retained_4h = value_at_or_after(records, "retained_co2_g_L", 4.0)
        ethanol_4h = value_at_or_after(records, "ethanol_g_L", 4.0)
        return 2.0 * retained_4h + 0.5 * co2_4h - 0.25 * ethanol_4h - 2.0 * max_stress_val
    if cond.environment_type == "rice_wine":
        return final_ethanol - 0.18 * final_sugar - 0.06 * final_carb - 2.0 * final_acetate - 3.0 * max_stress_val
    if cond.environment_type == "beer_wort":
        attenuation = 1.0 - final_sugar / max(1.0, float(records[0]["fermentable_sugar_g_L"]))
        return 40.0 * attenuation + 0.35 * final_ethanol - 2.5 * final_acetate - 0.20 * final_glycerol - 1.5 * abs(min_pH_val - 4.2)
    if cond.environment_type == "malt_mash_wort":
        wort_sugar_at_inoculation = value_at_or_after(records, "fermentable_sugar_g_L", max(cond.inoculation_time_h, cond.pre_saccharification_h))
        return 0.45 * wort_sugar_at_inoculation + 0.45 * final_ethanol - 0.12 * final_carb - 2.0 * max_stress_val
    return final_ethanol - final_sugar


def summarize(records: List[Dict[str, float | str]], cond: Condition, phys: StrainPhysiology) -> Dict[str, float | str]:
    first = records[0]
    final = records[-1]
    initial_fermentable = float(first["fermentable_sugar_g_L"])
    initial_carbohydrate = float(first["total_carbohydrate_g_L"])
    final_fermentable = float(final["fermentable_sugar_g_L"])
    final_carbohydrate = float(final["total_carbohydrate_g_L"])
    sugar_used = max(initial_fermentable - final_fermentable, 0.0)
    generated_or_released_sugar = max(0.0, max(float(r["fermentable_sugar_g_L"]) for r in records) - initial_fermentable)
    pHs = [float(r["pH"]) for r in records]
    stresses = [float(r["stress"]) for r in records]
    activities = [float(r["activity"]) for r in records]
    final_ethanol = float(final["ethanol_g_L"])
    final_co2 = float(final["co2_g_L"])

    summary: Dict[str, float | str] = {
        "condition_id": cond.condition_id,
        "environment_type": cond.environment_type,
        "strain": cond.strain,
        "design_axis": cond.design_axis,
        "hours": cond.hours,
        "dt_h": cond.dt_h,
        "record_every_h": cond.record_every_h,
        "fermentation_temp_C": cond.fermentation_temp_C,
        "pre_saccharification_h": cond.pre_saccharification_h,
        "pre_saccharification_temp_C": cond.pre_saccharification_temp_C,
        "inoculation_time_h": cond.inoculation_time_h,
        "initial_pH": cond.initial_pH,
        "kla_h": cond.kla_h,
        "amylase_activity": cond.amylase_activity,
        "glucoamylase_activity": cond.glucoamylase_activity,
        "protease_activity": cond.protease_activity,
        "initial_carbohydrate_g_L": initial_carbohydrate,
        "initial_fermentable_sugar_g_L": initial_fermentable,
        "max_fermentable_sugar_g_L": max(float(r["fermentable_sugar_g_L"]) for r in records),
        "released_fermentable_sugar_g_L": generated_or_released_sugar,
        "final_carbohydrate_g_L": final_carbohydrate,
        "final_fermentable_sugar_g_L": final_fermentable,
        "fermentable_sugar_used_g_L": sugar_used,
        "final_ethanol_g_L": final_ethanol,
        "final_co2_g_L": final_co2,
        "final_retained_co2_g_L": float(final["retained_co2_g_L"]),
        "final_biomass_g_L": float(final["biomass_g_L"]),
        "final_glycerol_g_L": float(final["glycerol_g_L"]),
        "final_acetate_g_L": float(final["acetate_g_L"]),
        "final_pH": float(final["pH"]),
        "min_pH": min(pHs),
        "max_stress": max(stresses),
        "mean_stress": mean(stresses),
        "mean_activity": mean(activities),
        "max_co2_rate_g_L_h": max_rate(records, "co2_g_L"),
        "max_ethanol_rate_g_L_h": max_rate(records, "ethanol_g_L"),
        "co2_2h_g_L": value_at_or_after(records, "co2_g_L", 2.0),
        "co2_4h_g_L": value_at_or_after(records, "co2_g_L", 4.0),
        "retained_co2_4h_g_L": value_at_or_after(records, "retained_co2_g_L", 4.0),
        "ethanol_24h_g_L": value_at_or_after(records, "ethanol_g_L", 24.0),
        "ethanol_72h_g_L": value_at_or_after(records, "ethanol_g_L", 72.0),
        "ethanol_120h_g_L": value_at_or_after(records, "ethanol_g_L", 120.0),
        "time_to_50pct_initial_fermentable_depletion_h": first_time_below(records, "fermentable_sugar_g_L", 0.50 * initial_fermentable) if initial_fermentable > EPS else float("nan"),
        "time_to_90pct_initial_fermentable_depletion_h": first_time_below(records, "fermentable_sugar_g_L", 0.10 * initial_fermentable) if initial_fermentable > EPS else float("nan"),
        "apparent_ethanol_yield_g_g_initial_fermentable_used": final_ethanol / sugar_used if sugar_used > EPS else float("nan"),
        "apparent_attenuation_fraction": 1.0 - final_fermentable / max(initial_fermentable, 1.0),
        "stuck_fermentation_flag": 1.0 if final_fermentable > max(20.0, 0.20 * max(initial_fermentable, 1.0)) else 0.0,
        "dominant_failure_mode": infer_failure_mode(records, cond, phys),
        "environment_objective_score": environment_objective(records, cond),
    }
    return summary


def write_csv(path: str, rows: Iterable[Dict[str, float | str]]) -> None:
    rows = list(rows)
    if not rows:
        raise ValueError(f"No rows to write for {path}")
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    fieldnames: List[str] = []
    for row in rows:
        for key in row.keys():
            if key not in fieldnames:
                fieldnames.append(key)
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def read_conditions_csv(path: str) -> List[Condition]:
    out: List[Condition] = []
    field_names = {f.name for f in fields(Condition)}
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            kwargs = {}
            for name in field_names:
                if name not in row or row[name] == "":
                    continue
                default_val = getattr(Condition("tmp", "tmp", "reference"), name)
                if isinstance(default_val, str):
                    kwargs[name] = row[name]
                else:
                    kwargs[name] = float(row[name])
            out.append(Condition(**kwargs))
    return out


def condition_to_row(cond: Condition) -> Dict[str, float | str]:
    return asdict(cond)


def make_default_design() -> List[Condition]:
    conditions: List[Condition] = []

    def add(cond: Condition) -> None:
        conditions.append(cond)

    # 1) Bread dough proxy: short proofing, CO2 productivity and retention matter.
    bread_temps = [22.0, 30.0, 35.0]
    bread_sugar_levels = [
        ("low_sugar", 12.0, 3.0, 120.0),
        ("standard_sugar", 25.0, 5.0, 110.0),
        ("sweet_dough", 55.0, 8.0, 90.0),
    ]
    for temp in bread_temps:
        for sugar_label, glucose, fructose, starch in bread_sugar_levels:
            for strain in ["baker_like", "fast_fermenter", "stress_tolerant"]:
                add(Condition(
                    condition_id=f"BREAD_{strain}_{sugar_label}_{int(temp)}C",
                    environment_type="bread_dough",
                    strain=strain,
                    hours=8.0,
                    dt_h=0.01,
                    record_every_h=0.10,
                    fermentation_temp_C=temp,
                    initial_pH=5.6,
                    buffer_capacity_mM_per_pH=120.0,
                    kla_h=0.6,
                    initial_oxygen_mg_L=4.0,
                    starch_g_L=starch,
                    dextrin_g_L=10.0,
                    protein_g_L=20.0,
                    glucose_g_L=glucose,
                    fructose_g_L=fructose,
                    nitrogen_g_L=0.22,
                    biomass_g_L=1.4,
                    amylase_activity=0.18,
                    glucoamylase_activity=0.06,
                    protease_activity=0.04,
                    k_starch_to_dextrin_h=0.012,
                    k_dextrin_to_maltose_h=0.015,
                    k_maltose_to_glucose_h=0.004,
                    enzyme_temp_opt_C=35.0,
                    enzyme_temp_width_C=12.0,
                    enzyme_pH_opt=5.5,
                    enzyme_pH_width=0.9,
                    co2_retention_fraction=0.65,
                    co2_escape_h=0.28,
                    design_axis=f"bread_temp_x_sugar_x_strain:{temp}_{sugar_label}",
                ))

    # 2) Rice wine proxy: simultaneous saccharification-fermentation.
    rice_temps = [18.0, 25.0, 30.0]
    rice_enzyme_levels = [
        ("weak_nuruk", 0.45, 0.40, 0.35),
        ("standard_nuruk", 0.75, 0.70, 0.55),
        ("strong_koji", 1.05, 0.95, 0.70),
    ]
    for temp in rice_temps:
        for enzyme_label, amy, gluco, prot in rice_enzyme_levels:
            for strain in ["rice_wine_like", "stress_tolerant", "fast_fermenter"]:
                add(Condition(
                    condition_id=f"RICE_{strain}_{enzyme_label}_{int(temp)}C",
                    environment_type="rice_wine",
                    strain=strain,
                    hours=192.0,
                    dt_h=0.05,
                    record_every_h=2.0,
                    fermentation_temp_C=temp,
                    initial_pH=4.4,
                    buffer_capacity_mM_per_pH=90.0,
                    kla_h=0.02,
                    initial_oxygen_mg_L=0.5,
                    starch_g_L=190.0,
                    dextrin_g_L=25.0,
                    protein_g_L=12.0,
                    glucose_g_L=12.0,
                    fructose_g_L=0.0,
                    nitrogen_g_L=0.18,
                    biomass_g_L=0.35,
                    amylase_activity=amy,
                    glucoamylase_activity=gluco,
                    protease_activity=prot,
                    k_starch_to_dextrin_h=0.035,
                    k_dextrin_to_maltose_h=0.055,
                    k_maltose_to_glucose_h=0.018,
                    k_maltotriose_to_glucose_h=0.010,
                    enzyme_temp_opt_C=32.0,
                    enzyme_temp_width_C=16.0,
                    enzyme_pH_opt=4.8,
                    enzyme_pH_width=1.1,
                    co2_retention_fraction=0.10,
                    co2_escape_h=1.20,
                    design_axis=f"rice_temp_x_saccharification_x_strain:{temp}_{enzyme_label}",
                ))

    # 3) Beer wort proxy: maltose/maltotriose-rich wort, oxygen at pitch, low ongoing O2.
    beer_styles = [
        ("lager_cool", 12.0, "lager_like"),
        ("ale_cool", 18.0, "ale_like"),
        ("ale_warm", 22.0, "ale_like"),
    ]
    wort_levels = [
        ("session_wort", 48.0, 13.0, 8.0, 3.0),
        ("standard_wort", 74.0, 20.0, 12.0, 4.0),
        ("high_gravity_wort", 112.0, 34.0, 18.0, 6.0),
    ]
    oxygen_levels = [("low_o2", 2.0, 0.02), ("aerated_pitch", 8.0, 0.08)]
    for style_label, temp, primary_strain in beer_styles:
        strain_options = [primary_strain, "stress_tolerant"]
        for wort_label, maltose, maltotriose, glucose, dextrin in wort_levels:
            for oxygen_label, initial_o2, kla in oxygen_levels:
                for strain in strain_options:
                    add(Condition(
                        condition_id=f"BEER_{strain}_{style_label}_{wort_label}_{oxygen_label}",
                        environment_type="beer_wort",
                        strain=strain,
                        hours=240.0,
                        dt_h=0.05,
                        record_every_h=2.0,
                        fermentation_temp_C=temp,
                        initial_pH=5.2,
                        buffer_capacity_mM_per_pH=110.0,
                        kla_h=kla,
                        initial_oxygen_mg_L=initial_o2,
                        dextrin_g_L=dextrin,
                        maltose_g_L=maltose,
                        maltotriose_g_L=maltotriose,
                        glucose_g_L=glucose,
                        protein_g_L=5.0,
                        nitrogen_g_L=0.28,
                        biomass_g_L=0.45,
                        amylase_activity=0.0,
                        glucoamylase_activity=0.0,
                        protease_activity=0.03,
                        k_protein_to_nitrogen_h=0.005,
                        co2_retention_fraction=0.05,
                        co2_escape_h=1.50,
                        design_axis=f"beer_style_x_wort_x_oxygen_x_strain:{style_label}_{wort_label}_{oxygen_label}",
                    ))

    # 4) Malt/mash-to-wort-to-fermentation proxy: pre-saccharification phase, then inoculation.
    mash_protocols = [
        ("low_mash_55C", 55.0, 2.0, 0.70, 0.45),
        ("standard_mash_62C", 62.0, 2.0, 0.95, 0.65),
        ("hot_mash_68C", 68.0, 2.0, 0.85, 0.50),
    ]
    ferment_temps = [("ale18", 18.0, "ale_like"), ("ale22", 22.0, "ale_like"), ("robust25", 25.0, "stress_tolerant")]
    grist_levels = [("standard_grist", 160.0, 18.0), ("heavy_grist", 230.0, 25.0)]
    for mash_label, mash_temp, mash_h, amy, gluco in mash_protocols:
        for ferm_label, ferment_temp, strain in ferment_temps:
            for grist_label, starch, protein in grist_levels:
                add(Condition(
                    condition_id=f"MALT_{strain}_{mash_label}_{ferm_label}_{grist_label}",
                    environment_type="malt_mash_wort",
                    strain=strain,
                    hours=168.0,
                    dt_h=0.02,
                    record_every_h=1.0,
                    fermentation_temp_C=ferment_temp,
                    pre_saccharification_h=mash_h,
                    pre_saccharification_temp_C=mash_temp,
                    inoculation_time_h=mash_h,
                    inoculum_biomass_g_L=0.45,
                    initial_pH=5.4,
                    buffer_capacity_mM_per_pH=120.0,
                    kla_h=0.04,
                    initial_oxygen_mg_L=6.0,
                    starch_g_L=starch,
                    dextrin_g_L=12.0,
                    protein_g_L=protein,
                    nitrogen_g_L=0.10,
                    biomass_g_L=0.0,
                    amylase_activity=amy,
                    glucoamylase_activity=gluco,
                    protease_activity=0.55,
                    k_starch_to_dextrin_h=1.10,
                    k_dextrin_to_maltose_h=1.15,
                    k_maltose_to_glucose_h=0.025,
                    k_maltotriose_to_glucose_h=0.010,
                    k_protein_to_nitrogen_h=0.035,
                    enzyme_temp_opt_C=62.0,
                    enzyme_temp_width_C=10.0,
                    enzyme_pH_opt=5.4,
                    enzyme_pH_width=0.75,
                    co2_retention_fraction=0.05,
                    co2_escape_h=1.50,
                    design_axis=f"malt_mash_x_ferment_x_grist:{mash_label}_{ferm_label}_{grist_label}",
                ))

    return conditions


def make_rankings(summary_rows: List[Dict[str, float | str]]) -> List[Dict[str, float | str]]:
    grouped: Dict[str, List[Dict[str, float | str]]] = {}
    for row in summary_rows:
        grouped.setdefault(str(row["environment_type"]), []).append(row)
    ranked: List[Dict[str, float | str]] = []
    for env, rows in grouped.items():
        rows_sorted = sorted(rows, key=lambda r: float(r["environment_objective_score"]), reverse=True)
        for rank, row in enumerate(rows_sorted, start=1):
            ranked.append({
                "environment_type": env,
                "rank": rank,
                "condition_id": row["condition_id"],
                "strain": row["strain"],
                "design_axis": row["design_axis"],
                "environment_objective_score": row["environment_objective_score"],
                "final_ethanol_g_L": row["final_ethanol_g_L"],
                "final_co2_g_L": row["final_co2_g_L"],
                "retained_co2_4h_g_L": row["retained_co2_4h_g_L"],
                "final_fermentable_sugar_g_L": row["final_fermentable_sugar_g_L"],
                "final_carbohydrate_g_L": row["final_carbohydrate_g_L"],
                "min_pH": row["min_pH"],
                "max_stress": row["max_stress"],
                "dominant_failure_mode": row["dominant_failure_mode"],
            })
    return ranked


def make_column_dictionary() -> List[Dict[str, str]]:
    descriptions = {
        "experiment_conditions.csv": "All condition parameters used to generate the benchmark sweep.",
        "timeseries.csv": "One row per recorded time point per condition.",
        "summary.csv": "One row per condition with final values and paper-ready metrics.",
        "environment_rankings.csv": "Conditions ranked within each environment by environment-specific objective score.",
        "manifest.csv": "Run metadata and file list.",
    }
    columns = [
        ("condition_id", "Unique identifier for one simulated experimental condition."),
        ("environment_type", "One of bread_dough, rice_wine, beer_wort, malt_mash_wort."),
        ("strain", "Yeast physiology preset used by the condition."),
        ("time_h", "Simulation time in hours."),
        ("phase", "pre_saccharification, pre_inoculation, or fermentation."),
        ("pH", "Buffered pH proxy computed from acid_mM and buffer capacity."),
        ("activity", "Composite yeast activity factor from substrate, N, ethanol, pH, and temperature."),
        ("stress", "Composite stress index from temperature, pH, ethanol, and osmotic stress."),
        ("enzyme_activity_factor", "Composite enzyme activity factor from temperature and pH."),
        ("fermentable_sugar_g_L", "maltose + maltotriose + glucose + fructose."),
        ("total_carbohydrate_g_L", "starch + dextrin + fermentable sugars."),
        ("retained_co2_g_L", "Bread-useful CO2 retention proxy, not a literal gas volume."),
        ("environment_objective_score", "Environment-specific scalar score for ranking conditions."),
        ("dominant_failure_mode", "Rule-based label such as completed, saccharification_limited, nitrogen_limited, yeast_stress."),
    ]
    rows: List[Dict[str, str]] = []
    for file_name, description in descriptions.items():
        rows.append({"file": file_name, "column": "__file__", "description": description})
    for column, description in columns:
        rows.append({"file": "multiple", "column": column, "description": description})
    for state in STATE_NAMES:
        rows.append({"file": "timeseries.csv", "column": state, "description": "State variable in g/L, mg/L for oxygen, or mM for acid_mM."})
    return rows


def run_benchmark(outdir: str, conditions: List[Condition], max_conditions: int | None = None) -> Dict[str, int | str]:
    os.makedirs(outdir, exist_ok=True)
    if max_conditions is not None:
        conditions = conditions[:max_conditions]

    condition_rows = [condition_to_row(c) for c in conditions]
    write_csv(os.path.join(outdir, "experiment_conditions.csv"), condition_rows)

    all_records: List[Dict[str, float | str]] = []
    summary_rows: List[Dict[str, float | str]] = []
    for i, cond in enumerate(conditions, start=1):
        records, summary = run_simulation(cond)
        all_records.extend(records)
        summary_rows.append(summary)
        if i % 20 == 0:
            print(f"simulated {i}/{len(conditions)} conditions", file=sys.stderr)

    ranking_rows = make_rankings(summary_rows)
    write_csv(os.path.join(outdir, "timeseries.csv"), all_records)
    write_csv(os.path.join(outdir, "summary.csv"), summary_rows)
    write_csv(os.path.join(outdir, "environment_rankings.csv"), ranking_rows)
    write_csv(os.path.join(outdir, "column_dictionary.csv"), make_column_dictionary())

    manifest = [
        {"key": "outdir", "value": outdir},
        {"key": "n_conditions", "value": len(conditions)},
        {"key": "n_timeseries_rows", "value": len(all_records)},
        {"key": "n_summary_rows", "value": len(summary_rows)},
        {"key": "script", "value": os.path.basename(__file__)},
        {"key": "model_note", "value": "Toy ODE-network benchmark; not calibrated for production fermentation."},
        {"key": "file", "value": "experiment_conditions.csv"},
        {"key": "file", "value": "timeseries.csv"},
        {"key": "file", "value": "summary.csv"},
        {"key": "file", "value": "environment_rankings.csv"},
        {"key": "file", "value": "column_dictionary.csv"},
    ]
    write_csv(os.path.join(outdir, "manifest.csv"), manifest)
    return {"outdir": outdir, "n_conditions": len(conditions), "n_timeseries_rows": len(all_records)}


def main() -> None:
    parser = argparse.ArgumentParser(description="Run a four-domain environment--yeast fermentation benchmark and write CSV outputs only.")
    parser.add_argument("--outdir", default="fermentation_benchmark_results", help="Output directory for CSV files.")
    parser.add_argument("--conditions", default=None, help="Optional custom experiment_conditions.csv to run instead of the default design.")
    parser.add_argument("--write-default-design-only", action="store_true", help="Only write the default experiment_conditions.csv, then stop.")
    parser.add_argument("--max-conditions", type=int, default=None, help="Debug option: run only the first N conditions.")
    args = parser.parse_args()

    if args.conditions:
        conditions = read_conditions_csv(args.conditions)
    else:
        conditions = make_default_design()

    os.makedirs(args.outdir, exist_ok=True)
    if args.write_default_design_only:
        write_csv(os.path.join(args.outdir, "experiment_conditions.csv"), [condition_to_row(c) for c in conditions])
        print(f"Wrote {len(conditions)} default conditions to {os.path.join(args.outdir, 'experiment_conditions.csv')}")
        return

    info = run_benchmark(args.outdir, conditions, max_conditions=args.max_conditions)
    print(f"Done. Wrote CSV outputs to: {info['outdir']}")
    print(f"Conditions: {info['n_conditions']}; time-series rows: {info['n_timeseries_rows']}")


if __name__ == "__main__":
    main()
