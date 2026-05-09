#!/usr/bin/env python3
"""
visualize_ode_network.py

Create paper-ready ODE-network diagrams for fermentation_benchmark_full.py.

Important: fermentation_benchmark_full.py does not store the ODE network as an
explicit graph object. The network is implicit in:

  1. STATE_NAMES                    -> dynamic state nodes
  2. rate_factors(), stress_index() -> regulatory / environmental factors
  3. rhs()                          -> ODE fluxes and dy/dt updates
  4. euler_step()                   -> synchronous numerical update

This script makes that implicit RHS structure explicit as nodes and directed
edges, then writes Graphviz DOT, SVG, PNG, and edge-list CSV files.

Run from the repository root:

    python visualize_ode_network.py

Outputs:

    ode_network_visualization/ode_network_core.dot
    ode_network_visualization/ode_network_core.svg
    ode_network_visualization/ode_network_core.png
    ode_network_visualization/ode_network_full.dot
    ode_network_visualization/ode_network_full.svg
    ode_network_visualization/ode_network_full.png
    ode_network_visualization/ode_network_full_edges.csv
    ode_network_visualization/ode_network_full_nodes.csv

The core graph is intended for the main paper figure. The full graph is useful
for supplementary material or debugging the model structure.
"""

from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Node:
    node_id: str
    label: str
    group: str
    kind: str = "state"  # state, factor, flux, source, sink


@dataclass(frozen=True)
class Edge:
    source: str
    target: str
    label: str
    relation: str
    sign: str = "+"  # +, -, +/-, mod
    module: str = ""


# ---------------------------------------------------------------------------
# Model-derived node sets
# ---------------------------------------------------------------------------

FALLBACK_STATE_NAMES = [
    "starch_g_L",
    "dextrin_g_L",
    "protein_g_L",
    "maltose_g_L",
    "maltotriose_g_L",
    "glucose_g_L",
    "fructose_g_L",
    "nitrogen_g_L",
    "oxygen_mg_L",
    "biomass_g_L",
    "ethanol_g_L",
    "co2_g_L",
    "retained_co2_g_L",
    "glycerol_g_L",
    "acetate_g_L",
    "acid_mM",
]


def load_state_names() -> List[str]:
    """Import STATE_NAMES from fermentation_benchmark_full.py if available."""
    try:
        import fermentation_benchmark_full as sim  # type: ignore
        return list(sim.STATE_NAMES)
    except Exception:
        return list(FALLBACK_STATE_NAMES)


def pretty_label(node_id: str) -> str:
    labels = {
        "starch_g_L": "Starch",
        "dextrin_g_L": "Dextrin",
        "protein_g_L": "Protein",
        "maltose_g_L": "Maltose",
        "maltotriose_g_L": "Maltotriose",
        "glucose_g_L": "Glucose",
        "fructose_g_L": "Fructose",
        "nitrogen_g_L": "Assimilable N",
        "oxygen_mg_L": "Dissolved O₂",
        "biomass_g_L": "Yeast biomass",
        "ethanol_g_L": "Ethanol",
        "co2_g_L": "CO₂",
        "retained_co2_g_L": "Retained CO₂",
        "glycerol_g_L": "Glycerol",
        "acetate_g_L": "Acetate",
        "acid_mM": "Organic acid",
        "temperature_C": "Temperature",
        "pH": "pH",
        "enzyme_activity_factor": "Enzyme activity",
        "fermentable_sugar_pool": "Fermentable sugar pool",
        "activity": "Yeast activity",
        "stress": "Stress index",
        "aerobic_fraction": "Aerobic fraction",
        "fermentative_fraction": "Fermentative fraction",
        "sugar_uptake_flux": "Sugar uptake flux",
        "oxygen_transfer_flux": "O₂ transfer flux",
        "growth_flux": "Growth flux",
        "death_flux": "Death flux",
        "gas_retention_flux": "Gas retention",
        "environment_oxygen": "External O₂",
        "ethanol_evaporation_sink": "Ethanol loss",
        "co2_escape_sink": "CO₂ escape",
    }
    return labels.get(node_id, node_id.replace("_g_L", "").replace("_mM", "").replace("_", " "))


def node_group(node_id: str) -> str:
    if node_id in {"starch_g_L", "dextrin_g_L", "protein_g_L", "maltose_g_L", "maltotriose_g_L", "glucose_g_L", "fructose_g_L", "fermentable_sugar_pool"}:
        return "substrate"
    if node_id in {"nitrogen_g_L", "oxygen_mg_L", "temperature_C", "pH", "environment_oxygen"}:
        return "environment"
    if node_id in {"enzyme_activity_factor"}:
        return "enzyme"
    if node_id in {"biomass_g_L", "activity", "stress", "aerobic_fraction", "fermentative_fraction", "sugar_uptake_flux", "growth_flux", "death_flux"}:
        return "yeast"
    if node_id in {"ethanol_g_L", "co2_g_L", "retained_co2_g_L", "glycerol_g_L", "acetate_g_L", "acid_mM", "gas_retention_flux", "ethanol_evaporation_sink", "co2_escape_sink"}:
        return "products"
    return "other"


def node_kind(node_id: str, state_names: Sequence[str]) -> str:
    if node_id in state_names:
        return "state"
    if node_id.endswith("_flux") or node_id in {"sugar_uptake_flux", "growth_flux", "death_flux", "gas_retention_flux"}:
        return "flux"
    if node_id in {"environment_oxygen"}:
        return "source"
    if node_id.endswith("_sink"):
        return "sink"
    return "factor"


# ---------------------------------------------------------------------------
# Explicit network extracted from rhs(), rate_factors(), stress_index()
# ---------------------------------------------------------------------------

def build_full_edges() -> List[Edge]:
    e: List[Edge] = []

    # Enzymatic saccharification / hydrolysis module from rhs().
    e += [
        Edge("temperature_C", "enzyme_activity_factor", "enzyme temp factor", "modulates", "mod", "enzyme"),
        Edge("pH", "enzyme_activity_factor", "enzyme pH factor", "modulates", "mod", "enzyme"),
        Edge("enzyme_activity_factor", "starch_g_L", "starch hydrolysis", "consumes", "-", "enzyme"),
        Edge("enzyme_activity_factor", "dextrin_g_L", "dextrin hydrolysis", "modulates", "mod", "enzyme"),
        Edge("enzyme_activity_factor", "protein_g_L", "proteolysis", "consumes", "-", "enzyme"),
        Edge("starch_g_L", "dextrin_g_L", "starch→dextrin", "conversion", "+", "enzyme"),
        Edge("dextrin_g_L", "maltose_g_L", "dextrin→maltose", "conversion", "+", "enzyme"),
        Edge("maltose_g_L", "glucose_g_L", "maltose→glucose", "conversion", "+", "enzyme"),
        Edge("maltotriose_g_L", "glucose_g_L", "maltotriose→glucose", "conversion", "+", "enzyme"),
        Edge("protein_g_L", "nitrogen_g_L", "protein→N", "conversion", "+", "enzyme"),
    ]

    # Fermentable sugar aggregation and uptake sharing.
    for sugar in ["maltose_g_L", "maltotriose_g_L", "glucose_g_L", "fructose_g_L"]:
        e.append(Edge(sugar, "fermentable_sugar_pool", "pool", "aggregate", "+", "yeast"))
        e.append(Edge(sugar, "sugar_uptake_flux", "uptake share", "consumed_by", "-", "yeast"))

    # Glucose repression of complex-sugar uptake.
    e += [
        Edge("glucose_g_L", "maltose_g_L", "glucose repression", "inhibits uptake", "-", "yeast"),
        Edge("glucose_g_L", "maltotriose_g_L", "glucose repression", "inhibits uptake", "-", "yeast"),
    ]

    # Rate factors from rate_factors() and stress_index().
    e += [
        Edge("fermentable_sugar_pool", "activity", "sugar limitation", "modulates", "+", "yeast"),
        Edge("nitrogen_g_L", "activity", "N limitation", "modulates", "+", "yeast"),
        Edge("ethanol_g_L", "activity", "ethanol inhibition", "inhibits", "-", "yeast"),
        Edge("temperature_C", "activity", "temp factor", "modulates", "mod", "yeast"),
        Edge("pH", "activity", "pH factor", "modulates", "mod", "yeast"),
        Edge("oxygen_mg_L", "aerobic_fraction", "O₂ availability", "modulates", "+", "yeast"),
        Edge("oxygen_mg_L", "fermentative_fraction", "low-O₂ fermentation", "modulates", "-", "yeast"),
        Edge("fermentable_sugar_pool", "fermentative_fraction", "Crabtree effect", "modulates", "+", "yeast"),
        Edge("temperature_C", "stress", "thermal stress", "modulates", "mod", "yeast"),
        Edge("pH", "stress", "pH stress", "modulates", "mod", "yeast"),
        Edge("ethanol_g_L", "stress", "ethanol stress", "increases", "+", "yeast"),
        Edge("fermentable_sugar_pool", "stress", "osmotic stress", "increases", "+", "yeast"),
        Edge("biomass_g_L", "sugar_uptake_flux", "biomass term", "drives", "+", "yeast"),
        Edge("activity", "sugar_uptake_flux", "q_s activity", "drives", "+", "yeast"),
    ]

    # Yeast-mediated fluxes from rhs().
    e += [
        Edge("sugar_uptake_flux", "growth_flux", "Yx/s", "drives", "+", "yeast"),
        Edge("aerobic_fraction", "growth_flux", "aerobic yield", "modulates", "+", "yeast"),
        Edge("growth_flux", "biomass_g_L", "growth", "produces", "+", "yeast"),
        Edge("stress", "death_flux", "stress death", "drives", "+", "yeast"),
        Edge("biomass_g_L", "death_flux", "death term", "drives", "+", "yeast"),
        Edge("death_flux", "biomass_g_L", "death", "decreases", "-", "yeast"),
        Edge("growth_flux", "nitrogen_g_L", "N assimilation", "consumes", "-", "yeast"),
        Edge("sugar_uptake_flux", "oxygen_mg_L", "aerobic O₂ use", "consumes", "-", "yeast"),
        Edge("aerobic_fraction", "oxygen_mg_L", "aerobic O₂ demand", "consumes", "-", "yeast"),
        Edge("sugar_uptake_flux", "ethanol_g_L", "ethanol yield", "produces", "+", "products"),
        Edge("fermentative_fraction", "ethanol_g_L", "fermentative yield", "modulates", "+", "products"),
        Edge("sugar_uptake_flux", "co2_g_L", "CO₂ yield", "produces", "+", "products"),
        Edge("fermentative_fraction", "co2_g_L", "fermentative CO₂", "modulates", "+", "products"),
        Edge("sugar_uptake_flux", "glycerol_g_L", "glycerol yield", "produces", "+", "products"),
        Edge("stress", "glycerol_g_L", "stress glycerol", "modulates", "+", "products"),
        Edge("sugar_uptake_flux", "acetate_g_L", "acetate yield", "produces", "+", "products"),
        Edge("aerobic_fraction", "acetate_g_L", "aerobic acetate", "modulates", "+", "products"),
        Edge("sugar_uptake_flux", "acid_mM", "acidification", "produces", "+", "products"),
        Edge("acetate_g_L", "acid_mM", "acetate acidity", "increases", "+", "products"),
    ]

    # Environmental feedback and gas handling.
    e += [
        Edge("environment_oxygen", "oxygen_transfer_flux", "kLa transfer", "drives", "+", "environment"),
        Edge("oxygen_transfer_flux", "oxygen_mg_L", "O₂ transfer", "produces", "+", "environment"),
        Edge("acid_mM", "pH", "buffered acidification", "decreases", "-", "environment"),
        Edge("co2_g_L", "gas_retention_flux", "retention fraction", "drives", "+", "products"),
        Edge("gas_retention_flux", "retained_co2_g_L", "retained CO₂", "produces", "+", "products"),
        Edge("retained_co2_g_L", "co2_escape_sink", "CO₂ escape", "loss", "-", "products"),
        Edge("ethanol_g_L", "ethanol_evaporation_sink", "evaporation", "loss", "-", "products"),
    ]
    return e


def build_core_edges() -> List[Edge]:
    """Smaller network for a main-text figure."""
    return [
        Edge("starch_g_L", "dextrin_g_L", "hydrolysis", "conversion", "+", "enzyme"),
        Edge("dextrin_g_L", "maltose_g_L", "saccharification", "conversion", "+", "enzyme"),
        Edge("maltose_g_L", "glucose_g_L", "glucoamylase", "conversion", "+", "enzyme"),
        Edge("protein_g_L", "nitrogen_g_L", "proteolysis", "conversion", "+", "enzyme"),
        Edge("temperature_C", "enzyme_activity_factor", "enzyme temp", "modulates", "mod", "enzyme"),
        Edge("pH", "enzyme_activity_factor", "enzyme pH", "modulates", "mod", "enzyme"),
        Edge("enzyme_activity_factor", "dextrin_g_L", "enzyme-mediated release", "modulates", "mod", "enzyme"),
        Edge("maltose_g_L", "fermentable_sugar_pool", "pool", "aggregate", "+", "yeast"),
        Edge("maltotriose_g_L", "fermentable_sugar_pool", "pool", "aggregate", "+", "yeast"),
        Edge("glucose_g_L", "fermentable_sugar_pool", "pool", "aggregate", "+", "yeast"),
        Edge("fructose_g_L", "fermentable_sugar_pool", "pool", "aggregate", "+", "yeast"),
        Edge("fermentable_sugar_pool", "sugar_uptake_flux", "substrate", "drives", "+", "yeast"),
        Edge("biomass_g_L", "sugar_uptake_flux", "biomass", "drives", "+", "yeast"),
        Edge("nitrogen_g_L", "activity", "N limitation", "modulates", "+", "yeast"),
        Edge("oxygen_mg_L", "aerobic_fraction", "O₂ availability", "modulates", "+", "yeast"),
        Edge("oxygen_mg_L", "fermentative_fraction", "low O₂", "modulates", "-", "yeast"),
        Edge("temperature_C", "activity", "temp", "modulates", "mod", "yeast"),
        Edge("pH", "activity", "pH", "modulates", "mod", "yeast"),
        Edge("ethanol_g_L", "activity", "ethanol inhibition", "inhibits", "-", "yeast"),
        Edge("activity", "sugar_uptake_flux", "q_s", "drives", "+", "yeast"),
        Edge("sugar_uptake_flux", "biomass_g_L", "growth", "produces", "+", "yeast"),
        Edge("sugar_uptake_flux", "ethanol_g_L", "ethanol", "produces", "+", "products"),
        Edge("sugar_uptake_flux", "co2_g_L", "CO₂", "produces", "+", "products"),
        Edge("co2_g_L", "retained_co2_g_L", "retention", "produces", "+", "products"),
        Edge("sugar_uptake_flux", "glycerol_g_L", "glycerol", "produces", "+", "products"),
        Edge("sugar_uptake_flux", "acetate_g_L", "acetate", "produces", "+", "products"),
        Edge("sugar_uptake_flux", "acid_mM", "acid", "produces", "+", "products"),
        Edge("acetate_g_L", "acid_mM", "acidity", "increases", "+", "products"),
        Edge("acid_mM", "pH", "acidification", "decreases", "-", "environment"),
        Edge("ethanol_g_L", "stress", "product stress", "increases", "+", "yeast"),
        Edge("stress", "biomass_g_L", "death", "decreases", "-", "yeast"),
    ]


def build_nodes(edges: Iterable[Edge], state_names: Sequence[str]) -> List[Node]:
    ids = set()
    for edge in edges:
        ids.add(edge.source)
        ids.add(edge.target)
    # Preserve a useful order.
    preferred = [
        "starch_g_L", "dextrin_g_L", "protein_g_L", "maltose_g_L", "maltotriose_g_L", "glucose_g_L", "fructose_g_L", "fermentable_sugar_pool",
        "temperature_C", "pH", "nitrogen_g_L", "oxygen_mg_L", "environment_oxygen", "enzyme_activity_factor",
        "biomass_g_L", "activity", "stress", "aerobic_fraction", "fermentative_fraction", "sugar_uptake_flux", "growth_flux", "death_flux", "oxygen_transfer_flux",
        "ethanol_g_L", "co2_g_L", "retained_co2_g_L", "glycerol_g_L", "acetate_g_L", "acid_mM", "gas_retention_flux", "ethanol_evaporation_sink", "co2_escape_sink",
    ]
    ordered = [x for x in preferred if x in ids] + sorted(ids - set(preferred))
    return [Node(x, pretty_label(x), node_group(x), node_kind(x, state_names)) for x in ordered]


# ---------------------------------------------------------------------------
# Writers
# ---------------------------------------------------------------------------

CLUSTER_LABELS = {
    "substrate": "Substrate pools",
    "environment": "Environment",
    "enzyme": "Enzymatic module",
    "yeast": "Yeast physiology / fluxes",
    "products": "Products and feedback",
    "other": "Other",
}


def dot_escape(text: str) -> str:
    return text.replace("\\", "\\\\").replace('"', '\\"')


def node_style(node: Node) -> str:
    if node.kind == "state":
        return 'shape=ellipse, style="rounded,filled", fillcolor="white", penwidth=1.4'
    if node.kind == "flux":
        return 'shape=diamond, style="filled", fillcolor="gray92", penwidth=1.2'
    if node.kind == "source":
        return 'shape=house, style="filled", fillcolor="gray95", penwidth=1.1'
    if node.kind == "sink":
        return 'shape=invhouse, style="filled", fillcolor="gray95", penwidth=1.1'
    return 'shape=box, style="rounded,filled", fillcolor="gray96", penwidth=1.1'


def edge_style(edge: Edge) -> str:
    # Grayscale/pattern-based so it survives journal printing.
    attrs = []
    if edge.sign == "-":
        attrs.extend(['style="dashed"', 'arrowhead="tee"'])
    elif edge.sign == "mod":
        attrs.extend(['style="dotted"', 'arrowhead="normal"'])
    else:
        attrs.extend(['style="solid"', 'arrowhead="normal"'])
    if edge.relation in {"aggregate"}:
        attrs.append('style="dotted"')
    attrs.append('penwidth=1.1')
    attrs.append(f'label="{dot_escape(edge.label)}"')
    return ", ".join(attrs)


def write_dot(path: Path, nodes: List[Node], edges: List[Edge], title: str) -> None:
    by_group: Dict[str, List[Node]] = {}
    for n in nodes:
        by_group.setdefault(n.group, []).append(n)

    lines: List[str] = []
    lines.append("digraph ODE_Network {")
    lines.append("  graph [rankdir=LR, splines=true, overlap=false, compound=true, fontsize=18, fontname=Helvetica, labelloc=t, label=\"%s\"];" % dot_escape(title))
    lines.append("  node [fontname=Helvetica, fontsize=11, margin=0.08];")
    lines.append("  edge [fontname=Helvetica, fontsize=8, arrowsize=0.7];")

    cluster_order = ["substrate", "enzyme", "environment", "yeast", "products", "other"]
    for idx, group in enumerate(cluster_order):
        ns = by_group.get(group, [])
        if not ns:
            continue
        lines.append(f"  subgraph cluster_{group} {{")
        lines.append(f"    label=\"{CLUSTER_LABELS.get(group, group)}\";")
        lines.append("    style=\"rounded\";")
        lines.append("    color=\"gray70\";")
        lines.append("    penwidth=1.0;")
        for n in ns:
            lines.append(f"    \"{dot_escape(n.node_id)}\" [label=\"{dot_escape(n.label)}\", {node_style(n)}];")
        lines.append("  }")

    for edge in edges:
        lines.append(f"  \"{dot_escape(edge.source)}\" -> \"{dot_escape(edge.target)}\" [{edge_style(edge)}];")
    lines.append("}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_csvs(outdir: Path, prefix: str, nodes: List[Node], edges: List[Edge]) -> None:
    with (outdir / f"{prefix}_nodes.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["node_id", "label", "group", "kind"])
        w.writeheader()
        for n in nodes:
            w.writerow({"node_id": n.node_id, "label": n.label, "group": n.group, "kind": n.kind})
    with (outdir / f"{prefix}_edges.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["source", "target", "label", "relation", "sign", "module"])
        w.writeheader()
        for e in edges:
            w.writerow({
                "source": e.source,
                "target": e.target,
                "label": e.label,
                "relation": e.relation,
                "sign": e.sign,
                "module": e.module,
            })


def render_dot(dot_path: Path, svg: bool = True, png: bool = True) -> None:
    dot = shutil.which("dot")
    if not dot:
        print(f"[warn] Graphviz 'dot' executable not found; wrote DOT only: {dot_path}")
        return
    if svg:
        subprocess.run([dot, "-Tsvg", str(dot_path), "-o", str(dot_path.with_suffix(".svg"))], check=True)
    if png:
        subprocess.run([dot, "-Tpng", str(dot_path), "-o", str(dot_path.with_suffix(".png"))], check=True)


def validate_state_coverage(nodes: List[Node], state_names: Sequence[str]) -> Tuple[List[str], List[str]]:
    node_ids = {n.node_id for n in nodes}
    missing = [s for s in state_names if s not in node_ids]
    extra_states = [n.node_id for n in nodes if n.kind == "state" and n.node_id not in state_names]
    return missing, extra_states


def build_and_write(mode: str, outdir: Path, state_names: Sequence[str]) -> None:
    if mode == "core":
        edges = build_core_edges()
        title = "Core ODE-Network Structure of the Fermentation Simulator"
    elif mode == "full":
        edges = build_full_edges()
        title = "Full ODE-Network Structure of the Fermentation Simulator"
    else:
        raise ValueError(mode)

    nodes = build_nodes(edges, state_names)
    prefix = f"ode_network_{mode}"
    write_csvs(outdir, prefix, nodes, edges)
    dot_path = outdir / f"{prefix}.dot"
    write_dot(dot_path, nodes, edges, title)
    render_dot(dot_path)

    missing, extra = validate_state_coverage(nodes, state_names)
    if mode == "full":
        if missing:
            print("[warn] Dynamic STATE_NAMES absent from full graph:", ", ".join(missing))
        if extra:
            print("[warn] Graph states not in imported STATE_NAMES:", ", ".join(extra))
    print(f"[ok] wrote {prefix}: {len(nodes)} nodes, {len(edges)} edges -> {outdir}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Visualize the implicit ODE network in fermentation_benchmark_full.py")
    parser.add_argument("--outdir", default="ode_network_visualization", help="Output directory")
    parser.add_argument("--mode", choices=["core", "full", "both"], default="both", help="Which network to draw")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    state_names = load_state_names()
    modes = ["core", "full"] if args.mode == "both" else [args.mode]
    for mode in modes:
        build_and_write(mode, outdir, state_names)

    readme = outdir / "README.txt"
    readme.write_text(
        "ODE-network visualization outputs\n"
        "=================================\n\n"
        "The ODE-network is implicit in fermentation_benchmark_full.py: STATE_NAMES defines dynamic nodes, "
        "rate_factors()/stress_index() define regulatory factors, rhs() defines fluxes, and euler_step() performs synchronous updates.\n\n"
        "Recommended manuscript figure: ode_network_core.svg or ode_network_core.png.\n"
        "Recommended supplementary figure: ode_network_full.svg or ode_network_full.png.\n"
        "Edge-list CSVs are included for transparent reporting.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
