#!/usr/bin/env python3
"""Render the PFCandidate truth validation from a harvested DQM file.

Per species and region: the cut-flow ladder of the track-found and track-missing
branches, the unconditional rung efficiencies versus pt and eta, the energy response per
branch and the foreign-energy share. Per PF type and region: the candidate class
fractions versus pt and eta. Species, regions and types are discovered from the file.

  makePFCandidateTruthPlots.py DQM_V0001_*.root --outputDir plots --sample "electron gun"
"""
import argparse, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning
plt.style.use(hep.style.CMS)

ap = argparse.ArgumentParser()
ap.add_argument("inFile")
ap.add_argument("--outputDir", default="pfcandidate_plots")
ap.add_argument("--sample", default="")
ap.add_argument("--minEntries", type=int, default=10, help="skip a species/region with fewer particles")
args = ap.parse_args()
os.makedirs(args.outputDir, exist_ok=True)

f = ROOT.TFile.Open(args.inFile)
def find_base(d, path="DQMData"):
    top = f.Get(path)
    for k in top.GetListOfKeys():
        if k.GetName().startswith("Run "):
            return path + "/" + k.GetName() + "/TruthInfo/Run summary/Offline/PFCandidates"
    raise SystemExit("no run folder in " + args.inFile)
base = find_base(f)
collections = [k.GetName() for k in f.Get(base).GetListOfKeys()]

RUNGS = ["trackFound", "trackInCandidate", "ecalCollected", "ecalLinked", "hcalCollected", "hcalLinked",
         "hgcalCollected", "hgcalLinked", "caloSameCandidate", "candidateFound", "merged", "partial", "clean", "pdgCorrect"]
CLASSES = ["matched", "merged", "split", "unmatched", "other", "unevaluated"]
COLORS = plt.cm.tab20.colors

def arrays(h):
    n = h.GetNbinsX()
    edges = np.array([h.GetXaxis().GetBinLowEdge(i) for i in range(1, n + 2)])
    vals = np.array([h.GetBinContent(i) for i in range(1, n + 1)])
    errs = np.array([h.GetBinError(i) for i in range(1, n + 1)])
    return edges, vals, errs

def label(ax, title):
    hep.cms.label("Private Work", data=False, rlabel=args.sample, ax=ax, fontsize=13)
    ax.text(0.0, 1.10, title, transform=ax.transAxes, fontsize=14, ha="left", va="bottom")

pages = []
def save(fig, name, section):
    fig.savefig(os.path.join(args.outputDir, name), dpi=100, bbox_inches="tight"); plt.close(fig)
    pages.append((section, name))

for coll in collections:
    truthDir = f.Get("%s/%s/truth" % (base, coll))
    species = [k.GetName() for k in truthDir.GetListOfKeys()]
    for s in species:
        for region in [k.GetName() for k in truthDir.Get(s).GetListOfKeys()]:
            d = truthDir.Get(s + "/" + region)
            pop = d.Get("num_simul_pt")
            if not pop or pop.GetEntries() < args.minEntries:
                continue
            section = "%s / %s / %s" % (coll, s, region)
            # ladders
            fig, axes = plt.subplots(1, 2, figsize=(18, 7))
            for ax, name, ttl in ((axes[0], "ladder_trackFound", "track found branch"), (axes[1], "ladder_trackMissing", "track missing / neutral branch")):
                h = d.Get(name)
                labels = [h.GetXaxis().GetBinLabel(b) for b in range(1, h.GetNbinsX() + 1)]
                vals = [h.GetBinContent(b) for b in range(1, h.GetNbinsX() + 1)]
                popn = vals[0] if vals and vals[0] > 0 else 1.
                ax.bar(range(len(vals)), [v / popn for v in vals], color="tab:blue")
                for i, v in enumerate(vals):
                    ax.text(i, v / popn + 0.01, "%d" % v, ha="center", fontsize=10)
                ax.set_xticks(range(len(labels))); ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=11)
                ax.set_ylim(0, 1.15); ax.set_ylabel("fraction of the branch population")
                label(ax, "%s, %s: %s" % (s, region, ttl))
            save(fig, "%s_%s_%s_ladder.png" % (coll, s, region), section)
            # rung efficiencies vs pt and eta
            for axis, xlabel, logx in (("pt", r"$p_T$ [GeV]", True), ("eta", r"$\eta$", False), ("caloeta", r"$\eta$ at calorimeter entrance", False)):
                fig, ax = plt.subplots(figsize=(11, 8))
                drawn = 0
                for i, rung in enumerate(RUNGS):
                    h = d.Get("eff_%s_vs_%s" % (rung, axis))
                    den = d.Get("num_%s_expected_%s" % (rung, axis))
                    if not h or not den or den.GetEntries() == 0:
                        continue
                    edges, vals, errs = arrays(h)
                    centers = 0.5 * (edges[1:] + edges[:-1])
                    mask = np.array([den.GetBinContent(b) > 0 for b in range(1, den.GetNbinsX() + 1)])
                    ax.errorbar(centers[mask], vals[mask], yerr=errs[mask], fmt="o-", ms=4, lw=1.2, color=COLORS[i % 20], label="%s (%d)" % (rung, den.GetEntries()))
                    drawn += 1
                if drawn == 0:
                    plt.close(fig); continue
                if logx: ax.set_xscale("log")
                ax.set_ylim(0, 1.15); ax.set_xlabel(xlabel); ax.set_ylabel("efficiency (each rung over the particles it applies to)")
                ax.legend(fontsize=10, ncol=2, loc="upper left", bbox_to_anchor=(1.01, 1.0), title="rung (denominator)")
                label(ax, "%s, %s: rung efficiencies" % (s, region))
                save(fig, "%s_%s_%s_rungs_%s.png" % (coll, s, region, axis), section)
            # response and foreign fraction
            fig, axes = plt.subplots(1, 2, figsize=(18, 7))
            for name, style in (("response_trackFound", "-"), ("response_trackMissing", "--")):
                h = d.Get(name)
                if h and h.GetEntries() > 0:
                    e, v, _ = arrays(h)
                    hep.histplot(v, e, ax=axes[0], histtype="step", linestyle=style, linewidth=1.8, label="%s (%d, mean %.2f)" % (name.replace("response_", ""), h.GetEntries(), h.GetMean()))
            axes[0].set_xlabel(r"$E_{cand}/E_{truth}$"); axes[0].set_ylabel("particles"); axes[0].legend(fontsize=11)
            label(axes[0], "%s, %s: energy response" % (s, region))
            h = d.Get("foreign_fraction")
            if h and h.GetEntries() > 0:
                e, v, _ = arrays(h)
                hep.histplot(v, e, ax=axes[1], histtype="step", linewidth=1.8, color="tab:red")
            axes[1].set_xlabel("foreign energy share of the candidate"); axes[1].set_ylabel("particles")
            label(axes[1], "%s, %s: contamination" % (s, region))
            save(fig, "%s_%s_%s_response.png" % (coll, s, region), section)

    recoDir = f.Get("%s/%s/reco" % (base, coll))
    for t in [k.GetName() for k in recoDir.GetListOfKeys()]:
        for region in [k.GetName() for k in recoDir.Get(t).GetListOfKeys()]:
            d = recoDir.Get(t + "/" + region)
            n = d.Get("num_reco_pt")
            if not n or n.GetEntries() < args.minEntries:
                continue
            section = "%s / candidates %s / %s" % (coll, t, region)
            for axis, xlabel, logx in (("pt", r"$p_T$ [GeV]", True), ("eta", r"$\eta$", False)):
                fig, ax = plt.subplots(figsize=(11, 8))
                for i, cls in enumerate(CLASSES):
                    h = d.Get("frac_%s_vs_%s" % (cls, axis)); num = d.Get("num_%s_%s" % (cls, axis))
                    if not h or not num or num.GetEntries() == 0:
                        continue
                    edges, vals, errs = arrays(h)
                    centers = 0.5 * (edges[1:] + edges[:-1])
                    mask = np.array([n.GetBinContent(b) > 0 for b in range(1, n.GetNbinsX() + 1)]) if axis == "pt" else np.array([d.Get("num_reco_eta").GetBinContent(b) > 0 for b in range(1, h.GetNbinsX() + 1)])
                    ax.errorbar(centers[mask], vals[mask], yerr=errs[mask], fmt="o-", ms=4, lw=1.2, color=COLORS[(2 * i) % 20], label="%s (%d)" % (cls, num.GetEntries()))
                if logx: ax.set_xscale("log")
                ax.set_ylim(0, 1.15); ax.set_xlabel(xlabel); ax.set_ylabel("fraction of candidates")
                ax.legend(fontsize=11, loc="upper left", bbox_to_anchor=(1.01, 1.0), title="class (candidates)")
                label(ax, "PF type %s, %s: candidate classes (%d)" % (t, region, n.GetEntries()))
                save(fig, "%s_reco_%s_%s_classes_%s.png" % (coll, t, region, axis), section)

with open(os.path.join(args.outputDir, "index.html"), "w") as html:
    html.write("<html><body style='font-family:sans-serif'><h1>PFCandidate truth validation</h1><p>%s</p>\n" % args.sample)
    current = None
    for section, name in pages:
        if section != current:
            html.write("<h2>%s</h2>\n" % section); current = section
        html.write("<a href='%s'><img src='%s' width='480'></a>\n" % (name, name))
    html.write("</body></html>\n")
print("wrote", len(pages), "plots to", args.outputDir)
