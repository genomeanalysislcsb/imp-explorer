"""Synthetic IMP3 output tree in the exact file formats of IMP3 3.0 (taken from a real sample), for tests and demos.

python3 tests/make_fixture.py OUTDIR [n_samples]  ->  OUTDIR/imp/<sample>/..., OUTDIR/metadata.tsv
"""
import gzip
import os
import random
import sys

GENERA = [("Bacteroides", "Bacteroidota", "Bacteroidia", "Bacteroidales", "Bacteroidaceae"),
          ("Faecalibacterium", "Bacillota_A", "Clostridia", "Oscillospirales", "Ruminococcaceae"),
          ("Blautia_A", "Bacillota_A", "Clostridia", "Lachnospirales", "Lachnospiraceae"),
          ("Akkermansia", "Verrucomicrobiota", "Verrucomicrobiae", "Verrucomicrobiales", "Akkermansiaceae"),
          ("Alistipes", "Bacteroidota", "Bacteroidia", "Bacteroidales", "Rikenellaceae"),
          ("Prevotella", "Bacteroidota", "Bacteroidia", "Bacteroidales", "Bacteroidaceae"),
          ("Ruminococcus_E", "Bacillota_A", "Clostridia", "Oscillospirales", "Acutalibacteraceae"),
          ("Methanobrevibacter_A", "Methanobacteriota", "Methanobacteria", "Methanobacteriales",
           "Methanobacteriaceae"),
          ("Escherichia", "Pseudomonadota", "Gammaproteobacteria", "Enterobacterales", "Enterobacteriaceae"),
          ("Bifidobacterium", "Actinomycetota", "Actinomycetes", "Actinomycetales", "Bifidobacteriaceae")]
KOS = ["K%05d" % i for i in (1803, 2111, 2112, 3043, 3046, 1689, 873, 1915, 2335, 3086, 1810, 2986, 2992, 878,
                             1006, 4043, 3406, 2470, 6147, 1897)]
MODULES = ["M%05d" % i for i in (1, 2, 3, 9, 11, 21, 26, 50, 63, 82)]


def w(path, content):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "wt") as fh:
        fh.write(content)


def fc_table(omic, g, counts):
    hdr = ('# Program:featureCounts v2.0.6; Command:"featureCounts" "-p" "-O" "--fraction" "-M" "--fraction" "-t" '
           '"CDS" "-g" "%s"\n' % g)
    hdr += "Geneid\tChr\tStart\tEnd\tStrand\tLength\tAssembly/%s.reads.paired.sorted.bam\tAssembly/%s.reads.single.sorted.bam\n" % (omic, omic)
    return hdr + "".join("%s\tc1;c2\t1;5\t9;20\t+;+\t900\t%.2f\t%.2f\n" % (k, v * 0.99, v * 0.01)
                         for k, v in counts.items())


def kreport(rows):
    return "".join("%.2f\t%d\t%d\t%s\t%s\t%s\n" % r for r in rows)


def make_sample(root, sid, group, rng, old_bracken, folder=None):
    b = os.path.join(root, folder or sid)
    shift = {"control": 1.0, "case_a": 1.8, "case_b": 0.6}[group]
    w(os.path.join(b, "status", "all.done"), "")
    # reads
    for o, n in (("mg", rng.randint(25, 40) * 10 ** 6), ("mt", rng.randint(35, 55) * 10 ** 6)):
        t = int(n * 0.99)
        f = int(t * (0.97 if o == "mg" else 0.6))
        w(os.path.join(b, "Stats", o, "%s.read_stats.tsv.gz" % o),
          "file\tformat\ttype\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len\tQ1\tQ2\tQ3\tsum_gap\tN50\tN50_num\tQ20(%%)\tQ30(%%)\tAvgQual\tGC(%%)\n"
          "%s.r1.fq\tFASTQ\tDNA\t%d\t0\t35\t148.7\t151\t150\t151\t151\t0\t151\t1\t98.1\t95.1\t27.5\t50.0\n"
          "%s.r1.trimmed.fq\tFASTQ\tDNA\t%d\t0\t40\t147.6\t150\t149\t150\t150\t0\t150\t1\t98.2\t%.2f\t27.8\t49.9\n"
          % (o, n, o, t, 95 + rng.random()))
        for rr in ("r1", "r2"):
            w(os.path.join(b, "Stats", o, "%s.%s_nonpareil.diversity.txt.gz" % (o, rr)),
              "model_a\t81.1\nmodel_b\t4.3\nnopareil_diversity\t%.3f\n" % (18 + rng.random() * 2))
        w(os.path.join(b, "Stats", o, "mgmt.assembly.contig_flagstat.txt.gz"),
          "%d\n%d\n%d\n0\n0\n0\n%d\n%d\n" % (2 * f * 1.01, 2 * f, 2 * f * 0.01, 2 * f * 0.96, 2 * f * 0.95))
    w(os.path.join(b, "Stats", "all_stats.tsv.gz"), "type\tstep\tnumber\n" + "".join(
        "%s\t%s.r1.%s\t%d\n" % (o, o, s, v) for o in ("mg", "mt")
        for s, v in (("fq", 3e7), ("trimmed.fq", 2.97e7), ("preprocessed.fq", 2.9e7 if o == "mg" else 1.8e7))))
    # assembly + bins
    contigs, c2b, bins = [], [], []
    for i in range(1, 401):
        contigs.append(("S%s_merged_C%010d" % (sid.lower(), i), rng.randint(500, 60000)))
    w(os.path.join(b, "Stats", "mgmt", "mgmt.assembly.length.txt.gz"), "".join("%s\t%d\n" % c for c in contigs))
    nb = rng.randint(6, 9)
    gen = rng.sample(GENERA, nb)
    for k, g in enumerate(gen):
        bins.append(("bin.%d" % k, g))
        for c in contigs[k * 30:(k + 1) * 30]:
            c2b.append((c[0], "bin.%d" % k))
    w(os.path.join(b, "Binning", "MetaWrap", "bins.contigs"), "".join("%s\t%s\n" % x for x in c2b))
    head = ["Bin", "Completeness", "Contamination", "Contig_N50", "Genome_Size", "GC_Content", "classification",
            "closest_genome_reference", "closest_genome_reference_radius", "closest_genome_taxonomy",
            "closest_genome_ani", "closest_genome_af"]
    lines = ["\t".join("%-10s " % h for h in head)]
    for bn, g in bins:
        cl = "d__%s;p__%s;c__%s;o__%s;f__%s;g__%s;s__%s sp%d" % (
            "Archaea" if g[0].startswith("Methano") else "Bacteria", g[1], g[2], g[3], g[4], g[0], g[0],
            rng.randint(1, 99))
        row = [bn, "%.2f" % rng.uniform(55, 100), "%.2f" % rng.uniform(0, 9), str(rng.randint(5000, 200000)),
               str(rng.randint(1500000, 4500000)), "%.2f" % rng.uniform(0.3, 0.6), cl, "GCF_%09d.1" % rng.randint(1, 10 ** 8),
               "95.0", cl, "%.2f" % rng.uniform(95, 99.9), "0.9"]
        lines.append("\t".join("%10s " % x for x in row))
    w(os.path.join(b, "Stats", "all_bin_stats.tsv"), "\n".join(lines) + "\n")
    for o in ("mg", "mt"):
        rc, dp = [], []
        for name, L in contigs:
            bi = next((x[1] for x in c2b if x[0] == name), None)
            fac = 1.0
            if bi is not None:
                g = gen[int(bi.split(".")[1])][0]
                fac = (shift if g in ("Akkermansia", "Methanobrevibacter_A") else 1.0) * (3 if o == "mt" and g == "Faecalibacterium" else 1)
            cnt = int(L / 150 * rng.uniform(1, 12) * fac)
            rc.append("%s\t0\t%d\t%d\n" % (name, L, cnt))
            dp.append("%s\t%.5f\n" % (name, cnt * 150 / L))
        w(os.path.join(b, "Stats", o, "mgmt.assembly.contig_read_count.txt.gz"), "".join(rc))
        w(os.path.join(b, "Stats", o, "mgmt.assembly.contig_depth.txt.gz"), "".join(dp))
    # taxonomy
    for o in ("mg", "mt"):
        ab = {g[0]: rng.lognormvariate(10, 1.2) * (shift if g[0] in ("Akkermansia", "Methanobrevibacter_A") else 1)
              for g in GENERA}
        tot = sum(ab.values())
        rk = {"G": [], "S": [], "P": []}
        phy = {}
        for g in GENERA:
            phy[g[1]] = phy.get(g[1], 0) + ab[g[0]]
        for i, g in enumerate(GENERA):
            v = int(ab[g[0]])
            rk["G"].append((100 * v / tot, v, v, "G", str(1000 + i), "            " + g[0]))
            rk["S"].append((100 * v / tot, v, v, "S", str(2000 + i), "              %s sp1" % g[0]))
        for i, (p, v) in enumerate(phy.items()):
            rk["P"].append((100 * v / tot, v, v, "P", str(3000 + i), "      " + p))
        for r, rows in rk.items():
            lvl = {"G": "genuses", "S": "species", "P": "phylums"}[r]
            if old_bracken:
                w(os.path.join(b, "Analysis", "taxonomy", "kraken", "%s.kraken_bracken_%s.report.gz" % (o, lvl)),
                  kreport([(0.5, 100, 100, "U", "0", "unclassified")] + rows))
            else:
                w(os.path.join(b, "Analysis", "taxonomy", "bracken", o, "%s.bracken" % r),
                  "name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\tnew_est_reads\tfraction_total_reads\n" +
                  "".join("%s\t%s\t%s\t%d\t0\t%d\t%.5f\n" % (x[5].strip(), x[4], r, x[1], x[1], x[0] / 100) for x in rows))
        if not old_bracken:                       # hierarchical mOTUs output (motus profile -p -k all)
            lines = ["# motus\n#consensus_taxonomy\t%s\n" % sid]
            for i, g in enumerate(GENERA):
                v = int(ab[g[0]] / 50)
                lines.append("k__Bacteria|p__%s|c__%s|o__%s|f__%s|g__%s\t%d\n" % (g[1], g[2], g[3], g[4], g[0], v))
                lines.append("k__Bacteria|p__%s|c__%s|o__%s|f__%s|g__%s|s__%s sp [ref_mOTU_v4_%05d]\t%d\n" % (
                    g[1], g[2], g[3], g[4], g[0], g[0], i, v))
            lines.append("unassigned\t99\n")
            w(os.path.join(b, "Analysis", "taxonomy", "mOTUs", "%s.mOTU.counts.tsv.gz" % o), "".join(lines))
        else:
            w(os.path.join(b, "Analysis", "taxonomy", "mOTUs", "%s.mOTU.counts.tsv.gz" % o),
              "# git tag version 3.1.0\n# call: motus profile\n#consensus_taxonomy\t%s\n" % sid +
              "".join("%s sp [ref_mOTU_v31_%05d]\t%d\n" % (g[0], i, int(ab[g[0]] / 50))
                      for i, g in enumerate(GENERA)) + "unassigned\t1234\n")
        # functions
        kos = {k: rng.lognormvariate(8, 1) * (2.5 if o == "mt" and k == "K01803" else 1) for k in KOS}
        if group == "case_a":
            kos["K00873"] *= 2
        w(os.path.join(b, "Analysis", "annotation", "%s.mantis.kegg_ko.counts.tsv.gz" % o), fc_table(o, "kegg_ko", kos))
        w(os.path.join(b, "Analysis", "annotation", "%s.mantis.kegg_module.counts.tsv.gz" % o),
          fc_table(o, "kegg_module", {m: rng.lognormvariate(9, 1) for m in MODULES}))


GROUPS = ("control", "case_a", "case_b")


def main(out, n=12, seed=7, folder_fmt="{sid}"):
    """n samples S01 ... in OUTDIR/imp/<folder_fmt>/, metadata in OUTDIR/metadata.tsv (sample_id, group, site, age,
    sex). folder_fmt formats the folder name from the sample id, e.g. "run_{sid_lower}" to test id translation."""
    rng = random.Random(seed)
    root = os.path.join(out, "imp")
    md = ["sample_id\tgroup\tsite\tage\tsex"]
    for i in range(n):
        group = GROUPS[i % 3]
        sid = "S%02d" % (i + 1)
        folder = folder_fmt.format(sid=sid, sid_lower=sid.lower(), i=i + 1)
        make_sample(root, sid, group, rng, old_bracken=(i % 2 == 0), folder=folder)
        md.append("%s\t%s\t%s\t%d\t%s" % (sid, group, "site_%d" % (1 + (i // 3) % 2), rng.randint(45, 75),
                                          rng.choice("FM")))
    w(os.path.join(out, "metadata.tsv"), "\n".join(md) + "\n")
    return root


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 12)
