import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

data = pd.read_csv("UMICH_SIGMAALDRICH_01172023.txt", sep="\t")

m = []
sms = []
for k in data["SMILES"]:
    try:
        sm = Chem.CanonSmiles(k)
        if sm not in sms:
            sms.append(sm)
            m.append(Chem.MolFromSmiles(sm))
    except:
        continue

acid_sp2 = Chem.MolFromSmarts("c[CX3](=O)[OX2H1]")
acid_sp3 = Chem.MolFromSmarts("C[CX3](=O)[OX2H1]")

aldehyde_sp2 = Chem.MolFromSmarts("[CX3H1](=O)[c]")
aldehyde_sp3 = Chem.MolFromSmarts("[CX3H1](=O)[#6!c]")

alcohol_sp2 = Chem.MolFromSmarts("[c][OX2H]")
alcohol_sp3 = Chem.MolFromSmarts("[OX2H][CX4!c;!$([CX4!c]([OX2H])[O,S,#7,#15])]")

amine_sp2 = Chem.MolFromSmarts("[NX3;H2;!$(NC=[!#6]);!$(NC#[!#6])][c]")
amine_sp3 = Chem.MolFromSmarts("[NX3;H2;!$(NC=[!#6]);!$(NC#[!#6])][#6!c]")

bromide_sp2 = Chem.MolFromSmarts("[Br][c]")
bromide_sp3 = Chem.MolFromSmarts("[Br][#6!c]")

iodide_sp2 = Chem.MolFromSmarts("[I][c]")
iodide_sp3 = Chem.MolFromSmarts("[I][#6!c]")

chloride_sp2 = Chem.MolFromSmarts("[Cl][c]")
chloride_sp3 = Chem.MolFromSmarts("[Cl][#6!c]")

boronate_sp2 = Chem.MolFromSmarts("[O]B(c)[O]")
boronate_sp3 = Chem.MolFromSmarts("[O]B([#6!c])[O]")

labels = ["acid_sp2", "acid_sp3", "aldehyde_sp2", "aldehyde_sp3", "alcohol_sp2", \
          "alcohol_sp3", "amine_sp2", "amine_sp3", "bromide_sp2", "bromide_sp3", \
          "iodide_sp2", "iodide_sp3", "chloride_sp2", "chloride_sp3", "boronate_sp2", \
          "boronate_sp3"]



ss_smarts = [acid_sp2, acid_sp3, aldehyde_sp2, aldehyde_sp3, alcohol_sp2, \
             alcohol_sp3, amine_sp2, amine_sp3, bromide_sp2, bromide_sp3, \
             iodide_sp2, iodide_sp3, chloride_sp2, chloride_sp3, boronate_sp2, \
             boronate_sp3]

bbs = {l:[] for l in labels}
bbs_counts = {label: 0 for label in labels}

for k in m:
    mol = k
    if not mol:
        continue
    for ss, smarts in zip(labels, ss_smarts):
        hits = mol.GetSubstructMatches(smarts)
        if hits:
            bbs[ss].append(mol)
            bbs_counts[ss] += len(hits)


for bb in bbs:
    print(bb, len(bbs[bb]))

d_acid_sp2 = AllChem.ReactionFromSmarts("[c:1][CX3](=O)[OX2H1]>>[2H][c:1]")
d_acid_sp3 = AllChem.ReactionFromSmarts("[C:1][CX3](=O)[OX2H1]>>[C:1][2H]")

d_aldehyde_sp2 = AllChem.ReactionFromSmarts("[CX3H1](=O)[c:1]>>[2H][c:1]")
d_aldehyde_sp3 = AllChem.ReactionFromSmarts("[CX3H1](=O)[#6!c:1]>>[2H][#6!c:1]")

d_alcohol_sp2 = AllChem.ReactionFromSmarts("[c:1][OX2H]>>[c:1][2H]")
d_alcohol_sp3 = AllChem.ReactionFromSmarts("[OX2H][CX4!c;!$([CX4!c]([OX2H])[O,S,#7,#15]):1]>>[CX4!c:1][2H]")

d_amine_sp2 = AllChem.ReactionFromSmarts("[NX3;H2;!$(NC=[!#6]);!$(NC#[!#6])][c:1]>>[2H][c:1]")
d_amine_sp3 = AllChem.ReactionFromSmarts("[NX3;H2;!$(NC=[!#6]);!$(NC#[!#6])][#6!c:1]>>[2H][#6!c:1]")

d_bromide_sp2 = AllChem.ReactionFromSmarts("[Br][c:1]>>[2H][c:1]")
d_bromide_sp3 = AllChem.ReactionFromSmarts("[Br][#6!c:1]>>[2H][#6!c:1]")

d_iodide_sp2 = AllChem.ReactionFromSmarts("[I][c:1]>>[2H][c:1]")
d_iodide_sp3 = AllChem.ReactionFromSmarts("[I][#6!c:1]>>[2H][c:1]")

d_chloride_sp2 = AllChem.ReactionFromSmarts("[Cl][c:1]>>[2H][c:1]")
d_chloride_sp3 = AllChem.ReactionFromSmarts("[Cl][#6!c:1]>>[2H][c:1]")

d_boronate_sp2 = AllChem.ReactionFromSmarts("[O]B([c:1])[O]>>[2H]([c:1])")
d_boronate_sp3 = AllChem.ReactionFromSmarts("[O]B([#6!c:1])[O]>>[2H]([c:1])")


def deut(sm, clas):
    try:
        if clas == "acid_sp2":
            j = d_acid_sp2.RunReactants((sm,))
        if clas == "acid_sp3":
            j = d_acid_sp3.RunReactants((sm,))
        if clas == "aldehyde_sp2":
            j = d_aldehyde_sp2.RunReactants((sm,))
        if clas == "aldehyde_sp3":
            j = d_aldehyde_sp3.RunReactants((sm,))
        if clas == "alcohol_sp2":
            j = d_alcohol_sp2.RunReactants((sm,))
        if clas == "alcohol_sp3":
            j = d_alcohol_sp3.RunReactants((sm,))
        if clas == "amine_sp2":
            j = d_amine_sp2.RunReactants((sm,))
        if clas == "amine_sp3":
            j = d_amine_sp3.RunReactants((sm,))
        if clas == "bromide_sp2":
            j = d_bromide_sp2.RunReactants((sm,))
        if clas == "bromide_sp3":
            j = d_bromide_sp3.RunReactants((sm,))
        if clas == "iodide_sp2":
            j = d_iodide_sp2.RunReactants((sm,))
        if clas == "iodide_sp3":
            j = d_iodide_sp3.RunReactants((sm,))
        if clas == "chloride_sp2":
            j = d_chloride_sp2.RunReactants((sm,))
        if clas == "chloride_sp3":
            j = d_chloride_sp3.RunReactants((sm,))
        if clas == "boronate_sp2":
            j = d_boronate_sp2.RunReactants((sm,))
        if clas == "boronate_sp3":
            j = d_boronate_sp3.RunReactants((sm,))
        return [Chem.MolToSmiles(k[0]) for k in j]
    except:
        print(Chem.MolToSmiles(sm), clas)

bbs_deut = {l:[] for l in labels}
bbs_deut_counts = {label: 0 for label in labels}


for j in bbs:
    print(j)
    print(len(bbs[j]))
    for k in bbs[j]:
        blocks = deut(k,j)
        bbs_deut_counts[j] += len(list(set(blocks)))
        for i in blocks:
            #if j == "aniline":
                #print(i)

            if i not in bbs_deut[j]:
                bbs_deut[j].append(i)

    print(len(bbs_deut[j]))

print(bbs_deut_counts)

labels_hal = ["acid_sp2",    "acid_sp3",   "aldehyde_sp2", "aldehyde_sp3", "alcohol_sp2",\
              "alcohol_sp3", "amine_sp2",  "amine_sp3",    "halide_sp2",  "halide_sp3",  \
              "boronate_sp2","boronate_sp3"]

hal_sp2 = []
hal_sp2.extend(bbs_deut["bromide_sp2"])
hal_sp2.extend(bbs_deut["iodide_sp2"])
hal_sp2.extend(bbs_deut["chloride_sp2"])
bbs_deut["halide_sp2"] = hal_sp2
hal_sp3 = []
hal_sp3.extend(bbs_deut["bromide_sp3"])
hal_sp3.extend(bbs_deut["iodide_sp3"])
hal_sp3.extend(bbs_deut["chloride_sp3"])
bbs_deut["halide_sp3"] = hal_sp3

mhal_sp2 = []
mhal_sp2.extend(bbs["bromide_sp2"])
mhal_sp2.extend(bbs["iodide_sp2"])
mhal_sp2.extend(bbs["chloride_sp2"])
bbs["halide_sp2"] = mhal_sp2
mhal_sp3 = []
mhal_sp3.extend(bbs["bromide_sp3"])
mhal_sp3.extend(bbs["iodide_sp3"])
mhal_sp3.extend(bbs["chloride_sp3"])
bbs["halide_sp3"] = mhal_sp3

bbs_deut_counts["halide_sp2"] = bbs_deut_counts["bromide_sp2"] + bbs_deut_counts["iodide_sp2"] + bbs_deut_counts["chloride_sp2"]
bbs_deut_counts["halide_sp3"] = bbs_deut_counts["bromide_sp3"] + bbs_deut_counts["iodide_sp3"] + bbs_deut_counts["chloride_sp3"]

for l in labels_hal:
    if "sp2" not in l:
        continue
    for l2 in labels_hal:
        if "sp2" not in l2:
            continue
        if l.split("_")[0] == l2.split("_")[0]:
            continue
        if l.split("_")[1] != l2.split("_")[1]:
            continue
        if "bromide" in l or "iodide" in l or "chloride" in l or "bromide" in l2 or "iodide" in l2 or "chloride" in l2:
            continue

from matplotlib import pyplot as plt
from matplotlib_venn import venn2, venn2_circles
import matplotlib
import math

font = {'family' : 'arial',
        'weight' : 'bold',
        'size'   : 5}
matplotlib.rc('font', **font)

n = 0
for k in labels_hal:
    if "sp2" in k:
        n = n + 1
x = 0
ystart = 1
fig,ax = plt.subplots(n,n,figsize=(4.25,3.5))

seen = []

for l in labels_hal:
    if "sp3" not in l:
        continue
    i = 0
    y = ystart
    for l2 in labels_hal:
        if l2 in seen:
            continue
        if y >= n:
            break
        if "sp3" not in l2:
            continue
        if l.split("_")[0] == l2.split("_")[0]:
            continue
        if l.split("_")[1] != l2.split("_")[1]:
            continue
        if "bromide" in l or "iodide" in l or "chloride" in l or "bromide" in l2 or "iodide" in l2 or "chloride" in l2:
            continue
        # print(x,y)
        # print(i, l, l2)
        samezies = 0
        for a1 in bbs_deut[l]:
            for a2 in bbs_deut[l2]:
                if a1 == a2:
                    samezies = samezies + 1
        v = venn2(subsets = (bbs_deut_counts[l], bbs_deut_counts[l2], samezies), set_labels = ("", "", ""), ax=ax[x][y])
        if len(bbs[l]) > 0:
            v.get_label_by_id("100").set_text(round(bbs_deut_counts[l] / 1000, 1))
        if len(bbs[l2]) > 0:
            v.get_label_by_id("010").set_text(round(bbs_deut_counts[l2] / 1000, 1))
        if samezies > 0:
            if round(samezies/1000,1) < 1:
                v.get_label_by_id("110").set_text(math.ceil(round(samezies/1000,1)))
            else:
                v.get_label_by_id("110").set_text(round(samezies/1000,1))

        print(x,y,l, l2)
        print(v)
        if x == 1 and y == 2:
            left = v.get_label_by_id("100")
            left.set_x(left._x - .20)
        if x == 1 and y == 3:
            left = v.get_label_by_id("100")
            left.set_x(left._x - .04)
        if x == 4 and y == 5:
            right = v.get_label_by_id("010")
            right.set_x(right._x + .20)
        if x == 0 and y == 1:
            right = v.get_label_by_id("010")
            right.set_x(right._x + .20)
        if x == 2 and y == 4:
            right = v.get_label_by_id("010")
            right.set_x(right._x + .08)
        v.patches[0].set_color("#0b1b82")
        v.patches[1].set_color('#f2f758')
        i = i + 1
        y = y + 1
    seen.append(l)
    ystart = ystart + 1
    x = x + 1

labels = ["Acid", "Aldehyde", "Alcohol", "Amine", "Halide", "Boronate"]
i = 0
for x in range(0,n):
    for y in range(0,n):
        if x == y:
            print(x,y)
            ax[x][y].set_xticks([])
            ax[x][y].set_yticks([])
            ax[x][y].axis('off')
            i = i + 1

ax[5][5].set_xlim([0,1])
ax[5][5].set_ylim([0,1])
seen = []

x = 0
xbreak = 1
ystart = 1
for l in labels_hal:
    if "sp2" not in l:
        continue
    i = 0
    y = ystart
    for l2 in labels_hal:
        if l2 in seen:
            continue
        if y >= n:
            break
        if "sp2" not in l2:
            continue
        if l.split("_")[0] == l2.split("_")[0]:
            continue
        if l.split("_")[1] != l2.split("_")[1]:
            continue
        if "bromide" in l or "iodide" in l or "chloride" in l or "bromide" in l2 or "iodide" in l2 or "chloride" in l2:
            continue
        # print(x,y)
        samezies = 0
        for a1 in bbs_deut[l]:
            for a2 in bbs_deut[l2]:
                if a1 == a2:
                    samezies = samezies + 1
        # print(l + " same " + l2)
        # print(len(aa[l]), samezies, len(aa[l2]))
        v = venn2(subsets = (bbs_deut_counts[l], bbs_deut_counts[l2], samezies), set_labels = ("", "", ""), ax=ax[y][x])
        if len(bbs[l]) > 0:
            v.get_label_by_id("100").set_text(round(bbs_deut_counts[l] / 1000, 1))
        if len(bbs[l2]) > 0:
            v.get_label_by_id("010").set_text(round(bbs_deut_counts[l2] / 1000, 1))
        if samezies > 0:
            if round(samezies/1000,1) < 1:
                v.get_label_by_id("110").set_text(math.ceil(round(samezies/1000,1)))
            else:
                v.get_label_by_id("110").set_text(round(samezies/1000,1))
        if y == 4 and x == 1:
            left = v.get_label_by_id("100")
            left.set_x(left._x - .05)
        if y == 5 and x == 4:
            right = v.get_label_by_id("010")
            right.set_x(right._x + .05)
        v.patches[0].set_color("#0b1b82")
        v.patches[1].set_color('#f2f758')
        i = i + 1
        y = y + 1
    seen.append(l)
    ystart = ystart + 1
    x = x + 1
plt.subplots_adjust(hspace = 0, wspace = 0)
plt.margins(0,0)
fig.savefig("venn_counts.png", dpi=900, bbox_inches="tight", pad_inches = 0)
plt.close()