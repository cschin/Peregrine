import random

random.seed(42)

rcmap = dict(zip("ACGT","TGCA"))

def rc_seq(seq):
    return "".join([rcmap[c] for c in seq[::-1]])

def sim_error(seq):
    out_seq = []
    for c in seq:
        if random.uniform(0, 1) < 0.01:
            c = random.choice( ('A','C','G','T', '', c+'A', c+'C', c+'G', c+'T') )
        out_seq.append(c)
    return "".join(out_seq)

seq = []
with open("./K12MG1655.fa") as f:
    for row in f:
        row = row.strip()
        if len(row) < 1:
            continue
        if ">" == row[0]:
            continue
        seq.append(row)

seq = "".join(seq)
seq = seq + seq[:40000]

rl = 15000
for j in range(8):
    read_count = 2 * len(seq) // rl
    sim_record = open(f"reads/reads_{j}.bed","w")
    read_file = open(f"reads/reads_{j}.fa","w")
    import random
    for i in range(read_count):
        rl2 = int(rl + random.gauss(0, 1500))
        s  = random.randint(0, len(seq)-40000)
        print(">{:02d}/{:06d}/{}_{}".format(j,i,0,rl2), file=read_file)
        seq_tmp = sim_error(seq[s:s+rl2])
        if random.randint(0,1) == 1:
            seq_tmp = rc_seq(seq_tmp)
        print(seq_tmp, file=read_file)
        print("{:02d}/{:06d}/{}_{}".format(j,i,0,rl2), s, s+rl2, sep="\t", file=sim_record)
    sim_record.close()
    read_file.close()

