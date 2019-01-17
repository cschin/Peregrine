

fn = "preads4falcon_mer"
mer_count = {}
with open(fn) as f:
    for row in f:
        row = row.strip()
        if row[0] == ">":
            continue
        row = row.split()
        mer_count.setdefault(row[2], 0)
        mer_count[row[2]] += 1

with open(fn) as f:
    for row in f:
        row = row.strip()
        if row[0] == ">":
            print(row)
            continue
        else:
            row = row.split()
            count = mer_count[row[2]]
            print(" ".join(row), count)

