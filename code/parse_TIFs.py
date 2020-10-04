import csv

with open('/tmp/output.tsv', 'wt') as out_file:
    tsv_writer = csv.writer(out_file, delimiter='\t')
    tsv_writer.writerow(['name', 'field'])
    tsv_writer.writerow(['Dijkstra', 'Computer Science'])
    tsv_writer.writerow(['Shelah', 'Math'])
    tsv_writer.writerow(['Aumann', 'Economic Sciences'])
f = open("./SteinmetzGilbert/nature12121-s2/S1_TIFs.txt", "r")

out_file = open("parsed_steinmetz_s1_tifs.txt", "w")
tsv_writer = csv.writer(out_file, delimiter='\t')

lines = f.readlines()

# header
tsv_writer.writerow(lines[0].split())

x = dict()

not_gene_names = ["ORF", "transcripts", "ORFs"]

# now the fun part
for line in lines[1:]:
    data = line.split()[:6]
    name = line.split()[-1]

    if name in not_gene_names:
        type_ = " ".join(line.split()[6:])
        stuff = data + [type_, "N/A"]
    else:
        type_ = " ".join(line.split()[6:-1])
        stuff = data + [type_, name]

    tsv_writer.writerow(stuff)
