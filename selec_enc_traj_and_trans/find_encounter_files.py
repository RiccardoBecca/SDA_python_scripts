import sys

total_solutes = int(sys.argv[1]) # p1 + p2 + crowders
subname = sys.argv[2]
files = sys.argv[3:]
to_transfer = []

for file in files:
    with open(file,"r") as fin:
        lines = fin.readlines()
    if len(lines) > total_solutes + 2: # +2 because of the header in encounter complex file
        to_transfer.append(file.replace(subname,''))
#print(to_transfer)
bash_string=''
for num in to_transfer:
    bash_string=bash_string+" "+num
print(bash_string)
