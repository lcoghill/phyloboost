from Bio import SeqIO


def parse_config():
    config_params = []
    config = open('private/config', 'r')
    for line in (line for line in config if not line.startswith('###')):
        line = line.rstrip('\n')
        line = line.split("=")
        config_params.append(line[1])

    return config_params


## Gather import configuration information from the config file
params = parse_config()  # retrieve the params from the config file



## Gather import configuration information from the config file
print "Please enter the name of your nucleotide file from Genbank."
print "The script will look in %s (specified input location) for the file." % params[7]
file_name = raw_input("Filename: ")
full_file_name = "".join([params[7], file_name])
print "Opening %s" % file_name

## Iterate through each record, parse out the important information and write that information to a csv file
f = open(full_file_name, 'r+')  # nt fasta file
csv_file = "".join([params[6], file_name, ".csv"])
csv = open(csv_file, "w+")  # output csv file

count = 0
for record in SeqIO.parse(full_file_name, "fasta"):
    seq_list = record.description.split("|")
    ## convert record to string
    csv_string = '%s, %s, %s, %s, %s' % (seq_list[1], seq_list[3], len(record.seq), seq_list[4][1:], record.seq)
    csv.write("%s\n" % csv_string)
    print "Record %s converted to CSV format" % count
    count = count + 1
f.close()
csv.close()