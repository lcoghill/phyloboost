### Read in a parse the tree file output from 'phylota-graph'

treeFile = open('tree_file.txt', 'r')

ci = [];

for line in treeFile:
    label = line.split('	', 1 ); #split the label id from the rest of the tree
    ci_ti = label[0].split('_', 1 ); #split the label into ti and ci
    ci.append (ci_ti[1]) #create a list of ci values
    
print ci

