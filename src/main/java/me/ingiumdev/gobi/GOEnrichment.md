```cat ~/IdeaProjects/gobi/data/GOEnrich/goa_human.gaf | grep -v "\!" | cut -f 4 | sort | uniq -c
391175
1399 NOT
14 NOT|colocalizes_with
9 NOT|contributes_to
1135 colocalizes_with
1131 contributes_to
```

--> We only want the files don't have anything in the third column
for ensembl, exclude the rows with nothing in hgnc

Don't only propagate the counts, propagate the gene ID too

breadth first search
continue searching up until you can't better your foudn LCA. If you found one at level 10, once you go above that, you
don't need to search anymore

