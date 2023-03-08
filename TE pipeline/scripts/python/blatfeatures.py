from datetime import datetime
starttime = datetime.now()

print('Running...')
print('Blat features...')

from difflib import SequenceMatcher

def fasta_reader(filename):
  from Bio.SeqIO.FastaIO import FastaIterator
  with open(filename) as handle:
    for record in FastaIterator(handle):
      yield record

rowseq = []
rowid = []
for entry in fasta_reader('blatmatch.fa.out'):
  rowseq.append(str(entry.seq))
  rowid.append(str(entry.id))

#view file of id and name of blat match results
with open('blatmatchf.txt') as i:
    linesb = i.readlines()

locfasta = []
queryb = []
matchb = []
startq = []
endq = []
sizeq = []
chrb = []
strand = []
startb = []
endb = []
locfastas = []
for i in range(len(linesb)):
    bseq = linesb[i]
    bt = [p for p, marker in enumerate(bseq) if marker == "\t"]
    locfasta.append(bseq[:bt[0]])
    queryb.append(bseq[bt[0]+1:bt[1]])
    matchb.append(bseq[bt[1]+1:bt[2]])
    startq.append(bseq[bt[2]+1:bt[3]])
    endq.append(bseq[bt[3]+1:bt[4]])
    sizeq.append(bseq[bt[4]+1:bt[5]])
    chrb.append(bseq[bt[5]+1:bt[6]])
    strand.append(bseq[bt[6]+1:bt[7]])
    startb.append(bseq[bt[7]+1:bt[8]])
    endb.append(bseq[bt[8]+1:-1])
    locfastas.append(locfasta[i].strip())

#view file of structural match results
with open('matchedta.txt') as i:
    linesm = i.readlines()

codem = []
chrm = []
startm = []
endm = []
strandm = []
classm = []
subclassm = []
idm = []
namem = []
methodm = []
motifm = []
tsd1m = []
tsd2m = []
matchtsdm = []
tir1m = []
tir2m = []
matchtirm = []
lenseqm = []
for i in range(len(linesm)):
    mseq = linesm[i]
    mt = [p for p, marker in enumerate(mseq) if marker == "\t"]
    codem.append(mseq[:mt[0]])
    chrm.append(mseq[mt[0]+1:mt[1]])
    startm.append(mseq[mt[1]+1:mt[2]])
    endm.append(mseq[mt[2]+1:mt[3]])
    strandm.append(mseq[mt[3]+1:mt[4]])
    classm.append(mseq[mt[4]+1:mt[5]])
    subclassm.append(mseq[mt[5]+1:mt[6]])
    idm.append(mseq[mt[6]+1:mt[7]])
    namem.append(mseq[mt[7]+1:mt[8]])
    methodm.append(mseq[mt[8]+1:mt[9]])
    motifm.append(mseq[mt[9]+1:mt[10]])
    mseq2 = mseq[mt[10]+1:-1]
    mt2 = [p for p, marker in enumerate(mseq2) if marker == "\t"]
    mu = [p for p, marker in enumerate(mseq2) if marker == "_"]
    tsd1m.append(mseq2[:mu[0]])
    tsd2m.append(mseq2[mu[0]+1:mu[1]])
    matchtsdm.append(mseq2[mu[1]+1:mt2[0]])
    tir1m.append(mseq2[mt2[0]+1:mu[2]])
    tir2m.append(mseq2[mu[2]+1:mu[3]])
    matchtirm.append(mseq2[mu[3]+1:mt2[1]])
    lenseqm.append(mseq2[mt2[1]+1:-1])


fblat = open('blatfeatures.txt', 'w')
fblat.write('QueryID\tChromosome\tStart\tEnd\tStrand\tClassification\tSubclassification\tQueryName\tMethod\tTSD\tTIR')

#run for all sequences
num = 0
for x in rowseq:
    seq = rowseq[num]
    fastaname = rowid[num]
    #clean sequence
    seq = ' '.join([i for i in seq if not i.isdigit()])
    seqnobreaks =  seq.replace("\n", " ")
    seqclean = seqnobreaks.replace(' ', '')
    seqcleanu = seqclean.lower()

    #get extra bp
    first = seqcleanu[:20]
    last = seqcleanu[-20:]
    middle = seqcleanu[20:-20]

    indexblat = locfastas.index(fastaname)
    query = queryb[indexblat]
    id = idm.index(query)
    tsd = tsd1m[id]
    lentsd = len(tsd)
    tir1 = tir1m[id]
    lentir1 = len(tir1)
    tir2 = tir2m[id]
    lentir2 = len(tir2)
    #blat tsd and tir
    tsd1b = first[-lentsd:]
    tsd2b = last[:lentsd]
    tir1b = middle[:lentir1]
    tir2f = middle[-lentir2:]
    tir2r = tir2f[::-1]
    tir2g = tir2r.replace('g', 'C')
    tir2c = tir2g.replace('c', 'G')
    tir2a = tir2c.replace('a', 'T')
    tir2t = tir2a.replace('t', 'A')
    tir2b = tir2t.lower()
    matchtsd = SequenceMatcher(None, tsd1b, tsd2b).ratio()
    matchtir = SequenceMatcher(None, tir1b, tir2b).ratio()
    matchtsd1 = matchtsd * 100
    matchtir1 = matchtir * 100
    matchtsdp = '{:.2f}'.format(matchtsd1)
    matchtirp = '{:.2f}'.format(matchtir1)
    if matchtsd >= 1:
        if matchtir >= 0.85:
            fblat.write('\n' + query + '\t' + chrb[num] + '\t' + startb[num] + '\t' + endb[num] + '\t' + strand[num] + '\t' + classm[id] + '\t' + subclassm[id] + '\t' + namem[id] + '\t' + 'Homology' + '\t' + tsd1b.upper() + '_' + tsd2b.upper() + '_' + str(matchtsdp) + '\t' + tir1b.upper() + '_' + tir2f.upper() + '_' + str(matchtirp))

    #for loop counter
    num = num + 1

fblat.close()

endtime = datetime.now()
print('Run time: {}'.format(endtime - starttime))
