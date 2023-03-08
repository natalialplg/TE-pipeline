#mainproject extraction of tsd/tir characteristics
from datetime import datetime
starttime = datetime.now()

print('Running...')
print('TSD and TIR features...')

import re
from difflib import SequenceMatcher

def fasta_reader(filename):
  from Bio.SeqIO.FastaIO import FastaIterator
  with open(filename) as handle:
    for record in FastaIterator(handle):
      yield record

rowseq = []
rowid = []
for entry in fasta_reader('outputnohseq.fa.out'):
  rowseq.append(str(entry.seq))
  rowid.append(str(entry.id))

#view file of EDTA results
with open('edta_seqe_filter_noh.txt') as i:
    linesc = i.readlines()

code = []
title = []
all = []
for i in range(len(linesc)):
    cseq = linesc[i]
    if i == 0:
        title.append(cseq)
    else:
        spc = [p for p, marker in enumerate(cseq) if marker == "\t"]
        code.append(cseq[spc[0]+1:spc[1]])
        all.append(cseq[spc[0]+1:-1])

with open('edta_seqe_filter_noh20.bed.txt') as i:
    lines = i.readlines()

scaffold = []
start = []
end = []
for i in range(len(lines)):
    tseq = lines[i]
    spt = [p for p, marker in enumerate(tseq) if marker == "\t"]
    scaffold.append(tseq[0:spt[0]])
    start.append(tseq[spt[0]+1:spt[1]])
    end.append(tseq[spt[1]+1:-1])

fgene = open('tsdtir.txt', 'w')
fgene.write('Number\tCode\tIdentifier\tScaffold\tStart\tEnd\tClassification\tSubclassification\tTSD\tlengthTSD\tTIR\tlengthTIR\tlengthSequence\tAutonomy')

fseq = open('tsdtirseq.bed.txt', 'w')

fteres = open('validtsdtir.txt', 'w')
fteres.write('Number\tCode\tIdentifier\tScaffold\tStart\tEnd\tClassification\tSubclassification\tTSD\tlengthTSD\tTIR\tlengthTIR\tlengthSequence\tAutonomy')

fteresedta = open('validtsdtiredta.txt', 'w')
fteresedta.write('Number\tCode\tScaffold\tProgram\tMainClassification\tStart\tEnd\tScore\tStrand\tPhase\tID\tParent\tName\tClassification\tSubclassification\tSequenceOntology\tIdentity\tMethod\tMotif\tTSD1\tTSD2\tmatchTSD\tTIR1\tTIR2')

valid = 0
dna = 0
ltr = 0
other = 0

#run for all sequences
num = 0
genenumber = []
allcounter = []
seqid = []
for x in rowseq:
    genecount = num + 1
    seq = rowseq[num]

    #clean sequence
    seq = ' '.join([i for i in seq if not i.isdigit()])
    seqnobreaks =  seq.replace("\n", " ")
    seqclean = seqnobreaks.replace(' ', '')
    seqcleanu = seqclean.lower()

    #get first # and last # of nucleotides
    first = seqcleanu[:39]
    last = seqcleanu[-40:]
    both = first + last

    #find unknown patterns (potential tsd)
    d = {}
    mincount = 2
    minlength = 2
    for sublength in range(minlength,int(len(seqcleanu)/mincount)):
        for i in range(0,len(both)-sublength):
            sub = both[i:i+sublength]
            cnt = both.count(sub)
            if cnt >= mincount and sub not in d:
                d[sub] = cnt

    #order patterns by length, longest to shortest
    patterns = []
    for i in set(d):
        patterns.append(i)
    patterns.sort(key = len, reverse = True)
    nonpatterns = [elem.replace('n', '') for elem in patterns]
    newpatterns = [x for x in nonpatterns if len(x) < 15]
    lenpatterns = len(newpatterns)
    cnt = 0
    tescaffold =[]
    testart = []
    teend = []
    classification = []
    subclassification = []
    autonomy = []
    tsd = []
    tir = []
    fullsequences = []
    sequencestsd = []
    sequences = []
    classificationyes = []
    subclassificationyes = []
    autonomyyes = []
    tirsecond = []
    tsdyes = []
    tiryes = []
    tir2yes =[]
    sequencestsdyes = []
    fullsequencesyes = []
    sequencesyes = []
    tescaffoldyes =[]
    testartyes = []
    teendyes = []
    lent = 0
    lentseq = 0
    lentsd = 0

    #run all patterns for tsd and tir
    equalcounter = 0
    while cnt < lenpatterns:
        #find position
        pattpos = []
        pattposleft = []
        pattposright = []
        for a in re.finditer(newpatterns[cnt], both):
            pattpos.append(a.start())
            pattpos.append(a.end())
        pattpos.remove(pattpos[len(pattpos) - 1])
        pattpos.remove(pattpos[0])

        for i in range(len(pattpos)):
            if 0 < pattpos[i] < 39:
                pattposleft.append(pattpos[i])
            if pattpos[i] > 39:
                pattposright.append(-pattpos[i] + 40)
        listpattpos = [(x,y) for x in pattposleft for y in pattposright]

        #get central sequence
        for i in range(len(listpattpos)):
            if len(newpatterns[cnt]) >= 2:
                eachpattpos = listpattpos[i]
                loc = eachpattpos[0]
                loc2 = eachpattpos[1]-1
                potteseq = seqcleanu[loc:loc2]
                startnum = int(start[genecount-1])+loc+1
                endnum = int(end[genecount-1])+loc2
                reverse = potteseq[::-1]
                reverseg = reverse.replace('g', 'C')
                reversec = reverseg.replace('c', 'G')
                reversea = reversec.replace('a', 'T')
                reverset = reversea.replace('t', 'A')
                reverseclean = reverset.lower()

                #find TIR and its length
                s = 1
                match1 = 0
                tir1 = potteseq[:s]
                tir2 = reverseclean[:s]
                for i in range(1,150):
                    s = s + 1
                    tir1 = potteseq[:s]
                    tir2 = reverseclean[:s]
                    match = SequenceMatcher(None, tir1, tir2).ratio()
                    if match < 0.85:
                        break
                    elif match >= 0.85:
                        #tir with tg-ca format
                        if tir1 == 'tg':
                            q = 3
                            p = 1
                            tir1 = potteseq[p:q]
                            for i in range(1,30):
                                tir1 = potteseq[p:q]
                                if tir1 == 'ca' or tir1 == 'ga' :
                                    tir2tg = potteseq[-q:]
                                    if tir2tg[:2] == 'tg' :
                                            if len(potteseq[:q]) < 15:
                                                lentwo = len(newpatterns[cnt]) + len(potteseq[:q])
                                                if lentwo > lent:
                                                    lenwofirst = len(potteseq)
                                                    lent = lentwo
                                                    if lenwofirst > lentseq:
                                                        equalcounter = equalcounter + 1
                                                        correcttsd = newpatterns[cnt]
                                                        correcttir = potteseq[:q]
                                                        lentsd = len(newpatterns[cnt])
                                                        if (seqcleanu.find(correcttsd+correcttir) != -1):
                                                            if (seqcleanu.find((potteseq[-len(correcttir):])+correcttsd) != -1):
                                                                if match >= match1:
                                                                    match1 = match
                                                                    tsd.append(correcttsd)
                                                                    tir.append(correcttir)
                                                                    tirsecond.append(potteseq[-len(correcttir):])
                                                                    tescaffold.append(scaffold[genecount-1])
                                                                    testart.append(startnum)
                                                                    teend.append(endnum)
                                                                    sequences.append(potteseq)
                                                                    if 3 <= len(correcttsd) <= 6:
                                                                        classification.append('LTR')
                                                                        subclassification.append('COPIA/GYPSY')
                                                                    else:
                                                                        classification.append('OTH1')
                                                                        subclassification.append('unk')
                                q = q + 1
                                p = p + 1
                        elif tir1 == 'cact':
                            m = 5
                            for i in range(1,30):
                                tir1 = potteseq[:m]
                                tir2 = reverseclean[:m]
                                if tir1 != tir2:
                                    break
                                else:
                                    lentwo = len(newpatterns[cnt]) + len(potteseq[:m])
                                    if lentwo > lent:
                                        lenwofirst = len(potteseq)
                                        lent = lentwo
                                        if lenwofirst > lentseq:
                                            equalcounter = equalcounter + 1
                                            correcttsd = newpatterns[cnt]
                                            correcttir = potteseq[:m]
                                            lentsd = len(newpatterns[cnt])
                                            if (seqcleanu.find(correcttsd+correcttir) != -1):
                                                if (seqcleanu.find((potteseq[-len(correcttir):])+correcttsd) != -1):
                                                    if match >= match1:
                                                        match1 = match
                                                        tsd.append(correcttsd)
                                                        tir.append(correcttir)
                                                        tirsecond.append(potteseq[-len(correcttir):])
                                                        tescaffold.append(scaffold[genecount-1])
                                                        testart.append(startnum)
                                                        teend.append(endnum)
                                                        sequences.append(potteseq)
                                                        classification.append('DNA')
                                                        if 2 <= len(correcttsd) <= 3:
                                                            subclassification.append('DTC/CACTA')
                                                        else:
                                                            subclassification.append('unk')
                                m = m + 1
                        # match perfect tsd and tir
                        else:
                            lentwo = len(newpatterns[cnt]) + len(tir1)
                            if lentwo > lent:
                                lenwofirst = len(potteseq)
                                lent = lentwo
                                if lenwofirst > lentseq:
                                    equalcounter = equalcounter + 1
                                    lentsd = len(newpatterns[cnt])
                                    correcttsd = newpatterns[cnt]
                                    correcttir = tir1
                                    if (seqcleanu.find(correcttsd+correcttir) != -1):
                                        if (seqcleanu.find((potteseq[-len(correcttir):])+correcttsd) != -1):
                                            if match >= match1:
                                                match1 = match
                                                tsd.append(correcttsd)
                                                tir.append(correcttir)
                                                tirsecond.append(potteseq[-len(correcttir):])
                                                tescaffold.append(scaffold[genecount-1])
                                                testart.append(startnum)
                                                teend.append(endnum)
                                                sequences.append(potteseq)
                                                if 2 <= len(correcttsd) <= 11:
                                                    #DTT Tc1 Mariner
                                                    if correcttsd == 'ta':
                                                        classification.append('DNA')
                                                        subclassification.append('DTT/Tc1-Mariner')
                                                    #DTH PIF-Harbringer
                                                    elif correcttsd == 'taa' or correcttsd == 'tta':
                                                        classification.append('DNA')
                                                        subclassification.append('DTH/PIF-Harbringer')
                                                    #DTA hAT
                                                    elif len(correcttsd) == 8:
                                                        classification.append('DNA')
                                                        subclassification.append('DTA/hAT')
                                                    #DTM Mutator
                                                    elif len(correcttsd) == 7:
                                                        classification.append('DNA')
                                                        subclassification.append('DTM/Mutator')
                                                    elif 9 <= len(correcttsd) <= 11:
                                                        classification.append('DNA')
                                                        subclassification.append('DTM/Mutator')
                                                    elif 3 <= len(correcttsd) <= 6:
                                                        if correcttir[:2] == 'tg':
                                                            secondtir = potteseq[-len(correcttir):]
                                                            if secondtir[-2:] == 'ca':
                                                                classification.append('LTR')
                                                                subclassification.append('COPIA/GYPSY')
                                                            elif secondtir[-2:] == 'ga':
                                                                classification.append('LTR')
                                                                subclassification.append('COPIA/GYPSY')
                                                            elif secondtir[-2:] == 'ta':
                                                                classification.append('LTR')
                                                                subclassification.append('COPIA/GYPSY')
                                                            elif len(potteseq) > 6000:
                                                                classification.append('LTR')
                                                                subclassification.append('COPIA/GYPSY')
                                                            else:
                                                                classification.append('LTR')
                                                                subclassification.append('unk')
                                                        else:
                                                            classification.append('DNA')
                                                            subclassification.append('unk')
                                                    else:
                                                        classification.append('DNA')
                                                        subclassification.append('unk')
                                                else:
                                                    classification.append('OTH2')
                                                    subclassification.append('unk')

        #while counter
        cnt = cnt + 1

    # for sequences without a defined tsd and tir structure
    if tsd == []:
        tsd.append('-')
        tir.append('-')
        tirsecond.append('-')
        tescaffold.append('-')
        testart.append('-')
        teend.append('-')
        sequences.append('-')
        classification.append('-')
        subclassification.append('-')

    tsdyes.append(tsd[-1])
    sequencesyes.append(sequences[-1])
    classificationyes.append(classification[-1])
    subclassificationyes.append(subclassification[-1])
    tiryes.append(tir[-1])
    tir2yes.append(tirsecond[-1])
    allcounter.append(genecount)
    seqid.append(rowid[num])
    sequencestsdyes.append(seqcleanu)
    tescaffoldyes.append(tescaffold[-1])
    testartyes.append(testart[-1])
    teendyes.append(teend[-1])

    for i in range(len(sequencesyes)):
        if subclassification[i] == '-':
            autonomyyes.append('-')
        elif len(potteseq) <= 600:
            autonomyyes.append('non-aut')
        else:
            autonomyyes.append('aut')

    tirrev = tirsecond[-1]
    reverse = tirrev[::-1]
    reverseg = reverse.replace('g', 'C')
    reversec = reverseg.replace('c', 'G')
    reversea = reversec.replace('a', 'T')
    reverset = reversea.replace('t', 'A')
    reverseclean = reverset.lower()
    matchtir = SequenceMatcher(None, tir[-1], reverseclean).ratio()
    matchtir1 = matchtir * 100
    matchtsd = SequenceMatcher(None, tsd[-1], tsd[-1]).ratio() * 100
    matchtsdp = '{:.2f}'.format(matchtsd)
    matchtirp = '{:.2f}'.format(matchtir1)
    #append file
    for i in range(len(tescaffoldyes)):
        fgene.write('\n')
        if tescaffoldyes[i] == '-':
            fgene.write(str(allcounter[genecount-1]) + '\t' + code[genecount-1] + '\t' + str(seqid[-1]) + '\t' + str(tescaffoldyes[-1]) + '\t' + str(testartyes[-1]) + '\t' + str(teendyes[-1]) + '\t' + classificationyes[i] + '\t' + subclassificationyes[i] + '\t' + tsdyes[i] + '_' + tsdyes[i] + '_' + matchtsdp + '\t' + '-' + '\t' + tiryes[i] + '_' + tir2yes[i] + '_' + matchtirp + '\t' + '-' + '\t' + '-' + '\t' + autonomyyes[i])
        else:
            fgene.write(str(allcounter[genecount-1]) + '\t' + code[genecount-1] + '\t' + str(seqid[-1]) + '\t' + str(tescaffoldyes[-1]) + '\t' + str(testartyes[-1]) + '\t' + str(teendyes[-1]) + '\t' + classificationyes[i] + '\t' + subclassificationyes[i] + '\t' + tsdyes[i] + '_' + tsdyes[i] + '_' + matchtsdp + '\t' + str(len(tsdyes[i])) + '\t' + tiryes[i] + '_' + tir2yes[i] + '_' + matchtirp + '\t' + str(len(tiryes[i])) + '\t' + str(len(sequencesyes[i])) + '\t' + autonomyyes[i])
            fteres.write('\n' + str(valid+1) + '\t' + code[genecount-1] + '\t' + str(seqid[-1]) + '\t' + str(tescaffoldyes[-1]) + '\t' + str(testartyes[-1]) + '\t' + str(teendyes[-1]) + '\t' + classificationyes[i] + '\t' + subclassificationyes[i] + '\t' + tsdyes[i] + '_' + tsdyes[i] + '_' + matchtsdp + '\t' + str(len(tsdyes[i])) + '\t' + tiryes[i] + '_' + tir2yes[i] + '_' + matchtirp + '\t' + str(len(tiryes[i])) + '\t' + str(len(sequencesyes[i])) + '\t' + autonomyyes[i])
            fteresedta.write('\n' + str(valid+1) + '\t' + all[genecount-1])
            valid = valid + 1
    for i in range(len(tescaffoldyes)):
        fseq.write(str(tescaffoldyes[-1]) + '\t' + str(testartyes[-1]) + '\t' + str(teendyes[-1]))
        fseq.write('\n')

    #for loop counter
    num = num + 1

fgene.close()
fseq.close()
fteres.close()
fteresedta.close()

endtime = datetime.now()
print('Run time: {}'.format(endtime - starttime))
