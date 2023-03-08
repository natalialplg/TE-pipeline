from datetime import datetime
start = datetime.now()

print('Running...')
print('Filter EDTA...')

#split data from EDTA into its corresponding sections
with open('seqextraction.txt') as i:
    lines = i.readlines()

gff = []
number = []
code = []
scaffold = []
program = []
mainclass = []
startedta = []
endedta = []
score = []
strand = []
phase = []
id = []
parent  = []
name = []
classedta = []
subclassedta = []
seqontology = []
identity = []
method = []
motif = []
tsd1edta = []
tsd2edta = []
matchtsd = []
tir1edta = []
tir2edta = []
codenumber = 10000
scaffoldnumber = 1
for i in range(len(lines)):
    edtaseq = lines[i]
    gff.append(edtaseq)
    sp = [p for p, marker in enumerate(edtaseq) if marker == "\t"]
    number.append(i+1)
    scaffold.append(edtaseq[0:sp[0]])
    codenumber = codenumber + 1
    if scaffoldnumber != scaffold[i]:
        scaffoldnumber = scaffold[i]
        codenumber = 10000
    code.append('GS' + scaffold[i] + 'TE' + str(codenumber))
    program.append(edtaseq[sp[0]+1:sp[1]])
    mainclass.append(edtaseq[sp[1]+1:sp[2]])
    startedta.append(edtaseq[sp[2]+1:sp[3]])
    endedta.append(edtaseq[sp[3]+1:sp[4]])
    score.append(edtaseq[sp[4]+1:sp[5]])
    strand.append(edtaseq[sp[5]+1:sp[6]])
    phase.append(edtaseq[sp[6]+1:sp[7]])
    more = edtaseq[sp[7]+1:]
    eq = [p for p, marker in enumerate(more) if marker == "="]
    pc = [p for p, marker in enumerate(more) if marker == ";"]
    sl = [p for p, marker in enumerate(more) if marker == "/"]
    id.append(more[eq[0]+1:pc[0]])
    if 'Parent' in more:
        parent.append(more[eq[1]+1:pc[1]])
        name.append(more[eq[2]+1:pc[2]])
        classedta.append(more[eq[3]+1:sl[0]])
        subclassedta.append(more[sl[0]+1:pc[3]])
        seqontology.append(more[eq[4]+1:pc[4]])
        identity.append(more[eq[5]+1:pc[5]])
        method.append(more[eq[6]+1:pc[6]])
        motif.append(more[eq[7]+1:pc[7]])
        tsd1edta.append(more[eq[8]+1:-1])
        tsd2edta.append('-')
        matchtsd.append('-')
        tir1edta.append('-')
        tir2edta.append('-')
    else:
        parent.append('-')
        name.append(more[eq[1]+1:pc[1]])
        classedta.append(more[eq[2]+1:sl[0]])
        subclassedta.append(more[sl[0]+1:pc[2]])
        seqontology.append(more[eq[3]+1:pc[3]])
        identity.append(more[eq[4]+1:pc[4]])
        firstmethod = more[eq[5]+1:]
        if identity[i] == 'NA':
            method.append(more[eq[5]+1:-1])
            motif.append('-')
            tsd1edta.append('-')
            tsd2edta.append('-')
            matchtsd.append('-')
            tir1edta.append('-')
            tir2edta.append('-')
        elif firstmethod[:1] == 'h':
            method.append(more[eq[5]+1:-1])
            motif.append('-')
            tsd1edta.append('-')
            tsd2edta.append('-')
            matchtsd.append('-')
            tir1edta.append('-')
            tir2edta.append('-')
        else:
            strucmore = more[eq[5]+1:]
            spc = [p for p, marker in enumerate(strucmore) if marker == ";"]
            seq = [p for p, marker in enumerate(strucmore) if marker == "="]
            method.append(strucmore[:spc[0]])
            if 'motif' in strucmore:
                motif.append(strucmore[seq[0]+1:spc[1]])
                tsd1edta.append(strucmore[seq[1]+1:-1])
                tsd2edta.append('-')
                matchtsd.append('-')
                tir1edta.append('-')
                tir2edta.append('-')
            else:
                struc2more = strucmore[seq[0]+1:]
                s2pc = [p for p, marker in enumerate(struc2more) if marker == ";"]
                s2eq = [p for p, marker in enumerate(struc2more) if marker == "="]
                us = [p for p, marker in enumerate(struc2more) if marker == "_"]
                motif.append('-')
                tsd1edta.append(struc2more[:us[0]])
                tsd2edta.append(struc2more[us[0]+1:us[1]])
                matchtsd.append(struc2more[us[1]+1:s2pc[0]])
                tir1edta.append(struc2more[s2eq[0]+1:us[2]])
                tir2edta.append(struc2more[us[2]+1:-1])

fedta = open('edta_seqe.txt', 'w')
fedta.write('Number\tCode\tScaffold\tProgram\tMainClassification\tStart\tEnd\tScore\tStrand\tPhase\tID\tParent\tName\tClassification\tSubclassification\tSequenceOntology\tIdentity\tMethod\tMotif\tTSD1\tTSD2\tmatchTSD\tTIR1\tTIR2')
for i in range(len(number)):
   fedta.write('\n' + str(i+1) + '\t' + code[i] + '\t' + scaffold[i] + '\t' + program[i] + '\t' + mainclass[i] + '\t' + startedta[i] + '\t' + endedta[i] + '\t' + score[i] + '\t' + strand[i] + '\t' + phase[i] + '\t' + id[i] + '\t' + parent[i] + '\t' + name[i] + '\t' + classedta[i] + '\t' + subclassedta[i] + '\t' + seqontology[i] + '\t' + identity[i] + '\t' + method[i] + '\t' + motif[i] + '\t' + tsd1edta[i] + '\t' + tsd2edta[i] + '\t' + matchtsd[i] + '\t' + tir1edta[i] + '\t' + tir2edta[i] +  '\t')
fedta.close()

# filter out helitrons
alledta = list(zip(number, code, scaffold, program, mainclass, startedta, endedta, score, strand, phase, id, parent, name, classedta, subclassedta, seqontology, identity, method, motif, tsd1edta, tsd2edta, matchtsd, tir1edta, tir2edta))
nohalledta = []
for i in range(len(number)):
    mainclassvalue = mainclass[i]
    if mainclassvalue != 'helitron':
        nohalledta.append(alledta[i])

nohnumber, nohcode, nohscaffold, nohprogram, nohmainclass, nohstartedta, nohendedta, nohscore, nohstrand, nohphase, nohid, nohparent, nohname, nohclassedta, nohsubclassedta, nohseqontology, nohidentity, nohmethod, nohmotif, nohtsd1edta, nohtsd2edta, nohmatchtsd, nohtir1edta, nohtir2edta = zip(*nohalledta)

fnohedta = open('edta_seqe_noh.txt', 'w')
fnohedta.write('Number\tCode\tScaffold\tProgram\tMainClassification\tStart\tEnd\tScore\tStrand\tPhase\tID\tParent\tName\tClassification\tSubclassification\tSequenceOntology\tIdentity\tMethod\tMotif\tTSD1\tTSD2\tmatchTSD\tTIR1\tTIR2')
for i in range(len(nohnumber)):
   fnohedta.write('\n' + str(i+1) + '\t' + nohcode[i]  + '\t' + nohscaffold[i] + '\t' + nohprogram[i] + '\t' + nohmainclass[i] + '\t' + nohstartedta[i] + '\t' + nohendedta[i] + '\t' + nohscore[i] + '\t' + nohstrand[i] + '\t' + nohphase[i] + '\t' + nohid[i] + '\t' + nohparent[i] + '\t' + nohname[i] + '\t' + nohclassedta[i] + '\t' + nohsubclassedta[i] + '\t' + nohseqontology[i] + '\t' + nohidentity[i] + '\t' + nohmethod[i] + '\t' + nohmotif[i] + '\t' + nohtsd1edta[i] + '\t' + nohtsd2edta[i] + '\t' + nohmatchtsd[i] + '\t' + nohtir1edta[i] + '\t' + nohtir2edta[i] +  '\t')
fnohedta.close()

#view file of EDTA results
columnsf = []
codef = []
scaffoldf = []
startj = []
startf = []
endj = []
endf = []
edtaf = []
with open('edta_seqe.txt') as i:
    linesf = i.readlines()

for j in range(len(linesf)):
    edtaseqf = linesf[j]
    if j == 0:
        columnsf.append(edtaseqf)
    else:
        spf = [p for p, marker in enumerate(edtaseqf) if marker == "\t"]
        codef.append(edtaseqf[spf[0]+1:spf[1]])
        scaffoldf.append(edtaseqf[spf[1]+1:spf[2]])
        startf.append(edtaseqf[spf[4]+1:spf[5]])
        endf.append(edtaseqf[spf[5]+1:spf[6]])
        edtaf.append(edtaseqf[spf[0]+1:-1])

filter1gff = []
filter1code = []
filter1scaffold = []
filter1start = []
filter1end = []
filter1edta = []
filter2gff = []
filter2code = []
filter2scaffold = []
filter2start = []
filter2end = []
filter2edta = []
filtergff = []
filtercode = []
filterscaffold = []
filterstart = []
filterend = []
filteredta = []

i = 0
while i <= (len(startf)-1):
    if i >= len(startf)-1:
        i = len(startf)-1
        filter1gff.append(gff[i])
        filter1code.append(codef[i])
        filter1scaffold.append(scaffoldf[i])
        filter1start.append(startf[i])
        filter1end.append(endf[i])
        filter1edta.append(edtaf[i])
        i = i + 2
    elif startf[i] == startf[i+1]:
        if endf[i] > endf[i+1]:
            filter1gff.append(gff[i])
            filter1code.append(codef[i])
            filter1scaffold.append(scaffoldf[i])
            filter1start.append(startf[i])
            filter1end.append(endf[i])
            filter1edta.append(edtaf[i])
            i = i + 2
        else:
            filter1gff.append(gff[i+1])
            filter1code.append(codef[i+1])
            filter1scaffold.append(scaffoldf[i+1])
            filter1start.append(startf[i+1])
            filter1end.append(endf[i+1])
            filter1edta.append(edtaf[i+1])
            i = i + 2
    else:
        filter1gff.append(gff[i])
        filter1code.append(codef[i])
        filter1scaffold.append(scaffoldf[i])
        filter1start.append(startf[i])
        filter1end.append(endf[i])
        filter1edta.append(edtaf[i])
        i = i + 1

i = 0
while i <= len(filter1start):
    if i >= len(filter1start)-1:
        i = len(filter1start)-1
        filter2gff.append(filter1gff[i])
        filter2code.append(filter1code[i])
        filter2scaffold.append(filter1scaffold[i])
        filter2start.append(filter1start[i])
        filter2end.append(filter1end[i])
        filter2edta.append(filter1edta[i])
        i = i + 2
    elif filter1end[i] == filter1end[i+1]:
        if filter1start[i] < filter1start[i+1]:
            filter2gff.append(filter1gff[i])
            filter2code.append(filter1code[i])
            filter2scaffold.append(filter1scaffold[i])
            filter2start.append(filter1start[i])
            filter2end.append(filter1end[i])
            filter2edta.append(filter1edta[i])
            i = i + 2
        else:
            filter2gff.append(filter1gff[i+1])
            filter2code.append(filter1code[i+1])
            filter2scaffold.append(filter1scaffold[i+1])
            filter2start.append(filter1start[i+1])
            filter2end.append(filter1end[i+1])
            filter2edta.append(filter1edta[i+1])
            i = i + 2
    else:
        filter2gff.append(filter1gff[i])
        filter2code.append(filter1code[i])
        filter2scaffold.append(filter1scaffold[i])
        filter2start.append(filter1start[i])
        filter2end.append(filter1end[i])
        filter2edta.append(filter1edta[i])
        i = i + 1

#file for small sequences
fsmalledta = open('edta_seqe_small.txt', 'w')
fsmalledta.write('Number\tCode\tScaffold\tProgram\tMainClassification\tStart\tEnd\tScore\tStrand\tPhase\tID\tParent\tName\tClassification\tSubclassification\tSequenceOntology\tIdentity\tMethod\tMotif\tTSD1\tTSD2\tmatchTSD\tTIR1\tTIR2')

#remove small sequences
for i in range(len(filter2start)):
    length = int(filter2end[i])-int(filter2start[i])+1
    if  length > 100:
        filtergff.append(filter2gff[i])
        filtercode.append(filter2code[i])
        filterscaffold.append(filter2scaffold[i])
        filterstart.append(filter2start[i])
        filterend.append(filter2end[i])
        filteredta.append(filter2edta[i])
    if 4 < length <=100:
        fsmalledta.write(str(i + 1) + '\t' + filter2edta[i] + '\n')
fsmalledta.close()

for i in range(len(filterscaffold)):
    startj.append(str(int(filterstart[i]) - 20))
    endj.append(str(int(filterend[i]) + 20))

#negative start values replaced with 1
for i in range(len(startj)):
    if int(startj[i]) < 0:
        startj[i] = 1

ffiltergff = open('edta_gff_filter.txt', 'w')
for i in range(len(filterscaffold)):
   ffiltergff.write(filtergff[i])
ffiltergff.close()

ffilterbed = open('edta_seqe_filter.bed.txt', 'w')
for i in range(len(filterscaffold)):
   ffilterbed.write(filterscaffold[i] + '\t' + str(startj[i]) + '\t' + str(endj[i]) + '\n')
ffilterbed.close()

ffilteredta = open('edta_seqe_filter.txt', 'w')
ffilteredta.write('Number\tCode\tScaffold\tProgram\tMainClassification\tStart\tEnd\tScore\tStrand\tPhase\tID\tParent\tName\tClassification\tSubclassification\tSequenceOntology\tIdentity\tMethod\tMotif\tTSD1\tTSD2\tmatchTSD\tTIR1\tTIR2')
for i in range(len(filterscaffold)):
   ffilteredta.write('\n' + str(i + 1) + '\t' + filteredta[i])
ffilteredta.close()

#remove helitrons of filtered sequences
nohfilteralledta = []
scaffoldnoh = []
startnoh = []
endnoh = []
for i in range(len(filterscaffold)):
    nohseq = filteredta[i]
    spnoh = [p for p, marker in enumerate(nohseq) if marker == "\t"]
    mainclassnoh = nohseq[spnoh[2]+1:spnoh[3]]
    if mainclassnoh != 'helitron':
        scaffoldnoh.append(nohseq[spnoh[0]+1:spnoh[1]])
        startnoh.append(nohseq[spnoh[3]+1:spnoh[4]])
        endnoh.append(nohseq[spnoh[4]+1:spnoh[5]])
        nohfilteralledta.append(nohseq)

fnohfilterbed = open('edta_seqe_filter_noh.bed.txt', 'w')
for i in range(len(startnoh)):
   fnohfilterbed.write(scaffoldnoh[i] + '\t' + str(startnoh[i]) + '\t' + str(endnoh[i]) + '\n')
fnohfilterbed.close()

fnohfilteredta = open('edta_seqe_filter_noh.txt', 'w')
fnohfilteredta.write('Number\tCode\tScaffold\tProgram\tMainClassification\tStart\tEnd\tScore\tStrand\tPhase\tID\tParent\tName\tClassification\tSubclassification\tSequenceOntology\tIdentity\tMethod\tMotif\tTSD1\tTSD2\tmatchTSD\tTIR1\tTIR2')
for i in range(len(nohfilteralledta)):
   fnohfilteredta.write('\n' + str(i+1) + '\t' + nohfilteralledta[i])
fnohfilteredta.close()

end= datetime.now()
print('Run time: {}'.format(end - start))
