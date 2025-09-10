import os
import numpy as np
from scipy.stats import chisquare
import matplotlib.ticker as mticker
import matplotlib.pyplot as plt


####### get overall stats for the graph   ######
####### create the Alleles.xls file that contains billalelic SNPs with read counts for each allele
####### includes results of chi-square contingency tests
####### input is P1.vcf


fout = open('Alleles.xls','w')
fn = open('raw','w')


listFreq = []
finalD = {'Diploid': 0, 'Triploid': 0,'Tetra_Penta': 0}
m = 0
liste = []
dico = {}
with open("P1.vcf", "r") as f:
    for line in f:
        if line[0] == '#':
            try:
                dico[line.split('=<ID=')[1].split(',')[0]] = int(line.split('length=')[1].split('>')[0])
            except:
                pass
        else:
            if 'INDEL' not in line:
                depth = float(line.split('DP=')[1].split(';')[0])
                MQ = float(line.split('MQ=')[1].split('\t')[0])
                if depth >= 60.0 and MQ >= 30.0:
                    GENO = line.split('\t')[-1].replace('\n','')
                    ref, alt = int(GENO.split(':')[-1].split(',')[0]), int(GENO.split(':')[-1].split(',')[1])
                    #MQ =  line.split('MQ=')[1].split('\t')[0]
                    qual =line.split('\t')[5]
                    DP = GENO.split(':')[-2]
                    if float(qual) >=30.0:
                        if GENO.split(':')[0] == '0/1':
                            if dico[line.split('\t')[0]] >= 1000:  #size of the scaffolds considered
                                m = m + 1
                                if alt >= ref :
                                    a = alt
                                    b = ref
                                else:
                                    a = ref
                                    b = alt

                                #freqMinor = float(a)/float(a+b)
                                listFreq.append(float(b)/float(a+b))
                                #listFreq.append(float(b) / float(a + b))

                                observed = [a,b]
                                expected = [float(a+b)/2, float(a+b)/2]
                                chi2_a, p_a = chisquare(f_obs=observed, f_exp=expected)

                                expected = [(float(a + b) / 4) * 3, float(a + b) / 4]
                                chi2_b, p_b = chisquare(f_obs=observed, f_exp=expected)

                                expected = [(float(a + b) / 5) * 4, float(a + b) / 5]
                                chi2_c, p_c = chisquare(f_obs=observed, f_exp=expected)

                                if p_a >= p_b and p_a > p_c :
                                    win = 'Diploid'
                                elif p_b > p_a or p_c >= p_a:
                                    win = 'Tetra_Penta'


                                finalD[win] = finalD[win] + 1
                                liste.append((line.split('\t')[0],p_a,p_b,p_c))
                                fout.write(line.split('\t')[0] + '\t' + line.split('\t')[1] + '\t' +  str(qual) + '\t' + str(MQ) + '\t' +
                                           str(DP) + '\t' + str(a) + '\t' + str(b) + '\t' + str(p_a) + '\t' + str(p_b) + '\t' + str(p_c) + '\t' + win + '\n')
                                fn.write(line)

fout.close()
fn.close()
print(finalD, m)
print('\n\n\n\n')



#####  generate histogram of allele frequencies for the entire genome  #####
fig = plt.figure(figsize=(6, 3))
ax1 = plt.subplot2grid((1,1), (0, 0),) #,colspan = 2)
ax1.hist(listFreq, bins=50, edgecolor='black',color = '#7f7f7f', alpha=0.85)  # bins = number of bars
ax1.set_xticks([0.0,0.1,0.2,0.3,0.4, 0.5])
#ax1.axvline(0.5, color = 'r',linestyle='--', linewidth=0.75)
ax1.set_ylabel('Number of SNP loci')
ax1.set_xlabel('Minor allele frequency')
plt.tight_layout()
plt.savefig('Figure_1.jpg',dpi=300,  format='jpg')



####### analyse only scaffolds > 100000 bp - generate part B of Figure 1  ######
####### generate a graph output for each scaffold > 100kbp

dicoS = {}
with open("P1.vcf", "r") as f:
    for line in f:
        if line[0] == '#':
            try:
                dicoS[line.split('=<ID=')[1].split(',')[0]] = int(line.split('length=')[1].split('>')[0])
            except:
                pass

fx = open('Alleles.xls')
ct = fx.readlines()
fx.close()

dico = {}
lscaf = []
for x in ct:
    if x.split('\t')[0] not in lscaf:
        lscaf.append(x.split('\t')[0])
    dico[x.split('\t')[0] + '|' + x.split('\t')[1]] = x

for sc in lscaf:

        if dicoS[sc] >= 100000:
            lsobjective = []
            j = 0
            while j < dicoS[sc]:
                if j % 1000 == 0:
                    lsobjective.append(j)
                j = j + 1
            lsx,lsy, lsz, lsw,pl1, pl2 = [], [], [], [],[], []
            os.system('samtools depth -r ' + sc + ' ../P1.masurcaV4.polishing/map1s.bam > ' + sc + '.depth')
            fout = open(sc + '.xls','w')
            ln = []
            with open(sc + ".depth", "r") as f:
                for line in f:
                    try:
                        info = dico[line.split('\t')[0] + '|' + line.split('\t')[1]]
                        val = float(info.split('\t')[6]) / float(info.split('\t')[4])
                        mn, mj = float(info.split('\t')[6]), float(info.split('\t')[5])
                        fout.write(line.replace('\n','') + '\t' + str(val) + '\n')
                        ln.append((line.split('\t')[1],line.split('\t')[2],val,mn,mj))
                    except:
                        fout.write(line)
                        ln.append((line.split('\t')[1],line.split('\t')[2]))

            fout.write('\n\n\n')
            fout.write('Location\tmean depth\tnb of SNPs\tMean nb reads minor allele\tChi2 1:1\tChi2 1:3\tChi2 1:4\n')
            i = 1
            dp, allele, alleleD, alleleM = [], [], [], []
            while i < len(ln):
                dp.append(int(ln[i][1]))
                try:
                    allele.append(float(ln[i][2]))
                    alleleD.append(float(ln[i][3]))
                    alleleM.append(float(ln[i][4]))
                except:
                    pass
                if i in lsobjective:
                    lsobjective.remove(i)
                    #print('*',i,allele,np.mean(allele))
                    if allele != []:
                        fq = np.mean(allele) * 100
                        minn = np.mean(alleleM)
                        maj = np.mean(alleleD)
                        observed = [minn, maj]
                        expected = [(minn + maj)/2, (minn + maj)/2]
                        chi2_a, p_a = chisquare(f_obs=observed, f_exp=expected)
                        if p_a >= 0.05:
                            pl1.append(0.6)
                        else:
                            pl1.append(np.nan)
                        unit = (minn + maj)/4
                        expected = [3*unit, unit]
                        chi2_c, p_c = chisquare(f_obs=observed, f_exp=expected)
                        if p_c >= 0.05:
                            pl2.append(0.6)
                        else:
                            pl2.append(np.nan)
                    else:
                        pl1.append(np.nan)
                        pl2.append(np.nan)
                        p_a = p_c = 'nan'
                    lsx.append(i)
                    if len(allele) != 0:
                        lsy.append(np.mean(allele))
                    else:
                        lsy.append(0.0)
                    lsz.append(np.mean(dp))
                    lsw.append(len(allele))
                    fout.write(str(i) + '\t' + str(np.mean(dp)) + '\t' + str(len(allele)) + '\t' + str(np.mean(allele)) + '\t' + str(p_a) + '\n')
                    dp, allele, alleleD, alleleM = [], [], [], []
                    i = i - 500
                i = i + 1
            #plotting part
            fig = plt.figure(figsize=(8, 5))
            ax1 = plt.subplot2grid((2, 1), (0, 0), )  # ,colspan = 2)
            ax2 = plt.subplot2grid((2, 1), (1, 0), )  # ,colspan = 2)
            C1 = "#000000"
            C2 = "#7f7f7f"
            C3 = "#000000"
            ax1.plot(lsx,lsy, color=C1)
            ax1.plot(lsx, pl1, marker = '+', linestyle="", color="k",markersize=4)
            ax11 = ax1.twinx()
            ax11.bar(lsx, lsw, width=750, alpha=0.85, color=C2, label="Bars")
            ax1.set_yticks([0.0,0.25,0.5])
            ax1.tick_params(axis="y", colors=C1)
            ax1.set_ylabel('Minor allele frequency',color=C1)
            ax11.set_ylabel('Number of SNPs',color=C2)
            ax11.yaxis.set_major_locator(mticker.MaxNLocator(nbins=4))
            ax11.tick_params(axis="y", colors= C1)
            ax2.plot(lsx,lsz, color=C3)
            ax2.yaxis.set_major_locator(mticker.MaxNLocator(nbins=4))
            ax2.set_ylabel('Reads depth',color=C3)
            ax2.tick_params(axis="y", colors=C3)
            plt.savefig('Figure_'+ sc.split('_')[0] + '_' + sc.split('_')[1] + '.jpg', dpi=300, format='jpg')
            fout.close()

