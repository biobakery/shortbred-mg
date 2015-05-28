#######################################################################################
# This file is provided under the Creative Commons Attribution 3.0 license.
#
# You are free to share, copy, distribute, transmit, or adapt this work
# PROVIDED THAT you attribute the work to the authors listed below.
# For more information, please see the following web page:
# http://creativecommons.org/licenses/by/3.0/
#
# This file is a component of the SflE Scientific workFLow Environment for reproducible
# research, authored by the Huttenhower lab at the Harvard School of Public Health
# (contact Curtis Huttenhower, chuttenh@hsph.harvard.edu).
#
# If you use this environment, the included scripts, or any related code in your work,
# please let us know, sign up for the SflE user's group (sfle-users@googlegroups.com),
# pass along any issues or feedback, and we'll let you know as soon as a formal citation
# is available.
#######################################################################################


#IMPORTANT NOTE:
#   GemReads.py adds "./" to the beginning of the reference folder for its input genomes.
#   I edited the version saved here to *DROP the "./". Without the change, it is difficult
#   to use sfle.
#
#   Line 119 of GemReads.py was:
#           path1='./'+ directory
#   Changed to:
#           path1=directory


import sfle
import sfleoo
import sys
import os
import string

Import( "*" )


oo = sfleoo.ooSfle(  fileDirOutput = fileDirOutput, fileDirTmp = fileDirTmp, fileDirInput=fileDirInput )


#--Programs--#
SimpleSim        = oo.fsrc("simplesim.py")
ScreenNucs      = oo.fsrc("ScreenNucs.py")
GemReads	      = oo.fsrc("gemsim" + os.sep + "GemReads.py")
Fastq2Fasta      = oo.fsrc("fastq2fasta.py")
CPUnix = "cp"

afnaGenomeStems = ["b.thetaiotaomicron.genome", "p.gingivalis.genome", "p.distasonis.genome",
"p.ruminicola.genome", "v.parvula.genome", "b.mallei.genome", "l.acidophilus.genome",
"e.eligens.genome", "s.aureus_n315.genome", "e.coli.genome", "b.longum.genome","f.nucleatum.genome",
"s.pneumoniae.genome", "c.difficile.genome", "n.meningitidis.genome", "c.jejuni.genome","p.multocida.genome",
"f.johnsoniae.genome","l.buccalis.genome","r.mucilaginosa.genome"]
dirGenomes = sfle.d(sfle.d(fileDirInput,"genomes"))
txtGenomes = os.path.abspath(oo.fin("GenomeNames.txt"))

afnaGenomes = []

for fnaGenomeStem in afnaGenomeStems:
	fnaGenome = os.path.abspath(sfle.d(fileDirInput,"genomes",fnaGenomeStem))
	afnaGenomes.append(fnaGenome)

afnaGenomes.sort()

astrRuns = ["05","10","25"]
c_strCLEANGENOMES = "Y"

dirWD = os.getcwd()

RemovePosGenes = oo.fsrc("RemovePosGenes.py")

for strDB in ["ARDB","VF"]:
	faaInputProts = os.path.abspath(sfle.d(fileDirInput,strDB+".input.faa"))

	for strMG in ["Illumina","454_mockYAT","454_deep","small"]:
	#for strMG in ["454_mockYAT"]:

		if strMG == "Illumina":
			iReads = 5000000
			iReadLen = "d" # "d" tells GemSim ("GemReads.py") to use the empirical model for length
			zipModel = os.path.abspath(oo.fin("ill100v5_s.gzip") ) # Illumina Model
			iPadLength = 100
		elif strMG== "454_mockYAT":
			iReads = 155890
			iReadLen = "d"
			zipModel = os.path.abspath(sfle.d(fileDirSrc,"gemsim","models","r454ti_s.gzip")) # 454 Model
			iPadLength = 450
		elif strMG=="454_deep":
			iReads = 1250000
			iReadLen = "d"
			zipModel = os.path.abspath(sfle.d(fileDirSrc,"gemsim","models","r454ti_s.gzip"))  # 454 Model
			iPadLength = 450
		elif strMG=="small":
			iReads = 4000000
			iReadLen = "d"
			zipModel = os.path.abspath(oo.fin("ill100v5_s.gzip") )  # 454 Model
			iPadLength = 100

		for strRun in astrRuns:
			if (strRun == "05"):
				iGenes = 150
				dSpike = int(strRun)/100.0
			if (strRun == "10"):
				iGenes = 500
				dSpike = int(strRun)/100.0
			elif (strRun == "25"):
				iGenes = 1000
				dSpike = int(strRun)/100.0

			afnaCleanGenomes = []
			afnaTmpGenomes = []

			strDBRun = strDB + strRun
			dirDBRun = sfle.d(fileDirTmp,strMG,strDBRun)

			#This is where all of the final padded genes and clean genomes are stored.
			dirGemSimRef   = sfle.d(dirDBRun,"final_genes")

			#Clean the KEGG bacterial genomes of matching ARDB or VF sequences, save the clean version in dirGemSimRef
			for fnaGenome in afnaGenomes:
				fnaGenomeStem = os.path.basename(fnaGenome)

				dirCleanGenomes = dirGemSimRef

				# This is genome that we use to "pad" the spiked sequences.
				fnaPadGenome = sfle.d(dirCleanGenomes,"b.longum.genome")

				fnaCleanGenome = sfle.d(dirCleanGenomes,fnaGenomeStem)
				afnaCleanGenomes.append(fnaCleanGenome)

				dirTmpGenomes = sfle.d(dirDBRun,"tmp_genomes")
				dirTmpProcessing = sfle.d(dirTmpGenomes,"Process" + fnaGenomeStem.replace(".genome",""))
				fnaTmpGenome = sfle.d(dirTmpGenomes,fnaGenomeStem)
				afnaTmpGenomes.append(fnaTmpGenome)

			    #python /n/data/users/jkaminski/sfle/input/metagenome/src/RemovePosGenes.py --prots ARDB.input.faa --genome c.difficile.genome --tmpgenome tmp.fna --clean testclean.fna
				oo.pipe([fnaGenome,faaInputProts],[fnaCleanGenome,fnaTmpGenome],RemovePosGenes,prots=faaInputProts, genome=fnaGenome,tmpgenome=fnaTmpGenome,clean=fnaCleanGenome,tmp=dirTmpProcessing,id=.90)
				Default(fnaCleanGenome)

			afnaCleanGenomes.sort()
			#print "There are", len(afnaCleanGenomes), " genomes."
			Default(afnaCleanGenomes)

			#####################################################################
			#NOTE: there are two different ARDB fna files:
			#ARDBnucs.fna = ~160 real nuc seqs,
			#ARDBbt.fna = 7,825 back translated sequences (using backtranseq from EMBOSS)
			#Please select the one below needed for your project.

			if (strDB == "ARDB" and strRun != "05"):
				fnaNucs = os.path.abspath(oo.fin("ARDBbtnucs.fna"))
			else:
				fnaNucs  = os.path.abspath(oo.fin(strDB + "nucs.fna"))
			######################################################################


            #--Tmp Files--#
			#Step 1 - Remove nuc seqs without corresponding prot.
			fnaScreened =sfle.d(dirDBRun,strDBRun+ "screened.fna")
			txtScreenLog =sfle.d(dirDBRun,strDBRun+ "screenlog.txt")
			#Step 2 - Make the padded fasta files. Save them in dirGemSimRef
			txtAbundance 	= sfle.d(dirDBRun,strDBRun+ "abund.txt")
			txtGS    		= sfle.d(dirDBRun,strDBRun+ "gs.txt")
			txtLogSim = sfle.d(dirGemSimRef,strDBRun+"simplesimlog.txt")
			#Step 3 - Call GemSim to make fastq synthetic MG's
			stderrSim =sfle.d(dirDBRun,strDBRun+ "sim.log")
			stubGemSim = sfle.d(dirDBRun,strDBRun)
			fastqSim =sfle.d(dirDBRun,strDBRun+ "_single.fastq")
         
			#Step 4 - Convert to fasta
			fastaSim =sfle.d(dirDBRun,strDBRun+ ".fasta")
			#Step 5 - Cut out excess text, reduce gene lables for spiked genes to ">USR_NAME_END"
			fastaFinal =sfle.d(dirDBRun,strDBRun+ "sim.fasta")
			#Step 6 - Count reads
			txtCount =sfle.d(dirDBRun,strDBRun+ "count.txt")



			#Step 1: Screen out nucleotide seqs that do not have a corresponding protein sequence
			oo.pipe([fnaNucs,faaInputProts],[fnaScreened,txtScreenLog],ScreenNucs,out=fnaScreened,nucs=fnaNucs,prots=faaInputProts,log=txtScreenLog)
			Default(fnaScreened)

            #Step 2: Padding - Make individual fasta files for each gene, make the abundance table and gold standard
			oo.pipe(fnaScreened,[txtAbundance,txtGS,txtLogSim], SimpleSim,nucs=fnaScreened,N=iGenes,muS=1,muG=1,gold=txtGS,genomes=txtGenomes,fastadir=dirGemSimRef,
			abund=txtAbundance,padgenome=fnaPadGenome,padlength=iPadLength, dirgenomes= os.path.abspath(dirGenomes),pctspike=dSpike,log=txtLogSim)
			Default(txtAbundance,txtGS,txtLogSim)

                    ### NEW STEP ### 

			#Step 3: Incorporate individual nucs into a synthetic metagenome, along  with the other genomes in "input/genomes", using GemReads.py
			oo.pipe([txtAbundance,afnaCleanGenomes],[stderrSim,fastqSim],GemReads,R=dirGemSimRef, n=iReads,l="d", m=zipModel,c="",q=64,o=stubGemSim,a = txtAbundance)
			Default(fastqSim)

			#Step 4: Convert fastq to fasta
			oo.pipe(fastqSim,fastaSim,Fastq2Fasta)

			#Step 5: Cut out excess text, reduce gene lables for spiked genes to ">USR_NAME_END"
			oo.pipe(fastaSim,fastaFinal,"sed",e="s/^>.*\(USR_.*_END\).*$/>\\1/g")
		 	Default(fastaFinal)

			#Step 6: Count up reads.
			oo.ex(fastaFinal,txtCount,"grep",args=[("-e",">"),("-c","")],outpipe= True)
			Default(txtCount)

