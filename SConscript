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



import sfle
import sfleoo
import sys
import os

Import( "*" )


oo = sfleoo.ooSfle(  fileDirOutput = fileDirOutput, fileDirTmp = fileDirTmp, fileDirInput=fileDirInput )
#--Constants--#
iReads = 5000000
iReadLen = 100
zipModel = os.path.abspath(oo.fin("ill100v5_s.gzip") ) # Illumina Model

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

afnaGenomes = []

for fnaGenomeStem in afnaGenomeStems:
	fnaGenome = os.path.abspath(sfle.d(fileDirInput,"genomes",fnaGenomeStem))
	afnaGenomes.append(fnaGenome)

afnaGenomes.sort()

astrRuns = ["05","10","25"]
c_strCLEANGENOMES = "Y"

dirWD = os.getcwd()

RemovePosGenes = oo.fsrc("RemovePosGenes.py")

#Please only run A) Section 1 alone, or B) Everything but Section 1.
#The reason for this is that GemSim takes the bacterial genomes as a folder, and not an array of files.
#I have not yet figured out how to tell the GemReads process to wait for the files in the folder to be updated
#before running.

for strDB in ["ARDB","VF"]:
	faaInputProts = os.path.abspath(sfle.d(fileDirInput,strDB+".input.faa"))

	afnaCleanGenomes = []
	afnaTmpGenomes = []


	for fnaGenome in afnaGenomes:
		fnaGenomeStem = os.path.basename(fnaGenome)

		#--TmpFiles --#
		dirCleanGenomes = sfle.d(fileDirTmp,strDB,"clean_genomes")
		fnaCleanGenome = sfle.d(dirCleanGenomes,fnaGenomeStem)
		afnaCleanGenomes.append(fnaCleanGenome)

		dirTmpGenomes = sfle.d(fileDirTmp,strDB,"tmp_genomes")
		dirTmpProcessing = sfle.d(dirTmpGenomes,"Process" + fnaGenomeStem.replace(".genome",""))
		fnaTmpGenome = sfle.d(dirTmpGenomes,fnaGenomeStem)
		afnaTmpGenomes.append(fnaTmpGenome)


	    #python /n/data/users/jkaminski/sfle/input/metagenome/src/RemovePosGenes.py --prots ARDB.input.faa --genome c.difficile.genome --tmpgenome tmp.fna --clean testclean.fna
		oo.pipe([fnaGenome,faaInputProts],[fnaCleanGenome,fnaTmpGenome],RemovePosGenes,prots=faaInputProts, genome=fnaGenome,tmpgenome=fnaTmpGenome,clean=fnaCleanGenome,tmp=dirTmpProcessing,id=.80)
		Default(fnaCleanGenome)



	afnaCleanGenomes.sort()
	#print "There are", len(afnaCleanGenomes), " genomes."
	Default(afnaCleanGenomes)


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


		txtGenomes = os.path.abspath(oo.fin("GenomeNames.txt"))

		if (c_strCLEANGENOMES == "Y"):
			fnaPadGenome = sfle.d(dirCleanGenomes,"b.longum.genome")
			dirGenomes = sfle.d(sfle.d(fileDirInput,"genomes"))

			# Please change this if you change dirCleanGenomes. GemSim will not take this
			# as absolute path, so I use the path relative to sfle.
			dirGemSimRef   = sfle.d("output","metagenome","tmp",strDB,"clean_genomes")

			#print "DIRECTORY TESTING"
			#print(dirWD)
			#print(str(dirCleanGenomes))
			#print(dirGemSimRef)


		else:
			fnaPadGenome = os.path.abspath(oo.fin("genomes" + os.sep + "b.longum.genome"))
			dirGenomes = sfle.d("input","metagenome","input","genomes")
			dirGemSimRef   = sfle.d("output","metagenome","tmp","GSRefGenomes")

		file_out2 = oo.ftmp("out.genome")
		strDBRun = strDB + strRun
		dirDBRun = sfle.d(fileDirTmp,strDBRun)
		print "The run is ",strDBRun

		#--Tmp Files--#
		fnaScreened =sfle.d(dirDBRun,strDBRun+ "screened.fna")
		txtScreenLog =sfle.d(dirDBRun,strDBRun+ "screenlog.txt")
		txtGS    =sfle.d(dirDBRun,strDBRun+ "gs.txt")
		stderrSim =sfle.d(dirDBRun,strDBRun+ "sim.log")

		fastaSim =sfle.d(dirDBRun,strDBRun+ ".fasta")
		fastaFinal =sfle.d(dirDBRun,strDBRun+ "sim.fasta")
		txtAbundance =sfle.d(dirDBRun,strDBRun+ "abund.txt")
		txtCount =sfle.d(dirDBRun,strDBRun+ "count.txt")


		# Gem is peculiar about its output. You have to specify a prefix instead of the full output file.
		stubGemSim = sfle.d(dirDBRun,strDBRun)
		fastqSim =sfle.d(dirDBRun,strDBRun+ "_single.fastq")
		print stubGemSim
		print fastqSim

	    #---Parameters---#


		#Screen out nucleotide seqs that do not have a corresponding protein sequwnce
		oo.pipe([fnaNucs,faaInputProts],[fnaScreened,txtScreenLog],ScreenNucs,out=fnaScreened,nucs=fnaNucs,prots=faaInputProts,log=txtScreenLog)



	    #Make individual fasta files for each gene, make the abundance table and gold standard
		oo.pipe(fnaScreened,[txtAbundance,txtGS], SimpleSim,nucs=fnaScreened,N=iGenes,muS=1,muG=1,gold=txtGS,genomes=txtGenomes,fastadir=dirGemSimRef,
		abund=txtAbundance,padgenome=fnaPadGenome,padlength=100, dirgenomes= os.path.abspath(dirGenomes),pctspike=dSpike)
		Default(txtAbundance)


	    #Incorporate individual nucs into a synthetic metagenome, along  with the other genomes in "input/genomes", using GemReads.py
		oo.pipe([txtAbundance,afnaCleanGenomes],[stderrSim,fastqSim],GemReads,R=dirGemSimRef, n=iReads,l="d", m=zipModel,c="",q=64,o=stubGemSim,a = txtAbundance)
		Default(fastqSim)

		oo.pipe(fastqSim,fastaSim,Fastq2Fasta)

	    #Cut out excess text, reduce gene lables for spiked genes to ">USR_NAME_END"
		oo.pipe(fastaSim,fastaFinal,"sed",e="s/^>.*\(USR_.*_END\).*$/>\\1/g")
	 	Default(fastaFinal)

		oo.ex(fastaFinal,txtCount,"grep",args=[("-e",">"),("-c","")],outpipe= True)
		#oo.ex(fastaFinal,txtCount,"grep",e="\">\"",c="",verbose=True)
		Default(txtCount)
