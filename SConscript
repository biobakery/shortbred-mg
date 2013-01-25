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



for InputDB in ["ARDB","VF"]:

    #<Insert code to create softlinks to my shuffled genomes>

    #--Input Files--#
	fnaNucs  = os.path.abspath(oo.fin(InputDB + "nucs.fna"))
	zipModel = os.path.abspath(oo.fin("ill100v5_s.gzip") )
	txtGenomes = os.path.abspath(oo.fin("GenomeNames.txt"))

    #--Tmp Files--#
	txtGS    = oo.ftmp(InputDB + "gs.txt")
	fnaShuffledGenome	= oo.ftmp("genomes" + os.sep + InputDB+ "shuf.fna")
	stderrSim = oo.ftmp(InputDB + "sim.log")
	fastqSim = oo.ftmp(InputDB + "_single.fastq")
	fastaSim = oo.ftmp(InputDB + ".fasta")
	fastaFinal = oo.ftmp(InputDB + "sim.fasta")
	txtAbundance = (oo.ftmp(InputDB + "abund.txt"))

    #--Programs--#
	SimpleSim        = oo.fsrc("simplesim.py")
	GemReads	      = oo.fsrc("gemsim" + os.sep + "GemReads.py")
	Fastq2Fasta      = oo.fsrc("fastq2fasta.py")

    #---Parameters---#
	stubGemSim = sfle.d(fileDirTmp,InputDB)

    #--Dirs------"
	dirGenomes = sfle.d("output","metagenome","tmp","genomes")


    #Make individual fasta files for each gene, make the abundance table and gold standard
	oo.pipe(fnaNucs,[txtAbundance,txtGS], SimpleSim,nucs=fnaNucs,N=150,muS=1,muG=6.25,gold=txtGS,genomes=txtGenomes,fastadir=dirGenomes,abund=txtAbundance)
 	Default(txtAbundance)


    #Incorporate it into a synthetic metagenome, with the other genomes in "input/genomes", using GemReads.py
	oo.pipe([txtAbundance],[stderrSim,fastqSim],GemReads,R=dirGenomes, n=5000000,l=100, m=zipModel,c="",q=64,o=stubGemSim,a = txtAbundance )

	oo.pipe(fastqSim,fastaSim,Fastq2Fasta)

    #Cut out excess text, reduce gene lables for spiked genes to ">USR_NAME_END"
	oo.pipe(fastaSim,fastaFinal,"sed",e="s/^>.*\(USR_.*_END\).*$/>\\1/g")
 	Default(fastaFinal)


