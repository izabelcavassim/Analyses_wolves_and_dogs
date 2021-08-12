# Workflow for the pyrho analyses of the wolves data

from gwf import Workflow
import glob
import os
import os.path
import itertools
#import pandas as pd

gwf = Workflow()
#bcftools query -l input.vcf

def merge_vcf_files(directory, species, size, outdir):
	inputs = []
	outputs=[f'{outdir}/{species}/BiSNP_6_bespokefixed_v4_minRGQ1minGQ20minDP10maxHet0.99missing2_5_v4_mergeGaTKfiltered_varnonvar_{size}{species}_jointcalled_allchr.vcf.gz']
	options = {
	'memory':'8g',
	'cores':8,
	'walltime':'23:59:59',
	}
	spec = f'''
	/u/local/Modules/default/init/modules.sh
	. "/u/local/apps/anaconda3/etc/profile.d/conda.sh"
	conda activate /u/home/m/mica20/.conda/envs/noncoding

	bcftools concat /u/scratch/c/cad17/wolf_recomb/dog_genomes_clares_cassdata_final_filtered_vcf_justSNPs_6files/{species}/BiSNP_6_bespokefixed_v4_minRGQ1minGQ20minDP10maxHet0.99missing2_5_v4_mergeGaTKfiltered_varnonvar_{size}{species}_jointcalled_chr{{1..38}}.vcf.gz -Oz -o {outdir}/{species}/BiSNP_6_bespokefixed_v4_minRGQ1minGQ20minDP10maxHet0.99missing2_5_v4_mergeGaTKfiltered_varnonvar_{size}{species}_jointcalled_allchr.vcf.gz
	'''
	print(spec)
	return inputs, outputs, options, spec

def demography_file_prep(species, size, inputdir, outdir):
	inputs = [f'{inputdir}/{species}/BiSNP_6_bespokefixed_v4_minRGQ1minGQ20minDP10maxHet0.99missing2_5_v4_mergeGaTKfiltered_varnonvar_{size}{species}_jointcalled_chr{i}.vcf.gz' for i in list(range(1, 38))]
	outputs=[f'{outdir}/{species}/BiSNP_6_bespokefixed_v4_minRGQ1minGQ20minDP10maxHet0.99missing2_5_v4_mergeGaTKfiltered_varnonvar_{size}{species}_jointcalled_chr{i}.smc.gz' for i in list(range(1, 38))]
	options = {
	'memory':'8g',
	'cores':8,
	'walltime':'23:59:59',
	#smc++ vcf2smc {input_vcf} {inputdir}/{species}/BiSNP_6_bespokefixed_v4_minRGQ1minGQ20minDP10maxHet0.99missing2_5_v4_mergeGaTKfiltered_varnonvar_{size}{species}_jointcalled_allchr.smc.gz

	}
	input_vcf = inputs[0]
	spec = f'''
	/u/local/Modules/default/init/modules.sh
	. "/u/local/apps/anaconda3/etc/profile.d/conda.sh"
	conda activate /u/home/m/mica20/.conda/envs/noncoding

	for i in {{1..38}};
    	do smc++ vcf2smc {inputdir}/{species}/BiSNP_6_bespokefixed_v4_minRGQ1minGQ20minDP10maxHet0.99missing2_5_v4_mergeGaTKfiltered_varnonvar_{size}{species}_jointcalled_chr$i.vcf.gz {outdir}/{species}/BiSNP_6_bespokefixed_v4_minRGQ1minGQ20minDP10maxHet0.99missing2_5_v4_mergeGaTKfiltered_varnonvar_{size}{species}_jointcalled_chr$i.smc.gz chr$i {species}:{species}1,{species}2,{species}3,{species}4,{species}5,{species}6,{species}7,{species}8,{species}9,{species}10,{species}11,{species}12,{species}13,{species}14,{species}15,{species}16,{species}17,{species}18 --ignore-missing
	done

	'''
	print(spec)
	return inputs, outputs, options, spec

def demography_estimation(species, inputdir, size):
	inputs=[f'{inputdir}/{species}/BiSNP_6_bespokefixed_v4_minRGQ1minGQ20minDP10maxHet0.99missing2_5_v4_mergeGaTKfiltered_varnonvar_{size}{species}_jointcalled_chr{i}.smc.gz' for i in list(range(1, 38))]
	outputs=[f'/u/home/m/mica20/project-kirk-bigdata/canid_project/scripts/output_{species}/model.final.json', f'/u/home/m/mica20/project-kirk-bigdata/canid_project/scripts/plot_{species}.pdf',f'/u/home/m/mica20/project-kirk-bigdata/canid_project/scripts/plo\
t_{species}.csv' ]
	options = {
	'memory':'8g',
	'cores':8,
	'walltime':'23:59:59',
	}
	spec = f'''
	/u/local/Modules/default/init/modules.sh
	. "/u/local/apps/anaconda3/etc/profile.d/conda.sh"
	conda activate /u/home/m/mica20/.conda/envs/noncoding

	smc++ estimate 4.5e-9 -o output_{species}/ {inputdir}/{species}/BiSNP_6_bespokefixed_v4_minRGQ1minGQ20minDP10maxHet0.99missing2_5_v4_mergeGaTKfiltered_varnonvar_{size}{species}_jointcalled_chr*
	
	smc++ plot plot_{species}.pdf output_{species}/model.final.json -c 
	'''
	print(spec)
	return inputs, outputs, options, spec


def vcf_stats(species, dir_files, sample):
	inputs = [f'{dir_files}/BiSNP_6_bespokefixed_v4_minRGQ1minGQ20minDP10maxHet0.99missing3_5_v4_mergeGaTKfiltered_varnonvar_All{sample}{species}_jointcalled_allchr.vcf.gz']
	outputs = [f'{sample}{species}_file.vchk', f'{sample}{species}_outdir/plot.py']
	options = {
	'memory':'8g',
	'cores':8,
	'walltime':'23:59:59',
	}
	spec = f'''
	/u/local/Modules/default/init/modules.sh
	. "/u/local/apps/anaconda3/etc/profile.d/conda.sh"
	conda activate /u/home/m/mica20/.conda/envs/noncoding


	bcftools stats {dir_files}/BiSNP_6_bespokefixed_v4_minRGQ1minGQ20minDP10maxHet0.99missing3_5_v4_mergeGaTKfiltered_varnonvar_All{sample}{species}_jointcalled_allchr.vcf.gz -s - > {sample}{species}_file.vchk
	plot-vcfstats -p {sample}{species}_outdir {sample}{species}_file.vchk
	'''
	print(spec)
	return inputs, outputs, options, spec

def spliting_vcf_files(species, dir_files, chrom, sample):
	''' Spliting VCF files per chromosome to run with pyrho '''
	inputs = [f'{dir_files}/{species}/BiSNP_6_bespokefixed_v4_minRGQ1minGQ20minDP10maxHet0.99missing2_5_v4_mergeGaTKfiltered_varnonvar_{sample}{species}_jointcalled_allchr.vcf.gz']
	outputs=[f'{dir_files}/{species}/BiSNP_6_bespokefixed_v4_minRGQ1minGQ20minDP10maxHet0.99missing2_5_v4_mergeGaTKfiltered_varnonvar_{sample}{species}_jointcalled_chr{chrom}.vcf.gz', f'{dir_files}/{species}/BiSNP_6_bespokefixed_v4_minRGQ1minGQ20minDP10maxHet0.99missing2_5_v4_mergeGaTKfiltered_varnonvar_{sample}{species}_jointcalled_chr{chrom}.vcf.gz.tbi']
	
	print(outputs)
	options = {
	'memory':'8g',
	'cores':8,
	'walltime':'23:59:59',
	}
	#bcftools view -S {individuals_file_directory}{species}_indivs.txt {inputs[0]} --force-samples | bgzip -c > {dir_files}/chroms/{species}/{sample}{species}_chr{chrom}.vcf.gz
	#vcftools  --gzvcf  {inputs[0]}  --chr chr{chrom}  --recode --recode-INFO-all --out  {dir_files}/chroms/{species}/{sample}{species}_chr{chrom}.vcf
	#bcftools view -r 1 file.vcf.gz

	individuals_file_directory = '/u/home/m/mica20/project-kirk-bigdata/canid_project/aw_pg_vcfs_and_demogs/'
	spec = f'''
	/u/local/Modules/default/init/modules.sh
	. "/u/local/apps/anaconda3/etc/profile.d/conda.sh"
	conda activate /u/home/m/mica20/.conda/envs/noncoding

	bcftools view -r chr{chrom} {inputs[0]} | bgzip -c > {dir_files}/{species}/BiSNP_6_bespokefixed_v4_minRGQ1minGQ20minDP10maxHet0.99missing2_5_v4_mergeGaTKfiltered_varnonvar_{sample}{species}_jointcalled_chr{chrom}.vcf.gz
	tabix -p vcf {dir_files}/{species}/BiSNP_6_bespokefixed_v4_minRGQ1minGQ20minDP10maxHet0.99missing2_5_v4_mergeGaTKfiltered_varnonvar_{sample}{species}_jointcalled_chr{chrom}.vcf.gz
	'''
	print(spec)
	return inputs, outputs, options, spec


def make_lookup_table(n, N, mu, smcpp_file, species, chrom, dir_files, results_dir):
	''' Lookuptable for pyrho ''' 
	inputs = [f'{dir_files}/Results/{smcpp_file}']
	outputs=[f'{results_dir}/logfile_{species}_n_{n}_N_{N}_lookuptable_{chrom}.txt', f'{results_dir}/{chrom}_{species}_n_{n}_N_{N}_lookuptable.hdf']
	options = {
	'memory':'50G',
	'cores':8,
	'walltime':'23:59:59',
	}
	print(outputs)
	spec = f'''
	/u/local/Modules/default/init/modules.sh
	. "/u/local/apps/anaconda3/etc/profile.d/conda.sh"
	conda activate /u/home/m/mica20/miniconda3/envs/canid_rec
	
	pyrho make_table -n {n} -N {N} --mu {mu} --logfile {results_dir}/logfile_{species}_n_{n}_N_{N}_lookuptable_{chrom}.txt --outfile {results_dir}/{chrom}_{species}_n_{n}_N_{N}_lookuptable.hdf \
	--approx --smcpp_file /u/home/m/mica20/project-kirk-bigdata/canid_project/scripts/Files_dog_recombination/{smcpp_file} --decimate_rel_tol 0.1
	'''.format(n=n, N=N, mu=mu, smcpp_file=smcpp_file, chrom=chrom, species=species, results_dir=results_dir)
	#print(spec)
	return inputs, outputs, options, spec


def get_hyperparam(n, N, mu, smcpp_file, species, chrom, sims, dir_files, results_dir):
	'''  Getting the hyperparam ''' 
	inputs = [f'{results_dir}/{chrom}_{species}_n_{n}_N_{N}_lookuptable.hdf', f'{dir_files}/Results/{smcpp_file}']
	outputs=[f'{results_dir}/{chrom}_{species}_n_{n}_N_{N}_hyperparam_results.txt']
	options = {
	'memory':'8g',
	'cores':1,
	'walltime':'23:59:59',
	}
	spec = f'''
	/u/local/Modules/default/init/modules.sh
	. "/u/local/apps/anaconda3/etc/profile.d/conda.sh"
	conda activate /u/home/m/mica20/miniconda3/envs/canid_rec

	pyrho hyperparam -n {n} --mu {mu} --blockpenalty 50,100 \
	--windowsize 25,50 --logfile {results_dir}/logfile_{species}_n_{n}_N_{N}_hyperparameter.txt  --tablefile {results_dir}/{chrom}_{species}_n_{n}_N_{N}_lookuptable.hdf \
	--num_sims {sims} \
	--smcpp_file {dir_files}/Results/{smcpp_file} --outfile {results_dir}/{chrom}_{species}_n_{n}_N_{N}_hyperparam_results.txt 
	'''.format(n=n, N=N, mu=mu, smcpp_file=smcpp_file, chrom=chrom, sims=sims, results_dir=results_dir)
	#print(spec)
	return inputs, outputs, options, spec

def optimize(n, N, mu, species, chrom, vcf, blockpenalty, windowsize, results_dir, dir_files):
	'''  Optmizing parameters (blockpenalty and windowsize) ''' 
	inputs = [f'{results_dir}/{chrom}_{species}_n_{n}_N_{N}_lookuptable.hdf', f'{vcf}']
	outputs=[f'{results_dir}/logfile_{species}_n_{n}_N_{N}_optimize_{chrom}.txt', f'{results_dir}/{chrom}_{species}_n_{n}_N_{N}.rmap']
	options = {
	'memory':'8g',
	'cores':1,
	'walltime':'23:59:59',
	}
	spec = f'''
	/u/local/Modules/default/init/modules.sh
	. "/u/local/apps/anaconda3/etc/profile.d/conda.sh"
	conda activate /u/home/m/mica20/miniconda3/envs/canid_rec

	pyrho optimize --tablefile {results_dir}/{chrom}_{species}_n_{n}_N_{N}_lookuptable.hdf \
	--logfile {results_dir}/logfile_{species}_n_{n}_N_{N}_optimize_{chrom}.txt \
	--vcffile {vcf} \
	--ploidy 2 \
	--outfile {results_dir}/{chrom}_{species}_n_{n}_N_{N}.rmap \
	--blockpenalty {blockpenalty} --windowsize {windowsize} \
	'''.format(n=n, N=N, mu=mu, chrom=chrom, blockpenalty=blockpenalty, windowsize=windowsize, results_dir=results_dir)
	print(spec)
	return inputs, outputs, options, spec

def compute_r2_distribution(n, lookuptable, species, chrom, results_dir, mafcut):
	''' Computing r '''
	inputs = [f'{results_dir}{chrom}_{species}_n_{n}_N_{N}_lookuptable.hdf']
	outputs=[f'{results_dir}{chrom}_{species}_expected_r2_{mafcut}.txt']
	options = {
	'memory':'8g',
	'cores':1,
	'walltime':'23:59:59',
	}
	spec = f'''
	/u/local/Modules/default/init/modules.sh
	. "/u/local/apps/anaconda3/etc/profile.d/conda.sh"
	conda activate /u/home/m/mica20/miniconda3/envs/canid_rec

	pyrho compute_r2 --quantiles .25,.5,.75 --compute_mean --samplesize {n} \
		--tablefile {results_dir}{lookuptable} \
		--outfile {results_dir}{chrom}_{species}_expected_r2_{mafcut}.txt
		--MAFcut {mafcut}
	'''.format(n=n, lookuptable=lookuptable, chrom=chrom, species=species)
	#print(spec)
	return inputs, outputs, options, spec

def converting_vcf_to_ldhat(chrm, sample, dir_files, species):
	''' Coverting vcf to ldhat '''
	inputs = [f'{dir_files}/{species}/BiSNP_6_bespokefixed_v4_minRGQ1minGQ20minDP10maxHet0.99missing2_5_v4_mergeGaTKfiltered_varnonvar_{sample}{species}_jointcalled_allchr.vcf.gz']
	outputs=[f'{dir_files}/chroms/LDhat/{species}/{sample}{species}_chr{chrm}.vcf.ldhat.sites', f'{dir_files}/chroms/LDhat/{species}/{sample}{species}_chr{chrm}.vcf.ldhat.locs']
	options = {
	'memory':'8g',
	'cores':1,
	'walltime':'23:59:59',
	}
	vcf_file = inputs[0]
	spec = f'''
	vcftools --gzvcf {vcf_file} --chr chr{chrm} --ldhat-geno --out {dir_files}/chroms/LDhat/{species}/{sample}{species}_chr{chrm}.vcf
	'''
	#print(spec)
	return inputs, outputs, options, spec

def compute_ldhat_lookuptable(rhomax, n, n_pts, theta, dir_files, sp):
	''' Computing LDhat lookup table '''
	inputs = []
	outputs=[f'{dir_files}/lookup_n{n}_theta_{theta}_{sp}.txt']
	options = {
	'memory':'50g',
	'cores':8,
	'walltime':'23:59:59',
	}

	spec = f'''
	### Make lookup table
	# run complete function to run a new table
	/u/home/m/mica20/project-kirk-bigdata/canid_project/software/LDhat/complete -n {n} -rhomax {rhomax} -n_pts {n_pts} -theta {theta}  
	# rename output file
	mv new_lk.txt {dir_files}/lookup_n{n}_theta_{theta}_{sp}.txt
	'''
	print(spec)
	return inputs, outputs, options, spec

def compute_rho_map(chrm, dir_files, sample, theta, species, n):
	inputs=[f'{dir_files}/chroms/LDhat/{species}/{sample}{species}_chr{chrm}.vcf.ldhat.sites', f'{dir_files}/chroms/LDhat/{species}/{sample}{species}_chr{chrm}.vcf.ldhat.locs'] #, f'{dir_files}/lookup_n{n}_theta_{theta}_{sp}.txt']
	outputs=[f'{dir_files}/Results/LDhat/{species}_{chrm}_acceptance_rates.txt', f'{dir_files}/Results/LDhat/{species}_{chrm}_hotspots.txt', f'{dir_files}/Results/LDhat/{species}_{chrm}_rates.txt', f'{dir_files}/Results/LDhat/{species}_{chrm}_summary.txt']
	options = {
	'memory':'8g',
	'cores':1,
	'walltime':'23:59:59',
	}

	spec = f'''
	
	###Make rho map
	# run rhomap
	/u/home/m/mica20/project-kirk-bigdata/canid_project/software/LDhat/rhomap -seq {dir_files}/chroms/LDhat/{species}/{sample}{species}_chr{chrm}.vcf.ldhat.sites -loc {dir_files}/chroms/LDhat/{species}/{sample}{species}_chr{chrm}.vcf.ldhat.locs -prefix {dir_files}/Results/LDhat/{species}_{chrm}_ -lk {dir_files}/lookup_n{n}_theta_{theta}_{sp}.txt -burn 100000 -its 10000000 -samp 3500
	'''
	print(spec)
	return inputs, outputs, options, spec


def compute_interval_map(chrm, dir_files, sample, theta, species, n):
	inputs=[f'{dir_files}/chroms/LDhat/{species}/{sample}{species}_chr{chrm}.vcf.ldhat.sites', f'{dir_files}/chroms/LDhat/{species}/{sample}{species}_chr{chrm}.vcf.ldhat.locs'] #, f'{dir_files}/lookup_n{n}_theta_{theta}_{sp}.txt']
	outputs=[f'{dir_files}/Results/LDhat/{species}_{chrm}_interval_rates.txt', f'{dir_files}/Results/LDhat/{species}_{chrm}_interval_bounds.txt']
	print(outputs)
	options = {
	'memory':'8g',
	'cores':1,
	'walltime':'50:59:59',
	}

	spec = f'''
	
	###Make rho map
	# run rhomap
	/u/home/m/mica20/project-kirk-bigdata/canid_project/software/LDhat/interval -seq {dir_files}/chroms/LDhat/{species}/{sample}{species}_chr{chrm}.vcf.ldhat.sites -loc {dir_files}/chroms/LDhat/{species}/{sample}{species}_chr{chrm}.vcf.ldhat.locs -prefix {dir_files}/Results/LDhat/{species}_{chrm}_interval_ -lk {dir_files}/lookup_n{n}_theta_{theta}_{sp}.txt -burn 100000 -its 10000000 -samp 3500
 	'''
	print(spec)
	return inputs, outputs, options, spec


def compute_stats(chrm, dir_files, sample, species):
	inputs =[f'{dir_files}/Results/LDhat/{species}_{chrm}_interval_rates.txt', f'{dir_files}/Results/LDhat/{species}_{chrm}_interval_bounds.txt']
	outputs = [f'{dir_files}/Results/LDhat/{species}_{chrm}_res.txt']
	options = {
	'memory':'8g',
	'cores':1,
	'walltime':'23:59:59',
	}
	spec = f'''
 	# run the summary stats
 	/u/home/m/mica20/project-kirk-bigdata/canid_project/software/LDhat/stat -input {dir_files}/Results/LDhat/{species}_{chrm}_interval_rates.txt -loc {dir_files}/chroms/LDhat/{species}/{sample}{species}_chr{chrm}.vcf.ldhat.locs -prefix {dir_files}/Results/LDhat/{species}_{chrm}_
	'''		
	print(spec)
	return inputs, outputs, options, spec

# Mutation rate in wolves 4.5 x 10-9
# Sample size in wolves 14, pg = 15
# aw_chrm38_mu4.5e-9_smcplot.csv


# Making a dictionary of breeds and sample sizes
metadata = {"AW":15, "BC":7, "EW":10, "IR": 10, "LB":10, "PG":15, "TM":9} # "AW" is 15, but is actually 14
metadata = {"BC":"Ten", "EW":"Ten", "IR":"Ten", "LB":"Ten", "TM":"Ten", "AW":"Fifteen", "PG":"Fifteen"} #, "MD":"group", "MW":"group"} # "AW" is 15, but is actually 14

#/u/scratch/c/cad17/wolf_recomb/dog_genomes_clares_cassdata_final_filtered_vcf_justSNPs_6files/MD/ ##### Missing ################################################
#/u/scratch/c/cad17/wolf_recomb/dog_genomes_clares_cassdata_final_filtered_vcf_justSNPs_6files/MW/ ##### Missing ################################################

dir_files = "/u/home/m/mica20/project-kirk-bigdata/canid_project/aw_pg_vcfs_and_demogs"
results_dir = "/u/scratch/c/cad17/wolf_recomb/Results"

dir_files = "/u/scratch/c/cad17/wolf_recomb"
chroms = range(1,39,1)
species = ["AW", "PG"]
metadata = {"AW":15, "BC":7, "EW":10, "IR": 10, "LB":10, "PG":15, "TM":9} # "AW" is 15, but is actually 14
metadata = {"BC":"Ten", "EW":"Ten", "IR":"Ten", "LB":"Ten", "TM":"Ten", "AW":"Fifteen", "PG":"Fifteen"} #, "MD":"group", "MW":"group"} # "AW" is 15, but is actually 14


#####################################################
# Subsetting files                                  #
#####################################################
dir_files = "/u/scratch/c/cad17/wolf_recomb/dog_genomes_clares_cassdata_final_filtered_vcf_justSNPs_6files"
for sp in species:
	for ch in chroms:
		gwf.target_from_template(f'Subsetting_vcf_{sp}_{ch}', spliting_vcf_files(dir_files=dir_files, chrom=ch, species=sp, sample=metadata[sp]))


######################################################
# COMPUTING DEMOGRAPHY.                              #
######################################################

smcpp = True
if smcpp == True:
	species = metadata.keys()
	vcf_dir =  "/u/scratch/c/cad17/wolf_recomb/dog_genomes_clares_cassdata_final_filtered_vcf_justSNPs_6files"
	out_dir = "/u/scratch/c/cad17/wolf_recomb/chroms/CONCATENATED"
	for sp in species:
		print(sp)
		vcf_dir_sp = vcf_dir
		size = metadata[sp]
		#gwf.target_from_template(f'Merging_vcf_{sp}', merge_vcf_files(directory=vcf_dir_sp, species=sp, size=size, outdir=out_dir))
		spliting_vcf_files
		gwf.target_from_template(f'Making_smcpp_input_{sp}', demography_file_prep(species=sp, size=size, inputdir=vcf_dir_sp, outdir=out_dir))
		gwf.target_from_template(f'Estimating_demography_{sp}',demography_estimation(species=sp, size=size, inputdir=out_dir))


######################################################
# COMPUTING Recombination                            #
######################################################


## !!!! here n is fixed, we need to change here once the other dog breeds are included 
gwf.target_from_template(f'Producing_lookuptable_AW_ldhat', compute_ldhat_lookuptable(n=15*2, rhomax=100, n_pts=101, theta=0.002, dir_files=dir_files, sp="AW"))
gwf.target_from_template(f'Producing_lookuptable_PG_ldhat', compute_ldhat_lookuptable(n=15*2, rhomax=100, n_pts=101, theta=0.002, dir_files=dir_files, sp="PG"))


## Estimating stats for each vcf file
#gwf.target_from_template(f'Producing_VCF_stats_AW', vcf_stats(dir_files=dir_files, sample=15, species='AW'))
#gwf.target_from_template(f'Producing_VCF_stats_PG', vcf_stats(dir_files=dir_files, sample=15, species='PG'))

species = ["AW", "PG"]
#for loop here
for sp in species:
	for ch in chroms:
		gwf.target_from_template(f'Subsetting_vcf_{sp}_{ch}_lhdat', converting_vcf_to_ldhat(dir_files=dir_files, chrm=ch, species=sp, sample=metadata[sp]))
		#gwf.target_from_template(f'Computing_rho_chr_{sp}_{ch}_ldhat', compute_rho_map(chrm=ch, dir_files=dir_files, sample=metadata[sp], theta=0.002, species=sp, n=15*2))
		gwf.target_from_template(f'Computing_rho_interval_chr_{sp}_{ch}_ldhat', compute_interval_map(chrm=ch, dir_files=dir_files, sample=metadata[sp], theta=0.002, species=sp, n=15*2))
		gwf.target_from_template(f'Computing_interval_stats_{sp}_{ch}_ldhat', compute_stats(chrm=ch, dir_files=dir_files, sample=metadata[sp], species=sp))	


